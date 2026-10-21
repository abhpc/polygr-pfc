#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "pfc_dcu.h"

#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <sched.h>
#include <unistd.h>
#include <vector>

struct PfcDcu {
	MPI_Comm world, group, comm;
	int world_rank, rank, ranks, W, H, Wc, Wp;
	int group_rank, group_size, local_real_count, local_spectrum_count;
	bool leader;
	bool shared_transport;
	MPI_Win shared_window;
	ptrdiff_t lH, lW;
	size_t real_count, row_count, spectrum_count, transfer_count, send_extent, recv_extent;
	double *q, *saved, *A, *B;
	hipfftDoubleComplex *P, *Q, *rows, *send, *recv;
	double *host_send, *host_recv;
	int *x_offset, *x_width, *y_offset, *global_x;
	std::vector<int> group_real_count, group_real_displacement;
	std::vector<int> group_spectrum_count, group_spectrum_displacement;
	std::vector<int> send_count, send_displacement, recv_count, recv_displacement;
	std::vector<double*> shared_recv;
	std::vector<size_t> shared_bytes, forward_shared_offset, backward_shared_offset;
	hipfftHandle x_forward, x_backward, y_transform;
	bool has_x, has_y;
	unsigned long long exchanges;
	double ready_download_seconds, mpi_seconds, upload_seconds, sync_seconds;
};

static void fail(MPI_Comm comm, const char *message) {
	int rank = -1;
	MPI_Comm_rank(comm, &rank);
	std::fprintf(stderr, "PFC_DCU_ERROR rank=%d %s\n", rank, message);
	std::fflush(stderr);
	MPI_Abort(comm, 1);
	std::abort();
}

static void hip_check(PfcDcu *ctx, hipError_t status, const char *operation) {
	if(status != hipSuccess) {
		char message[512];
		std::snprintf(message, sizeof(message), "%s: %s", operation, hipGetErrorString(status));
		fail(ctx->world, message);
	}
}

static void fft_check(PfcDcu *ctx, hipfftResult status, const char *operation) {
	if(status != HIPFFT_SUCCESS) {
		char message[256];
		std::snprintf(message, sizeof(message), "%s: hipFFT status %d", operation, (int)status);
		fail(ctx->world, message);
	}
}

static void mpi_check(PfcDcu *ctx, int status, const char *operation) {
	if(status != MPI_SUCCESS) fail(ctx->world, operation);
}

static size_t product(PfcDcu *ctx, size_t a, size_t b) {
	if(b != 0 && a > std::numeric_limits<size_t>::max()/b) fail(ctx->world, "array size overflow");
	return a*b;
}

static size_t padded_block(PfcDcu *ctx, size_t count) {
	if(count > std::numeric_limits<size_t>::max()-3) fail(ctx->world, "transpose padding overflow");
	return (count+3)&~size_t(3);
}

static int bind_rank(PfcDcu *ctx, MPI_Comm node, const char *pci, const cpu_set_t &assigned) {
	int node_rank, node_size;
	mpi_check(ctx, MPI_Comm_rank(node, &node_rank), "affinity MPI_Comm_rank");
	mpi_check(ctx, MPI_Comm_size(node, &node_size), "affinity MPI_Comm_size");
	cpu_set_t preferred, allowed, used;
	CPU_ZERO(&preferred); CPU_ZERO(&allowed); CPU_ZERO(&used);
	char path[256], topology[8*CPU_SETSIZE];
	std::snprintf(path, sizeof(path), "/sys/bus/pci/devices/%s/local_cpulist", pci);
	for(char *p = path; *p; ++p) if(*p >= 'A' && *p <= 'F') *p += 'a'-'A';
	FILE *file = std::fopen(path, "r");
	if(!file) fail(ctx->world, path);
	if(!std::fgets(topology, sizeof(topology), file) || std::ferror(file) || std::fclose(file)) fail(ctx->world, path);
	const char *p = topology;
	while(*p && *p != '\n') {
		char *end;
		errno = 0;
		long first = std::strtol(p, &end, 10), last = first;
		if(errno || end == p || first < 0 || first >= CPU_SETSIZE) fail(ctx->world, "invalid GPU local CPU list");
		p = end;
		if(*p == '-') {
			errno = 0;
			last = std::strtol(p+1, &end, 10);
			if(errno || end == p+1 || last < first || last >= CPU_SETSIZE) fail(ctx->world, "invalid GPU local CPU range");
			p = end;
		}
		for(long cpu = first; cpu <= last; ++cpu) CPU_SET(cpu, &preferred);
		if(*p == ',') ++p;
		else if(*p && *p != '\n') fail(ctx->world, "invalid GPU local CPU separator");
	}
	std::vector<cpu_set_t> assignments(node_size), preferences(node_size);
	mpi_check(ctx, MPI_Allgather(&assigned, sizeof(assigned), MPI_BYTE, assignments.data(), sizeof(assigned), MPI_BYTE, node), "MPI_Allgather assigned CPUs");
	mpi_check(ctx, MPI_Allgather(&preferred, sizeof(preferred), MPI_BYTE, preferences.data(), sizeof(preferred), MPI_BYTE, node), "MPI_Allgather GPU CPU preferences");
	for(int r = 0; r < node_size; ++r) for(int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
		if(CPU_ISSET(cpu, &assignments[r])) CPU_SET(cpu, &allowed);
	int selected = -1;
	for(int r = 0; r < node_size; ++r) {
		int candidate = -1;
		for(int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
			if(CPU_ISSET(cpu, &allowed) && !CPU_ISSET(cpu, &used) && CPU_ISSET(cpu, &preferences[r])) { candidate = cpu; break; }
		if(candidate < 0) for(int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
			if(CPU_ISSET(cpu, &allowed) && !CPU_ISSET(cpu, &used)) { candidate = cpu; break; }
		if(candidate < 0) fail(ctx->world, "not enough distinct assigned CPUs for scientific MPI ranks");
		CPU_SET(candidate, &used);
		if(r == node_rank) selected = candidate;
	}
	cpu_set_t singleton, actual;
	CPU_ZERO(&singleton); CPU_SET(selected, &singleton);
	CPU_ZERO(&actual);
	if(sched_setaffinity(0, sizeof(singleton), &singleton) || sched_getaffinity(0, sizeof(actual), &actual) ||
		CPU_COUNT(&actual) != 1 || !CPU_ISSET(selected, &actual)) fail(ctx->world, "could not bind rank within its node assigned CPU union");
	return selected;
}

template<class T> static void device_allocate(PfcDcu *ctx, T **pointer, size_t count) {
	hip_check(ctx, hipMalloc((void**)pointer, product(ctx, std::max(count, size_t(1)), sizeof(T))), "hipMalloc");
}

static unsigned int blocks(size_t count) {
	return (unsigned int)std::min(size_t(65535), (count+255)/256);
}

static void configure_shared(PfcDcu *ctx, const std::vector<int> &heights, const std::vector<int> &widths) {
	long system_page = sysconf(_SC_PAGESIZE);
	if(system_page < 64 || system_page%64) fail(ctx->world, "unsupported host registration page size");
	size_t page = (size_t)system_page;
	ctx->shared_recv.resize(ctx->ranks); ctx->shared_bytes.resize(ctx->ranks);
	ctx->forward_shared_offset.resize(ctx->ranks); ctx->backward_shared_offset.resize(ctx->ranks);
	for(int r = 0; r < ctx->ranks; ++r) {
		size_t outgoing = 0, incoming = 0;
		for(int k = 0; k < ctx->ranks; ++k) {
			size_t send = padded_block(ctx, product(ctx, heights[r], widths[k]));
			size_t recv = padded_block(ctx, product(ctx, widths[r], heights[k]));
			if(send > (size_t)INT_MAX/2 || recv > (size_t)INT_MAX/2 ||
				outgoing > (size_t)INT_MAX/2-send || incoming > (size_t)INT_MAX/2-recv)
				fail(ctx->world, "shared transpose extent overflow");
			if(k == ctx->rank) {
				ctx->forward_shared_offset[r] = incoming;
				ctx->backward_shared_offset[r] = outgoing;
			}
			outgoing += send; incoming += recv;
		}
		size_t capacity = std::max(outgoing, incoming);
		if(ctx->forward_shared_offset[r] > capacity || (size_t)ctx->send_count[r]/2 > capacity-ctx->forward_shared_offset[r] ||
			ctx->backward_shared_offset[r] > capacity || (size_t)ctx->recv_count[r]/2 > capacity-ctx->backward_shared_offset[r])
			fail(ctx->world, "shared transpose exceeds peer receive allocation");
		size_t bytes = product(ctx, std::max(capacity, size_t(1)), sizeof(hipfftDoubleComplex));
		if(bytes > std::numeric_limits<size_t>::max()-(page-1)) fail(ctx->world, "shared registration rounding overflow");
		ctx->shared_bytes[r] = (bytes+page-1)/page*page;
	}
	size_t bytes = ctx->shared_bytes[ctx->rank];
	if(bytes > std::numeric_limits<size_t>::max()-(page-1)) fail(ctx->world, "shared allocation size overflow");
	MPI_Aint allocation = (MPI_Aint)(bytes+page-1);
	if(allocation < 0 || (size_t)allocation != bytes+page-1) fail(ctx->world, "shared MPI allocation size overflow");
	MPI_Info info;
	mpi_check(ctx, MPI_Info_create(&info), "MPI_Info_create shared buffer");
	mpi_check(ctx, MPI_Info_set(info, "alloc_shared_noncontig", "true"), "MPI_Info_set shared buffer");
	void *owned;
	mpi_check(ctx, MPI_Win_allocate_shared(allocation, 1, info, ctx->comm, &owned, &ctx->shared_window), "MPI_Win_allocate_shared");
	mpi_check(ctx, MPI_Info_free(&info), "MPI_Info_free shared buffer");
	mpi_check(ctx, MPI_Win_set_errhandler(ctx->shared_window, MPI_ERRORS_RETURN), "MPI_Win_set_errhandler");
	for(int r = 0; r < ctx->ranks; ++r) {
		MPI_Aint size;
		int unit;
		void *base;
		mpi_check(ctx, MPI_Win_shared_query(ctx->shared_window, r, &size, &unit, &base), "MPI_Win_shared_query");
		if(size < 0 || unit != 1 || !base) fail(ctx->world, "invalid shared buffer mapping");
		std::uintptr_t address = reinterpret_cast<std::uintptr_t>(base);
		if(address > std::numeric_limits<std::uintptr_t>::max()-(page-1)) fail(ctx->world, "shared pointer alignment overflow");
		std::uintptr_t aligned = (address+page-1)/page*page;
		if((size_t)size < aligned-address || ctx->shared_bytes[r] > (size_t)size-(aligned-address)) fail(ctx->world, "shared buffer registration extent overflow");
		ctx->shared_recv[r] = reinterpret_cast<double*>(aligned);
	}
	ctx->host_recv = ctx->shared_recv[ctx->rank];
	std::memset(ctx->host_recv, 0, bytes);
	mpi_check(ctx, MPI_Win_lock_all(0, ctx->shared_window), "MPI_Win_lock_all shared buffer");
	mpi_check(ctx, MPI_Win_sync(ctx->shared_window), "MPI_Win_sync first touch");
	mpi_check(ctx, MPI_Barrier(ctx->comm), "shared first touch MPI_Barrier");
	mpi_check(ctx, MPI_Win_sync(ctx->shared_window), "MPI_Win_sync first touch acquire");
	for(int r = 0; r < ctx->ranks; ++r)
		hip_check(ctx, hipHostRegister(ctx->shared_recv[r], ctx->shared_bytes[r], hipHostRegisterPortable), "hipHostRegister shared receive buffer");
	mpi_check(ctx, MPI_Barrier(ctx->comm), "shared registration ready MPI_Barrier");
}

__global__ static void pack_rows(const hipfftDoubleComplex *rows, hipfftDoubleComplex *packed,
	size_t count, int Wc, const int *offset, const int *width) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		int x = i%Wc;
		size_t y = i/Wc;
		packed[offset[x]+y*width[x]] = rows[i];
	}
}

__global__ static void unpack_spectrum(const hipfftDoubleComplex *packed, hipfftDoubleComplex *spectrum,
	size_t count, int H, const int *offset) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		int y = i%H;
		size_t x = i/H;
		spectrum[i] = packed[offset[y]+x];
	}
}

__global__ static void pack_spectrum(const hipfftDoubleComplex *spectrum, hipfftDoubleComplex *packed,
	size_t count, int H, const int *offset) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		int y = i%H;
		size_t x = i/H;
		packed[offset[y]+x] = spectrum[i];
	}
}

__global__ static void unpack_rows(const hipfftDoubleComplex *packed, hipfftDoubleComplex *rows,
	size_t count, int Wc, const int *offset, const int *width) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		int x = i%Wc;
		size_t y = i/Wc;
		rows[i] = packed[offset[x]+y*width[x]];
	}
}

__global__ static void scale_spectrum(hipfftDoubleComplex *data, size_t count, double scale) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		data[i].x *= scale;
		data[i].y *= scale;
	}
}

__global__ static void nonlinear(double *q, size_t count, int W, int Wp, double gamma, double delta) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		if(i%Wp < (size_t)W) {
			double p = q[i];
			double p2 = p*p;
			q[i] = gamma*p2+delta*p*p2;
		}
	}
}

__global__ static void advance(hipfftDoubleComplex *P, hipfftDoubleComplex *Q,
	const double *A, const double *B, size_t count) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		P[i].x = A[i]*P[i].x+B[i]*Q[i].x;
		P[i].y = A[i]*P[i].y+B[i]*Q[i].y;
		Q[i] = P[i];
	}
}

__global__ static void diagnostic_operator(hipfftDoubleComplex *Q, size_t count,
	int H, const int *global_x, double divWH, double dkx, double dky) {
	for(size_t i = blockIdx.x*blockDim.x+threadIdx.x; i < count; i += (size_t)blockDim.x*gridDim.x) {
		int h = i%H;
		int gw = global_x[i/H];
		double kx2 = gw*dkx;
		kx2 *= kx2;
		double ky = (h < H/2 ? h : h-H)*dky;
		double k2 = kx2+ky*ky;
		double d2 = 1.0-k2;
		d2 *= d2;
		Q[i].x *= divWH;
		Q[i].y *= divWH;
		Q[i].x *= d2;
		Q[i].y *= d2;
	}
}

static void launch_check(PfcDcu *ctx) {
	hip_check(ctx, hipGetLastError(), "HIP kernel launch");
}

static void exchange(PfcDcu *ctx, bool inverse) {
	const std::vector<int> &counts = inverse ? ctx->recv_count : ctx->send_count;
	const std::vector<int> &displacements = inverse ? ctx->recv_displacement : ctx->send_displacement;
	const std::vector<int> &incoming_displacements = inverse ? ctx->send_displacement : ctx->recv_displacement;
	size_t outgoing = inverse ? ctx->recv_extent : ctx->send_extent;
	size_t incoming = inverse ? ctx->send_extent : ctx->recv_extent;
	size_t self_begin = (size_t)incoming_displacements[ctx->rank]/2;
	size_t self_extent = padded_block(ctx, (size_t)counts[ctx->rank]/2);
	if(self_begin > incoming || self_extent > incoming-self_begin) fail(ctx->world, "self transpose receive extent overflow");
	size_t self_end = self_begin+self_extent;
	// Wall times include preceding GPU work and transport peer waits, not isolated bandwidth.
	double start = MPI_Wtime();
	if(ctx->shared_transport) {
		const std::vector<size_t> &offsets = inverse ? ctx->backward_shared_offset : ctx->forward_shared_offset;
		for(int r = 0; r < ctx->ranks; ++r) if(counts[r]) {
			if(r == ctx->rank)
				hip_check(ctx, hipMemcpyAsync(ctx->recv+self_begin, ctx->send+displacements[r]/2,
					(size_t)counts[r]*sizeof(double), hipMemcpyDeviceToDevice, 0), "copy self transpose on device");
			else hip_check(ctx, hipMemcpyAsync(ctx->shared_recv[r]+2*offsets[r], ctx->send+displacements[r]/2,
				(size_t)counts[r]*sizeof(double), hipMemcpyDeviceToHost, 0), "push shared host transpose");
		}
		hip_check(ctx, hipDeviceSynchronize(), "complete shared host pushes");
	} else hip_check(ctx, hipMemcpy(ctx->host_send, ctx->send, outgoing*sizeof(hipfftDoubleComplex), hipMemcpyDeviceToHost), "transpose staging download");
	ctx->ready_download_seconds += MPI_Wtime()-start;
	start = MPI_Wtime();
	// The MPI bucket covers shared-memory visibility/reuse waits or Alltoallv.
	if(ctx->shared_transport) {
		// Publish DMA writes before receivers read; the epoch also covers direct shared access.
		mpi_check(ctx, MPI_Win_sync(ctx->shared_window), "MPI_Win_sync transpose publish");
		mpi_check(ctx, MPI_Barrier(ctx->comm), "shared transpose publish MPI_Barrier");
		mpi_check(ctx, MPI_Win_sync(ctx->shared_window), "MPI_Win_sync transpose acquire");
	} else mpi_check(ctx, MPI_Alltoallv(ctx->host_send, counts.data(), displacements.data(), MPI_DOUBLE,
		ctx->host_recv, inverse ? ctx->send_count.data() : ctx->recv_count.data(),
		inverse ? ctx->send_displacement.data() : ctx->recv_displacement.data(), MPI_DOUBLE, ctx->comm), "transpose MPI_Alltoallv");
	ctx->mpi_seconds += MPI_Wtime()-start;
	start = MPI_Wtime();
	if(ctx->shared_transport) {
		// The self block is already on the device; upload the two surrounding spans only.
		if(self_begin) hip_check(ctx, hipMemcpyAsync(ctx->recv, ctx->host_recv,
			self_begin*sizeof(hipfftDoubleComplex), hipMemcpyHostToDevice, 0), "transpose prefix upload");
		if(self_end < incoming) hip_check(ctx, hipMemcpyAsync(ctx->recv+self_end, ctx->host_recv+2*self_end,
			(incoming-self_end)*sizeof(hipfftDoubleComplex), hipMemcpyHostToDevice, 0), "transpose suffix upload");
		hip_check(ctx, hipDeviceSynchronize(), "complete transpose surrounding uploads");
	} else hip_check(ctx, hipMemcpy(ctx->recv, ctx->host_recv, incoming*sizeof(hipfftDoubleComplex), hipMemcpyHostToDevice), "transpose staging upload");
	ctx->upload_seconds += MPI_Wtime()-start;
	if(ctx->shared_transport) {
		start = MPI_Wtime();
		mpi_check(ctx, MPI_Barrier(ctx->comm), "shared transpose receive reuse MPI_Barrier");
		ctx->mpi_seconds += MPI_Wtime()-start;
	}
	ctx->exchanges++;
}

static void forward(PfcDcu *ctx, double *real, hipfftDoubleComplex *spectrum) {
	if(ctx->has_x) fft_check(ctx, hipfftExecD2Z(ctx->x_forward, real, ctx->rows), "x forward");
	if(ctx->row_count) {
		pack_rows<<<blocks(ctx->row_count), 256>>>(ctx->rows, ctx->send, ctx->row_count,
			ctx->Wc, ctx->x_offset, ctx->x_width);
		launch_check(ctx);
	}
	exchange(ctx, false);
	if(ctx->spectrum_count) {
		unpack_spectrum<<<blocks(ctx->spectrum_count), 256>>>(ctx->recv, spectrum, ctx->spectrum_count, ctx->H, ctx->y_offset);
		launch_check(ctx);
	}
	if(ctx->has_y) fft_check(ctx, hipfftExecZ2Z(ctx->y_transform, spectrum, spectrum, HIPFFT_FORWARD), "y forward");
}

static void backward(PfcDcu *ctx, hipfftDoubleComplex *spectrum, double *real) {
	if(ctx->has_y) fft_check(ctx, hipfftExecZ2Z(ctx->y_transform, spectrum, spectrum, HIPFFT_BACKWARD), "y backward");
	if(ctx->spectrum_count) {
		pack_spectrum<<<blocks(ctx->spectrum_count), 256>>>(spectrum, ctx->send, ctx->spectrum_count, ctx->H, ctx->y_offset);
		launch_check(ctx);
	}
	exchange(ctx, true);
	if(ctx->row_count) {
		unpack_rows<<<blocks(ctx->row_count), 256>>>(ctx->recv, ctx->rows, ctx->row_count,
			ctx->Wc, ctx->x_offset, ctx->x_width);
		launch_check(ctx);
	}
	if(ctx->has_x) fft_check(ctx, hipfftExecZ2D(ctx->x_backward, ctx->rows, real), "x backward");
}

extern "C" PfcDcu *pfc_dcu_create(int W, int H, ptrdiff_t lH, ptrdiff_t lh0,
	ptrdiff_t lW, ptrdiff_t lw0, MPI_Comm comm) {
	PfcDcu *ctx = new(std::nothrow) PfcDcu{};
	if(!ctx) fail(comm, "context allocation failed");
	ctx->world = comm;
	try {
		mpi_check(ctx, MPI_Comm_dup(comm, &ctx->world), "MPI_Comm_dup world");
		mpi_check(ctx, MPI_Comm_set_errhandler(ctx->world, MPI_ERRORS_RETURN), "MPI_Comm_set_errhandler");
		mpi_check(ctx, MPI_Comm_rank(ctx->world, &ctx->world_rank), "MPI_Comm_rank world");
		cpu_set_t original_affinity;
		CPU_ZERO(&original_affinity);
		if(sched_getaffinity(0, sizeof(original_affinity), &original_affinity)) fail(ctx->world, "sched_getaffinity failed");
		if(W <= 0 || W > INT_MAX-2 || H <= 0 || lH < 0 || lH > H || lW < 0 || lW > W/2+1 ||
			lh0 < 0 || lh0 > H-lH || lw0 < 0 || lw0 > W/2+1-lW) fail(ctx->world, "invalid FFTW slab dimensions");
		ctx->W = W; ctx->H = H; ctx->Wc = W/2+1; ctx->Wp = 2*ctx->Wc;
		size_t local_real = product(ctx, lH, ctx->Wp), local_spectrum = product(ctx, lW, H);
		if(local_real > INT_MAX || local_spectrum > INT_MAX) fail(ctx->world, "group transfer exceeds MPI int count range");
		ctx->local_real_count = (int)local_real;
		ctx->local_spectrum_count = (int)local_spectrum;

		MPI_Comm shared;
		mpi_check(ctx, MPI_Comm_split_type(ctx->world, MPI_COMM_TYPE_SHARED, ctx->world_rank, MPI_INFO_NULL, &shared), "MPI_Comm_split_type");
		int local_rank = 0;
		mpi_check(ctx, MPI_Comm_rank(shared, &local_rank), "shared MPI_Comm_rank");
		int card = local_rank%4;
		struct DeviceInfo { int visible, device, leader; char pci[64], name[256]; } node_devices[4]{}, info{};
		// Initialize the other leaders' HIP helper threads only after their CPU binding.
		if(local_rank == 0) {
			int visible;
			hip_check(ctx, hipGetDeviceCount(&visible), "node hipGetDeviceCount");
			if(visible < 4) fail(ctx->world, "node rank 0 must see all four allocated DCUs");
			hip_check(ctx, hipSetDevice(0), "node hipSetDevice");
			for(int d = 0; d < 4; ++d) {
				node_devices[d].visible = visible; node_devices[d].device = d;
				hipDeviceProp_t properties;
				hip_check(ctx, hipGetDeviceProperties(&properties, d), "node hipGetDeviceProperties");
				hip_check(ctx, hipDeviceGetPCIBusId(node_devices[d].pci, sizeof(node_devices[d].pci), d), "node hipDeviceGetPCIBusId");
				std::snprintf(node_devices[d].name, sizeof(node_devices[d].name), "%s", properties.name);
				for(char *p = node_devices[d].name; *p; ++p) if(*p == ' ' || *p == '\t' || *p == '\n') *p = '_';
			}
		}
		mpi_check(ctx, MPI_Bcast(node_devices, sizeof(node_devices), MPI_BYTE, 0, shared), "MPI_Bcast node device identities");
		int cpu = bind_rank(ctx, shared, node_devices[card].pci, original_affinity);
		mpi_check(ctx, MPI_Comm_split(shared, card, ctx->world_rank, &ctx->group), "MPI_Comm_split card group");
		mpi_check(ctx, MPI_Comm_free(&shared), "MPI_Comm_free shared");
		mpi_check(ctx, MPI_Comm_rank(ctx->group, &ctx->group_rank), "card MPI_Comm_rank");
		mpi_check(ctx, MPI_Comm_size(ctx->group, &ctx->group_size), "card MPI_Comm_size");
		ctx->leader = ctx->group_rank == 0;
		mpi_check(ctx, MPI_Comm_split(ctx->world, ctx->leader ? 0 : MPI_UNDEFINED,
			ctx->world_rank, &ctx->comm), "MPI_Comm_split leaders");

		std::vector<long long> members(ctx->leader ? (size_t)4*ctx->group_size : 0);
		long long original[] = {(long long)lH, (long long)lh0, (long long)lW, (long long)lw0};
		mpi_check(ctx, MPI_Gather(original, 4, MPI_LONG_LONG, members.data(), 4, MPI_LONG_LONG, 0, ctx->group), "MPI_Gather card slabs");
		std::vector<int> local_y, local_x;
		if(ctx->leader) {
			mpi_check(ctx, MPI_Comm_rank(ctx->comm, &ctx->rank), "leaders MPI_Comm_rank");
			mpi_check(ctx, MPI_Comm_size(ctx->comm, &ctx->ranks), "leaders MPI_Comm_size");
			ctx->group_real_count.resize(ctx->group_size); ctx->group_real_displacement.resize(ctx->group_size);
			ctx->group_spectrum_count.resize(ctx->group_size); ctx->group_spectrum_displacement.resize(ctx->group_size);
			size_t real_offset = 0, spectrum_offset = 0;
			for(int m = 0; m < ctx->group_size; ++m) {
				long long mH = members[4*m], mh0 = members[4*m+1], mW = members[4*m+2], mw0 = members[4*m+3];
				if(mH < 0 || mH > H || mh0 < 0 || mh0 > H-mH || mW < 0 || mW > ctx->Wc || mw0 < 0 || mw0 > ctx->Wc-mW)
					fail(ctx->world, "invalid card member FFTW slab");
				size_t real_count = product(ctx, mH, ctx->Wp), spectrum_count = product(ctx, mW, H);
				if(real_count > INT_MAX || spectrum_count > INT_MAX || real_offset+real_count > INT_MAX || spectrum_offset+spectrum_count > INT_MAX)
					fail(ctx->world, "card aggregate exceeds MPI int count range");
				ctx->group_real_count[m] = (int)real_count;
				ctx->group_real_displacement[m] = (int)real_offset;
				ctx->group_spectrum_count[m] = (int)spectrum_count;
				ctx->group_spectrum_displacement[m] = (int)spectrum_offset;
				real_offset += real_count; spectrum_offset += spectrum_count;
				for(long long y = 0; y < mH; ++y) local_y.push_back((int)(mh0+y));
				for(long long x = 0; x < mW; ++x) local_x.push_back((int)(mw0+x));
			}
			if(local_y.size() > (size_t)H || local_x.size() > (size_t)ctx->Wc) fail(ctx->world, "invalid card aggregate dimensions");
			ctx->lH = local_y.size(); ctx->lW = local_x.size();
			ctx->real_count = real_offset;
			ctx->row_count = product(ctx, ctx->lH, ctx->Wc);
			ctx->spectrum_count = spectrum_offset;
			ctx->transfer_count = std::max(ctx->row_count, ctx->spectrum_count);
			if(ctx->transfer_count > (size_t)INT_MAX/2) fail(ctx->world, "MPI transpose exceeds int count/displacement range");
			hip_check(ctx, hipGetDeviceCount(&info.visible), "hipGetDeviceCount");
			if(info.visible < 1) fail(ctx->world, "no scheduler-visible DCU");
			info.device = card%info.visible;
			info.leader = ctx->world_rank;
			hip_check(ctx, hipSetDevice(info.device), "hipSetDevice");
			hip_check(ctx, hipDeviceGetPCIBusId(info.pci, sizeof(info.pci), info.device), "hipDeviceGetPCIBusId");
			if(std::strcmp(info.pci, node_devices[card].pci)) fail(ctx->world, "leader device visibility disagrees with node PCI assignment");
			std::snprintf(info.name, sizeof(info.name), "%s", node_devices[card].name);
			MPI_Comm node_leaders;
			mpi_check(ctx, MPI_Comm_split_type(ctx->comm, MPI_COMM_TYPE_SHARED, ctx->rank, MPI_INFO_NULL, &node_leaders), "MPI_Comm_split_type node leaders");
			int node_count;
			mpi_check(ctx, MPI_Comm_size(node_leaders, &node_count), "node leaders MPI_Comm_size");
			ctx->shared_transport = node_count == ctx->ranks;
			std::vector<char> identities(sizeof(info.pci)*node_count);
			mpi_check(ctx, MPI_Allgather(info.pci, sizeof(info.pci), MPI_CHAR, identities.data(), sizeof(info.pci), MPI_CHAR, node_leaders), "MPI_Allgather physical devices");
			for(int a = 0; a < node_count; ++a) for(int b = 0; b < a; ++b)
				if(!std::strcmp(identities.data()+a*sizeof(info.pci), identities.data()+b*sizeof(info.pci))) fail(ctx->world, "card leaders mapped to the same physical device");
			mpi_check(ctx, MPI_Comm_free(&node_leaders), "MPI_Comm_free node leaders");
		}
		mpi_check(ctx, MPI_Bcast(&info, sizeof(info), MPI_BYTE, 0, ctx->group), "MPI_Bcast card identity");
		char host[MPI_MAX_PROCESSOR_NAME+1];
		int host_length;
		mpi_check(ctx, MPI_Get_processor_name(host, &host_length), "MPI_Get_processor_name");
		host[host_length] = '\0';
		std::printf("PFC_DEVICE rank=%d local_rank=%d host=%s device=%d pci=%s visible=%d name=%s role=%s leader_rank=%d cpu=%d\n",
			ctx->world_rank, local_rank, host, info.device, info.pci, info.visible, info.name, ctx->leader ? "leader" : "member", info.leader, cpu);
		std::printf("PFC_SLAB rank=%d W=%d H=%d lH=%td lh0=%td lW=%td lw0=%td alloc=%zu\n",
			ctx->world_rank, W, H, lH, lh0, lW, lw0, std::max(product(ctx, lH, ctx->Wc), local_spectrum));
		std::fflush(stdout);
		if(!ctx->leader) return ctx;

		// Card batches concatenate original slabs; explicit indices retain their global coordinates.
		std::vector<int> dimensions((size_t)2*ctx->ranks), heights(ctx->ranks), widths(ctx->ranks), h_start(ctx->ranks), w_start(ctx->ranks);
		int aggregate[] = {(int)ctx->lH, (int)ctx->lW};
		mpi_check(ctx, MPI_Allgather(aggregate, 2, MPI_INT, dimensions.data(), 2, MPI_INT, ctx->comm), "MPI_Allgather leader dimensions");
		int total_h = 0, total_w = 0;
		for(int r = 0; r < ctx->ranks; ++r) {
			heights[r] = dimensions[2*r]; widths[r] = dimensions[2*r+1];
			if(heights[r] < 0 || heights[r] > H-total_h || widths[r] < 0 || widths[r] > ctx->Wc-total_w) fail(ctx->world, "invalid leader dimensions");
			h_start[r] = total_h; w_start[r] = total_w;
			total_h += heights[r]; total_w += widths[r];
		}
		if(total_h != H || total_w != ctx->Wc) fail(ctx->world, "incomplete leader dimensions");
		std::vector<int> all_y(H), all_x(ctx->Wc);
		mpi_check(ctx, MPI_Allgatherv(local_y.data(), aggregate[0], MPI_INT, all_y.data(), heights.data(), h_start.data(), MPI_INT, ctx->comm), "MPI_Allgatherv global y indices");
		mpi_check(ctx, MPI_Allgatherv(local_x.data(), aggregate[1], MPI_INT, all_x.data(), widths.data(), w_start.data(), MPI_INT, ctx->comm), "MPI_Allgatherv global x indices");
		ctx->send_count.resize(ctx->ranks); ctx->send_displacement.resize(ctx->ranks);
		ctx->recv_count.resize(ctx->ranks); ctx->recv_displacement.resize(ctx->ranks);
		std::vector<int> x_offset(ctx->Wc, -1), x_width(ctx->Wc, -1), y_offset(H, -1);
		size_t send_offset = 0, recv_offset = 0, send_samples = 0, recv_samples = 0;
		for(int r = 0; r < ctx->ranks; ++r) {
			size_t send_count = product(ctx, ctx->lH, widths[r]), recv_count = product(ctx, ctx->lW, heights[r]);
			// Isolate GPU writers on separate host cache lines, including non-divisible slabs.
			size_t send_padded = padded_block(ctx, send_count), recv_padded = padded_block(ctx, recv_count);
			if(send_padded > (size_t)INT_MAX/2 || recv_padded > (size_t)INT_MAX/2 ||
				send_offset > (size_t)INT_MAX/2-send_padded || recv_offset > (size_t)INT_MAX/2-recv_padded)
				fail(ctx->world, "MPI transpose count/displacement overflow");
			ctx->send_count[r] = (int)(2*send_count); ctx->send_displacement[r] = (int)(2*send_offset);
			ctx->recv_count[r] = (int)(2*recv_count); ctx->recv_displacement[r] = (int)(2*recv_offset);
			for(int x = 0; x < widths[r]; ++x) {
				int gx = all_x[w_start[r]+x];
				if(gx < 0 || gx >= ctx->Wc || x_offset[gx] != -1) fail(ctx->world, "overlapping FFTW frequency slabs");
				x_offset[gx] = (int)(send_offset+x); x_width[gx] = widths[r];
			}
			for(int y = 0; y < heights[r]; ++y) {
				int gy = all_y[h_start[r]+y];
				if(gy < 0 || gy >= H || y_offset[gy] != -1) fail(ctx->world, "overlapping FFTW real slabs");
				y_offset[gy] = (int)(recv_offset+(size_t)y*ctx->lW);
			}
			send_offset += send_padded; recv_offset += recv_padded;
			send_samples += send_count; recv_samples += recv_count;
		}
		if(send_samples != ctx->row_count || recv_samples != ctx->spectrum_count ||
			std::find(x_offset.begin(), x_offset.end(), -1) != x_offset.end() ||
			std::find(y_offset.begin(), y_offset.end(), -1) != y_offset.end()) fail(ctx->world, "incomplete FFTW slab coverage");
		ctx->send_extent = send_offset; ctx->recv_extent = recv_offset;
		ctx->transfer_count = std::max(send_offset, recv_offset);

		device_allocate(ctx, &ctx->q, ctx->real_count);
		device_allocate(ctx, &ctx->saved, ctx->real_count);
		device_allocate(ctx, &ctx->A, ctx->spectrum_count);
		device_allocate(ctx, &ctx->B, ctx->spectrum_count);
		device_allocate(ctx, &ctx->P, ctx->spectrum_count);
		device_allocate(ctx, &ctx->Q, ctx->spectrum_count);
		device_allocate(ctx, &ctx->rows, ctx->row_count);
		device_allocate(ctx, &ctx->send, ctx->transfer_count);
		device_allocate(ctx, &ctx->recv, ctx->transfer_count);
		device_allocate(ctx, &ctx->x_offset, ctx->Wc);
		device_allocate(ctx, &ctx->x_width, ctx->Wc);
		device_allocate(ctx, &ctx->y_offset, H);
		device_allocate(ctx, &ctx->global_x, ctx->lW);
		size_t staging_bytes = product(ctx, std::max(ctx->transfer_count, size_t(1)), sizeof(hipfftDoubleComplex));
		hip_check(ctx, hipHostMalloc((void**)&ctx->host_send, staging_bytes), "hipHostMalloc send");
		if(ctx->shared_transport) configure_shared(ctx, heights, widths);
		else hip_check(ctx, hipHostMalloc((void**)&ctx->host_recv, staging_bytes), "hipHostMalloc recv");
		hip_check(ctx, hipMemcpy(ctx->x_offset, x_offset.data(), ctx->Wc*sizeof(int), hipMemcpyHostToDevice), "x offset upload");
		hip_check(ctx, hipMemcpy(ctx->x_width, x_width.data(), ctx->Wc*sizeof(int), hipMemcpyHostToDevice), "x width upload");
		hip_check(ctx, hipMemcpy(ctx->y_offset, y_offset.data(), H*sizeof(int), hipMemcpyHostToDevice), "y offset upload");
		if(ctx->lW) hip_check(ctx, hipMemcpy(ctx->global_x, local_x.data(), ctx->lW*sizeof(int), hipMemcpyHostToDevice), "global x index upload");
		hip_check(ctx, hipMemset(ctx->q, 0, std::max(ctx->real_count, size_t(1))*sizeof(double)), "initialize real padding");
		if(ctx->lH) {
			int n[] = {W}, real_embed[] = {ctx->Wp}, complex_embed[] = {ctx->Wc};
			fft_check(ctx, hipfftPlanMany(&ctx->x_forward, 1, n, real_embed, 1, ctx->Wp,
				complex_embed, 1, ctx->Wc, HIPFFT_D2Z, (int)ctx->lH), "plan x forward");
			fft_check(ctx, hipfftPlanMany(&ctx->x_backward, 1, n, complex_embed, 1, ctx->Wc,
				real_embed, 1, ctx->Wp, HIPFFT_Z2D, (int)ctx->lH), "plan x backward");
			ctx->has_x = true;
		}
		if(ctx->lW) {
			int n[] = {H}, embed[] = {H};
			fft_check(ctx, hipfftPlanMany(&ctx->y_transform, 1, n, embed, 1, H,
				embed, 1, H, HIPFFT_Z2Z, (int)ctx->lW), "plan y complex");
			ctx->has_y = true;
		}
		std::printf("PFC_EXEC leader_rank=%d members=%d leaders=%d lH=%td lW=%td\n", ctx->world_rank, ctx->group_size, ctx->ranks, ctx->lH, ctx->lW);
		std::printf("PFC_TRANSFER leader_rank=%d mode=%s\n", ctx->world_rank, ctx->shared_transport ? "host-shared" : "host-mpi");
		std::fflush(stdout);
		pfc_dcu_sync(ctx);
		return ctx;
	}
	catch(const std::bad_alloc&) { fail(ctx->world, "host allocation failed"); }
	return NULL;
}

static void scatter_real(PfcDcu *ctx, double *host) {
	mpi_check(ctx, MPI_Scatterv(ctx->host_send, ctx->group_real_count.data(), ctx->group_real_displacement.data(), MPI_DOUBLE,
		host, ctx->local_real_count, MPI_DOUBLE, 0, ctx->group), "MPI_Scatterv original real slabs");
}

extern "C" void pfc_dcu_upload(PfcDcu *ctx, const double *host) {
	mpi_check(ctx, MPI_Gatherv(host, ctx->local_real_count, MPI_DOUBLE, ctx->host_send,
		ctx->group_real_count.data(), ctx->group_real_displacement.data(), MPI_DOUBLE, 0, ctx->group), "MPI_Gatherv initial density");
	if(ctx->leader) hip_check(ctx, hipMemcpy(ctx->q, ctx->host_send, ctx->real_count*sizeof(double), hipMemcpyHostToDevice), "initial density upload");
}

extern "C" void pfc_dcu_download(PfcDcu *ctx, double *host) {
	if(ctx->leader) hip_check(ctx, hipMemcpy(ctx->host_send, ctx->q, ctx->real_count*sizeof(double), hipMemcpyDeviceToHost), "density download");
	scatter_real(ctx, host);
}

extern "C" void pfc_dcu_upload_ab(PfcDcu *ctx, const double *A, const double *B) {
	mpi_check(ctx, MPI_Gatherv(A, ctx->local_spectrum_count, MPI_DOUBLE, ctx->host_send,
		ctx->group_spectrum_count.data(), ctx->group_spectrum_displacement.data(), MPI_DOUBLE, 0, ctx->group), "MPI_Gatherv A");
	if(ctx->leader) hip_check(ctx, hipMemcpy(ctx->A, ctx->host_send, ctx->spectrum_count*sizeof(double), hipMemcpyHostToDevice), "A upload");
	mpi_check(ctx, MPI_Gatherv(B, ctx->local_spectrum_count, MPI_DOUBLE, ctx->host_send,
		ctx->group_spectrum_count.data(), ctx->group_spectrum_displacement.data(), MPI_DOUBLE, 0, ctx->group), "MPI_Gatherv B");
	if(ctx->leader) hip_check(ctx, hipMemcpy(ctx->B, ctx->host_send, ctx->spectrum_count*sizeof(double), hipMemcpyHostToDevice), "B upload");
}

extern "C" int pfc_dcu_is_leader(const PfcDcu *ctx) {
	return ctx->leader ? 1 : 0;
}

extern "C" void pfc_dcu_step(PfcDcu *ctx, int t, double gamma, double delta) {
	if(!ctx->leader) return;
	if(t%100 == 0) {
		forward(ctx, ctx->q, ctx->P);
		if(ctx->spectrum_count) {
			scale_spectrum<<<blocks(ctx->spectrum_count), 256>>>(ctx->P, ctx->spectrum_count, 1.0/ctx->W/ctx->H);
			launch_check(ctx);
		}
	}
	if(ctx->real_count) {
		nonlinear<<<blocks(ctx->real_count), 256>>>(ctx->q, ctx->real_count, ctx->W, ctx->Wp, gamma, delta);
		launch_check(ctx);
	}
	forward(ctx, ctx->q, ctx->Q);
	if(ctx->spectrum_count) {
		advance<<<blocks(ctx->spectrum_count), 256>>>(ctx->P, ctx->Q, ctx->A, ctx->B, ctx->spectrum_count);
		launch_check(ctx);
	}
	backward(ctx, ctx->Q, ctx->q);
}

extern "C" void pfc_dcu_diagnostic(PfcDcu *ctx, double *host_density,
	double *host_operator, double dx, double dy) {
	if(ctx->leader) {
		size_t bytes = ctx->real_count*sizeof(double);
		hip_check(ctx, hipMemcpy(ctx->saved, ctx->q, bytes, hipMemcpyDeviceToDevice), "save diagnostic density");
		hip_check(ctx, hipMemcpy(ctx->host_send, ctx->q, bytes, hipMemcpyDeviceToHost), "diagnostic density download");
	}
	scatter_real(ctx, host_density);
	if(ctx->leader) {
		forward(ctx, ctx->q, ctx->Q);
		if(ctx->spectrum_count) {
			double divWH = 1.0/ctx->W/ctx->H;
			double dkx = 2.0*3.14159265358979323846/(dx*ctx->W);
			double dky = 2.0*3.14159265358979323846/(dy*ctx->H);
			diagnostic_operator<<<blocks(ctx->spectrum_count), 256>>>(ctx->Q, ctx->spectrum_count,
				ctx->H, ctx->global_x, divWH, dkx, dky);
			launch_check(ctx);
		}
		backward(ctx, ctx->Q, ctx->q);
		hip_check(ctx, hipMemcpy(ctx->host_send, ctx->q, ctx->real_count*sizeof(double), hipMemcpyDeviceToHost), "diagnostic operator download");
	}
	scatter_real(ctx, host_operator);
	if(ctx->leader) {
		hip_check(ctx, hipMemcpy(ctx->q, ctx->saved, ctx->real_count*sizeof(double), hipMemcpyDeviceToDevice), "restore diagnostic density");
		forward(ctx, ctx->q, ctx->P);
		if(ctx->spectrum_count) {
			scale_spectrum<<<blocks(ctx->spectrum_count), 256>>>(ctx->P, ctx->spectrum_count, 1.0/ctx->W/ctx->H);
			launch_check(ctx);
		}
		pfc_dcu_sync(ctx);
	}
}

extern "C" void pfc_dcu_sync(PfcDcu *ctx) {
	if(ctx->leader) {
		double start = MPI_Wtime();
		hip_check(ctx, hipDeviceSynchronize(), "hipDeviceSynchronize");
		ctx->sync_seconds += MPI_Wtime()-start;
	}
}

extern "C" void pfc_dcu_destroy(PfcDcu *ctx) {
	if(!ctx) return;
	if(ctx->leader) {
		pfc_dcu_sync(ctx);
		if(ctx->shared_transport) {
			mpi_check(ctx, MPI_Barrier(ctx->comm), "shared DMA completion MPI_Barrier");
			for(int r = 0; r < ctx->ranks; ++r)
				hip_check(ctx, hipHostUnregister(ctx->shared_recv[r]), "hipHostUnregister shared receive buffer");
			mpi_check(ctx, MPI_Barrier(ctx->comm), "shared unregister MPI_Barrier");
			mpi_check(ctx, MPI_Win_unlock_all(ctx->shared_window), "MPI_Win_unlock_all shared buffer");
			mpi_check(ctx, MPI_Win_free(&ctx->shared_window), "MPI_Win_free shared buffer");
		}
		if(ctx->has_x) {
			fft_check(ctx, hipfftDestroy(ctx->x_forward), "destroy x forward");
			fft_check(ctx, hipfftDestroy(ctx->x_backward), "destroy x backward");
		}
		if(ctx->has_y) fft_check(ctx, hipfftDestroy(ctx->y_transform), "destroy y complex");
		hip_check(ctx, hipFree(ctx->q), "free density");
		hip_check(ctx, hipFree(ctx->saved), "free saved density");
		hip_check(ctx, hipFree(ctx->A), "free A");
		hip_check(ctx, hipFree(ctx->B), "free B");
		hip_check(ctx, hipFree(ctx->P), "free P");
		hip_check(ctx, hipFree(ctx->Q), "free Q");
		hip_check(ctx, hipFree(ctx->rows), "free rows");
		hip_check(ctx, hipFree(ctx->send), "free send");
		hip_check(ctx, hipFree(ctx->recv), "free receive");
		hip_check(ctx, hipFree(ctx->x_offset), "free x offsets");
		hip_check(ctx, hipFree(ctx->x_width), "free x widths");
		hip_check(ctx, hipFree(ctx->y_offset), "free y offsets");
		hip_check(ctx, hipFree(ctx->global_x), "free global x indices");
		hip_check(ctx, hipHostFree(ctx->host_send), "free pinned send");
		if(!ctx->shared_transport) hip_check(ctx, hipHostFree(ctx->host_recv), "free pinned receive");
		std::printf("PFC_PROFILE leader_rank=%d exchanges=%llu ready_download_seconds=%.9f mpi_seconds=%.9f upload_seconds=%.9f sync_seconds=%.9f\n",
			ctx->world_rank, ctx->exchanges, ctx->ready_download_seconds, ctx->mpi_seconds, ctx->upload_seconds, ctx->sync_seconds);
		std::fflush(stdout);
		mpi_check(ctx, MPI_Comm_free(&ctx->comm), "MPI_Comm_free leaders");
	}
	mpi_check(ctx, MPI_Comm_free(&ctx->group), "MPI_Comm_free card group");
	if(MPI_Comm_free(&ctx->world) != MPI_SUCCESS) fail(MPI_COMM_WORLD, "MPI_Comm_free world");
	delete ctx;
}
