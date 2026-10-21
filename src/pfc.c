/*
Author: Petri Hirvonen, petri.hirvonen@aalto.fi
December 7, 2017
*/

#define _DEFAULT_SOURCE
#include <ctype.h>
#include <complex.h>
#include <dirent.h>
#include <errno.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef PFC_DCU
#include "pfc_dcu.h"
#endif

#define pi 3.14159265358979323846264338327
#define PFC_PATH 4096

static void fail(const char *message) {
	fprintf(stderr, "PFC error: %s (%s)\n", message, strerror(errno));
	MPI_Abort(MPI_COMM_WORLD, 1);
	abort();
}

#ifdef PFC_DCU
// glibc rand/srand use this same 128-byte generator; HIP must not consume its state.
static struct random_data physical_rng;
static int32_t physical_rng_state[32];
static int physical_rand(void) {
	int32_t value;
	if(random_r(&physical_rng, &value)) fail("Physical RNG sampling");
	return value;
}
static void physical_srand(unsigned int seed) {
	if(srandom_r(seed, &physical_rng)) fail("Physical RNG seed");
}
#define rand physical_rand
#define srand physical_srand
#endif

static const char *basename_of(const char *path) {
	const char *slash = strrchr(path, '/');
	return slash ? slash+1 : path;
}

static void safe_prefix(char *result, const char *path) {
	if(strlen(path) >= PFC_PATH || !*basename_of(path)) fail("Invalid output prefix");
	strcpy(result, path);
	char *p = result+(basename_of(path)-path);
	for(; *p; p++) {
		unsigned char c = (unsigned char)*p;
		if(!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
		     (c >= '0' && c <= '9') || c == '_' || c == '-')) *p = '_';
	}
}

static void make_path(char *result, const char *prefix, const char *suffix) {
	if(snprintf(result, PFC_PATH, "%s%s", prefix, suffix) >= PFC_PATH) fail("Path too long");
}

// Only the documented old spelling may fall back, in the same directory.
static FILE *open_restart(const char *filename) {
	FILE *file = fopen(filename, "r");
	if(file || errno != ENOENT) return file;
	const char *base = basename_of(filename);
	const char *marker = NULL;
	for(const char *p = base; (p = strstr(p, "-t-")); p++) marker = p;
	if(!marker || marker == base) return NULL;
	const char *digits = marker+3, *end = digits;
	while(*end >= '0' && *end <= '9') end++;
	if(end == digits || strcmp(end, ".dat")) return NULL;
	while(digits+1 < end && *digits == '0') digits++;
	char prefix[PFC_PATH], safe[PFC_PATH], mapped[PFC_PATH];
	size_t length = (size_t)(marker-filename);
	if(length >= sizeof(prefix)) fail("Restart path too long");
	memcpy(prefix, filename, length);
	prefix[length] = '\0';
	safe_prefix(safe, prefix);
	if(snprintf(mapped, sizeof(mapped), "%s-%s", safe, digits) >= (int)sizeof(mapped)) fail("Restart path too long");
	file = fopen(mapped, "r");
	if(file) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(!rank) printf("PFC_RESTART requested=%s actual=%s\n", filename, mapped);
	}
	return file;
}

static void reserve_output(const char *prefix, int rank) {
	if(rank == 0) {
		char directory[PFC_PATH], filename[PFC_PATH];
		const char *base = basename_of(prefix);
		size_t directory_length = (size_t)(base-prefix);
		if(directory_length) {
			memcpy(directory, prefix, directory_length);
			directory[directory_length] = '\0';
		} else strcpy(directory, ".");
		DIR *dir = opendir(directory);
		if(!dir) fail(directory);
		struct dirent *entry;
		while((entry = readdir(dir))) {
			if(strncmp(entry->d_name, base, strlen(base)) || entry->d_name[strlen(base)] != '-') continue;
			const char *p = entry->d_name+strlen(base)+1, *start = p;
			while(*p >= '0' && *p <= '9') p++;
			if(p != start && !strcmp(p, ".dat")) fail("Output data already exists for this prefix");
		}
		if(closedir(dir)) fail(directory);
		make_path(filename, prefix, ".out");
		FILE *file = fopen(filename, "wx");
		if(!file) fail(filename);
		if(fclose(file)) fail(filename);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

// contains stuff related to output
struct Output {
	char name[PFC_PATH];	// safe output prefix; input path stays unchanged
	int T_print;		// interval for printing output
	int T_write;		// interval for saving state
	int *written_steps;
	size_t written_count;
};

// contains stuff related to data arrays
struct Arrays {
	int W;				// system width
	int H;				// system height
	ptrdiff_t lH;		// local system height
	ptrdiff_t lW;		// local system width
	ptrdiff_t lh0;		// local vertical start index
	ptrdiff_t lw0;		// local horizontal start index
	
	double *A;		// operator for linear part, e^{-k^2 \hat{\mathcal{L}} \Delta t}
	double *B;		// operator for nonlinear part ...
	double *p;		// array for \psi_C
	double *q;		// another array
	
	fftw_plan p_P;	// FFTW plan for F(p)
	fftw_plan q_Q;	// F(q)
	fftw_plan Q_q;	// F^-1(Q)
#ifdef PFC_DCU
	PfcDcu *dcu;
#endif
};

// contains model parameters
struct Model {
	double alpha;
	double beta;
	double gamma;
	double delta;
	double c;
};

// contains stuff related to relaxation
struct Relaxation {
	time_t t0;			// start time
	int id;				// rank of process
	int ID;				// number of processes
	int t;				// time step count
	double f;			// average free energy density
	double p;			// average density
	double d;			// sampling step size for calculation box optimization
	int T;				// total number of iterations
	double dx;			// x-discretization
	double dy;			// y-
	double dt;			// time step
	int T_optimize;		// interval for calculation box optimization
	int block;
};

// seeds processes' random number generators
void seed_rngs(struct Relaxation *relaxation, FILE *input) {
	int seed;								// seed for next process
	// seed for 0th process from input file
	if(fscanf(input, " %d", &seed) != 1) {
		printf("Invalid input for random number generator seed!\n");
		fail("Invalid solver input or restart data");
	}
	if(relaxation->id == 0) {
		if(seed == 0) srand(time(NULL));	// random seed
		else srand(seed);					// user-specified seed
		seed = rand();						// sample new seed for next process
	}
	else {
		MPI_Recv(&seed, 1, MPI_INT, relaxation->id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	// receive from previous
		srand(seed);																			// seed
		seed = rand();																			// new seed for next
	}
	if(relaxation->id != relaxation->ID-1) MPI_Send(&seed, 1, MPI_INT, relaxation->id+1, 0, MPI_COMM_WORLD);	// send to next
}

// configures arrays and FFTW plans
void configure_arrays(struct Arrays *arrays, FILE *input) {
	if(fscanf(input, " %d %d", &arrays->W, &arrays->H) != 2 || arrays->W <= 0 || arrays->W > INT_MAX-2 || arrays->H <= 0) {
		printf("Invalid input for array dimensions!\n");
		fail("Invalid solver input or restart data");
	}

	ptrdiff_t lWHh = fftw_mpi_local_size_2d_transposed(arrays->H, arrays->W/2+1, MPI_COMM_WORLD, &arrays->lH, &arrays->lh0, &arrays->lW, &arrays->lw0);
	if(lWHh < 0 || lWHh > PTRDIFF_MAX/2) fail("Invalid FFTW allocation size");
	ptrdiff_t lWHp = 2*lWHh;	// preserve FFTW's allocation capacity and slab boundaries
	if(lWHh == 0) lWHh = 1, lWHp = 2;
	if((size_t)lWHp > SIZE_MAX/sizeof(double)) fail("Host allocation byte size overflow");

	arrays->A = (double*)fftw_malloc(lWHh*sizeof(double));	// allocate arrays
	arrays->B = (double*)fftw_malloc(lWHh*sizeof(double));
	
	arrays->p = (double*)fftw_malloc(lWHp*sizeof(double));
	arrays->q = (double*)fftw_malloc(lWHp*sizeof(double));
	if(!arrays->A || !arrays->B || !arrays->p || !arrays->q) fail("Host array allocation");
	memset(arrays->p, 0, lWHp*sizeof(double));
	memset(arrays->q, 0, lWHp*sizeof(double));
#ifdef PFC_DCU
	arrays->dcu = pfc_dcu_create(arrays->W, arrays->H, arrays->lH, arrays->lh0, arrays->lW, arrays->lw0, MPI_COMM_WORLD);
#else
	arrays->p_P = fftw_mpi_plan_dft_r2c_2d(arrays->H, arrays->W, arrays->p, (fftw_complex*)arrays->p, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
	arrays->q_Q = fftw_mpi_plan_dft_r2c_2d(arrays->H, arrays->W, arrays->q, (fftw_complex*)arrays->q, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
	arrays->Q_q = fftw_mpi_plan_dft_c2r_2d(arrays->H, arrays->W, (fftw_complex*)arrays->q, arrays->q, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);
	if(!arrays->p_P || !arrays->q_Q || !arrays->Q_q) fail("FFTW planning");
#endif
}

// configures output
void configure_output(struct Output *output, FILE *input) {
	// read intervals for printing and saving from input file
	if(fscanf(input, " %d %d", &output->T_print, &output->T_write) != 2 || output->T_print <= 0 || output->T_write <= 0) {
		printf("Invalid input for output configuration!\n");
		fail("Invalid solver input or restart data");
	}
}

// configures model
void configure_model(struct Model *model, FILE *input) {
	// read model parameters from input file
	if(fscanf(input, " %lf %lf %lf %lf %lf", &model->alpha, &model->beta, &model->gamma, &model->delta, &model->c) != 5) {
		printf("Invalid input for model parameters!\n");
		fail("Invalid solver input or restart data");
	}
	if(!isfinite(model->alpha) || !isfinite(model->beta) || !isfinite(model->gamma) ||
	   !isfinite(model->delta) || !isfinite(model->c)) fail("Non-finite model parameter");
}

// configures relaxation
void configure_relaxation(struct Relaxation *relaxation, FILE *input) {
	// read from input file
	if(fscanf(input, " %d %lf %lf %lf %d", &relaxation->T, &relaxation->dx, &relaxation->dy, &relaxation->dt, &relaxation->T_optimize) != 5 ||
	   relaxation->T < 0 || relaxation->T_optimize < 0 ||
	   !isfinite(relaxation->dx) || !isfinite(relaxation->dy) || !isfinite(relaxation->dt) ||
	   relaxation->dx <= 0 || relaxation->dy <= 0 || relaxation->dt <= 0) {
		printf("Invalid input for relaxation settings!\n");
		fail("Invalid solver input or restart data");
	}
}

// initializes the density field with random noise
// p and A are the average density and noise amplitude
void random_state(struct Arrays *arrays, double p, double A) {
	int Wp = 2*(arrays->W/2+1);
	int w, h, k;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			// average + amplitude * random from U(-0.5, 0.5)
			arrays->q[k+w] = p+A*((double)rand()/RAND_MAX-0.5);
		}
	}
}

// returns the one-mode approximation for close-packed lattice
// x0, y0 give rotation (theta) axis; u, v give translation; l0 dimensionless lattice constant
double OMA(double x0, double y0, double u, double v, double l0, double theta) {
	double qx, qy;
	qx = 2.0*pi/l0;						// wave numbers
	qy = qx/sqrt(3.0);
	double cosine = cos(theta);
	double sine = sin(theta);
	double x = cosine*x0-sine*y0+u;		// rotation matrix applied
	double y = sine*x0+cosine*y0+v;
	return cos(qx*x)*cos(qy*y)+0.5*cos(2.0*qy*y);
}

// initializes the density field with random crystallites
// dx, dy and l0 control the lenght scale
// p and A are the average density and noise amplitude
// N is number of crystallites and R their radius
void polycrystalline_state(struct Arrays *arrays, double dx, double dy, double l0, double p, double A, int N, double R) {
	double crystallites[3*N];	// array for crystallite coordinates and orientations
	int w, h, gh, k, n, id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if(id == 0) {						   // 0th process samples crystallites
		for(n = 0; n < N; n++) {
			crystallites[3*n] = (double)rand()/RAND_MAX*arrays->W;
			crystallites[3*n+1] = (double)rand()/RAND_MAX*arrays->H;
			crystallites[3*n+2] = 2.0*pi*rand()/RAND_MAX;
		}
	}
	// broadcasted to other processes
	MPI_Bcast(crystallites, 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int Wp = 2*(arrays->W/2+1);
	// x-, y- and squared distances from crystallite center, Gaussian envelope
	double u, v, r2, G;
	double div = 0.5/(R*R*l0*l0);
	int min;	// index for closest crystallite
	double r2min;	// squared distance for closest crystallite
	double umin = 0.0;	// u and v for closest crystallite
	double vmin = 0.0;
	for(h = 0; h < arrays->lH; h++) {
		gh = arrays->lh0+h;	// global vertical index
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			arrays->q[k+w] = p;	// (1) set to average density
			min = -1;
			r2min = 1.0e100;
			for(n = 0; n < N; n++) {
				u = w-crystallites[3*n];
				// closest periodic image searched
				if(u > 0.5*arrays->W) u -= arrays->W;
				else if(u < -0.5*arrays->W) u += arrays->W;
				u *= dx;
				v = gh-crystallites[3*n+1];
				if(v > 0.5*arrays->H) v -= arrays->H;
				else if(v < -0.5*arrays->H) v += arrays->H;
				v *= dy;
				r2 = u*u+v*v;
				if(r2 < r2min) {	// closest crystallite so far?
					min = n;
					r2min = r2;
					umin = u;
					vmin = v;
				}
			}
			G = exp(-r2min*div);
			// (2) superimpose crystallites
			arrays->q[k+w] += A*G*OMA(umin, vmin, 0.0, 0.0, l0, crystallites[3*min+2]);
		}
	}
}

// initializes the density field with data from a data file
void read_state(struct Arrays *arrays, FILE *file, double p, double A) {
	int Wp = 2*(arrays->W/2+1);
	int w, h, k;
	double p0 = 0.0;
	double Amin = 1.0e100;
	double Amax = -1.0e100;
	int rank, ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ranks);
	size_t local_count = (size_t)arrays->W*arrays->lH;
	if(local_count > INT_MAX) fail("Restart slab exceeds MPI count range");
	long long local[] = {(long long)arrays->lH, (long long)arrays->lh0};
	long long *slabs = rank == 0 ? malloc((size_t)2*ranks*sizeof(long long)) : NULL;
	if(rank == 0 && !slabs) fail("Restart slab metadata allocation");
	if(MPI_Gather(local, 2, MPI_LONG_LONG, slabs, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fail("Restart slab metadata communication");
	char *text = NULL;
	size_t *starts = NULL;
	int *lengths = NULL;
	int parallel = 0;
	if(rank == 0) {
		if((size_t)arrays->W > SIZE_MAX/(size_t)arrays->H/sizeof(double)) fail("Restart field byte size overflow");
		size_t total = (size_t)arrays->W*arrays->H;
		for(int r = 0; r < ranks; r++) {
			if(slabs[2*r] < 0 || slabs[2*r] > arrays->H || slabs[2*r+1] < 0 || slabs[2*r+1] > arrays->H-slabs[2*r])
				fail("Invalid restart slab");
			if((size_t)arrays->W*slabs[2*r] > INT_MAX) fail("Invalid restart slab extent");
		}
		if(fseek(file, 0, SEEK_END) == 0) {
			long end = ftell(file);
			if(end < 0 || (uintmax_t)end > SIZE_MAX-1) fail("Restart text byte size overflow");
			size_t bytes = (size_t)end;
			if(fseek(file, 0, SEEK_SET)) fail("Restart file rewind");
			// shortcut: root buffers restart ASCII, or FP64 on fallback; stream if it exceeds root memory.
			text = malloc(bytes+1);
			starts = malloc((size_t)ranks*sizeof(size_t));
			lengths = malloc((size_t)ranks*sizeof(int));
			if(!text || !starts || !lengths) fail("Restart text allocation");
			if(fread(text, 1, bytes, file) != bytes || ferror(file)) fail("Restart text read");
			text[bytes] = '\0';
			if(fseek(file, 0, SEEK_SET)) fail("Restart file rewind");
			size_t cursor = 0, covered = 0;
			parallel = 1;
			for(int r = 0; r < ranks && parallel; r++) {
				size_t count = (size_t)arrays->W*slabs[2*r];
				if(count && (size_t)arrays->W*slabs[2*r+1] != covered) { parallel = 0; break; }
				starts[r] = cursor;
				for(size_t i = 0; i < count; i++) {
					while(cursor < bytes && isspace((unsigned char)text[cursor])) cursor++;
					if(cursor == bytes) { parallel = 0; break; }
					while(cursor < bytes && !isspace((unsigned char)text[cursor])) cursor++;
				}
				if(cursor-starts[r] > INT_MAX) parallel = 0;
				else lengths[r] = (int)(cursor-starts[r]);
				covered += count;
			}
			if(covered != total) parallel = 0;
		} else clearerr(file);
	}
	if(MPI_Bcast(&parallel, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fail("Restart parser selection communication");
	if(parallel) {
		int bytes;
		if(MPI_Scatter(lengths, 1, MPI_INT, &bytes, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			fail("Restart text lengths communication");
		char *chunk = rank == 0 ? text+starts[0] : malloc((size_t)bytes+1);
		if(!chunk) fail("Restart local text allocation");
		if(rank == 0) {
			for(int r = 1; r < ranks; r++)
				if(MPI_Send(text+starts[r], lengths[r], MPI_BYTE, r, 74, MPI_COMM_WORLD) != MPI_SUCCESS)
					fail("Restart text distribution");
		} else if(MPI_Recv(chunk, bytes, MPI_BYTE, 0, 74, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
			fail("Restart text reception");
		int valid = 1;
		if(local_count) {
			FILE *input = fmemopen(chunk, (size_t)bytes, "r");
			if(!input) fail("Restart text stream");
			for(size_t i = 0; i < local_count; i++)
				if(fscanf(input, "%lf", &arrays->p[i]) != 1 || !isfinite(arrays->p[i])) { valid = 0; break; }
			int next;
			do { next = fgetc(input); } while(next != EOF && isspace((unsigned char)next));
			if(next != EOF) valid = 0;
			if(ferror(input) || fclose(input)) fail("Restart local text stream read");
		}
		if(rank != 0) free(chunk);
		// scanf also accepts adjacent numbers such as 1-2; retain the original parser for those files.
		if(MPI_Allreduce(&valid, &parallel, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
			fail("Restart parser validation communication");
	}
	free(text);
	free(starts);
	free(lengths);
	if(rank == 0) printf("PFC_RESTART_PARSE mode=%s\n", parallel ? "parallel" : "root-fallback");
	if(!parallel && rank == 0) {
		size_t total = (size_t)arrays->W*arrays->H;
		double *field = malloc(total*sizeof(double));
		if(!field) fail("Restart field allocation");
		for(size_t i = 0; i < total; i++)
			if(fscanf(file, "%lf", &field[i]) != 1 || !isfinite(field[i])) fail("Invalid restart density");
		for(int r = 0; r < ranks; r++) {
			if(slabs[2*r] < 0 || slabs[2*r+1] < 0) fail("Invalid restart slab");
			size_t count = (size_t)arrays->W*slabs[2*r];
			size_t offset = (size_t)arrays->W*slabs[2*r+1];
			if(count > INT_MAX || offset > total || count > total-offset) fail("Invalid restart slab extent");
			if(r == 0) memcpy(arrays->p, field+offset, count*sizeof(double));
			else if(MPI_Send(field+offset, (int)count, MPI_DOUBLE, r, 73, MPI_COMM_WORLD) != MPI_SUCCESS)
				fail("Restart slab distribution");
		}
		free(field);
	} else if(!parallel && MPI_Recv(arrays->p, (int)local_count, MPI_DOUBLE, 0, 73, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
		fail("Restart slab reception");
	free(slabs);
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			arrays->q[k+w] = arrays->p[(size_t)h*arrays->W+w];
			p0 += arrays->q[k+w];
			if(arrays->q[k+w] < Amin) Amin = arrays->q[k+w];
			if(arrays->q[k+w] > Amax) Amax = arrays->q[k+w];
		}
	}
	if(arrays->lH == 0) return;
	if(!(Amax > Amin)) fail("Constant restart slab cannot be rescaled by the original algorithm");
	p0 /= arrays->W*arrays->lH;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			arrays->q[k+w] = (arrays->q[k+w] - p0)/(Amax - Amin)*A + p;
		}
	}
}

// initializes the density field
void initialize_system(struct Arrays *arrays, FILE *input) {
	int init;
	if(fscanf(input, " %d", &init) != 1) {
		printf("Invalid initialization type!\n");
		fail("Invalid solver input or restart data");
	}
	if(init == 0) {	// random initialization
		double p, A;
		if(fscanf(input, " %lf %lf", &p, &A) != 2) {
			printf("Invalid initialization parameters!\n");
			fail("Invalid solver input or restart data");
		}
		if(!isfinite(p) || !isfinite(A)) fail("Non-finite initial density parameter");
		random_state(arrays, p, A);
	}
	else if(init == 1) {	// polycrystalline initialization
		int N;
		double dx, dy, l0, p, A, R;
		if(fscanf(input, " %lf %lf %lf %lf %lf %d %lf", &dx, &dy, &l0, &p, &A, &N, &R) != 7) {
			printf("Invalid initialization parameters!\n");
			fail("Invalid solver input or restart data");
		}
		if(N <= 0 || N > INT_MAX/3 || R <= 0 || l0 <= 0 || dx <= 0 || dy <= 0 ||
		   !isfinite(p) || !isfinite(A) || !isfinite(R) || !isfinite(l0) || !isfinite(dx) || !isfinite(dy))
			fail("Invalid crystallite parameters");
		polycrystalline_state(arrays, dx, dy, l0, p, A, N, R);
	}
	else if(init == 2) {			// initial state read from a data file
		char filename[PFC_PATH];
		double p, A;
		if(fscanf(input, " %4095s %lf %lf", filename, &p, &A) != 3) {
			printf("Invalid initialization parameters!\n");
			fail("Invalid solver input or restart data");
		}
		if(!isfinite(p) || !isfinite(A)) fail("Non-finite restart parameter");
		FILE *file = open_restart(filename);
		if(file == NULL) {
			printf("Data file not found!\n");
			fail("Invalid solver input or restart data");
		}
		read_state(arrays, file, p, A);
		if(fclose(file)) fail(filename);
	}
	else fail("Unknown initialization type");
#ifdef PFC_DCU
	pfc_dcu_upload(arrays->dcu, arrays->q);
#endif
}

// saves state into a file
void write_state(struct Arrays *arrays, struct Relaxation *relaxation, struct Output *output) {
	char filename[PFC_PATH], suffix[64];
	snprintf(suffix, sizeof(suffix), "-%d.dat", relaxation->t);
	make_path(filename, output->name, suffix);
#ifdef PFC_DCU
	pfc_dcu_download(arrays->dcu, arrays->q);
#endif
	int owned = 0;
	for(size_t j = 0; j < output->written_count; j++) if(output->written_steps[j] == relaxation->t) owned = 1;
	if(!owned) {
		int *steps = realloc(output->written_steps, (output->written_count+1)*sizeof(int));
		if(!steps) fail("Output ownership allocation");
		output->written_steps = steps;
		output->written_steps[output->written_count++] = relaxation->t;
	}
	if((size_t)arrays->lH > SIZE_MAX/(size_t)arrays->W) fail("Output sample count overflow");
	size_t samples = (size_t)arrays->lH*arrays->W;
	if(samples > (SIZE_MAX-1)/32) fail("Output formatting size overflow");
	size_t capacity = 32*samples+1, used = 0;
	char *text = malloc(capacity);
	if(!text) fail("Output formatting allocation");
	int Wp = 2*(arrays->W/2+1);
	// Format the same %e values in parallel, then retain the original rank write order.
	for(ptrdiff_t h = 0; h < arrays->lH; h++) {
		for(int w = 0; w < arrays->W; w++) {
			int length = snprintf(text+used, capacity-used, "%e\n", arrays->q[Wp*h+w]);
			if(length < 0 || (size_t)length >= capacity-used) fail("Output formatting failed");
			used += length;
		}
	}
	for(int i = 0; i < relaxation->ID; i++) {
		if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) fail("Output MPI_Barrier");
		if(relaxation->id == i) {
			FILE *file = fopen(filename, i == 0 ? (owned ? "w" : "wx") : "a");
			if(!file) fail(filename);
			if(fwrite(text, 1, used, file) != used || ferror(file) || fclose(file)) fail(filename);
		}
	}
	free(text);
}

// prints output, format:
// time step count, elapsed wall-clock time (s), dx (dimensionless*), dy (*), average free energy density (*), average density (*)
void print(struct Relaxation *relaxation, struct Output *output) {
	if(relaxation->id == 0) {			// only 0th process
		printf("%d %d %.9lf %.9lf %.9lf %.9lf\n", relaxation->t, (int)(time(NULL)-relaxation->t0), relaxation->dx, relaxation->dy, relaxation->f, relaxation->p);
		char filename[PFC_PATH];
		make_path(filename, output->name, ".out");
		FILE *file;
		file = fopen(filename, "a");	// need only append, empty file was already generated
		if(!file) fail(filename);
		fprintf(file, "%d %d %.9lf %.9lf %.9lf %.9lf\n", relaxation->t, (int)(time(NULL)-relaxation->t0), relaxation->dx, relaxation->dy, relaxation->f, relaxation->p);
		if(ferror(file) || fclose(file)) fail(filename);
	}
}

// updates operators for linear and nonlinear parts
void update_AB(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	int w, h, gw, k;
	double ky, kx2, k2;									// (squared) k-vector (components)
	double divWH = 1.0/arrays->W/arrays->H;				// 1/(W*H)
	double dkx = 2.0*pi/(relaxation->dx*arrays->W);		// discretization in k-space
	double dky = 2.0*pi/(relaxation->dy*arrays->H);
	double d2, l, expl;									// (1+\nabla^2)^2, variables for intermediate results
	double c = model->c;
	for(w = 0; w < arrays->lW; w++) {
		gw = arrays->lw0+w;								// global x-index
		kx2 = gw*dkx;									// kx
		kx2 *= kx2;										// kx^2
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			if(h < arrays->H/2) ky = h*dky;				// upper half negative ky-values
			else ky = (h-arrays->H)*dky;
			k2 = kx2+ky*ky;								// k^2
			d2 = (1.0-k2);
			d2 *= d2;
			l = model->alpha+model->beta*d2;			// alpha+beta*(1+\nabla^2)^2 in k-space
			expl = exp(-((1.0-c)+c*k2)*l*relaxation->dt);			// e^{-k^2 \hat{\mathcal{L}} \Delta t}
			arrays->A[k+h] = expl;
			if(l == 0.0) arrays->B[k+h] = -((1.0-c)+c*k2)*relaxation->dt;		// avoid divide by zero
			else arrays->B[k+h] = (expl-1.0)/l;
			arrays->B[k+h] *= divWH;					// descaling
		}
	}
#ifdef PFC_DCU
	pfc_dcu_upload_ab(arrays->dcu, arrays->A, arrays->B);
#endif
}

// scales p arrays (in k-space) by 1/(W*H) ((I)FFTs cause scaling of data by sqrt(W*H))
void scale_P(struct Arrays *arrays) {
	fftw_complex *P;			// complex data pointers
	P = (fftw_complex*)&arrays->p[0];
	double divWH = 1.0/arrays->W/arrays->H;
	int i;
	int lA = arrays->lW*arrays->H;
	for(i = 0; i < lA; i++) P[i] *= divWH;
}

// scales q arrays (in k-space) by 1/(W*H)
void scale_Q(struct Arrays *arrays) {
	fftw_complex *Q;
	Q = (fftw_complex*)&arrays->q[0];
	double divWH = 1.0/arrays->W/arrays->H;
	int i;
	int lA = arrays->lW*arrays->H;
	for(i = 0; i < lA; i++) Q[i] *= divWH;
}

// computes average free energy density (f) and average density (p) (f & p -> fp)
void fp(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	int Wp = 2*(arrays->W/2+1);
	int w, h, k;
#ifdef PFC_DCU
	pfc_dcu_diagnostic(arrays->dcu, arrays->p, arrays->q, relaxation->dx, relaxation->dy);
#else
	double dkx = 2.0*pi/(relaxation->dx*arrays->W);
	double dky = 2.0*pi/(relaxation->dy*arrays->H);
	memcpy(arrays->p, arrays->q, Wp*arrays->lH*sizeof(double));		// save current state to p array
	fftw_execute(arrays->q_Q);
	scale_Q(arrays);
	double kx2, ky, k2, d2;
	fftw_complex *Q = (fftw_complex*)&arrays->q[0];
	int gw;
	for(w = 0; w < arrays->lW; w++) {
		gw = arrays->lw0+w;
		kx2 = gw*dkx;
		kx2 *= kx2;
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			if(h < arrays->H/2) ky = h*dky;
			else ky = (h-arrays->H)*dky;
			k2 = kx2+ky*ky;
			d2 = (1.0-k2);
			d2 *= d2;
			Q[k+h] *= d2;		// compute (1-\nabla^2)^2 \psi_C in k-space
		}
	}
	fftw_execute(arrays->Q_q);
#endif
	relaxation->f = 0.0;		// reset variables for average (free energy) density
	relaxation->p = 0.0;
	double div3 = 1.0/3.0;
	double p, q, p2;				// variables to simplify expressions
	double alpha = model->alpha;	// variables for model parameters to simplify expressions
	double beta = model->beta;
	double gamma = model->gamma;
	double delta = model->delta;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			p = arrays->p[k+w];
			q = arrays->q[k+w];
			p2 = p*p;
			// average free energy density
			relaxation->f += 0.5*( alpha*p2 + beta*p*q ) + div3*gamma*p*p2 + 0.25*delta*p2*p2;
			// restore q arrays
			arrays->q[k+w] = p;
			// average densities
			relaxation->p += p;
		}
	}
	// communication between processes, note that can't reduce as average because chunks may have different sizes!
	MPI_Allreduce(MPI_IN_PLACE, &relaxation->f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &relaxation->p, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double divWH = 1.0/arrays->W/arrays->H;
	relaxation->f *= divWH;
	relaxation->p *= divWH;
#ifndef PFC_DCU
	fftw_execute(arrays->p_P);
	scale_P(arrays);
#endif
}

// performs one iteration of the semi-implicit spectral method
void step(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
#ifdef PFC_DCU
	pfc_dcu_step(arrays->dcu, relaxation->t, model->gamma, model->delta);
#else
	int w, h, k;
	int Wp = 2*(arrays->W/2+1);
	if(relaxation->t%100 == 0) {	// avoids numerical instability (I guess there's some accumulating numerical error, don't know really...)
		memcpy(arrays->p, arrays->q, Wp*arrays->lH*sizeof(double));
		fftw_execute(arrays->p_P);	// update k-space solution fully from real-space solution
		scale_P(arrays);
	}
	double q, q2;
	double gamma = model->gamma;
	double delta = model->delta;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			q = arrays->q[k+w];
			q2 = q*q;
			// compute nonlinear part (C)
			arrays->q[k+w] = gamma*q2 + delta*q*q2;
		}
	}
	fftw_execute(arrays->q_Q);
	fftw_complex *P = (fftw_complex*)&arrays->p[0];
	fftw_complex *Q = (fftw_complex*)&arrays->q[0];
	for(w = 0; w < arrays->lW; w++) {
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			P[k+h] = arrays->A[k+h]*P[k+h] + arrays->B[k+h]*Q[k+h];
			// copy p's to q's
			Q[k+h] = P[k+h];
		}
	}
	fftw_execute(arrays->Q_q);
#endif
}

// optimizes calculation box size for system
// samples slightly different box sizes by varying dx and dy and interpolates optimum quadratically
void optimize(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	// sample sizes
	double dxs[] = {relaxation->dx, relaxation->dx-relaxation->d, relaxation->dx+relaxation->d, relaxation->dx, relaxation->dx};
	double dys[] = {relaxation->dy, relaxation->dy, relaxation->dy, relaxation->dy-relaxation->d, relaxation->dy+relaxation->d};
	double fs[5];	// average free energy densities of sizes sampled
	int i;
	for(i = 0; i < 5; i++) {		// sample sizes
		relaxation->dx = dxs[i];
		relaxation->dy = dys[i];
		fp(arrays, model, relaxation);				// compute energy (and density (not needed here though))
		fs[i] = relaxation->f;	// save average free energy density
	}
	// interpolate new dx and dy (minimum of surface f = a+b*dx+c*dy+d*dx^2+e*dy^2)
	relaxation->dx = (relaxation->d*(fs[1]-fs[2])+2.0*dxs[0]*(-2.0*fs[0]+fs[1]+fs[2]))/(2.0*(fs[1]-2.0*fs[0]+fs[2]));
	relaxation->dy = (relaxation->d*(fs[3]-fs[4])+2.0*dys[0]*(-2.0*fs[0]+fs[3]+fs[4]))/(2.0*(fs[3]-2.0*fs[0]+fs[4]));
	// check that change in discretization is acceptable
	double l0 = 4.0*pi/sqrt(3.0);						// approximate dimensionless lenght scale (lattice constant)
	double dw = arrays->W*(relaxation->dx-dxs[0]);		// change in horizontal system size
	double dh = arrays->H*(relaxation->dy-dys[0]);		// ... vertical ...
	double dr = sqrt(dw*dw+dh*dh);						// "change vector"
	double x = 0.25*l0;									// limit to 1/4 of lattice constant (to ensure stability)
	if(dr > x) {	// if the change in system dimensions exceeds 1/4 of the lattice constant ...
		x /= dr;
		relaxation->dx = x*relaxation->dx+(1.0-x)*dxs[0];	// ... truncate the change to 1/4 of the lattice constant
		relaxation->dy = x*relaxation->dy+(1.0-x)*dys[0];
	}
	// update A and B
	update_AB(arrays, model, relaxation);	// dx and dy changed -> need to update operators
	// update sampling step size (tries to keep it in the same ballpark with dr (for, hopefully, more accurate optimization))
	double ddx = relaxation->dx-dxs[0];
	double ddy = relaxation->dy-dys[0];
	double ddr = sqrt(ddx*ddx+ddy*ddy);					// discretization change
	if(ddr < relaxation->d) relaxation->d *= 0.5;		// if change vector < d, halve d
	else relaxation->d *= 2.0;							// otherwise double
	if(relaxation->d < 1.0e-6) relaxation->d *= 2.0;	// can cause numerical issues if d gets too small
}

// relaxes the system for T time steps
void relax(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation, struct Output *output) {
	size_t count = (size_t)relaxation->T+1;
	if(count > 715827882) fail("Too many steps for buffered MPI timing");
	double *timing = calloc(3*count, sizeof(double));
	if(!timing) fail("Step timing allocation");
	double *elapsed = timing, *cumulative = timing+count, *optimization = timing+2*count;
	double block_start = MPI_Wtime();
	const char *phase = basename_of(output->name);
	int timed_rank = 1;
#ifdef PFC_DCU
	timed_rank = pfc_dcu_is_leader(arrays->dcu);
#endif
	for(relaxation->t = 0; relaxation->t <= relaxation->T; relaxation->t++) {
		if(relaxation->T_optimize > 0 && relaxation->t > 0 && relaxation->t%relaxation->T_optimize == 0) {
			double start = MPI_Wtime();
			optimize(arrays, model, relaxation);
#ifdef PFC_DCU
			pfc_dcu_sync(arrays->dcu);
#endif
			if(timed_rank) optimization[relaxation->t] = MPI_Wtime()-start;
		}
		if(relaxation->t%output->T_print == 0) {	// print output?
			fp(arrays, model, relaxation);
			print(relaxation, output);
		}
		if(relaxation->t%output->T_write == 0) write_state(arrays, relaxation, output);		// write out state?
		if(relaxation->t < relaxation->T) {
			double start = MPI_Wtime();
			step(arrays, model, relaxation);
#ifdef PFC_DCU
			pfc_dcu_sync(arrays->dcu);
#endif
			if(timed_rank) {
				elapsed[relaxation->t] = MPI_Wtime()-start;
				cumulative[relaxation->t] = MPI_Wtime()-block_start;
			}
		}
	}
	// One block reduction observes every completed step without changing step communication.
	MPI_Allreduce(MPI_IN_PLACE, timing, (int)(3*count), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(relaxation->id == 0) {
		for(int t = 0; t <= relaxation->T; t++) {
			if(relaxation->T_optimize > 0 && t > 0 && t%relaxation->T_optimize == 0)
				printf("PFC_TIMING kind=optimize phase=%s block=%d step=%d seconds=%.9f\n", phase, relaxation->block, t, optimization[t]);
			if(t < relaxation->T)
				printf("PFC_TIMING kind=step phase=%s block=%d step=%d seconds=%.9f cumulative_seconds=%.9f\n", phase, relaxation->block, t+1, elapsed[t], cumulative[t]);
		}
		fflush(stdout);
	}
	free(timing);
}

// frees allocated arrays
void clear_arrays(struct Arrays *arrays) {
#ifdef PFC_DCU
	if(arrays->dcu) pfc_dcu_destroy(arrays->dcu);
#else
	if(arrays->p_P) fftw_destroy_plan(arrays->p_P);
	if(arrays->q_Q) fftw_destroy_plan(arrays->q_Q);
	if(arrays->Q_q) fftw_destroy_plan(arrays->Q_q);
#endif
	fftw_free(arrays->p);
	fftw_free(arrays->q);
	fftw_free(arrays->A);
	fftw_free(arrays->B);
}

int main(int argc, char **argv) {
	// init MPI
	MPI_Init(&argc, &argv);
	double invocation_start = MPI_Wtime();
	fftw_mpi_init();
#ifdef PFC_DCU
	if(initstate_r(1, (char*)physical_rng_state, sizeof(physical_rng_state), &physical_rng))
		fail("Physical RNG state initialization");
#endif
	// create structs
	struct Arrays arrays = {0};
	struct Model model = {0};
	struct Relaxation relaxation = {0};
	struct Output output = {0};
	relaxation.t0 = time(NULL);
	relaxation.t = 0;
	relaxation.d = 0.0001;
	MPI_Comm_rank(MPI_COMM_WORLD, &relaxation.id);
	MPI_Comm_size(MPI_COMM_WORLD, &relaxation.ID);
	if(argc != 2) fail("Usage: polygrpfc INPUT_PREFIX");
	safe_prefix(output.name, argv[1]);
	// input stream
	char filename[PFC_PATH];
	make_path(filename, argv[1], ".in");
	FILE *input = fopen(filename, "r");
	if(input == NULL) fail(filename);
	reserve_output(output.name, relaxation.id);
	if(relaxation.id == 0 && strcmp(output.name, argv[1]))
		printf("PFC_PREFIX requested=%s actual=%s\n", argv[1], output.name);
	// read input
	char label;
	char line[1024];
	int initialized = 0, model_configured = 0;
	while(fscanf(input, " %c", &label) != EOF) {
		if(label == '#' && fgets(line, 1024, input)) {}			// comment
		else if(label == '\n') {}								// empty line
		else if(label == 'S') {									// seed random number generators
			seed_rngs(&relaxation, input);
		}
		else if(label == 'O') {									// output
			configure_output(&output, input);
		}
		else if(label == 'A') {									// set up arrays and FFTW plans
			if(arrays.q) { clear_arrays(&arrays); memset(&arrays, 0, sizeof(arrays)); }
			configure_arrays(&arrays, input);
			initialized = 0;
		}
		else if(label == 'I') {									// initialization
			if(!arrays.q) fail("Arrays must precede initialization");
			initialize_system(&arrays, input);
			initialized = 1;
		}
		else if(label == 'M') {									// model
			configure_model(&model, input);
			model_configured = 1;
		}
		else if(label == 'R') {									// relaxation
			if(!initialized || !model_configured || output.T_print <= 0 || output.T_write <= 0)
				fail("Relaxation needs initialized arrays, model and output intervals");
			configure_relaxation(&relaxation, input);
			update_AB(&arrays, &model, &relaxation);
			if(relaxation.block++ == 0) {
#ifdef PFC_DCU
				pfc_dcu_sync(arrays.dcu);
#endif
				double initialization = MPI_Wtime()-invocation_start;
				MPI_Allreduce(MPI_IN_PLACE, &initialization, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				if(relaxation.id == 0) printf("PFC_TIMING kind=init phase=%s seconds=%.9f\n", basename_of(output.name), initialization);
			}
			relax(&arrays, &model, &relaxation, &output);
		}
		else {
			printf("Invalid label: %c!\n", label);				// bad input
			fail("Invalid input label");
		}
	}
	if(ferror(input) || fclose(input)) fail(filename);
	
	clear_arrays(&arrays);										// clean-up
	free(output.written_steps);
	double total = MPI_Wtime()-invocation_start;
	MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(relaxation.id == 0) {
		printf("PFC_TIMING kind=total phase=%s seconds=%.9f\n", basename_of(output.name), total);
		fflush(stdout);
	}
	MPI_Finalize();

	return 0;

}
