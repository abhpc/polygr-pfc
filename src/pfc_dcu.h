#ifndef PFC_DCU_H
#define PFC_DCU_H

#include <stddef.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PfcDcu PfcDcu;

/* FP64 throughout; failures abort all ranks in the supplied MPI communicator.
 * Density buffers remain caller-owned: lH rows of 2*(W/2+1) doubles.
 * A/B buffers use the original lW*H transposed slab. All ranks call transfer
 * functions; step/sync execute on card leaders without a per-step barrier.
 * Transfers complete before returning; device q/P persist until destroy.
 * Slabs are the original fftw_mpi_local_size_2d_transposed distribution. */
PfcDcu *pfc_dcu_create(int W, int H, ptrdiff_t lH, ptrdiff_t lh0,
	ptrdiff_t lW, ptrdiff_t lw0, MPI_Comm comm);
void pfc_dcu_destroy(PfcDcu *ctx);
void pfc_dcu_upload(PfcDcu *ctx, const double *host);
void pfc_dcu_download(PfcDcu *ctx, double *host);
void pfc_dcu_upload_ab(PfcDcu *ctx, const double *A, const double *B);
int pfc_dcu_is_leader(const PfcDcu *ctx);
void pfc_dcu_step(PfcDcu *ctx, int t, double gamma, double delta);
/* Padded density and (1+laplacian)^2 density; restores q and normalized P. */
void pfc_dcu_diagnostic(PfcDcu *ctx, double *host_density,
	double *host_operator, double dx, double dy);
void pfc_dcu_sync(PfcDcu *ctx);

#ifdef __cplusplus
}
#endif

#endif
