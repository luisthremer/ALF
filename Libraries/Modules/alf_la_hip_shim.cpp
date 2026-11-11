// ALF GPU=amd (ALF_LA_HIP): C shim of the alf_la HIP backend (module alf_la_hip_mod).
// BLAS through hipBLAS, solvers through hipSOLVER, on the null stream. Every call finishes (hipDeviceSynchronize) before it
// returns: the OpenMP target regions of the offload runtime (libgomp GCN plugin) run on their own queues, so a library call
// has to be ordered against them explicitly. All pointers are device pointers (from use_device_addr).
#include <hip/hip_runtime_api.h>
#include <hipblas/hipblas.h>
#include <hipsolver/hipsolver.h>
#include <cstdint>
#include <cstdio>

static hipblasHandle_t hb = nullptr;
static hipsolverHandle_t hs = nullptr;
static int* dinfo = nullptr;

static hipblasOperation_t bop(int o) { return o == 0 ? HIPBLAS_OP_N : (o == 1 ? HIPBLAS_OP_T : HIPBLAS_OP_C); }
static hipsolverOperation_t sop(int o) { return o == 0 ? HIPSOLVER_OP_N : (o == 1 ? HIPSOLVER_OP_T : HIPSOLVER_OP_C); }
static hipDoubleComplex* Z(void* p) { return static_cast<hipDoubleComplex*>(p); }
static hipDoubleComplex mkz(double r, double i) { return make_hipDoubleComplex(r, i); }

extern "C" {

int alfhp_init(int dev, char* name, int nname) {
  { hipError_t e = hipSetDevice(dev); if (e != hipSuccess) { fprintf(stderr, "alfhp_init: hipSetDevice(%d): %s\n", dev, hipGetErrorString(e)); return 1; } }
  hipDeviceProp_t pr;
  if (hipGetDeviceProperties(&pr, dev) != hipSuccess) return 2;
  snprintf(name, nname, "%s (%s)", pr.name, pr.gcnArchName);
  if (hipblasCreate(&hb) != HIPBLAS_STATUS_SUCCESS) return 3;
  if (hipsolverCreate(&hs) != HIPSOLVER_STATUS_SUCCESS) return 4;
  if (hipMalloc((void**)&dinfo, sizeof(int)) != hipSuccess) return 5;
  return 0;
}
void alfhp_finalize() {
  if (dinfo) hipFree(dinfo);
  if (hs) hipsolverDestroy(hs);
  if (hb) hipblasDestroy(hb);
  dinfo = nullptr; hs = nullptr; hb = nullptr;
}
int alfhp_wait() { return hipDeviceSynchronize() == hipSuccess ? 0 : 1; }

#define SYNC_RET(st) do { if ((st) != 0) return 1000 + (int)(st); return hipDeviceSynchronize() == hipSuccess ? 0 : 999; } while (0)

int alfhp_zgemm(int ta, int tb, int m, int n, int k, double ar, double ai, void* A, int lda, void* B, int ldb,
                double br, double bi, void* C, int ldc) {
  hipDoubleComplex al = mkz(ar, ai), be = mkz(br, bi);
  int st = hipblasZgemm(hb, bop(ta), bop(tb), m, n, k, &al, Z(A), lda, Z(B), ldb, &be, Z(C), ldc);
  SYNC_RET(st);
}
int alfhp_dgemm(int m, int n, int k, void* A, int lda, void* B, int ldb, void* C, int ldc) {
  double one = 1.0, zero = 0.0;
  int st = hipblasDgemm(hb, HIPBLAS_OP_N, HIPBLAS_OP_N, m, n, k, &one, (double*)A, lda, (double*)B, ldb, &zero, (double*)C, ldc);
  SYNC_RET(st);
}
int alfhp_zgemm_batched(int ta, int tb, int m, int n, int k, double ar, double ai, void* A, int lda, int64_t sA, void* B,
                        int ldb, int64_t sB, double br, double bi, void* C, int ldc, int64_t sC, int nb) {
  hipDoubleComplex al = mkz(ar, ai), be = mkz(br, bi);
  int st = hipblasZgemmStridedBatched(hb, bop(ta), bop(tb), m, n, k, &al, Z(A), lda, sA, Z(B), ldb, sB, &be, Z(C), ldc, sC, nb);
  SYNC_RET(st);
}
int alfhp_zgeru(int m, int n, double ar, double ai, void* x, int incx, void* y, int incy, void* A, int lda) {
  hipDoubleComplex al = mkz(ar, ai);
  int st = hipblasZgeru(hb, m, n, &al, Z(x), incx, Z(y), incy, Z(A), lda);
  SYNC_RET(st);
}
int alfhp_ztrsm(int side, int uplo, int op, int diag, int m, int n, double ar, double ai, void* A, int lda, void* B, int ldb) {
  hipDoubleComplex al = mkz(ar, ai);
  hipblasSideMode_t s = side == 0 ? HIPBLAS_SIDE_LEFT : HIPBLAS_SIDE_RIGHT;
  hipblasFillMode_t u = uplo == 1 ? HIPBLAS_FILL_MODE_UPPER : HIPBLAS_FILL_MODE_LOWER;
  hipblasDiagType_t d = diag == 1 ? HIPBLAS_DIAG_UNIT : HIPBLAS_DIAG_NON_UNIT;
  int st = hipblasZtrsm(hb, s, u, bop(op), d, m, n, &al, Z(A), lda, Z(B), ldb);
  SYNC_RET(st);
}

// ---- hipSOLVER: work sizes in complex(8) elements; *_ returns the device info (0 = ok) or >= 1000 on an API failure ----
static int fin(hipsolverStatus_t st) {
  if (st != HIPSOLVER_STATUS_SUCCESS) return 1000 + (int)st;
  if (hipDeviceSynchronize() != hipSuccess) return 999;
  int h = 0;
  if (hipMemcpy(&h, dinfo, sizeof(int), hipMemcpyDeviceToHost) != hipSuccess) return 998;
  return h;
}
int alfhp_zgeqrf_ws(int64_t m, int64_t n, int64_t lda, int64_t* ws) {
  int l = 0;
  if (hipsolverZgeqrf_bufferSize(hs, (int)m, (int)n, nullptr, (int)lda, &l) != HIPSOLVER_STATUS_SUCCESS) return 1;
  *ws = l; return 0;
}
int alfhp_zgeqrf(int64_t m, int64_t n, void* A, int64_t lda, void* tau, void* work, int64_t lwork) {
  return fin(hipsolverZgeqrf(hs, (int)m, (int)n, Z(A), (int)lda, Z(tau), Z(work), (int)lwork, dinfo));
}
int alfhp_zungqr_ws(int64_t m, int64_t n, int64_t k, int64_t lda, int64_t* ws) {
  int l = 0;
  if (hipsolverZungqr_bufferSize(hs, (int)m, (int)n, (int)k, nullptr, (int)lda, nullptr, &l) != HIPSOLVER_STATUS_SUCCESS) return 1;
  *ws = l; return 0;
}
int alfhp_zungqr(int64_t m, int64_t n, int64_t k, void* A, int64_t lda, void* tau, void* work, int64_t lwork) {
  return fin(hipsolverZungqr(hs, (int)m, (int)n, (int)k, Z(A), (int)lda, Z(tau), Z(work), (int)lwork, dinfo));
}
int alfhp_zunmqr_ws(int side, int trans, int64_t m, int64_t n, int64_t k, int64_t lda, int64_t ldc, int64_t* ws) {
  int l = 0;
  hipsolverSideMode_t s = side == 0 ? HIPSOLVER_SIDE_LEFT : HIPSOLVER_SIDE_RIGHT;
  if (hipsolverZunmqr_bufferSize(hs, s, sop(trans), (int)m, (int)n, (int)k, nullptr, (int)lda, nullptr, nullptr, (int)ldc, &l) !=
      HIPSOLVER_STATUS_SUCCESS) return 1;
  *ws = l; return 0;
}
int alfhp_zunmqr(int side, int trans, int64_t m, int64_t n, int64_t k, void* A, int64_t lda, void* tau, void* C, int64_t ldc,
                 void* work, int64_t lwork) {
  hipsolverSideMode_t s = side == 0 ? HIPSOLVER_SIDE_LEFT : HIPSOLVER_SIDE_RIGHT;
  return fin(hipsolverZunmqr(hs, s, sop(trans), (int)m, (int)n, (int)k, Z(A), (int)lda, Z(tau), Z(C), (int)ldc, Z(work), (int)lwork, dinfo));
}
int alfhp_zgetrf_ws(int64_t m, int64_t n, int64_t lda, int64_t* ws) {
  int l = 0;
  if (hipsolverZgetrf_bufferSize(hs, (int)m, (int)n, nullptr, (int)lda, &l) != HIPSOLVER_STATUS_SUCCESS) return 1;
  *ws = l; return 0;
}
int alfhp_zgetrf(int64_t m, int64_t n, void* A, int64_t lda, void* ipiv, void* work, int64_t lwork) {
  return fin(hipsolverZgetrf(hs, (int)m, (int)n, Z(A), (int)lda, Z(work), (int)lwork, (int*)ipiv, dinfo));
}
int alfhp_zgetrs_ws(int trans, int64_t n, int64_t nrhs, int64_t lda, int64_t ldb, int64_t* ws) {
  int l = 0;
  if (hipsolverZgetrs_bufferSize(hs, sop(trans), (int)n, (int)nrhs, nullptr, (int)lda, nullptr, nullptr, (int)ldb, &l) !=
      HIPSOLVER_STATUS_SUCCESS) return 1;
  *ws = l; return 0;
}
int alfhp_zgetrs(int trans, int64_t n, int64_t nrhs, void* A, int64_t lda, void* ipiv, void* B, int64_t ldb, void* work, int64_t lwork) {
  return fin(hipsolverZgetrs(hs, sop(trans), (int)n, (int)nrhs, Z(A), (int)lda, (int*)ipiv, Z(B), (int)ldb, Z(work), (int)lwork, dinfo));
}

}  // extern "C"
