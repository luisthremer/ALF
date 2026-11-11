// Gate A' (Phase 3b): SYCL submission on the sycl::queue that the OpenMP runtime hands out through omp interop
// (exactly as alf_la_om_mod::om_init does). Small kernels + oneMKL SYCL gemm, called one op at a time from Fortran.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <cstdio>
using cd = std::complex<double>;
static sycl::queue *Q = nullptr;
static double *V; static cd *A, *B, *C; static int N, K;
#define GUARD(stmt) try { stmt; } catch (sycl::exception &e) { std::fprintf(stderr, "SB sycl exception: %s (code %d)\n", e.what(), e.code().value()); std::fflush(stderr); std::abort(); } catch (std::exception &e) { std::fprintf(stderr, "SB exception: %s\n", e.what()); std::fflush(stderr); std::abort(); }
extern "C" {
int sb_init(void *q) { Q = static_cast<sycl::queue *>(q); return Q->is_in_order() ? 1 : 0; }
void sb_set(double *v, cd *a, cd *b, cd *c, int n, int k) { V = v; A = a; B = b; C = c; N = n; K = k; }
void sb_kern(double c, double d) {
  double *v = V;
  GUARD(Q->parallel_for(sycl::range<1>(N), [=](sycl::id<1> i) { int jj = (int)i[0] + 1; v[i] = v[i] * c + d * (1.0 + 1.e-4 * jj); }));
}
void sb_gemm() {
  GUARD(oneapi::mkl::blas::column_major::gemm(*Q, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans, N, K, K,
                                        cd(1, 0), A, N, B, K, cd(0, 0), C, N));
}
void sb_copy() {
  double *v = V; cd *a = A, *c = C; int n = N;
  GUARD(Q->parallel_for(sycl::range<2>(K, N), [=](sycl::id<2> id) {
    size_t jj = id[0], ii = id[1];
    cd x = c[jj * n + ii];
    a[jj * n + ii] = cd(x.real() * 0.9 + 1.e-3 * v[ii], x.imag() * 0.9);
  }));
}
void sb_wait() { GUARD(Q->wait_and_throw()); }
}
