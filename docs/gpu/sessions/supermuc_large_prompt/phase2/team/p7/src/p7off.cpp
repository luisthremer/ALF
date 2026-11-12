// P7 measurement library (LD_PRELOAD, application only): interposes the oneMKL OpenMP-offload VARIANTS that the Fortran "!$omp dispatch" resolves to
//   mkl_blas_zgemm_omp_offload_lp64, mkl_blas_ztrsm_omp_offload_lp64, mkl_blas_zgemm_batch_strided_omp_offload_lp64
// (signature = the BLAS one + the appended interop object, see $MKLROOT/include/mkl_blas_omp_variant.h).  With ALFP7_PROF=1 the call is NOT forwarded; instead
// the same oneMKL SYCL USM routine (oneapi::mkl::blas::column_major::gemm / trsm / gemm_batch, the routine the variant calls) is run on a separate in-order
// PROFILING queue with the same arguments (the arguments are already device pointers at this point), after waiting for the interop queue and the
// ALF shim queue, and waited for.  Device time = event command_end - command_start of its last kernel; host time = submit + wait.  Printed at exit:
//   ALFP7_OFF <routine shape> calls <n> dev_s <s> host_s <s>
// Without ALFP7_PROF the call is forwarded to the real variant and only counted.  The result is the same call, fully serialised.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <omp.h>
#include <dlfcn.h>
#include <execinfo.h>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <string>
using cd = std::complex<double>;
using namespace oneapi::mkl;
namespace {
struct Acc { long long n = 0; double dev = 0, host = 0; };
std::map<std::string, Acc> g_acc;
std::mutex g_mx;
sycl::queue *g_pq = nullptr, *g_shimq = nullptr;
int g_on = -1;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
bool on() { if (g_on < 0) { const char *e = std::getenv("ALFP7_PROF"); g_on = (e && e[0] == '1') ? 1 : 0; } return g_on == 1; }
struct Report {
  ~Report() { for (auto &kv : g_acc) std::fprintf(stderr, "ALFP7_OFF %-34s calls %9lld dev_s %10.4f host_s %10.4f\n", kv.first.c_str(), kv.second.n, kv.second.dev, kv.second.host); }
} g_rep;
transpose tr(char c) { return (c == 'N' || c == 'n') ? transpose::nontrans : ((c == 'T' || c == 't') ? transpose::trans : transpose::conjtrans); }
// call site: the first return address inside the main executable (text below 0x10000000, ALF.out is not PIE); resolved offline with nm -n
std::string callsite() {
  void *bt[24]; int n = backtrace(bt, 24); std::string r; int k = 0;
  for (int i = 1; i < n && k < 4; ++i) if ((unsigned long)bt[i] < 0x10000000ul) { char b[24]; std::snprintf(b, sizeof b, "%s%lx", k ? "," : " @", (unsigned long)bt[i]); r += b; ++k; }
  return k ? r : " @0";
}
template <class F> void prof(void *interop, const std::string &key0, F &&call) {
  const std::string key = key0 + callsite();
  int rc = 0;
  sycl::queue *iq = nullptr;
  if (interop) iq = static_cast<sycl::queue *>(omp_get_interop_ptr((omp_interop_t)interop, omp_ipr_targetsync, &rc));
  std::lock_guard<std::mutex> l(g_mx);
  if (!g_shimq) { auto f = (void *(*)())dlsym(RTLD_DEFAULT, "alfsh_get_queue"); if (f) g_shimq = static_cast<sycl::queue *>(f()); }
  sycl::queue *base = g_shimq ? g_shimq : iq;
  if (!g_pq) g_pq = new sycl::queue(base->get_context(), base->get_device(), sycl::property_list{sycl::property::queue::in_order{}, sycl::property::queue::enable_profiling{}});
  if (iq && rc == 0) iq->wait();
  if (g_shimq) g_shimq->wait();
  double t0 = now();
  sycl::event e = call(*g_pq);
  e.wait();
  double h = now() - t0;
  double d = 1e-9 * (double)(e.get_profiling_info<sycl::info::event_profiling::command_end>() - e.get_profiling_info<sycl::info::event_profiling::command_start>());
  auto &a = g_acc[key]; a.n++; a.dev += d; a.host += h;
}
}  // namespace
extern "C" {
void mkl_blas_zgemm_omp_offload_lp64(const char *ta, const char *tb, const int *m, const int *n, const int *k, const cd *alpha, const cd *a, const int *lda, const cd *b,
                                     const int *ldb, const cd *beta, cd *c, const int *ldc, void *interop) {
  if (!on()) {
    static auto real = (void (*)(const char *, const char *, const int *, const int *, const int *, const cd *, const cd *, const int *, const cd *, const int *, const cd *, cd *, const int *, void *))dlsym(RTLD_NEXT, "mkl_blas_zgemm_omp_offload_lp64");
    real(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, interop); return;
  }
  char key[96]; std::snprintf(key, sizeof key, "zgemm %dx%dx%d t%c%c", *m, *n, *k, *ta, *tb);
  prof(interop, key, [&](sycl::queue &q) { return blas::column_major::gemm(q, tr(*ta), tr(*tb), *m, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc); });
}
void mkl_blas_ztrsm_omp_offload_lp64(const char *sd, const char *ul, const char *ta, const char *dg, const int *m, const int *n, const cd *alpha, const cd *a, const int *lda, cd *b,
                                     const int *ldb, void *interop) {
  if (!on()) {
    static auto real = (void (*)(const char *, const char *, const char *, const char *, const int *, const int *, const cd *, const cd *, const int *, cd *, const int *, void *))dlsym(RTLD_NEXT, "mkl_blas_ztrsm_omp_offload_lp64");
    real(sd, ul, ta, dg, m, n, alpha, a, lda, b, ldb, interop); return;
  }
  char key[96]; std::snprintf(key, sizeof key, "ztrsm %dx%d %c%c%c%c", *m, *n, *sd, *ul, *ta, *dg);
  prof(interop, key, [&](sycl::queue &q) {
    return blas::column_major::trsm(q, (*sd == 'L' || *sd == 'l') ? side::left : side::right, (*ul == 'U' || *ul == 'u') ? uplo::upper : uplo::lower, tr(*ta),
                                    (*dg == 'N' || *dg == 'n') ? diag::nonunit : diag::unit, *m, *n, *alpha, a, *lda, b, *ldb);
  });
}
void mkl_blas_zgemm_batch_strided_omp_offload_lp64(const char *ta, const char *tb, const int *m, const int *n, const int *k, const cd *alpha, const cd *a, const int *lda, const long *sa,
                                                   const cd *b, const int *ldb, const long *sb, const cd *beta, cd *c, const int *ldc, const long *sc, const int *bs, void *interop) {
  if (!on()) {
    static auto real = (void (*)(const char *, const char *, const int *, const int *, const int *, const cd *, const cd *, const int *, const long *, const cd *, const int *, const long *, const cd *, cd *, const int *, const long *, const int *, void *))dlsym(RTLD_NEXT, "mkl_blas_zgemm_batch_strided_omp_offload_lp64");
    real(ta, tb, m, n, k, alpha, a, lda, sa, b, ldb, sb, beta, c, ldc, sc, bs, interop); return;
  }
  char key[96]; std::snprintf(key, sizeof key, "zgemm_batch_strided %dx%dx%d b%d t%c%c", *m, *n, *k, *bs, *ta, *tb);
  prof(interop, key, [&](sycl::queue &q) { return blas::column_major::gemm_batch(q, tr(*ta), tr(*tb), *m, *n, *k, *alpha, a, *lda, *sa, b, *ldb, *sb, *beta, c, *ldc, *sc, *bs); });
}
}
