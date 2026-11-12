// P7 measurement library (LD_PRELOAD, application only): interposes the exported USM complex<double> ZGEMM and ZTRSM of libmkl_sycl_blas (the
// oneMKL SYCL API; the Fortran "!$omp dispatch" calls and the block path both end up there; MKL_VERBOSE lists them as
// oneapi::mkl::blas::column_major::gemm[complex<double>]).  With ALFP7_PROF=1 each call is run on a separate in-order PROFILING queue (same context
// and device) after a wait on the caller's queue, and waited for; its device time (event command_end - command_start, the last kernel of the call)
// and host time (submit + wait) are accumulated per shape and printed at exit:
//   ALFP7_MKL <routine shape> calls <n> dev_s <s> host_s <s>
// Calls on a queue that already profiles (the block library's queue) are passed through.  Fully serialised; numerics unchanged.
// Without ALFP7_PROF the calls are passed through unchanged.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <dlfcn.h>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <string>
#include <vector>
using cd = std::complex<double>;
using namespace oneapi::mkl;
namespace {
struct Acc { long long n = 0; double dev = 0, host = 0; };
std::map<std::string, Acc> g_acc;
std::mutex g_mx;
sycl::queue *g_pq = nullptr;
int g_on = -1;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
bool on() { if (g_on < 0) { const char *e = std::getenv("ALFP7_PROF"); g_on = (e && e[0] == '1') ? 1 : 0; } return g_on == 1; }
long long g_seen = 0, g_passq = 0, g_profq = 0;
struct Report {
  ~Report() { std::fprintf(stderr, "ALFP7_MKL_SEEN calls %lld on_profiling_queue %lld\n", g_seen, g_profq); for (auto &kv : g_acc) std::fprintf(stderr, "ALFP7_MKL %-30s calls %9lld dev_s %10.4f host_s %10.4f\n", kv.first.c_str(), kv.second.n, kv.second.dev, kv.second.host); }
} g_rep;
template <class F> sycl::event prof(sycl::queue &q, const std::string &key, F &&call) {
  std::lock_guard<std::mutex> l(g_mx);
  g_seen++;
  if (!on()) return call(q);
  // a queue that already profiles (the block library's queue, or the interop queue under LIBOMPTARGET_PLUGIN_PROFILE=T): run on it
  bool pq_given = q.has_property<sycl::property::queue::enable_profiling>();
  if (pq_given) g_profq++;
  sycl::queue *use = &q;
  if (!pq_given) {
    if (!g_pq) g_pq = new sycl::queue(q.get_context(), q.get_device(), sycl::property_list{sycl::property::queue::in_order{}, sycl::property::queue::enable_profiling{}});
    q.wait();
    use = g_pq;
  } else g_passq++;
  double t0 = now();
  sycl::event e = call(*use);
  e.wait();
  double h = now() - t0;
  double d = 1e-9 * (double)(e.get_profiling_info<sycl::info::event_profiling::command_end>() - e.get_profiling_info<sycl::info::event_profiling::command_start>());
  auto &a = g_acc[key]; a.n++; a.dev += d; a.host += h;
  return e;
}
std::string gk(long m, long n, long k, transpose a, transpose b) { char s[96]; std::snprintf(s, sizeof s, "zgemm %ldx%ldx%ld t%d%d", m, n, k, (int)a, (int)b); return s; }
std::string tk(side sd, uplo ul, transpose t, long m, long n) { char s[96]; std::snprintf(s, sizeof s, "ztrsm %ldx%ld s%d u%d t%d", m, n, (int)sd, (int)ul, (int)t); return s; }
using VP = value_or_pointer<cd>;
using EV = std::vector<sycl::event>;
}  // namespace

#define SYM_GEMM "_ZN6oneapi3mkl4blas12column_major4gemmERN4sycl3_V15queueENS0_9transposeES7_lllNS0_16value_or_pointerISt7complexIdEEEPKSA_lSD_lSB_PSA_lRKSt6vectorINS4_5eventESaISG_EE"
#define SYM_GEMMM "_ZN6oneapi3mkl4blas12column_major4gemmERN4sycl3_V15queueENS0_9transposeES7_lllNS0_16value_or_pointerISt7complexIdEEEPKSA_lSD_lSB_PSA_lNS1_12compute_modeERKSt6vectorINS4_5eventESaISH_EE"
#define SYM_TRSM "_ZN6oneapi3mkl4blas12column_major4trsmERN4sycl3_V15queueENS0_4sideENS0_4uploENS0_9transposeENS0_4diagEllNS0_16value_or_pointerISt7complexIdEEEPKSD_lPSD_lRKSt6vectorINS4_5eventESaISJ_EE"
#define SYM_TRSMM "_ZN6oneapi3mkl4blas12column_major4trsmERN4sycl3_V15queueENS0_4sideENS0_4uploENS0_9transposeENS0_4diagEllNS0_16value_or_pointerISt7complexIdEEEPKSD_lPSD_lNS1_12compute_modeERKSt6vectorINS4_5eventESaISK_EE"

sycl::event p7_gemm(sycl::queue &, transpose, transpose, long, long, long, VP, const cd *, long, const cd *, long, VP, cd *, long, const EV &) __asm__(SYM_GEMM);
sycl::event p7_gemm(sycl::queue &q, transpose ta, transpose tb, long m, long n, long k, VP al, const cd *A, long lda, const cd *B, long ldb, VP be, cd *C, long ldc, const EV &dep) {
  using Fn = sycl::event (*)(sycl::queue &, transpose, transpose, long, long, long, VP, const cd *, long, const cd *, long, VP, cd *, long, const EV &);
  static Fn real = (Fn)dlsym(RTLD_NEXT, SYM_GEMM);
  return prof(q, gk(m, n, k, ta, tb), [&](sycl::queue &qq) { return real(qq, ta, tb, m, n, k, al, A, lda, B, ldb, be, C, ldc, dep); });
}
sycl::event p7_gemmm(sycl::queue &, transpose, transpose, long, long, long, VP, const cd *, long, const cd *, long, VP, cd *, long, blas::compute_mode, const EV &) __asm__(SYM_GEMMM);
sycl::event p7_gemmm(sycl::queue &q, transpose ta, transpose tb, long m, long n, long k, VP al, const cd *A, long lda, const cd *B, long ldb, VP be, cd *C, long ldc, blas::compute_mode cm, const EV &dep) {
  using Fn = sycl::event (*)(sycl::queue &, transpose, transpose, long, long, long, VP, const cd *, long, const cd *, long, VP, cd *, long, blas::compute_mode, const EV &);
  static Fn real = (Fn)dlsym(RTLD_NEXT, SYM_GEMMM);
  return prof(q, gk(m, n, k, ta, tb), [&](sycl::queue &qq) { return real(qq, ta, tb, m, n, k, al, A, lda, B, ldb, be, C, ldc, cm, dep); });
}
sycl::event p7_trsm(sycl::queue &, side, uplo, transpose, diag, long, long, VP, const cd *, long, cd *, long, const EV &) __asm__(SYM_TRSM);
sycl::event p7_trsm(sycl::queue &q, side sd, uplo ul, transpose t, diag dg, long m, long n, VP al, const cd *A, long lda, cd *B, long ldb, const EV &dep) {
  using Fn = sycl::event (*)(sycl::queue &, side, uplo, transpose, diag, long, long, VP, const cd *, long, cd *, long, const EV &);
  static Fn real = (Fn)dlsym(RTLD_NEXT, SYM_TRSM);
  return prof(q, tk(sd, ul, t, m, n), [&](sycl::queue &qq) { return real(qq, sd, ul, t, dg, m, n, al, A, lda, B, ldb, dep); });
}
sycl::event p7_trsmm(sycl::queue &, side, uplo, transpose, diag, long, long, VP, const cd *, long, cd *, long, blas::compute_mode, const EV &) __asm__(SYM_TRSMM);
sycl::event p7_trsmm(sycl::queue &q, side sd, uplo ul, transpose t, diag dg, long m, long n, VP al, const cd *A, long lda, cd *B, long ldb, blas::compute_mode cm, const EV &dep) {
  using Fn = sycl::event (*)(sycl::queue &, side, uplo, transpose, diag, long, long, VP, const cd *, long, cd *, long, blas::compute_mode, const EV &);
  static Fn real = (Fn)dlsym(RTLD_NEXT, SYM_TRSMM);
  return prof(q, tk(sd, ul, t, m, n), [&](sycl::queue &qq) { return real(qq, sd, ul, t, dg, m, n, al, A, lda, B, ldb, cm, dep); });
}
