// P7 measurement build of libalfsb.so (copy of Libraries/Modules/alf_la_sycl_blk.cpp; same kernels, same expressions; NOT part of the source tree).
// LD_PRELOAD it in front of the production libalfsb.so.  With ALFP7_PROF=1 every operation of the asynchronous block path is run on a
// SEPARATE in-order queue PQ (same context and device, enable_profiling) after a wait on the interop queue, and waited for; its device time
// (command_end - command_start of the event, 1 ns ticks) and its host time (submit + wait) are accumulated per operation and printed at exit:
//   ALFP7_BLK <op> calls <n> dev_s <s> host_s <s> mean_global_items <w>
// Without ALFP7_PROF the code is the original (async on the interop queue).  The numerical result is the same (fully serialised).
// gemm: oneMKL returns the event of its last kernel only, so dev_s is a lower bound for gemm; host_s is the submit-to-complete time.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <chrono>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <string>
using cd = std::complex<double>;
namespace {
sycl::queue *Q = nullptr, *PQ = nullptr;
struct Acc { long long n = 0; double dev = 0, host = 0, wi = 0; };
std::map<std::string, Acc> g_acc;
std::mutex g_mx;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int g_sub = -1;
bool sub_on() { if (g_sub < 0) { const char *e = std::getenv("ALFP7_SUBMIT"); g_sub = (e && e[0] == '1') ? 1 : 0; } return g_sub == 1; }
std::map<std::string, Acc> g_sacc;
struct Report {
  ~Report() {
    for (auto &kv : g_sacc)
      std::fprintf(stderr, "ALFP7_SUBMIT %-10s calls %9lld submit_s %10.4f\n", kv.first.c_str(), kv.second.n, kv.second.host);
    if (!PQ) return;
    for (auto &kv : g_acc)
      std::fprintf(stderr, "ALFP7_BLK %-8s calls %9lld dev_s %10.4f host_s %10.4f mean_global_items %10.1f\n", kv.first.c_str(), kv.second.n, kv.second.dev,
                   kv.second.host, kv.second.n ? kv.second.wi / kv.second.n : 0.0);
  }
} g_report;
template <class F> int guarded(const char *what, F &&f) {
  try { f(); return 0; }
  catch (sycl::exception &e) { std::fprintf(stderr, "alfsb: %s: SYCL exception: %s\n", what, e.what()); }
  catch (std::exception &e) { std::fprintf(stderr, "alfsb: %s: exception: %s\n", what, e.what()); }
  return -1;
}
// submit(queue&) -> sycl::event
template <class F> int op(const char *name, double items, F &&submit) {
  return guarded(name, [&] {
    if (!PQ) {
      if (sub_on()) { double t0 = now(); submit(*Q); double h = now() - t0; std::lock_guard<std::mutex> l(g_mx); auto &a = g_sacc[name]; a.n++; a.host += h; }
      else submit(*Q);
      return;
    }
    Q->wait();
    double t0 = now();
    sycl::event e = submit(*PQ);
    e.wait();
    double h = now() - t0;
    double d = 1e-9 * (double)(e.get_profiling_info<sycl::info::event_profiling::command_end>() - e.get_profiling_info<sycl::info::event_profiling::command_start>());
    std::lock_guard<std::mutex> l(g_mx);
    auto &a = g_acc[name]; a.n++; a.dev += d; a.host += h; a.wi += items;
  });
}
}  // namespace

extern "C" {
int alfsb_init(void *q) {
  if (!q) return -1;
  Q = static_cast<sycl::queue *>(q);
  const char *e = std::getenv("ALFP7_PROF");
  if (e && e[0] == '1') {
    PQ = new sycl::queue(Q->get_context(), Q->get_device(), sycl::property_list{sycl::property::queue::in_order{}, sycl::property::queue::enable_profiling{}});
    std::fprintf(stderr, "ALFP7: block path serialised on a profiling queue (interop queue in_order=%d)\n", (int)Q->is_in_order());
  }
  return Q->is_in_order() ? 0 : 1;
}
void *alfsb_alloc_host(int64_t bytes) { try { return sycl::malloc_host(bytes, *Q); } catch (...) { return nullptr; } }
void alfsb_free_host(void *p) { if (p && Q) sycl::free(p, *Q); }
int alfsb_wait() { return guarded("wait", [&] { Q->wait_and_throw(); if (PQ) PQ->wait_and_throw(); }); }
int alfsb_copy(void *dst, const void *src, int64_t bytes) {
  return op("copy", (double)bytes, [&](sycl::queue &q) { return q.memcpy(dst, src, (size_t)bytes); });
}
int alfsb_gather(int nfl, int nfa, int nfb, int ld, int b, int nd, cd *gbk, const int *bsk, const cd *Gk) {
  return op("gather", (double)nfl * b * b, [&](sycl::queue &q) {
    return q.parallel_for(sycl::range<3>(nfl, b, b), [=](sycl::id<3> id) {
      long fl = id[0], c = id[1], a = id[2];
      long nf = (fl == 1) ? nfb : nfa;
      gbk[a + (long)nd * (c + (long)nd * fl)] = Gk[(bsk[a + (long)nd * fl] - 1) + (long)ld * ((bsk[c + (long)nd * fl] - 1) + (long)ld * (nf - 1))];
    });
  });
}
int alfsb_scale(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd *Gk, const int *bsk, const cd *cbk) {
  return op("scale", (double)nfl * b * ld, [&](sycl::queue &q) {
    return q.parallel_for(sycl::range<3>(nfl, b, ld), [=](sycl::id<3> id) {
      long fl = id[0], a = id[1] + 1, j = id[2] + 1;
      long nf = (fl == 1) ? nfb : nfa;
      long sa = bsk[(a - 1) + (long)nd * fl];
      long cj = 0;
      for (long c = 1; c <= b; ++c)
        if (bsk[(c - 1) + (long)nd * fl] == j) cj = c;
      cd *Gf = Gk + (long)ld * ld * (nf - 1);
      cd v = cbk[(a - 1) + (long)ldc * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cbk[(nd + cj - 1) + (long)ldc * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cbk[(nd + a - 1) + (long)ldc * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
int alfsb_fgather(const cd *G, int ld, int nrow, const int *S, int k, int nf, cd *C, cd *Qt, int *slot) {
  return op("fgather", (double)k * nrow, [&](sycl::queue &q) {
    return q.parallel_for(sycl::range<2>(k, nrow), [=](sycl::id<2> id) {
      long m = id[0], i = id[1];
      long s = S[m];
      C[i + (long)ld * m] = G[i + (long)ld * (s - 1)];
      cd qv = -G[(s - 1) + (long)ld * i];
      if (i + 1 == s) { qv = qv + 1.0; slot[i + (long)ld * (nf - 1)] = 0; }
      Qt[i + (long)ld * m] = qv;
    });
  });
}
int alfsb_gemm(int ta, int tb, int m, int n, int k, double ar, double ai, const cd *A, int lda, const cd *B, int ldb, double br, double bi, cd *C, int ldc) {
  auto t = [](int x) { return x == 0 ? oneapi::mkl::transpose::nontrans : (x == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans); };
  char nm[32]; std::snprintf(nm, sizeof nm, "gemm m%d", m);
  return op(nm, 8.0 * m * n * k, [&](sycl::queue &q) {
    return oneapi::mkl::blas::column_major::gemm(q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, B, ldb, cd(br, bi), C, ldc);
  });
}
}
