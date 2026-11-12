// Phase 3b step 2: the device sequence of one submatrix block (blk_begin gather + blk_end scale + 2 x flush) in SYCL, on the
// OpenMP interop queue. Kernels are line by line translations of k_blk_gather, k_blk_scale, k_sub_flush_gather.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <cstdio>
using cd = std::complex<double>;
using oneapi::mkl::transpose;
static sycl::queue *Q = nullptr;
#define GUARD(stmt) try { stmt; } catch (sycl::exception &e) { std::fprintf(stderr, "BB sycl exception: %s\n", e.what()); std::fflush(stderr); std::abort(); }
extern "C" {
int bb_init(void *q) { Q = static_cast<sycl::queue *>(q); return Q->is_in_order() ? 1 : 0; }
void *bb_alloc_host(long bytes) { return sycl::malloc_host(bytes, *Q); }
void bb_wait() { GUARD(Q->wait_and_throw()); }
void bb_memcpy(void *dst, const void *src, long bytes) { GUARD(Q->memcpy(dst, src, bytes)); }
// gbk(nd,nd,2) = Gk(ld,ld,nfg)(bsk(a,fl), bsk(c,fl), nf)  ; nf = fl (flavours 1,2)
void bb_gather(int nfl, int ld, int b, int nd, cd *gbk, const int *bsk, const cd *Gk) {
  GUARD(Q->parallel_for(sycl::range<3>(nfl, b, b), [=](sycl::id<3> id) {
    long fl = id[0], c = id[1], a = id[2];
    long nf = fl;
    gbk[a + (long)nd * (c + (long)nd * fl)] = Gk[(bsk[a + (long)nd * fl] - 1) + (long)ld * ((bsk[c + (long)nd * fl] - 1) + (long)ld * nf)];
  }));
}
void bb_scale(int nfl, int ld, int b, int nd, int ldc, cd *Gk, const int *bsk, const cd *cbk) {
  GUARD(Q->parallel_for(sycl::range<3>(nfl, b, ld), [=](sycl::id<3> id) {
    long fl = id[0], a1 = id[1] + 1, j = id[2] + 1;
    long nf = fl;
    long sa = bsk[(a1 - 1) + (long)nd * fl];
    long cj = 0;
    for (long c = 1; c <= b; ++c) if (bsk[(c - 1) + (long)nd * fl] == j) cj = c;
    cd v = cbk[(a1 - 1) + (long)ldc * fl] * Gk[(sa - 1) + (long)ld * ((j - 1) + (long)ld * nf)];
    if (cj > 0) v = cbk[(nd + cj - 1) + (long)ldc * fl] * v;
    Gk[(sa - 1) + (long)ld * ((j - 1) + (long)ld * nf)] = v;
    if (cj == 0) Gk[(j - 1) + (long)ld * ((sa - 1) + (long)ld * nf)] = cbk[(nd + a1 - 1) + (long)ldc * fl] * Gk[(j - 1) + (long)ld * ((sa - 1) + (long)ld * nf)];
  }));
}
// C(:,m) = G(:,S(m)); Qt(:,m) = e_S(m) - G(S(m),:)^T; slot(i,nf)=0 at i==S(m);  G is G_d(1,1,nf) (ld x ld), slot (ld,nfg)
void bb_fgather(const cd *G, int ld, int nrow, const int *S, int k, int nf, cd *C, cd *Qt, int *slot) {
  GUARD(Q->parallel_for(sycl::range<2>(k, nrow), [=](sycl::id<2> id) {
    long m = id[0], i = id[1];
    long s = S[m];
    C[i + (long)ld * m] = G[i + (long)ld * (s - 1)];
    cd q = -G[(s - 1) + (long)ld * i];
    if (i + 1 == s) { q = q + cd(1.0, 0.0); slot[i + (long)ld * (nf - 1)] = 0; }
    Qt[i + (long)ld * m] = q;
  }));
}
void bb_gemm(int ta, int tb, int m, int n, int k, double ar, double ai, const cd *A, int lda, const cd *B, int ldb, double br, double bi, cd *C, int ldc) {
  auto t = [](int x) { return x == 0 ? transpose::nontrans : (x == 1 ? transpose::trans : transpose::conjtrans); };
  GUARD(oneapi::mkl::blas::column_major::gemm(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, B, ldb, cd(br, bi), C, ldc));
}
}
