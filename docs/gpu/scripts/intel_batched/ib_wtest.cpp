// Bitwise check: multi-walker kernels of alf_la_sycl_blk.cpp (alfsb_*_w) vs W single-walker calls on the same data.
// Build: icpx -fsycl -O2 -qmkl=sequential -I<repo>/Libraries/Modules ib_wtest.cpp -lmkl_sycl_blas
// IB_PLANT=1 perturbs one walker's input after the reference run: the check must then FAIL (positive control).
#include "alf_la_sycl_blk.cpp"
#include <cstring>
#include <vector>
static void fillr(std::vector<cd> &v, unsigned s) {
  for (size_t i = 0; i < v.size(); ++i) {
    s = s * 1103515245u + 12345u;
    double x = (s >> 8) / 16777216.0;
    s = s * 1103515245u + 12345u;
    double y = (s >> 8) / 16777216.0;
    v[i] = cd(x - 0.5, y - 0.5);
  }
}
int main() {
#ifdef IB_OMPINIT
  // oneAPI 2025.3 with the UR Level-Zero v1 adapter: a pure SYCL program fails at platform init (UR_RESULT_ERROR_UNINITIALIZED);
  // ALF initialises the device through libomptarget first. Do the same.
  {
    int x = 0;
#pragma omp target map(tofrom : x)
    x = 1;
  }
#endif
  sycl::queue q(sycl::gpu_selector_v, sycl::property::queue::in_order());
  if (alfsb_init(&q) != 0) return 2;
  const int W = 5, ld = 576, nd = 32, b = 32, nfl = 2, k = 24;
  const int64_t sG = (int64_t)ld * ld * 2, sgb = (int64_t)nd * nd * 2, sbs = (int64_t)nd * 2, ldc = 2 * nd + k * k,
                scb = ldc * 2, sC = (int64_t)ld * k, sS = k, ssl = (int64_t)ld * 2;
  std::vector<cd> hG(sG * W), hcb(scb * W), hC(sC * W);
  fillr(hG, 1); fillr(hcb, 2); fillr(hC, 3);
  std::vector<int> hbs(sbs * W), hS(sS * W), hkw(W);
  for (int w = 0; w < W; ++w) {
    for (int i = 0; i < sbs; ++i) hbs[w * sbs + i] = 1 + (int)((i * 7 + w * 13) % ld);
    for (int i = 0; i < k; ++i) hS[w * sS + i] = 1 + (int)((i * 11 + w * 5) % ld);
    hkw[w] = (w == 3) ? k - 5 : k;  // one walker with fewer pending columns
  }
  auto dz = [&](size_t n) { return sycl::malloc_device<cd>(n, q); };
  auto di = [&](size_t n) { return sycl::malloc_device<int>(n, q); };
  cd *G1 = dz(sG * W), *G2 = dz(sG * W), *gb1 = dz(sgb * W), *gb2 = dz(sgb * W), *cb = dz(scb * W), *C1 = dz(sC * W), *C2 = dz(sC * W),
     *Q1 = dz(sC * W), *Q2 = dz(sC * W), *R1 = dz(sG * W), *R2 = dz(sG * W);
  int *bs = di(sbs * W), *S = di(sS * W), *kw = di(W), *sl1 = di(ssl * W), *sl2 = di(ssl * W);
  q.memcpy(G1, hG.data(), 16 * sG * W); q.memcpy(G2, hG.data(), 16 * sG * W);
  q.memcpy(R1, hG.data(), 16 * sG * W); q.memcpy(R2, hG.data(), 16 * sG * W);
  q.memcpy(cb, hcb.data(), 16 * scb * W); q.memcpy(bs, hbs.data(), 4 * sbs * W); q.memcpy(S, hS.data(), 4 * sS * W);
  q.memcpy(kw, hkw.data(), 4 * W);
  q.memset(C1, 0, 16 * sC * W); q.memset(C2, 0, 16 * sC * W); q.memset(Q1, 0, 16 * sC * W); q.memset(Q2, 0, 16 * sC * W);
  q.memset(sl1, 1, 4 * ssl * W); q.memset(sl2, 1, 4 * ssl * W);
  q.wait();
  // reference: single-walker calls
  for (int w = 0; w < W; ++w) {
    alfsb_gather(nfl, 2, 1, ld, b, nd, gb1 + sgb * w, bs + sbs * w, G1 + sG * w);
    alfsb_scale(nfl, 2, 1, ld, b, nd, ldc, G1 + sG * w, bs + sbs * w, cb + scb * w);
    alfsb_fgather(G1 + sG * w, ld, ld, S + sS * w, hkw[w], 2, C1 + sC * w, Q1 + sC * w, sl1 + ssl * w);
    alfsb_gemm(0, 1, ld, ld, hkw[w], -1.0, 0.0, C1 + sC * w, ld, Q1 + sC * w, ld, 1.0, 0.0, R1 + sG * w, ld);
  }
  q.wait();
  if (getenv("IB_PLANT")) q.memcpy(G2 + sG * 2 + 1000, hG.data() + 7, 16).wait();
  alfsb_gather_w(W, nfl, 2, 1, ld, b, nd, gb2, sgb, bs, sbs, G2, sG);
  alfsb_scale_w(W, nfl, 2, 1, ld, b, nd, ldc, G2, sG, bs, sbs, cb, scb);
  alfsb_fgather_w(W, G2, sG, ld, ld, S, sS, kw, k, 2, C2, sC, Q2, sC, sl2, ssl);
  // gemm: walkers have different k here, so mode 0 per walker (the lockstep code pads to kmax; padding columns are zero)
  for (int w = 0; w < W; ++w)
    alfsb_gemm_w(0, 1, 0, 1, ld, ld, hkw[w], -1.0, 0.0, C2 + sC * w, ld, 0, Q2 + sC * w, ld, 0, 1.0, 0.0, R2 + sG * w, ld, 0);
  q.wait();
  auto cmpz = [&](const char *nm, cd *a, cd *c, size_t n) {
    std::vector<cd> x(n), y(n);
    q.memcpy(x.data(), a, 16 * n); q.memcpy(y.data(), c, 16 * n); q.wait();
    size_t bad = 0;
    for (size_t i = 0; i < n; ++i) bad += std::memcmp(&x[i], &y[i], 16) != 0;
    std::printf("%-8s n=%zu differing=%zu %s\n", nm, n, bad, bad ? "FAIL" : "BITWISE");
    return bad == 0;
  };
  auto cmpi = [&](const char *nm, int *a, int *c, size_t n) {
    std::vector<int> x(n), y(n);
    q.memcpy(x.data(), a, 4 * n); q.memcpy(y.data(), c, 4 * n); q.wait();
    size_t bad = 0;
    for (size_t i = 0; i < n; ++i) bad += x[i] != y[i];
    std::printf("%-8s n=%zu differing=%zu %s\n", nm, n, bad, bad ? "FAIL" : "BITWISE");
    return bad == 0;
  };
  bool ok = cmpz("gather", gb1, gb2, sgb * W) & cmpz("scale", G1, G2, sG * W) & cmpz("C", C1, C2, sC * W) & cmpz("Qt", Q1, Q2, sC * W) &
            cmpi("slot", sl1, sl2, ssl * W) & cmpz("gemm", R1, R2, sG * W);
  // gemm_batch (mode 1, equal k for all walkers) vs W single gemms: report the max difference (expected nonzero, Phase 6)
  q.memcpy(R1, hG.data(), 16 * sG * W); q.memcpy(R2, hG.data(), 16 * sG * W); q.wait();
  alfsb_gemm_w(0, W, 0, 1, ld, ld, k, -1.0, 0.0, C1, ld, sC, Q1, ld, sC, 1.0, 0.0, R1, ld, sG);
  alfsb_gemm_w(1, W, 0, 1, ld, ld, k, -1.0, 0.0, C1, ld, sC, Q1, ld, sC, 1.0, 0.0, R2, ld, sG);
  q.wait();
  {
    std::vector<cd> x(sG * W), y(sG * W);
    q.memcpy(x.data(), R1, 16 * sG * W); q.memcpy(y.data(), R2, 16 * sG * W); q.wait();
    double md = 0; size_t nd2 = 0;
    for (size_t i = 0; i < x.size(); ++i) { md = std::max(md, std::abs(x[i] - y[i])); nd2 += std::memcmp(&x[i], &y[i], 16) != 0; }
    std::printf("gemm_batch(mode1) vs loop(mode0): differing=%zu maxabs=%.3e (information only)\n", nd2, md);
  }
  std::printf("RESULT %s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
