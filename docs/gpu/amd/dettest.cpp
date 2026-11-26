// Standalone check of k_det_phase (alf_hip_kernels.hip.cpp): synthetic LU diagonal whose raw product underflows.
// Reference: long double product of the unit-modulus factors a_i/|a_i|, times (-1)^{#pivots}, conjugated (cs(3) = conj(det)/|det|).
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <complex>
extern "C" int alfhk_det_phase(int n, int ld, void* A, void* ipiv, void* cs);
typedef std::complex<long double> ld_c;
static int run(const char* name, int n, double mag_lo, double mag_hi, int nlo, unsigned seed, bool pivots) {
  std::vector<double> A(2UL * n * n, 0.0);
  std::vector<int> piv(n);
  srand(seed);
  ld_c ref(1, 0); int cnt = 0;
  for (int i = 0; i < n; i++) {
    double th = 6.283185307179586 * rand() / RAND_MAX;
    double m = (i < nlo) ? mag_lo : mag_hi;      // first nlo entries tiny, the rest huge: raw running product underflows
    A[2UL * (i + (long)i * n)] = m * cos(th); A[2UL * (i + (long)i * n) + 1] = m * sin(th);
    piv[i] = (pivots && (rand() % 3 == 0)) ? ((i + 1) % n) + 1 : i + 1;
    if (piv[i] != i + 1) cnt++;
    ref *= ld_c(cosl(th), sinl(th));
  }
  if (cnt & 1) ref = -ref;
  ref = std::conj(ref);
  double *dA, *dcs; int* dp;
  hipMalloc(&dA, A.size() * 8); hipMalloc(&dp, n * 4); hipMalloc(&dcs, 3 * 16);
  hipMemcpy(dA, A.data(), A.size() * 8, hipMemcpyHostToDevice); hipMemcpy(dp, piv.data(), n * 4, hipMemcpyHostToDevice);
  hipMemset(dcs, 0, 48);
  alfhk_det_phase(n, n, dA, dp, dcs);
  hipDeviceSynchronize();
  double cs[6]; hipMemcpy(cs, dcs, 48, hipMemcpyDeviceToHost);
  std::complex<double> got(cs[4], cs[5]);
  double err = std::abs(std::complex<long double>(got.real(), got.imag()) - ref);
  bool ok = std::isfinite(err) && err < 1e-10 && std::fabs(std::abs(got) - 1.0) < 1e-10;
  printf("%-34s n=%5d  phase=(% .6f,% .6f) |phase|=%.6f  ref=(% .6f,% .6f)  err=%.2e  %s\n", name, n, got.real(), got.imag(), std::abs(got),
         (double)ref.real(), (double)ref.imag(), err, ok ? "OK" : "FAIL");
  hipFree(dA); hipFree(dp); hipFree(dcs);
  return ok ? 0 : 1;
}
int main() {
  int bad = 0;
  bad += run("unit magnitudes (control)", 3136, 1.0, 1.0, 0, 1, true);
  bad += run("mild spread", 3136, 0.5, 2.0, 1568, 2, true);
  bad += run("underflow in per-lane product", 3136, 1e-14, 1e14, 1568, 3, true);
  bad += run("underflow, no pivots, N=6272", 6272, 1e-7, 1e7, 3136, 4, false);
  printf("failures: %d\n", bad);
  return bad;
}
