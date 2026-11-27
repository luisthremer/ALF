// ALF GPU=amd: small device kernels written in HIP instead of OpenMP target regions.
// A libgomp-GCN target region with a parallel loop costs ~150 us to launch; a HIP kernel launch + hipDeviceSynchronize ~15 us.
// Every function takes device pointers (from use_device_addr in alf_hipk_mod / linalg_hip_mod) and ends in fin() = alf_hip_fin(): the work is synchronised
// before it returns, or, with the deferred synchronisation (below), only marked as pending and synchronised by alf_hip_presync() before
// the next OpenMP kernel launch or device copy, because the OpenMP target regions run on their own queues. Index conventions follow
// the Fortran arrays (column major).
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>

struct cd {
  double r, i;
  __device__ cd() : r(0), i(0) {}
  __device__ cd(double a, double b) : r(a), i(b) {}
};
__device__ inline cd operator+(cd a, cd b) { return cd(a.r + b.r, a.i + b.i); }
__device__ inline cd operator-(cd a, cd b) { return cd(a.r - b.r, a.i - b.i); }
__device__ inline cd operator-(cd a) { return cd(-a.r, -a.i); }
__device__ inline cd operator*(cd a, cd b) { return cd(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
__device__ inline cd operator*(double s, cd a) { return cd(s * a.r, s * a.i); }
__device__ inline cd operator*(cd a, double s) { return cd(s * a.r, s * a.i); }
__device__ inline cd operator/(cd a, double s) { return cd(a.r / s, a.i / s); }
__device__ inline cd& operator+=(cd& a, cd b) { a.r += b.r; a.i += b.i; return a; }

// ---- deferred synchronisation -------------------------------------------------------------------------------------
// All HIP kernels, hipBLAS calls and the QR run on stream 0 and are therefore ordered among themselves. Only the OpenMP
// target regions (libgomp queue) and the copies of the libgomp plugin need the HIP work to be finished: the patched GCN
// plugin calls alf_hip_presync() before every OpenMP kernel launch and before every device copy/free
// (alf_gcn_set_presync, registered in alfhk_register_presync). So the wrappers only mark the work as pending.
// ALF_HIP_ASYNC=0 restores a synchronisation after every call (debugging, A/B timing).
extern "C" {
int alf_hip_pending = 0;
static int g_async = -1;
static int g_registered = 0;
int alf_hip_presync() {
  if (!alf_hip_pending) return 0;
  alf_hip_pending = 0;
  hipError_t e = hipDeviceSynchronize();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: sync error %s\n", hipGetErrorString(e)); return 2; }
  return 0;
}
int alfhk_register_presync();
int alf_hip_fin() {
  hipError_t e = hipGetLastError();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: launch error %s\n", hipGetErrorString(e)); return 1; }
  // the mode is fixed at the first call: register the plugin hook here too (HIP kernels may run before linalg_hip_mod's init)
  if (g_async < 0 && !g_registered) alfhk_register_presync();
  if (g_async < 0) { const char* v = getenv("ALF_HIP_ASYNC"); g_async = (v && v[0] == '0') ? 0 : g_registered;   // without the plugin hook every call must synchronise
    printf(" alf_la: HIP deferred synchronisation: %s\n", g_async ? "ON" : (g_registered ? "OFF (ALF_HIP_ASYNC=0)" : "OFF (libgomp plugin has no presync hook)")); }
  if (g_async) { alf_hip_pending = 1; return 0; }
  e = hipDeviceSynchronize();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: sync error %s\n", hipGetErrorString(e)); return 2; }
  return 0;
}
// hook the presync into the (patched) libgomp GCN plugin; returns 1 if the plugin has the hook, else 0 (then every call syncs)
int alfhk_register_presync() {
  void* h = dlopen("libgomp-plugin-gcn.so.1", RTLD_NOLOAD | RTLD_LAZY);
  if (!h) return 0;
  void (*set)(int (*)()) = (void (*)(int (*)()))dlsym(h, "alf_gcn_set_presync");
  if (!set) return 0;
  set(alf_hip_presync);
  g_registered = 1;
  return 1;
}
}
static int fin() { return alf_hip_fin(); }
static inline unsigned nb(long n, int t) { return (unsigned)((n + t - 1) / t); }
#define LAUNCH(k, grid, blk, ...) hipLaunchKernelGGL(k, dim3(grid), dim3(blk), 0, 0, __VA_ARGS__)


// ---- linalg_hip_mod helpers: copies, identity, column pre-pivoting of la_zqrpp --------------------------------
__global__ void k_copyz(long n, const cd* a, cd* b) {
  long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) b[i] = a[i];
}
__global__ void k_ident(int n, int ld, cd* U) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; U[i + (long)j * ld] = (i == j) ? cd(1, 0) : cd(0, 0);
}
// nrm(j) = 2-norm of column j (scaled); NaN if the column holds a NaN or Inf; one block per column
__global__ void k_colnorms(int m, int n, const cd* A, int lda, double* nrm) {
  __shared__ double sm[256]; __shared__ int sf[256];
  int j = blockIdx.x; const cd* a = A + (long)j * lda;
  double mx = 0; int nf = 0;
  for (int i = threadIdx.x; i < m; i += 256) {
    double re = a[i].r, im = a[i].i;
    if (((__double_as_longlong(re) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL) || ((__double_as_longlong(im) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL)) nf++;
    mx = fmax(mx, fmax(fabs(re), fabs(im)));
  }
  sm[threadIdx.x] = mx; sf[threadIdx.x] = nf; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) { sm[threadIdx.x] = fmax(sm[threadIdx.x], sm[threadIdx.x + o]); sf[threadIdx.x] += sf[threadIdx.x + o]; } __syncthreads(); }
  mx = sm[0]; nf = sf[0]; __syncthreads();
  double s = 0;
  if (nf == 0 && mx > 0.0) for (int i = threadIdx.x; i < m; i += 256) { double v = a[i].r / mx; s += v * v; v = a[i].i / mx; s += v * v; }
  sm[threadIdx.x] = s; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) sm[threadIdx.x] += sm[threadIdx.x + o]; __syncthreads(); }
  if (threadIdx.x == 0) nrm[j] = (nf > 0) ? __longlong_as_double(0x7FF8000000000000LL) : (mx > 0.0 ? mx * sqrt(sm[0]) : 0.0);
}
__global__ void k_gather_cols(int m, int n, const cd* A, int lda, const int* ipvt, cd* scr) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)m * n) return;
  int i = t % m, j = t / m; scr[t] = A[i + (long)(ipvt[j] - 1) * lda];
}
__global__ void k_scatter_cols(int m, int n, const cd* scr, cd* A, int lda) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)m * n) return;
  int i = t % m, j = t / m; A[i + (long)j * lda] = scr[t];
}

extern "C" {
int alfhk_copyz(long n, void* a, void* b) {
  if (n <= 0) return 0;
  LAUNCH(k_copyz, nb(n, 256), 256, n, (const cd*)a, (cd*)b);
  return fin();
}
int alfhk_ident(int n, int ld, void* U) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_ident, nb(t,256), 256, n, ld, (cd*)U); return fin(); }
int alfhk_colnorms(int m, int n, void* A, int lda, void* nrm) { if (n<=0) return 0; LAUNCH(k_colnorms, n, 256, m, n, (const cd*)A, lda, (double*)nrm); return fin(); }
static cd* gscr = nullptr; static long gscr_n = 0;
int alfhk_permute_cols(int m, int n, void* A, int lda, void* ipvt) {
  long t=(long)m*n; if (t<=0) return 0;
  if (gscr_n < t) { if (gscr) (void)hipFree(gscr); gscr = nullptr; gscr_n = 0; if (hipMalloc((void**)&gscr, t * sizeof(cd)) != hipSuccess) return 3; gscr_n = t; }
  LAUNCH(k_gather_cols, nb(t,256), 256, m, n, (const cd*)A, lda, (const int*)ipvt, gscr);
  LAUNCH(k_scatter_cols, nb(t,256), 256, m, n, (const cd*)gscr, (cd*)A, lda);
  return fin();
}
}  // extern "C"

// ---- gpu_stab_mod kernels (P2.4 step 2 amd, commit 4) ----------------------------------------------------------------
// Each one reproduces the OpenMP kernel of the same name in gpu_stab_mod (CM contracts: S1 cs(6), LU phase cs(4/5),
// tracked CGR phase cs(7) -> cs(3), NaN-or-Inf tests).  Index conventions: Fortran column major, cs(k) = cs[k-1].
__device__ inline cd operator/(cd a, cd b) {   // Smith's algorithm (gfortran's complex division without -ffast-math)
  if (fabs(b.r) >= fabs(b.i)) { double q = b.i / b.r, d = b.r + b.i * q; return cd((a.r + a.i * q) / d, (a.i - a.r * q) / d); }
  double q = b.r / b.i, d = b.r * q + b.i;
  return cd((a.r * q + a.i) / d, (a.i * q - a.r) / d);
}
__device__ inline cd conjg(cd a) { return cd(a.r, -a.i); }
__device__ inline double cabsd(cd a) { return hypot(a.r, a.i); }
__device__ inline bool nonfinite(double x) { return (__double_as_longlong(x) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL; }
// complex product over the 64 lanes of one wave (butterfly; every lane holds the total afterwards)
__device__ inline cd wave_prod(cd x) {
  for (int o = 32; o > 0; o >>= 1) { cd y(__shfl_xor(x.r, o, 64), __shfl_xor(x.i, o, 64)); x = x * y; }
  return x;
}
__device__ inline int wave_sum_i(int c) { for (int o = 32; o > 0; o >>= 1) c += __shfl_xor(c, o, 64); return c; }
// det of the ZGEQRF reflector with scalar z: -(u*u), u = z/|z|; 1 for z = 0 (refl_det of gpu_stab_mod)
__device__ inline cd refl_det(cd z) {
  if (z.r == 0.0 && z.i == 0.0) return cd(1.0, 0.0);
  cd u = z / cabsd(z);
  return -(u * u);
}
#define IJ long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return; int i = t % n, j = t / n;
__global__ void k_setd(int n, cd* D) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) D[i] = cd(1, 0); }
__global__ void k_setl(int n, double* L) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) L[i] = 0.0; }
__global__ void k_scale_cols(int n, int ld, cd* U, const cd* D) { IJ U[i + (long)j * ld] = U[i + (long)j * ld] * D[j]; }
__global__ void k_diagabs(int n, int ld, const cd* U, cd* D) {
  int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) D[i] = cd(cabsd(U[i + (long)i * ld]), 0.0);
}
__global__ void k_divrows(int n, int ld, cd* U, const cd* D) { IJ if (j >= i) U[i + (long)j * ld] = U[i + (long)j * ld] / D[i].r; }
// cs(1) = Phase = prod U(i,i) (* sgn, conjugated if cj), cs(2) = 1/Phase, cs(6) = prod refl_det(tau(i)) * Phase
__global__ void k_phase(int n, int ld, const cd* U, int sgn, int cj, const cd* tau, cd* cs) {
  cd ph(1.0, 0.0), dq(1.0, 0.0);
  for (int i = threadIdx.x; i < n; i += 64) { ph = ph * U[i + (long)i * ld]; dq = dq * refl_det(tau[i]); }
  ph = wave_prod(ph); dq = wave_prod(dq);
  if (threadIdx.x == 0) {
    if (sgn < 0) ph = -ph;
    if (cj) ph = conjg(ph);
    cs[0] = ph; cs[1] = cd(1.0, 0.0) / ph; cs[5] = dq * ph;
  }
}
__global__ void k_scale_row1(int n, int ld, cd* U, const cd* cs) { int j = blockIdx.x * blockDim.x + threadIdx.x; if (j < n) U[(long)j * ld] = cs[1] * U[(long)j * ld]; }
__global__ void k_scale_col1(int n, int ld, cd* U, const cd* cs) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) U[i] = cs[0] * U[i]; }
__global__ void k_triu(int n, int ld, const cd* U, int ldr, cd* R) { IJ R[i + (long)j * ldr] = (i <= j) ? U[i + (long)j * ld] : cd(0, 0); }
// mode 0: W(i,j) = X(perm(i),j)  1: W(i,j) = X(i,perm(j))  2: W(perm(i),j) = X(i,j)  3: W(i,perm(j)) = X(i,j)
__global__ void k_perm(int mode, int n, const cd* X, int ldx, const int* perm, cd* W, int ldw) {
  IJ
  if (mode == 0) W[i + (long)j * ldw] = X[(perm[i] - 1) + (long)j * ldx];
  else if (mode == 1) W[i + (long)j * ldw] = X[i + (long)(perm[j] - 1) * ldx];
  else if (mode == 2) W[(perm[i] - 1) + (long)j * ldw] = X[i + (long)j * ldx];
  else W[i + (long)(perm[j] - 1) * ldw] = X[i + (long)j * ldx];
}
__global__ void k_perm_vec(int n, const double* L, const int* perm, double* W) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) W[i] = L[perm[i] - 1]; }
__global__ void k_ct(int n, const cd* A, int lda, cd* B, int ldb) { IJ B[i + (long)j * ldb] = conjg(A[j + (long)i * lda]); }
__global__ void k_rowdiv(int n, int lda, cd* A, const cd* D) { IJ A[i + (long)j * lda] = A[i + (long)j * lda] / D[i]; }
__global__ void k_colinv(int n, int lda, cd* A, const cd* D) { IJ double sv = 1.0 / D[j].r; A[i + (long)j * lda] = A[i + (long)j * lda] * sv; }
__global__ void k_tpup(int n, cd* T, const cd* R, const cd* Dr, const cd* Dl) {
  IJ cd v = (Dr[i] * T[i + (long)j * n]) * Dl[j]; T[i + (long)j * n] = v + R[i + (long)j * n];
}
// LU phase: cs(4) = prod A(i,i)/|A(i,i)| with the pivot sign (0 on a zero pivot), cs(5) = raw prod A(i,i) with the pivot sign
__global__ void k_det_phase(int n, int lda, const cd* A, const int* ipiv, cd* cs) {
  cd raw(1.0, 0.0), p(1.0, 0.0); int sw = 0, zero = 0;
  for (int i = threadIdx.x; i < n; i += 64) {
    cd a = A[i + (long)i * lda]; double am = cabsd(a);
    raw = raw * a;
    if (am == 0.0) zero = 1; else p = p * (a / am);
    if (ipiv[i] != i + 1) sw++;
  }
  raw = wave_prod(raw); p = wave_prod(p); sw = wave_sum_i(sw); zero = wave_sum_i(zero);
  if (threadIdx.x == 0) {
    if (sw & 1) { raw = -raw; p = -p; }
    cs[3] = zero ? cd(0.0, 0.0) : p; cs[4] = raw;
  }
}
// CGR phase: cs(3) = ph * sgn * prod (R_ii phase, conjugated for nvar = 2) * reflector factors; ph = cs(7) (dev = 1) or ph0
__global__ void k_cgr_phase(int n, int ld, const cd* T, const cd* tau, cd* cs, int nvar, int sgn, int bad, double p0r, double p0i, int dev) {
  cd acc(1.0, 0.0);
  for (int i = threadIdx.x; i < n; i += 64) {
    cd z = tau[i], tii = T[i + (long)i * ld]; double at = cabsd(tii); cd f;
    if (nvar == 1) f = tii / at; else { f = conjg(tii) / at; z = conjg(z); }
    if (z.r != 0.0 || z.i != 0.0) { double x = cabsd(z); cd w = cd(1.0, 0.0) - (2.0 * (z / x)) * (z.r / x); f = f * (w / cabsd(w)); }
    acc = acc * f;
  }
  acc = wave_prod(acc);
  if (threadIdx.x == 0) {
    cd ph = dev ? cs[6] : cd(p0r, p0i);
    if (sgn < 0) ph = -ph;
    ph = ph * acc;
    if (bad == 1) ph = -ph;
    cs[2] = ph;
  }
}
// scale test: r(1) = number of NaN/Inf scales (+1 if bad), r(2) = Re D(1), r(3) = Re D(n)
__global__ void k_testscale(int n, const cd* D, double* r, int bad) {
  __shared__ int s[256];
  int c = 0;
  for (int i = threadIdx.x; i < n; i += 256) if (nonfinite(D[i].r) || nonfinite(D[i].i)) c++;
  s[threadIdx.x] = c; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) s[threadIdx.x] += s[threadIdx.x + o]; __syncthreads(); }
  if (threadIdx.x == 0) { r[0] = (double)(s[0] + bad); r[1] = D[0].r; r[2] = D[n - 1].r; }
}
// precision check: per column (serial in i) sum and max of |A-B| over finite pairs, NaN/Inf counts; then the column sums in order
__global__ void k_prec1(int n, const cd* A, const cd* B, double* pw, double* pm, double* pa, double* pb) {
  int j = blockIdx.x * blockDim.x + threadIdx.x; if (j >= n) return;
  double cs = 0, xm = 0, na = 0, nb = 0;
  for (int i = 0; i < n; i++) {
    cd a = A[i + (long)j * n], b = B[i + (long)j * n];
    bool fa = nonfinite(a.r) || nonfinite(a.i), fb = nonfinite(b.r) || nonfinite(b.i);
    if (fa) na += 1;
    if (fb) nb += 1;
    if (!(fa || fb)) { double d = cabsd(a - b); cs += d; xm = fmax(xm, d); }
  }
  pw[j] = cs; pm[j] = xm; pa[j] = na; pb[j] = nb;
}
__global__ void k_prec2(int n, const double* pw, const double* pm, const double* pa, const double* pb, double* r) {
  double xs = 0, xm = 0, na = 0, nb = 0;
  for (int t = 0; t < n; t++) { xs += pw[t]; xm = fmax(xm, pm[t]); na += pa[t]; nb += pb[t]; }
  r[3] = xm; r[4] = xs; r[5] = na; r[6] = nb;
}

extern "C" {
int alfhk_setd(int n, void* D) { if (n<=0) return 0; LAUNCH(k_setd, nb(n,256), 256, n, (cd*)D); return fin(); }
int alfhk_setl(int n, void* L) { if (n<=0) return 0; LAUNCH(k_setl, nb(n,256), 256, n, (double*)L); return fin(); }
int alfhk_scale_cols(int n, int ld, void* U, void* D) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_scale_cols, nb(t,256), 256, n, ld, (cd*)U, (const cd*)D); return fin(); }
int alfhk_diagabs(int n, int ld, void* U, void* D) { if (n<=0) return 0; LAUNCH(k_diagabs, nb(n,256), 256, n, ld, (const cd*)U, (cd*)D); return fin(); }
int alfhk_divrows(int n, int ld, void* U, void* D) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_divrows, nb(t,256), 256, n, ld, (cd*)U, (const cd*)D); return fin(); }
int alfhk_phase(int n, int ld, void* U, int sgn, int cj, void* tau, void* cs) { LAUNCH(k_phase, 1, 64, n, ld, (const cd*)U, sgn, cj, (const cd*)tau, (cd*)cs); return fin(); }
int alfhk_scale_row1(int n, int ld, void* U, void* cs) { if (n<=0) return 0; LAUNCH(k_scale_row1, nb(n,256), 256, n, ld, (cd*)U, (const cd*)cs); return fin(); }
int alfhk_scale_col1(int n, int ld, void* U, void* cs) { if (n<=0) return 0; LAUNCH(k_scale_col1, nb(n,256), 256, n, ld, (cd*)U, (const cd*)cs); return fin(); }
int alfhk_triu(int n, int ld, void* U, int ldr, void* R) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_triu, nb(t,256), 256, n, ld, (const cd*)U, ldr, (cd*)R); return fin(); }
int alfhk_perm(int mode, int n, void* X, int ldx, void* perm, void* W, int ldw) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_perm, nb(t,256), 256, mode, n, (const cd*)X, ldx, (const int*)perm, (cd*)W, ldw); return fin(); }
int alfhk_perm_vec(int n, void* L, void* perm, void* W) { if (n<=0) return 0; LAUNCH(k_perm_vec, nb(n,256), 256, n, (const double*)L, (const int*)perm, (double*)W); return fin(); }
int alfhk_ct(int n, void* A, int lda, void* B, int ldb) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_ct, nb(t,256), 256, n, (const cd*)A, lda, (cd*)B, ldb); return fin(); }
int alfhk_rowdiv(int n, int lda, void* A, void* D) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_rowdiv, nb(t,256), 256, n, lda, (cd*)A, (const cd*)D); return fin(); }
int alfhk_colinv(int n, int lda, void* A, void* D) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_colinv, nb(t,256), 256, n, lda, (cd*)A, (const cd*)D); return fin(); }
int alfhk_tpup(int n, void* T, void* R, void* Dr, void* Dl) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_tpup, nb(t,256), 256, n, (cd*)T, (const cd*)R, (const cd*)Dr, (const cd*)Dl); return fin(); }
int alfhk_det_phase(int n, int ld, void* A, void* ipiv, void* cs) { LAUNCH(k_det_phase, 1, 64, n, ld, (const cd*)A, (const int*)ipiv, (cd*)cs); return fin(); }
int alfhk_cgr_phase(int n, int ld, void* T, void* tau, void* cs, int nvar, int sgn, int bad, double p0r, double p0i, int dev) {
  LAUNCH(k_cgr_phase, 1, 64, n, ld, (const cd*)T, (const cd*)tau, (cd*)cs, nvar, sgn, bad, p0r, p0i, dev); return fin(); }
int alfhk_testscale(int n, void* D, void* r, int bad) { if (n<=0) return 0; LAUNCH(k_testscale, 1, 256, n, (const cd*)D, (double*)r, bad); return fin(); }
static double* pscr = nullptr; static long pscr_n = 0;
int alfhk_prec(int n, void* A, void* B, void* r) {
  if (pscr_n < n) { if (pscr) (void)hipFree(pscr); pscr = nullptr; pscr_n = 0; if (hipMalloc((void**)&pscr, 4L * n * sizeof(double)) != hipSuccess) return 3; pscr_n = n; }
  double *pw = pscr, *pm = pscr + n, *pa = pscr + 2L * n, *pb = pscr + 3L * n;
  LAUNCH(k_prec1, nb(n, 64), 64, n, (const cd*)A, (const cd*)B, pw, pm, pa, pb);
  LAUNCH(k_prec2, 1, 1, n, (const double*)pw, (const double*)pm, (const double*)pa, (const double*)pb, (double*)r);
  return fin();
}
}  // extern "C"

// ---- gpu_wrap_mod and gpu_lowrank_mod kernels (P2.4 step 2 amd, commit 5) ----------------------------------------------
// Same loop bodies as the OpenMP kernels of the same name (one thread per loop iteration of the collapsed loops).
// M(P(l), j) = W(l, j)
__global__ void k_scatter_rows(int nd, int ncol, const int* P, cd* M, int ld, const cd* W, int ldw) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)nd * ncol) return;
  int l = t % nd, j = t / nd;
  M[(P[l] - 1) + (long)j * ld] = W[l + (long)j * ldw];
}
// k_vslice_R: M(r,j) = (product over the active Op_V factors of row r, in order (desc: reversed)) * M(r,j)
__global__ void k_vslice_r(int nd, int ncol, int ld, int nf, int t, int cj, int desc, int ldptr, int nop, int ldn,
                           const int* vr_ptr, const int* vr_ent, const int* V_meta, const int* V_type, const cd* nsig,
                           const cd* V_exp, cd* M) {
  long tt = (long)blockIdx.x * blockDim.x + threadIdx.x; if (tt >= (long)nd * ncol) return;
  int r = tt % nd, j = tt / nd;
  int pa = vr_ptr[r + (long)(nf - 1) * ldptr], pb = vr_ptr[r + 1 + (long)(nf - 1) * ldptr] - 1;
  if (pb < pa) return;
  cd x = M[r + (long)j * ld];
  for (int q = 0; q <= pb - pa; q++) {
    int p = desc ? pb - q : pa + q;
    int n = vr_ent[2 * (p - 1)], i = vr_ent[2 * (p - 1) + 1];
    const int* mt = V_meta + 4L * ((n - 1) + (long)(nf - 1) * nop);
    if (mt[3] != 1) {
      int sp = (int)rint(nsig[(n - 1) + (long)(t - 1) * ldn].r);     // nint (round half away from zero; sp is integral)
      cd z = V_exp[mt[2] + i + (long)(sp + V_type[n - 1]) * mt[0] - 1];
      if (cj) z.i = -z.i;
      x = z * x;
    }
  }
  M[r + (long)j * ld] = x;
}
// k_blk_gather: gbk(a,c,fl) = Gk(bsk(a,fl), bsk(c,fl), nf(fl)),  gbk(nd,nd,2), bsk(nd,2), Gk(ld,ld,*)
__global__ void k_blk_gather(int nfl, int nfa, int nfb, int ld, int b, int nd, const cd* G, const int* bsk, cd* gbk) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)nfl * b * b) return;
  int a = t % b, c = (t / b) % b, fl = (int)(t / ((long)b * b)) + 1;
  int nf = (fl == 2) ? nfb : nfa;
  const int* bs = bsk + (long)(fl - 1) * nd;
  gbk[a + (long)c * nd + (long)(fl - 1) * nd * nd] = G[(bs[a] - 1) + (long)(bs[c] - 1) * ld + (long)(nf - 1) * ld * ld];
}
// k_blk_scale: row factor cbk(a) on row bsk(a) (and column factor cbk(nd+c) if the column is block site c); column factor on
// the columns bsk(a) for the rows outside the block.  One thread per (fl, a, j), as the OpenMP collapse(3) loop.
__global__ void k_blk_scale(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd* G, const int* bsk, const cd* cbk) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)nfl * b * ld) return;
  int j = t % ld, a = (t / ld) % b, fl = (int)(t / ((long)ld * b)) + 1;
  int nf = (fl == 2) ? nfb : nfa;
  const int* bs = bsk + (long)(fl - 1) * nd;
  const cd* cb = cbk + (long)(fl - 1) * ldc;
  cd* Gf = G + (long)(nf - 1) * ld * ld;
  int sa = bs[a] - 1, cj = 0;
  for (int c = 0; c < b; c++) if (bs[c] - 1 == j) cj = c + 1;
  cd v = cb[a] * Gf[sa + (long)j * ld];
  if (cj > 0) v = cb[nd + cj - 1] * v;
  Gf[sa + (long)j * ld] = v;
  if (cj == 0) Gf[j + (long)sa * ld] = cb[nd + a] * Gf[j + (long)sa * ld];
}
// k_sub_flush_gather: C(i,m) = G(i,S(m)); Qt(i,m) = -G(S(m),i) (+1 and slot(i,nf) = 0 on i = S(m))
__global__ void k_sub_flush_gather(int ld, int nrow, int k, int nf, int ldslot, const cd* G, const int* S, cd* C, cd* Qt, int* slot) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)k * nrow) return;
  int i = t % nrow, m = t / nrow, s = S[m] - 1;
  C[i + (long)m * ld] = G[i + (long)s * ld];
  cd q = -G[s + (long)i * ld];
  if (i == s) { q.r += 1.0; slot[i + (long)(nf - 1) * ldslot] = 0; }
  Qt[i + (long)m * ld] = q;
}
extern "C" {
int alfhk_scatter_rows(int nd, int ncol, void* P, void* M, int ld, void* W, int ldw) {
  long tot = (long)nd * ncol; if (tot <= 0) return 0;
  LAUNCH(k_scatter_rows, nb(tot, 256), 256, nd, ncol, (const int*)P, (cd*)M, ld, (const cd*)W, ldw); return fin(); }
int alfhk_vslice_r(int nd, int ncol, int ld, int nf, int t, int cj, int desc, int ldptr, int nop, int ldn, void* vr_ptr, void* vr_ent,
                   void* V_meta, void* V_type, void* nsig, void* V_exp, void* M) {
  long tot = (long)nd * ncol; if (tot <= 0) return 0;
  LAUNCH(k_vslice_r, nb(tot, 256), 256, nd, ncol, ld, nf, t, cj, desc, ldptr, nop, ldn, (const int*)vr_ptr, (const int*)vr_ent,
         (const int*)V_meta, (const int*)V_type, (const cd*)nsig, (const cd*)V_exp, (cd*)M); return fin(); }
int alfhk_blk_gather(int nfl, int nfa, int nfb, int ld, int b, int nd, void* G, void* bsk, void* gbk) {
  long tot = (long)nfl * b * b; if (tot <= 0) return 0;
  LAUNCH(k_blk_gather, nb(tot, 128), 128, nfl, nfa, nfb, ld, b, nd, (const cd*)G, (const int*)bsk, (cd*)gbk); return fin(); }
int alfhk_blk_scale(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, void* G, void* bsk, void* cbk) {
  long tot = (long)nfl * b * ld; if (tot <= 0) return 0;
  LAUNCH(k_blk_scale, nb(tot, 256), 256, nfl, nfa, nfb, ld, b, nd, ldc, (cd*)G, (const int*)bsk, (const cd*)cbk); return fin(); }
int alfhk_sub_flush_gather(int ld, int nrow, int k, int nf, int ldslot, void* G, void* S, void* C, void* Qt, void* slot) {
  long tot = (long)k * nrow; if (tot <= 0) return 0;
  LAUNCH(k_sub_flush_gather, nb(tot, 256), 256, ld, nrow, k, nf, ldslot, (const cd*)G, (const int*)S, (cd*)C, (cd*)Qt, (int*)slot);
  return fin(); }
}  // extern "C"

// ---- gpu_obs_mod kernels (P2.4 step 2 amd, commit 6): equal-time correlations (k_corr) and the scalars Zkin, Zpot, Zrho of
// measure_on; one block per output element, block tree reduction (summation order differs from the OpenMP reductions).
// From gpu-port; the H7 weights come in through uo(no) (host array o_uo).
#define RB 128
// One block per output element m = (imj,no,no2) (gather over the lattice sites I) for m < MT; blocks MT..MT+Norb-1 do the acc0 sums (kinds 3, 4).
__global__ void k_corr(int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld,
                       const cd* G, const int* jc, const int* inv, cd* acc, cd* acc0) {
  __shared__ double red[6][RB];
  long MT = (long)N * Norb * Norb;
  long m = blockIdx.x;
  bool zero_part = (m >= MT);
  int imj, no, no2;
  if (!zero_part) { imj = m % N; no = (m / N) % Norb; no2 = m / ((long)N * Norb); }
  else { no = (int)(m - MT); imj = 0; no2 = 0; }
  cd s1, s2, s3;
  const long ld2 = (long)ld * ld;
  for (int I = threadIdx.x; I < N; I += blockDim.x) {
    if (!zero_part) {
      int J = jc[I + (long)imj * N] - 1;
      int i1 = inv[I + (long)no * N] - 1, j1 = inv[J + (long)no2 * N] - 1;
      cd one(1.0, 0.0);
      if (kind == 1) {
        cd z;
        for (int nf = 0; nf < NFL; nf++) { cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one; z += a; }
        s1 += z * ns;
      } else if (kind == 2) {
        cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
        s1 += (a * G[i1 + (long)j1 * ld]) * ns;
      } else if (kind == 3) {
        cd zi, zj, z;
        for (int nf = 0; nf < NFL; nf++) {
          zi += one - G[i1 + (long)i1 * ld + nf * ld2];
          zj += one - G[j1 + (long)j1 * ld + nf * ld2];
          cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one;
          z += a * G[i1 + (long)j1 * ld + nf * ld2];
        }
        s1 += ((zi * ns) * (zj * ns)) + z * ns;
      } else {
        cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
        cd b = -G[j1 + (long)i1 * ld + ld2]; if (i1 == j1) b = b + one;
        cd zxy = a * G[i1 + (long)j1 * ld + ld2] + b * G[i1 + (long)j1 * ld];
        cd zi = (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
        cd zj = (one - G[j1 + (long)j1 * ld + ld2]) - (one - G[j1 + (long)j1 * ld]);
        cd zz = a * G[i1 + (long)j1 * ld] + b * G[i1 + (long)j1 * ld + ld2] + zi * zj;
        s1 += zz; s2 += zxy; s3 += (2.0 * zxy + zz) / 3.0;
      }
    } else {
      int i1 = inv[I + (long)no * N] - 1;
      cd one(1.0, 0.0);
      if (kind == 3) {
        cd zi;
        for (int nf = 0; nf < NFL; nf++) zi += one - G[i1 + (long)i1 * ld + nf * ld2];
        s1 += zi * ns;
      } else {
        s1 += (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
      }
    }
  }
  red[0][threadIdx.x] = s1.r; red[1][threadIdx.x] = s1.i;
  red[2][threadIdx.x] = s2.r; red[3][threadIdx.x] = s2.i;
  red[4][threadIdx.x] = s3.r; red[5][threadIdx.x] = s3.i;
  __syncthreads();
  for (int st = RB / 2; st > 0; st >>= 1) {
    if (threadIdx.x < st) for (int q = 0; q < 6; q++) red[q][threadIdx.x] += red[q][threadIdx.x + st];
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    cd w(wr, wi);
    cd r1(red[0][0], red[1][0]), r2(red[2][0], red[3][0]), r3(red[4][0], red[5][0]);
    if (!zero_part) {
      acc[m + (long)(k - 1) * MT] += r1 * w;
      if (kind == 4) { acc[m + (long)k * MT] += r2 * w; acc[m + (long)(k + 1) * MT] += r3 * w; }
    } else {
      acc0[no + (long)(k - 1) * Norb] += r1 * w;
      if (kind == 4) acc0[no + (long)(k + 1) * Norb] += (r1 * w) / sq3;
    }
  }
}

// Zkin (nk terms), Zpot, Zrho in one block; out[0..5] = re/im of Zkin(before N_SUN), Zpot, Zrho
__global__ void k_scalars(int nk, const int* hi, const int* hj, const int* hf, const cd* hc, int N, int Norb, int NFL, int dec,
                          int ld, const cd* G, const int* inv, const double* uo, double* out) {
  __shared__ double red[6][RB];
  const long ld2 = (long)ld * ld;
  cd zk, zp, zr;
  cd one(1.0, 0.0);
  for (int k = threadIdx.x; k < nk; k += blockDim.x) {
    int i1 = hi[k] - 1, j1 = hj[k] - 1, nf = hf[k] - 1;
    cd zc = -G[j1 + (long)i1 * ld + nf * ld2];
    if (i1 == j1) zc = zc + one;
    zk += hc[k] * zc;
  }
  for (int t = threadIdx.x; t < N * Norb; t += blockDim.x) {
    int I = t % N, no = t / N;
    int i1 = inv[I + (long)no * N] - 1;
    cd zc = one - G[i1 + (long)i1 * ld];
    cd zt = one - G[i1 + (long)i1 * ld + (long)(dec - 1) * ld2];
    zp += (zc * zt) * uo[no];
    for (int nf = 0; nf < NFL; nf++) zr += one - G[i1 + (long)i1 * ld + nf * ld2];
  }
  red[0][threadIdx.x] = zk.r; red[1][threadIdx.x] = zk.i;
  red[2][threadIdx.x] = zp.r; red[3][threadIdx.x] = zp.i;
  red[4][threadIdx.x] = zr.r; red[5][threadIdx.x] = zr.i;
  __syncthreads();
  for (int st = RB / 2; st > 0; st >>= 1) {
    if (threadIdx.x < st) for (int q = 0; q < 6; q++) red[q][threadIdx.x] += red[q][threadIdx.x + st];
    __syncthreads();
  }
  if (threadIdx.x == 0) for (int q = 0; q < 6; q++) out[q] = red[q][0];
}
extern "C" {
int alfhk_corr(int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld,
               void* G, void* jc, void* inv, void* acc, void* acc0) {
  long MT = (long)N * Norb * Norb;
  long blocks = MT + ((kind == 3 || kind == 4) ? Norb : 0);
  LAUNCH(k_corr, blocks, RB, kind, k, N, Norb, NFL, ns, sq3, wr, wi, ld, (const cd*)G, (const int*)jc, (const int*)inv, (cd*)acc, (cd*)acc0);
  return fin();
}
static double* dres = nullptr;
int alfhk_scalars(int nk, void* hi, void* hj, void* hf, void* hc, int N, int Norb, int NFL, int dec, int ld, void* G, void* inv,
                  void* uo, double* hres) {
  if (!dres && hipMalloc((void**)&dres, 6 * sizeof(double)) != hipSuccess) return 3;
  LAUNCH(k_scalars, 1, RB, nk, (const int*)hi, (const int*)hj, (const int*)hf, (const cd*)hc, N, Norb, NFL, dec, ld, (const cd*)G,
         (const int*)inv, (const double*)uo, dres);
  int r = fin();
  if (r) return r;
  return hipMemcpy(hres, dres, 6 * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess ? 0 : 4;
}
}  // extern "C"

// ---- checkerboard group product of gpu_wrap_mod (k_cb_left_g / k_cb_right_g with el_left_d / el_right_d) ----------------
// Same arithmetic as the OpenMP element bodies (real elements: a*Re, a*Im sums in the order l = 1..nd; complex elements: the
// complex products and sums in the same order), same thread decomposition. One HIP launch per group instead of one libgomp
// target region (~150 us); the elements of a group have disjoint site lists, so the threads write disjoint entries.
// meta(4,e): 1 real(1)/complex, 2 nd, 3 offset into P, 4 offset into R/C; all Fortran 1-based, e 1-based.
__device__ inline void el_left(cd* M, long ld, int j1, int j2, int e, int k, const int* meta, const int* P, const double* R, const cd* C) {
  const int* mt = meta + 4L * (e - 1);
  int nd = mt[1], po = mt[2]; long mo = mt[3] + (long)(k - 1) * nd * nd;
  cd w[8];
  if (mt[0] == 1) {
    if (nd == 1) {
      int r = P[po] - 1; double a = R[mo];
      for (int j = j1; j <= j2; j++) { cd x = M[r + (long)(j - 1) * ld]; M[r + (long)(j - 1) * ld] = cd(a * x.r, a * x.i); }
    } else if (nd == 2) {
      int r1 = P[po] - 1, r2 = P[po + 1] - 1;
      double a11 = R[mo], a21 = R[mo + 1], a12 = R[mo + 2], a22 = R[mo + 3];
      for (int j = j1; j <= j2; j++) {
        cd w1 = M[r1 + (long)(j - 1) * ld], w2 = M[r2 + (long)(j - 1) * ld];
        M[r1 + (long)(j - 1) * ld] = cd(a11 * w1.r + a12 * w2.r, a11 * w1.i + a12 * w2.i);
        M[r2 + (long)(j - 1) * ld] = cd(a21 * w1.r + a22 * w2.r, a21 * w1.i + a22 * w2.i);
      }
    } else {
      for (int j = j1; j <= j2; j++) {
        for (int l = 0; l < nd; l++) w[l] = M[(P[po + l] - 1) + (long)(j - 1) * ld];
        for (int i = 0; i < nd; i++) {
          double a = R[mo + i]; double zr = a * w[0].r, zi = a * w[0].i;
          for (int l = 1; l < nd; l++) { a = R[mo + i + (long)l * nd]; zr = zr + a * w[l].r; zi = zi + a * w[l].i; }
          M[(P[po + i] - 1) + (long)(j - 1) * ld] = cd(zr, zi);
        }
      }
    }
  } else {
    if (nd == 1) {
      int r = P[po] - 1; cd c11 = C[mo];
      for (int j = j1; j <= j2; j++) M[r + (long)(j - 1) * ld] = c11 * M[r + (long)(j - 1) * ld];
    } else if (nd == 2) {
      int r1 = P[po] - 1, r2 = P[po + 1] - 1;
      cd c11 = C[mo], c21 = C[mo + 1], c12 = C[mo + 2], c22 = C[mo + 3];
      for (int j = j1; j <= j2; j++) {
        cd w1 = M[r1 + (long)(j - 1) * ld], w2 = M[r2 + (long)(j - 1) * ld];
        M[r1 + (long)(j - 1) * ld] = c11 * w1 + c12 * w2;
        M[r2 + (long)(j - 1) * ld] = c21 * w1 + c22 * w2;
      }
    } else {
      for (int j = j1; j <= j2; j++) {
        for (int l = 0; l < nd; l++) w[l] = M[(P[po + l] - 1) + (long)(j - 1) * ld];
        for (int i = 0; i < nd; i++) {
          cd z = C[mo + i] * w[0];
          for (int l = 1; l < nd; l++) z = z + C[mo + i + (long)l * nd] * w[l];
          M[(P[po + i] - 1) + (long)(j - 1) * ld] = z;
        }
      }
    }
  }
}
__device__ inline void el_right(cd* M, long ld, int i1, int i2, int istep, int e, int k, const int* meta, const int* P, const double* R, const cd* C) {
  const int* mt = meta + 4L * (e - 1);
  int nd = mt[1], po = mt[2]; long mo = mt[3] + (long)(k - 1) * nd * nd;
  cd w[8];
  if (mt[0] == 1) {
    if (nd == 1) {
      long c = P[po] - 1; double a = R[mo];
      for (int i = i1; i <= i2; i += istep) { cd x = M[(i - 1) + c * ld]; M[(i - 1) + c * ld] = cd(a * x.r, a * x.i); }
    } else if (nd == 2) {
      long c1 = P[po] - 1, c2 = P[po + 1] - 1;
      double a11 = R[mo], a21 = R[mo + 1], a12 = R[mo + 2], a22 = R[mo + 3];
      for (int i = i1; i <= i2; i += istep) {
        cd w1 = M[(i - 1) + c1 * ld], w2 = M[(i - 1) + c2 * ld];
        M[(i - 1) + c1 * ld] = cd(a11 * w1.r + a21 * w2.r, a11 * w1.i + a21 * w2.i);
        M[(i - 1) + c2 * ld] = cd(a12 * w1.r + a22 * w2.r, a12 * w1.i + a22 * w2.i);
      }
    } else {
      for (int i = i1; i <= i2; i += istep) {
        for (int l = 0; l < nd; l++) w[l] = M[(i - 1) + (long)(P[po + l] - 1) * ld];
        for (int j = 0; j < nd; j++) {
          double a = R[mo + (long)j * nd]; double zr = a * w[0].r, zi = a * w[0].i;
          for (int l = 1; l < nd; l++) { a = R[mo + l + (long)j * nd]; zr = zr + a * w[l].r; zi = zi + a * w[l].i; }
          M[(i - 1) + (long)(P[po + j] - 1) * ld] = cd(zr, zi);
        }
      }
    }
  } else {
    if (nd == 1) {
      long c = P[po] - 1; cd c11 = C[mo];
      for (int i = i1; i <= i2; i += istep) M[(i - 1) + c * ld] = M[(i - 1) + c * ld] * c11;
    } else if (nd == 2) {
      long c1 = P[po] - 1, c2 = P[po + 1] - 1;
      cd c11 = C[mo], c21 = C[mo + 1], c12 = C[mo + 2], c22 = C[mo + 3];
      for (int i = i1; i <= i2; i += istep) {
        cd w1 = M[(i - 1) + c1 * ld], w2 = M[(i - 1) + c2 * ld];
        M[(i - 1) + c1 * ld] = w1 * c11 + w2 * c21;
        M[(i - 1) + c2 * ld] = w1 * c12 + w2 * c22;
      }
    } else {
      for (int i = i1; i <= i2; i += istep) {
        for (int l = 0; l < nd; l++) w[l] = M[(i - 1) + (long)(P[po + l] - 1) * ld];
        for (int j = 0; j < nd; j++) {
          cd z = w[0] * C[mo + (long)j * nd];
          for (int l = 1; l < nd; l++) z = z + w[l] * C[mo + l + (long)j * nd];
          M[(i - 1) + (long)(P[po + j] - 1) * ld] = z;
        }
      }
    }
  }
}
// left: threads over (column tile jt, element e); right: threads over (element e, row tile it), rows it, it+ntile, ...
__global__ void k_cb_left(cd* M, long ld, int ncol, int e1, int nel, int k, int tl, int ntile, const int* meta, const int* P, const double* R, const cd* C) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)ntile * nel) return;
  int jt = (int)(t / nel) + 1, e = e1 + (int)(t % nel);
  int j1 = (jt - 1) * tl + 1, j2 = min(ncol, jt * tl);
  el_left(M, ld, j1, j2, e, k, meta, P, R, C);
}
__global__ void k_cb_right(cd* M, long ld, int nrow, int e1, int nel, int k, int ntile, const int* meta, const int* P, const double* R, const cd* C) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)ntile * nel) return;
  int e = e1 + (int)(t / ntile), it = (int)(t % ntile) + 1;
  el_right(M, ld, it, nrow, ntile, e, k, meta, P, R, C);
}
extern "C" {
int alfhk_cb_left(void* M, int ld, int ncol, int e1, int e2, int k, int tl, void* meta, void* P, void* R, void* C) {
  int nel = e2 - e1 + 1; if (nel <= 0 || ncol <= 0) return 0;
  int ntile = (ncol + tl - 1) / tl; long tot = (long)ntile * nel;
  LAUNCH(k_cb_left, nb(tot, 128), 128, (cd*)M, (long)ld, ncol, e1, nel, k, tl, ntile, (const int*)meta, (const int*)P, (const double*)R, (const cd*)C);
  return fin();
}
int alfhk_cb_right(void* M, int ld, int nrow, int e1, int e2, int k, int tl, void* meta, void* P, void* R, void* C) {
  int nel = e2 - e1 + 1; if (nel <= 0 || nrow <= 0) return 0;
  int ntile = (nrow + tl - 1) / tl; long tot = (long)ntile * nel;
  LAUNCH(k_cb_right, nb(tot, 128), 128, (cd*)M, (long)ld, nrow, e1, nel, k, ntile, (const int*)meta, (const int*)P, (const double*)R, (const cd*)C);
  return fin();
}
}  // extern "C"
