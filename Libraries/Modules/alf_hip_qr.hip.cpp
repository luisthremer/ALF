// ALF GPU=amd: blocked Householder QR (zgeqrf semantics) and application of Q (zunmqr semantics) in HIP.
//
// Why: rocSOLVER's zgeqrf/zunmqr at n = 576 launch ~2300 tiny kernels each (8.8 ms per call, ~0.12 TFLOP/s on an MI210).
// The stabilisation of the DQMC calls one of each per QR, ~217 times per sweep at L = 24, i.e. ~44% of the sweep time.
//
// Design (LAPACK conventions: H_j = I - tau_j v_j v_j^H, v_j(j) = 1, v_j(i < j) = 0, R in the upper triangle, v below;
//         zlarfg for the reflector, so the output agrees with zgeqrf up to rounding):
//   * k_qr_panel: ONE workgroup factors a panel of NB columns completely in registers (row i of the panel lives in thread
//     i mod T, register slot i / T); the two reductions per column (norm, dot products with the remaining panel columns)
//     are wave-shuffle reductions plus one LDS exchange and one barrier each.
//   * k_apply_left: applies a range of reflectors to a set of independent vectors (columns of C). One WAVE (64 lanes) owns
//     NV columns, so every reduction is a pure shuffle reduction without barrier. It is used for the trailing update of
//     the QR (NB reflectors, Q^H) and for the whole unmqr (k reflectors) -- no T factor, no GEMM, one launch each.
//   * unmqr from the right is C Q^(H) = (Q^(H)^H applied from the left to C^H)^H: transposes through a scratch buffer.
// Supported sizes: m <= 2304 (else the caller falls back to hipSOLVER). All kernels run on the default stream; the
// functions end in alf_hip_fin(): the work is synchronised before they return, or, with the deferred synchronisation (default when the
// patched libgomp plugin has the presync hook), only marked as pending and synchronised by alf_hip_presync() before the next OpenMP
// kernel launch or device copy (the OpenMP target regions run on their own queue).
#include <hip/hip_runtime.h>
#include <cstdio>

struct qc {
  double r, i;
  __device__ qc() : r(0), i(0) {}
  __device__ qc(double a, double b) : r(a), i(b) {}
};
__device__ inline qc operator+(qc a, qc b) { return qc(a.r + b.r, a.i + b.i); }
__device__ inline qc operator-(qc a, qc b) { return qc(a.r - b.r, a.i - b.i); }
__device__ inline qc operator*(qc a, qc b) { return qc(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
__device__ inline qc qconj(qc a) { return qc(a.r, -a.i); }
// wave64 sum with DPP (no LDS crossbar, unlike __shfl_xor): xor-1, xor-2, half-row mirror, row mirror, then the row sums are
// combined with the row_bcast15/31 patterns (rows 1 and 3, then 2 and 3, masked-off lanes add 0) and read from lane 63;
// the result is wave-uniform (a scalar), every lane returns the total
template <int CTRL>
__device__ inline double dpp_d(double x) {
  union { double d; int i[2]; } u, v;
  u.d = x;
  v.i[0] = __builtin_amdgcn_mov_dpp(u.i[0], CTRL, 0xf, 0xf, false);
  v.i[1] = __builtin_amdgcn_mov_dpp(u.i[1], CTRL, 0xf, 0xf, false);
  return v.d;
}
template <int CTRL, int ROWMASK>
__device__ inline double dpp_m(double x) {   // masked rows receive the DPP source, all others 0
  union { double d; int i[2]; } u, v;
  u.d = x;
  v.i[0] = __builtin_amdgcn_update_dpp(0, u.i[0], CTRL, ROWMASK, 0xf, false);
  v.i[1] = __builtin_amdgcn_update_dpp(0, u.i[1], CTRL, ROWMASK, 0xf, false);
  return v.d;
}
__device__ inline double wsum(double x) {
  x += dpp_d<0xB1>(x);        // quad_perm [1,0,3,2]
  x += dpp_d<0x4E>(x);        // quad_perm [2,3,0,1]
  x += dpp_d<0x141>(x);       // row_half_mirror
  x += dpp_d<0x140>(x);       // row_mirror
  x += dpp_m<0x142, 0xa>(x);  // row_bcast:15 -> rows 1, 3
  x += dpp_m<0x143, 0xc>(x);  // row_bcast:31 -> rows 2, 3
  union { double d; int i[2]; } u, v;
  u.d = x;
  v.i[0] = __builtin_amdgcn_readlane(u.i[0], 63);
  v.i[1] = __builtin_amdgcn_readlane(u.i[1], 63);
  return v.d;
}

extern "C" int alf_hip_fin();
static int qfin() { return alf_hip_fin(); }

// ---------------------------------------------------------------------------------------------------------------------
// Panel factorisation. P = A(j0,j0) (mp rows, pw <= NB columns, leading dimension ldp), tau = tau(j0)
// Per column: (1) the norm of the column below the diagonal (one value: DPP wave sum + LDS exchange of the waves, one
// barrier) and (2) the dot products of the new reflector with the remaining panel columns (up to 2(NB-1) real values).
// For (2) a per-value wave reduction costs ~1000 instructions per column, so the partial sums go through the LDS instead:
// every thread writes its values, thread (v, segment) sums 32 threads' entries of value v, the pairs/waves are combined in
// two more LDS steps (three barriers, ~100 LDS operations per thread).
template <int NB, int T, int R>
__global__ __launch_bounds__(T) void k_qr_panel(int mp, int pw, qc* P, long ldp, qc* tau) {
  constexpr int NW = T / 64;
  constexpr int NVMAX = 2 * (NB - 1);
  constexpr int STR = NVMAX + 1;       // odd stride: conflict-free writes
  __shared__ double red[2][NW][2];
  __shared__ double alS[2][2];
  __shared__ double Pl[T * STR];
  __shared__ double red2[NW][32];
  __shared__ double totS[32];
  const int t = threadIdx.x, w = t >> 6, ln = t & 63;
  qc a[R][NB];
#pragma unroll
  for (int r = 0; r < R; r++) {
    int i = t + r * T;
#pragma unroll
    for (int c = 0; c < NB; c++) a[r][c] = (i < mp && c < pw) ? P[i + (long)c * ldp] : qc(0, 0);
  }
  int rb = 0;
#pragma unroll
  for (int c = 0; c < NB; c++) {
    if (c < pw) {   // (no 'break': it stops the full unrolling and puts a[][] into scratch memory)
      // ---- reflector of column c (zlarfg) ----
      double sq = 0;
#pragma unroll
      for (int r = 0; r < R; r++) {
        int i = t + r * T;
        if (i > c) sq += a[r][c].r * a[r][c].r + a[r][c].i * a[r][c].i;     // rows >= mp hold zeros
      }
      sq = wsum(sq);
      if (ln == 0) red[rb][w][0] = sq;
      if (t == c) { alS[rb][0] = a[0][c].r; alS[rb][1] = a[0][c].i; }
      __syncthreads();
      double xn2 = 0;
#pragma unroll
      for (int q = 0; q < NW; q++) xn2 += red[rb][q][0];
      const double alr = alS[rb][0], ali = alS[rb][1];
      rb ^= 1;
      qc tc(0, 0), sc(1, 0);
      double beta = alr;
      if (!(xn2 == 0.0 && ali == 0.0)) {
        double nrm = sqrt(alr * alr + ali * ali + xn2);
        beta = (alr >= 0.0) ? -nrm : nrm;
        tc = qc((beta - alr) / beta, -ali / beta);
        double dr = alr - beta, d2 = dr * dr + ali * ali;
        sc = qc(dr / d2, -ali / d2);
      }
#pragma unroll
      for (int r = 0; r < R; r++) {
        int i = t + r * T;
        if (i > c) a[r][c] = a[r][c] * sc;
      }
      if (t == c) { a[0][c] = qc(beta, 0.0); tau[c] = tc; }
      // ---- apply H^H = I - conj(tau) v v^H to the remaining panel columns ----
      if (c + 1 < pw) {
        qc vv[R];     // the reflector, masked: 0 above the diagonal, 1 on it (no predicates in the loops below)
#pragma unroll
        for (int r = 0; r < R; r++) {
          int i = t + r * T;
          vv[r] = (i == c) ? qc(1, 0) : ((i > c) ? a[r][c] : qc(0, 0));
        }
#pragma unroll
        for (int k = c + 1; k < NB; k++) {
          if (k < pw) {
            double xr = 0, xi = 0;
#pragma unroll
            for (int r = 0; r < R; r++) {
              xr += vv[r].r * a[r][k].r + vv[r].i * a[r][k].i;
              xi += vv[r].r * a[r][k].i - vv[r].i * a[r][k].r;
            }
            Pl[t * STR + 2 * (k - c - 1)] = xr;
            Pl[t * STR + 2 * (k - c - 1) + 1] = xi;
          }
        }
        __syncthreads();
        const int nv = 2 * (pw - 1 - c);
        {
          const int v = t & 31, sg = t >> 5;
          double part = 0;
          if (v < nv) {
#pragma unroll 8
            for (int q = 0; q < 32; q++) part += Pl[(sg * 32 + q) * STR + v];
          }
          part += __shfl_xor(part, 32, 64);
          if (ln < 32) red2[w][v] = part;
        }
        __syncthreads();
        if (t < 32) {
          double tot = 0;
#pragma unroll
          for (int q = 0; q < NW; q++) tot += red2[q][t];
          totS[t] = tot;
        }
        __syncthreads();
        const qc ct = qconj(tc);
#pragma unroll
        for (int k = c + 1; k < NB; k++) {
          if (k < pw) {
            qc s(totS[2 * (k - c - 1)], totS[2 * (k - c - 1) + 1]);
            qc coef = ct * s;
#pragma unroll
            for (int r = 0; r < R; r++) a[r][k] = a[r][k] - vv[r] * coef;
          }
        }
      }
    }
  }
#pragma unroll
  for (int r = 0; r < R; r++) {
    int i = t + r * T;
#pragma unroll
    for (int c = 0; c < NB; c++)
      if (i < mp && c < pw) P[i + (long)c * ldp] = a[r][c];
  }
}

// ---------------------------------------------------------------------------------------------------------------------
// C(:, 0..nc-1) <- (product of reflectors j = 0..k-1 of V) applied from the left, one wave per NV columns.
//   asc != 0: j = 0,1,..,k-1 else j = k-1,..,0;  conjtau != 0: use conj(tau_j).   (Left,N: desc, tau; Left,C: asc, conj(tau))
// V(i,j) for i > j are the vectors (unit diagonal implied, upper part ignored); nq = number of rows of C and V.
template <int Q, int NV>
__global__ __launch_bounds__(256) void k_apply_left(int nq, int k, int nc, int asc, int conjtau, const qc* V, long ldv,
                                                    const qc* tau, qc* C, long ldc) {
  const int lane = threadIdx.x & 63;
  const int wave = (blockIdx.x * blockDim.x + threadIdx.x) >> 6;
  const int c0 = wave * NV;
  if (c0 >= nc) return;
  qc y[NV][Q];
#pragma unroll
  for (int v = 0; v < NV; v++)
#pragma unroll
    for (int q = 0; q < Q; q++) {
      int i = lane + 64 * q;
      y[v][q] = (c0 + v < nc && i < nq) ? C[i + (long)(c0 + v) * ldc] : qc(0, 0);
    }
  for (int s = 0; s < k; s++) {
    const int j = asc ? s : k - 1 - s;
    qc t = tau[j];
    if (t.r == 0.0 && t.i == 0.0) continue;
    if (conjtau) t = qconj(t);
    qc vv[Q];
#pragma unroll
    for (int q = 0; q < Q; q++) {
      int i = lane + 64 * q;
      vv[q] = (i > j && i < nq) ? V[i + (long)j * ldv] : (i == j ? qc(1, 0) : qc(0, 0));
    }
    qc dot[NV];
#pragma unroll
    for (int v = 0; v < NV; v++) {
      double xr = 0, xi = 0;
#pragma unroll
      for (int q = 0; q < Q; q++) {  // conj(vv) * y
        xr += vv[q].r * y[v][q].r + vv[q].i * y[v][q].i;
        xi += vv[q].r * y[v][q].i - vv[q].i * y[v][q].r;
      }
      dot[v] = qc(wsum(xr), wsum(xi));
    }
#pragma unroll
    for (int v = 0; v < NV; v++) {
      qc coef = t * dot[v];
#pragma unroll
      for (int q = 0; q < Q; q++) y[v][q] = y[v][q] - vv[q] * coef;
    }
  }
#pragma unroll
  for (int v = 0; v < NV; v++)
#pragma unroll
    for (int q = 0; q < Q; q++) {
      int i = lane + 64 * q;
      if (c0 + v < nc && i < nq) C[i + (long)(c0 + v) * ldc] = y[v][q];
    }
}

// B(j,i) = conj(A(i,j)),  A m x n (lda), B n x m (ldb)
__global__ void k_ctranspose(int m, int n, const qc* A, long lda, qc* B, long ldb) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)m * n) return;
  int i = t % m, j = t / m;
  B[j + (long)i * ldb] = qconj(A[i + (long)j * lda]);
}

#define QLAUNCH(kern, grid, blk, ...) hipLaunchKernelGGL(kern, dim3(grid), dim3(blk), 0, 0, __VA_ARGS__)

template <int Q, int NV>
static void apply_inst(int nq, int k, int nc, int asc, int conjtau, const qc* V, long ldv, const qc* tau, qc* C, long ldc) {
  long waves = (nc + NV - 1) / NV;
  QLAUNCH((k_apply_left<Q, NV>), (waves + 3) / 4, 256, nq, k, nc, asc, conjtau, V, ldv, tau, C, ldc);
}
// returns false if nq is too large
static bool apply_left(int nq, int k, int nc, int asc, int conjtau, const qc* V, long ldv, const qc* tau, qc* C, long ldc) {
  if (nc <= 0 || k <= 0) return true;
  int Q = (nq + 63) / 64;
#define AP(QQ, NVV) apply_inst<QQ, NVV>(nq, k, nc, asc, conjtau, V, ldv, tau, C, ldc); return true
  if (Q <= 1) { AP(1, 4); }
  if (Q <= 2) { AP(2, 4); }
  if (Q <= 3) { AP(3, 4); }
  if (Q <= 4) { AP(4, 4); }
  if (Q <= 5) { AP(5, 4); }
  if (Q <= 6) { AP(6, 4); }
  if (Q <= 9) { AP(9, 4); }
  if (Q <= 12) { AP(12, 2); }
  if (Q <= 16) { AP(16, 2); }
  if (Q <= 24) { AP(24, 1); }
  if (Q <= 36) { AP(36, 1); }
#undef AP
  return false;
}

template <int NB, int T, int R>
static void panel_inst(int mp, int pw, qc* P, long ldp, qc* tau) {
  QLAUNCH((k_qr_panel<NB, T, R>), 1, T, mp, pw, P, ldp, tau);
}

extern "C" {
// A (m x n, lda) <- QR factors in LAPACK zgeqrf format, tau(min(m,n)). Return 0 ok, 3 = size not supported (A untouched).
int alfhq_geqrf(int m, int n, void* Ad, int lda, void* taud) {
  qc* A = (qc*)Ad;
  qc* tau = (qc*)taud;
  int k = m < n ? m : n;
  if (k <= 0) return 0;
  if (m > 2304) return 3;
  const int NB = (m <= 1024) ? 16 : 8;
  for (int j0 = 0; j0 < k; j0 += NB) {
    int pw = (k - j0 < NB) ? (k - j0) : NB, mp = m - j0;
    qc* P = A + j0 + (long)j0 * lda;
    if (NB == 16) {
      int R = (mp + 255) / 256;
      if (R <= 1) panel_inst<16, 256, 1>(mp, pw, P, lda, tau + j0);
      else if (R == 2) panel_inst<16, 256, 2>(mp, pw, P, lda, tau + j0);
      else if (R == 3) panel_inst<16, 256, 3>(mp, pw, P, lda, tau + j0);
      else panel_inst<16, 256, 4>(mp, pw, P, lda, tau + j0);
    } else {
      int R = (mp + 511) / 512;
      if (R <= 1) panel_inst<8, 512, 1>(mp, pw, P, lda, tau + j0);
      else if (R == 2) panel_inst<8, 512, 2>(mp, pw, P, lda, tau + j0);
      else if (R == 3) panel_inst<8, 512, 3>(mp, pw, P, lda, tau + j0);
      else if (R == 4) panel_inst<8, 512, 4>(mp, pw, P, lda, tau + j0);
      else panel_inst<8, 512, 5>(mp, pw, P, lda, tau + j0);
    }
    int nc = n - j0 - pw;
    if (nc > 0) apply_left(mp, pw, nc, 1, 1, P, lda, tau + j0, A + j0 + (long)(j0 + pw) * lda, lda);
  }
  return qfin();
}

static qc* qscr = nullptr;
static long qscr_n = 0;

// C <- Q C, Q^H C (side 0) or C Q, C Q^H (side 1);  trans 0 = N, else C.  A holds k reflectors (nq x k, lda), nq = m or n.
int alfhq_unmqr(int side, int trans, int m, int n, int k, void* Ad, int lda, void* taud, void* Cd, int ldc) {
  const qc* A = (const qc*)Ad;
  const qc* tau = (const qc*)taud;
  qc* C = (qc*)Cd;
  if (m <= 0 || n <= 0 || k <= 0) return 0;
  int nq = side == 0 ? m : n;
  if (nq > 2304) return 3;
  if (side == 0) {
    // Left,N: descending with tau; Left,C: ascending with conj(tau)
    apply_left(m, k, n, trans != 0, trans != 0, A, lda, tau, C, ldc);
    return qfin();
  }
  // right: C Q = (Q^H C^H)^H,  C Q^H = (Q C^H)^H
  long need = (long)m * n;
  if (qscr_n < need) {
    if (qscr) hipFree(qscr);
    qscr = nullptr; qscr_n = 0;
    if (hipMalloc((void**)&qscr, need * sizeof(qc)) != hipSuccess) return 4;
    qscr_n = need;
  }
  QLAUNCH(k_ctranspose, (need + 255) / 256, 256, m, n, (const qc*)C, (long)ldc, qscr, (long)n);   // qscr(n x m) = C^H
  int tl = (trans == 0);   // right,N -> left,C on C^H ; right,C -> left,N
  apply_left(n, k, m, tl, tl, A, lda, tau, qscr, n);
  QLAUNCH(k_ctranspose, (need + 255) / 256, 256, n, m, (const qc*)qscr, (long)n, C, (long)ldc);
  return qfin();
}
}  // extern "C"
