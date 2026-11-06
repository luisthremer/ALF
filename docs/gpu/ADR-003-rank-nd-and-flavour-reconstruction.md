# ADR-003: non-diagonal / rank > 1 vertices (PARKING #15) and flavour reconstruction (PARKING #16) on the GPU

Status: PROPOSAL, 2026-11-06, branch `gpu-port-delayed-submatrix-rank-nd` (from `gpu-port-delayed-submatrix`). No source edited yet.
Written from reading the code only; nothing has been built or run for this note (the laptop has no GPU; `nividia`/`intel` not touched).

## 0. Findings that change the brief

1. **Spin_Peierls is not opened by #15 and #16.** Its J-vertices are field **type 4** (`Hamiltonian_Spin_Peierls_smod.F90:460`, two HS
   fields per vertex), which `gpu_supported` refuses (PARKING #29). Its `Calc_Fl` block is dead code: line 312 terminates the run
   when `.not. SU2_Symm`, and the `Calc_Fl` allocation sits behind the same condition. So Spin_Peierls needs #15 **plus type 4**.
   Type 4 is small (see 4.6) but it is a third item.
2. **#16 alone opens Plain_Vanilla** (diagonal rank-1, type 2) whenever it sets `Calc_Fl` (half filling `Ham_chem = 0`, or `U < 0`;
   `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90:252`). `Hamiltonian_Hubbard_smod` has the block commented out (line 289).
3. **`test_vs_ed` has no reference for any of the four models.** `testsuite/test_vs_ed/test_specs.yaml` has only `Hubbard` (finite T,
   T=0, MALA). The required "test_vs_ed for tV, Kondo, Spin_Peierls, Plain_Vanilla" needs new, independent ED references (small
   Python ED for 4-site chains / 2x2). Plain_Vanilla at half filling is the Hubbard model, so the existing Hubbard ED energy
   (-1.47261997, 4-site chain, beta 2) applies to it. The new ED codes need their own validation (CPU vs ED first). See open question Q2.
4. **The slice-level kernels are the bigger part of #15**, not only the per-site update. `gpu_Vslice_R/L`, `gpu_Vwrap_all`, `gpu_Vop_R/L`
   (stabilization, tau pass, wrap hook; `gpu_wrap_mod.F90`) apply, per matrix row, the factors of all operators containing that row.
   That is valid only for diagonal operators. Non-diagonal vertices couple rows P(1..N) through `Op%U`.

## 1. What the CPU does (the specification)

Per vertex n and flavour (Wrapgr_mod.F90:198-270, Operator_mod.F90 `Op_Wrapup`):
1. `Op_Wrapup(G, Op, f_old, N_Type=1)`: `G <- D(sp) U^H G U D(-sp)` on rows/columns P (diagonal operator: scalings by `E_exp(i,+-sp)`).
2. `Upgrade2` (or `Upgrade_lowrank`) acts **in the eigenbasis**: r = `Op%N_non_zero` (<= N), `Delta_m = exp(Z1 E_m) - 1`,
   `Mat = I + (I - G_PP) Delta`, ratio = det(Mat) * exp(Z1 alpha); accept: G <- G - G(:,P) W (E_P^T - G(P,:)) with a r x r matrix W
   (r = 1: W = Delta / Mat; for r > 1 the closed form of W is written down in `docs/algorithm-notes.md` sections 10-12 on the CPU branch and is
   **checked numerically against the recursive x_v/y_v + ZGEMM code of `Upgrade2` in A0 before it is used**; I have not verified the
   ordering of Delta and Mat^-1 for r > 1 here).
3. `Op_Wrapup(..., N_Type=2)`: `G <- U G U^H` on P (no-op for diagonal).
4. `weight_reconstruction(Ratio)` after the flavour loop if `reconstruction_needed`.

Low-rank CPU code (`Lowrank_update_mod.F90`): eligibility = support N <= 8, any structure; submatrix with non-diagonal content is
flushed at once (K_ND = 1) because of accuracy (documented there); delayed has no such loss. Hence **for non-diagonal models delayed
is the scheme that matters; submatrix gives nothing on them.**

## 2. What the device code assumes today

| Place | Assumption |
|---|---|
| `gpu_supported` (gpu_state_mod:479) | refuses `.not. diag`, `N_non_zero > 1`, type /= 1,2, `reconstruction_needed` |
| `build_V` mirror | `V_meta(N, poff, eoff, skip)`, `V_P`, `V_exp` = real table `E_exp(i,sp)`; **no `U`, no `E`, no N_non_zero, no diag flag** |
| `gpu_update_mod:site_update` | rank 1 only: one `G(p,p)` D->H per flavour, `k_upd_prep/geru` (rank-1 x y^T); N > 1 diagonal handled by `gpu_Vop_wrap` (two kernels), `N_Type 2` omitted |
| `gpu_lowrank_mod` | fast paths (`k_dly_*`, `k_sub_*`, block mode `blk_*`) only for `V_meta(1,n,nf) == 1`; other vertices flush first (`gpu_lr_site_ok`) |
| `gpu_wrap_mod` slice kernels | per-row product of diagonal factors; ordering by operator index |
| `gpu_obs_mod` | Hubbard only (N_FL = 1, or 2 with N_SUN = 1); other models measure on the host after D->H of G |
| `gpu_tau_mod` | refuses `reconstruction_needed`; host shadows + host `OBSERT` for non-Hubbard models |

Already in place and reusable: `N_FL_eff` / `Calc_Fl_map` indexing in stabilization, tau pass and update loops (udv arrays are
`N_FL_eff` deep); G_d is allocated for all `N_FL` flavours; `la_gemm` backend (cuBLAS / oneMKL); `gpu_Vop_wrap` for N > 1.

## 3. Stage B: flavour reconstruction (smaller, do first)

Touch points (all small, no new device kernel needed because every consumer of the reconstructed block is a host routine):
1. `gpu_supported`: remove the `reconstruction_needed` refusal. Keep refusing the cases still unsupported, written as explicit
   checks: `Calc_Fl` together with the Hubbard **device** observables (`gpu_obs_setup` must refuse when `reconstruction_needed`, because
   it reads G_d for all flavours), and with global/Langevin moves (already refused).
2. `gpu_update_mod:site_update`: `if (reconstruction_needed) call ham%weight_reconstruction(Ratio)` before `Product(Ratio)`
   (as `upgrade_mod.F90:208`). Ratio is host data; `Ratio` is initialised to 1 for the non-computed flavours; block mode and the delayed
   paths all go through `site_update`, so one call site.
3. `gpu_tau_mod`: drop `reconstruction_needed` from `tau_off_refuse`; add `ham%GR_reconstruction` / `GRT_reconstruction` on the host shadows
   (G00_T, GTT_T, GT0_T, G0T_T) after the Symm step and before `OBSERT`, at `nt = 0` (`tau_begin`) and in `tau_obs`, exactly as `tau_m_mod.F90:119-130`.
4. `main.F90` already calls `GR_reconstruction` on the host GR after `gpu_g_host(GR)` (lines 1085-1157); `weight_reconstruction(Phase_array)`
   is host. Verify only: `gpu_g_host` must not copy the non-computed flavour block over the reconstructed one, and `g_in` must not
   rely on that block.
Risk: low. Test: kernel level = GR / GT0 / G0T of the non-computed flavour equal to the CPU after wrap, stabilization and a sweep segment
(fixed fields); observable level = Plain_Vanilla at `Ham_chem = 0` vs ED (Hubbard value) and vs CPU. Plant: skip the weight call (fails the
energy check at about 1e-2), swap `nf_calc`/`nf_reconst` in the reconstruction (fails the Green-function comparison).

## 4. Stage A: non-diagonal and rank > 1 vertices

### 4.1 Data layout (plain contiguous arrays)
Extend `build_V` only (no derived types): per (n,nf) add to `V_meta`: `N_non_zero`, `diag` flag, `uoff`; new arrays
`V_U(uoff + i + (j-1)*N)` (complex, the full N x N `Op%U`, columns ordered as `Op_set` leaves them: non-zero eigenvalues first) and
`V_E(eoff2 + i)` (real eigenvalues, needed for `Delta_m` on the host - the host keeps using `Op_V%E` there, so `V_E` is only needed if a
device-side ratio is added later; **do not add it now**). `V_exp` (table `E_exp(i,sp)`) is already the factor `Op_Wrapup` uses for the
non-diagonal case too (`VH1(:,i) = U(:,i) * E_Exp(i,-sp)`), so it is reused unchanged. For types 1 and 2 only. Offsets stay int32.

### 4.2 General wrap (Op_Wrapup N_Type 1 and 2, rows P and columns P of G)
Two kernels per step, same structure as `k_vop_wrap_rows/cols`, but with the dense N x N block:
columns: `G(:,P) <- G(:,P) * (U D(-sp))` for all rows (parallel over rows, N x N matvec per row in registers; N <= 8 for eligible, any N allowed);
rows: `G(P,:) <- (D(sp) U^H) G(P,:)` (parallel over columns). Order columns then rows as ZSLGEMM does ('r' then 'l'); the corner block
differs from ALF at the rounding level only. N_Type 2: `G <- U G U^H` with two more kernels. Used in `site_wrap` (replacing the
diagonal-only branch) and in `gpu_wrap_mod` `gpu_Vop_wrap`/hooks.

### 4.3 Rank-r ratio and update (standard scheme first)
- Ratio: gather the r x r block `G(P,P)` into a small device buffer (one kernel), one D->H of `r*r` elements per flavour (was 1), the
  ratio code of `Upgrade2` on the host (the same lines, determinant branches for r = 1, 2, general), weight reconstruction, `ranf_wrap`
  in ALF's order.
- Accept: host builds `W = Delta Mat^-1` (r x r, `la`-free small LU, same rounding class as `Upgrade2`), uploads W (r*r complex, one
  H->D, like the 2 reals today); device: gather `C = G(:,P)` (Ndim x r), `R = E_P^T - G(P,:)` (r x Ndim), `G <- G - C (W R)` by two
  `la_gemm` calls (`Ndim x r x r`, then `Ndim x Ndim x r`) or one fused kernel for small r. r = 1 keeps today's `k_upd_prep/geru` path
  unchanged (bit-identical to the present GPU code for Hubbard).
- Field flip on the host, one scalar H->D, as now.
Eligibility for the existing fast paths stays `N = 1`; every other vertex goes `flush -> general wrap -> rank-r update -> general wrap
(N_Type 2)`. That is the **standard path for general vertices** and it works under all three `Sequential_scheme` values at once, because
`gpu_lr_site_ok` already flushes before non-eligible vertices. Delayed/submatrix then give no speed-up on tV/Kondo until stage A3, but
they are correct.

### 4.4 Slice kernels (stabilization, tau pass)
`gpu_Vslice_R/L`, `gpu_Vwrap_all`: for a slice with non-diagonal operators, build once (at setup, like `hop_gptr`) **groups = maximal runs of
consecutive operators with pairwise disjoint P**; the diagonal operators stay in the existing fused per-row kernel; non-diagonal groups use the
4.2 kernels, one launch per group, groups applied in operator order (ordering of overlapping operators is preserved by construction). Mixed
slices (Kondo, Spin_Peierls: diagonal U-vertices plus bond vertices) are applied in operator order: runs of diagonal operators fused as now,
non-diagonal runs per group. Needs `Op_mmultL/R` semantics for non-diagonal operators (`Operator_mod.F90:558-720`): `M(P,:) <- U D(sp) U^H M(P,:)` (and the transposed/conjugate variants and `sign`); to be read out and unit-tested against the host routine.
Expect a launch-count regression for non-diagonal models (about one launch per group per slice per matrix); acceptable for stage A, to be measured.

### 4.5 Refusal logic
Replace the three blanket refusals by checks that say exactly what is unsupported: field type not in {1,2} (until 4.6), `N > N_MAX_DEV`
(proposed 8, the CPU low-rank limit; larger operators refused with the reason), non-diagonal operator with `g_t`, and keep everything
else. `GPU=none` build: none of this code is compiled (`#ifdef ALF_GPU`); the CPU files are untouched, so bit-identity with upstream is by
construction and checked by the existing `acc_gpunone.sh`.

### 4.6 Type 4 (needed for Spin_Peierls; optional, owner's decision)
`Fields_Phi` type 4: `Phi = Phi_st(Nint(Re f),2) * sqrt(1 + Im f)`-style function of the complex field (`Fields_mod.F90:127`), plus `Gama`.
The device wrap kernels read the field from `nsigma_d` (already complex) and use the table for types 1, 2 only. Type 4 needs `exp(phi g E_i)`
evaluated on the device (no table: the argument is continuous in the imaginary part) and the host ratio code is unchanged (`nsigma%Phi`).
The flip proposal (`nsigma%flip`) stays on the host. No MALA/Langevin (type 3) is involved.

### 4.7 Delayed / submatrix for rank > 1 (performance, after correctness)
A3 (delayed, rank r): extend `k_dly_*` to append r columns per accept (CPU: `Lowrank_update_mod` section 9-11) and the wrap hook for a
general wrap (`X(P,:) <- L X(P,:)`, `Y(P,:) <- R^T Y(P,:)`). A4 (submatrix): not proposed - the CPU itself flushes non-diagonal
content immediately (K_ND = 1), so there is nothing to gain. Decision point after stage A2 with measured launch counts.

## 5. Staged plan and order

Order decided from the code: **B first** (4 touch points, no new kernels, opens Plain_Vanilla), then **A**.

| Stage | Content | Opens | New kernels |
|---|---|---|---|
| B | flavour reconstruction (section 3) | Plain_Vanilla | none |
| A0 | mirror (`V_U`, flags), `gpu_supported` rewritten (4.5); kernel self-tests of Op_Wrapup N_Type 1/2 and Op_mmultL/R vs host | - | wrap/mult kernels for dense N x N |
| A1 | standard path for general vertices: wrap, rank-r ratio/accept (4.2, 4.3), slice kernels (4.4) | tV, Kondo with `standard`, `delayed`, `submatrix` (correct, no speed-up) | gather, small gemm |
| A2 | acceptance: ED / CPU statistics on A30 and B580, `Nwrap` 4 for tV | tV, Kondo validated | - |
| A3 | type 4 (4.6) | Spin_Peierls | on-the-fly exp(phi g E) |
| A4 | delayed rank r (4.7), only if A2 shows a need | speed | extend `k_dly_*` |

Not started: MPI, multi-chain, tuning. Both compilers each stage; production optimisation level; nvfortran `-O2` `blk_end` pin kept
(any new routine in `gpu_lowrank_mod.F90` is tested at -O2/-O3 on nvfortran; new kernels are added to `gpu_update_mod`/`gpu_wrap_mod`, not to `blk_end`'s neighbourhood).

## 6. Verification plan (pre-registration comes next, before any run)
- One acceptance file per stage in `testsuite/gpu_reference/acceptance/` (names to be allocated: 401 B, 402 A1/A2 unless the owner prefers
  other numbers), written and committed before any result, with the positive controls of VERIFICATION.md section 5.
- Kernel level (fixed fields, relative 1e-10): G after one wrap, after stabilization, after a sweep segment, per model and per
  `Sequential_scheme`; flavour block reconstructed for Plain_Vanilla; `G` from the general wrap vs `Op_Wrapup`; rank-r update vs `Upgrade2`.
- Observable level: ED energy within 3 sigma, chi^2 reported; CPU reference with >= 10x statistics; no bin-by-bin comparison across compilers;
  schemes bin-for-bin where they were before (Hubbard unchanged: r = 1 path untouched).
- Planted defects (each must fail its gate): drop `N_Type 2` rotation; use `U` instead of `U^H`; wrong `V_exp` sign index; swap the column order of `U`;
  skip the weight reconstruction; stale `nsigma_d` after accept; rank-r update with W transposed.
- Dev/acceptance split, `env -i`, provenance stamps, T0/T1/T2 tiers, verifier and reviewer roles exactly as in VERIFICATION.md.

## 7. Open questions for the owner
Q1. Spin_Peierls needs type 4 (A3). Include it in this work, or stop after tV and Kondo?
Q2. ED references for tV, Kondo, Spin_Peierls do not exist (only Hubbard has one). Write independent small ED codes (Python, 4-site chain / 2x2) and validate them against the CPU
    build first, or accept CPU-vs-GPU statistics with >= 10x CPU statistics as the observable gate for these models (ED stays only for Plain_Vanilla via the Hubbard value)? My recommendation: write the ED references, since the CPU result itself is then also tested.
Q3. For tV/Kondo the delayed/submatrix schemes give no speed-up after A1 (all their vertices are non-diagonal or N > 1 and flush). Do A4 (delayed rank r) now or after A2?
Q4. Branch hygiene: this worktree's auto-created branch `claude/gpu-port-delayed-submatrix-rank-nd-353e76` is not used; work is on `gpu-port-delayed-submatrix-rank-nd`
    (created from `gpu-port-delayed-submatrix`, 2ac9d670). OK to leave the other one for deletion?
