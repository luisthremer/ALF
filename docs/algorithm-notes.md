# Delayed and submatrix local updates in ALF — algorithm notes

Status: derivation for Checkpoint 1. Numerically checked (Sec. 6) and cross-checked against an
independent re-derivation (Sec. 7): no open mismatches.

References (full citations in `CLAUDE.md`):
[SX-D] Sun & Xu, PRB 109, 235140 (2024), arXiv:2308.12005 (delay update in DQMC).
[SX-S] Sun & Xu, SciPost Phys. 18, 055 (2025), arXiv:2404.09989 (submatrix update in DQMC).
[G11] Gull et al., PRB 83, 075122 (2011), arXiv:1010.3690 (submatrix updates, CT-AUX).
[N09] Nukala et al., PRB 80, 195111 (2009) (submatrix-type fast update, Hirsch–Fye).
[A08] Alvarez et al., SC'08 (delayed updates, Hirsch–Fye/DCA).

## 0. ALF conventions used here

* Per flavour `nf` (colour index suppressed, N_SUN-fold degenerate):
  `G = (1 + B(τ,0) B(β,τ))^{-1}`, `G_{xy} = <c_x c_y^†>`  (`model.tex`, Eq. `Green_eq`).
* One time slice: `B_τ = e^{V_M} ⋯ e^{V_1} e^{-Δτ T}`  (`model.tex`, Eq. `Btau.eq`).
  Vertex `n` = `Op_V(n,nf)`: sites `P(1:N)`, matrix `O = U diag(E) U^†`, coupling `g`
  (or `g_t(τ)`), shift `alpha`, field `φ = nsigma%Phi(n,τ)`, `e^{V_n} = P^T e^{g φ O} P`.
  `N_non_zero` eigenvalues are non-zero, stored first (`E(1..N_non_zero)`).
* Weight (`Op_phase`, `Upgrade2` doc):
  `W(C) = [ Π_nf ( Π_{n,τ} e^{g α φ} ) det(1 + B_nf) ]^{N_SUN} Π_{n,τ} γ(φ)` (× `e^{-S_0}`).
* G is carried through the slice **vertex by vertex** (`Wrapgr_mod.F90`):
  * `WRAPGRUP` (τ → τ+1): hopping; then for `n = 1..M`:
    `Op_Wrapup(type 1, old φ)` = rotate rows/cols `P` into the eigenbasis and propagate
    `G ← e^{V_n} G e^{-V_n}`; propose/update at `n`; `Op_Wrapup(type 2)` = rotate back.
  * `WRAPGRDO` (τ → τ−1): for `n = M..1`: `Op_Wrapdo(type 2)` = rotate; propose/update at `n`;
    `Op_Wrapdo(type 1, new φ)`: `G ← e^{-V_n} G e^{V_n}`; hopping at the end.
  * For a **diagonal** operator (`Op%diag`, `U = 1`) type 2 is a no-op and type 1 is
    `G(P(i),:) *= t_i`, `G(:,P(i)) /= t_i` with `t_i = e^{± g φ E_i}`: a diagonal similarity
    transform acting only on rows/columns `P`.

In the update step the current G is in the eigenbasis of vertex `n` on rows/cols `P`.
A rank-1 vertex has `N_non_zero = 1`; its active index is `p = P(1)`, eigenvalue `E_1`.

## 1. Acceptance ratio (rank-1 vertex)

Field change `φ → φ'`. Per computed flavour `f` (`Calc_Fl_map`):

    Z1_f = g_f (φ' − φ),     Δ_f = e^{Z1_f E_1} − 1,

and in the vertex frame `e^{V'} = (1 + Δ_f e_p e_p^T) e^{V}`, so `B' = (1 + Δ_f e_p e_p^T) B`.
With `G = (1+B)^{-1}`:

    1 + B' = [1 + Δ_f e_p e_p^T (1 − G)] (1 + B)
    det(1+B')/det(1+B) = 1 + Δ_f (1 − G_pp) =: d_f            (Sylvester)

`Upgrade2` forms exactly this as `Mat(1,1) = myexp − Z·G(P,P)` with `myexp = 1+Δ`, `Z = Δ`.
Full ratio and Metropolis weight, **as `Upgrade2` computes them**:

    R_f      = d_f · e^{Z1_f α_f}
    (ham%weight_reconstruction fills non-computed flavours if N_FL_eff < N_FL)
    R        = (Π_f R_f)^{N_SUN} · γ(φ')/γ(φ)
    Weight   = S0_ratio · T0_proposal_ratio · | Re(Phase·R) / Re(Phase) |
    accept  ⇔ Weight > ranf_wrap()
    on accept: Phase ← Phase · R / sqrt(R · conj(R))

The only matrix data needed is **one complex number per computed flavour: `G_pp` of the current G.**
The new schemes compute `G_pp` differently (Secs. 3, 4) and must then evaluate the scalar
expressions above in the same order and form as `Upgrade2` does. Random numbers are consumed in
the same order as in the standard update: `nsigma%flip` (types 2, 3), then the `T0_proposal` draw, then
the Metropolis draw.

## 2. Rank-1 G update (standard; reference)

`G' = (1 + B')^{-1} = G [1 + Δ e_p e_p^T (1 − G)]^{-1}`; Sherman–Morrison:

    G' = G − G(:,p) · (Δ/d) · (e_p^T − G(p,:))

This is `Upgrade2`'s `x_v(p) = u/(1+v u) = Δ/d`, `y_v = e_p − G(p,:)`, `ZGERU` with `−x_v(p)`.

## 3. Delayed update

After `j` accepted rank-1 flips since the last flush (per flavour):

    G_j = G0 − X_j Y_j^T,       X_j, Y_j ∈ C^{N×j}.

* Proposal at `p`: `G_j(p,p) = G0(p,p) − X(p,1:j)·Y(p,1:j)`  — O(j).
* Accept (Δ, d from Sec. 1):

      c = G_j(:,p) = G0(:,p) − X Y(p,:)^T      (ZGEMV, O(Nj))
      r = G_j(p,:) = G0(p,:) − X(p,:) Y^T      (ZGEMV, O(Nj))
      X_{j+1} = [X, (Δ/d) c],   Y_{j+1} = [Y, e_p − r^T]

  This is Sec. 2 applied to `G_j`, so it is exact at every step.
* Flush when `j = n_d` or at a forced flush point: `G0 ← G0 − X Y^T` (ZGEMM, O(N² j)), `j ← 0`.
* **Wraps between vertices.** A vertex wrap is `G ← L G R`, where L and R act only on rows/cols P. To keep
  `G = G0 − X Y^T`: `G0 ← L G0 R`, `X ← L X`, `Y ← R^T Y`. For a diagonal operator:
  `X(P(i),:) *= t_i`, `Y(P(i),:) /= t_i`  — O(j) per wrap. With no flush per vertex, the identity
  holds exactly across all vertices of a slice.
* Repeated sites in a block are allowed: each step is an exact rank-1 update of the current G.
* Rank-N_P extension (not built now): each accepted vertex adds `N_P` columns to X, Y and the
  ratio needs the `N_P×N_P` block `G_j(P,P)`; Upgrade2's Woodbury algebra carries over.
  Storage is therefore organised as column counts, not "number of accepted flips".

Cost per slice (N_acc accepted): ratio O(M·j̄); new columns O(N_acc N j̄) (level 2);
flushes O(N_acc N²) (level 3). [SX-D] Table I: extra overhead `O(β n_d N²)` per sweep at level 1/2.

## 4. Submatrix update

Accepted flips at **distinct** active sites `S = (s_1..s_k)` since the last flush, `D = diag(Δ_1..Δ_k)`.
For distinct sites the factors `(1 + Δ_j e_{s_j} e_{s_j}^T)` commute and multiply to
`1 + E_S D E_S^T` (`E_S`: N×k selector), so `G_k = G0 [1 + E_S D E_S^T (1−G0)]^{-1}`. Woodbury gives

    G_k = G0 − G0(:,S) Γ̃^{-1} D (1 − G0)(S,:),      Γ̃ = I_k + D (I − G0)_{SS}

and `det(1+B_k)/det(1+B_0) = det Γ̃`. (Γ̃ = D Γ with Γ = D^{-1} + (I−G0)_SS; we store Γ̃^{-1},
which avoids dividing by Δ.)

* Proposal at `p ∉ S`: `(1−G0)(S,p) = −G0(S,p)`, so

      a = Γ̃^{-1} (D G0(S,p))                   (O(k²))
      G_k(p,p) = G0(p,p) + G0(p,S) · a

  and `d = 1 + Δ(1 − G_k(p,p))` is then used exactly as in Sec. 1.
* Accept: border Γ̃ with `w = −D G0(S,p)`, `q = −Δ G0(p,S)`, `σ = 1 + Δ(1 − G0(p,p))`.
  Schur complement `s = σ − q Γ̃^{-1} w = d`: **the Schur complement is the
  fermionic ratio factor**, so the ratio and the inverse update share their arithmetic. Then

      Γ̃_{k+1}^{-1} = [ Γ̃^{-1} + (Γ̃^{-1}w)(qΓ̃^{-1})/s ,  −Γ̃^{-1}w/s ]
                     [ −(qΓ̃^{-1})/s                   ,   1/s       ]      (O(k²), ZGERU)

  with `Γ̃^{-1}w = −a` (already computed for the ratio) and `qΓ̃^{-1}` one extra ZGEMV.
* Flush: `M = Γ̃^{-1} D (E_S^T − G0(S,:))` (k×N, ZGEMM), `G0 ← G0 − G0(:,S) M` (ZGEMM, O(N² k)),
  using a copy of `G0(:,S)`.
* **Wraps between vertices (diagonal operators only).** Let a wrap be `G ← T G T^{-1}`, with T
  diagonal and non-trivial only on P. Read G0 entries lazily (always from the current G0) and
  apply the wrap to G0 as ALF already does. Then the identity above survives if
  `Γ̃ → T_S Γ̃ T_S^{-1}` (D commutes with T_S), i.e.

      for each s_m ∈ S ∩ P with factor t:  Γ̃^{-1}(m,:) *= t,  Γ̃^{-1}(:,m) /= t     (O(k))

  * Up-sweep: the wrap at `n` precedes the update at `n`. For one vertex per site per slice
    (Hubbard), wraps never touch S, so no rescaling occurs.
  * Down-sweep: the un-wrap with the **new** field follows the update at `n`. It always touches
    the site just added, so the rescaling is required. The numerical check confirms this; without it the down-sweep fails.
* Restrictions and forced flushes:
  1. **Non-diagonal operators** (`U ≠ 1`) mix S and non-S rows. Flush before any wrap of a
     non-diagonal vertex. In the proposed design, non-diagonal vertices are not submatrix-eligible.
  2. **Same site twice in one block** (`p ∈ S`): `(1+Δ_2 ee^T)(1+Δ_1 ee^T) ≠ 1 + (Δ_1+Δ_2)ee^T`,
     and `(1−G0)(S,p) ≠ −G0(S,p)`. Flush before proposing at `p ∈ S`.
* Cost per slice: ratio O(M k̄²) (level 2, cache-resident k×k); flush O(N² k) + O(N k²) (level 3).
  [SX-S] Table I: overhead `O(β n_d² N)` versus delayed `O(β n_d N²)`.

## 5. Multiple flavours, N_SUN, phases

* X, Y (delayed) and S, D, Γ̃^{-1} (submatrix) are kept per computed flavour. Acceptance is joint:
  on accept, every computed flavour with `N_non_zero = 1` is extended. A flavour with
  `N_non_zero = 0` has `d_f = 1` and no update, as in `Upgrade2`.
* N_SUN enters only through the scalar power `(Π_f R_f)^{N_SUN}`. The colour copies share G.
* The phase is a scalar update; `Op_phase` and the recomputation at stabilization (`main.F90`) are
  untouched. Complex `g` (density-channel HS, `Mz=.F.`) and `e^{Z1 α}` enter only through the
  same scalar expressions.

## 6. Numerical check

Script: N=24, complex g, 4-valued φ, random distinct site order, fixed accept mask,
up- and down-sweep. It compares standard rank-1, delayed, submatrix and the exact `(1+A)^{-1}`:

    up   nd=1000 delayed  : max|dR/R|=1.65e-15 max|dG|/max|G|=3.93e-15
    up   nd=1000 submatrix: max|dR/R|=1.06e-14 max|dG|/max|G|=1.11e-14
    up   nd=   5 delayed  : max|dR/R|=1.70e-15 max|dG|/max|G|=1.48e-15
    up   nd=   5 submatrix: max|dR/R|=2.53e-15 max|dG|/max|G|=3.21e-15
    down nd=1000 delayed  : max|dR/R|=2.14e-14 max|dG|/max|G|=9.43e-15
    down nd=1000 submatrix: max|dR/R|=2.28e-14 max|dG|/max|G|=1.45e-14
    down nd=   5 delayed  : max|dR/R|=1.53e-14 max|dG|/max|G|=7.09e-15
    down nd=   5 submatrix: max|dR/R|=2.69e-14 max|dG|/max|G|=9.10e-15
    standard vs exact: 6.2e-15 (up), 4.5e-15 (down)

The script also asserts `Δ·s = R` for the Γ-form Schur complement at every accepted step.
The script lives in `docs/scripts/check_formulas.py`.

## 7. Comparison with the independent derivation

A subagent re-derived everything from the papers and the ALF source without seeing these notes
(`docs/independent-derivation.md`, its own numpy check `docs/scripts/check_updates_independent.py`
with output `check_updates_independent.txt`: N=6, 3 slices, complex g, ALF-style type-1/2 wraps, compared with the exact inverse).

| Item | These notes | Independent | Status |
|---|---|---|---|
| Ratio, weight, phase | §1 | same (`R_f = 1+Δ(1−G_qq)`, `^N_SUN`, `e^{Z1α}`, γ, `\|Re(Phase R)/Re(Phase)\|`) | agree |
| Rank-1 update | `G − G(:,p)(Δ/d)(e_p − G(p,:))` | same | agree |
| Delayed | `X ← [X,(Δ/d)c]`, `Y ← [Y, e_p − r]` | `x = c`, `y = (Δ/R)(e_p − r)`; the scalar sits in the other factor | equivalent |
| Wrap on pending factors | `X ← LX`, `Y ← R^T Y` | `X ← RX`, `Y ← R^{-T}Y` (same, their R is a similarity) | agree |
| Γ | `Γ̃ = I + D(I−G0)_SS` | `Γ = D^{-1} + (I−G0)_SS` = [SX-S] Eq. (29) | `Γ̃ = DΓ`, equivalent |
| Schur complement | `s̃ = d` (ratio factor) | `s = R/Δ` | equivalent |
| Final update | `G0 − G0(:,S)Γ̃^{-1}D(1−G0)(S,:)` | `G0 − G0(:,S)Γ^{-1}(E_S^T − G0(S,:))` | identical |
| Frame rule | rescale `Γ̃^{-1}` row/col of `s ∈ S∩P` | same; needed in the down-sweep (error 0.16 without it) | agree |
| Same site twice | flush (Sun–Xu bordering invalid) | Sun–Xu Γ fails (error 1.7); a corrected border using `(1−G0)(S,q)` with the δ-term works | agree it is a hazard; **we flush** (simpler, never occurs in Hubbard) |
| Non-diagonal ops | not eligible for submatrix | generalised-vector form works but loses the speed-up; use delayed | agree |
| Δ = 0 | `Γ̃` form has no `1/Δ` | must exclude Δ=0 from `Γ` | `Γ̃` form avoids the issue |
| Transposes | `ZGERU`/`'T'` | plain, never conjugate | agree |

**[SX-S] reconciliation** (from the PDF; the text extraction had garbled it): Eq. (29)
`Γ = blockdiag(I + D^{(m)-1}) − P G0 P^T`, Eq. (30) `E = I − P blockdiag(D(I+D)^{-1}) P^T`, Eq. (28)
`G = (G0 + G0 P Γ^{-1} P G0) E`. This equals our final update via
`G0(:,S)(1 + Γ^{-1} G0_SS) = G0(:,S) Γ^{-1}(1 + D^{-1})`. [SX-S] apply the whole slice `e^{V}` at once. ALF's per-vertex G
differs by a diagonal similarity (diagonal ops), which leaves the ratios unchanged.
**[G11]** works with `A = G^{-1}` with column updates. ALF's G corresponds to their `G^T` with `γ ↔ Δ`; their
`Γ_pq = G0(p,q) − δ_pq(1+γ)/γ = −Γ_SX^T`, and ratio `−βγ = Δ s = R`.

No unresolved mismatches.

## 8. Open points / uncertain

* None in the formulas. Implementation choices (diagonal-only eligibility, flush on a repeated site) are in the design spec.
* Recommended block sizes from the papers: [SX-D] `k·n_d = 2^λ`, `λ = max(6, [log2(N/20)])`.
  [SX-S] `λ = [log2(N/12)]`. These are tuned for their platforms; we scan in Phase 5.
