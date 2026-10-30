# Delayed and submatrix local updates in ALF: an independent derivation

Sources: Sun & Xu arXiv:2404.09989 (submatrix, Sec. II), Sun & Xu arXiv:2308.12005 (delayed), Gull et al. arXiv:1010.3690 (CT-AUX submatrix), and the ALF source files `Prog/upgrade_mod.F90` (Upgrade2), `Prog/Wrapgr_mod.F90`, `Prog/Operator_mod.F90` and `Documentation/model.tex`. I read the ALF code myself; everything below is in ALF notation. A numpy check is in `check_updates.py` in this directory, and its output is in `check_output.txt`.

Notation: `^T` means plain transpose throughout. Never use a conjugate transpose: ALF uses ZGERU and `'T'`, and G is complex. `e_q` is a unit vector. `S` is the ordered list of indices accepted in the current block. `1` is the identity.

---

## 0. Conventions, checked against the source

* `G = (1 + A)^{-1}` with `A = B(τ,0) B(β,τ)` split at the current operator position. `G_xy = <c_x c_y^†>`. There is one G per flavor `nf`; colour enters only through the power `N_SUN`.
* Order of operators: WRAPGRUP applies `Op_V(1)` first after the hopping, so in code order `B_τ = e^{V_M} … e^{V_1} e^{-Δτ T}`. (`model.tex` Btau.eq writes `∏_{k=1}^{M_V}` without stating the order; the code defines it. This does not matter for commuting diagonal operators.)
* `Op_Wrapup` with N_type=1 and field φ: `G ← D U† G U D^{-1}` on rows and columns P, with `D = diag(e^{g φ E})`. With N_type=2: `G ← U G U†`. Together these give `G ← e^{V} G e^{-V}` with `e^{V} = U D U†`.
* `Op_Wrapdo` with N_type=2: `G ← U† G U`. With N_type=1 and field φ: `G ← U D^{-1} G D U†`.
* After a type-1 wrap (up-sweep) or a type-2 wrapdo (down-sweep), G is in the eigenbasis of operator n. Row and column `P(i)` then belong to eigenvector `i`, and the non-zero eigenvalues come first (`Op_set`). For `Op%diag`, U = 1 and type 2 does nothing.
* At the moment of the proposal, `A = e^{V_n(φ)} M` in the site basis and `A_rot = D(φ) M'` in the rotated frame, where M and M' do not depend on φ_n. Changing `φ → φ'` therefore gives
  `A_rot → (1 + Δ) A_rot`, with `Δ = diag(Δ_m)` on rows `P(m)`, `Δ_m = exp(Z1·E_m) − 1`, and `Z1 = g_loc (φ' − φ)` (`g_loc = g`, or `g_t(nt)`).
  This holds in both sweeps: in the up-sweep the type-1 wrap used the old field; in the down-sweep e^{V_n} is still on the left of A.

Two identities give everything that follows. Using `A G = 1 − G`:

```
1 + (1+Δ)A = [1 + Δ(1−G)] (1+A)        ⇒   det ratio = det[1 + Δ(1−G)]
G' = G [1 + Δ(1−G)]^{-1}
```

---

## 1. Acceptance ratio for a rank-1 vertex (Upgrade2, mode "FINAL")

For one non-zero eigenvalue per flavor, `q = Op_V(n,nf)%P(1)` in the rotated frame:

```
Δ_nf  = exp(Z1_nf · E_nf(1)) − 1,            Z1_nf = g_loc,nf (φ' − φ)
R_nf  = 1 + Δ_nf (1 − G_nf(q,q))             (Upgrade2: Mat(1,1) = myexp − Z·G(P1,P1))
Ratio(nf) = R_nf · exp(Z1_nf · alpha_nf)
```

* For more than one non-zero eigenvalue, Upgrade2 uses `Mat = 1 + (1−G_QQ)Δ` (columns scaled by Δ) and takes its determinant. This equals `det[1 + Δ(1−G)_QQ]`.
* If `N_non_zero = 0` for a flavor, `R_nf = 1`.
* Flavors not computed (`Calc_Fl_map`/`N_FL_eff`) are filled by `ham%weight_reconstruction(Ratio)`.

Total ratio and Metropolis weight:

```
Ratiotot = [ ∏_nf Ratio(nf) ]^{N_SUN} · γ(φ')/γ(φ)        (nsigma_new%Gama / nsigma%Gama(n,nt))
Weight   = S0_ratio · T0_proposal_ratio · | Re(Phase·Ratiotot) / Re(Phase) |
accept  ⇔ Weight > ranf_wrap()
on accept: Phase ← Phase · Ratiotot/|Ratiotot| ;  nsigma%f(n,nt) ← φ'
```

* The exp(Z1·alpha) factor sits inside the N_SUN power. This matches the weight `[∏ e^{gαφ} det]^{N_SUN} ∏γ` and `Op_phase`.
* `T0_proposal`, drawn with ranf before Upgrade2 is called, decides whether a proposal happens at all. A delayed or submatrix version must consume random numbers in the same order to reproduce the chain. Roundoff can still flip rare decisions.
* **Only one element of the current G is needed:** `G_nf(q,q)` per computed flavor, in the frame Upgrade2 sees (rotated for non-diagonal operators). For a diagonal operator this is the plain diagonal element `G(P1,P1)` of the wrapped G.

## 2. Rank-1 G update after acceptance

Sherman–Morrison on `G' = G[1 + Δ e_q (e_q − G_{q,:})^T]^{-1}` gives

```
G' = G − G(:,q) · (Δ/R) · (e_q^T − G(q,:))
```

This matches Upgrade2 term by term:

* `u(P1) = Δ`
* `v = e_P − G(P,:)^T`
* `x_v(P1) = Δ/(1 + v(P1)Δ) = Δ/R`
* `xp = G(:,P1)`
* `ZGERU(−x_v(P1), xp, y_v = v)`

The ZGEMM branch for `size(P) > 1` does the same with `x_v(P(1:r),:)`. For r > 1 the recursion over n builds `Δ(1 + (1−G)_QQ Δ)^{-1}` column by column.

## 3. Delayed update

Keep `G = G0 − X Y^T`, with X and Y of size N×k per flavor. Before any acceptance, k = 0.

Needed element:

```
G(q,q) = G0(q,q) − Σ_j X(q,j) Y(q,j)
```

On acceptance at q, with Δ and R from Sec. 1, append one column to each factor:

```
x_{k+1} = G(:,q)          = G0(:,q) − X · Y(q,:)^T
y_{k+1} = (Δ/R)(e_q − G(q,:)^T) = (Δ/R)(e_q − G0(q,:)^T + Y · X(q,:)^T)
```

This is exactly Sun & Xu 2308.12005 Eqs. (22)–(26) with `U = x`, `S V = −y^T`. For rank r per vertex, append r columns: `x = G(:,Q)` and `y^T = [1 + Δ(1−G)_QQ]^{-1} Δ (E_Q^T − G(Q,:))`.

**Flush** when k reaches n_d, at the end of the slice, and before anything non-local. Non-local means the hopping wrap, measurements (`Obser`), `Wrapgr_Random_update`/`PlaceGR`, and stabilization:

```
G0 ← G0 − X Y^T      (one ZGEMM, N×N×k)
```

**Per-vertex wraps.** Each ALF wrap is a local similarity `G ← R G R^{-1}`, where R differs from 1 only on rows and columns P. The identity stays exact if the wrap is applied to all three stored objects:

```
G0 ← R G0 R^{-1}      (as ALF does now, O(N|P|) for diagonal R)
X  ← R X              (rows P only, O(k|P|²))
Y  ← R^{-T} Y         (rows P only)
```

The four concrete cases:

* Diagonal operator with factor w on row p: `X(p,:) *= w`, `Y(p,:) /= w`.
* Type-1 wrapup, non-diagonal: `X(P,:) ← D U† X(P,:)`, `Y(P,:) ← D^{-1} U^T Y(P,:)`.
* Type-2 wrapup: `X(P,:) ← U X(P,:)`, `Y(P,:) ← conj(U) Y(P,:)`.
* Wrapdo: use the analogous inverse matrices. The type-1 wrapdo uses the **new** field.

Forgetting to transform X and Y (rows P) is a silent error. The delayed scheme is exact for any operator: non-diagonal, overlapping, or repeated sites. This was checked numerically.

## 4. Submatrix update

### 4.1 Reference-frame derivation (diagonal operators, or no wraps in between)

Suppose all accepted updates in the block act in one fixed frame as left multiplications `(1 + e_{p_j} Δ_j e_{p_j}^T)`. A diagonal similarity commutes with a diagonal Δ, so per-vertex diagonal wraps do not change this. Write `S = (p_1..p_k)` and `Δ_S = diag(Δ_j)`, distinct sites for now:

```
G_k = G0 [1 + P_S Δ_S P_S^T (1−G0)]^{-1}
    = G0 − G0(:,S) Γ^{-1} (E_S^T − G0(S,:))                     (Woodbury)
Γ   = Δ_S^{-1} + (1 − G0)_{SS}       i.e.  Γ_ij = δ_ij (1/Δ_i + 1) − G0(p_i,p_j)
det ratio of the whole block = det(Δ_S) det(Γ)
```

Here `E_S^T` is the k×N matrix whose rows are `e_{p_j}^T`. The Δ-weighted alternative `G0(:,S) Δ_S [1 + (1−G0)_SS Δ_S]^{-1}(...)` avoids `1/Δ`.

**Needed element for a new proposal at q:**

```
G_k(q,q) = G0(q,q) − G0(q,S) Γ^{-1} (δ_{S,q} − G0(S,q))
         = G0(q,q) + G0(q,S) Γ^{-1} G0(S,q)           if q ∉ S
R = 1 + Δ_q (1 − G_k(q,q))   ×   all the other factors of Sec. 1
```

**Bordered recursion on acceptance:**

```
Γ_{k+1} = [[Γ_k, b],[c, d]],
b = (1−G0)(S,q) = δ_{S,q} − G0(S,q),   c = −G0(q,S),   d = 1/Δ_q + 1 − G0(q,q)
s = d − c Γ_k^{-1} b                     (Schur complement)
Γ_{k+1}^{-1} = [[Γ^{-1} + Γ^{-1}b c Γ^{-1}/s,  −Γ^{-1}b/s],[−cΓ^{-1}/s,  1/s]]
```

**Schur relation:** `s = 1/Δ_q + 1 − G_k(q,q) = R_q/Δ_q`, so `R_q = Δ_q · s`. Equivalently `det Γ_{k+1}/det Γ_k = s`. In code, compute R directly and use `s = R/Δ` for the border: `G_k(q,q)` already gives `Γ^{-1}b` and `cΓ^{-1}`.

**Final update (flush):** with `Z = G0(:,S) Γ^{-1}` (N×k ZGEMM),

```
G ← G0 − Z (E_S^T − G0(S,:))
```

Copy `G0(:,S)` and `G0(S,:)` before overwriting G0.

### 4.2 Two vertices on the same site in one block

The bordering above silently assumes `q ∉ S`; the Sun–Xu Eqs. (35)–(37) have no δ term. If `q = p_i` for some earlier i, the products `(1+eΔ_1e^T)(1+eΔ_2e^T)` do not add Δ's. I derived the general product of rank-1 factors `∏_{j=k..1}(1 + u_j Δ_j v_j^T) = 1 + U L V^T` with `L = (Δ^{-1} − T)^{-1}`, where T is strictly lower triangular with `T_ij = v_i^T u_j`. This gives

```
Γ = Δ^{-1} − T + V^T (1 − G0) U .
```

For unit vectors this reads `Γ_ij = δ_ij/Δ_i + θ(i ≤ j) δ_{p_i p_j} − G0(p_i,p_j)`. In words: **the new border column uses (1−G0)(S,q), including the δ term; the border row uses −G0(q,S), with no δ.** `G_k(q,q)` must also include the δ term (first line of the formula in 4.1).

* With this Θ term, repeated sites are exact (numerically verified).
* The Sun–Xu form without it is wrong: errors of O(1) in the test.
* Alternatives are a flush before any vertex whose site is already in S, or Gull's "Γ → Γ + uv^T" modification.

### 4.3 ALF per-vertex wraps, diagonal operators

There are two consistent implementations.

**(A) Frozen G0.** Do not wrap G0 inside a block.

* Diagonal elements and Δ are invariant under diagonal similarity, so `G_cur(q,q) = G_ref(q,q)`. All ratios use the frozen G0, and Γ is built from frozen entries.
* Accumulate the diagonal wrap factors `w_i`, exactly the scalings `Op_Wrapup`/`Op_Wrapdo` would have applied: old field in the up-sweep, new field in the down-sweep, times `g_t(nt)` if allocated.
* At flush: `G ← diag(w) G_ref diag(w)^{-1}`.
* No Γ rescaling is ever needed.

**(B) Eager wraps on G0,** as ALF does today. Then Γ must live in the current frame, `Γ_cur = W_S Γ_ref W_S^{-1}`. That is, `Γ_cur = Δ^{-1} + Θ + (−G0_cur)_SS`; the Δ and Θ parts are invariant.

Rule: whenever a wrap multiplies row p of G0 by w (and column p by 1/w), for every index i with `p_i = p`:

```
Γ^{-1}(i,:) *= w ,   Γ^{-1}(:,i) /= w      (the same rule holds for Γ)
```

* **Up-sweep, distinct sites:** the type-1 wrap at vertex n hits `p_n` *before* it can join S. Later vertices hit other sites. **No rescaling is needed** (verified).
* **Down-sweep:** the post-update `Op_Wrapdo` type 1 with the **new** field hits `p_n`, which has just joined S if accepted. **Rescale the last row of Γ^{-1} by `w = e^{−g φ_new E}` and the last column by `1/w`.** If the proposal was rejected, nothing is rescaled. Omitting this gives errors of 0.16 in G and 4% in ratios in the test.
* **Up-sweep with repeated sites:** the wrap at the second vertex on p hits `p ∈ S`, so rescaling is needed there too.

### 4.4 What breaks for non-diagonal operators

* Type-1 and type-2 wraps rotate rows and columns P with `U`/`U†`. Stored site-indexed entries `G0(q,S)`, `G0(S,q)` and Γ are no longer just rescaled.
* The update in the reference frame is `R^{-1}(1+Δ)R`, which is no longer a site-diagonal Δ.
* With overlapping P blocks the factors do not commute either.
* The **generalized** form still works: carry `u_j = R e_q` and `v_j = R^{-T} e_q` through all later wraps (O(k|P|²) per vertex), as X and Y are carried in Sec. 3. Use `Γ = Δ^{-1} − T + V^T(1−G0)U`. This Γ is similarity-invariant, so it needs no rescaling. Use border column `V^T(1−G0)e_q`, border row `−e_q^T G0 U`, `G(q,q) = G0(q,q) − G0(q,:)UΓ^{-1}V^T(1−G0)e_q`, and flush `G = G0 − G0UΓ^{-1}V^T(1−G0)`. I verified this numerically for overlapping bond projectors.
* The supports of u and v grow under overlapping wraps, so most of the submatrix saving is lost. The practical recommendation is: **submatrix for diagonal operators only; delayed or flush otherwise.**

### 4.5 Other pitfalls

* Diagonal operators with `N > 1` get `N_non_zero = N` in Op_set, even when some `E = 0`. Indices with `Δ = 0` must not enter S, because Γ contains `1/Δ`. Use the Δ-weighted form if Δ can be tiny.
* An accepted proposal with `Δ = 0` (for example `g_t(nt) = 0`) should simply be skipped.
* Each flavor needs its own S, Γ^{-1} and G0. A flavor with `Δ_nf = 0` contributes `R_nf = 1` and does not grow its S.
* Flush before hopping, measurements, global-τ updates (`Wrapgr_Random_update`), and stabilization.

## 5. Reconciliation with the papers

**Sun & Xu 2404.09989.** They use the same convention as ALF: G = (I + B B)^{-1}, `e^{V(s')} = (I+Δ)e^{V(s)}`, and ratio `det[I + Δ(I−G)]` in Eq. (7). Their `S = I + 𝒱D` with `𝒱 = (1−G)_xx` is exactly ALF's `Mat`.

* **Eq. (29):** `Γ^{(i)} = blockdiag(I + D^{(m)−1}) − P G^{(0)} P^T`, i.e. `Γ = Δ^{-1} + 1 − G0_SS`. For distinct sites this is identical to my Γ. Their `D^{(m)}` is the k×k Δ of step m; k = 1 for Hubbard.
* **Eq. (30):** `E = I − P blockdiag(D(I+D)^{-1}) P^T` is N×N and acts from the right, dividing columns S by `(1+Δ_j)`.
* **Eq. (28):** `G = (G0 + G0(:,S)Γ^{-1}G0(S,:)) E`. I showed it equals my `G0 − G0(:,S)Γ^{-1}(E_S^T − G0(S,:))`, using `G0(:,S)(1 + Γ^{-1}G0_SS) = G0(:,S)Γ^{-1}(1 + Δ^{-1})`.
* **Eqs. (34)–(37):** `W = −G0(S_old,x)`, `Q = −G0(x,S_old)`, `Σ = I + D^{-1} − G0_xx`. These are my b, c, d for `x ∉ S`.
* **Eq. (38):** the block inverse; the Schur complement is `(Σ − QΓ^{-1}W) = det S/det D`.
* They apply the whole `e^{V(s)}` of a slice at once. ALF's per-vertex G differs from theirs by a diagonal similarity for diagonal operators, so ratios are identical.
* Their bordering is valid only for distinct sites (see 4.2).

**Sun & Xu 2308.12005.** Eqs. (22)–(27) are my Sec. 3: `G^{(i)} = G^{(0)} + Σ U S V`.

**Gull et al. 1010.3690.** Their Eq. (16) is a column update of `A = G^{-1}`. ALF's update is a row update of `G^{-1}`: `G'^{-1}_{pj} = G^{-1}_{pj} + Δ(G^{-1}_{pj} − δ_pj)`. Hence Gull's G corresponds to ALF's `G^T` with `γ ↔ Δ`.

* Their Eq. (25) `Γ_pq = G0(p,q) − δ_pq(1+γ)/γ` equals `−Γ_SunXu^T`.
* Their row factor `D^{-1} = 1/(1+γ)` in Eqs. (24) and (26) becomes Sun–Xu's column factor E.
* Their ratio `−βγ` in Eqs. (36) and (39), with β the Schur complement of Γ_Gull (equal to −s), becomes `Δ·s = R`.
* Their remark that changing the same spin twice needs `Γ → Γ + uv^T` is consistent with 4.2.

## 6. Numerical check

`check_updates.py` uses N=6, L=3, complex g, and emulates ALF's type-1/type-2 wraps with old and new fields in both sweeps. It compares each method against sequential Upgrade2-like rank-1 updates and against an exact recomputation of `(1+B(τ,0)B(β,τ))^{-1}` after every slice. All methods replay the same accept/reject decisions.

| Case | Correct methods | max\|G − G_exact\| | Failing variants |
|---|---|---|---|
| Diagonal operators, distinct sites | delayed, subA, subB, subGen | ≤ 6e-15 | subB without rescaling: 0.16 (fine after the up-sweep alone) |
| Diagonal operators, repeated sites | delayed, subA, subB (with Θ), subGen | ≤ 1.3e-14 | Sun–Xu Γ without Θ: 1.7; subB without rescaling: 0.77 |
| Non-diagonal overlapping bond projectors | delayed, subGen | ≤ 8e-15 | — |

Across all three cases:

* The exact determinant ratio matches the sequential `R` to 2.5e-15.
* All correct methods reproduce the sequential ratios to 1e-15.
* `|Δ·s − R|/|R|` ≤ 3.5e-16.

UNCERTAIN, not covered here: performance and stability over long blocks (the conditioning of Γ for large n_d), and the global-τ (`mode = "INTERMEDIATE"`) path. My recommendation is to flush before that path rather than mix it in.
