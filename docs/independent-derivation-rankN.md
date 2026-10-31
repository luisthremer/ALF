# Rank-r fast, delayed and submatrix updates in ALF conventions (independent derivation)

Sources used: `Prog/upgrade_mod.F90` (Upgrade2), `Prog/Operator_mod.F90` (Op_set, Op_Wrapup, Op_Wrapdo),
`Libraries/Modules/Mat_subroutines_mod.F90` (ZSLGEMM), `Prog/Wrapgr_mod.F90` (WRAPGRUP/WRAPGRDO), `Prog/Hop_mod.F90`,
arXiv:2308.12005 (Sun & Xu, delay update), arXiv:2404.09989 (Sun & Xu, submatrix), arXiv:1010.3690 (Gull et al.).
arXiv:2602.03656 was not consulted. `Lowrank_update_mod.F90` was deliberately not read.

Notation: everything is per flavour. G = (1+A)^{-1}, G_xy = <c_x c_y^†>, where A is the cyclic product of
time-slice matrices cut at the current position. `^T` is the plain transpose, `^†` the conjugate transpose,
`conj` is elementwise conjugation. For an index list P = (P_1..P_n) let E_P be the N×n matrix with columns
e_{P_i}. So G_{:,P} = G E_P, G_{P,:} = E_P^T G and G_{PP} = E_P^T G E_P. Write 𝟙 for the identity.

## 0. What ALF does (read off the source)

**ZSLGEMM(side, op, n, ., ., A, P, Mat)**, generic path and the n ≤ 8 special cases:
`side='l'`: Mat(P,:) ← op(A)·Mat(P,:); `side='r'`: Mat(:,P) ← Mat(:,P)·op(A); op ∈ {N: A, T: A^T, C: A^†}.

**Op_set**: O (Hermitian, N_op×N_op on P) = U diag(E) U^†. For a non-diagonal operator the non-zero eigenvalues
come first (r = N_non_zero) and the zero ones are filled from the back. Column 1 of U is divided by det U, which
is a phase for unitary U, so U stays unitary and U^† = U^{-1}. For a diagonal operator, U = 𝟙, E = diag O and r = N_op.
E_exp(n, ±s) = exp(±g φ_s E_n) for n ≤ r and 1 otherwise.

**Wraps** (block on P, identity elsewhere; Ĺ means the N×N embedding), with e^{±} ≡ diag(exp(±g φ E)):

| wrap | ALF calls | L (rows P) | R (cols P) |
|---|---|---|---|
| type 1, up (old φ) | 'r','n' with U e^{-}, then 'l','T' with conj(U) e^{+} | e^{+} U^† | U e^{-} |
| type 2, up | 'l','n' U ; 'r','c' U | U | U^† |
| type 2, down | 'r','n' U ; 'l','c' U | U^† | U |
| type 1, down (new φ) | 'l','n' U e^{-} ; 'r','T' conj(U) e^{+} | U e^{-} | e^{+} U^† |

In all four cases R = L^{-1}, so every wrap is a similarity G ← L̂ G L̂^{-1}, i.e. A ← L̂ A L̂^{-1}. For diagonal
operators type 2 is the identity, and type 1 rescales rows P by e^{±gφE} and columns by the inverse factor.
The hopping wraps are also similarities. Up (`HOP_MOD_mmthr`, `mmthl_m1`): G ← B_T G B_T^{-1}. Down (`mmthl`,
`mmthr_m1`): G ← B_T^{-1} G B_T. Here B_T = Π_nc e^{-Δτ T_nc}.

**Per-vertex sequence.** Up: type-1 wrap with the old field, then the update, then the type-2 wrap. Down: type-2
wrap, then the update, then the type-1 wrap with the (possibly new) field. Write V_n = U e^{gφE} U^† on P.
If A = Y·V_n, the type-1 up wrap gives A ← e^{gφE}U^† Y U. The vertex exponential is now the leftmost diagonal
factor in the eigenframe, so φ→φ' is A ← F A with F = 𝟙 + Σ_{m≤r} Δ_m e_{Q_m} e_{Q_m}^T, Q = (P_1..P_r) and
Δ_m = e^{g(φ'-φ)E_m} − 1. In matrix form, F − 𝟙 = E_Q Δ E_Q^T with Δ = diag(Δ_m). The same holds after the
type-2 down wrap. The numpy check confirms this frame bookkeeping against products built from scratch, to 1e-14.

## 1. Exact rank-r ratio and update

Since 𝟙+B = G^{-1}, 𝟙 + F B = [𝟙 + (F−𝟙)(𝟙−G)] G^{-1}. Hence

- ratio: R = det[𝟙_N + E_QΔE_Q^T(𝟙−G)] = det[𝟙_r + Δ(𝟙−G_QQ)] = det[𝟙_r + (𝟙−G_QQ)Δ] (Sylvester)
- Upgrade2 builds `Mat(n,m) = δ_nm(1+Δ_m) − Δ_m G(P_n,P_m)`, which is 𝟙 + (𝟙−G_QQ)Δ. Then
  Ratio(nf) = det(Mat)·exp(Z1·α), and the total ratio is Π_nf Ratio(nf)^{N_SUN}·γ'/γ. The extra factors do not
  depend on the scheme.

**Update.** G' = G[𝟙 + u v^T]^{-1} with u = E_Q Δ (N×r) and v^T = (𝟙−G)_{Q,:} = E_Q^T − G_{Q,:}. This is ALF's
`u(P_n,n) = Δ_n`, `v(i,n) = δ_{i,P_n} − G(P_n,i)`. Since v^T u = (𝟙−G_QQ)Δ, we have 𝟙 + v^T u = Mat and

  **G' = G − G_{:,Q} Δ Mat^{-1} (E_Q^T − G_{Q,:})**,   with Δ Mat^{-1} = [𝟙 + Δ(𝟙−G_QQ)]^{-1} Δ.

ALF's sequential Sherman–Morrison loop builds x_v, y_v with x_v y_v^T = u Mat^{-1} v^T. This holds with the
plain transpose (`zgemv 'T'`, `ZGEMM 'N','T'`). x_v vanishes outside the rows Q, so G x_v = G_{:,Q} x_v(Q,:)
(`Zarr`, `grarr`). The check asserts this identity at every accept. Cost: ratio O(r³), update O(N² r).

## 2. Delayed update, G = G0 − X Y^T

State: G0 (N×N, what ALF holds in GR) and X, Y (N×k). Store Yt = Y^T (k×N).

- **Proposal**: G_QQ = G0_QQ − X_{Q,:} Yt_{:,Q} costs O(k r²). Then Mat = 𝟙 + (𝟙−G_QQ)Δ and the ratio is det Mat.
- **Accept**: form G_{:,Q} = G0_{:,Q} − X Yt_{:,Q} and G_{Q,:} = G0_{Q,:} − X_{Q,:} Yt, at O(N k r). Append
  x_new = G_{:,Q} Δ Mat^{-1} (N×r) and yt_new = E_Q^T − G_{Q,:} (r×N), that is y_new = E_Q − (G_{Q,:})^T.
  This is equivalent to Eqs. (22)-(24) of 2308.12005 up to the sign convention.
- **Flush** when k ≥ n_d: G0 ← G0 − X Yt, O(N² k).
- **Wrap** G ← L̂ G R̂ applied by ALF to G0: G_new = L̂G0R̂ − (L̂X)(Y^T R̂), so
  **X_{P,:} ← L X_{P,:}, Yt_{:,P} ← Yt_{:,P} R**, equivalently **Y_{P,:} ← R^T Y_{P,:}** (plain transpose).
  This does not need R = L^{-1}. Cost is O(n_P² k). Per wrap:
  - type 1 up: X_P ← e^{+}U^† X_P, Y_P ← e^{-} U^T Y_P
  - type 2 up: X_P ← U X_P, Y_P ← conj(U) Y_P
  - type 2 down: X_P ← U^† X_P, Y_P ← U^T Y_P
  - type 1 down (new φ'): X_P ← U e^{-}(φ') X_P, Y_P ← conj(U) e^{+}(φ') Y_P

  The diagonal-operator versions are the same with U = 𝟙. A negative control in the check shows that using
  R^† instead of R^T fails with an O(1) error.

## 3. Submatrix update (growing k×k matrix, lazy reads of G0)

**State**: G0, an ordered list S of k "pending" sites, and a k×k matrix M with

  **G = G0 − G0_{:,S} M (E_S^T − G0_{S,:})**   (M is stored, not its inverse).

This is the general form of a product of updates supported on S: if 𝟙+B = (𝟙+K)(𝟙+B0) with K = E_S K_S E_S^T,
then M = K_S[𝟙 + (𝟙−G0_SS)K_S]^{-1}. It stays well defined when K_S is singular, as with zero rows from padding or
Δ_m = 0. The inverse form Γ = M^{-1} (Gull; Γ_{ik×ik} of 2404.09989) needs K_S to be invertible.

**Proposal** for a vertex at Q (r sites, which may overlap S): read the lazy blocks G0_{Q,S}, G0_{S,Q} and G0_{QQ}
(O(kr) entries). Let I_{SQ} = E_S^T E_Q be the 0/1 coincidence matrix. Then

  t = M (I_{SQ} − G0_{S,Q}) (k×r),  G_QQ = G0_QQ − G0_{Q,S} t,  Mat = 𝟙_r + (𝟙 − G_QQ)Δ,  ratio = det Mat.

G_QQ is the exact current block, so the ratio equals Upgrade2's by construction. In the check the relative
difference is 7e-15. Cost: O(k² r + k r² + r³).

**Accept.** Let S' = S ∪ (Q∖S) with k' = |S'|. Let J_Q (k'×r) place Q inside S' and M̃ = M padded with zeros to
k'×k'. A repeated site simply maps to its existing position. With G0_{:,Q} = G0_{:,S'}J_Q and
(𝟙−G0)_{Q,:} = J_Q^T (𝟙−G0)_{S',:}:

  G_{:,Q} = G0_{:,S'} a,   a = J_Q − M̃ (J_Q − G0_{S',Q})   (k'×r)
  (𝟙−G)_{Q,:} = b^T (E_{S'}^T − G0_{S',:}),   b^T = J_Q^T + G0_{Q,S'} M̃   (r×k')

Inserting these into the rank-r formula of §1 gives

  **M' = M̃ + a Δ Mat^{-1} b^T**   (cost O(k'² r)).

This one formula covers new sites, repeated sites (Q∩S ≠ ∅) and multi-site rank-r vertices. For Q∩S = ∅ it
reduces to the bordered block inverse of 2404.09989 Eq. (38). In that case M' = [[M + M G0_SQ X G0_QS M,
M G0_SQ X], [X G0_QS M, X]] with X = (Δ^{-1} + 𝟙 − G_QQ)^{-1}, and the Schur complement is Δ^{-1}+𝟙−G_QQ, so
det Mat = det Δ · det(Schur). The bordered form cannot express overlap. In the Γ form overlap gives K'_{S'} =
(𝟙+Δ̃)(𝟙+K̃) − 1, which needs a re-inversion.

**Wrap** G ← L̂ G L̂^{-1} with support P_w (the full Op%P, including zero-eigenvalue sites). ALF applies it to G0.
Because L̂ G0 E_S = G0_new L̂E_S and E_S^T(𝟙−G0)L̂^{-1} = E_S^T L̂^{-1}(𝟙−G0_new):

  G = G0_new − G0_new (L̂E_S) M (E_S^T L̂^{-1}) (𝟙 − G0_new).

This needs R = L^{-1}, which holds for every ALF wrap. The cases:

- **P_w ∩ S = ∅**: L̂E_S = E_S, so M is unchanged. Cost: a membership test.
- **Diagonal L** (type-1 wrap of a diagonal operator), any overlap: M_{i,:} ← L_ii M_{i,:} and M_{:,i} ← M_{:,i} R_ii
  for i ∈ S∩P_w. No padding is needed. Cost O(k |S∩P_w|).
- **P_w ⊆ S, non-diagonal**: L̂E_S = E_S L_S, where L_S is the k×k identity with the block L at the positions of
  P_w in S. Then **M ← L_S M R_S**, i.e. M_{pos,:} ← L M_{pos,:} and M_{:,pos} ← M_{:,pos} R. Cost O(k n_w²).
- **Partial overlap, non-diagonal**: L̂E_S has non-zero rows on P_w∖S, so the state leaves the span of
  G0_{:,S}. Two options:
  1. **Pad**: append the sites P_w∖S to S with zero rows and columns in M. This is exact, because the padded
     sites contribute nothing. Then apply the full-overlap rule. Cost O(k'n_w²), and k grows by |P_w∖S|.
  2. Convert or flush (see §4), then continue with an empty S.

  For r < N_op the type-2 wrap right after an accept always partially overlaps, on the zero-eigenvalue sites
  P_{r+1..N_op}, so padding happens routinely. Skipping it gives an O(1) error in the negative control.
- **Flush**: G = G0 − G0_{:,S} M (E_S^T − G0_{S,:}), cost O(N² k). Do it when k > k_max or at the end.

## 4. Hopping wraps, and conversion from submatrix to delayed form

**Delayed state across a hopping wrap** (no flush), using §2 with dense L̂ and R̂:
- up, G ← B_T G B_T^{-1}: **X ← B_T X, Yt ← Yt B_T^{-1}** (Y ← B_T^{-T} Y)
- down, G ← B_T^{-1} G B_T: **X ← B_T^{-1} X, Yt ← Yt B_T** (Y ← B_T^T Y)

X transforms like the left side of G and Yt (k×N) like the right side. The existing routines (`mmthr` on X,
`mmthl_m1` on Yt; `mmthr_m1` on X, `mmthl` on Yt) apply if they accept non-square operands; I did not check
this. Cost is O(N² k) dense, or O(N k·n_checkerboard) with checkerboard factors.

**Conversion of a submatrix state to delayed form**: X_add = G0_{:,S} M (N×k, O(N k²)) and
Yt_add = E_S^T − G0_{S,:}. Append these to (X, Yt) and reset S = ∅. After that the state carries through a
hopping wrap as above. A submatrix state cannot cross a dense hopping wrap directly: L̂E_S is dense, so the
partial-overlap padding would put every site into S.

**Hybrid used in the check**: a delayed base Gb = G0 − X Yt with the submatrix layer
G = Gb − Gb_{:,S} M (E_S^T − Gb_{S,:}) on top. Lazy reads of Gb cost O(k_d) per entry. Wraps transform the base
by §2 and M by §3. At hopping wraps, or when |S| > k_max, the layer is converted into base columns. The base is
flushed when its column count reaches n_d. The whole state survives a full up+down sweep over 2 slices without
a dense flush.

**Cost summary** (k = pending size, k_d = delayed columns, n_w = wrap support):

| step | standard | delayed | submatrix |
|---|---|---|---|
| ratio | O(r³) | O(k r²) | O(k² r + k r²) |
| accept | O(N² r) | O(N k r) | O(k² r) |
| wrap (vertex) | O(N n_w) on G | + O(n_w² k) on X,Y | + 0 / O(k n_w) / O(k n_w²), plus padding |
| flush | — | O(N² k) | O(N² k), or conversion O(N k²) |
| hop wrap | O(N³) on G | + O(N² k) on X,Y | convert, then as delayed |

## 5. Numerical check

`check_indep.py` uses N = 10, one flavour, 2 time slices, a random complex Hermitian hopping and a random
"rest-of-chain" matrix. There are 10 vertices with complex g and continuous fields:
- diagonal: 1-site and 2-site
- non-diagonal: bonds (r=2); a 3-site rank-2 operator (r=2 < N_op=3); [[1,1],[1,1]] (r=1 < 2); a 4-site full-rank
  operator with unsorted P; a disjoint bond
- the bond (4,5) twice and site 3 three times (repeated sites)

The check runs a full up sweep then a down sweep with the exact ALF wrap sequence and a fixed accept mask
(32 of 40 accepted). After every vertex it compares each scheme's ratio with det(𝟙+FA)/det(𝟙+A) and its G with
(𝟙+A)^{-1}. Output (`python3 -I check_indep.py`):

```
operators (P, diag, r, N_op): [([3], True, 1, 1), ([3, 4], True, 2, 2), ([4, 5], False, 2, 2), ([5, 6, 7], False, 2, 3), ([1, 2], False, 1, 2), ([4, 5], False, 2, 2), ([7, 2, 8, 3], False, 4, 4), ([3], True, 1, 1), ([0, 9], False, 2, 2), ([6, 1], False, 2, 2)]
vertex visits: 40, accepted: 32
frame check (tracked A vs product from scratch): after up 5.26e-15, after down 1.16e-14
scheme                     max|G-G_exact|   max|ratio/exact-1|  final |G-G_scratch|
std (Upgrade2)                   1.37e-14             7.34e-15             7.38e-15
dly(nd=5)                        1.30e-14             6.00e-15             7.92e-15  [flushes=10]
dly(nd=inf)                      1.29e-14             6.89e-15             8.51e-15  [flushes=0]
sub(kmax=6,nd=8)                 1.41e-14             6.89e-15             8.64e-15  [pads=13, conversions=9, base flushes=5]
sub(kmax=inf,nd=inf)             1.29e-14             7.33e-15             8.39e-15  [pads=11, conversions=3, base flushes=0]
sub(flush@hop)                   1.24e-14             6.89e-15             7.88e-15  [pads=11, conversions=3, base flushes=4]
ALL CHECKS PASSED
negative controls (max|G-G_exact| of the broken scheme; must be O(1e-3) or larger):
  Y <- R^dag Y instead of R^T Y      : dly(nd=inf) 1.41e+00
  no padding on partial overlap      : sub(kmax=inf,nd=inf) 5.91e-01
  hop: Y^T not multiplied by B^-1    : dly(nd=inf) 1.45e-01
```

## Uncertain or worth double-checking

1. **Where X, Y are transformed in ALF.** In `WRAPGRUP` and `WRAPGRDO`, `Lowrank_wrap_hook` is called only after
   the **type-1** wrap. Nothing is called after the type-2 wrap: up, after the update; down, before the update.
   The type-2 wrap also changes G0 for non-diagonal operators, so X, Y (or M) must also be transformed by
   (U, U^†) up and (U^†, U) down. Otherwise the hook must fold the previous type-2 transform into the next one.
   I did not read `Lowrank_update_mod.F90`, so I cannot tell which is the case; it should be verified. In the
   down direction the hook is called after the type-1 wrap that uses the new field, which matches §2.
2. Zero eigenvalues in the `g_t_alloc` path: the wrap uses exp(φ g_t E(I)) with the stored E(I), which may be
   ~1e-15 rather than exactly 0. This is harmless for the schemes, because wraps are exact similarities and F
   touches only Q = P(1..r).
3. With padding, S can grow quickly for chains of overlapping non-diagonal bonds. k_max should then count padded
   sites. Alternatively, convert to delayed form on partial overlap, which costs O(N k²) and keeps exactness.
4. Multi-flavour: everything is per flavour, but S, M and the padding differ between flavours when N_non_zero or
   U differ between flavours.
5. I assumed `HOP_MOD_*` and `Op_Wrap*` can be applied to non-square X (N×k) and Yt (k×N) operands. If they
   cannot, explicit small GEMMs are needed.
