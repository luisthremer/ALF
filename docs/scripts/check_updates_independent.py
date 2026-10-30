"""
Independent numerical check of delayed / submatrix update formulas in ALF conventions.

G = (1 + B(tau,0) B(beta,tau))^{-1};  B_tau = e^{V_M} ... e^{V_1} e^{-dtau T}.
Up-sweep  (WRAPGRUP): hop;  for n=1..M: type-1 wrap (old field), propose/update, type-2 wrap.
Down-sweep(WRAPGRDO): for n=M..1: type-2 wrapdo, propose/update, type-1 wrapdo (NEW field); hop.
Proposal at eigen-index q = P[0] in the rotated frame, Delta = exp(g (phi'-phi) E[0]) - 1,
rank-1 per vertex (one non-zero eigenvalue).

Methods compared (same accept/reject decisions replayed):
  seq      : Upgrade2-like rank-1 update of the full G
  delayed  : G = G0 - X Y^T, wraps applied to G0, X (rows), Y (rows, inverse-transpose)
  subA     : submatrix, G0 frozen at block start, diag wraps accumulated in w, applied at flush (diag ops only)
  subB     : submatrix, G0 wrapped eagerly, Gamma^{-1} rescaled when a wrap hits a site in S (diag ops only)
  subBnaive: as subB but WITHOUT rescaling (expected to fail in down-sweep)
  subXu    : as subB but Sun-Xu Gamma without the repeated-site Theta term (fails if sites repeat in a block)
  subGen   : generalized submatrix with vectors u_j, v_j carried along by all wraps (works for any op)
"""
import numpy as np

rng = np.random.default_rng(12345)
SCHUR_DEV = [0.0]


def herm_exp(H, c):
    E, U = np.linalg.eigh(H)
    return (U * np.exp(c * E)) @ U.conj().T


class Op:
    def __init__(self, P, O, g):
        self.P = np.array(P)
        self.O = np.array(O, dtype=complex)
        self.g = g
        E, U = np.linalg.eigh(self.O)
        # non-zero eigenvalues first (as in ALF Op_set); we only use rank-1 ops: exactly one non-zero
        order = np.argsort(np.abs(E) < 1e-12, kind="stable")
        self.E = E[order]
        self.U = U[:, order]
        self.nnz = int(np.sum(np.abs(self.E) > 1e-12))
        assert self.nnz == 1
        self.diag = np.allclose(self.O, np.diag(np.diag(self.O)))
        if self.diag:
            self.E = np.real(np.diag(self.O)).copy()
            self.U = np.eye(len(P), dtype=complex)


def embed(N, P, blk):
    M = np.eye(N, dtype=complex)
    M[np.ix_(P, P)] = blk
    return M


def expV(N, op, phi):
    return embed(N, op.P, op.U @ np.diag(np.exp(op.g * phi * op.E)) @ op.U.conj().T)


class Model:
    def __init__(self, N, L, dtau, ops):
        self.N, self.L, self.ops = N, L, ops
        T = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)) * 0.3
        T = T + T.conj().T
        self.Bhop = herm_exp(T, -dtau)
        self.Bhop_inv = herm_exp(T, dtau)

    def Bslice(self, phis):
        B = self.Bhop.copy()
        for n, op in enumerate(self.ops):
            B = expV(self.N, op, phis[n]) @ B
        return B

    def G_exact(self, fields, tau):
        N = self.N
        left = np.eye(N, dtype=complex)
        for t in range(tau):
            left = self.Bslice(fields[t]) @ left
        right = np.eye(N, dtype=complex)
        for t in range(tau, self.L):
            right = self.Bslice(fields[t]) @ right
        return np.linalg.inv(np.eye(N) + left @ right)

    def detw(self, fields):
        B = np.eye(self.N, dtype=complex)
        for t in range(self.L):
            B = self.Bslice(fields[t]) @ B
        return np.linalg.det(np.eye(self.N) + B)


# ---------------------------------------------------------------- wrap matrices (R: G <- R G R^{-1})
def R_up1(N, op, phi):   # Op_Wrapup type 1: e^{g phi E} U^dag  on block P
    return embed(N, op.P, np.diag(np.exp(op.g * phi * op.E)) @ op.U.conj().T)


def R_up2(N, op):        # Op_Wrapup type 2: U
    return embed(N, op.P, op.U)


def R_do2(N, op):        # Op_Wrapdo type 2: U^dag
    return embed(N, op.P, op.U.conj().T)


def R_do1(N, op, phi):   # Op_Wrapdo type 1: U e^{-g phi E}
    return embed(N, op.P, op.U @ np.diag(np.exp(-op.g * phi * op.E)))


def sim(R, G):
    return R @ G @ np.linalg.inv(R)


# ---------------------------------------------------------------- engines
class Seq:
    def __init__(s, G): s.G = G.copy()
    def wrap(s, R): s.G = sim(R, s.G)
    def Gqq(s, q): return s.G[q, q]
    def accept(s, q, D, Rat):
        N = s.G.shape[0]
        e = np.zeros(N); e[q] = 1
        s.G = s.G - np.outer(s.G[:, q], (D / Rat) * (e - s.G[q, :]))
    def flush(s): pass
    def full(s): return s.G


class Delayed:
    def __init__(s, G, nd): s.G0 = G.copy(); s.X = []; s.Y = []; s.nd = nd
    def wrap(s, R):
        Ri = np.linalg.inv(R)
        s.G0 = R @ s.G0 @ Ri
        s.X = [R @ x for x in s.X]
        s.Y = [Ri.T @ y for y in s.Y]
    def _col(s, q):
        c = s.G0[:, q].copy()
        for x, y in zip(s.X, s.Y): c -= x * y[q]
        return c
    def _row(s, q):
        r = s.G0[q, :].copy()
        for x, y in zip(s.X, s.Y): r -= x[q] * y
        return r
    def Gqq(s, q): return s._col(q)[q]
    def accept(s, q, D, Rat):
        N = s.G0.shape[0]
        e = np.zeros(N); e[q] = 1
        s.X.append(s._col(q)); s.Y.append((D / Rat) * (e - s._row(q)))
        if len(s.X) == s.nd: s.flush()
    def flush(s):
        for x, y in zip(s.X, s.Y): s.G0 = s.G0 - np.outer(x, y)
        s.X, s.Y = [], []
    def full(s):
        G = s.G0.copy()
        for x, y in zip(s.X, s.Y): G -= np.outer(x, y)
        return G


def border(Ginv, col, row, d):
    """inverse of [[Gamma, col],[row, d]] given Ginv = Gamma^{-1}; returns (new inverse, Schur complement)"""
    k = Ginv.shape[0]
    if k == 0:
        return np.array([[1.0 / d]], dtype=complex), d
    a = Ginv @ col
    b = row @ Ginv
    sch = d - row @ a
    new = np.zeros((k + 1, k + 1), dtype=complex)
    new[:k, :k] = Ginv + np.outer(a, b) / sch
    new[:k, k] = -a / sch
    new[k, :k] = -b / sch
    new[k, k] = 1.0 / sch
    return new, sch


class SubB:
    """G0 wrapped eagerly; diag ops; Gamma^{-1} in current frame (rescaled if rescale=True)."""
    def __init__(s, G, nd, rescale=True, theta=True):
        s.G0 = G.copy(); s.S = []; s.D = []; s.Ginv = np.zeros((0, 0), complex)
        s.nd, s.rescale, s.theta = nd, rescale, theta
    def wrap(s, R):
        if not np.allclose(R, np.diag(np.diag(R))):   # hopping: only allowed right after a flush
            assert not s.S
            s.G0 = sim(R, s.G0); return
        w = np.diag(R)
        s.G0 = (w[:, None] * s.G0) / w[None, :]
        if s.rescale and s.S:
            wS = w[s.S]
            s.Ginv = (wS[:, None] * s.Ginv) / wS[None, :]
    def _vecs(s, q):
        S = s.S
        colG = s.G0[S, q]                       # G0_{S,q}
        rowG = s.G0[q, S]                       # G0_{q,S}
        delta = np.array([1.0 if p == q else 0.0 for p in S])
        return colG, rowG, delta
    def Gqq(s, q):
        if not s.S: return s.G0[q, q]
        colG, rowG, delta = s._vecs(q)
        # G_cur(q,q) = G0_qq - G0_{q,S} Gamma^{-1} (delta_{S,q} - G0_{S,q})
        return s.G0[q, q] - rowG @ s.Ginv @ (delta - colG)
    def accept(s, q, D, Rat):
        colG, rowG, delta = s._vecs(q) if s.S else (np.zeros(0), np.zeros(0), np.zeros(0))
        col = (delta if s.theta else 0 * delta) - colG    # (1-G0)_{S,q} incl. Theta term
        row = -rowG                                           # -G0_{q,S}
        d = 1.0 / D + 1.0 - s.G0[q, q]
        s.Ginv, sch = border(s.Ginv, col, row, d)
        if s.rescale and s.theta: SCHUR_DEV[0] = max(SCHUR_DEV[0], abs(sch * D - Rat) / abs(Rat))
        s.S.append(q); s.D.append(D)
        if len(s.S) == s.nd: s.flush()
    def flush(s):
        if not s.S: return
        N = s.G0.shape[0]; S = s.S
        ES = np.zeros((len(S), N)); ES[np.arange(len(S)), S] = 1
        s.G0 = s.G0 - s.G0[:, S] @ s.Ginv @ (ES - s.G0[S, :])
        s.S, s.D, s.Ginv = [], [], np.zeros((0, 0), complex)
    def full(s):
        if not s.S: return s.G0.copy()
        N = s.G0.shape[0]; S = s.S
        ES = np.zeros((len(S), N)); ES[np.arange(len(S)), S] = 1
        return s.G0 - s.G0[:, S] @ s.Ginv @ (ES - s.G0[S, :])


class SubA(SubB):
    """G0 frozen at block start, diag wraps accumulated in w and applied at flush."""
    def __init__(s, G, nd):
        super().__init__(G, nd)
        s.w = np.ones(G.shape[0], complex)
    def wrap(s, R):
        if not np.allclose(R, np.diag(np.diag(R))):
            assert not s.S and np.allclose(s.w, 1)
            s.G0 = sim(R, s.G0); return
        s.w = s.w * np.diag(R)
    def flush(s):
        SubB.flush(s)
        s.G0 = (s.w[:, None] * s.G0) / s.w[None, :]
        s.w = np.ones_like(s.w)
    def full(s):
        G = SubB.full(s)
        return (s.w[:, None] * G) / s.w[None, :]


class SubGen:
    """Generalized submatrix: vectors u_j, v_j (current frame), Gamma = Delta^{-1} - T + V^T (1-G0) U."""
    def __init__(s, G, nd):
        s.G0 = G.copy(); s.U = []; s.V = []; s.Ginv = np.zeros((0, 0), complex); s.nd = nd
    def wrap(s, R):
        Ri = np.linalg.inv(R)
        s.G0 = R @ s.G0 @ Ri
        s.U = [R @ u for u in s.U]
        s.V = [Ri.T @ v for v in s.V]
    def _UV(s):
        return np.array(s.U).T, np.array(s.V).T   # N x k
    def Gqq(s, q):
        if not s.U: return s.G0[q, q]
        U, V = s._UV()
        N = s.G0.shape[0]; e = np.zeros(N); e[q] = 1
        return s.G0[q, q] - s.G0[q, :] @ U @ s.Ginv @ (V.T @ (e - s.G0[:, q]))
    def accept(s, q, D, Rat):
        N = s.G0.shape[0]; e = np.zeros(N, complex); e[q] = 1
        if s.U:
            U, V = s._UV()
            col = V.T @ (e - s.G0 @ e)          # v_i^T (1-G0) u_new   (i < new: keep V^T U term)
            row = -(e @ s.G0 @ U)               # -v_new^T G0 u_j      (j < new: V^T U term cancelled by -T)
        else:
            col = row = np.zeros(0)
        d = 1.0 / D + 1.0 - s.G0[q, q]
        s.Ginv, sch = border(s.Ginv, col, row, d)
        SCHUR_DEV[0] = max(SCHUR_DEV[0], abs(sch * D - Rat) / abs(Rat))
        s.U.append(e.copy()); s.V.append(e.copy())
        if len(s.U) == s.nd: s.flush()
    def flush(s):
        s.G0 = s.full(); s.U, s.V = [], []; s.Ginv = np.zeros((0, 0), complex)
    def full(s):
        if not s.U: return s.G0.copy()
        U, V = s._UV()
        N = s.G0.shape[0]
        return s.G0 - s.G0 @ U @ s.Ginv @ (V.T @ (np.eye(N) - s.G0))


# ---------------------------------------------------------------- driver
def run(model, fields0, engines, nd, check_exact_ratio=True, label=""):
    N, L, ops = model.N, model.L, model.ops
    fields = [f.copy() for f in fields0]
    G = model.G_exact(fields, 0)
    eng = {k: f(G) for k, f in engines.items()}
    maxdev = {k: 0.0 for k in eng}
    ratdev = {k: 0.0 for k in eng}
    exdev = 0.0

    def all_wrap(R):
        for e in eng.values(): e.wrap(R)

    def propose(t, n, op):
        nonlocal exdev
        q = op.P[0]
        phi, phin = fields[t][n], -fields[t][n]
        D = np.exp(op.g * (phin - phi) * op.E[0]) - 1.0
        rats = {k: 1.0 + D * (1.0 - e.Gqq(q)) for k, e in eng.items()}
        R0 = rats["seq"]
        if check_exact_ratio:
            f2 = [f.copy() for f in fields]; f2[t][n] = phin
            Rex = model.detw(f2) / model.detw(fields)
            exdev = max(exdev, abs(Rex - R0) / abs(R0))
        for k in eng: ratdev[k] = max(ratdev[k], abs(rats[k] - R0) / abs(R0))
        acc = rng.random() < min(1.0, abs(R0)) * 0.9 + 0.1   # accept often; decisions shared
        if acc:
            for k, e in eng.items(): e.accept(q, D, rats[k] if k != "seq" else R0)
            fields[t][n] = phin
        return acc

    def compare(Gex):
        for k, e in eng.items():
            maxdev[k] = max(maxdev[k], np.max(np.abs(e.full() - Gex)))

    # up-sweep
    for t in range(L):
        for e in eng.values(): e.flush()
        all_wrap(model.Bhop)
        for n, op in enumerate(ops):
            all_wrap(R_up1(N, op, fields[t][n]))
            propose(t, n, op)
            all_wrap(R_up2(N, op))
        for e in eng.values(): e.flush()
        compare(model.G_exact(fields, t + 1))
    print('   after up-sweep only: ' + ', '.join(f'{k}={v:.1e}' for k, v in maxdev.items()))
    # down-sweep
    for t in reversed(range(L)):
        for e in eng.values(): e.flush()
        for n in reversed(range(len(ops))):
            op = ops[n]
            all_wrap(R_do2(N, op))
            propose(t, n, op)
            all_wrap(R_do1(N, op, fields[t][n]))
        for e in eng.values(): e.flush()
        all_wrap(model.Bhop_inv)
        compare(model.G_exact(fields, t))
    print(f"--- {label}  (nd={nd})   exact-vs-seq ratio rel.dev = {exdev:.2e}")
    for k in eng:
        print(f"   {k:10s}  max|G-G_exact| = {maxdev[k]:.2e}   max rel ratio dev vs seq = {ratdev[k]:.2e}")
    print(f"   max |Delta*Schur - R|/|R| over correct submatrix engines = {SCHUR_DEV[0]:.2e}"); SCHUR_DEV[0] = 0.0


def make_fields(L, M):
    return [rng.choice([-1.0, 1.0], size=M) for _ in range(L)]


N, L, dtau = 6, 3, 0.1
g = 0.35 + 0.2j

# 1) diagonal rank-1 ops on distinct sites
ops = [Op([i], [[1.0]], g) for i in range(N)]
m = Model(N, L, dtau, ops)
nd = 4
eng = {"seq": Seq, "delayed": lambda G: Delayed(G, nd), "subA": lambda G: SubA(G, nd),
       "subB": lambda G: SubB(G, nd), "subBnaive": lambda G: SubB(G, nd, rescale=False),
       "subXu": lambda G: SubB(G, nd, theta=False), "subGen": lambda G: SubGen(G, nd)}
run(m, make_fields(L, len(ops)), eng, nd, label="diag ops, distinct sites")

# 2) diagonal ops, repeated sites within one slice (two vertices on sites 1 and 3), different E
ops = [Op([i], [[1.0]], g) for i in range(N)] + [Op([1], [[-0.7]], g), Op([3], [[1.3]], 0.5)]
m = Model(N, L, dtau, ops)
nd = 8
eng = {"seq": Seq, "delayed": lambda G: Delayed(G, nd), "subA": lambda G: SubA(G, nd),
       "subB": lambda G: SubB(G, nd), "subBnaive": lambda G: SubB(G, nd, rescale=False),
       "subXu": lambda G: SubB(G, nd, theta=False), "subGen": lambda G: SubGen(G, nd)}
run(m, make_fields(L, len(ops)), eng, nd, label="diag ops, repeated sites in block")

# 3) non-diagonal rank-1 ops (bond projectors), overlapping bonds, plus diagonal ones
proj = np.array([[0.5, 0.5], [0.5, 0.5]])
ops = [Op([0, 1], proj, g), Op([1, 2], proj, g), Op([3], [[1.0]], g), Op([2, 4], proj * 1.0, 0.4), Op([4, 5], proj, g)]
m = Model(N, L, dtau, ops)
nd = 5
eng = {"seq": Seq, "delayed": lambda G: Delayed(G, nd), "subGen": lambda G: SubGen(G, nd)}
run(m, make_fields(L, len(ops)), eng, nd, label="non-diag rank-1 ops (overlapping bonds)")
