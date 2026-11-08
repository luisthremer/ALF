#!/usr/bin/env python3
"""Autonomous GPU pipeline: correlation ratio R_chi^AFM(U, L) of the half-filled pi-flux Hubbard model (finite-size crossing,
Uc, nu, eta_phi), including the equilibration study that fixes Nskip.

  rchi_pipeline.py run      [--work DIR] [--bin ALF.out] [--quick] [--deadline-hours H]   whole pipeline (resumable)
  rchi_pipeline.py analyze  [--work DIR] [--quick]                                         analysis + plots only

Stages (each skipped when its result file exists):
  1. warm   equilibration study: NSweep=1 runs at a few (L,U), plot warmup.png, nskip.json
  2. prod   production rounds: runs of NSweep=5 sweeps per bin, continued from confout_0 round after round; after every
            U point (all L) the plots are refreshed; stops at the deadline or after MAX_ROUNDS
  3. final  analysis.png, summary.md, summary.json

Model (Hubbard_smod): square lattice L x L, pi flux (N_Phi = L^2/2), Mz = .F. (SU(2) symmetric charge decoupling),
N_FL = 1, N_SUN = 2, Beta = L, Dtau = 0.1, Symm, Checkerboard, Sequential_scheme = "submatrix" on the GPU.
R_chi = 1 - S(Q+dq)/S(Q), S = SpinZ_eq (k space), Q = (pi,pi), dq = (2pi/L,0),(0,2pi/L) (+- , averaged).
"""
import os, sys, re, json, time, math, shutil, subprocess, argparse, gzip
import numpy as np

# ------------------------------------------------------------------ settings
LS_FULL = [12, 16, 20, 24, 28, 32, 36]
US_FULL = [5.4, 5.3, 5.5, 5.2, 5.6, 5.35, 5.45, 5.0, 6.0, 5.1, 5.8, 5.25, 5.55, 5.7]       # priority order
LS_QUICK = [4, 6, 8, 10]
US_QUICK = [5.2, 5.4, 5.6, 5.8]
WARM_FULL = {12: ([5.0, 5.4, 6.0], 400), 24: ([5.0, 5.4, 6.0], 300), 36: ([5.4, 6.0], 150)}   # L: (U list, number of 1-sweep bins)
WARM_QUICK = {4: ([5.4], 40), 8: ([5.4], 40)}
BIN_SWEEPS = 5            # sweeps per production bin
BINS_FIRST = 10           # production bins (after Nskip) in round 1
BINS_MORE = 20            # production bins added per further round
MAX_ROUNDS = int(os.environ.get("RCHI_MAX_ROUNDS", "12"))
MIN_FREE_GB = 5.5
KEEP_RE = re.compile(r"^(SpinZ_eq|Ener_scal)$")

def ndelay(L):
    return 32 if L <= 28 else 64

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

# ------------------------------------------------------------------ parameters / runs
def params_text(L, U, nsweep, nbin):
    return f"""&VAR_ham_name
ham_name = "Hubbard"
/
&VAR_lattice
L1 = {L}
L2 = {L}
Lattice_type = "Square"
Model = "Hubbard"
/
&VAR_Model_Generic
Checkerboard = .T.
Symm = .T.
N_SUN = 2
N_FL = 1
Phi_X = 0.d0
Phi_Y = 0.d0
Bulk = .T.
N_Phi = {L*L//2}
Dtau = 0.1d0
Beta = {L}.d0
Projector = .F.
Theta = 10.d0
/
&VAR_QMC
Nwrap = 10
NSweep = {nsweep}
NBin = {nbin}
Ltau = 0
LOBS_ST = 0
LOBS_EN = 0
CPU_MAX = 0.0
Propose_S0 = .F.
Global_moves = .F.
N_Global = 1
Global_tau_moves = .F.
N_Global_tau = 1
Sequential = .T.
Nt_sequential_start = 0
Nt_sequential_end = -1
Langevin = .F.
HMC = .F.
Sequential_scheme = "submatrix"
N_delay = {ndelay(L)}
/
&VAR_errors
n_skip = 1
N_rebin = 1
N_Cov = 0
Extended_Zone = .F.
N_BZ_Zones = 1
/
&VAR_Hubbard
Mz = .F.
Continuous = .F.
ham_T = 1.d0
ham_chem = 0.d0
ham_U = {U:.6f}d0
ham_h0 = 0.d0
/
! slash terminates namelist statement - DO NOT REMOVE
"""

class Ctx:
    def __init__(self, a):
        self.work = os.path.abspath(os.path.expanduser(a.work))
        self.bin = os.path.abspath(os.path.expanduser(a.bin))
        self.quick = a.quick
        self.deadline = time.time() + a.deadline_hours * 3600
        self.Ls = LS_QUICK if a.quick else LS_FULL
        self.Us = US_QUICK if a.quick else US_FULL
        self.warm = WARM_QUICK if a.quick else WARM_FULL
        self.cpu = a.cpu
        os.makedirs(self.work + "/logs", exist_ok=True)
        self.logf = open(self.work + "/logs/pipeline.log", "a")
    def log(self, *m):
        s = f"[{now()}] " + " ".join(str(x) for x in m)
        print(s, flush=True); self.logf.write(s + "\n"); self.logf.flush()
    def free_gb(self):
        return shutil.disk_usage(self.work).free / 1e9
    def late(self):
        return time.time() > self.deadline

def seed_for(L, U, tag):
    return 100003 + 7919 * L + int(round(U * 1000)) * 31 + tag * 1000003

def nbins_done(rundir):
    p = rundir + "/Ener_scal"
    if not os.path.exists(p):
        return 0
    return sum(1 for l in open(p) if l.strip())

def run_segment(ctx, rundir, L, U, nsweep, nbin_new, tag):
    """Run ALF in rundir for nbin_new further bins (continues from confout_0 when present).  Returns True on success."""
    os.makedirs(rundir, exist_ok=True)
    first = not os.path.exists(rundir + "/confout_0")
    open(rundir + "/parameters", "w").write(params_text(L, U, nsweep, nbin_new))
    if first:
        lines = open(ctx.work + "/seeds_template").read().split("\n")
        lines[0] = f"{seed_for(L, U, tag):12d}"
        open(rundir + "/seeds", "w").write("\n".join(lines))
    else:
        if os.path.exists(rundir + "/confin_0"):
            os.remove(rundir + "/confin_0")
        os.rename(rundir + "/confout_0", rundir + "/confin_0")
    if os.path.exists(rundir + "/RUNNING"):
        os.remove(rundir + "/RUNNING")
    env = dict(os.environ, TMPDIR=ctx.work + "/tmp", OMP_NUM_THREADS="1")
    if not ctx.cpu:
        env["OMP_TARGET_OFFLOAD"] = "MANDATORY"
    n0 = nbins_done(rundir)
    t0 = time.time()
    with open(rundir + "/run.log", "a") as lg:
        r = subprocess.run([ctx.bin], cwd=rundir, stdout=lg, stderr=subprocess.STDOUT, env=env)
    dt = time.time() - t0
    n1 = nbins_done(rundir)
    ok = (r.returncode == 0 and n1 == n0 + nbin_new and os.path.exists(rundir + "/confout_0"))
    if ok and first and not ctx.cpu:
        info = open(rundir + "/info").read()
        if "GPU path (gpu_state): ACTIVE" not in info or "device equal-time measurement ACTIVE" not in info:
            ctx.log("ABORT: GPU path or device measurement not active in", rundir)
            sys.exit(3)
    if ok:
        # the unused equal-time files (Green_eq, Den_eq) are only archived: gzip the segment, ALF starts new files next time
        for f in ("Green_eq", "Den_eq"):
            if os.path.exists(f"{rundir}/{f}"):
                with open(f"{rundir}/{f}", "rb") as fi, gzip.open(f"{rundir}/{f}.seg{n0:05d}.gz", "wb", 6) as fo:
                    shutil.copyfileobj(fi, fo)
                os.remove(f"{rundir}/{f}")
    ctx.log(f"{'ok ' if ok else 'FAIL'} L={L} U={U} bins {n0}->{n1} ({nsweep} sweeps/bin) rc={r.returncode} {dt:.0f}s free={ctx.free_gb():.1f}GB")
    return ok

# ------------------------------------------------------------------ parsing
NUM = re.compile(r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?")
def numbers(path):
    t = open(path).read().replace("D", "E").replace("(", " ").replace(")", " ").replace(",", " ")
    return np.array(t.split(), dtype=float)

_kcache = {}
def read_sq(path, L):
    """SpinZ_eq (k space): returns kidx map and array S[nbin, L*L] (real part), index = (nx*L + ny), n = round(k L/2pi) mod L."""
    x = numbers(path)
    nk = L * L
    blk = 5 + 4 * nk
    if len(x) % blk:
        raise ValueError(f"{path}: {len(x)} numbers, not a multiple of {blk}")
    nb = len(x) // blk
    x = x.reshape(nb, blk)[:, 5:].reshape(nb, nk, 4)
    key = L
    if key not in _kcache:
        n = np.rint(x[0, :, :2] * L / (2 * np.pi)).astype(int) % L
        idx = np.full((L, L), -1, int)
        for i, (a, b) in enumerate(n):
            idx[a, b] = i
        if (idx < 0).any():
            raise ValueError("k grid is not a full L x L grid of multiples of 2pi/L")
        _kcache[key] = idx
    return _kcache[key], x[:, :, 2]

def read_scal(path):
    x = numbers(path)
    return x.reshape(-1, 4)[:, 1]          # (count, re, im, sign) per bin

def run_series(rundir, L):
    idx, S = read_sq(rundir + "/SpinZ_eq", L)
    h = L // 2
    SQ = S[:, idx[h, h]]
    SD = 0.25 * (S[:, idx[(h + 1) % L, h]] + S[:, idx[(h - 1) % L, h]] + S[:, idx[h, (h + 1) % L]] + S[:, idx[h, (h - 1) % L]])
    E = read_scal(rundir + "/Ener_scal")
    n = min(len(SQ), len(E))
    return dict(SQ=SQ[:n], SD=SD[:n], E=E[:n], Smean=S[:n].mean(axis=1))

# ------------------------------------------------------------------ stage 1: equilibration study
def moving_average(x, w):
    c = np.cumsum(np.insert(x, 0, 0.0))
    return (c[w:] - c[:-w]) / w

def estimate_nskip(x, w=10, nsig=3.0):
    """Sweeps to discard: first index after which the moving average (window w) of x stays within nsig standard deviations
    (of the moving average over the second half of the series) of the second-half mean."""
    n = len(x)
    if n < 4 * w:
        return n
    tail = x[n // 2:]
    mu = tail.mean()
    ma_tail = moving_average(tail, w)
    sig = max(ma_tail.std(), 1e-300)
    ma = moving_average(x, w)
    bad = np.where(np.abs(ma - mu) > nsig * sig)[0]
    return 0 if len(bad) == 0 else int(bad[-1] + w)

def stage_warm(ctx):
    res = ctx.work + "/nskip.json"
    if os.path.exists(res):
        return json.load(open(res))
    ctx.log("stage warm: equilibration study")
    data = {}
    for L, (Us, nb) in ctx.warm.items():
        for U in Us:
            d = f"{ctx.work}/warm/L{L}_U{U:.3f}"
            if not os.path.exists(d + "/confout_0"):
                if os.path.exists(d):
                    shutil.rmtree(d)
                if not run_segment(ctx, d, L, U, 1, nb, 9):
                    ctx.log("warm run failed", L, U); sys.exit(2)
            data[(L, U)] = run_series(d, L)
    perL = {}
    for (L, U), s in data.items():
        ns = max(estimate_nskip(s["SQ"]), estimate_nskip(s["E"]))
        perL[L] = max(perL.get(L, 0), ns)
        ctx.log(f"warm: L={L} U={U}: Nskip(S(Q))={estimate_nskip(s['SQ'])} Nskip(E)={estimate_nskip(s['E'])} of {len(s['SQ'])} sweeps")
    Lw = sorted(perL)
    nskip = {}
    for L in ctx.Ls:
        v = float(np.interp(L, Lw, [perL[l] for l in Lw]))
        nskip[str(L)] = int(max(20, math.ceil(1.5 * v + 5)))       # safety factor 1.5
    json.dump(dict(per_L_measured=perL, nskip_sweeps=nskip, rule="3-sigma moving average (w=10) of S(Q), E; x1.5+5, min 20, linear in L"),
              open(res, "w"), indent=1)
    plot_warm(ctx, data, perL, nskip)
    ctx.log("nskip (sweeps):", nskip)
    return json.load(open(res))

def plot_warm(ctx, data, perL, nskip):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Lw = sorted(perL)
    fig, ax = plt.subplots(2, len(Lw), figsize=(4.6 * len(Lw), 7), squeeze=False)
    cols = {}
    for (L, U), s in sorted(data.items()):
        j = Lw.index(L)
        for i, (key, lab) in enumerate((("SQ", r"$S(Q)/\langle S(Q)\rangle_{tail}$"), ("E", r"$E/\langle E\rangle_{tail}$"))):
            x = s[key]; n = len(x); mu = x[n // 2:].mean()
            c = cols.setdefault(U, f"C{len(cols)}")
            ax[i][j].plot(np.arange(1, n + 1), x / mu, ".", ms=2, color=c, alpha=0.35)
            w = min(10, n // 4)
            ax[i][j].plot(np.arange(w, n + 1), moving_average(x / mu, w), "-", color=c, lw=1.6, label=f"U={U}")
            ax[i][j].axhline(1.0, color="k", lw=0.5)
        ax[0][j].set_title(f"L={L}")
    for j, L in enumerate(Lw):
        for i in range(2):
            ax[i][j].axvline(nskip[str(L)] if str(L) in nskip else perL[L], color="r", ls="--")
            ax[i][j].set_xlabel("sweep"); ax[i][j].set_ylabel(["S(Q) / tail mean", "E / tail mean"][i])
            ax[i][j].legend(fontsize=8)
            lo, hi = ax[i][j].get_ylim(); ax[i][j].set_ylim(max(lo, 0.5), min(hi, 1.5))
    fig.suptitle(r"Equilibration from a random start (1 sweep per bin); dashed: Nskip used in production")
    fig.tight_layout(); fig.savefig(ctx.work + "/warmup.png", dpi=130); plt.close(fig)

# ------------------------------------------------------------------ stage 2: production
def rundir(ctx, L, U):
    return f"{ctx.work}/prod/L{L}_U{U:.3f}"

def stage_prod(ctx, nskip):
    state_f = ctx.work + "/state.json"
    st = json.load(open(state_f)) if os.path.exists(state_f) else {"round": 0}
    for rnd in range(max(1, st["round"] + 1) if False else 1, MAX_ROUNDS + 1):
        if st.get("done_round", 0) >= rnd:
            continue
        for U in ctx.Us:
            for L in ctx.Ls:
                if ctx.late():
                    ctx.log("deadline reached, stopping production"); return
                if ctx.free_gb() < MIN_FREE_GB:
                    ctx.log(f"STOP: free disk {ctx.free_gb():.1f} GB < {MIN_FREE_GB}"); return
                nsk = math.ceil(nskip[str(L)] / BIN_SWEEPS)
                target = nsk + BINS_FIRST + BINS_MORE * (rnd - 1)
                d = rundir(ctx, L, U)
                have = nbins_done(d)
                if have >= target:
                    continue
                for attempt in range(2):
                    if run_segment(ctx, d, L, U, BIN_SWEEPS, target - have, 1):
                        break
                    ctx.log("retry" if attempt == 0 else "giving up on this run for now")
                    if os.path.exists(d + "/confin_0") and not os.path.exists(d + "/confout_0"):
                        os.rename(d + "/confin_0", d + "/confout_0")
                    if attempt == 0:
                        have = nbins_done(d)
                        if have >= target: break
            try:
                final_analysis(ctx, nskip, tag=f"round{rnd}_after_U{U}")
            except Exception as e:
                ctx.log("analysis failed (continuing):", repr(e))
        st["done_round"] = rnd
        json.dump(st, open(state_f, "w"))
        ctx.log(f"round {rnd} complete")

# ------------------------------------------------------------------ stage 3: analysis
def collect(ctx, nskip, jmax=20):
    """data[(L,U)] = dict(n, SQ, SD, jk_SQ(J), jk_SD(J)) after removing the first Nskip sweeps (in bins)."""
    runs = {}
    for U in sorted(ctx.Us):
        for L in ctx.Ls:
            d = rundir(ctx, L, U)
            if not os.path.exists(d + "/SpinZ_eq"):
                continue
            s = run_series(d, L)
            k = math.ceil(nskip[str(L)] / BIN_SWEEPS)
            for key in s:
                s[key] = s[key][k:]
            if len(s["SQ"]) >= 4:
                runs[(L, U)] = s
    if not runs:
        raise RuntimeError("no production data yet")
    J = int(min(jmax, min(len(s["SQ"]) for s in runs.values())))
    out = {}
    for key, s in runs.items():
        n = len(s["SQ"])
        edges = np.linspace(0, n, J + 1).astype(int)
        jk = {}
        for q in ("SQ", "SD"):
            tot = s[q].sum()
            bs = np.array([s[q][edges[j]:edges[j + 1]].sum() for j in range(J)])
            nb = np.diff(edges)
            jk[q] = (tot - bs) / (n - nb)
        out[key] = dict(n=n, SQ=s["SQ"].mean(), SD=s["SD"].mean(), jkSQ=jk["SQ"], jkSD=jk["SD"], Smean=s["Smean"].mean())
    return out, J

def crossing_analysis(Ls, Us_by_L, tabs_R, tabs_S, u0=5.4):
    """tabs_R[L] = dict U->R, tabs_S[L] = dict U->S(Q) for ONE sample.  Returns dict pair-> (Uc, invnu, eta)."""
    from scipy.interpolate import CubicSpline
    from scipy.optimize import brentq
    spl = {}
    for L in Ls:
        us = sorted(tabs_R.get(L, {}))
        if len(us) >= 4:
            spl[L] = (CubicSpline(us, [tabs_R[L][u] for u in us]), CubicSpline(us, [np.log(tabs_S[L][u]) for u in us]), us[0], us[-1])
    res = {}
    for L1, L2 in zip(Ls[:-1], Ls[1:]):          # neighbouring sizes (L, L+4 for the production set)
        if True:
            if L1 not in spl or L2 not in spl:
                continue
            f1, g1, a1, b1 = spl[L1]; f2, g2, a2, b2 = spl[L2]
            lo, hi = max(a1, a2), min(b1, b2)
            if hi <= lo:
                continue
            grid = np.linspace(lo, hi, 2001)
            d = f1(grid) - f2(grid)
            sc = np.where(((d[:-1] <= 0) & (d[1:] > 0)) | ((d[:-1] >= 0) & (d[1:] < 0)))[0]
            if len(sc) == 0:
                continue
            roots = [brentq(lambda u: f1(u) - f2(u), grid[k], grid[k + 1]) for k in sc]
            uc = min(roots, key=lambda r: abs(r - u0))
            d1, d2 = f1.derivative()(uc), f2.derivative()(uc)
            invnu = math.log(d2 / d1) / math.log(L2 / L1) if d1 > 0 and d2 > 0 else float("nan")
            slope = (g2(uc) - g1(uc)) / math.log(L2 / L1)               # S(Q) ~ L^(1-eta)
            res[(L1, L2)] = (uc, invnu, 1.0 - slope)
    return res

def fit_ext(x_L, y):
    """y(L) = a + b / L^c; returns a (falls back to c=1 for < 4 points / failed fits)."""
    from scipy.optimize import curve_fit
    x_L = np.asarray(x_L, float); y = np.asarray(y, float)
    f = lambda L, a, b, c: a + b / L ** c
    if len(x_L) >= 4:
        try:
            p, _ = curve_fit(f, x_L, y, p0=[y[np.argmax(x_L)], (y[0] - y[-1]) * x_L[0], 1.0], bounds=([-10, -1e4, 0.2], [10, 1e4, 6.0]), maxfev=20000)
            return p
        except Exception:
            pass
    if len(x_L) >= 2:
        b, a = np.polyfit(1.0 / x_L, y, 1)
        return np.array([a, b, 1.0])
    return np.array([y[0], 0.0, 1.0])

def full_estimates(data, Ls, J, ctx):
    """Sample 0 = full data, samples 1..J = jackknife."""
    Us = sorted({u for (_, u) in data})
    out = []
    for s in range(J + 1):
        R, S = {}, {}
        for (L, U), d in data.items():
            sq = d["SQ"] if s == 0 else d["jkSQ"][s - 1]
            sd = d["SD"] if s == 0 else d["jkSD"][s - 1]
            R.setdefault(L, {})[U] = 1.0 - sd / sq
            S.setdefault(L, {})[U] = sq
        pairs = crossing_analysis(sorted(R), None, R, S)
        out.append((R, pairs))
    return out

def jk_err(vals):
    v = np.asarray(vals, float)
    J = len(v)
    return math.sqrt((J - 1) / J * np.sum((v - v.mean()) ** 2))

def final_analysis(ctx, nskip, tag="final", plot=True):
    data, J = collect(ctx, nskip)
    Ls = sorted({L for (L, _) in data}); Us = sorted({u for (_, u) in data})
    ests = full_estimates(data, Ls, J, ctx)
    R0, pairs0 = ests[0]
    # R(L,U) with jackknife errors
    Rtab = {}
    for (L, U) in data:
        jk = [ests[s][0][L][U] for s in range(1, J + 1)]
        Rtab[(L, U)] = (R0[L][U], jk_err(jk))
    summ = dict(tag=tag, J=J, Ls=Ls, Us=Us, R={f"{L}_{U}": v for (L, U), v in Rtab.items()},
                bins={f"{L}_{U}": d["n"] for (L, U), d in data.items()}, pairs={})
    series = {}
    for pr, (uc, inu, eta) in pairs0.items():
        jk = [ests[s][1].get(pr) for s in range(1, J + 1)]
        jk = [x for x in jk if x is not None]
        err = [jk_err([x[i] for x in jk]) if len(jk) >= 3 else float("nan") for i in range(3)]
        summ["pairs"][f"{pr[0]}_{pr[1]}"] = dict(Uc=uc, invnu=inu, eta=eta, err=err, n_jk=len(jk))
        series[pr] = (uc, inu, eta, err)
    # extrapolations (jackknife of the whole fit)
    def extrap(sel):
        prs = sorted(series)
        if len(prs) < 2:
            return None
        x = [p[0] for p in prs]
        f0 = fit_ext(x, [pairs0[p][sel] for p in prs])
        jkv = []
        for s in range(1, J + 1):
            pj = ests[s][1]
            pp = [p for p in prs if p in pj]
            if len(pp) < 2:
                continue
            jkv.append(fit_ext([p[0] for p in pp], [pj[p][sel] for p in pp])[0])
        return f0, (jk_err(jkv) if len(jkv) >= 3 else float("nan"))
    ex = {k: extrap(i) for i, k in enumerate(("Uc", "invnu", "eta"))}
    if ex["Uc"]:
        summ["Uc_inf"] = (float(ex["Uc"][0][0]), ex["Uc"][1])
    if ex["invnu"]:
        a = float(ex["invnu"][0][0]); e = ex["invnu"][1]
        summ["nu"] = (1.0 / a, e / a ** 2)
    if ex["eta"]:
        summ["eta_phi"] = (float(ex["eta"][0][0]), ex["eta"][1])
    json.dump(summ, open(f"{ctx.work}/summary.json", "w"), indent=1, default=float)
    write_md(ctx, summ, Rtab, series)
    if plot:
        make_plot(ctx, summ, Rtab, series, ex)
    ctx.log(f"analysis ({tag}): {len(data)} runs, J={J}, Uc={summ.get('Uc_inf')}, nu={summ.get('nu')}, eta={summ.get('eta_phi')}")
    return summ

def write_md(ctx, summ, Rtab, series):
    L = ["# R_chi^AFM crossing analysis (" + summ["tag"] + ")", "",
         f"jackknife blocks J={summ['J']}; Uc(L), 1/nu(L), eta(L) from crossings of neighbouring sizes (L, L+4); 1/L = 1/L_small", "",
         "| pair | Uc | 1/nu | eta_phi |", "|---|---|---|---|"]
    for k, v in summ["pairs"].items():
        L.append(f"| {k} | {v['Uc']:.4f}({v['err'][0]*1e4:.0f}) | {v['invnu']:.3f}({v['err'][1]*1e3:.0f}) | {v['eta']:.3f}({v['err'][2]*1e3:.0f}) |")
    for k in ("Uc_inf", "nu", "eta_phi"):
        if k in summ:
            L.append(f"\n{k} = {summ[k][0]:.4f} +- {summ[k][1]:.4f}")
    L += ["", "| L | U | bins used | R +- err |", "|---|---|---|---|"]
    for (l, u), (r, e) in sorted(Rtab.items()):
        L.append(f"| {l} | {u} | {summ['bins'][f'{l}_{u}']} | {r:.5f} +- {e:.5f} |")
    open(ctx.work + "/summary.md", "w").write("\n".join(L) + "\n")

def make_plot(ctx, summ, Rtab, series, ex):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mk = ["s", "o", "^", "v", "D", "p", "<", ">", "h"]
    cl = ["#2b5ca8", "#4daf4a", "#f5a030", "#e6452f", "#7b5aa6", "#444444", "#999999", "#17becf", "#8c564b"]
    fig = plt.figure(figsize=(8, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.28, wspace=0.28)
    a = fig.add_subplot(gs[0, :])
    Ls = summ["Ls"]
    for i, L in enumerate(Ls):
        pts = sorted((u, v[0], v[1]) for (l, u), v in Rtab.items() if l == L)
        if not pts:
            continue
        u, r, e = map(np.array, zip(*pts))
        a.errorbar(u, r, e, marker=mk[i % len(mk)], mfc="none", color=cl[i % len(cl)], ms=7, lw=1.4, capsize=2, label=f"L={L}")
    a.set_xlabel("U"); a.set_ylabel(r"$R_\chi^{\rm AFM}$"); a.legend(loc="lower right", fontsize=9, ncol=1)
    a.text(0.01, 0.97, "a", transform=a.transAxes, fontweight="bold", va="top")
    if series:
        ins = a.inset_axes([0.12, 0.62, 0.32, 0.33])
        prs = sorted(series)
        x = np.array([1.0 / p[0] for p in prs]); y = np.array([series[p][0] for p in prs]); e = np.array([series[p][3][0] for p in prs])
        ins.errorbar(x, y, e, fmt="s", mfc="none", color="k", capsize=2)
        if ex["Uc"]:
            p = ex["Uc"][0]; xs = np.linspace(1e-4, x.max() * 1.05, 100)
            ins.plot(xs, p[0] + p[1] / (1 / xs) ** p[2], "k-", lw=1)
            ins.errorbar([0], [p[0]], [ex["Uc"][1]], fmt="s", mfc="none", color="k", capsize=2, label=f"$U_c$={p[0]:.2f}({ex['Uc'][1]*100:.0f})")
            ins.legend(fontsize=7, loc="lower right")
        ins.set_xlabel("1/L", fontsize=8); ins.set_ylabel(r"$U_c(L)$", fontsize=8); ins.tick_params(labelsize=7)
    for pos, key, ylab, lab in ((gs[1, 0], "invnu", r"$1/\nu(L)$", "b"), (gs[1, 1], "eta", r"$\eta_\phi(L)$", "c")):
        b = fig.add_subplot(pos)
        sel = 1 if key == "invnu" else 2
        prs = sorted(series)
        if prs:
            x = np.array([1.0 / p[0] for p in prs]); y = np.array([series[p][sel] for p in prs]); e = np.array([series[p][3][sel] for p in prs])
            b.errorbar(x, y, e, fmt="s", mfc="none", color="k", capsize=2)
            f = ex["invnu" if key == "invnu" else "eta"]
            if f:
                p = f[0]; xs = np.linspace(1e-4, x.max() * 1.05, 100)
                b.plot(xs, p[0] + p[1] / (1 / xs) ** p[2], "k-", lw=1, label=r"$a+b/L^c$")
                val = (1 / p[0], f[1] / p[0] ** 2) if key == "invnu" else (p[0], f[1])
                b.errorbar([0], [p[0]], [f[1]], fmt="s", mfc="none", color="k", capsize=2,
                           label=(r"$\nu$" if key == "invnu" else r"$\eta_\phi$") + f" = {val[0]:.3f}({min(val[1]*1e3, 999):.0f})")
            b.legend(fontsize=9)
        b.set_xlabel("1/L"); b.set_ylabel(ylab); b.text(-0.12, 1.02, lab, transform=b.transAxes, fontweight="bold", va="bottom")
        b.set_xlim(-0.003, None)
    fig.suptitle(f"pi-flux Hubbard, beta = L, GPU submatrix DQMC  ({summ['tag']}, J={summ['J']} jackknife blocks)", fontsize=10)
    fig.savefig(ctx.work + "/Rchi_crossing.png", dpi=150, bbox_inches="tight")
    fig.savefig(ctx.work + "/Rchi_crossing.pdf", bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "analyze"])
    ap.add_argument("--work", default="~/alf-gpu/work/rchi")
    ap.add_argument("--bin", default="~/alf-gpu/work/fg1/bin/ALF-gpu-prod.out")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--cpu", action="store_true", help="CPU binary (no GPU checks); test only")
    ap.add_argument("--deadline-hours", type=float, default=72.0)
    ap.add_argument("--seeds", default="~/alf-gpu/work/fg1/src_gpu/Scripts_and_Parameters_files/Start/seeds")
    a = ap.parse_args()
    ctx = Ctx(a)
    os.makedirs(ctx.work + "/tmp", exist_ok=True)
    if not os.path.exists(ctx.work + "/seeds_template"):
        shutil.copy(os.path.expanduser(a.seeds), ctx.work + "/seeds_template")
    if a.cmd == "analyze":
        nskip = json.load(open(ctx.work + "/nskip.json"))["nskip_sweeps"]
        final_analysis(ctx, nskip); return
    ctx.log("=== pipeline start", "quick" if a.quick else "full", "work", ctx.work, "bin", ctx.bin)
    sh = ctx.work + "/bin"; os.makedirs(sh, exist_ok=True)
    frozen = sh + "/ALF.out"
    if not os.path.exists(frozen):
        shutil.copy(ctx.bin, frozen)
    ctx.bin = frozen
    nskip = stage_warm(ctx)["nskip_sweeps"]
    stage_prod(ctx, nskip)
    final_analysis(ctx, nskip, tag="final")
    ctx.log("=== pipeline finished; results:", ctx.work + "/Rchi_crossing.png", ctx.work + "/summary.md", ctx.work + "/warmup.png")
    open(ctx.work + "/DONE", "w").write(now() + "\n")

if __name__ == "__main__":
    main()
