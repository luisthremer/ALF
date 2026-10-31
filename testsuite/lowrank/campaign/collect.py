#!/usr/bin/env python3
"""Collect a campaign's runs into one table with full provenance (runs locally, stdlib only).

Usage: python3 -I collect.py CAMPAIGN_DIR [--out FILE.json] [--strict [--allow KEY ...]]

One row per run of campaign.json, from the latest attempt runs/<id>/<latest>/: run.prov, the job's job.prov, info
(CPU Time, Acceptance, Precision, commit, START/FIN TIME, Sequential scheme / N_delay in every format) and timers.txt
(T1 format, see parse_timers; "absent" if the binary was built without -DALF_TIMERS; warnings for a restart block,
sweep != sum of the nine categories, flush_unattributed != 0, or sweep + bin_overhead off info's CPU Time by > 1 %). Exit 1 if an OK run lacks a required key, a .prov
file repeats a key, a run appears twice, or the job script's sha256 differs from campaign.json. Runs whose status is
not OK (FAILED, INVALID pinning, or info contradicting the request) are excluded from the statistics and listed with the
reason; runs never started are listed as missing. --strict also fails on any required field reading "unknown" or
"unavailable(...)", and on cluster profiles (not 'local') on i_mpi_root=unset and on omp_affinity "none printed" with
THREADS > 1; --allow KEY (repeatable; job./run./info. prefix optional) accepts a field a node legitimately lacks.
Output: CAMPAIGN_DIR/collected.json (default) and a .tsv next to it.
"""
import datetime
import hashlib
import json
import os
import re
import sys

REQUIRED_JOB = [
    "campaign_id", "job_name", "slurm_job_id", "cluster", "node", "hostname", "profile", "node_list", "partition",
    "exclusive", "ear", "cpu_model", "cpus", "sockets", "cores_per_socket", "threads_per_core", "numa_nodes",
    "mem_total_mib", "governor_job_start", "freq_driver_job_start", "cur_freq_job_start", "max_freq_job_start",
    "boost_job_start", "mpi_version", "i_mpi_root", "modules_job", "pin_env", "omp_num_threads", "mkl_num_threads",
    "mkl_dynamic", "omp_display_affinity", "mkl_verbose_in_runs", "ranks", "threads", "job_script_sha256",
    "profile_sha256", "provenance_sh_sha256", "rank_wrap_sha256", "binary", "binary_sha256", "build_info",
    "build_info_sha256", "build_compiler", "build_compiler_version", "build_flags_prog", "build_alf_lib",
    "build_modules", "build_git_commit", "fakeintel", "fakeintel_lib", "fakeintel_lib_sha256", "mkl_probe_rc",
    "mkl_probe_line", "mkl_version", "mkl_isa", "mkl_threading", "fakeintel_maps_probe"]
REQUIRED_RUN = [
    "run_id", "attempt", "scheme_requested", "n_delay_requested", "job_prov", "mkl_verbose_env", "governor_start",
    "cur_freq_start", "max_freq_start", "boost_start", "run_start", "governor_end", "cur_freq_end", "wall_s",
    "wall_resolution",
    "launch_rc", "max_rss_kb", "max_rss_source", "rank_hosts", "rank_affinity_preinit", "fakeintel_maps",
    "ld_preload", "pin_map", "pin_check", "omp_affinity", "status"]
REQUIRED_INFO = ["commit", "start_time", "fin_time", "cpu_time_s", "acceptance", "pgreen_mean", "pgreen_max",
                 "pphase_max", "beta_info", "dtau_info", "sweeps_info", "bins_info", "scheme_used", "n_delay_mode",
                 "n_delay_used"]
# fields that may legitimately read unavailable(...) under --strict: the build tree's git commit (a binary built in a
# git-archive copy has none; the authoritative commit is info's "This executable represents commit" line, required)
ALLOWED_UNAVAILABLE = {"build_git_commit"}


class CollectError(Exception):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_prov(path):
    d = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            if "=" not in line:
                raise CollectError("%s:%d: not KEY=VALUE: %r" % (path, n, line))
            k, v = line.split("=", 1)
            if k in d:
                raise CollectError("%s:%d: duplicate key %r" % (path, n, k))
            d[k] = v
    return d


def _num(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse_ndelay(text):
    """Sequential scheme / N_delay of an info file -> (scheme_used, mode, n_used, n0, ndim). Tolerant of all formats:
    'Sequential scheme: X, N_delay = n' | '..., N_delay = auto (static n0 = a, Ndim = b)' | '..., N_delay = auto+tune
    (static n0 = a, Ndim = b)' followed by 'Auto N_delay: candidates ...: ... s/proposal, chosen n' or
    'Auto N_delay: tuning skipped (<reason>), N_delay = n0'. No scheme line: the standard scheme."""
    lines = re.findall(r"Sequential scheme:\s*([A-Za-z_]+)\s*,\s*N_delay\s*=\s*([^\n]*?)\s*$", text, re.M)
    if not lines:
        if re.search(r"Auto N_delay:", text):
            raise CollectError("'Auto N_delay:' line without a 'Sequential scheme:' line")
        return "STANDARD", "n/a", "n/a", "n/a", "n/a"
    if len(set(lines)) > 1:
        raise CollectError("conflicting 'Sequential scheme:' lines: %s" % sorted(set(lines)))
    scheme, rest = lines[0][0].upper(), lines[0][1]
    m = re.match(r"^(-?\d+)\b", rest)
    if m:
        if re.search(r"Auto N_delay:", text):
            raise CollectError("'Auto N_delay:' line with a fixed N_delay")
        return scheme, "fixed", int(m.group(1)), "n/a", "n/a"
    m = re.match(r"^auto(\s*\+\s*tune)?\s*\(\s*static\s+n0\s*=\s*(\d+)\s*,\s*Ndim\s*=\s*(\d+)\s*\)", rest, re.I)
    if not m:
        raise CollectError("unparsed N_delay value %r" % rest)
    tune, n0, ndim = bool(m.group(1)), int(m.group(2)), int(m.group(3))
    chosen = re.findall(r"Auto N_delay:[^\n]*?chosen\s+(\d+)", text)
    skipped = re.findall(r"Auto N_delay:\s*tuning skipped\s*\(([^)\n]*)\)\s*,\s*N_delay\s*=\s*(\d+)", text)
    if tune:
        if len(chosen) + len(skipped) != 1:
            raise CollectError("auto+tune: expected exactly one 'Auto N_delay:' result line, found %d"
                               % (len(chosen) + len(skipped)))
        if chosen:
            return scheme, "auto+tune", int(chosen[0]), n0, ndim
        return scheme, "auto+tune (skipped: %s)" % skipped[0][0].strip(), int(skipped[0][1]), n0, ndim
    if chosen or skipped:
        raise CollectError("'Auto N_delay:' result line without auto+tune")
    return scheme, "auto", n0, n0, ndim


def parse_info(path):
    text = open(path, errors="replace").read()

    def one(rx, conv=str, group=1):
        m = re.search(rx, text, re.M)
        return conv(m.group(group)) if m else None

    d = {
        "commit": one(r"This executable represents commit\s+(\S+)"),
        "start_time": one(r"START TIME:\s*(\S+)"),
        "fin_time": one(r"FIN TIME:\s*(\S+)"),
        "cpu_time_s": one(r"^\s*CPU Time\s*:\s*(\S+)", _num),
        "acceptance": one(r"^\s*Acceptance\s*:\s*(\S+)", _num),
        "pgreen_mean": one(r"^\s*Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)", _num),
        "pgreen_max": one(r"^\s*Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)", _num, 2),
        "pphase_max": one(r"^\s*Precision Phase, Max\s*:\s*(\S+)", _num),
        "beta_info": one(r"^\s*Beta\s*:\s*(\S+)", _num),
        "dtau_info": one(r"^\s*dtau,Ltrot_eff\s*:\s*(\S+)", _num),
        "sweeps_info": one(r"^\s*Sweeps\s*:\s*(\d+)", int),
        "bins_info": one(r"^\s*Bins\s*:\s*(\d+)", int),
        "model_info": one(r"^\s*Model is\s*:\s*(\S+)") or "n/a (not printed by this model)",
    }
    s, mode, used, n0, ndim = parse_ndelay(text)
    d.update(scheme_used=s, n_delay_mode=mode, n_delay_used=used, n_delay_n0=n0, n_delay_ndim=ndim)
    if d["start_time"] and d["fin_time"] and d["cpu_time_s"] is not None:
        t0 = datetime.datetime.fromisoformat(d["start_time"])
        t1 = datetime.datetime.fromisoformat(d["fin_time"])
        d["info_wall_minus_cpu_s"] = round((t1 - t0).total_seconds() - d["cpu_time_s"], 3)
    return d


TIMER_CATEGORIES = ("update", "flush", "vwrap", "hop", "stab", "meas", "global", "langevin_hmc", "other")
TIMER_TOTAL = "sweep"                 # inclusive body of the sweep loop = sum of the nine exclusive categories
TIMER_EXTRA = ("bin_overhead",)       # per-bin wall time outside the sweeps (not part of sweep)


def parse_timers(path):
    """timers.txt of a -DALF_TIMERS binary (T1, Prog/Timers_mod.F90; testsuite/lowrank/README.md):
    '#' comment lines; header TIMER_RATE <count_rate>, TIMER_CLOCK <name>, TIMER_RANKS <n>, TIMER_TICK <s>; per bin
    TIMER_SWEEPS <bin> <n>, TIMER <bin> <id> <rank0> <max> <mean> (seconds, ES format) for the nine exclusive categories,
    sweep and bin_overhead, COUNTER <bin> <name> <rank0> <max> <mean> (rank0/max integers, mean real; flush_rank_mean all
    real). A second header (TIMER_RATE after bin data: a restart in the same directory, the file is appended) starts a
    new block; the last block is kept and a warning recorded. Checks (warnings, not errors): sweep = sum of the nine per
    bin (rank0 and mean, within 1e-6 relative + 2 ticks per category), flush_unattributed = 0."""
    blocks, cur = [], None
    for n, line in enumerate(open(path), 1):
        w = line.split()
        if not w or w[0].startswith("#"):
            continue
        try:
            if w[0] == "TIMER_RATE" and len(w) == 2:
                cur = {"timers": {}, "counters": {}, "sweeps": {}, "rate": _num(w[1]), "clock": None, "ranks": None,
                       "tick": None, "first_line": n}
                blocks.append(cur)
                continue
            if cur is None:
                raise CollectError("%s:%d: data before the TIMER_RATE header: %r" % (path, n, line.rstrip()))
            if w[0] == "TIMER_CLOCK" and len(w) == 2:
                cur["clock"] = w[1]
            elif w[0] == "TIMER_RANKS" and len(w) == 2:
                cur["ranks"] = int(w[1])
            elif w[0] == "TIMER_TICK" and len(w) == 2:
                cur["tick"] = _num(w[1])
            elif w[0] in ("TIMER", "COUNTER") and len(w) == 6:
                key = "timers" if w[0] == "TIMER" else "counters"
                b = cur[key].setdefault(w[1], {})
                if w[2] in b:
                    raise CollectError("%s:%d: duplicate %s %s in bin %s" % (path, n, w[0], w[2], w[1]))
                b[w[2]] = {"rank0": _num(w[3]), "max": _num(w[4]), "mean": _num(w[5])}
            elif w[0] == "TIMER_SWEEPS" and len(w) == 3:
                if w[1] in cur["sweeps"]:
                    raise CollectError("%s:%d: duplicate TIMER_SWEEPS for bin %s" % (path, n, w[1]))
                cur["sweeps"][w[1]] = int(w[2])
            else:
                raise CollectError("%s:%d: unknown line %r" % (path, n, line.rstrip()))
        except ValueError:
            raise CollectError("%s:%d: bad number in %r" % (path, n, line.rstrip()))
    if not blocks:
        raise CollectError("%s: no TIMER_RATE header" % path)
    t = blocks[-1]
    t["warnings"] = []
    if len(blocks) > 1:
        t["warnings"].append("%s has %d header blocks (restart in the same directory?): kept the last one (from line %d)"
                             % (os.path.basename(path), len(blocks), t["first_line"]))
    tick = t["tick"] or 0.0
    for b in sorted(t["timers"], key=int):
        e = t["timers"][b]
        missing = [k for k in TIMER_CATEGORIES + (TIMER_TOTAL,) + TIMER_EXTRA if k not in e]
        if missing:
            raise CollectError("%s: bin %s lacks TIMER %s" % (path, b, ",".join(missing)))
        if b not in t["sweeps"]:
            raise CollectError("%s: bin %s has no TIMER_SWEEPS line" % (path, b))
        for col in ("rank0", "mean"):
            tot = sum(e[k][col] for k in TIMER_CATEGORIES)
            if abs(tot - e[TIMER_TOTAL][col]) > 1e-6 * abs(e[TIMER_TOTAL][col]) + 2 * len(TIMER_CATEGORIES) * tick + 1e-12:
                t["warnings"].append("bin %s: sum of the nine categories %.9g != sweep %.9g (%s)"
                                     % (b, tot, e[TIMER_TOTAL][col], col))
    for b, c in sorted(t["counters"].items(), key=lambda x: int(x[0])):
        u = c.get("flush_unattributed")
        if u is not None and (u["rank0"] != 0 or u["max"] != 0):
            t["warnings"].append("bin %s: flush_unattributed = %g (max %g), must be 0" % (b, u["rank0"], u["max"]))
    t["flush_unattributed_max"] = max([c["flush_unattributed"]["max"] for c in t["counters"].values()
                                       if "flush_unattributed" in c] or [0])
    return t


def strict_problems(rid, job, run, info, threads, allow=()):
    """--strict findings for one OK run (list of strings). allow: field names (with or without job./run./info.)."""
    def allowed(scope, k):
        return k in ALLOWED_UNAVAILABLE or k in allow or ("%s.%s" % (scope, k)) in allow
    out = []
    for scope, keys, d in (("run", REQUIRED_RUN, run), ("job", REQUIRED_JOB, job), ("info", REQUIRED_INFO, info)):
        for k in keys:
            if bad_value(d[k]) and not allowed(scope, k):
                out.append("%s %s.%s=%s" % (rid, scope, k, d[k]))
    if job.get("profile") != "local":
        if job.get("i_mpi_root") == "unset" and not allowed("job", "i_mpi_root"):
            out.append("%s job.i_mpi_root=unset (Intel MPI expected on the cluster)" % rid)
        if threads > 1 and run.get("omp_affinity", "").startswith("none printed") and not allowed("run", "omp_affinity"):
            out.append("%s run.omp_affinity=%s with THREADS=%d (no OpenMP affinity record)" % (rid, run["omp_affinity"], threads))
    return out


def bad_value(v):
    s = str(v)
    return s == "" or "unknown" in s.lower() or s.startswith("unavailable")


def collect(cdir, strict=False, allow=()):
    cj_path = os.path.join(cdir, "campaign.json")
    camp = json.load(open(cj_path))
    cid = camp["campaign_id"]
    job_sha = {j["job"]: j["script_sha256"] for j in camp["jobs"]}
    rows, excluded, missing, strict_bad, warnings = [], [], [], [], []
    seen = set()
    for req in camp["runs"]:
        rid = req["run_id"]
        if rid in seen:
            raise CollectError("run %s listed twice in campaign.json" % rid)
        seen.add(rid)
        rdir = os.path.join(cdir, "runs", rid)
        status_line = open(os.path.join(rdir, "status")).read().strip() if os.path.exists(os.path.join(rdir, "status")) else ""
        attempts = sorted(d for d in os.listdir(rdir) if re.match(r"^a\d+$", d)) if os.path.isdir(rdir) else []
        base = {"campaign_id": cid, "run_id": rid, "status": status_line or "never started", "attempts": len(attempts)}
        word = status_line.split()[0] if status_line else ""
        if word not in ("OK", "FAILED", "INVALID"):          # never started, NOT_RUN, RUNNING (node lost)
            missing.append(base)
            continue
        latest = open(os.path.join(rdir, "latest")).read().strip()
        adir = os.path.join(rdir, latest)
        if word in ("FAILED", "INVALID"):
            base["log"] = os.path.relpath(os.path.join(adir, "stdout.log"), cdir)
            excluded.append(base)
            continue
        run = read_prov(os.path.join(adir, "run.prov"))
        job = read_prov(os.path.join(cdir, run["job_prov"]))
        for k in REQUIRED_RUN:
            if k not in run:
                raise CollectError("%s: run.prov lacks required key %s" % (rid, k))
        for k in REQUIRED_JOB:
            if k not in job:
                raise CollectError("%s: %s lacks required key %s" % (rid, run["job_prov"], k))
        if job["job_script_sha256"] != job_sha[req["job"]]:
            raise CollectError("%s: job script sha256 %s differs from campaign.json %s"
                               % (rid, job["job_script_sha256"], job_sha[req["job"]]))
        if job["campaign_id"] != cid or run["run_id"] != rid:
            raise CollectError("%s: provenance belongs to %s/%s" % (rid, job["campaign_id"], run["run_id"]))
        info_path = os.path.join(adir, "info")
        info = parse_info(info_path)
        for k in REQUIRED_INFO:
            if info.get(k) is None:
                raise CollectError("%s: info lacks %s" % (rid, k))
        row = dict(base)
        row.update({k: v for k, v in req.items() if k not in ("keys",)})
        row["keys"] = ",".join("%s=%s" % kv for kv in sorted(req["keys"].items()))
        row.update({"job." + k: v for k, v in job.items()})
        row.update({"run." + k: v for k, v in run.items()})
        row.update(info)
        row["sec_per_sweep"] = info["cpu_time_s"] / (req["nsweep"] * req["nbin"])
        # setup = launch wall time (job script, around the launcher) - CPU Time (info: sweep phase only, the timer is
        # restarted after the initial Green function): MPI start-up, Hamiltonian setup, initial stabilisation, output
        row["setup_s"] = round(float(run["wall_s"]) - info["cpu_time_s"], 3)
        # sacct MaxRSS is job level (all steps so far): only an upper bound for this run
        row["rss_upper_bound"] = run["max_rss_source"].startswith("sacct")
        row["info_path"] = os.path.relpath(info_path, cdir)
        row["info_sha256"] = sha256(info_path)
        row["stdout_sha256"] = sha256(os.path.join(adir, "stdout.log"))
        row["campaign_json_sha256"] = sha256(cj_path)
        tpath = os.path.join(adir, "timers.txt")
        row["timers_present"] = "present" if os.path.exists(tpath) else "absent"
        row["timers"] = parse_timers(tpath) if os.path.exists(tpath) else None
        if row["timers"]:
            for w in row["timers"]["warnings"]:
                warnings.append("%s: timers: %s" % (rid, w))
            tsum = sum(e[TIMER_TOTAL]["rank0"] + e["bin_overhead"]["rank0"] for e in row["timers"]["timers"].values())
            row["timers_vs_cpu_time"] = tsum / info["cpu_time_s"] if info["cpu_time_s"] else None
            if row["timers_vs_cpu_time"] is not None and abs(row["timers_vs_cpu_time"] - 1.0) > 0.01:
                warnings.append("%s: timers: sum of sweep + bin_overhead over bins = %.4g x info's CPU Time (expected "
                                "within 1 %%)" % (rid, row["timers_vs_cpu_time"]))
        # consistency of the run with its request: excluded (listed), not fatal
        problems = []
        if info["scheme_used"] != req["scheme"].upper():
            problems.append("scheme used %s != requested %s" % (info["scheme_used"], req["scheme"]))
        nd = req["n_delay_requested"]
        if req["scheme"] != "standard":
            want = {0: "auto", -1: "auto+tune"}.get(nd, "fixed")
            if not info["n_delay_mode"].startswith(want) or (want == "fixed" and info["n_delay_used"] != nd):
                problems.append("N_delay %s/%s in info, requested %s" % (info["n_delay_mode"], info["n_delay_used"], nd))
        if info["sweeps_info"] != req["nsweep"] or info["bins_info"] != req["nbin"]:
            problems.append("sweeps/bins %s/%s in info, requested %s/%s"
                            % (info["sweeps_info"], info["bins_info"], req["nsweep"], req["nbin"]))
        if abs(info["beta_info"] - req["beta"]) > 1e-9 or abs(info["dtau_info"] - req["dtau"]) > 1e-9:
            problems.append("beta/dtau %s/%s in info, requested %s/%s"
                            % (info["beta_info"], info["dtau_info"], req["beta"], req["dtau"]))
        if int(job["ranks"]) != req["ranks"] or job["omp_num_threads"] != str(req["threads"]) \
                or job["mkl_num_threads"] != str(req["threads"]):
            problems.append("layout in job.prov differs from the request")
        bc = job.get("build_git_commit", "")
        if bad_value(bc):
            row["build_commit_match"] = "n/a (no git commit in build_info)"
        else:
            row["build_commit_match"] = "yes" if bc.startswith(str(info["commit"])) else "NO"
            if row["build_commit_match"] == "NO":
                warnings.append("%s: build_info git commit %s != info commit %s (build_info.sh not run right after "
                                "the build?)" % (rid, bc[:12], info["commit"]))
        if not job["exclusive"].startswith("yes"):
            warnings.append("%s: node not exclusive according to scontrol (%s)" % (rid, job["exclusive"]))
        if problems:
            base.update(status="INVALID: " + "; ".join(problems), log=row["info_path"])
            excluded.append(base)
            continue
        if strict:
            strict_bad += strict_problems(rid, job, run, info, req["threads"], allow)
        rows.append(row)
    if strict_bad:
        raise CollectError("--strict: fields without a real value:\n  " + "\n  ".join(strict_bad))
    return {"campaign_id": cid, "campaign_dir": os.path.abspath(cdir), "campaign_json_sha256": sha256(cj_path),
            "collector_sha256": sha256(os.path.abspath(__file__)),
            "collected_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "rows": rows, "excluded": excluded, "missing": missing, "warnings": warnings}


TSV_COLS = ["run_id", "point", "L", "layout", "scheme", "n_delay_requested", "n_delay_mode", "n_delay_used",
            "repetition", "sec_per_sweep", "acceptance", "setup_s", "run.max_rss_kb", "commit", "job.node",
            "job.cpu_model", "job.mkl_isa", "job.fakeintel", "timers_present", "info_sha256"]


def main(argv):
    strict = "--strict" in argv
    allow = []
    while "--allow" in argv:
        i = argv.index("--allow")
        if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
            sys.stderr.write("collect.py: --allow needs a KEY\n" + __doc__)
            return 2
        allow.append(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    args = [a for a in argv[1:] if not a.startswith("--")]
    out = None
    if "--out" in argv:
        i = argv.index("--out")
        if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
            sys.stderr.write("collect.py: --out needs a FILE\n" + __doc__)
            return 2
        out = argv[i + 1]
        args.remove(out)
    if len(args) != 1:
        sys.exit(__doc__)
    cdir = args[0]
    try:
        res = collect(cdir, strict, set(allow))
        res["strict"] = strict
        res["strict_allow"] = sorted(allow)
    except (CollectError, OSError, KeyError) as e:
        sys.stderr.write("collect.py: ERROR: %s\n" % e)
        return 1
    out = out or os.path.join(cdir, "collected.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    with open(re.sub(r"\.json$", "", out) + ".tsv", "w") as f:
        f.write("\t".join(TSV_COLS) + "\n")
        for r in res["rows"]:
            f.write("\t".join(str(r.get(c, "")) for c in TSV_COLS) + "\n")
    print("campaign %s: %d OK rows, %d excluded, %d missing -> %s"
          % (res["campaign_id"], len(res["rows"]), len(res["excluded"]), len(res["missing"]), out))
    for w in res["warnings"]:
        print("  WARNING  " + w)
    for e in res["excluded"]:
        print("  EXCLUDED %s: %s (attempts %d, log %s)" % (e["run_id"], e["status"], e["attempts"], e.get("log", "-")))
    for m in res["missing"]:
        print("  MISSING  %s: %s (attempts %d)" % (m["run_id"], m["status"], m["attempts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
