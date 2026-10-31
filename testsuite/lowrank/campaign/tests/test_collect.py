#!/usr/bin/env python3
"""Unit checks of collect.py's parsers on fixtures: every N_delay info format, timers.txt (real T1 sample), duplicate keys.
Run: python3 -I testsuite/lowrank/campaign/tests/test_collect.py   (prints one line per check, exit 0 iff all pass)"""
import os
import sys
import tempfile

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import collect  # noqa: E402
import report  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
fails = 0


def check(name, cond, detail=""):
    global fails
    print("%s %s %s" % ("PASS" if cond else "FAIL", name, detail))
    fails += 0 if cond else 1


def info(name):
    return collect.parse_info(os.path.join(FIX, name))


def raises(fn, *a):
    try:
        fn(*a)
    except collect.CollectError as e:
        return str(e)
    return None


d = info("info_fixed.txt")
check("fixed", (d["scheme_used"], d["n_delay_mode"], d["n_delay_used"]) == ("DELAYED", "fixed", 32), str(d["n_delay_used"]))
check("fixed numbers", d["cpu_time_s"] == 1.5 and d["commit"] == "26b354c7" and abs(d["info_wall_minus_cpu_s"] - 1.0) < 1e-9
      and d["sweeps_info"] == 20 and d["bins_info"] == 3 and d["beta_info"] == 2.0, str(d))
d = info("info_auto.txt")
check("auto", (d["scheme_used"], d["n_delay_mode"], d["n_delay_used"], d["n_delay_n0"], d["n_delay_ndim"])
      == ("SUBMATRIX", "auto", 24, 24, 144))
d = info("info_autotune.txt")
check("auto+tune chosen", (d["n_delay_mode"], d["n_delay_used"], d["n_delay_n0"]) == ("auto+tune", 24, 24))
d = info("info_autotune_skipped.txt")
check("auto+tune skipped", (d["n_delay_mode"], d["n_delay_used"]) == ("auto+tune (skipped: cross-slice)", 24))
check("auto+tune without result line is an error",
      raises(collect.parse_info, os.path.join(FIX, "info_autotune_noresult.txt")) is not None)
d = info("info_standard.txt")
check("standard (no scheme line)", (d["scheme_used"], d["n_delay_mode"], d["n_delay_used"]) == ("STANDARD", "n/a", "n/a"))
txt = open(os.path.join(FIX, "info_fixed.txt")).read()
check("chosen line without auto+tune is an error",
      raises(collect.parse_ndelay, txt + "Auto N_delay: candidates 8 16: 1E-6 2E-6 s/proposal, chosen 8\n") is not None)
check("conflicting scheme lines are an error",
      raises(collect.parse_ndelay, txt + "Sequential scheme: SUBMATRIX, N_delay = 32\n") is not None)
check("repeated identical scheme lines (one per Temp_ dir / rank) are fine",
      collect.parse_ndelay(txt + "Sequential scheme: DELAYED, N_delay = 32\n")[2] == 32)
check("tolerant spacing", collect.parse_ndelay("Sequential scheme:delayed ,N_delay=auto+tune(static n0=16,Ndim=64)\n"
                                               "Auto N_delay: candidates 8 16 32: 1 2 3 s/proposal, chosen 32\n")
      == ("DELAYED", "auto+tune", 32, 16, 64))

t = collect.parse_timers(os.path.join(FIX, "timers.txt"))
check("timers header (T1 format, # comments skipped)", t["rate"] == 1e9 and t["clock"] == "system_clock"
      and t["ranks"] == 1 and t["tick"] == 1e-6 and t["sweeps"] == {"1": 140, "2": 140, "3": 140}, str(t["sweeps"]))
check("timers values in seconds", t["timers"]["2"]["stab"]["max"] == 10.664966 and t["timers"]["1"]["sweep"]["rank0"] == 17.570315)
check("timers counters (integer rank0/max, real mean)", t["counters"]["2"]["accepts"]["max"] == 349337
      and abs(t["counters"]["1"]["flush_rank_mean"]["mean"] - 26.51377337) < 1e-9)
check("real sample: sweep = sum of the nine, flush_unattributed = 0: no warnings", t["warnings"] == [], str(t["warnings"]))
val, cnt = report.split_per_sweep(t)
e2, e3 = t["timers"]["2"], t["timers"]["3"]
check("timing split drops bin 1 (bins 2, 3 / 280 sweeps)",
      abs(val["sweep"] - (e2["sweep"]["mean"] + e3["sweep"]["mean"]) / 280) < 1e-12
      and abs(val["stab"] - (e2["stab"]["mean"] + e3["stab"]["mean"]) / 280) < 1e-12, str(val))
check("nine categories add up to sweep", abs(sum(val[c] for c in report.CATEGORIES) - val["sweep"]) < 1e-9)
check("bin_overhead separate, not in sweep", "bin_overhead" in val and "bin_overhead" not in report.CATEGORIES)
check("counters per sweep; flush_rank_mean averaged over bins", abs(cnt["proposals"] - 5760.0) < 1e-9
      and abs(cnt["flush_rank_mean"] - (26.44889461 + 26.49889461) / 2) < 1e-6, str(cnt))
t = collect.parse_timers(os.path.join(FIX, "timers_restart.txt"))
check("restart: two header blocks, last kept, warning", t["sweeps"]["1"] == 140 and len(t["warnings"]) == 1
      and "2 header blocks" in t["warnings"][0], str(t["warnings"]))
t = collect.parse_timers(os.path.join(FIX, "timers_unattributed.txt"))
check("flush_unattributed != 0 is flagged", t["flush_unattributed_max"] == 3
      and any("flush_unattributed = 3" in w for w in t["warnings"]), str(t["warnings"]))

with tempfile.TemporaryDirectory() as tmp:
    p = os.path.join(tmp, "x.prov")
    open(p, "w").write("a=1\nb=2\na=3\n")
    check("duplicate key in .prov is an error", raises(collect.read_prov, p) is not None)
    p2 = os.path.join(tmp, "t.txt")
    open(p2, "w").write("TIMER_RATE 1000000000\nTIMER 2 sweep 1 2 3\nTIMER 2 sweep 1 2 3\n")
    check("duplicate timer line is an error", raises(collect.parse_timers, p2) is not None)

# --strict rules (strict_problems): generic unavailable values, cluster-only rules, --allow
job = {k: "x" for k in collect.REQUIRED_JOB}
run = {k: "x" for k in collect.REQUIRED_RUN}
inf = {k: 1 for k in collect.REQUIRED_INFO}
job.update(profile="supermuc", i_mpi_root="unset")
run.update(omp_affinity="none printed (no OpenMP parallel region was entered)")
p = collect.strict_problems("r", job, run, inf, 4)
check("cluster: i_mpi_root=unset and omp_affinity none with THREADS>1 are bad", len(p) == 2, str(p))
check("cluster: omp_affinity none is fine with THREADS=1", len(collect.strict_problems("r", job, run, inf, 1)) == 1)
check("--allow job.i_mpi_root / omp_affinity", collect.strict_problems("r", job, run, inf, 4,
                                                                     {"job.i_mpi_root", "omp_affinity"}) == [])
job["profile"] = "local"
check("local profile: cluster rules off", collect.strict_problems("r", job, run, inf, 4) == [])
job["mem_total_mib"] = "unavailable(no free)"
check("unavailable value is bad", len(collect.strict_problems("r", job, run, inf, 4)) == 1)
check("--allow mem_total_mib", collect.strict_problems("r", job, run, inf, 4, {"mem_total_mib"}) == [])
job["build_git_commit"] = "unavailable(archive)"
check("build_git_commit may be unavailable", collect.strict_problems("r", job, run, inf, 4, {"mem_total_mib"}) == [])

print("test_collect: %s (%d failures)" % ("ALL PASS" if not fails else "FAILED", fails))
sys.exit(1 if fails else 0)
