#!/usr/bin/env python3
"""Generate a Slurm benchmark campaign (bash-only job scripts with provenance) from a campaign description.

Usage: python3 -I make_campaign.py DESCRIPTION.py [--root DIR] [--dry-run]

DESCRIPTION.py is a python file that defines a dict CAMPAIGN (see README.md and tests/campaign_test.py):
  id, profile (name in profiles/ or a path), partition, binary (path on the target machine), build_info (default
  <binary>.build_info.txt), binary_sha256 (optional), fakeintel (on|off, default: profile), repetitions, seed,
  max_jobs_in_queue, max_attempts, template (optional parameter template dir), wall_safety, points = [
    {name, model, keys {KEY: VALUE} (make_params.py syntax), L [list], schemes [standard|delayed|submatrix],
     n_delay [ints >= -1; 0 = auto, -1 = auto+tune], nsweep, nbin (int or {L: int}), layouts ["RxT", ...],
     partial_node (allow ranks*threads != cores per node), n_fl, orbitals, sec_per_sweep_est (number or {L: number}),
     setup_s_est (number or {L: number}), mem_gb_per_rank (optional override of the estimate)} ]

Writes <root>/<id>/ (default root: <ALF>/runs/campaigns): campaign.json (resolved description, run list, sha256 of
every script), campaign.env, jobs.list, jobs/<job>.sbatch, runs/<run_id>/input/{parameters,seeds}, scripts/
(profile.sh, provenance.sh, rank_wrap.sh), submit_all.sh, logs/. One job = one point x L x layout x repetition with all
schemes (standard once, every other scheme x n_delay), in an order rotated by the repetition index.

Refuses (exit 1): a profile containing the placeholder marker '[[' anywhere, or a description with a '[[...]]'
placeholder in any string value; fewer than 3 repetitions for a non-local profile; ranks*threads != cores per node
unless partial_node; memory or wall-time estimate above the node / partition limit; an existing campaign directory
with a different description (re-running with the same description is a no-op; nothing is ever deleted).
"""
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ALF = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MAKE_PARAMS = os.path.join(HERE, "..", "make_params.py")
SCHEMES = ("standard", "delayed", "submatrix")
PROFILE_VARS = ("PROFILE_NAME", "CORES_PER_NODE", "NODE_MEM_GB", "PARTITIONS", "ACCOUNT", "SBATCH_EXTRA",
                "FAKEINTEL_DEFAULT")


def die(msg):
    sys.stderr.write("make_campaign.py: " + msg + "\n")
    sys.exit(1)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fortran_float(v):
    return float(str(v).lower().replace("d", "e"))


def per_l(value, L, what, point):
    if isinstance(value, dict):
        if L not in value:
            die("point %s: %s has no entry for L=%d" % (point, what, L))
        return value[L]
    if value is None:
        die("point %s: %s missing (calibrate with a pilot run)" % (point, what))
    return value


def placeholder_check(path, text):
    bad = [(i + 1, l.rstrip()) for i, l in enumerate(text.splitlines()) if "[[" in l]
    if bad:
        die("%s still contains placeholders '[[...]]' (fill them first):\n" % path +
            "\n".join("  line %d: %s" % b for b in bad))


def read_profile(path):
    with open(path) as f:
        placeholder_check(path, f.read())
    script = '. "$1" >/dev/null 2>&1 || exit 3; for v in %s; do eval "printf \'%%s\\037\' \\"\\${$v}\\""; done' % \
        " ".join(PROFILE_VARS)
    out = subprocess.run(["/bin/bash", "-c", script, "_", path], capture_output=True, text=True)
    if out.returncode != 0:
        die("sourcing profile %s failed (rc=%d)" % (path, out.returncode))
    vals = dict(zip(PROFILE_VARS, out.stdout.split("\x1f")))
    prof = {"name": vals["PROFILE_NAME"], "cores_per_node": int(vals["CORES_PER_NODE"]),
            "node_mem_gb": float(vals["NODE_MEM_GB"]), "account": vals["ACCOUNT"],
            "sbatch_extra": [l.strip() for l in vals["SBATCH_EXTRA"].splitlines() if l.strip()],
            "fakeintel_default": vals["FAKEINTEL_DEFAULT"] or "off", "partitions": {}}
    for item in vals["PARTITIONS"].split():
        name, minutes = item.split(":")
        prof["partitions"][name] = int(minutes)
    return prof


PLACEHOLDER = re.compile(r"\[\[[^\]]*\]\]")


def placeholder_values(obj, where="CAMPAIGN"):
    """(location, value) of every string value (also dict keys) containing a [[...]] placeholder."""
    out = []
    if isinstance(obj, str):
        if PLACEHOLDER.search(obj):
            out.append((where, obj))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out += placeholder_values(k, where + " key") + placeholder_values(v, "%s[%r]" % (where, k))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out += placeholder_values(v, "%s[%d]" % (where, i))
    return out


def load_description(path):
    with open(path) as f:
        text = f.read()
    ns = {}
    exec(compile(text, path, "exec"), ns)
    if "CAMPAIGN" not in ns:
        die("%s defines no CAMPAIGN dict" % path)
    bad = placeholder_values(ns["CAMPAIGN"])
    if bad:
        die("%s still contains placeholders '[[...]]':\n" % path + "\n".join("  %s = %r" % b for b in bad))
    return ns["CAMPAIGN"], text


def mem_estimate_gb(ndim, n_fl, beta, dtau, nwrap):
    """GR, work matrices (~12 N x N) plus the UDV stack (2 N x N per stabilisation point) per flavour, complex double."""
    nstab = int(math.ceil(beta / dtau / nwrap)) + 1
    return (16.0 * ndim * ndim * n_fl * (12 + 2 * nstab) + 200e6) / 1e9


def nd_tag(nd):
    return "ndm%d" % (-nd) if nd < 0 else "nd%d" % nd


def build(desc, prof, profile_path):
    cid = desc["id"]
    if not re.match(r"^[A-Za-z0-9][A-Za-z0-9_-]*$", cid):
        die("campaign id %r: use letters, digits, '_' and '-'" % cid)
    part = desc.get("partition")
    if part not in prof["partitions"]:
        die("partition %r not in profile %s (%s)" % (part, prof["name"], " ".join(prof["partitions"])))
    max_wall_min = prof["partitions"][part]
    fakeintel = desc.get("fakeintel", prof["fakeintel_default"])
    if fakeintel not in ("on", "off"):
        die("fakeintel must be on or off")
    reps = int(desc["repetitions"])
    if reps < 3 and prof["name"] != "local":
        die("repetitions = %d: at least 3 are required for a cluster campaign (median and spread)" % reps)
    seed = int(desc["seed"])
    safety = float(desc.get("wall_safety", 1.5))
    prefix = "c_" + cid
    jobs, runs = [], []
    for p in desc["points"]:
        name = p["name"]
        if not re.match(r"^[A-Za-z0-9]+$", name):
            die("point name %r: letters and digits only" % name)
        keys = dict(p.get("keys", {}))
        for k in ("Beta", "Dtau"):
            if k not in keys:
                die("point %s: keys must set %s (provenance and estimates)" % (name, k))
        beta, dtau = fortran_float(keys["Beta"]), fortran_float(keys["Dtau"])
        nwrap = int(keys.get("Nwrap", 10))
        for s in p["schemes"]:
            if s not in SCHEMES:
                die("point %s: unknown scheme %r" % (name, s))
        nds = [int(n) for n in p.get("n_delay", [])]
        if any(n < -1 for n in nds):
            die("point %s: n_delay must be >= -1 (0 = auto, -1 = auto+tune)" % name)
        if any(s != "standard" for s in p["schemes"]) and not nds:
            die("point %s: low-rank schemes need an n_delay list" % name)
        for L in p["L"]:
            L = int(L)
            nsweep = int(per_l(p.get("nsweep"), L, "nsweep", name))
            nbin = int(per_l(p.get("nbin"), L, "nbin", name))
            sps = float(per_l(p.get("sec_per_sweep_est"), L, "sec_per_sweep_est", name))
            setup = float(per_l(p.get("setup_s_est"), L, "setup_s_est", name))
            ndim = L * L * int(p.get("orbitals", 1))
            n_fl = int(p.get("n_fl", 2))
            mem = float(p["mem_gb_per_rank"]) if "mem_gb_per_rank" in p else mem_estimate_gb(ndim, n_fl, beta, dtau, nwrap)
            entries = []
            if "standard" in p["schemes"]:
                entries.append(("standard", "n/a"))
            for s in p["schemes"]:
                if s != "standard":
                    entries.extend((s, nd) for nd in nds)
            for lay in p["layouts"]:
                m = re.match(r"^(\d+)x(\d+)$", lay)
                if not m:
                    die("point %s: layout %r is not RANKSxTHREADS" % (name, lay))
                ranks, threads = int(m.group(1)), int(m.group(2))
                if ranks * threads != prof["cores_per_node"] and not p.get("partial_node", False):
                    die("point %s layout %s: ranks x threads = %d != %d cores per node (set partial_node to allow)"
                        % (name, lay, ranks * threads, prof["cores_per_node"]))
                if ranks * threads > prof["cores_per_node"]:
                    die("point %s layout %s: %d cores > %d per node" % (name, lay, ranks * threads, prof["cores_per_node"]))
                if ranks * mem > 0.9 * prof["node_mem_gb"]:
                    die("point %s L=%d layout %s: memory estimate %d x %.2f GB > 90%% of %.0f GB node memory"
                        % (name, L, lay, ranks, mem, prof["node_mem_gb"]))
                est_s = len(entries) * (setup + nsweep * nbin * sps) + 120.0
                wall_min = int(math.ceil(est_s * safety / 60.0)) + 5
                if wall_min > max_wall_min:
                    die("point %s L=%d layout %s: wall estimate %d min > %d min of partition %s"
                        % (name, L, lay, wall_min, max_wall_min, part))
                for rep in range(1, reps + 1):
                    jname = "%s_%s_L%d_%s_r%d" % (prefix, name, L, lay, rep)
                    shift = (rep - 1) % len(entries)
                    order = entries[shift:] + entries[:shift]
                    jruns = []
                    for s, nd in order:
                        rid = "%s_L%d_%s_r%d_%s" % (name, L, lay, rep, s) + ("" if s == "standard" else "_" + nd_tag(nd))
                        jruns.append(rid)
                        runs.append({"run_id": rid, "job": jname, "point": name, "model": p.get("model", "?"),
                                     "L": L, "ndim": ndim, "beta": beta, "dtau": dtau, "scheme": s,
                                     "n_delay_requested": nd, "layout": lay, "ranks": ranks, "threads": threads,
                                     "repetition": rep, "seed": seed + rep - 1, "nsweep": nsweep, "nbin": nbin,
                                     "keys": keys})
                    jobs.append({"job": jname, "point": name, "L": L, "layout": lay, "ranks": ranks,
                                 "threads": threads, "repetition": rep, "runs": jruns, "wall_min": wall_min,
                                 "est_s": est_s, "mem_gb_per_rank": mem})
    ids = [r["run_id"] for r in runs]
    if len(set(ids)) != len(ids):
        die("duplicate run ids (repeated L, layout or n_delay in a point?)")
    return {"prefix": prefix, "fakeintel": fakeintel, "max_wall_min": max_wall_min, "jobs": jobs, "runs": runs}


def hhmm(minutes):
    return "%d:%02d:00" % (minutes // 60, minutes % 60)


def main(argv):
    args, root, dry = [], os.path.join(ALF, "runs", "campaigns"), False
    i = 1
    while i < len(argv):
        if argv[i] == "--root":
            root = argv[i + 1]
            i += 1
        elif argv[i] == "--dry-run":
            dry = True
        else:
            args.append(argv[i])
        i += 1
    if len(args) != 1:
        sys.exit(__doc__)
    desc_path = os.path.abspath(args[0])
    desc, desc_text = load_description(desc_path)
    profile_path = desc["profile"]
    if os.sep not in profile_path:
        profile_path = os.path.join(HERE, "profiles", profile_path + ".sh")
    prof = read_profile(profile_path)
    plan = build(desc, prof, profile_path)
    binary = desc["binary"]
    build_info = desc.get("build_info", binary + ".build_info.txt")
    bin_sha = desc.get("binary_sha256", "")
    if not (os.path.isabs(binary) and os.path.isabs(build_info)):
        die("binary and build_info must be absolute paths on the target machine (jobs run in other directories)")
    if prof["name"] == "local":
        if not os.path.isfile(binary) or not os.path.isfile(build_info):
            die("local profile: binary %s or build info %s missing (run build_info.sh)" % (binary, build_info))
        recorded = [l.split("=", 1)[1].strip() for l in open(build_info) if l.startswith("binary_sha256=")]
        if not recorded or recorded[0] != sha256(binary):
            die("local profile: sha256 of %s does not match %s" % (binary, build_info))
        bin_sha = bin_sha or recorded[0]

    scripts = {"make_campaign.py": os.path.abspath(__file__), "job_template.sh": os.path.join(HERE, "job_template.sh"),
               "provenance.sh": os.path.join(HERE, "provenance.sh"), "rank_wrap.sh": os.path.join(HERE, "rank_wrap.sh"),
               "submit_all.sh": os.path.join(HERE, "submit_all.sh"), "profile.sh": profile_path,
               "make_params.py": MAKE_PARAMS, "description": desc_path}
    resolved = {"campaign_id": desc["id"], "profile": prof, "partition": desc["partition"], "binary": binary,
                "build_info": build_info, "binary_sha256": bin_sha, "fakeintel": plan["fakeintel"],
                "repetitions": desc["repetitions"], "max_jobs_in_queue": int(desc.get("max_jobs_in_queue", 10)),
                "max_attempts": int(desc.get("max_attempts", 2)), "job_prefix": plan["prefix"],
                "script_sha256": {k: sha256(v) for k, v in sorted(scripts.items())},
                "description_text": desc_text, "jobs": plan["jobs"], "runs": plan["runs"]}
    core_h = sum(j["est_s"] for j in plan["jobs"]) * prof["cores_per_node"] / 3600.0
    core_h_lim = sum(j["wall_min"] for j in plan["jobs"]) * prof["cores_per_node"] / 60.0
    summary = ("campaign %s: %d jobs, %d runs, profile %s, partition %s; estimate %.1f core-h (wall limits: %.1f core-h, "
               "whole nodes of %d cores)" % (desc["id"], len(plan["jobs"]), len(plan["runs"]), prof["name"],
                                             desc["partition"], core_h, core_h_lim, prof["cores_per_node"]))
    if dry:
        print(summary)
        return 0

    # job scripts (in memory first: their sha256 is part of campaign.json, which decides idempotency)
    with open(os.path.join(HERE, "job_template.sh")) as f:
        tmpl = f.read()
    sbatch_lines = ["#SBATCH --partition=%s" % desc["partition"]] if prof["name"] != "local" else []
    if prof["account"]:
        sbatch_lines.append("#SBATCH --account=%s" % prof["account"])
    sbatch_lines += ["#SBATCH " + l for l in prof["sbatch_extra"]]
    byid = {r["run_id"]: r for r in plan["runs"]}
    job_text = {}
    for j in plan["jobs"]:
        runlines = "\n".join("%s %s %s" % (rid, byid[rid]["scheme"], byid[rid]["n_delay_requested"]) for rid in j["runs"])
        subst = {"JOB_NAME": j["job"], "RANKS": str(j["ranks"]), "THREADS": str(j["threads"]),
                 "WALL": hhmm(j["wall_min"]), "CAMPAIGN_ID": desc["id"],
                 "FAKEINTEL": plan["fakeintel"], "BINARY": binary, "BINARY_SHA256": bin_sha, "BUILD_INFO": build_info,
                 "MAX_ATTEMPTS": str(resolved["max_attempts"]), "RUNS": runlines}
        text = tmpl.replace("@@SBATCH_LINES@@\n", "".join(l + "\n" for l in sbatch_lines))
        for k, v in subst.items():
            text = text.replace("@@%s@@" % k, v)
        if "@@" in text:
            die("job template: unsubstituted token in " + j["job"])
        job_text[j["job"]] = text
        j["script_sha256"] = hashlib.sha256(text.encode()).hexdigest()

    C = os.path.join(os.path.abspath(root), desc["id"])
    cj = os.path.join(C, "campaign.json")
    new_json = json.dumps(resolved, indent=1, sort_keys=True) + "\n"
    if os.path.exists(C):
        if os.path.exists(cj) and open(cj).read() == new_json:
            print("campaign %s already generated in %s (unchanged)" % (desc["id"], C))
            return 0
        die("%s exists with a different description or scripts: use a new campaign id (nothing is overwritten)" % C)
    for d in ("jobs", "runs", "scripts", "logs"):
        os.makedirs(os.path.join(C, d))
    for name, src in (("provenance.sh", scripts["provenance.sh"]), ("rank_wrap.sh", scripts["rank_wrap.sh"]),
                      ("profile.sh", profile_path)):
        shutil.copyfile(src, os.path.join(C, "scripts", name))
    shutil.copyfile(scripts["submit_all.sh"], os.path.join(C, "submit_all.sh"))
    shutil.copyfile(desc_path, os.path.join(C, "description.py"))
    for f in ("scripts/rank_wrap.sh", "submit_all.sh"):
        os.chmod(os.path.join(C, f), 0o755)

    template = desc.get("template")
    for r in plan["runs"]:
        inp = os.path.join(C, "runs", r["run_id"], "input")
        kv = ["L1=%d" % r["L"], "L2=%d" % r["L"]] + ["%s=%s" % kv for kv in r["keys"].items()]
        kv += ["NSweep=%d" % r["nsweep"], "NBin=%d" % r["nbin"], 'Sequential_scheme="%s"' % r["scheme"],
               "SEED=%d" % r["seed"]]
        if r["scheme"] != "standard":          # standard ignores N_delay: the template value stays
            kv.append("N_delay=%d" % r["n_delay_requested"])
        cmd = [sys.executable, "-I", MAKE_PARAMS] + (["--template", template] if template else []) + [inp] + kv
        subprocess.run(cmd, check=True)
    for name, text in job_text.items():
        path = os.path.join(C, "jobs", name + ".sbatch")
        with open(path, "w") as f:
            f.write(text)
        os.chmod(path, 0o755)
    with open(os.path.join(C, "jobs.list"), "w") as f:
        for j in plan["jobs"]:
            f.write("%s %s\n" % (j["job"], " ".join(j["runs"])))
    with open(os.path.join(C, "campaign.env"), "w") as f:
        f.write("CAMPAIGN_ID=%s\nJOB_PREFIX=%s\nMAX_JOBS_IN_QUEUE=%d\nMAX_ATTEMPTS=%d\n"
                % (desc["id"], plan["prefix"], resolved["max_jobs_in_queue"], resolved["max_attempts"]))
    with open(cj, "w") as f:
        f.write(new_json)
    print(summary)
    print("generated %s; submit with %s/submit_all.sh" % (C, C))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
