#!/usr/bin/env python3
"""GitHub-style network graph of this git repo as a self-updating HTML page.

  netgraph.py build  [--repo PATH] [--out DIR] [--hide-remote NAME ...] [--only PATTERN ...]
  netgraph.py watch  [same options]                                     rebuild whenever a ref changes

The page (DIR/index.html, default manager_tools/out/) re-reads its data every 2 s, so a
browser tab stays current as long as `watch` runs. Defaults for --hide-remote and --only
can be set in netgraph.json next to this script ({"hide_remotes": [...], "only": [...], "row_tags": [...], "show_unnamed": false}). Use ./manager_tools/netgraph.sh to start both.
"""
import argparse, fnmatch, json, os, subprocess, sys, time, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))


def git(repo, *args):
    return subprocess.run(["git", "-C", repo, *args], check=True,
                          capture_output=True, text=True, errors="replace").stdout


def matches(name, kind, only):
    """Branch filter: a pattern matches the name, or the name without its remote prefix."""
    if not only:
        return True
    short = name.split("/", 1)[1] if kind == "remote" and "/" in name else name
    return any(fnmatch.fnmatchcase(name, p) or fnmatch.fnmatchcase(short, p) for p in only)


def collect(repo, hide=(), only=()):
    allrefs = []  # (name, kind, hash, full refname)
    for line in git(repo, "for-each-ref", "--format=%(refname)%1f%(objectname)%1f%(*objectname)",
                    "refs/heads", "refs/remotes", "refs/tags").splitlines():
        full, obj, peeled = line.split("\x1f")
        if full.endswith("/HEAD") or any(full.startswith(f"refs/remotes/{r}/") for r in hide):
            continue
        if full.startswith("refs/heads/"):
            allrefs.append((full[11:], "local", peeled or obj, full))
        elif full.startswith("refs/remotes/"):
            allrefs.append((full[13:], "remote", peeled or obj, full))
        else:
            allrefs.append((full[10:], "tag", peeled or obj, full))
    if only:  # only the selected branches and their history; tags only label what is shown
        allrefs = [r for r in allrefs if r[1] == "tag" or matches(r[0], r[1], only)]
        starts = [r[3] for r in allrefs if r[1] != "tag"]
    else:
        starts = ["--branches", *[f"--exclude={r}/*" for r in hide], "--remotes", "--tags", "HEAD"]
    raw = git(repo, "log", *starts, "--date-order", "--reverse",
              "--pretty=format:%H%x1f%P%x1f%ct%x1f%an%x1f%s%x1e") if starts else ""
    commits, index = [], {}
    for rec in raw.split("\x1e"):
        rec = rec.strip("\n")
        if not rec:
            continue
        h, p, t, a, s = rec.split("\x1f")
        index[h] = len(commits)
        commits.append({"h": h, "p": p.split() if p else [], "t": int(t), "a": a, "s": s})

    refs = [(n, k, h) for n, k, h, _ in allrefs if h in index]
    try:
        head = git(repo, "rev-parse", "HEAD").strip()
        head_ref = git(repo, "symbolic-ref", "-q", "--short", "HEAD").strip()
    except subprocess.CalledProcessError:
        head, head_ref = "", ""
    try:
        default = git(repo, "symbolic-ref", "--short", "refs/remotes/origin/HEAD").strip()
    except subprocess.CalledProcessError:
        default = ""
    names = {n for n, _, _ in refs}
    if default not in names:
        default = next((n for n in ("master", "main", "fork/master") if n in names), "")
    if only and not (head_ref and matches(head_ref, "local", only)):
        head = ""  # checked-out branch is filtered out: open at the newest commit instead
    return commits, index, refs, head, head_ref, default


def layout(commits, index, refs, default, head_ref="", row_tags=()):
    """Give every commit a row. Each branch claims its first-parent chain; leftovers get packed rows."""
    n = len(commits)
    row = [-1] * n
    rows = []  # {"name", "kind", "refs": [...]}

    def prio(r):
        name, kind, h = r
        return (0 if name == default else 1 if name == head_ref else 2, {"local": 0, "remote": 1}.get(kind, 2), -index[h], name)

    for name, kind, h in sorted((r for r in refs if r[1] != "tag"), key=prio):
        i = index[h]
        if row[i] != -1:
            continue
        rid = len(rows)
        rows.append({"name": name, "kind": kind})
        while i is not None and row[i] == -1:
            row[i] = rid
            ps = commits[i]["p"]
            i = index.get(ps[0]) if ps else None

    # leftover commits (only reachable via merges or tags): chains, packed into shared rows
    children = [[] for _ in range(n)]
    for i, c in enumerate(commits):
        for p in c["p"]:
            if p in index:
                children[index[p]].append(i)
    tag_at = {}
    for name, kind, h in refs:
        if kind == "tag" and any(fnmatch.fnmatchcase(name, p) for p in row_tags):
            tag_at.setdefault(index[h], name)
    anon_start = []  # per anonymous row: leftmost x used (chains arrive right to left)
    for i in range(n - 1, -1, -1):
        if row[i] != -1:
            continue
        chain, j = [], i
        while j is not None and row[j] == -1:
            chain.append(j)
            row[j] = -2
            ps = commits[j]["p"]
            j = index.get(ps[0]) if ps else None
        lo, hi = min(chain), max(chain)
        tag = tag_at.get(chain[0])
        if tag:  # merged branch that was archived as a tag: its own row, named after the tag
            rid = len(rows)
            rows.append({"name": tag, "kind": "tag"})
            for j in chain:
                row[j] = rid
            continue
        lo_ext = index.get(commits[lo]["p"][0], lo) if commits[lo]["p"] else lo
        hi_ext = max([hi] + children[hi])
        slot = next((k for k, st in enumerate(anon_start) if hi_ext < st - 1), None)
        if slot is None:
            slot = len(anon_start)
            anon_start.append(n + 1)
            rows.append({"name": "", "kind": "anon", "slot": slot})
        anon_start[slot] = min(anon_start[slot], lo_ext)
        rid = next(r for r, x in enumerate(rows) if x["kind"] == "anon" and x.get("slot") == slot)
        for j in chain:
            row[j] = rid

    # order rows: default first, then by where each row starts (staircase like GitHub)
    start = [n] * len(rows)
    for i, r in enumerate(row):
        start[r] = min(start[r], i)
    order = sorted(range(len(rows)), key=lambda r: (rows[r]["name"] != default, {"anon": 2, "tag": 1}.get(rows[r]["kind"], 0), start[r]))
    remap = {old: new for new, old in enumerate(order)}
    rows = [rows[o] for o in order]
    for r in rows:
        r.pop("slot", None)
    return [remap[r] for r in row], rows


def build(repo, out, hide=(), only=(), row_tags=(), show_unnamed=True):
    repo = os.path.abspath(repo)
    commits, index, refs, head, head_ref, default = collect(repo, hide, only)
    row, rows = layout(commits, index, refs, default, head_ref, row_tags)
    labels = {}
    for name, kind, h in refs:
        labels.setdefault(index[h], []).append([name, kind])
    data = {
        "repo": repo_name(repo),
        "generated": int(time.time()),
        "head": index.get(head, -1), "headRef": head_ref, "default": default,
        "rows": rows,
        # compact commit records: [hash, row, time, author, subject, [parent indices]]
        "c": [[c["h"][:10], row[i], c["t"], c["a"], c["s"], [index[p] for p in c["p"] if p in index]]
              for i, c in enumerate(commits)],
        "labels": labels,
        "showUnnamed": show_unnamed,
    }
    blob = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    ver = hashlib.sha1(blob.encode()).hexdigest()[:12]
    os.makedirs(out, exist_ok=True)

    def write(name, text):
        tmp = os.path.join(out, name + ".tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, os.path.join(out, name))

    write("data.js", f"window.NG_DATA={blob};window.NG_DATA.version={json.dumps(ver)};")
    write("version.js", f"window.NG_VERSION={json.dumps(ver)};")
    with open(os.path.join(HERE, "page.html"), encoding="utf-8") as f:
        write("index.html", f.read())
    return os.path.join(out, "index.html"), len(commits), len(rows)


def repo_name(repo):
    common = os.path.abspath(os.path.join(repo, git(repo, "rev-parse", "--git-common-dir").strip()))
    if os.path.basename(common) == ".git":
        return os.path.basename(os.path.dirname(common))
    return os.path.basename(git(repo, "rev-parse", "--show-toplevel").strip())


def ref_state(repo):
    return git(repo, "for-each-ref", "--format=%(objectname) %(refname)") + git(repo, "rev-parse", "HEAD", "--symbolic-full-name", "HEAD")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["build", "watch"])
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    cfg = {}
    cfg_path = os.path.join(HERE, "netgraph.json")
    if os.path.exists(cfg_path):
        with open(cfg_path) as f:
            cfg = json.load(f)
    ap.add_argument("--hide-remote", nargs="*", default=cfg.get("hide_remotes", []), metavar="NAME",
                    help="remotes to leave out, e.g. origin")
    ap.add_argument("--only", nargs="*", default=cfg.get("only", []), metavar="PATTERN",
                    help="show only branches matching these patterns (e.g. master 'cm/*'); default: all")
    ap.add_argument("--row-tags", nargs="*", default=cfg.get("row_tags", []), metavar="PATTERN",
                    help="merged, deleted branches whose tip carries a matching tag get their own row named after it")
    ap.add_argument("--hide-unnamed", action="store_true", default=not cfg.get("show_unnamed", True),
                    help="start with the rows of merged/deleted branches hidden (toggle in the page)")
    ap.add_argument("--interval", type=float, default=2.0)
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    t = time.time()
    page, nc, nr = build(a.repo, out, a.hide_remote, a.only, a.row_tags, not a.hide_unnamed)
    print(f"{page}  ({nc} commits, {nr} rows, {time.time() - t:.1f}s)", flush=True)
    if a.cmd == "build":
        return
    print("watching refs, Ctrl-C to stop", flush=True)
    state = ref_state(a.repo)
    try:
        while True:
            time.sleep(a.interval)
            try:
                now = ref_state(a.repo)
                if now != state:
                    state = now
                    build(a.repo, out, a.hide_remote, a.only, a.row_tags, not a.hide_unnamed)
                    print(time.strftime("%H:%M:%S"), "rebuilt", flush=True)
            except subprocess.CalledProcessError as e:  # e.g. mid-rebase; try again next tick
                print("git error, retrying:", e.stderr.strip()[:200], flush=True)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
