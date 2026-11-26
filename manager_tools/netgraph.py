#!/usr/bin/env python3
"""GitHub-style network graph of this git repo as a self-updating HTML page.

  netgraph.py build  [--repo PATH] [--out DIR] [--hide-remote NAME ...]   write the page once
  netgraph.py watch  [same options]                                     rebuild whenever a ref changes

The page (DIR/index.html, default manager_tools/out/) re-reads its data every 2 s, so a
browser tab stays current as long as `watch` runs. Use ./manager_tools/netgraph.sh to start both.
"""
import argparse, json, os, subprocess, sys, time, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))


def git(repo, *args):
    return subprocess.run(["git", "-C", repo, *args], check=True,
                          capture_output=True, text=True, errors="replace").stdout


def collect(repo, hide=()):
    excl = [f"--exclude={r}/*" for r in hide]
    raw = git(repo, "log", "--branches", *excl, "--remotes", "--tags", "HEAD", "--date-order", "--reverse",
              "--pretty=format:%H%x1f%P%x1f%ct%x1f%an%x1f%s%x1e")
    commits, index = [], {}
    for rec in raw.split("\x1e"):
        rec = rec.strip("\n")
        if not rec:
            continue
        h, p, t, a, s = rec.split("\x1f")
        index[h] = len(commits)
        commits.append({"h": h, "p": p.split() if p else [], "t": int(t), "a": a, "s": s})

    refs = []  # (name, kind, hash)
    for line in git(repo, "for-each-ref", "--format=%(refname)%1f%(objectname)%1f%(*objectname)",
                    "refs/heads", "refs/remotes", "refs/tags").splitlines():
        name, obj, peeled = line.split("\x1f")
        h = peeled or obj
        if name.endswith("/HEAD") or h not in index or any(name.startswith(f"refs/remotes/{r}/") for r in hide):
            continue
        if name.startswith("refs/heads/"):
            refs.append((name[11:], "local", h))
        elif name.startswith("refs/remotes/"):
            refs.append((name[13:], "remote", h))
        else:
            refs.append((name[10:], "tag", h))
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
    return commits, index, refs, head, head_ref, default


def layout(commits, index, refs, default, head_ref=""):
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
    order = sorted(range(len(rows)), key=lambda r: (rows[r]["name"] != default, rows[r]["kind"] == "anon", start[r]))
    remap = {old: new for new, old in enumerate(order)}
    rows = [rows[o] for o in order]
    for r in rows:
        r.pop("slot", None)
    return [remap[r] for r in row], rows


def build(repo, out, hide=()):
    repo = os.path.abspath(repo)
    commits, index, refs, head, head_ref, default = collect(repo, hide)
    row, rows = layout(commits, index, refs, default, head_ref)
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
    ap.add_argument("--hide-remote", nargs="*", default=[], metavar="NAME", help="remotes to leave out, e.g. origin")
    ap.add_argument("--interval", type=float, default=2.0)
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    t = time.time()
    page, nc, nr = build(a.repo, out, a.hide_remote)
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
                    build(a.repo, out, a.hide_remote)
                    print(time.strftime("%H:%M:%S"), "rebuilt", flush=True)
            except subprocess.CalledProcessError as e:  # e.g. mid-rebase; try again next tick
                print("git error, retrying:", e.stderr.strip()[:200], flush=True)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
