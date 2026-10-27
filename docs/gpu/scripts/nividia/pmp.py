"""Poor man's sampling profiler, run inside gdb (gdb-oneapi on nividia; perf/VTune are blocked).

usage: PMP_OUT=stacks.txt PMP_INTERVAL=0.05 gdb-oneapi -batch -x pmp.py --args ./ALF.out
gdb launches the program as its child (allowed with yama ptrace_scope=1), an external shell loop stops it
with SIGINT (not passed to the program) every PMP_INTERVAL seconds (a Python thread cannot run while
gdb.execute("continue") holds the GIL), and the main thread's stack is
written as one folded line "main;...;leaf" (inline frames included, needs -g). Count per unique
stack, same format as flamegraph.pl input. A SIGINT that arrives while the program is already stopped
is delivered right after "continue"; such stops (< interval/4 after resuming) are skipped, not counted.
"""
import collections
import os
import subprocess
import time

import gdb

interval = float(os.environ.get("PMP_INTERVAL", "0.05"))
out_path = os.environ.get("PMP_OUT", "stacks.txt")

for cmd in ("set pagination off", "set confirm off", "set print thread-events off",
            "set print inferior-events off", "handle SIGINT stop print nopass",
            "handle SIGPIPE SIGUSR1 SIGUSR2 SIGCHLD SIGALRM nostop noprint pass"):
    gdb.execute(cmd, to_string=True)



def stack():
    names = []
    f = gdb.newest_frame()
    while f is not None:
        names.append(f.name() or f"?{f.pc():#x}")
        if names[-1] in ("main", "MAIN__"):
            break
        try:
            f = f.older()
        except gdb.error:
            names.append("<unwind-error>")
            break
    return ";".join(reversed(names))


counts = collections.Counter()
gdb.execute("starti", to_string=True)
pid = gdb.selected_inferior().pid
ticker = subprocess.Popen(["bash", "-c", f"while kill -0 {pid} 2>/dev/null; do sleep {interval}; "
                           f"kill -INT {pid} 2>/dev/null; done"])
t0 = time.time()
skipped = 0
while True:
    resumed = time.time()
    try:
        gdb.execute("continue", to_string=True)
    except gdb.error:
        break
    if gdb.selected_inferior().pid == 0:
        break
    if time.time() - resumed < interval / 4:
        skipped += 1
        continue
    try:
        gdb.execute("thread 1", to_string=True)
        counts[stack()] += 1
    except gdb.error as e:
        counts[f"<error {e}>"] += 1
ticker.wait()
with open(out_path, "w", encoding="utf-8") as fh:
    for s, n in counts.most_common():
        fh.write(f"{s} {n}\n")
print(f"pmp: {sum(counts.values())} samples ({skipped} skipped) in {time.time() - t0:.1f} s -> {out_path}")
