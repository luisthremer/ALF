#!/bin/bash
# Ticket 207: precision table of runs/BIN/CONFIG/info.  usage: t207_prectab.sh CONFIG BIN [BIN ...]
W=$HOME/alf-gpu/work/t207
C=$1; shift
printf '%-10s %-12s %-12s %-12s %-12s %-10s %-9s\n' bin GreenMean GreenMax tauMean tauMax PhaseMax CPUtime
for b in "$@"; do
  f=$W/runs/$b/$C/info
  python3 - "$f" "$b" <<'PY'
import re, sys
t = open(sys.argv[1]).read().replace("\n", " ")
def g(label, n):
    m = re.search(label + r"[^:]*:\s+" + r"\s+".join([r"([-+0-9.eEdD]+)"] * n), t)
    return m.groups() if m else ("nan",) * n
gm, gx = g("Precision Green", 2)
tm, tx = g("Precision tau", 2)
(ph,) = g("Precision Phase", 1)
(cp,) = g("CPU Time", 1)
f = lambda x: "%.3e" % float(x.replace("D", "E"))
print("%-10s %-12s %-12s %-12s %-12s %-10s %-9s" % (sys.argv[2], f(gm), f(gx), f(tm), f(tx), f(ph), "%.0fs" % float(cp)))
PY
done
