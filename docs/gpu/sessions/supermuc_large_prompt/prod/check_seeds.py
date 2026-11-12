#!/usr/bin/env python3
"""Check that all seeds of one or more ALF seeds files are positive integers and pairwise distinct (across all files).

usage: check_seeds.py SEEDFILE [SEEDFILE ...]       exit 1 on a non-integer, non-positive or repeated seed (python3 stdlib only)
Use it on every RUNDIR/seeds of a campaign (different run directories must not share a seed: identical seeds = identical random streams).
"""
import sys


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    seen, bad, total = {}, [], 0
    for f in sys.argv[1:]:
        n = 0
        with open(f) as fh:
            for ln, line in enumerate(fh, 1):
                for tok in line.split():
                    n += 1
                    try:
                        s = int(tok)
                    except ValueError:
                        bad.append("%s:%d: not an integer: %r" % (f, ln, tok))
                        continue
                    if s <= 0:
                        bad.append("%s:%d: not positive: %d" % (f, ln, s))
                    if s in seen:
                        bad.append("%s:%d: seed %d already used at %s" % (f, ln, s, seen[s]))
                    else:
                        seen[s] = "%s:%d" % (f, ln)
        print("%s: %d seed(s)" % (f, n))
        total += n
    for b in bad:
        print("VIOLATION " + b)
    print("%s: %d seeds in %d file(s), %d distinct, %d problem(s)" % ("FAIL" if bad else "PASS", total, len(sys.argv) - 1, len(seen), len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
