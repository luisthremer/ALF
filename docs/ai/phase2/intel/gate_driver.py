#!/usr/bin/env python3
"""Runs gate.py cases on the 8 tiles of one node (one worker per tile, ZE_AFFINITY_MASK = tile).
usage: gate_driver.py CASEFILE OUTROOT   (CASEFILE lines: name exe start key=value ...)"""
import os, subprocess, sys, threading, queue
cases, root = sys.argv[1:3]
here = os.path.dirname(os.path.abspath(__file__))
q = queue.Queue()
for l in open(cases):
    if l.strip() and not l.startswith("#"): q.put(l.split())
lock = threading.Lock()
def worker(t):
    while True:
        try: c = q.get_nowait()
        except queue.Empty: return
        env = dict(os.environ, ZE_AFFINITY_MASK=str(t))
        r = subprocess.run([sys.executable, os.path.join(here, "gate.py"), os.path.join(root, c[0])] + c[1:], env=env,
                           capture_output=True, text=True)
        with lock: print((r.stdout.strip() or ("FAIL %s %s" % (c[0], r.stderr.strip()[-300:]))), flush=True)
ts = [threading.Thread(target=worker, args=(t,)) for t in range(8)]
[t.start() for t in ts]; [t.join() for t in ts]
