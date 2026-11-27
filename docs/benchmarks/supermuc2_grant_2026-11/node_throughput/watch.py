#!/usr/bin/env python3
"""Log every change of a file's mtime (ns) and size while a process lives.

usage: watch.py FILE PID > watch.txt
Output lines: mtime_ns size wall_ns. With NSweep = 1, rank 0 writes one bin per sweep; one bin is a short burst of
HDF5 writes, so analyse.py groups changes closer than 50 ms into one bin. Poll every 2 ms (negligible cpu).
"""
import os
import sys
import time

path, pid = sys.argv[1], int(sys.argv[2])
last = None
while True:
    try:
        os.kill(pid, 0)
    except OSError:
        break
    try:
        st = os.stat(path)
        cur = (st.st_mtime_ns, st.st_size)
        if cur != last:
            print(cur[0], cur[1], time.time_ns(), flush=True)
            last = cur
    except FileNotFoundError:
        pass
    time.sleep(0.002)
