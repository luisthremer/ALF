#!/bin/bash
# usage: q.sh submit NAME script-file   (drops the script into the daemon's queue)
#        q.sh out NAME [lines]           (prints queue/NAME.out)
#        q.sh wait NAME [seconds]        (waits for queue/NAME.done)
F='grep -v -E "module|Unload|WARNING: Module|ERROR: Unload|EAR|SLURM_EAR|^$|spack/release"'
case $1 in
  submit) timeout 60 ssh supermuc2 "cat > \$SCRATCH/large_prompt/team/p5/queue/$2.tmp && mv \$SCRATCH/large_prompt/team/p5/queue/$2.tmp \$SCRATCH/large_prompt/team/p5/queue/$2.sh" < $3 2>&1 | eval $F ;;
  out) timeout 60 ssh supermuc2 "tail -n ${3:-200} \$SCRATCH/large_prompt/team/p5/queue/$2.out" 2>&1 | eval $F ;;
  wait) for i in $(seq 1 $(( ${3:-600} / 10 ))); do timeout 60 ssh supermuc2 "test -f \$SCRATCH/large_prompt/team/p5/queue/$2.done" 2>/dev/null && { echo done; exit 0; }; sleep 10; done; echo notyet ;;
esac
