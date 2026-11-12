source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P
for m in "$@"; do echo "######## $m"; ./bin/qr2_test $m ${NB:-16} 2>&1 | grep -v -E "^warning|^\s*$" ; done
