#!/bin/bash
# Unit test of recover.sh restore/backup logic with a fake h5ls and a cp that fails on purpose (disk full).
# Usage: bash recover_unit_test.sh   (runs in a temporary directory, needs no HDF5)
# Added 2026-11-16 after the quota crash: restore must not fall back to an older source when a COPY fails,
# and rc_backup must keep the old .segbak until the new one is complete.
R=$(cd "$(dirname "$0")" && pwd)/recover.sh
S=$(mktemp -d); mkdir -p $S/bin
cat > $S/bin/h5ls <<'H5'
#!/bin/bash
f=${@: -1}; c=$(cat "$f" 2>/dev/null) || exit 1
case "$c" in bins=*) echo "/Ener_scal/sign          Dataset {${c#bins=}/Inf}";; conf) echo ok;; *) exit 1;; esac
H5
chmod +x $S/bin/h5ls; export PATH=$S/bin:$PATH
cd $S
pass=0; fail=0; ok() { if eval "$2"; then echo "PASS $1"; pass=$((pass+1)); else echo "FAIL $1"; fail=$((fail+1)); fi; }
setup() {  # one set Temp_0, 2 walkers
  rm -rf run; mkdir -p run/Temp_0/.snap/A run/Temp_0/.snap/B run/Temp_0/.segbak; cd run
  echo "mpi_per_parameter_set = 2" > parameters; : > RUNNING
  echo BAD > Temp_0/data.h5; echo conf > Temp_0/confout_0.h5; echo conf > Temp_0/confout_1.h5
  echo bins=100 > Temp_0/.snap/A/data.h5; echo conf > Temp_0/.snap/A/confout_0.h5; echo conf > Temp_0/.snap/A/confout_1.h5
  echo "bins=100 time=2 seq=2 walkers=2" > Temp_0/.snap/A/COMPLETE
  echo bins=95 > Temp_0/.snap/B/data.h5; echo conf > Temp_0/.snap/B/confout_0.h5; echo conf > Temp_0/.snap/B/confout_1.h5
  echo "bins=95 time=1 seq=1 walkers=2" > Temp_0/.snap/B/COMPLETE
  echo bins=50 > Temp_0/.segbak/data.h5; echo conf > Temp_0/.segbak/confin_0.h5; echo conf > Temp_0/.segbak/confin_1.h5; : > Temp_0/.segbak/COMPLETE
  source $R
}
cd $S
# T1: every copy of data.h5 fails (disk full) -> FATAL, no fallback, live + sources unchanged
setup; cp() { case "$*" in *data.h5*) return 1;; *) command cp "$@";; esac; }
rc_prepare > ../t1.log 2>&1; rc=$?; unset -f cp
ok "T1 prepare refuses" '[ $rc -ne 0 ]'
ok "T1 no fallback logged" 'grep -q "NOT falling back" ../t1.log && ! grep -q "restored from" ../t1.log'
ok "T1 live untouched" '[ "$(cat Temp_0/data.h5)" = BAD ]'
ok "T1 segbak untouched" '[ "$(cat Temp_0/.segbak/data.h5)" = bins=50 ] && [ -f Temp_0/.segbak/COMPLETE ]'
ok "T1 no rtmp left" '! ls Temp_0/*.rtmp >/dev/null 2>&1'
cd ..
# T2: snapshot A invalid (marker mismatch) -> next valid source B used
setup; echo "bins=999 time=2 seq=2 walkers=2" > Temp_0/.snap/A/COMPLETE
rc_prepare > ../t2.log 2>&1; rc=$?
ok "T2 prepare ok" '[ $rc -eq 0 ]'
ok "T2 restored from B" '[ "$(cat Temp_0/data.h5)" = bins=95 ] && grep -q "snapshot-B" ../t2.log'
ok "T2 new segbak = restored state" '[ "$(cat Temp_0/.segbak/data.h5)" = bins=95 ] && [ -f Temp_0/.segbak/COMPLETE ] && [ ! -e Temp_0/.segbak.new ]'
cd ..
# T3: normal restore from A, then backup copy fails -> old segbak kept, .segbak.new removed
setup; cp() { case "$*" in *".segbak.new"*) return 1;; *) command cp "$@";; esac; }
rc_prepare > ../t3.log 2>&1; rc=$?; unset -f cp
ok "T3 prepare refuses" '[ $rc -ne 0 ]'
ok "T3 live restored from A" '[ "$(cat Temp_0/data.h5)" = bins=100 ]'
ok "T3 old segbak kept" '[ "$(cat Temp_0/.segbak/data.h5)" = bins=50 ] && [ -f Temp_0/.segbak/COMPLETE ]'
ok "T3 no .segbak.new" '[ ! -e Temp_0/.segbak.new ]'
cd ..
# T4: no RUNNING (clean start) -> backup replaces segbak
setup; rm RUNNING; echo bins=120 > Temp_0/data.h5; rm Temp_0/confout_*; echo conf > Temp_0/confin_0.h5; echo conf > Temp_0/confin_1.h5
rc_prepare > ../t4.log 2>&1; rc=$?
ok "T4 ok + segbak replaced" '[ $rc -eq 0 ] && [ "$(cat Temp_0/.segbak/data.h5)" = bins=120 ] && [ -f Temp_0/.segbak/COMPLETE ] && [ $(ls Temp_0/.segbak | wc -l) -eq 4 ]'
cd ..
echo "pass=$pass fail=$fail"
