#!/bin/bash
# harness for the run.sh / chain.sh fixes of 207-3f (no GPU, no real gpu.lock: fake HOME with its own gpu.lock); run on nividia
F=$HOME/alf-gpu/work/t207-3f
R=$F/tmp/rt; FH=$R/fh
export TMPDIR=$F/tmp MINFREE_MB=1
rm -rf -- "$R"; mkdir -p $R/bin $R/runs/dev1/hub12x12_cb_p $R/logs $R/stubproc $FH/alf-gpu
ln -s $HOME/alf-gpu/env $FH/alf-gpu/env; ln -s $HOME/alf-gpu/venv $FH/alf-gpu/venv
printf '1\n2\n3\n' > $R/runs/dev1/hub12x12_cb_p/seeds; echo "ListOfParameters" > $R/runs/dev1/hub12x12_cb_p/parameters
printf '#!/bin/bash\necho fake run\n' > $R/bin/ALF_gnu_fake.out; chmod +x $R/bin/ALF_gnu_fake.out; sha256sum $R/bin/ALF_gnu_fake.out > $R/bin/ALF_gnu_fake.out.sha256
printf '#!/bin/bash\necho "stub proc $1 $2" >> %s/logs/stub_calls.txt; echo "{}" > $1/qref_proc.json\n' $R > $R/stubproc/t207_d5_proc.sh
mk() { # TAG: a run dir as left by a disk wait (inputs + stamp + flag)
  d=$R/runs/$1/hub12x12_cb_p; mkdir -p $R/runs/$1; cp -r $R/runs/dev1/hub12x12_cb_p $d; echo "stamp" > $d/stamp.txt; : > $d/.diskwait; mkdir -p $d/qref; }
run() { # label, then env assignments before the script
  HOME=$FH T207_3W=$R "$@"; }
chk() { d=$R/runs/$1/hub12x12_cb_p; echo "   -> flag present: $([ -e $d/.diskwait ] && echo YES || echo no); provenance rc: $(grep -o '"rc": [0-9]*' $d/provenance.json 2>/dev/null || echo none)"; }

echo "== T1 flag + missing binary (new run.sh)"
mk t1; HOME=$FH T207_3W=$R SET=dev1 TAG=t1 MINFREE_MB=1 bash $F/t207_3_run.sh gnu_missing hub12x12_cb_p; echo "   rc=$?"; chk t1
echo "== T2 flag + run.log present"
mk t2; echo old > $R/runs/t2/hub12x12_cb_p/run.log
HOME=$FH T207_3W=$R SET=dev1 TAG=t2 MINFREE_MB=1 bash $F/t207_3_run.sh gnu_fake hub12x12_cb_p; echo "   rc=$?"; chk t2; echo "   run.log intact: $(cat $R/runs/t2/hub12x12_cb_p/run.log)"
echo "== T3 flag + provenance.json present (must stay untouched)"
mk t3; echo '{"rc": 0, "old": true}' > $R/runs/t3/hub12x12_cb_p/provenance.json; a=$(sha256sum < $R/runs/t3/hub12x12_cb_p/provenance.json)
HOME=$FH T207_3W=$R SET=dev1 TAG=t3 MINFREE_MB=1 bash $F/t207_3_run.sh gnu_fake hub12x12_cb_p; echo "   rc=$?"; chk t3; [ "$a" = "$(sha256sum < $R/runs/t3/hub12x12_cb_p/provenance.json)" ] && echo "   provenance unchanged"
echo "== T4 flag + seeds differ from the source set"
mk t4; echo 9 >> $R/runs/t4/hub12x12_cb_p/seeds
HOME=$FH T207_3W=$R SET=dev1 TAG=t4 MINFREE_MB=1 bash $F/t207_3_run.sh gnu_fake hub12x12_cb_p; echo "   rc=$?"; chk t4
echo "== T4b flag + parameters differ"
mk t4b; echo x >> $R/runs/t4b/hub12x12_cb_p/parameters
HOME=$FH T207_3W=$R SET=dev1 TAG=t4b MINFREE_MB=1 bash $F/t207_3_run.sh gnu_fake hub12x12_cb_p; echo "   rc=$?"; chk t4b
echo "== T5 valid resume (control: must still work)"
mk t5
HOME=$FH T207_3W=$R SET=dev1 TAG=t5 MINFREE_MB=1 bash $F/t207_3_run.sh gnu_fake hub12x12_cb_p | tail -1; chk t5; cat $R/runs/t5/hub12x12_cb_p/run.log

echo "== T6 chain, HOLDLOCK=1, flag + missing binary: OLD scripts (t207-3e), timeout 20 s"
O=$R/old; mkdir -p $O/logs $O/bin $O/runs; cp $HOME/alf-gpu/work/t207-3e/t207_3_{chain,run}.sh $O/; cp -r $R/runs/dev1 $O/runs/; cp $R/bin/* $O/bin/
d=$O/runs/f_def_1em9_dev1/hub12x12_cb_p; mkdir -p $O/runs/f_def_1em9_dev1; cp -r $R/runs/dev1/hub12x12_cb_p $d; echo stamp > $d/stamp.txt; : > $d/.diskwait
( cd $O; HOME=$FH T207_3W=$O HOLDLOCK=1 PROC_DIR=$R/stubproc BD=gnu_missing timeout 20 bash $O/t207_3_chain.sh flt def dev1 1e-9 > /dev/null 2>&1; echo "   old chain rc=$? (124 = still looping at the timeout)"; tail -3 $O/logs/chain_flt.txt | cut -c1-150 )
echo "== T6 same with the NEW scripts"
d=$R/runs/f_def_1em9_dev1/hub12x12_cb_p; mkdir -p $R/runs/f_def_1em9_dev1; cp -r $R/runs/dev1/hub12x12_cb_p $d; echo stamp > $d/stamp.txt; : > $d/.diskwait
cp $F/t207_3_chain.sh $F/t207_3_run.sh $R/
( cd $R; HOME=$FH T207_3W=$R HOLDLOCK=1 PROC_DIR=$R/stubproc BD=gnu_missing timeout 60 bash $R/t207_3_chain.sh flt def dev1 1e-9 > /dev/null 2>&1; echo "   new chain rc=$?"; cat $R/logs/chain_flt.txt | cut -c1-170 ); chk f_def_1em9_dev1
echo "== T7 time guard: DISKWAIT_S=-700 (now - tw0 > DISKWAIT_S + 600 at once), binary present"
mv $R/logs/chain_flt.txt $R/logs/chain_flt_T6.txt
( cd $R; HOME=$FH T207_3W=$R HOLDLOCK=1 DISKWAIT_S=-700 PROC_DIR=$R/stubproc BD=gnu_fake timeout 60 bash $R/t207_3_chain.sh flt def dev2 1e-9 > /dev/null 2>&1; echo "   rc=$?"; cat $R/logs/chain_flt.txt | cut -c1-170 ); echo "   run dir created: $([ -e $R/runs/f_def_1em9_dev2 ] && echo yes || echo no)"
echo "== T8 pr skips an already processed run (non-HOLDLOCK path): first call processes, second call skips"
rm -f $R/logs/chain_flt.txt $R/logs/stub_calls.txt
( cd $R; HOME=$FH T207_3W=$R PROC_DIR=$R/stubproc BD=gnu_fake bash $R/t207_3_chain.sh flt def dev1 2e-9 > /dev/null 2>&1; HOME=$FH T207_3W=$R PROC_DIR=$R/stubproc BD=gnu_fake bash $R/t207_3_chain.sh flt def dev1 2e-9 > /dev/null 2>&1 )
echo "   stub proc calls: $(wc -l < $R/logs/stub_calls.txt) (want 1)"; grep -E "already processed|refusing|rc != 0|no run" $R/logs/chain_flt.txt | cut -c1-150
