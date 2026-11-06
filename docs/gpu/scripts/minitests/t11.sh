cd ~/alf-gpu/work/dlsub
for L in 20 32; do
  for nd in 8 16 32 64 128; do
    printf '%s\n' "s/^L1 .*/L1 = $L/" "s/^L2 .*/L2 = $L/" 's/^NSweep .*/NSweep = 1/' 's/^NBin .*/NBin = 1/' "s/^N_delay .*/N_delay = $nd/" > sed_t11.sed
    out=$(VAR=srcT SCHEME=submatrix TMO=300 bash rc.sh nv t11_${L}_$nd sed_t11.sed | head -1)
    wall=$(echo "$out" | sed 's/.*wall=\([0-9.]*\).*/\1/')
    up=$(grep -h "GPU time in sequential" r_t11_${L}_$nd/info | sed 's/.*: *\([0-9.]*\) s.*/\1/')
    nb=$(grep -h "low-rank update counts" r_t11_${L}_$nd/info | sed 's/.*blocks \([0-9]*\).*/\1/')
    echo "L=$L N_delay=$nd wall=${wall}s update_time=${up}s blocks=$nb"
  done
done
