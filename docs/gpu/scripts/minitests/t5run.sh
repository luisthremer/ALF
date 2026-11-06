cd ~/alf-gpu/work/dlsub
M=$1
base='s/^L1 .*/L1 = 4/
s/^L2 .*/L2 = 4/
s/^NSweep .*/NSweep = 1/
s/^NBin .*/NBin = 2/'
cfg() { # id | description | extra sed lines (newline separated, via printf) | optional ham
  id=$1; desc=$2; extra=$3
  { echo "$base"; printf '%b\n' "$extra"; } > sed_t5_$id.sed
  line="$id|$desc|"
  for s in standard submatrix; do
    out=$(SCHEME=$s TMO=120 bash rc.sh $M t5_${id}_$s sed_t5_$id.sed ${4:-Hubbard} | head -1)
    rc=$(echo "$out" | sed 's/.*rc=\([0-9]*\).*/\1/')
    act=$(echo "$out" | grep -o "ACTIVE\|INACTIVE" | head -1)
    line="$line $s:rc=$rc:${act:--}"
  done
  d=$(python3 cmp.py r_t5_${id}_standard/Ener_scal r_t5_${id}_submatrix/Ener_scal 2>&1 | tail -1)
  echo "$line | max rel diff std-vs-sub Ener_scal: $d"
}
cfg 01 "dense hopping" 's/^Checkerboard .*/Checkerboard = .F./'
cfg 02 "Honeycomb 4x4" 's/^Lattice_type .*/Lattice_type = "Honeycomb"/'
cfg 03 "Bilayer_square 4x4" 's/^Lattice_type .*/Lattice_type = "Bilayer_square"/'
cfg 04 "Bilayer_honeycomb 3x3" 's/^Lattice_type .*/Lattice_type = "Bilayer_honeycomb"/\ns/^L1 .*/L1 = 3/\ns/^L2 .*/L2 = 3/'
cfg 05 "N_leg_ladder 4x2" 's/^Lattice_type .*/Lattice_type = "N_leg_ladder"/\ns/^L2 .*/L2 = 2/'
cfg 06 "Symm=F" 's/^Symm .*/Symm = .F./'
cfg 07 "Ltau=0" 's/^Ltau .*/Ltau = 0/'
cfg 08 "Mz=F (N_FL=1)" 's/^Mz .*/Mz = .F./'
cfg 09 "Mz=F N_SUN=4" 's/^Mz .*/Mz = .F./\ns/^N_SUN .*/N_SUN = 4/'
cfg 10 "Continuous HS (type 3)" 's/^Continuous .*/Continuous = .T./'
cfg 11 "Nwrap=5" 's/^Nwrap .*/Nwrap = 5/'
cfg 12 "Nwrap=20" 's/^Nwrap .*/Nwrap = 20/'
cfg 13 "U=0" 's/^ham_U .*/ham_U = 0.d0/'
cfg 14 "chem pot 0.5" 's/^ham_chem .*/ham_chem = 0.5d0/'
cfg 15 "Projector" 's/^Projector .*/Projector = .T./\ns/^Theta .*/Theta = 2.d0/\ns/^Ltau .*/Ltau = 0/'
cfg 16 "Propose_S0" 's/^Propose_S0 .*/Propose_S0 = .T./'
cfg 17 "N_delay=1" 's/^N_delay .*/N_delay = 1/'
cfg 18 "N_delay=64" 's/^N_delay .*/N_delay = 64/'
cfg 19 "Dtau=0.05 Beta=2" 's/^Dtau .*/Dtau = 0.05d0/\ns/^Beta .*/Beta = 2.d0/'
cfg 20 "partial sequential range" 's/^Nt_sequential_start .*/Nt_sequential_start = 5/\ns/^Nt_sequential_end .*/Nt_sequential_end = 30/'
cfg 21 "cross-slice delay" 's/^Lowrank_cross_slice .*/Lowrank_cross_slice = .T./'
cfg 22 "Hubbard_Plain_Vanilla" '' Hubbard_Plain_Vanilla
s/^N_Global_tau .*/N_Global_tau = 1/'
cfg 24 "Langevin" 's/^Langevin .*/Langevin = .T./\ns/^Continuous .*/Continuous = .T./'
