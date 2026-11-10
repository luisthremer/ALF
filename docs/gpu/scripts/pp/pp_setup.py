# usage: pp_setup.py DIR NMPI_PER_SET NSWEEP NBIN [NDELAY]   (run with the venv python on nividia)
import sys, os, re, shutil, subprocess, importlib.util
d, nper, nsweep, nbin = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
nd = int(sys.argv[5]) if len(sys.argv) > 5 else 32
s = importlib.util.spec_from_file_location('p', os.path.expanduser('~/alf-gpu/work/rchi_pipeline.py')); m = importlib.util.module_from_spec(s)
try: s.loader.exec_module(m)
except SystemExit: pass
t = m.params_text(24, 5.0, nsweep, nbin)
t = re.sub(r"ham_U = [^\n]*", "ham_U = {3.0d0, 5.25d0, 8.0d0}", t)
t = t.replace('ham_name = "Hubbard"', 'ham_name = {"Hubbard"}')
t = re.sub(r"N_delay = \d+", "N_delay = %d" % nd, t)
t = t.replace("&VAR_errors", "&VAR_TEMP\nN_exchange_steps = 0\nN_Tempering_frequency = 100000\nmpi_per_parameter_set = %d\nTempering_calc_det = .F.\n/\n&VAR_errors" % nper)
if len(sys.argv) > 6: t = re.sub(r"CPU_MAX = [^\n]*", "CPU_MAX = " + sys.argv[6], t)
os.makedirs(d, exist_ok=True)
open(os.path.join(d, 'parameters'), 'w').write(t)
shutil.copy(os.path.expanduser('~/alf-gpu/work/rchi/seeds_template'), os.path.join(d, 'seeds'))
shutil.copy(os.path.expanduser('~/alf-gpu/work/pp/prepare_sim.py'), d)
subprocess.run([sys.executable, 'prepare_sim.py'], cwd=d, check=True)
print(sorted(os.listdir(d)))
