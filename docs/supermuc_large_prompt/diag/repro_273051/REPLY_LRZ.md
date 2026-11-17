# Reply to LRZ incident 273051 (draft for the owner; send it yourself)

**Subject:** Re: [273051] Faulty GPU tile i20r02c06s02: one-node reproducer

*Transparency note: this message was co-authored with Claude Opus 5.5 (Anthropic), an AI assistant. It assembled the reproducer and drafted the text. All job IDs and numbers come from our own job logs.*

Dear Mr. Block,

Thank you for the swap. Yes, here is a one-node reproducer. It is the same input that failed on the original node.

**Location:** `/hppfs/scratch/0C/di38zuw/large_prompt/diag/repro_273051/`, also packed as `/hppfs/scratch/0C/di38zuw/large_prompt/diag/repro_273051.tar.gz`. Both are only readable with admin rights, since my scratch directory is private. The tarball (12 MB) is also attached to this mail.

**How to run it:** `sbatch -w <node> run_repro.sh [NREP]`, from that directory or a copy of it. It needs no other files.
- Account pn93la is set in the script; please change it to yours.
- One repetition takes about 5 minutes on a healthy node. The default is 3 repetitions.

**What it does:**
- Our application ALF runs with 8 MPI ranks, one per tile (flat hierarchy, rank r uses device r).
- It runs a fixed input with fixed random seeds.
- It uses OpenMP offload with Level Zero and oneMKL, built with intel/2025.3.0, intel-mkl/2025.3.0 and intel-mpi/2021.15.0.
- The binary and its source are included.

**How to read the result:** the script prints one line per repetition.
- **Healthy:** `rc=0 ... acceptance=0.318327091942149 verdict=OK`. This result is bit-identical on every healthy node we used, including a check today on i20r06c06s10 (job 562002, 2 repetitions, also identical in Precision Green to the reference).
- **Faulty:** `rc=3 ... verdict=FAULT | Abort(3) on node 1 | ... exceeding the threshold of 10 for G difference!`. "node 1" means MPI rank 1, i.e. device 1. On the original node this happened in **4 of 4 runs** (jobs 560331 and 560620), each time about 3 minutes into the run. With tile 1 hidden (`ZE_AFFINITY_MASK=0,2,3,4,5,6,7`, job 560623) the same node ran clean.
- The `reference/` folder holds the healthy result and the log of the failed run.

**A lighter test, with a caveat:** `stress/` contains a small standalone test. It repeats zgemm, zgetrf, zgeqrf, zunmqr and ztrsm at N=1024 and checks every result bit for bit. It catches an injected error on all 8 tiles. However, it was never run on the faulty tile, so we cannot say whether it reproduces the fault. Please use the ALF reproducer as the primary test.

Kind regards,
Luis Thremer
