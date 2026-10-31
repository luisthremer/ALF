## Default N_delay

**Chosen default: N_delay = 32** (a single value; the optimum does not move with N). At L=12 (62 accepts per slice) the
best points are N_delay=16 and 32 for both schemes, within 1-2 % of each other, which is inside the run-to-run noise (about
2 % between the three standard L=12 runs); at L=20 (172 accepts per slice) 32 is the best for both schemes (delayed 1.172,
submatrix 1.112 s/sweep, N_delay=64 within 1 % for submatrix). Beyond 64 the time rises again (L=20: N_delay=256 is 33 %
slower than 32 for delayed and 44 % slower for submatrix), so a large value is not free in time even though the per-slice
cap limits the memory: the pending-vertex work per proposal grows with the number of pending vertices, so the plateau
expected from the cap is not seen; the cap only prevents N_delay > vertices per slice from mattering. Below 16 the flush
count is too high (N_delay=4 is 1.3-1.5x slower than 32 at L=20).

Precision against N_delay: `Precision Green Max` is 5.9e-7 (L=12) and 5.6e-7 (L=20) for every run of both schemes and for standard
(ratio 1.00 to three digits); it is set by one event of the trajectory (identical for all runs) and carries no information
about N_delay. The mean is more sensitive: the submatrix Green mean grows slightly with N_delay (L=20: 8.75e-12 at 4 to
9.12e-12 at 256, +4 %; L=12: 1.90e-11 to 1.95e-11, +2.6 %), the delayed mean less (8.75e-12 to 8.87e-12, +1.4 %); at the chosen
N_delay=32 both are within 0.1-0.4 % of standard (8.82e-12 vs 8.82e-12 at L=20 for submatrix, 8.75e-12 for delayed).
`Precision Phase Max` scatters between 0.9e-14 and 1.6e-14 without trend (standard: 0.9e-14 at L=12, 1.4e-14 at L=20).
The ratio of 1.9 between submatrix and standard Green precision max seen on the density set at Checkpoint 4 does not
appear for this Mz Hubbard model at any N_delay in {4,...,256}.
