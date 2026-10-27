# ALF GPU port – documentation index

| File | Purpose | Status |
|---|---|---|
| [HARD.md](HARD.md) | What makes this port hard, ranked by impact, evidence + strategy | Phase 0 ✅ |
| [PLAN.md](PLAN.md) | Architecture, phases, ticket list, risks, owner decisions | Phase 1 ✅ |
| [MACHINES.md](MACHINES.md) | Probed hardware/toolchains (`nividia`, `intel`, SuperMUC-NG Ph2) | Phase 0–1 |
| [BASELINE.md](BASELINE.md) | CPU builds, unit tests, test_vs_ed, Hubbard reference runs, validation lessons | Phase 0 |
| [PROFILE.md](PROFILE.md) | Measured CPU time per kernel (vs R1 §2), call counts, profiling method | Phase 0 |
| [ADR-001](ADR-001-programming-model.md), [ADR-002](ADR-002-alf-la-and-data-model.md) | Architecture decisions (programming model; `alf_la` and device data model) | accepted, Phase 1 ✅ |
| [BENCHMARKS.md](BENCHMARKS.md) | Time per sweep CPU vs GPU, per-kernel breakdown, crossover | per phase |
| [FEATURES.md](FEATURES.md) | ALF feature × {CPU, NVIDIA, Intel} support table, future work | maintained |
| [DISCREPANCIES.md](DISCREPANCIES.md) | Research/documentation claims contradicted by code or measurement | maintained |
| [tickets/](tickets/) | `NNN-short-name.md`, template `000-template.md` | maintained |
