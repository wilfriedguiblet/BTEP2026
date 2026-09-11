# Validation of the curated Quadron replay

Date: 2026-09-11. All tests used existing local data or small synthetic fixtures. No original biological input, source result, paper table, or manuscript figure was replaced.

## Actual numerical reruns

1. **Saved annotated slopes -> M0/M1:** completed under `results/quadron_paper_reproduction` in BTEP.
2. **Prepared beta/metadata + manifest + merged windows -> three slope tables -> M0/M1:** completed under `results/quadron_paper_reproduction_prepared` in BTEP.

Both routes used the original four analysis scripts with explicit minimal-model flags and a 1,000-bp sequence flank. Both checked all columns of six final model tables against the immutable, hash-checked paper references. The prepared route additionally checked every column in each regenerated slope/annotation table. [environment_validated.csv](environment_validated.csv) records the current successful environment; per-run code/config copies and full commands are retained in each run's provenance directory.

Portable diagnostic copies are retained in [verification/](verification/). They are verification artifacts, not replacement paper results. Numerical equivalence is the recorded success criterion; regenerated image bytes are not required to match the original figures.

## Functional and protection tests

- Offline processed-GEO fixture imports matching sample/metadata order and beta values.
- Synthetic 100-bp window expansion merges overlapping classes into a shared neighborhood with both indicators retained.
- First-stage slopes, SEs, and means match direct OLS, including a probe with a missing beta value.
- `--slopes-only` emits only the slope checkpoint and no excluded regression families.
- The comparator accepts row reordering and harmless floating-point differences, but rejects changed counts, formulas, tiny P values, missing estimates, and missing/duplicate rows.
- Path tests reject nonexistent-parent traversal, report-file symlinks, report-parent symlinks that escape outputs, model-table report destinations, and dangling symlinks. Reference sentinel hashes remain unchanged.
- The pre-existing manuscript artifact validator passes after curation; the two original PNGs, six paper tables, citations and links remain valid.

## What remains unverified

The source-level mode was fixture-tested but not run against the original full class-separated Quadron BEDs, which were not available locally. Original Quadron scoring/classification provenance, raw-IDAT normalization/probe-QC history, and independent correctness of the historical coordinate/strand conventions remain unresolved. Successful replay from retained checkpoints demonstrates numerical reproducibility of the presented results, not those upstream provenance facts, a dependence-robust inference model, or independent biological replication.

One independent provenance audit and one independent code review were completed. The review identified two path-protection gaps; both were repaired and covered by regression tests. No additional analysis variants or undocumented parameter searches were introduced to obtain agreement.