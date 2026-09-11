# Replication curation record

Date: 2026-09-11. Scope: the paper's final two Quadron minimal models only.
Budget: one provenance audit and one independent code review, with focused executable checks. No open-ended analysis search.

## Queue

- Complete: identify the four upstream/final analysis scripts and render/validation helpers.
- Complete: comment the exact formulas, statistical assumptions, sequence/coordinate handling, and sample joins.
- Complete: test a slope-only entry point against direct `lm` and exclude older first-stage regression side products.
- Complete: add an explicit configuration, numerical parity gate, and protected-output runner.
- Complete: execute both models from saved checkpoints and compare every result column.
- Complete: validate shared-manifest mapped coordinates and reconstruct all three slope checkpoints from prepared inputs, then reproduce both models again.
- Complete: comment rendering helpers, document the replay routes, review code, and report remaining gaps.
- Complete: preserve portable verification summaries; all curated R scripts parse and the focused input/comparison/path tests pass.

## Provenance audit

Read-only independent audit found no historical command for either final prefix. The recovered Nextflow call was a help/configuration invocation, not evidence of an executed model. Existing checkpoints include all three prepared beta/metadata pairs, three annotated slope files, and merged/unmerged Quadron windows. Original Quadron class-separated BEDs and a verified original class cutoff/version were not recovered locally. The default sequence flank is 1000 bp, but no historical command proved that setting. Upstream manifest identity and the coordinate convention remain unresolved.

## Candidate replay routes

Generation/review provenance: supervisor and independent provenance audit, this curation pass.

1. Saved slope + unmerged-window checkpoint replay: preferred, because it tests the published two models without reconstructing uncertain upstream choices.
2. Prepared beta + merged/unmerged windows: secondary, because manifest identity and sample-dependent first-stage computations must also match.
3. GEO processed matrices + original Quadron BEDs: optional source-level route, gated on explicitly supplied original inputs; not equivalent to reproducing raw-IDAT normalization or the Quadron training/prediction process.

Pairwise decisions: route 1 over route 2 (fewer unresolved upstream inputs); route 2 over route 3 (local checkpoints available, original BED provenance absent). These rank replay feasibility, not scientific hypotheses or evidence strength. None is a replacement for independent biological replication.

Editorial Elo ledger (initial 1000, K=32): after route 1 defeats route 2, scores are 1016/984/1000; after route 2 defeats route 3, scores are approximately 1016/1000.7/983.3. Following the successful prepared-input reconstruction, route 2 is retained as a verified, more extensive companion to the faster route 1. These scores express workflow preference only and are not evidence probabilities.

The selected local hypothesis is that the existing script's minimal branch, explicit Quadron checkpoints and sequence flank of 1000 can reproduce all six paper tables. A table-by-table numerical comparison can disconfirm it. A mismatch must be reported; expected tables and model definitions must never be changed to force agreement.

## Executed evidence

The checkpoint-based `models` run and the prepared-beta `prepared` run both completed. Both matched all columns of all six manuscript model tables. In the checkpoint replay, the largest parsed numerical difference was zero in every table. The prepared run verified its three full regenerated slope/annotation tables before final modeling. The shared manifest contained 657,261 unique probe IDs and 381,179 mapped standard-chromosome entries; every saved checkpoint coordinate matched the mapped subset. Unmapped entries are allowed and left to the historical exclusion rule rather than guessed or repaired.

## Code review and repair

The single independent code review identified (1) unresolved `..` traversal through a nonexistent ancestor and (2) a verification report symlink that could target protected data. Both were repaired locally: physical ancestor resolution rejects ambiguous traversal, full report destinations and model/reference trees are checked, and reports are written by temporary-file rename. The new `test_paths.R` fixtures reproduce the unsafe path patterns and verify unchanged sentinel hashes. Re-verification of the numerical results passed after these changes.

Upstream parser and slope-only additions were tested with offline synthetic fixtures; comments otherwise preserve existing scientific calculations. The earlier notebook, obsolete analysis variants, and source data were not edited. Tests of the prepared route add real evidence for that route but do not recover original Quadron threshold or raw-array preprocessing provenance.

## Completion and stopping reason

[collect_validation.R](collect_validation.R) packaged the two completed replay summaries and three upstream parity reports under [verification/](verification/). The checkpoint route's maximum parsed numerical error was zero; the prepared route's largest final-table difference was approximately 9.95e-14, below the declared tolerances. The current successful package versions are recorded separately from historical-environment claims.

Completed gates: script parsing; local processed-GEO import; 100-bp mixed-window annotation; first-stage slope/SE/mean parity with missing-beta handling; strict comparison failure cases; filesystem protection cases; real three-tissue checkpoint replay; real prepared-input reconstruction with full slope checkpoint comparison; six-table final model parity in both routes; unchanged-paper artifact validation. No original source/result/figure was overwritten. One inline evidence-packaging command stalled because of terminal input corruption; it was stopped and replaced by the successful saved packaging script. The scientific reruns had already completed and were unaffected.

The one provenance audit and one independent code review are complete. Stop on the successful executable gates and stable replay scope. The optional full processed-source route remains unexecuted on original class BEDs; original Quadron labels, raw-IDAT preprocessing, and independent coordinate/strand correctness are disclosed follow-up provenance requirements, not inferred from successful numerical replay.