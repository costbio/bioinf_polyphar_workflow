# Dock stage: GPU engine, replicates, and remote execution

**Date:** 2026-11-03
**Status:** design, not implemented
**Scope:** stage 5 of the panel redesign (`docs/superpowers/specs/2026-11-03-panel-selection-redesign-design.md`)

The current `dock` stage is built end to end around a **deterministic CPU engine running
on the same machine as the rest of the pipeline**. Three things break that at once:
the engine becomes Uni-Dock on GPU, scores become means over 5 replicates, and the GPU
lives on a different machine from the run folder. This is a rewrite of the stage's
execution model, not a parameter change.

## What is actually changing

| | now | after |
|---|---|---|
| engine | smina 2020.12.10 (Vina 1.1.2), CPU | Uni-Dock 1.2.0, CUDA |
| determinism | guaranteed; a test pins it | **impossible**; 95-98% of ligands move between identical runs |
| a score | one number | mean over 5 replicates, with SE |
| where it runs | in-container, same host as the run folder | a different machine with a GPU |
| work unit | target x 50-ligand chunk | **target x replicate**, whole library |

## Work unit: (target, replicate), not ligand chunks

Chunking by ligand was right for CPU smina, where each ligand is an independent
~200 core-second job. It is wrong for Uni-Dock, which batches ligands onto the GPU and is
**less** efficient on small inputs — the authors recommend at least 10,000 ligands per
invocation, and overhead dominates below a few hundred.

So one unit = **one target, one replicate, the whole prepared library**. For a subtype:

- 164 targets x 5 replicates = **820 units**
- each unit: 1,960 ligands x 0.176 s = **~345 s (~6 minutes)**
- per subtype: ~78.6 GPU-hours; all seven: ~550

Six minutes is a good granularity for a resumable remote job: short enough that losing
one to a dropped session is cheap, long enough that per-unit overhead is negligible.

## Execution model: session-sized bundles, pushed and returned

POSEIDON has no GPU. The GPU hosts are Colab Pro+ (primary), TRUBA, or a rented
container. Colab's constraints drive the design: VMs are **ephemeral**, there is **no
inbound network**, sessions are **attended** and capped at ~24 hours, and the client
hangs while the kernel keeps running.

That rules out a pull-based worker (no inbound connectivity, no persistent host) and a
shared filesystem. The DoGSiteScorer file-queue precedent does not transfer either: that
queue works because the container and the host worker see the same directory.

**The model is therefore: export a bundle, run it remotely, return results, merge.**

```
POSEIDON                                   GPU host (ephemeral)
--------                                   --------------------
polyphar dock-export  ──► bundle.tgz ────►  unpack
  receptors + pockets                       polyphar dock-worker
  prepared library                            for each outstanding unit:
  units.tsv (outstanding only)                  run Uni-Dock
  run_key.json                                  write r<i>.csv + .done
                       ◄──── results.tgz ◄──  pack only new results
polyphar dock-merge
  verify + merge into 07_docking/
  recompute what is outstanding
```

Each round trip is **session-sized**: `dock-export --max-units N` takes only as much work
as the session can plausibly finish, so a dropped session loses at most the units in
flight. The next export recomputes what is outstanding, so progress is monotone with no
coordination between sessions.

Sizes are small enough that this is not a constraint: ~100 MB out (receptors, pockets,
library for the selected targets), ~400 MB back per subtype of results.

**Bundle contents must be self-describing.** `run_key.json` carries the run identity, the
panel hash and the engine settings, so `dock-merge` can refuse results produced under
different settings rather than silently mixing them.

## Identity, and what `.done` can still mean

Today a chunk is done when `chunks/<chunk>.done` holds a job key covering the receptor,
pocket, chunk, seed, exhaustiveness, num_modes, addH and the **smina version**. That key
asserts something strong: re-running would produce the same numbers.

**With a non-deterministic engine that assertion is false**, and pretending otherwise is
the main hazard in this rewrite. The consequences:

- **A unit's identity includes its replicate index.** The key is
  `(receptor sha, pocket sha, library sha, replicate i, seed_i, search_mode, spacing,
  engine version)`.
- **Replicates are not interchangeable re-runs.** Replicate 2 is a distinct observation,
  not a repeat of replicate 1. Re-running a completed unit does not reproduce it, it
  produces a *different* sample.
- **Therefore resume must skip, never redo.** A `.done` marker certifies "this
  (target, replicate) was run and its results are recorded" — not "these are the scores".
  The result file is the only record of the scores.
- **A corrupted or partial unit is discarded and re-run as the same replicate index**,
  and that is the one case where a repeat is legitimate. It must be logged, because the
  discarded sample was a real draw and dropping it is a (small) selection effect.

Seeds: replicate `i` uses `base_seed + i` (42..46 by default), recorded per unit. The
engine's seeding is known broken upstream (issue #133), so this is provenance, not a
reproducibility claim — and the design must not imply otherwise.

## Results schema

Per target, keep both the raw draws and the summary:

```
07_docking/<label>/replicates/r0.csv .. r4.csv    raw: ID, Ligand, best_score
07_docking/<label>/chunks/r<i>.done               unit markers (job key)
07_docking/df_results_<label>.csv                 merged summary
```

`df_results_<label>.csv` columns: `Ligand`, `ID`, `Pocket`, `n`, `score_mean`,
`score_sd`, `score_se`, `score_min`, `score_max`.

`RESULT_COLUMNS = ["Pose", "SMINA_Score", "Ligand", "ID", "Pocket"]` is retired. The
engine-specific name goes, and a single score column cannot express a distribution.

**Keep the raw per-replicate scores.** They cost little, they let the analysis recompute
any summary, and they are the only way to revisit the jitter question later — which is
still open, since the root cause (floating-point non-associativity against VRAM-driven
batch composition) was never settled.

### Poses

Five replicates x 164 targets x 1,960 ligands of poses is far more than is useful. **Keep
poses for replicate 0 only**, and make even that optional (`--poses none|first|all`).

Replicate 0's poses are what the **pose-in-pocket validation** consumes — the outstanding
test of whether a best-scoring pose actually sits in the selected cavity, which is the
assumption the whole pocket-similarity premise rests on. So `first` is the default.

## Provenance the manifest must record

Per unit, not just per stage: engine version, `unidock-tools` version, CUDA runtime,
**GPU model** (throughput and plausibly scores vary by device), `search_mode`,
`--batch_size`, replicate index and seed, prepared-library sha, receptor and pocket sha,
and the host path used (Colab / TRUBA / container image digest or lockfile hash).

GPU model matters enough to be a first-class field: if replicates for one target are
split across an L4 and an A100, that is a potential systematic difference inside what is
being averaged. **The merge step should warn when a target's replicates span GPU models.**

## Host checks and preflight

- `run`'s preflight checks `smina --version` and will exit 2 without it. It must become
  engine-aware, and must not demand smina when the engine is Uni-Dock.
- `polyphar dock` on a host with no GPU must **fail fast with guidance**, not start and
  crash. POSEIDON is exactly that host, so this is the common case, not an edge case.
- `polyphar run` must stop cleanly before `dock` on a GPU-less host, having completed
  everything up to it, rather than treating it as a failure.
- The worker must verify the bundle's `run_key.json` against its own settings before
  doing any work.

## ETA

The existing calibration is a single seconds-per-ligand figure from smina and does not
apply. Uni-Dock needs its own: measured 0.176 s/ligand on an L4, and **flat across a 5x
range of box volume**, which is a real simplification over CPU smina where box volume was
the dominant cost driver. So a per-(target, replicate) estimate of
`n_ligands x rate(GPU model)` is adequate, with the rate calibrated per device class.

## Tests

**Retire:**
- `test_cpu_and_chunk_invariance` — pins a property the engine no longer has. Deleting it
  is correct, not a concession; leaving it would force the rewrite to fake determinism.
- `test_reproduces_2023_targets` — the 2023 comparison is a non-goal.
- Any test asserting a single `SMINA_Score` column.

**Add:**
- replicate index is part of the unit key; two replicates of one target are distinct units
- resume skips a completed replicate and never re-runs it
- a corrupted unit is discarded, re-run under the same replicate index, and logged
- the merged summary's `n`, mean and SE are correct over known inputs
- `dock-merge` refuses results whose `run_key.json` disagrees with the run
- a merge spanning GPU models warns
- a GPU-less host fails fast with a usable message
- poses are written for replicate 0 only under the default

Determinism tests become **statistical**: with a stubbed engine returning fixed values,
the summary must be exact; with one returning jittered values, the SE must track the
input distribution. Do not assert exact scores against the real engine anywhere.

## Open questions

1. **How are replicates distributed across sessions?** All 5 replicates of one target
   together (same GPU, cheaper to merge, but a lost session loses a whole target), or
   replicate-major across targets (spreads device variation, more round trips)? The GPU
   model warning above argues for target-major; resilience argues the other way.
2. **Is 5 replicates enough?** Chosen by decision, not measurement. The SE-against-N curve
   is an outstanding test, and the threshold-flipping problem — `min_targets` counts
   targets below a cutoff, so jitter near the cutoff flips classifications in a way
   averaging the score does not fix — has no power statement yet.
3. **Does the score offset survive conformer preparation?** The +0.326 kcal/mol
   smina-to-Uni-Dock offset was measured on *unprepared* ligands. Stage 4 changes them.
4. **TRUBA has never been sized** — no partition, no queue wait, no validated Apptainer
   build. It is the named fallback for when attended Colab sessions become impractical,
   and at ~26 attended session-days for the full production run, that is likely.

## Implementation order

1. Split the stage: `dock-export`, `dock-worker`, `dock-merge`, keeping `dock` as the
   single-host convenience path for a machine that has both a GPU and the run folder.
2. Replicate-aware identity, keys and resume.
3. The Uni-Dock engine adapter and its ETA.
4. Results schema and the merge summary.
5. Host checks and preflight.
6. Bundle format and `run_key.json` verification.

Steps 1 and 2 carry the design risk; the rest is mechanical.
