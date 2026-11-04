# Docking on Colab Pro+

POSEIDON has no GPU. These four scripts run `polyphar dock-worker` on a Colab GPU VM
through [`google-colab-cli`](https://pypi.org/project/google-colab-cli/), and are the
production form of the harness that validated Uni-Dock on an L4 on **2026-11-03**:

| Measured there | |
|---|---|
| Throughput | 0.157–0.180 s/ligand across boxes of 3,679–20,112 Å³ (flat over a 5.5× range) |
| Agreement with smina | n=580, Spearman rho=0.951, Uni-Dock +0.326 kcal/mol less negative |
| Determinism (issue #133) | 138–142 of 145 ligands move between identical same-seed runs, max 1.202 kcal/mol |

The determinism result is why a unit is `(target, replicate)` and why `.done` certifies
that a unit *was drawn*, never that redrawing it would agree.

## Why four short execs instead of one long one

The kernel websocket drops. It dropped during the validation run while the work was still
going. So nothing long-running is ever a child of the kernel: `launch.py` starts the worker
under `nohup` and returns immediately, and `poll.py` reads progress afterwards. Those polls
are also the activity that stops Colab reclaiming the VM.

## The round trip

On POSEIDON, in the run folder:

    export POLYPHAR_VERSION=1.1
    polyphar-docker dock-export --bundle 07_docking/export/<name> \
        --max-units N --engine unidock --engine-version 'Uni-Dock 1.2.0'
    tar czf bundle.tgz -C 07_docking/export <name>

Size `N` to what one session can finish. At ~0.18 s/ligand and 1,960 ligands a unit is
about six minutes, so ~200 units is a 20-hour session. **Do not export a second bundle
before merging the first** — both would contain the same units, both would be worked, and
the second merge's draws are dropped.

Then upload the bundle and, so the worker can reach you, the Pushover credentials:

    colab upload bundle.tgz /content/bundle.tgz -s <session>
    colab upload ~/.polyphar/pushover.env /content/pushover.env -s <session>

The credentials go as a **file, not through `colab exec --env`**: colab-cli records request
bodies in `~/.config/colab-cli/colab.log` at DEBUG level, and a token passed as an argument
can end up there. Without the file the worker runs silent — it still calls `notify()`, but
`notify()` finds no credentials, logs one warning and skips, and that log dies with the VM.
On a 20-hour session that is the difference between hearing it finish and discovering hours
later that it died.

Then:

    export POLYPHAR_COMMIT=$(git rev-parse HEAD)     # must be the commit that exported
    colab exec -f bin/colab/setup.py  --timeout 3600 --env POLYPHAR_COMMIT=$POLYPHAR_COMMIT
    colab exec -f bin/colab/launch.py --timeout 600
    colab exec -f bin/colab/poll.py   --timeout 300  # repeat
    colab exec -f bin/colab/pack.py   --timeout 600

Download `/content/results.tgz`, then back on POSEIDON:

    tar xzf results.tgz -C 07_docking/export/<name>
    polyphar-docker dock-merge --bundle 07_docking/export/<name>

`dock-merge` verifies each unit's key before folding it in, and recomputes what is still
outstanding. Only `results/` travels home: the receptors and the library came from the run
folder and merging re-reads them there.

## The commit pin is not optional

`POLYPHAR_COMMIT` must be the commit that produced the bundle. Unit keys record the code
version, so a worker running different code produces draws `dock-merge` cannot vouch for.
`setup.py` refuses to run without it.

## Watch the unit balance

Colab Pro+ had 1,598.9 units at the end of the validation session, against an estimated
1,446 for the full seven-subtype, three-replicate run. That has no margin. `setup.py`
prints the GPU it was given; an A100 burns units faster than an L4 for work that is not
compute-bound, and this work is not.
