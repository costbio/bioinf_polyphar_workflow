#!/usr/bin/env python3
"""Calibrate a DoGSite3 volume cutoff on DoGSite3's own scale against the 2023 DoGSiteScorer 2
selection (drugScore >= 0.5), and compare the two selections pocket by pocket.

    # 1. run a tool on the 2023 inputs (pdbfixer chains + AlphaFold models), as the pipeline does
    python scripts/calibrate_dogsite3_volume.py run --inputs <dir> --out <dir> \
        --dogsite /opt/ZBH/dogsite3/dogsite3 --tool dogsite3 --nprocs 60 [--limit N --seed S]
    # 2. analyse
    python scripts/calibrate_dogsite3_volume.py analyze --ds2-drugscore DrugScore.csv \
        --ds2-pockets <2023 pockets dir> --ds3 <dogsite3 out dir> --json result.json
    # 3. optional: compare a rerun of DoGSiteScorer 2 with the 2023 DrugScore.csv
    python scripts/calibrate_dogsite3_volume.py compare-ds2 --ds2-drugscore DrugScore.csv \
        --rerun <dogsite2 out dir>

Inputs layout (`run --inputs`): <inputs>/pdbfixer/<UniProt>/<pdb>_<chain>.pdb and
<inputs>/af/<UniProt>/AF-<acc>-F1.pdb (the 2023 HCC net2 files). Outputs go to
<out>/<UniProt>/<stem> with the pipeline's command lines (bin/polyphar-dogsite-worker's
dogsite_cmd: dogsite3 `-p <model> -o <out> -s -d --writeSiteResiduesPDB -v 1`; dogsite2
`-s -i -y -d -w 4`). A model whose _desc.txt exists is skipped (resumable).

Models: a 2023 DrugScore.csv row belongs to (UniProt ID, Protein) with Protein = 2OGQ_A for a
PDB chain and AF-<acc> for an AlphaFold model; a DoGSite3 output to (UniProt, upper-cased
<pdb>_<chain>) or (UniProt, AF-<acc>). Only models with results of BOTH tools count.

Pockets AND subpockets are counted alike on both sides (every _desc.txt row: P_<n> and
P_<n>_<m>; both tools were run with -s), because the pipeline selects among all rows: 2023
kept any row with drugScore >= 0.5, the DoGSite3 stage keeps any row with volume >= cutoff.
A main-pockets-only variant (P_<n> rows) is reported as a sensitivity check.

Cutoff: f = (DS2 rows with drugScore >= 0.5) / (all DS2 rows) on the common models; the DS3
cutoff is the volume v such that the same fraction f of all DS3 rows on those models has
volume >= v, i.e. the (1 - f) quantile (numpy "higher" so that >= keeps at least f), rounded
to the nearest 10 A^3 ("proposed").

Agreement: residue sets (chain, residue number + insertion code) of each pocket's residue file
(DS2 <...>_res_P_*.pdb of 2023, DS3 <...>_P_*_res.pdb; both tools ran on the same files, so
numbering agrees). Pockets match when the Jaccard index of their residue sets is >= 0.5.
Per model, precision = DS3-selected rows matching some DS2-selected row / DS3-selected rows,
recall = DS2-selected rows matched by some DS3-selected row / DS2-selected rows (pooled over
models: micro averages); "matchable" = DS2-selected rows matched by ANY DS3 row (upper bound
of recall at any cutoff). DS2 rows whose residue file is missing are left out of matching.
"""
import argparse
import concurrent.futures
import importlib.machinery
import importlib.util
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
JACCARD = 0.5
DS2_SELECT = 0.5
SANITY = [("P53350", "2OGQ_A"), ("P51955", "2XKC_A"), ("O60218", "4GA8_A"),
          ("P06493", "6GU2_A")]


def _worker():
    path = REPO / "bin" / "polyphar-dogsite-worker"
    loader = importlib.machinery.SourceFileLoader("polyphar_dogsite_worker", str(path))
    mod = importlib.util.module_from_spec(importlib.util.spec_from_loader(loader.name, loader))
    loader.exec_module(mod)
    return mod


# --- run ---------------------------------------------------------------------------------

def _models(inputs: Path):
    out = []
    for sub in ("pdbfixer", "af"):
        for p in sorted((inputs / sub).glob("*/*.pdb")):
            with open(p, "rb") as fh:
                head = fh.read(6)
            if head.startswith((b"HEADER", b"ATOM", b"REMARK", b"CRYST1", b"MODEL")):
                out.append(p)
    return out


def cmd_run(args):
    w = _worker()
    models = _models(Path(args.inputs))
    if args.limit:
        random.Random(args.seed).shuffle(models)
        models = sorted(models[:args.limit])
    outdir = Path(args.out).resolve()   # absolute: the tool runs with cwd = its out dir
    log = outdir / "times.tsv"
    outdir.mkdir(parents=True, exist_ok=True)

    def one(model: Path):
        acc, stem = model.parent.name, model.stem
        d = outdir / acc
        d.mkdir(exist_ok=True)
        out = w.out_for(str(d / stem), args.tool)
        if os.path.isfile(out + "_desc.txt"):
            return None
        t0 = time.monotonic()
        try:
            r = subprocess.run(w.dogsite_cmd(args.dogsite, str(model.resolve()), out, args.tool),
                               cwd=d, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True, errors="replace", timeout=args.timeout)
            rc, err = r.returncode, r.stderr or ""
        except subprocess.TimeoutExpired:
            rc, err = "timeout", ""
        return f"{acc}\t{stem}\t{rc}\t{time.monotonic() - t0:.2f}\t" \
               f"{' '.join(err.split())[-200:]}"

    done = 0
    with concurrent.futures.ThreadPoolExecutor(args.nprocs) as ex, open(log, "a") as fh:
        for line in ex.map(one, models):
            done += 1
            if line:
                fh.write(line + "\n")
                fh.flush()
    print(f"{done} models, log {log}")


# --- analysis ------------------------------------------------------------------------------

def residues(path: Path) -> frozenset:
    out = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 27:
                out.add((line[21], line[22:27].strip()))
    return frozenset(out)


def jaccard(a, b) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0


def _ds2_table(drugscore, pockets_dir: Path):
    import pandas as pd

    df = pd.read_csv(drugscore, index_col=0, dtype={"Pocket": str, "Protein": str})
    df["model"] = list(zip(df["UniProt ID"], df["Protein"].str.upper()))
    listing = {}
    files = []
    for acc, fname in zip(df["UniProt ID"], df["Pocket Filename"]):
        if acc not in listing:
            d = pockets_dir / acc
            listing[acc] = {n.lower(): d / n for n in os.listdir(d)} if d.is_dir() else {}
        files.append(listing[acc].get(str(fname).lower()))
    df["file"] = files
    return df


def _ds3_table(ds3: Path):
    import pandas as pd

    frames, seen = [], {}
    for desc in sorted(ds3.glob("*/*_desc.txt")):
        acc, stem = desc.parent.name, desc.name[:-len("_desc.txt")]
        protein = stem[:-3] if stem.startswith("AF-") and stem.endswith("-F1") else stem.upper()
        key = (acc, protein)
        seen.setdefault(key, []).append(stem)
        d = pd.read_csv(desc, sep="\t", dtype={"name": str})
        d = d.rename(columns={"name": "Pocket"})
        d["UniProt ID"], d["Protein"], d["stem"] = acc, protein, stem
        d["file"] = [desc.parent / f"{stem}_{p}_res.pdb" for p in d["Pocket"]]
        frames.append(d[["UniProt ID", "Protein", "stem", "Pocket", "volume", "file"]])
    df = pd.concat(frames, ignore_index=True)
    ambiguous = {k for k, v in seen.items() if len(v) > 1}   # chains differing only by case
    df["model"] = list(zip(df["UniProt ID"], df["Protein"]))
    return df[~df["model"].isin(ambiguous)].reset_index(drop=True), sorted(ambiguous)


def _quantile_cutoff(volumes, f):
    import numpy as np

    v = np.quantile(np.asarray(volumes, float), 1 - f, method="higher")
    return float(v), float(round(v / 10.0) * 10)


def _match(ds2, ds3, cutoff):
    """Pooled precision / recall / matchable of DS3 (volume >= cutoff) vs DS2 (drugScore >=
    0.5) selections; residue sets are cached on the frames (column `res`)."""
    tp3 = n3 = rec2 = n2 = matchable = 0
    per_model = []
    g3 = dict(tuple(ds3.groupby("model")))
    for model, a in ds2.groupby("model"):
        b = g3.get(model)
        if b is None:
            continue
        s2 = [r for r, s in zip(a["res"], a["drugScore"]) if s >= DS2_SELECT and r is not None]
        all3 = list(b["res"])
        s3 = [r for r, v in zip(b["res"], b["volume"]) if v >= cutoff]
        hit3 = sum(any(jaccard(x, y) >= JACCARD for y in s2) for x in s3)
        hit2 = sum(any(jaccard(y, x) >= JACCARD for x in s3) for y in s2)
        can2 = sum(any(jaccard(y, x) >= JACCARD for x in all3) for y in s2)
        tp3, n3, rec2, n2, matchable = tp3 + hit3, n3 + len(s3), rec2 + hit2, n2 + len(s2), \
            matchable + can2
        per_model.append({"model": f"{model[0]} {model[1]}", "ds2_selected": len(s2),
                          "ds3_selected": len(s3), "ds3_hits": hit3, "ds2_recalled": hit2})
    return {"cutoff": cutoff, "ds3_selected": n3, "ds2_selected_with_residues": n2,
            "precision": tp3 / n3 if n3 else None, "recall": rec2 / n2 if n2 else None,
            "matchable": matchable / n2 if n2 else None,
            "models_both_selected_nonzero": sum(1 for m in per_model
                                                if m["ds2_selected"] and m["ds3_selected"])}


def _sanity(ds2, ds3, cutoff):
    out = []
    for acc, protein in SANITY:
        a = ds2[ds2["model"] == (acc, protein)]
        b = ds3[ds3["model"] == (acc, protein)]
        if a.empty or b.empty:
            out.append({"model": f"{acc} {protein}", "note": "not in both"})
            continue
        top = a.sort_values("drugScore", ascending=False).iloc[0]
        best = None
        for _, r in b.iterrows():
            j = jaccard(top["res"], r["res"]) if top["res"] is not None else 0.0
            if best is None or j > best[0]:
                best = (j, r)
        j, r = best
        out.append({"model": f"{acc} {protein}",
                    "ds2_top": {"pocket": top["Pocket"], "drugScore": round(float(top["drugScore"]), 3),
                                "volume": round(float(top["volume"]), 1)},
                    "ds2_max_volume": round(float(a["volume"].max()), 1),
                    "ds3_best_match": {"pocket": r["Pocket"], "volume": round(float(r["volume"]), 1),
                                       "jaccard": round(j, 2),
                                       "selected_at_cutoff": bool(r["volume"] >= cutoff)},
                    "ds3_max_volume": round(float(b["volume"].max()), 1),
                    "ds2_selected": int((a["drugScore"] >= DS2_SELECT).sum()),
                    "ds3_selected_at_cutoff": int((b["volume"] >= cutoff).sum())})
    return out


def cmd_analyze(args):
    import numpy as np

    ds2 = _ds2_table(args.ds2_drugscore, Path(args.ds2_pockets))
    ds3, ambiguous = _ds3_table(Path(args.ds3))
    common = set(ds2["model"]) & set(ds3["model"])
    ds2 = ds2[ds2["model"].isin(common)].reset_index(drop=True)
    ds3 = ds3[ds3["model"].isin(common)].reset_index(drop=True)
    ds2["res"] = [residues(p) if p is not None and p.is_file() else None for p in ds2["file"]]
    ds3["res"] = [residues(p) if p.is_file() else None for p in ds3["file"]]

    res = {"models_ds2": int(len(set(_ds2_table(args.ds2_drugscore,
                                                 Path(args.ds2_pockets))["model"]))),
           "models_ds3": int(len(set(_ds3_table(Path(args.ds3))[0]["model"]))),
           "models_common": len(common), "ambiguous_ds3_models": [list(a) for a in ambiguous],
           "jaccard_match": JACCARD}
    variants = {}
    for name, keep2, keep3 in (
            ("all rows (pockets + subpockets)", np.ones(len(ds2), bool), np.ones(len(ds3), bool)),
            ("main pockets only (P_<n>)", ds2["Pocket"].str.count("_") == 1,
             ds3["Pocket"].str.count("_") == 1)):
        a, b = ds2[keep2], ds3[keep3]
        f = float((a["drugScore"] >= DS2_SELECT).mean())
        raw, proposed = _quantile_cutoff(b["volume"], f)
        v = {"ds2_rows": int(len(a)), "ds3_rows": int(len(b)),
             "ds2_selected": int((a["drugScore"] >= DS2_SELECT).sum()), "f": f,
             "cutoff_quantile": raw, "cutoff_proposed": proposed,
             "ds3_selected_at_proposed": int((b["volume"] >= proposed).sum()),
             "models_with_selected_ds2": int(a[a["drugScore"] >= DS2_SELECT]["model"].nunique()),
             "models_with_selected_ds3": int(b[b["volume"] >= proposed]["model"].nunique()),
             "models_with_selected_ds3_at_300": int(b[b["volume"] >= 300]["model"].nunique()),
             "ds3_selected_at_300": int((b["volume"] >= 300).sum())}
        v["agreement"] = [_match(a, b, c) for c in sorted({proposed, raw, 300.0,
                                                          *args.extra_cutoffs})]
        variants[name] = v
    res["variants"] = variants
    m2 = ds2.groupby("model")["volume"].max()
    m3 = ds3.groupby("model")["volume"].max()
    ratio = (m3 / m2.reindex(m3.index)).dropna()
    q = ratio.quantile([0.25, 0.5, 0.75])
    res["largest_volume_ratio_ds3_over_ds2"] = {"n": int(len(ratio)), "q25": float(q[0.25]),
                                                "median": float(q[0.5]), "q75": float(q[0.75])}
    res["ds2_rows_without_residue_file"] = int(ds2["res"].isna().sum())
    res["ds3_rows_without_residue_file"] = int(ds3["res"].isna().sum())
    proposed = variants["all rows (pockets + subpockets)"]["cutoff_proposed"]
    res["sanity"] = _sanity(ds2, ds3, proposed)
    text = json.dumps(res, indent=2, default=str)
    if args.json:
        Path(args.json).write_text(text + "\n")
    print(text)


# --- DoGSiteScorer 2 rerun vs 2023 ---------------------------------------------------------

def cmd_compare_ds2(args):
    import pandas as pd

    old = pd.read_csv(args.ds2_drugscore, index_col=0, dtype={"Pocket": str, "Protein": str})
    old["model"] = list(zip(old["UniProt ID"], old["Protein"].str.upper()))
    rows = []
    for desc in sorted(Path(args.rerun).glob("*/*.pdb_desc.txt")):
        acc, stem = desc.parent.name, desc.name[:-len(".pdb_desc.txt")]
        protein = stem[:-3] if stem.startswith("AF-") and stem.endswith("-F1") else stem.upper()
        new = pd.read_csv(desc, sep="\t", dtype={"name": str}).set_index("name")
        o = old[old["model"] == (acc, protein)].set_index("Pocket")
        same = sorted(new.index) == sorted(o.index)
        diff = None
        if same and len(o):
            diff = float((new["drugScore"] - o.loc[new.index, "drugScore"]).abs().max())
        rows.append({"model": f"{acc} {protein}", "in_2023": len(o) > 0, "pockets_new": len(new),
                     "pockets_2023": len(o), "same_names": same, "max_abs_drugscore_diff": diff})
    df = pd.DataFrame(rows)
    in23 = df[df["in_2023"]]
    summary = {"models_rerun": len(df), "models_in_2023": int(len(in23)),
               "same_pocket_names": int(in23["same_names"].sum()),
               "drugscore_within_1e-4": int((in23["max_abs_drugscore_diff"] <= 1e-4).sum()),
               "max_abs_drugscore_diff": (float(in23["max_abs_drugscore_diff"].max())
                                          if in23["max_abs_drugscore_diff"].notna().any()
                                          else None),
               "differing_models": in23[~in23["same_names"] |
                                        (in23["max_abs_drugscore_diff"] > 1e-4)]
                                   ["model"].tolist()[:30]}
    text = json.dumps(summary, indent=2)
    if args.json:
        Path(args.json).write_text(text + "\n")
    print(text)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--inputs", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--dogsite", required=True)
    r.add_argument("--tool", choices=("dogsite3", "dogsite2"), default="dogsite3")
    r.add_argument("--nprocs", type=int, default=os.cpu_count())
    r.add_argument("--timeout", type=float, default=6 * 3600)
    r.add_argument("--limit", type=int, default=0, help="random subset of N models")
    r.add_argument("--seed", type=int, default=42)
    a = sub.add_parser("analyze")
    a.add_argument("--ds2-drugscore", required=True)
    a.add_argument("--ds2-pockets", required=True)
    a.add_argument("--ds3", required=True)
    a.add_argument("--json")
    a.add_argument("--extra-cutoffs", type=float, nargs="*", default=[100, 150, 200, 250])
    c = sub.add_parser("compare-ds2")
    c.add_argument("--ds2-drugscore", required=True)
    c.add_argument("--rerun", required=True)
    c.add_argument("--json")
    args = p.parse_args(argv)
    {"run": cmd_run, "analyze": cmd_analyze, "compare-ds2": cmd_compare_ds2}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
