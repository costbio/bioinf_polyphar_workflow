"""Stage 07 (library): prepare the docking ligand library once, deterministically.

Input: the raw e-Drug3D SDF (`--library`, default inputs/e-Drug3D_2083.sdf — the same
option `dock` takes, and the same file).
Outputs in 07_docking/library/: library_prepared.sdf (the dockable records, minimised),
ligands.csv (one row per SOURCE record, kept or not, with the numbers the filter used),
failed.tsv, manifest.json (written last = done).

Why this is a stage of its own, and why it runs before `dock`:
preparation is deterministic and identical for every target, so doing it inside `dock`
would repeat it per target (and, with replicates, per replicate) and make the docking
inputs depend on when a target was docked. Prepared once and hashed, the library becomes a
fingerprinted input of `dock`: `dock` picks this file up automatically (see
`prepared_state` / `prepared_library` below, used by `dock._library_path`), so its existing
`library_md5` fingerprint covers the preparation, and nothing re-runs `uniprot` or
`pockets`.

Two different questions are answered by two different checks, on purpose. "Must the
library be prepared again?" is `fingerprint()` / `units()`: it hashes the source library
AND every preparation parameter, the RDKit version included, because a different RDKit
would produce different coordinates. "Which file should `dock` read?" is
`prepared_state()`: it asks only about the file's own identity (manifest present, the
prepared SDF still at its recorded sha256, the source still at the sha256 it was prepared
from). Using the first for the second would mean an rdkit-less `dock` host silently docks
the RAW library — and, because that changes `library_md5`, re-docks every target and
overwrites prepared results with unprepared ones. `prepared_state` also reports WHY it
declines, and `dock` warns loudly (stderr + notification) and records
`library_prepared: false` in its manifest rather than failing.

What it does, and why
---------------------
**Dockability filter: exclude `rotatable bonds > 15` OR `MW > 900`** (both bounds
exclusive: rotB == 15 and MW == 900.0 are KEPT). e-Drug3D's inclusion limit is MW <= 2000,
so it holds peptides far outside Vina's domain (afamelanotide at 64 rotatable bonds,
bacitracin at 1,421 Da); including them stalled a calibration run for hours. On
e-Drug3D_2083 this excludes 123 of 2,083 (5.9%) and keeps 1,960.

**The rotatable-bond definition is non-standard and is pinned here on purpose.** It is
`docs/superpowers/specs/2026-11-03-panel-selection-spikes/ligsize.py`, reproduced verbatim
in `bridges()` / `profile()` below: single, non-ring, non-terminal bonds between heavy
atoms, where non-ring is decided by **bridge detection on the heavy-atom bond graph**, and
**amide C-N bonds are NOT excluded**. RDKit's `CalcNumRotatableBonds` and `obabel` both
use the conventional Lipinski definition, which excludes amides and would silently change
which ligands are dropped: on e-Drug3D_2083 RDKit's count differs for 1,316 of 2,083
records and crosses the 15 boundary for 46 of them (SAQUINAVIR 16 -> 13, DIPYRIDAMOLE
12 -> 16, VANCOMYCIN 14 -> 23). The molecular weight is likewise the spike's: the sum of
its own element-mass table over every atom in the molfile (hydrogens included), unknown
elements counted as carbon — NOT RDKit's `MolWt`. That one is pinned for reproducibility,
not because it changes the result: swapping in RDKit's `MolWt` crosses the 900 boundary
for **0 of 2,083** records (largest difference 3.02 Da, a deuterium isotope). Both are
computed from the molfile text, with no RDKit involved, so the cutoffs mean exactly what
the design document measured.

**Conformers: minimise in place with MMFF94**, keeping the library's own graph and
stereochemistry. Vina-family docking samples torsions but not ring pucker or stereo
centres, so the input conformer is fixed for the whole run; starting from the library's
coordinates and relaxing them preserves the deposited stereochemistry, which regenerating
every ligand would not. "In place" means *same molecule, same starting point*, NOT a small
perturbation: over the 1,960 kept records the largest per-atom displacement has median
1.61 A, p99 8.27 A and maximum 12.63 A. The bond graph is preserved in all 1,960. One
conformer per ligand. The conformer is marked 3D on output: e-Drug3D's records are tagged
2D although their Z coordinates are non-zero, and that tag is metadata, not geometry.

**MMFF94 has no chirality restraint**, so on a badly strained deposited geometry it can
walk a near-planar stereocentre through its plane. The signed chirality volume of every
tetrahedral stereocentre is therefore measured before and after minimisation; a record
where any of them changes sign keeps its deposited geometry (the minimised conformer is
discarded), is counted as `stereo_inversions` and logged to failed.tsv. On e-Drug3D_2083
exactly one record does this: TRABECTEDIN (ID 1808), whose atom 18 goes S -> R (chirality
volume +1.861 -> -0.643) because its deposited geometry is strained (MMFF94 energy 647.5
-> 405.6 kcal/mol). Without the guard an approved antitumour agent would be docked as the
wrong diastereomer.

**Implausible formal charges are flagged, not swallowed.** A C/N/O atom carrying
|formal charge| > 2 is almost always a malformed `M  CHG` line, not chemistry: ISOSORBIDE
DINITRATE (ID 271) declares `M  CHG  5   1  15  13   1  14   1  17  -1  18  -1`, i.e. five
entries of which only the four trailing ones are real, so atom 1 — a carbon — is read as
formal charge +15. RDKit sanitises it and MMFF94 parameterises it, and the prepared record
would carry `[C@H+15]` into docking, where derived partial charges are nonsense. Such a
record is **passed through flagged** rather than dropped (dropping it would silently change
the library), counted as `implausible_formal_charge` and logged. It is the only kept record
of e-Drug3D_2083 with |net charge| > 4.

**Protonation (for Methods).** Ligands keep their deposited protonation, while receptors
are protonated at pH 7.0 by pdbfixer (`fix`). The asymmetry is deliberate: the library's
charge states are curated per drug, and re-deriving them at a nominal pH would overwrite
that curation. In the PREPARED library 430 of the 1,960 kept records carry a non-zero net
formal charge (-4: 24, -3: 3, -2: 63, -1: 276, +1: 55, +2: 8, and +15: the one corrupt
record above). Counted on the SOURCE molfiles it is 453, because the 23 half-written nitro
groups each read as net -1 there; the nitro repair restores them to neutral, so 430 is the
figure that describes what is docked.

A record RDKit cannot sanitise is first repaired, not regenerated: every sanitisation
failure in e-Drug3D_2083 (23 of 2,083 — metronidazole, nifedipine, clonazepam,
chloramphenicol, ...) is a **neutral five-valent nitro group**, `C-N(=O)=O`, which RDKit
rejects and which `charge_separate_hypervalent_nitrogen` rewrites as `C-[N+](=O)[O-]`
without moving an atom. ETKDG cannot fix those at all (it raises AtomValenceException on
all 23), so the design document's "regenerate with ETKDG for records failing sanitisation"
would have cost 23 approved drugs their minimisation. ETKDG remains the fallback for a
record that still fails after the repair, and every record taking either path is logged to
failed.tsv (counted as `charge_separated` / `etkdg`, not as `failed`: the ligand is
prepared either way). A record that cannot be prepared at all is written through
**verbatim** rather than dropped — dropping it would silently change the library — and is
counted as `failed` and logged. A molecule with no MMFF94 parameters keeps its input
geometry and is logged the same way; on e-Drug3D_2083 those are the 6 **boron-containing**
drugs (BORTEZOMIB, TAVABOROLE, IXAZOMIB, IXAZOMIB CITRATE, CRISABOROLE, VABORBACTAM) —
the library contains no metal-containing record at all.

RDKit is a hard requirement of this stage, checked once in `units()` / `run()`: without it
every record would "fail to prepare" and the stage would write an unprepared library while
reporting success.

Determinism: records are written in library order, every RDKit step is seeded or
parameter-free (MMFF94 minimisation has no RNG; ETKDG uses a fixed seed), and the original
`> <field>` data block of each record is copied through byte for byte, so `dock` still
finds `ID`, `name`, `cas` and `status`. Running the stage twice gives a byte-identical
file, and `--nprocs` changes nothing.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import json
import os
from pathlib import Path

from polyphar import pool
from polyphar.common import WorkUnits, append_failure, sha256_file, write_manifest

STAGE = "ligands"
DIRNAME = "library"                 # inside 07_docking (see RunLayout.docking)
OUTPUT = "library_prepared.sdf"
TABLE = "ligands.csv"
MAX_ROTATABLE_BONDS = 15            # exclude > this
MAX_MW = 900.0                      # exclude > this (Da)
ETKDG_SEED = 42                     # only used for records that fail sanitisation
MMFF_VARIANT = "MMFF94"
MMFF_MAX_ITERS = 2000
# |formal charge| above this on a C/N/O atom is a malformed `M  CHG` line, not chemistry
# (ISOSORBIDE DINITRATE, ID 271, reads as carbon +15). Flagged, counted, never silent.
MAX_ABS_FORMAL_CHARGE = 2
ORGANIC_Z = (6, 7, 8)

COLUMNS = ["record_index", "ID", "name", "heavy_atoms", "mw", "rotatable_bonds",
           "kept", "exclusion", "method", "mmff_converged", "mmff_status", "net_charge",
           "note"]

# mmff_status, kept apart from mmff_converged on purpose: the boolean cannot say whether
# MMFF94 ran at all (no parameters / a discarded conformer are not "did not converge").
MMFF_CONVERGED = "converged"
MMFF_NOT_CONVERGED = "not_converged"
MMFF_NO_PARAMETERS = "no_parameters"
MMFF_STEREO_INVERSION = "stereo_inversion"

RDKIT_REQUIRED = (
    "the ligands stage needs RDKit and this interpreter cannot import it ({why}). "
    "Without RDKit every record would be reported as 'could not prepare, source kept' and "
    "the stage would write an UNPREPARED library while claiming success. Run the stage "
    "inside the polyphar image, or with an interpreter that has rdkit.")

# Element masses of the spike (ligsize.py). Kept as it is on purpose: the 123/1,960 split
# the design document reports was measured with exactly this table.
MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998,
        "CL": 35.45, "BR": 79.904, "I": 126.90, "H": 1.008, "NA": 22.99, "K": 39.098,
        "CA": 40.078, "MG": 24.305, "FE": 55.845, "ZN": 65.38, "SE": 78.971, "B": 10.81,
        "SI": 28.085, "AS": 74.922, "LI": 6.94, "AL": 26.982, "GD": 157.25, "PT": 195.08,
        "TC": 98.0, "AU": 196.97, "AG": 107.87, "MN": 54.938, "CU": 63.546, "CR": 51.996,
        "CO": 58.933, "SR": 87.62, "BA": 137.33, "SM": 150.36, "RA": 226.0, "Y": 88.906,
        "IN": 114.82, "TL": 204.38, "BI": 208.98, "SB": 121.76, "HG": 200.59, "ST": 0.0}
UNKNOWN_ELEMENT_MASS = 12.0


def add_arguments(parser) -> None:
    # Defined exactly as `dock` defines it: `polyphar run` merges the stages' options and
    # rejects two stages declaring one option differently.
    from polyphar.stages.dock import DEFAULT_LIBRARY

    parser.add_argument("--library", default=DEFAULT_LIBRARY,
                        help=f"ligand SDF, relative to the run folder (default {DEFAULT_LIBRARY})")


# --- the pinned size / flexibility profile (no RDKit) ---------------------------------

def bridges(n: int, edges) -> set:
    """Indices of the bridge edges of the graph on nodes 1..n (iterative Tarjan).

    Verbatim from ligsize.py. A bond that is not a bridge is in a ring, so this is the
    'non-ring' half of the rotatable-bond definition."""
    adj = collections.defaultdict(list)
    for i, (a, b) in enumerate(edges):
        adj[a].append((b, i))
        adj[b].append((a, i))
    disc, low, out, t = {}, {}, set(), [0]
    for s in range(1, n + 1):
        if s in disc:
            continue
        stack = [(s, -1, iter(adj[s]))]
        disc[s] = low[s] = t[0]
        t[0] += 1
        while stack:
            u, pe, it = stack[-1]
            advanced = False
            for v, ei in it:
                if ei == pe:
                    continue
                if v not in disc:
                    disc[v] = low[v] = t[0]
                    t[0] += 1
                    stack.append((v, ei, iter(adj[v])))
                    advanced = True
                    break
                low[u] = min(low[u], disc[v])
            if not advanced:
                stack.pop()
                if stack:
                    p, ppe, _ = stack[-1]
                    low[p] = min(low[p], low[u])
                    if low[u] > disc[p]:
                        out.add(pe)
    return out


def profile(lines) -> dict:
    """{heavy_atoms, mw, rotatable_bonds} of one SDF record, from its molfile text.

    `lines` are the record's lines (line endings may be present), starting at the title
    line. Rotatable bonds are the spike's definition: single, non-ring (bridge),
    non-terminal heavy-atom bonds, amides included. Raises ValueError on a molfile whose
    counts line cannot be read."""
    lines = [ln.rstrip("\r\n") for ln in lines]
    if len(lines) < 4:
        raise ValueError("record too short to hold a molfile")
    try:
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    except ValueError:
        raise ValueError(f"unreadable counts line: {lines[3]!r}") from None

    el, mw = [], 0.0
    for ln in lines[4:4 + na]:
        e = (ln[31:34].strip() or "C").upper()
        el.append(e)
        mw += MASS.get(e, UNKNOWN_ELEMENT_MASS)
    heavy_idx = {i + 1 for i, e in enumerate(el) if e != "H"}

    edges, orders = [], []
    for ln in lines[4 + na:4 + na + nb]:
        try:
            a1, a2, o = int(ln[0:3]), int(ln[3:6]), int(ln[6:9])
        except ValueError:
            continue
        if a1 in heavy_idx and a2 in heavy_idx:
            edges.append((a1, a2))
            orders.append(o)

    br = bridges(na, edges)
    deg = collections.Counter()
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    rot = sum(1 for i, (a, b) in enumerate(edges)
              if orders[i] == 1 and i in br and deg[a] > 1 and deg[b] > 1)
    return {"heavy_atoms": len(heavy_idx), "mw": mw, "rotatable_bonds": rot}


def exclusion(prof: dict) -> str:
    """"" when the record is dockable, else why it is not (both bounds exclusive)."""
    reasons = []
    if prof["rotatable_bonds"] > MAX_ROTATABLE_BONDS:
        reasons.append(f"rotatable_bonds={prof['rotatable_bonds']}>{MAX_ROTATABLE_BONDS}")
    if prof["mw"] > MAX_MW:
        reasons.append(f"mw={prof['mw']:.1f}>{MAX_MW:g}")
    return "; ".join(reasons)


# --- conformer preparation (RDKit) ----------------------------------------------------

def _molblock(lines) -> str:
    """The molfile part of a record (title line .. 'M  END')."""
    out = []
    for ln in lines:
        out.append(ln)
        if ln.rstrip("\r\n") == "M  END":
            break
    return "".join(out)


def _data_block(lines) -> str:
    """The record's `> <field>` data items, verbatim, without the $$$$ line."""
    out, started = [], False
    for ln in lines:
        if ln.rstrip("\r\n") == "$$$$":
            break
        if not started and ln.startswith(">"):
            started = True
        if started:
            out.append(ln)
    return "".join(out)


def charge_separate_hypervalent_nitrogen(mol) -> int:
    """Complete the charge-separated form of an over-valent neutral nitrogen, in place.

    Every sanitisation failure in e-Drug3D_2083 is a nitro group whose charges are only
    half written: the molfile carries `M  CHG` for the anionic oxygen but leaves the
    nitrogen neutral, so RDKit sees a four-valent neutral N and refuses the record
    (metronidazole, nifedipine, chloramphenicol, ... — 23 of 2,083). Setting the nitrogen
    to +1 completes the standard `[N+](=O)[O-]` form. A fully neutral `N(=O)=O`, which
    RDKit's molfile parser normalises by itself, is handled too: one N=O becomes
    [N+]-[O-].

    Only the valence model changes — no atom moves — so the deposited conformer and its
    stereochemistry survive. Returns the number of nitrogens repaired.

    The gate deliberately refuses to judge an aromatic nitrogen. `GetBondTypeAsDouble()`
    scores an aromatic bond 1.5, so an aromatic ring N bearing an anionic O (the HOBt-anion
    motif) sums to 4.0 although its true valence is 3, and a naive `> 3` test would wrongly
    set it to +1. e-Drug3D_2083 has no bond-type-4 bond at all, so nothing here is reached
    by it today; the guard is against the next library."""
    from rdkit import Chem

    fixed = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 0:
            continue
        bonds = list(atom.GetBonds())
        if any(b.GetBondType() == Chem.BondType.AROMATIC or b.GetIsAromatic()
               for b in bonds):
            continue                      # 1.5 per bond: the sum below means nothing
        valence = sum(b.GetBondTypeAsDouble() for b in bonds) + atom.GetNumExplicitHs()
        if valence <= 3 or valence != int(valence):
            continue                      # non-integral: not a Kekule valence either
        anionic = [b for b in atom.GetBonds()
                   if b.GetOtherAtom(atom).GetAtomicNum() == 8
                   and b.GetOtherAtom(atom).GetFormalCharge() == -1]
        if anionic:                      # the oxygen is already [O-]: complete the pair
            atom.SetFormalCharge(1)
            fixed += 1
            continue
        for bond in atom.GetBonds():     # fully neutral nitro: separate the charges
            other = bond.GetOtherAtom(atom)
            if (bond.GetBondType() == Chem.BondType.DOUBLE and other.GetAtomicNum() == 8
                    and other.GetDegree() == 1 and other.GetFormalCharge() == 0):
                bond.SetBondType(Chem.BondType.SINGLE)
                other.SetFormalCharge(-1)
                atom.SetFormalCharge(1)
                fixed += 1
                break
    return fixed


def _sanitised(molblock: str):
    """(mol, needs_etkdg, note).

    A molblock RDKit sanitises directly gives (mol, False, ""). One it does not is parsed
    unsanitised and repaired by charge_separate_hypervalent_nitrogen, which is what every
    sanitisation failure in e-Drug3D_2083 actually is (23 of 2,083: metronidazole,
    nifedipine, chloramphenicol, ...), and keeps its deposited geometry. Only a record
    that still fails then is partially sanitised and handed to ETKDG, because its
    coordinates cannot be trusted once the valence model is broken."""
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=True)
    if mol is not None:
        return mol, False, ""
    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("RDKit cannot parse the molblock")
    mol.UpdatePropertyCache(strict=False)
    fixed = charge_separate_hypervalent_nitrogen(mol)
    if fixed:
        mol.UpdatePropertyCache(strict=False)
        if not Chem.SanitizeMol(mol, catchErrors=True):
            return mol, False, (f"{fixed} neutral hypervalent N rewritten "
                                f"charge-separated; geometry kept")
    # Last ditch: sanitise as far as RDKit can (valence, kekulisation and aromaticity
    # perception dropped) and let ETKDG rebuild the conformer. A molecule whose valence
    # model is still wrong makes ETKDG raise, and the caller then keeps the source record.
    ops = (Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
           ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    if Chem.SanitizeMol(mol, sanitizeOps=ops, catchErrors=True):
        raise ValueError("RDKit cannot sanitise the molblock")
    return mol, True, ""


def _embed(mol) -> None:
    """A fresh ETKDG conformer, seeded. Only for records that failed sanitisation."""
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.randomSeed = ETKDG_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("ETKDG could not embed the molecule")


def chirality_volumes(mol) -> dict:
    """{atom index: signed chirality volume} at every tetrahedral stereocentre.

    The volume is the scalar triple product of the vectors FROM THE CENTRAL ATOM to its
    three lowest-indexed neighbours. That origin is not a detail: it is the convention
    RDKit's own `AssignChiralTypesFrom3D` uses, and it is what tracks an umbrella
    inversion, where the centre pops through the plane of its substituents. The signed
    volume of the neighbours' own tetrahedron does NOT: on TRABECTEDIN it keeps its sign
    (+4.93 -> +6.53) while RDKit's perceived tag goes CCW -> CW and this measure goes
    +1.861 -> -0.643.

    Neighbours are ordered by atom index, which `MolToMolBlock` preserves, so the two
    conformers of one molecule are directly comparable; the SIGN is the handedness and an
    inversion flips it. Stereocentres are perceived from the 3D coordinates on a copy, so
    nothing on `mol` changes."""
    from rdkit import Chem

    probe = Chem.Mol(mol)
    try:
        Chem.AssignStereochemistryFrom3D(probe)
    except Exception:
        return {}
    conf = mol.GetConformer()
    out = {}
    for atom in probe.GetAtoms():
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        nbrs = sorted(n.GetIdx() for n in atom.GetNeighbors())
        if len(nbrs) < 3:
            continue
        origin = conf.GetAtomPosition(atom.GetIdx())
        a, b, c = (conf.GetAtomPosition(i) - origin for i in nbrs[:3])
        out[atom.GetIdx()] = float(a.CrossProduct(b).DotProduct(c))
    return out


def inverted_stereocentres(before: dict, after: dict) -> list:
    """Atom indices whose signed chirality volume changed sign between two conformers."""
    return sorted(i for i, v in before.items()
                  if i in after and v * after[i] < 0.0)


def implausible_charges(mol) -> list:
    """(index, symbol, charge) of C/N/O atoms carrying |formal charge| > MAX_ABS...

    A malformed `M  CHG` line, not chemistry: ISOSORBIDE DINITRATE (ID 271) declares five
    entries and lists four, so its first field reads as atom 1 (a carbon) at +15. RDKit
    sanitises and MMFF94 parameterises it happily, so nothing downstream would notice."""
    return [(a.GetIdx(), a.GetSymbol(), a.GetFormalCharge()) for a in mol.GetAtoms()
            if a.GetAtomicNum() in ORGANIC_Z
            and abs(a.GetFormalCharge()) > MAX_ABS_FORMAL_CHARGE]


def _minimise(mol) -> tuple[str, str]:
    """MMFF94 in place. (status, note).

    `status` is one of MMFF_CONVERGED / MMFF_NOT_CONVERGED / MMFF_NO_PARAMETERS /
    MMFF_STEREO_INVERSION. No parameters leaves the geometry alone. So does a minimisation
    that inverts a stereocentre: MMFF94 carries no chirality restraint, so from a strained
    deposited geometry it can push a near-planar centre through its plane (TRABECTEDIN,
    ID 1808 — the one record of e-Drug3D_2083 where it happens). The minimised conformer
    is then discarded and the deposited one restored, and the caller logs and counts it."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if not AllChem.MMFFHasAllMoleculeParams(mol):
        return MMFF_NO_PARAMETERS, f"no {MMFF_VARIANT} parameters; input geometry kept"
    before = chirality_volumes(mol)
    deposited = Chem.Conformer(mol.GetConformer())
    code = AllChem.MMFFOptimizeMolecule(mol, mmffVariant=MMFF_VARIANT,
                                        maxIters=MMFF_MAX_ITERS)
    flipped = inverted_stereocentres(before, chirality_volumes(mol))
    if flipped:
        mol.RemoveAllConformers()
        mol.AddConformer(deposited, assignId=True)
        return MMFF_STEREO_INVERSION, (
            f"{MMFF_VARIANT} inverted the stereocentre(s) at atom(s) "
            f"{', '.join(str(i) for i in flipped)}; minimised conformer discarded, "
            f"deposited geometry kept")
    if code == 0:
        return MMFF_CONVERGED, ""
    return (MMFF_NOT_CONVERGED,
            f"{MMFF_VARIANT} did not converge in {MMFF_MAX_ITERS} iterations")


def prepare_record(lines) -> dict:
    """Prepare one kept record.

    Returns {text, method, mmff_converged, mmff_status, net_charge, bad_charges, note}.
    `text` is the record ready for docking: the minimised molfile followed by the source
    record's data items verbatim and `$$$$`. Raises ValueError when the record cannot be
    prepared at all (the caller then keeps the source record)."""
    from rdkit import Chem

    mol, needs_etkdg, note = _sanitised(_molblock(lines))
    method = f"{MMFF_VARIANT.lower()}"
    if needs_etkdg:
        _embed(mol)
        method = f"etkdg+{MMFF_VARIANT.lower()}"
        note = "sanitisation failed; conformer regenerated with ETKDG"
    if mol.GetNumConformers() == 0:
        raise ValueError("no conformer")
    status, mmff_note = _minimise(mol)
    if mmff_note:
        note = f"{note}; {mmff_note}" if note else mmff_note
    if status in (MMFF_NO_PARAMETERS, MMFF_STEREO_INVERSION):
        # the written geometry is the deposited one, not an MMFF94 one: say so in `method`
        # so the `minimised` count means "this conformer came out of MMFF94".
        method = method.replace(MMFF_VARIANT.lower(), "none")
    bad_charges = implausible_charges(mol)
    if bad_charges:
        where = ", ".join(f"atom {i} ({sym}) {chg:+d}" for i, sym, chg in bad_charges)
        flag = (f"implausible formal charge(s) carried through from the molfile: {where}; "
                f"likely a malformed `M  CHG` line — docking partial charges will be wrong")
        note = f"{note}; {flag}" if note else flag
    # e-Drug3D tags its records 2D although the Z coordinates are non-zero: a metadata
    # error, not a geometry one. Docking engines read the flag, so correct it.
    mol.GetConformer().Set3D(True)
    block = Chem.MolToMolBlock(mol, includeStereo=True, kekulize=True)
    if not block.endswith("\n"):
        block += "\n"
    return {"text": block + _data_block(lines) + "$$$$\n", "method": method,
            "mmff_converged": status == MMFF_CONVERGED if status in (
                MMFF_CONVERGED, MMFF_NOT_CONVERGED) else "",
            "mmff_status": status, "bad_charges": bad_charges,
            "net_charge": Chem.GetFormalCharge(mol), "note": note}


def _silence_rdkit_warnings() -> None:
    """Drop RDKit's per-record warnings (every e-Drug3D record is 'tagged as 2D, but at
    least one Z coordinate is not zero'): 2,083 of them bury the stage's own output.
    Errors still come through, and every record's outcome is in ligands.csv anyway."""
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.info")
    except ImportError:
        pass


def _prepare_one(index: int, lines: list) -> dict:
    """One record: profile, filter, prepare. Never raises for a record that cannot be
    prepared — the source text is returned instead, flagged.

    A missing RDKit is NOT a per-record failure and is not caught here: `units()` / `run()`
    refuse to start without it, so reaching this function means the import works."""
    _silence_rdkit_warnings()
    prof = profile(lines)
    why = exclusion(prof)
    row = {"record_index": index, **prof, "kept": not why, "exclusion": why,
           "method": "", "mmff_converged": "", "mmff_status": "", "net_charge": "",
           "note": ""}
    if why:
        return {"row": row, "text": None, "failed": None, "logged": None,
                "bad_charges": []}
    try:
        prepared = prepare_record(lines)
    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        row.update(method="source", note=f"not prepared: {reason}; source record kept")
        return {"row": row, "text": "".join(lines), "failed": reason, "logged": reason,
                "bad_charges": []}
    row.update(method=prepared["method"], mmff_converged=prepared["mmff_converged"],
               mmff_status=prepared["mmff_status"], net_charge=prepared["net_charge"],
               note=prepared["note"])
    return {"row": row, "text": prepared["text"], "failed": None,
            "logged": prepared["note"] or None, "bad_charges": prepared["bad_charges"]}


def _worker(index: int, lines: list) -> dict:
    """`_prepare_one` inside pool.task(), like every other stage's pool worker (dock's
    _unit_worker, fix.py, similarity.py): failures cross the process boundary as
    TaskError and SIGTERM ends the worker instead of returning it to the idle loop."""
    return pool.task(lambda: _prepare_one(index, lines))


# --- stage ----------------------------------------------------------------------------

def stage_dir(layout) -> Path:
    return layout.docking / DIRNAME


def _source(layout, args, copy: bool = True) -> Path:
    """The raw library, never the prepared one."""
    from polyphar.stages.dock import _raw_library_path

    return _raw_library_path(layout, args, copy=copy)


def _rdkit_version() -> str:
    try:
        import rdkit

        return rdkit.__version__
    except Exception:
        return "unavailable"


def _versions() -> dict:
    return {"rdkit": _rdkit_version(), "mmff_variant": MMFF_VARIANT,
            "mmff_max_iters": MMFF_MAX_ITERS, "etkdg_seed": ETKDG_SEED,
            "max_rotatable_bonds": MAX_ROTATABLE_BONDS, "max_mw": MAX_MW,
            "max_abs_formal_charge": MAX_ABS_FORMAL_CHARGE,
            # the stereo guard decides which geometry is written, so it is a preparation
            # parameter: changing it must re-prepare the library.
            "stereo_inversion_guard": "signed chirality volume, central-atom origin",
            "rotatable_bond_definition":
                "single, non-ring (bridge), non-terminal heavy-atom bonds, amides included "
                "(ligsize.py); NOT RDKit/Lipinski"}


def fingerprint(source: Path) -> str:
    """What the prepared library depends on: the source library and every preparation
    parameter (not --nprocs, which cannot change the result)."""
    payload = json.dumps({"source_sha256": sha256_file(source), **_versions()},
                         sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _manifest(layout) -> dict | None:
    try:
        data = json.loads((stage_dir(layout) / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _finished(layout, source: Path) -> dict | None:
    """The manifest's counts when the prepared library is current for `source`, else None.

    This is the "must it be prepared AGAIN?" question, so it compares the full
    `fingerprint()` — the RDKit version and every preparation parameter included. It is
    NOT the question `dock` asks (see `prepared_state`)."""
    manifest = _manifest(layout)
    if manifest is None or manifest.get("fingerprint") != fingerprint(source):
        return None
    for name in (OUTPUT, TABLE):
        if not (stage_dir(layout) / name).is_file():
            return None
    counts = manifest.get("counts")
    return dict(counts) if isinstance(counts, dict) else {}


def _recorded_source_sha(manifest: dict) -> str | None:
    """The sha256 of the library this manifest was prepared from."""
    recorded = manifest.get("source_sha256")
    if isinstance(recorded, str) and recorded:
        return recorded
    for item in manifest.get("inputs") or []:    # manifests written before the key existed
        if isinstance(item, dict) and item.get("sha256"):
            return str(item["sha256"])
    return None


def prepared_state(layout, source) -> tuple:
    """(prepared library path or None, why it cannot be used or None).

    The question `dock` asks — *which file do I dock?* — and nothing else. It is decided
    on the files' own identity: a `ligands` manifest, the prepared SDF present and still
    at its recorded sha256, and the source library still at the sha256 it was prepared
    from. Deliberately NO RDKit version and NO preparation-parameter check: those belong
    to `fingerprint()` / `units()`, which decide whether to prepare again. Hashing the
    RDKit version here made an rdkit-less (or differently-versioned) `dock` process fall
    back to the RAW library, which changes `library_md5`, re-docks every target and
    overwrites prepared results with unprepared ones.

    The reason is None when there is simply nothing prepared here (no manifest: not a
    warning), and a message when a manifest IS present but the prepared library still
    cannot be used — the states `dock` must shout about."""
    source = Path(source)
    manifest = _manifest(layout)
    if manifest is None:
        return None, None
    out = stage_dir(layout) / OUTPUT
    if not source.is_file():
        return None, f"the source library {source} does not exist"
    if not out.is_file():
        return None, (f"the ligands manifest is there but {out.name} is not "
                      f"(expected at {out}): the prepared library was deleted or moved")
    recorded_out = manifest.get("library_sha256")
    recorded_src = _recorded_source_sha(manifest)
    if not recorded_out or not recorded_src:
        return None, (f"the ligands manifest at {stage_dir(layout)} records no "
                      f"library_sha256/source_sha256: it is too old or truncated")
    try:
        actual_out, actual_src = sha256_file(out), sha256_file(source)
    except OSError as exc:
        return None, f"cannot hash the libraries: {exc}"
    if actual_out != recorded_out:
        return None, (f"{out.name} changed after it was prepared (sha256 "
                      f"{actual_out[:12]}, manifest says {recorded_out[:12]})")
    if actual_src != recorded_src:
        return None, (f"the source library {source.name} changed after the prepared one "
                      f"was built (sha256 {actual_src[:12]}, manifest says "
                      f"{recorded_src[:12]}): rerun `polyphar ligands`")
    return out, None


def prepared_library(layout, source) -> Path | None:
    """The prepared library for `source`, or None. See `prepared_state`."""
    try:
        return prepared_state(layout, source)[0]
    except OSError:
        return None


def rdkit_problem() -> str | None:
    """Why this stage cannot run in this interpreter, or None.

    A missing RDKit is one missing dependency, not 2,083 per-record failures: caught per
    record it would classify every ligand as "could not prepare, source kept", write an
    unprepared library and report success, and `dock` would pick that up as prepared."""
    try:
        import rdkit                      # noqa: F401
        from rdkit import Chem            # noqa: F401
        from rdkit.Chem import AllChem    # noqa: F401
    except Exception as exc:
        return RDKIT_REQUIRED.format(why=f"{type(exc).__name__}: {exc}")
    return None


def _require_rdkit() -> None:
    problem = rdkit_problem()
    if problem:
        raise RuntimeError(problem)


def _records(source: Path) -> list:
    """The source records' line lists, in library order (IDs validated by dock)."""
    from polyphar.stages.dock import _library_records

    return [rec.lines for rec in _library_records(source)]


def _fields(lines) -> dict:
    from polyphar.stages.dock import _records as split

    recs = split("".join(lines))
    return recs[0].fields if recs else {}


def units(layout, args) -> WorkUnits:
    """Records to prepare; 0 once the prepared library is current for this source."""
    source = _source(layout, args, copy=False)
    if _finished(layout, source) is not None:
        return WorkUnits(0, "ligands")
    _require_rdkit()   # there is work to do, and it cannot be done without RDKit
    records = _records(source)
    return WorkUnits(len(records), "ligands",
                     [float(sum(len(ln) for ln in r)) for r in records])


def run(layout, args, progress) -> dict:
    source = _source(layout, args)
    done = _finished(layout, source)
    if done is not None:   # current for this source: touch nothing
        return done
    _require_rdkit()

    d = stage_dir(layout)
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # the outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")

    records = _records(source)
    fields = [_fields(lines) for lines in records]
    sizes = [float(sum(len(ln) for ln in r)) for r in records]
    results: dict[int, dict] = {}
    failures: list[tuple[str, str]] = []

    def on_result(item, result):
        results[item] = result
        progress.advance(1, size=sizes[item])

    def on_failure(item, why):
        failures.append((f"record {item}", why))
        progress.advance(1, size=sizes[item])

    jobs = [(i, _worker, (i, list(lines))) for i, lines in enumerate(records)]
    try:
        pool.run_pool(jobs, max(1, getattr(args, "nprocs", 1)),
                      on_result=on_result, on_failure=on_failure)
    finally:
        for item, reason in sorted(failures):
            append_failure(d, item, STAGE, reason)

    # Written in library order, so neither --nprocs nor the pool's completion order can
    # change a byte of the output.
    rows, kept_text, logged = [], [], []
    n_failed = n_etkdg = n_no_mmff = n_charge_separated = 0
    n_stereo = n_bad_charge = 0
    for i, lines in enumerate(records):
        got = results.get(i)
        ident = fields[i].get("ID", "") or f"record_{i}"
        name = fields[i].get("name", "")
        if got is None:     # the worker itself died: keep the source record, flag it
            rows.append({"record_index": i, "ID": ident, "name": name, "heavy_atoms": "",
                         "mw": "", "rotatable_bonds": "", "kept": True, "exclusion": "",
                         "method": "source", "mmff_converged": "", "mmff_status": "",
                         "net_charge": "",
                         "note": "worker failed; source record kept"})
            kept_text.append("".join(lines))
            n_failed += 1
            continue
        row = got["row"]
        rows.append({**row, "ID": ident, "name": name,
                     "mw": round(row["mw"], 3) if row["mw"] != "" else ""})
        if got["text"] is not None:
            kept_text.append(got["text"])
        if got["failed"]:
            n_failed += 1
        if "etkdg" in (row.get("method") or ""):
            n_etkdg += 1
        if row.get("mmff_status") == MMFF_NO_PARAMETERS:
            n_no_mmff += 1
        if row.get("mmff_status") == MMFF_STEREO_INVERSION:
            n_stereo += 1
        if "charge-separated" in (row.get("note") or ""):
            n_charge_separated += 1
        if got["bad_charges"]:
            n_bad_charge += 1
        if got["logged"]:
            logged.append((f"{ident} (record {i})", got["logged"]))

    for item, reason in logged:
        append_failure(d, item, STAGE, reason)

    _atomic_write(d / OUTPUT, "".join(kept_text))
    _write_table(d / TABLE, rows)

    kept = sum(1 for r in rows if r["kept"])
    counts = {
        "records": len(records), "kept": kept, "excluded": len(records) - kept,
        "excluded_rotatable_bonds": sum(1 for r in rows
                                        if "rotatable_bonds=" in (r["exclusion"] or "")),
        "excluded_mw": sum(1 for r in rows if "mw=" in (r["exclusion"] or "")),
        # `minimised` means "this conformer came out of MMFF94": a record whose minimised
        # conformer was discarded (stereo inversion) or never built (no parameters) keeps
        # its deposited geometry and is NOT counted here.
        "minimised": sum(1 for r in rows if r["kept"] and r["method"].endswith(
            MMFF_VARIANT.lower())),
        "etkdg": n_etkdg, "charge_separated": n_charge_separated,
        "no_mmff_parameters": n_no_mmff, "stereo_inversions": n_stereo,
        "implausible_formal_charge": n_bad_charge, "failed": n_failed,
    }
    # The manifest is written last: its presence (with this fingerprint) marks the stage
    # as finished, and library_sha256 is the hash `dock` fingerprints the library by.
    write_manifest(d, STAGE,
                   args={"nprocs": getattr(args, "nprocs", 1),
                         "library": str(source)},
                   inputs=[source], versions=_versions(), counts=counts, timing={},
                   extra={"fingerprint": fingerprint(source),
                          "library_sha256": sha256_file(d / OUTPUT),
                          # what `prepared_state` (and so `dock`) checks: the identity of
                          # the two files, independent of the RDKit that made them.
                          "source_sha256": sha256_file(source),
                          "outputs": {name: (d / name).stat().st_size
                                      for name in (OUTPUT, TABLE)}})
    return counts


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
        fh.write(text)
    os.replace(tmp, path)


def _write_table(path: Path, rows) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in COLUMNS})
    os.replace(tmp, path)
