"""Tests for the ligands stage: the pinned rotatable-bond definition, the dockability
filter and its boundaries, conformer preparation and its fallbacks, determinism, resume,
and the wiring that lets `dock` pick the prepared library up.

The profile half (rotatable bonds, MW, the filter) is pure Python and always runs; the
preparation half needs RDKit and is marked `rdkit`.
"""
import argparse
import csv
import json
import types

import pytest

from polyphar.stages import ligands

# --- fixtures: hand-written V2000 molfiles ------------------------------------------


def _molfile(title, atoms, bonds, fields=(), charges=(), end=True):
    """atoms = [(element, x, y, z)], bonds = [(a1, a2, order)] (1-based), charges =
    [(atom, charge)] written as an `M  CHG` line. One SDF record."""
    out = [f"{title}\n", "  polyphar test\n", "\n",
           f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0  0  0  0  0999 V2000\n"]
    for el, x, y, z in atoms:
        out.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {el:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n")
    for a, b, order in bonds:
        out.append(f"{a:>3}{b:>3}{order:>3}  0  0  0  0\n")
    for atom, charge in charges:
        out.append(f"M  CHG  1 {atom:>3} {charge:>3}\n")
    out.append("M  END\n")
    for key, value in fields:
        out.append(f"> <{key}>\n{value}\n\n")
    if end:
        out.append("$$$$\n")
    return "".join(out)


def _chain(elements, title="x", fields=()):
    """A straight single-bonded chain, 1.5 A apart along x."""
    atoms = [(el, 1.5 * i, 0.0, 0.0) for i, el in enumerate(elements)]
    bonds = [(i, i + 1, 1) for i in range(1, len(elements))]
    return _molfile(title, atoms, bonds, fields)


N_METHYLACETAMIDE = _molfile(
    "N-methylacetamide",
    [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("O", 2.1, 1.1, 0.0),
     ("N", 2.1, -1.2, 0.0), ("C", 3.5, -1.4, 0.0)],
    [(1, 2, 1), (2, 3, 2), (2, 4, 1), (4, 5, 1)],
    [("ID", "AMIDE1"), ("name", "N-METHYLACETAMIDE")])

CYCLOHEXANE = _molfile(
    "cyclohexane",
    [("C", 1.5, 0.0, 0.0), ("C", 0.75, 1.3, 0.0), ("C", -0.75, 1.3, 0.0),
     ("C", -1.5, 0.0, 0.0), ("C", -0.75, -1.3, 0.0), ("C", 0.75, -1.3, 0.0)],
    [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1), (6, 1, 1)],
    [("ID", "RING1"), ("name", "CYCLOHEXANE")])

BUTANE = _chain(["C", "C", "C", "C"], "butane",
                [("ID", "BUT1"), ("name", "BUTANE")])

# A half-charged nitro group -- the anionic oxygen is marked, the nitrogen is left
# neutral and so has valence 4. This is exactly how e-Drug3D writes its nitro drugs and
# is the only sanitisation failure the real library has (23 of 2,083, checked against
# metronidazole record 301). RDKit normalises a FULLY neutral N(=O)=O by itself, so that
# form would not reproduce the failure.
NITROMETHANE = _molfile(
    "nitromethane",
    [("C", 0.0, 0.0, 0.0), ("N", 1.5, 0.0, 0.0), ("O", 2.1, 1.1, 0.0),
     ("O", 2.1, -1.1, 0.0)],
    [(1, 2, 1), (2, 3, 2), (2, 4, 1)],
    [("ID", "NITRO1"), ("name", "NITROMETHANE")],
    charges=[(4, -1)])


def _lines(record_text):
    return record_text.splitlines(keepends=True)


# --- the pinned rotatable-bond definition -------------------------------------------

def test_rotatable_bonds_counts_the_amide_bond_rdkit_would_not():
    """The whole point of pinning the definition: ligsize.py counts the amide C-N, the
    conventional RDKit/Lipinski definition does not. If this ever flips, the dockability
    filter silently selects a different library."""
    assert ligands.profile(_lines(N_METHYLACETAMIDE))["rotatable_bonds"] == 1


@pytest.mark.rdkit
def test_the_conventional_definition_really_disagrees_on_the_amide():
    """Guards the claim above against the actual RDKit descriptor."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromMolBlock(ligands._molblock(_lines(N_METHYLACETAMIDE)))
    assert rdMolDescriptors.CalcNumRotatableBonds(mol) == 0
    assert ligands.profile(_lines(N_METHYLACETAMIDE))["rotatable_bonds"] == 1


@pytest.mark.parametrize("record,expected", [
    (BUTANE, 1),             # only the middle bond: the two end bonds are terminal
    (CYCLOHEXANE, 0),        # every bond is in the ring (no bridges)
    (N_METHYLACETAMIDE, 1),  # the amide C-N
])
def test_rotatable_bond_count_matches_the_spike(record, expected):
    assert ligands.profile(_lines(record))["rotatable_bonds"] == expected


@pytest.mark.repo_files
def test_rotatable_bonds_match_ligsize_on_every_fixture():
    """Run the spike's own code over the same records and require identical numbers."""
    import importlib.util
    from pathlib import Path

    spike = (Path(__file__).resolve().parents[1] / "docs" / "superpowers" / "specs"
             / "2026-11-03-panel-selection-spikes" / "ligsize.py")
    spec = importlib.util.spec_from_file_location("ligsize", spike)
    module = importlib.util.module_from_spec(spec)
    # ligsize.py runs its report at import time; only its bridges() is reusable, and the
    # rest of the definition is re-stated here exactly as the spike states it.
    source = spike.read_text().split("mols = []")[0]
    exec(compile(source, str(spike), "exec"), module.__dict__)  # noqa: S102

    import collections

    for record in (BUTANE, CYCLOHEXANE, N_METHYLACETAMIDE, NITROMETHANE):
        b = [ln for ln in _lines(record) if ln.rstrip("\r\n") != "$$$$"]
        na, nb = int(b[3][0:3]), int(b[3][3:6])
        el = [(ln[31:34].strip() or "C").upper() for ln in b[4:4 + na]]
        heavy_idx = {i + 1 for i, e in enumerate(el) if e != "H"}
        edges, orders = [], []
        for ln in b[4 + na:4 + na + nb]:
            a1, a2, o = int(ln[0:3]), int(ln[3:6]), int(ln[6:9])
            if a1 in heavy_idx and a2 in heavy_idx:
                edges.append((a1, a2))
                orders.append(o)
        br = module.bridges(na, edges)
        deg = collections.Counter()
        for a, c in edges:
            deg[a] += 1
            deg[c] += 1
        rot = sum(1 for i, (a, c) in enumerate(edges)
                  if orders[i] == 1 and i in br and deg[a] > 1 and deg[c] > 1)
        mw = sum(module.MASS.get(e, 12.0) for e in el)
        got = ligands.profile(_lines(record))
        assert got["rotatable_bonds"] == rot
        assert got["mw"] == pytest.approx(mw)
        assert ligands.bridges(na, edges) == br


def test_molecular_weight_uses_the_spike_table_over_every_atom():
    # a C-H chain: the table's own masses, hydrogens included
    record = _chain(["C", "C", "H"])
    assert ligands.profile(_lines(record))["mw"] == pytest.approx(2 * 12.011 + 1.008)
    assert ligands.profile(_lines(record))["heavy_atoms"] == 2


def test_unknown_elements_are_counted_as_carbon():
    record = _chain(["C", "XX"])
    assert ligands.profile(_lines(record))["mw"] == pytest.approx(12.011 + 12.0)


# --- the dockability filter ---------------------------------------------------------

@pytest.mark.parametrize("rot,mw,excluded", [
    (15, 900.0, False),    # both bounds are exclusive: the boundary is KEPT
    (16, 900.0, True),
    (15, 900.1, True),
    (0, 0.0, False),
    (64, 1646.0, True),    # afamelanotide
])
def test_filter_boundaries_are_exclusive(rot, mw, excluded):
    why = ligands.exclusion({"rotatable_bonds": rot, "mw": mw})
    assert bool(why) is excluded


def test_filter_names_both_reasons_when_both_apply():
    why = ligands.exclusion({"rotatable_bonds": 35, "mw": 1421.0})
    assert "rotatable_bonds=35>15" in why and "mw=1421.0>900" in why


def test_filter_constants_are_the_designed_ones():
    assert ligands.MAX_ROTATABLE_BONDS == 15
    assert ligands.MAX_MW == 900.0


# --- preparation --------------------------------------------------------------------

@pytest.mark.rdkit
def test_prepare_keeps_the_data_fields_verbatim_and_marks_the_conformer_3d():
    got = ligands.prepare_record(_lines(BUTANE))
    assert "> <ID>\nBUT1\n" in got["text"]
    assert "> <name>\nBUTANE\n" in got["text"]
    assert got["text"].endswith("$$$$\n")
    assert got["method"] == "mmff94"
    # the counts/dimensionality line of the written molfile says 3D
    assert "3D" in got["text"].splitlines()[1]


@pytest.mark.rdkit
def test_prepare_moves_atoms_but_keeps_the_atom_count():
    before = ligands.profile(_lines(BUTANE))
    after = ligands.profile(ligands.prepare_record(_lines(BUTANE))["text"].splitlines(
        keepends=True))
    assert after["heavy_atoms"] == before["heavy_atoms"]
    assert after["rotatable_bonds"] == before["rotatable_bonds"]


@pytest.mark.rdkit
def test_a_neutral_nitro_group_is_repaired_not_regenerated():
    """The real sanitisation failure in e-Drug3D. ETKDG cannot embed these at all, so
    the repair must happen before the fallback or 23 approved drugs lose minimisation."""
    from rdkit import Chem

    assert Chem.MolFromMolBlock(ligands._molblock(_lines(NITROMETHANE))) is None
    got = ligands.prepare_record(_lines(NITROMETHANE))
    assert got["method"] == "mmff94"            # geometry kept, not re-embedded
    assert "charge-separated" in got["note"]
    assert "> <ID>\nNITRO1\n" in got["text"]


@pytest.mark.rdkit
def test_etkdg_is_the_fallback_for_a_record_that_cannot_be_sanitised(monkeypatch):
    """A record the repair does not rescue gets a fresh ETKDG conformer instead of its
    deposited one, and the method says so (the stage logs exactly these to failed.tsv)."""
    real = ligands._sanitised

    def unsanitisable(molblock):
        mol, _needs, _note = real(molblock)
        return mol, True, ""       # as if sanitisation had failed

    monkeypatch.setattr(ligands, "_sanitised", unsanitisable)
    before = ligands.profile(_lines(BUTANE))
    got = ligands.prepare_record(_lines(BUTANE))
    assert got["method"] == "etkdg+mmff94"
    assert "ETKDG" in got["note"]
    # same molecule, different (regenerated) coordinates
    after = ligands.profile(got["text"].splitlines(keepends=True))
    assert after["heavy_atoms"] == before["heavy_atoms"]
    assert got["text"].splitlines()[4] != BUTANE.splitlines()[4]


# --- MMFF94 has no chirality restraint ----------------------------------------------

# A tetrahedral carbon with four different substituents, at the corners of a regular
# tetrahedron. `mirror=True` negates z, i.e. is the other enantiomer. (The real case is
# TRABECTEDIN, ID 1808 -- see the regression test at the end.)
def _stereo_molfile(mirror=False):
    s = -0.87 if mirror else 0.87
    return _molfile(
        "stereocentre",
        [("C", 0.0, 0.0, 0.0), ("C", 0.87, 0.87, s), ("N", 0.87, -0.87, -s),
         ("O", -0.87, 0.87, -s), ("F", -0.87, -0.87, s)],
        [(1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1)],
        [("ID", "STEREO1"), ("name", "A STEREOCENTRE")])


@pytest.mark.rdkit
def test_the_chirality_volume_has_a_sign_and_the_mirror_image_flips_it():
    from rdkit import Chem

    up = Chem.MolFromMolBlock(ligands._molblock(_lines(_stereo_molfile())),
                              removeHs=False)
    down = Chem.MolFromMolBlock(ligands._molblock(_lines(_stereo_molfile(mirror=True))),
                                removeHs=False)
    vol_up, vol_down = ligands.chirality_volumes(up), ligands.chirality_volumes(down)
    assert vol_up and set(vol_up) == set(vol_down)
    assert ligands.inverted_stereocentres(vol_up, vol_down) == sorted(vol_up)
    assert ligands.inverted_stereocentres(vol_up, vol_up) == []


@pytest.mark.rdkit
def test_a_minimisation_that_inverts_a_stereocentre_is_discarded(monkeypatch):
    """The guard itself: if MMFF94 flips a signed chirality volume, the minimised
    conformer is thrown away, the deposited geometry is kept, and the record says so."""
    from rdkit.Chem import AllChem

    def mirror(mol, **kwargs):      # stand in for a minimisation that inverts everything
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x, p.y, -p.z))
        return 0

    monkeypatch.setattr(AllChem, "MMFFOptimizeMolecule", mirror)
    lines = _lines(_stereo_molfile())
    got = ligands.prepare_record(lines)
    assert got["mmff_status"] == ligands.MMFF_STEREO_INVERSION
    assert got["mmff_converged"] == ""          # MMFF94 did not produce this geometry
    assert "inverted the stereocentre" in got["note"]
    assert "none" in got["method"] and "mmff94" not in got["method"]
    # the written coordinates are the deposited ones, not the mirrored ones
    from rdkit import Chem as _Chem
    kept = _Chem.MolFromMolBlock(got["text"].split("$$$$")[0], removeHs=False)
    source = _Chem.MolFromMolBlock(ligands._molblock(lines), removeHs=False)
    assert ligands.inverted_stereocentres(ligands.chirality_volumes(source),
                                          ligands.chirality_volumes(kept)) == []


@pytest.mark.rdkit
def test_an_ordinary_record_is_not_flagged_as_a_stereo_inversion(tmp_path):
    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE, _stereo_molfile()])
    counts = ligands.run(layout, _args(tmp_path), _Progress())
    assert counts["stereo_inversions"] == 0
    assert counts["minimised"] == 3


# --- implausible formal charges ------------------------------------------------------

# ISOSORBIDE DINITRATE (ID 271) in miniature: a charge line that declares more entries
# than it lists, so atom 1 -- a carbon -- reads as formal charge +15.
BAD_CHARGE = _molfile("bad charge", [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0),
                                     ("C", 3.0, 0.0, 0.0)],
                      [(1, 2, 1), (2, 3, 1)],
                      [("ID", "BADCHG1"), ("name", "CORRUPT CHARGE")],
                      charges=[(1, 15)])


@pytest.mark.rdkit
def test_an_implausible_formal_charge_is_flagged_counted_and_passed_through(tmp_path):
    """Dropping the record would silently change the library, so it is kept -- but it is
    never silent: counted, noted in ligands.csv and logged to failed.tsv."""
    layout = _layout(tmp_path, [BAD_CHARGE, BUTANE])
    counts = ligands.run(layout, _args(tmp_path), _Progress())
    assert counts["implausible_formal_charge"] == 1
    assert counts["kept"] == 2 and counts["failed"] == 0
    assert "BADCHG1" in (ligands.stage_dir(layout) / ligands.OUTPUT).read_text()
    log = (ligands.stage_dir(layout) / "failed.tsv").read_text()
    assert "BADCHG1" in log and "implausible formal charge" in log
    with open(ligands.stage_dir(layout) / ligands.TABLE, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert "implausible formal charge" in rows[0]["note"]
    assert rows[0]["net_charge"] == "15" and rows[1]["net_charge"] == "0"


@pytest.mark.rdkit
def test_a_normal_ionic_drug_is_not_flagged():
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(ligands._molblock(_lines(NITROMETHANE)), removeHs=False,
                               sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    ligands.charge_separate_hypervalent_nitrogen(mol)
    assert ligands.implausible_charges(mol) == []   # +1 / -1 is chemistry, not corruption


@pytest.mark.rdkit
def test_the_nitro_repair_refuses_to_judge_an_aromatic_nitrogen():
    """GetBondTypeAsDouble() scores an aromatic bond 1.5, so an aromatic ring N bearing an
    anionic O (the HOBt-anion motif) sums to 4.0 although its real valence is 3. The naive
    `> 3` gate set it to +1. e-Drug3D_2083 has no bond-type-4 bond, so this is a guard for
    the next library, not a fix for this one."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("c1ccc2c(c1)n(=O)nn2", sanitize=False)
    if mol is None:
        pytest.skip("rdkit would not parse the probe")
    mol.UpdatePropertyCache(strict=False)
    aromatic_n = [a for a in mol.GetAtoms()
                  if a.GetAtomicNum() == 7
                  and any(b.GetBondType() == Chem.BondType.AROMATIC for b in a.GetBonds())]
    assert aromatic_n, "the probe has no aromatic nitrogen"
    before = [a.GetFormalCharge() for a in aromatic_n]
    ligands.charge_separate_hypervalent_nitrogen(mol)
    assert [a.GetFormalCharge() for a in aromatic_n] == before


@pytest.mark.rdkit
def test_mmff_status_separates_no_parameters_from_not_converged(monkeypatch):
    from rdkit.Chem import AllChem

    monkeypatch.setattr(AllChem, "MMFFHasAllMoleculeParams", lambda mol: False)
    got = ligands.prepare_record(_lines(BUTANE))
    assert got["mmff_status"] == ligands.MMFF_NO_PARAMETERS
    assert got["mmff_converged"] == ""          # NOT False: MMFF94 never ran

    monkeypatch.undo()
    monkeypatch.setattr(AllChem, "MMFFOptimizeMolecule", lambda mol, **kw: 1)
    got = ligands.prepare_record(_lines(BUTANE))
    assert got["mmff_status"] == ligands.MMFF_NOT_CONVERGED
    assert got["mmff_converged"] is False


@pytest.mark.rdkit
def test_charge_separation_moves_no_atom():
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(ligands._molblock(_lines(NITROMETHANE)),
                               removeHs=False, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    before = [tuple(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    assert ligands.charge_separate_hypervalent_nitrogen(mol) == 1
    assert mol.GetAtomWithIdx(1).GetFormalCharge() == 1
    after = [tuple(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    assert before == after


# --- the stage ----------------------------------------------------------------------

def _args(run_dir, nprocs=1, library="inputs/lib.sdf"):
    return argparse.Namespace(run_dir=str(run_dir), nprocs=nprocs, library=library,
                              dry_run=False, no_notify=True)


def _layout(tmp_path, records):
    from polyphar.common import RunLayout

    layout = RunLayout(tmp_path)
    layout.inputs.mkdir(parents=True, exist_ok=True)
    (layout.inputs / "lib.sdf").write_text("".join(records))
    return layout


class _Progress:
    def __init__(self):
        self.n = 0

    def advance(self, n=1, size=None):
        self.n += n


BIG = _molfile("huge", [("C", 1.5 * i, 0.0, 0.0) for i in range(80)],
               [(i, i + 1, 1) for i in range(1, 80)],
               [("ID", "BIG1"), ("name", "A PEPTIDE")])   # 77 rotatable bonds, 961 Da

# RDKit cannot parse an unknown element symbol, but the pure-Python profile can (it
# counts it as carbon), so this record is dockable by the filter and unpreparable.
UNPARSABLE = _chain(["C", "C", "XX"], "broken",
                    [("ID", "BROKEN1"), ("name", "NOT A MOLECULE")])


@pytest.mark.rdkit
def test_run_filters_prepares_and_tabulates(tmp_path):
    layout = _layout(tmp_path, [BUTANE, BIG, CYCLOHEXANE])
    counts = ligands.run(layout, _args(tmp_path), _Progress())

    assert counts["records"] == 3 and counts["kept"] == 2 and counts["excluded"] == 1
    assert counts["excluded_rotatable_bonds"] == 1 and counts["excluded_mw"] == 1
    assert counts["failed"] == 0

    out = (ligands.stage_dir(layout) / ligands.OUTPUT).read_text()
    assert out.count("$$$$") == 2
    assert "BIG1" not in out and "BUT1" in out and "RING1" in out

    with open(ligands.stage_dir(layout) / ligands.TABLE, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert [r["ID"] for r in rows] == ["BUT1", "BIG1", "RING1"]   # library order
    assert [r["kept"] for r in rows] == ["True", "False", "True"]
    assert "rotatable_bonds=77>15" in rows[1]["exclusion"]


@pytest.mark.rdkit
def test_run_is_deterministic_and_nprocs_does_not_change_it(tmp_path):
    records = [BUTANE, CYCLOHEXANE, NITROMETHANE, BIG]
    one = tmp_path / "one"
    two = tmp_path / "two"
    a, b = _layout(one, records), _layout(two, records)
    ligands.run(a, _args(one, nprocs=1), _Progress())
    ligands.run(b, _args(two, nprocs=3), _Progress())
    first = (ligands.stage_dir(a) / ligands.OUTPUT).read_bytes()
    second = (ligands.stage_dir(b) / ligands.OUTPUT).read_bytes()
    assert first == second and first

    # and running it again in place reproduces the same bytes after the outputs are wiped
    (ligands.stage_dir(a) / "manifest.json").unlink()
    ligands.run(a, _args(one, nprocs=2), _Progress())
    assert (ligands.stage_dir(a) / ligands.OUTPUT).read_bytes() == first


@pytest.mark.rdkit
def test_resume_reports_zero_units_and_changes_nothing(tmp_path):
    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    assert ligands.units(layout, args).count == 2
    counts = ligands.run(layout, args, _Progress())

    out = ligands.stage_dir(layout) / ligands.OUTPUT
    before = (out.read_bytes(), out.stat().st_mtime_ns)
    assert ligands.units(layout, args).count == 0
    progress = _Progress()
    assert ligands.run(layout, args, progress) == counts
    assert progress.n == 0
    assert (out.read_bytes(), out.stat().st_mtime_ns) == before


# --- which library does `dock` actually read? (states a-e) ---------------------------
#
# The question "must the library be prepared again?" (ligands.fingerprint / units, which
# hashes the rdkit version) and "which file should dock read?" (ligands.prepared_state,
# which must not) used to share one answer. Where they disagree, dock silently docked the
# RAW library, which changes library_md5, re-docks every target and os.replace()s prepared
# results with unprepared ones. These tests pin dock._library_path in every state.

def _prepared(layout):
    return ligands.stage_dir(layout) / ligands.OUTPUT


@pytest.mark.rdkit
def test_state_a_no_ligands_stage_dock_reads_the_raw_library(tmp_path):
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    raw = layout.inputs / "lib.sdf"
    assert dock._library_path(layout, args, copy=False) == raw
    # nothing was prepared here, so there is nothing to warn about
    assert ligands.prepared_state(layout, raw) == (None, None)


@pytest.mark.rdkit
def test_state_b_a_finished_stage_makes_dock_read_the_prepared_library(tmp_path):
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    assert dock._library_path(layout, args, copy=False) == _prepared(layout)
    assert ligands.prepared_state(layout, layout.inputs / "lib.sdf")[1] is None
    assert dock._library_is_prepared(layout, args, _prepared(layout)) is True


@pytest.mark.rdkit
def test_state_c_a_deleted_prepared_file_falls_back_to_raw_and_warns(tmp_path, capsys):
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    _prepared(layout).unlink()                       # manifest stays, the SDF is gone
    raw = layout.inputs / "lib.sdf"
    assert dock._library_path(layout, args, copy=False) == raw
    assert "was deleted or moved" in ligands.prepared_state(layout, raw)[1]
    assert "WARNING" in capsys.readouterr().err
    assert dock._library_is_prepared(layout, args, raw) is False


@pytest.mark.rdkit
def test_state_d_a_changed_source_library_falls_back_to_raw_and_warns(tmp_path, capsys):
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    assert ligands.units(layout, args).count == 0
    raw = layout.inputs / "lib.sdf"
    raw.write_text(BUTANE + CYCLOHEXANE + NITROMETHANE)
    assert ligands.units(layout, args).count == 3    # it must be prepared again
    assert ligands.prepared_library(layout, raw) is None
    assert "changed after the prepared one" in ligands.prepared_state(layout, raw)[1]
    assert dock._library_path(layout, args, copy=False) == raw
    assert "WARNING" in capsys.readouterr().err


@pytest.mark.rdkit
def test_state_d2_a_tampered_prepared_file_falls_back_to_raw(tmp_path, capsys):
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    out = _prepared(layout)
    out.write_text(out.read_text() + "\n")           # no longer its recorded sha256
    assert dock._library_path(layout, args, copy=False) == layout.inputs / "lib.sdf"
    assert "changed after it was prepared" in ligands.prepared_state(
        layout, layout.inputs / "lib.sdf")[1]
    assert "WARNING" in capsys.readouterr().err


@pytest.mark.rdkit
def test_state_e_a_different_rdkit_in_dock_still_reads_the_prepared_library(
        tmp_path, monkeypatch):
    """The reproducible one: prepare with one interpreter, dock from another whose rdkit
    differs (or is missing). The files are intact, so dock MUST use the prepared library —
    anything else re-docks every target with an unprepared one."""
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    manifest = json.loads((ligands.stage_dir(layout) / "manifest.json").read_text())
    # the dock process reports a different rdkit than the one that prepared the library
    monkeypatch.setattr(ligands, "_rdkit_version", lambda: "a-different-version")

    raw = layout.inputs / "lib.sdf"
    assert ligands.fingerprint(raw) != manifest["fingerprint"]   # units() WOULD rerun
    assert ligands._finished(layout, raw) is None
    assert ligands.prepared_state(layout, raw) == (_prepared(layout), None)
    assert dock._library_path(layout, args, copy=False) == _prepared(layout)


def test_prepared_state_does_not_depend_on_rdkit_being_importable(tmp_path, monkeypatch):
    """Runs on an interpreter WITHOUT rdkit too: a hand-made manifest is enough, because
    prepared_state asks only about file identity."""
    from polyphar.common import sha256_file
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE])
    raw = layout.inputs / "lib.sdf"
    d = ligands.stage_dir(layout)
    d.mkdir(parents=True, exist_ok=True)
    out = d / ligands.OUTPUT
    out.write_text(BUTANE)
    (d / "manifest.json").write_text(json.dumps({
        "stage": "ligands", "counts": {"kept": 1},
        "fingerprint": "whatever-another-rdkit-produced",
        "library_sha256": sha256_file(out), "source_sha256": sha256_file(raw)}))
    monkeypatch.setattr(ligands, "_rdkit_version", lambda: "unavailable")
    assert ligands.prepared_state(layout, raw) == (out, None)
    assert dock._library_path(layout, _args(tmp_path), copy=False) == out
    assert dock._library_is_prepared(layout, _args(tmp_path), out) is True


def test_an_old_manifest_without_source_sha256_uses_the_recorded_input(tmp_path):
    from polyphar.common import sha256_file

    layout = _layout(tmp_path, [BUTANE])
    raw = layout.inputs / "lib.sdf"
    d = ligands.stage_dir(layout)
    d.mkdir(parents=True, exist_ok=True)
    out = d / ligands.OUTPUT
    out.write_text(BUTANE)
    (d / "manifest.json").write_text(json.dumps({
        "stage": "ligands", "library_sha256": sha256_file(out),
        "inputs": [{"path": str(raw), "sha256": sha256_file(raw)}]}))
    assert ligands.prepared_state(layout, raw) == (out, None)


# --- RDKit is a dependency, not 2,083 per-record failures ----------------------------

def test_the_stage_refuses_to_run_without_rdkit(tmp_path, monkeypatch):
    """Without rdkit the blanket per-record `except` used to classify every record as
    'could not prepare, source kept', write an UNPREPARED library and report success --
    and `dock` picked that up as prepared. `ligands` is not a FRACTION_STAGE, so
    --max-failed-fraction never caught it either."""
    monkeypatch.setattr(ligands, "rdkit_problem",
                        lambda: ligands.RDKIT_REQUIRED.format(why="ModuleNotFoundError"))
    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    for call in (lambda: ligands.units(layout, args), lambda: ligands.run(
            layout, args, _Progress())):
        with pytest.raises(RuntimeError, match="needs RDKit"):
            call()
    assert not (ligands.stage_dir(layout) / ligands.OUTPUT).exists()


def test_rdkit_problem_tells_the_truth_about_this_interpreter():
    problem = ligands.rdkit_problem()
    try:
        import rdkit  # noqa: F401
    except ImportError:
        assert problem and "needs RDKit" in problem
    else:
        assert problem is None


@pytest.mark.rdkit
def test_a_finished_stage_does_not_need_rdkit_to_report_zero_units(tmp_path, monkeypatch):
    """Resuming a run folder whose ligands stage is done must not demand the dependency."""
    layout = _layout(tmp_path, [BUTANE])
    args = _args(tmp_path)
    ligands.run(layout, args, _Progress())
    monkeypatch.setattr(ligands, "rdkit_problem", lambda: "no rdkit here")
    assert ligands.units(layout, args).count == 0
    assert ligands.run(layout, args, _Progress())["kept"] == 1


def test_preflight_requires_the_library_and_rdkit_for_ligands(tmp_path, monkeypatch):
    from polyphar.common import RunLayout
    from polyphar.stages import run as run_stage

    layout = RunLayout(tmp_path)
    layout.inputs.mkdir(parents=True, exist_ok=True)
    args = _args(tmp_path)
    # the library check used to run only for `dock`
    problems = run_stage.preflight_problems(layout, args, ["ligands"])
    assert any("not found" in p for p in problems), problems
    (layout.inputs / "lib.sdf").write_text(BUTANE)
    monkeypatch.setattr(ligands, "rdkit_problem", lambda: "rdkit is missing")
    assert run_stage.preflight_problems(layout, args, ["ligands"]) == [
        "ligands: rdkit is missing"]
    monkeypatch.setattr(ligands, "rdkit_problem", lambda: None)
    assert run_stage.preflight_problems(layout, args, ["ligands"]) == []


@pytest.mark.rdkit
def test_an_unpreparable_record_is_kept_verbatim_and_logged(tmp_path):
    """Dropping it would silently change the library, so it is written through unchanged
    and counted as failed."""
    layout = _layout(tmp_path, [UNPARSABLE, BUTANE])
    counts = ligands.run(layout, _args(tmp_path), _Progress())
    assert counts["failed"] == 1 and counts["kept"] == 2
    out = (ligands.stage_dir(layout) / ligands.OUTPUT).read_text()
    assert out.startswith(UNPARSABLE)          # verbatim, first, still in library order
    log = (ligands.stage_dir(layout) / "failed.tsv").read_text()
    assert "BROKEN1" in log and "cannot parse" in log


@pytest.mark.rdkit
def test_worker_never_raises_for_an_unpreparable_record():
    got = ligands._worker(0, _lines(UNPARSABLE))
    assert got["failed"] and got["text"] == UNPARSABLE
    assert got["row"]["method"] == "source"


@pytest.mark.rdkit
def test_sanitisation_fallbacks_are_logged_to_failed_tsv(tmp_path):
    layout = _layout(tmp_path, [NITROMETHANE])
    counts = ligands.run(layout, _args(tmp_path), _Progress())
    assert counts["charge_separated"] == 1 and counts["failed"] == 0
    log = (ligands.stage_dir(layout) / "failed.tsv").read_text()
    assert "NITRO1" in log and "charge-separated" in log


@pytest.mark.rdkit
def test_manifest_records_the_hash_dock_fingerprints(tmp_path):
    from polyphar.common import sha256_file

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    ligands.run(layout, _args(tmp_path), _Progress())
    manifest = json.loads((ligands.stage_dir(layout) / "manifest.json").read_text())
    out = ligands.stage_dir(layout) / ligands.OUTPUT
    assert manifest["library_sha256"] == sha256_file(out)
    assert manifest["fingerprint"] == ligands.fingerprint(layout.inputs / "lib.sdf")
    assert manifest["versions"]["max_rotatable_bonds"] == 15
    assert "NOT RDKit/Lipinski" in manifest["versions"]["rotatable_bond_definition"]


@pytest.mark.rdkit
def test_dock_picks_the_prepared_library_up(tmp_path):
    """The only coupling to dock: once this stage has run on the source library, dock's
    library IS the prepared file, so dock's library_md5 fingerprint covers preparation."""
    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE, CYCLOHEXANE])
    args = _args(tmp_path)
    assert dock._library_path(layout, args, copy=False) == layout.inputs / "lib.sdf"
    ligands.run(layout, args, _Progress())
    got = dock._library_path(layout, args, copy=False)
    assert got == ligands.stage_dir(layout) / ligands.OUTPUT
    assert dock._raw_library_path(layout, args, copy=False) == layout.inputs / "lib.sdf"
    # and the prepared library is still a library dock can read
    assert list(dock.read_library(got)["ID"]) == ["BUT1", "RING1"]


def test_prepared_library_is_ignored_when_the_stage_has_not_run(tmp_path):
    layout = _layout(tmp_path, [BUTANE])
    assert ligands.prepared_library(layout, layout.inputs / "lib.sdf") is None


# --- pipeline wiring ----------------------------------------------------------------

def test_ligands_is_in_the_pipeline_order_before_dock():
    from polyphar.benchmark import ORDER, STAGE_DIRS

    assert "ligands" in ORDER
    assert ORDER.index("ligands") < ORDER.index("dock")
    assert STAGE_DIRS["ligands"] == "07_docking/library"


def test_the_prepared_library_survives_docks_orphan_sweep(tmp_path):
    """Behavioural, not a membership assertion: 07_docking/<label> dirs that are not
    current targets are MOVED to obsolete/, and 07_docking/library/ must not be among
    them. The counter-factual at the end shows the guard is what prevents it."""
    import types

    from polyphar.stages import dock

    layout = _layout(tmp_path, [BUTANE])
    d = layout.docking
    for sub in ("library", "ligands", "obsolete", "T_P_1", "GONE_P_9"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    (d / "library" / ligands.OUTPUT).write_text("PREPARED\n")
    (d / "library" / "manifest.json").write_text("{}")
    (d / "ligands" / "c50_00000.sdf").write_text("chunk\n")
    (d / "df_results_T_P_1.csv").write_text("a\n")
    (d / "df_results_GONE_P_9.csv").write_text("b\n")
    # the retired replicate-0 schema: an orphan under EVERY label, current or not
    (d / "df_results_T_P_1_fixed.csv").write_text("legacy\n")
    plan = dock._Plan(targets=[types.SimpleNamespace(label="T_P_1")], library_sha="",
                      n_ligands=0, usable=[], missing=[], units=[], opts={})

    orphans = dock._orphans(layout, plan)
    assert sorted(p.name for p in orphans) == [
        "GONE_P_9", "df_results_GONE_P_9.csv", "df_results_T_P_1_fixed.csv"]
    dock._obsolete(layout, orphans)
    assert (d / "library" / ligands.OUTPUT).read_text() == "PREPARED\n"
    assert (d / "library" / "manifest.json").is_file()
    assert (d / "ligands" / "c50_00000.sdf").is_file()
    assert not (d / "GONE_P_9").exists()

    # counter-factual: without the RESERVED_DIRS entry the prepared library IS an orphan
    without = tuple(x for x in dock.RESERVED_DIRS if x != ligands.DIRNAME)
    saved, dock.RESERVED_DIRS = dock.RESERVED_DIRS, without
    try:
        assert ligands.DIRNAME in {p.name for p in dock._orphans(layout, plan)}
    finally:
        dock.RESERVED_DIRS = saved


def test_cli_wires_the_subcommand():
    from polyphar.cli import STAGES, build_parser

    assert STAGES["ligands"][0] == "polyphar.stages.ligands"
    args = build_parser().parse_args(["ligands", "--run-dir", "/tmp/x"])
    assert args.command == "ligands" and args.library == "inputs/e-Drug3D_2083.sdf"


def test_library_option_matches_dock_so_run_can_merge_them():
    from polyphar.stages import dock

    a, b = argparse.ArgumentParser(), argparse.ArgumentParser()
    ligands.add_arguments(a)
    dock.add_arguments(b)
    mine = [x for x in a._actions if "--library" in x.option_strings][0]
    theirs = [x for x in b._actions if "--library" in x.option_strings][0]
    assert (mine.default, mine.type, mine.nargs) == (theirs.default, theirs.type,
                                                     theirs.nargs)


def test_projected_units_knows_the_stage():
    from polyphar.benchmark import PRIORS, projected_units

    units = projected_units("ligands", dict(PRIORS))
    assert units is not None and units.cost_basis == "ligands"


# --- the real library ---------------------------------------------------------------

@pytest.mark.regression
def test_real_edrug3d_split_is_123_excluded_1960_kept():
    """The numbers the design document reports, against the real library (ADEP_DATA)."""
    import os
    from pathlib import Path

    from polyphar.stages.dock import _library_records

    sdf = Path(os.environ["ADEP_DATA"]) / "docking_screen_data" / "libraries" \
        / "e-Drug3D_2083.sdf"
    if not sdf.is_file():
        pytest.skip(f"{sdf} not present")
    records = _library_records(sdf)
    assert len(records) == 2083
    excluded = [r for r in records if ligands.exclusion(ligands.profile(r.lines))]
    assert len(excluded) == 123
    assert len(records) - len(excluded) == 1960


@pytest.mark.regression
@pytest.mark.rdkit
def test_real_trabectedin_is_caught_and_keeps_its_deposited_geometry():
    """TRABECTEDIN (ID 1808, record 1807) is the one record of e-Drug3D_2083 where MMFF94
    inverts a stereocentre (atom 18, S -> R). Without the guard an approved antitumour
    agent would be docked as the wrong diastereomer."""
    import os
    from pathlib import Path

    from polyphar.stages.dock import _library_records

    sdf = Path(os.environ["ADEP_DATA"]) / "docking_screen_data" / "libraries" \
        / "e-Drug3D_2083.sdf"
    if not sdf.is_file():
        pytest.skip(f"{sdf} not present")
    record = next(r for r in _library_records(sdf) if r.fields.get("ID") == "1808")
    assert "TRABECTEDIN" in record.fields.get("name", "").upper()
    got = ligands.prepare_record(list(record.lines))
    assert got["mmff_status"] == ligands.MMFF_STEREO_INVERSION
    assert "18" in got["note"] and "deposited geometry kept" in got["note"]

    from rdkit import Chem

    before = Chem.MolFromMolBlock(ligands._molblock(list(record.lines)), removeHs=False)
    after = Chem.MolFromMolBlock(got["text"].split("$$$$")[0], removeHs=False)
    assert ligands.inverted_stereocentres(ligands.chirality_volumes(before),
                                          ligands.chirality_volumes(after)) == []


@pytest.mark.regression
@pytest.mark.rdkit
def test_real_isosorbide_dinitrate_is_flagged_not_dropped():
    """ID 271: `M  CHG  5   1  15  13   1  14   1  17  -1  18  -1` declares five entries
    and lists four, so atom 1 (a carbon) reads as +15. RDKit sanitises it and MMFF94
    parameterises it, so nothing downstream would notice."""
    import os
    from pathlib import Path

    from polyphar.stages.dock import _library_records

    sdf = Path(os.environ["ADEP_DATA"]) / "docking_screen_data" / "libraries" \
        / "e-Drug3D_2083.sdf"
    if not sdf.is_file():
        pytest.skip(f"{sdf} not present")
    record = next(r for r in _library_records(sdf) if r.fields.get("ID") == "271")
    assert "ISOSORBIDE" in record.fields.get("name", "").upper()
    got = ligands.prepare_record(list(record.lines))
    assert got["bad_charges"] and got["bad_charges"][0][1] == "C"
    assert got["bad_charges"][0][2] == 15
    assert "implausible formal charge" in got["note"]
    assert got["text"].count("$$$$") == 1       # passed through, not dropped


@pytest.mark.repo_files
def test_the_docker_wrapper_knows_every_stage():
    """bin/polyphar-docker's stage_index() knew only the 2023 eight stages. An unknown
    name produced an empty string, `${from:-0}` / `${to:-N}` substituted the full range,
    and `--from ligands` therefore demanded $DOGSITE_BIN for a stage that touches neither
    DoGSiteScorer nor pockets."""
    import re
    import subprocess
    from pathlib import Path

    from polyphar.benchmark import ORDER

    wrapper = Path(__file__).resolve().parents[1] / "bin" / "polyphar-docker"
    body = re.search(r"stage_index\(\) \{.*?\n\}", wrapper.read_text(), re.S).group(0)
    for i, stage in enumerate(ORDER):
        out = subprocess.run(["bash", "-c", f"{body}\nstage_index {stage}"],
                             capture_output=True, text=True, check=True)
        assert out.stdout.strip() == str(i), f"{stage}: {out.stdout!r}, expected {i}"
    # the pockets bound the wrapper tests against, and the last index it defaults --to to
    text = wrapper.read_text()
    assert f"pockets_index={ORDER.index('pockets')}" in text
    assert f"last_index={len(ORDER) - 1}" in text


def test_stage_module_exposes_the_stage_contract():
    for name in ("units", "run", "add_arguments"):
        assert callable(getattr(ligands, name))
    assert isinstance(ligands.STAGE, str)
    assert isinstance(ligands.units, types.FunctionType)


# --- elements the scoring function cannot type ----------------------------------------

def test_profile_reports_the_elements_present():
    record = _molfile("X", [("C", 0, 0, 0), ("B", 1.5, 0, 0), ("O", 3.0, 0, 0)],
                      [(1, 2, 1), (2, 3, 1)])
    assert ligands.profile(_lines(record))["elements"] == ["B", "C", "O"]


def test_a_boron_ligand_is_excluded_because_vina_cannot_type_it():
    """Uni-Dock answers a boron ligand by writing no pose at all, with no error. Six such
    drugs (bortezomib, ixazomib, tavaborole, ...) were missing from every replicate on
    net1, so no target could ever be complete."""
    prof = {"heavy_atoms": 3, "mw": 100.0, "rotatable_bonds": 0, "elements": ["B", "C", "O"]}
    why = ligands.exclusion(prof)
    assert "cannot type" in why and "B" in why


def test_the_ordinary_drug_elements_are_all_dockable():
    prof = {"heavy_atoms": 20, "mw": 300.0, "rotatable_bonds": 3,
            "elements": sorted(ligands.DOCKABLE_ELEMENTS)}
    assert ligands.exclusion(prof) == ""


def test_a_record_with_no_element_information_is_not_excluded_for_elements():
    """profile() defaults a blank element column to carbon; exclusion() must not invent a
    reason when `elements` is absent entirely (older callers, fixtures)."""
    assert ligands.exclusion({"heavy_atoms": 5, "mw": 80.0, "rotatable_bonds": 0}) == ""
