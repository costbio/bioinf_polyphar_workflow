"""The three defects that reached live data in November 2026.

Each was trivial once seen, and each cost real work: receptors paired with another
structure's pocket (9 of 37 net1 targets docked into empty space), six boron drugs the
engine silently never poses (no target could ever complete), and a pocket of UNK residues
that made PocketMatch spin for 2.5 hours.

They are pinned here against the MECHANISM rather than the instance, so a different
structure, a different element or a different unmodelled residue is caught too. The
mechanism each shares: an upstream tool answered in a way the code read as success.
"""

import pytest

from polyphar import geometry, modelling
from polyphar.stages import ligands, similarity


def _pdb(path, residues, origin=(0.0, 0.0, 0.0)):
    """residues = [(resname, [atom names])], laid out from `origin`."""
    lines, i = [], 0
    ox, oy, oz = origin
    for r, (name, atoms) in enumerate(residues, start=1):
        for atom in atoms:
            i += 1
            lines.append(f"ATOM  {i:5d}  {atom:<3} {name} A{r:4d}    "
                         f"{ox + r:8.3f}{oy + i * 0.5:8.3f}{oz + r:8.3f}"
                         f"  1.00  0.00           {atom[0]}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n")
    return path


LEU = ("LEU", ["N", "CA", "C", "O", "CB", "CG", "CD1"])
UNK = ("UNK", ["N", "CA", "C", "O", "CB"])


# --- 1. a pocket from another structure ------------------------------------------------

def test_a_box_the_receptor_does_not_occupy_is_refused(tmp_path):
    """83 of net1's 479 target rows paired a pocket with a receptor from a different
    structure. Where the frames differed the box landed off the protein and every ligand
    was scored in vacuum -- three of those targets produced a FULL set of poses and would
    have entered the analysis as real affinities."""
    receptor = _pdb(tmp_path / "r.pdb", [LEU] * 30)
    stranger = _pdb(tmp_path / "other.pdb", [LEU] * 8, origin=(500.0, 0.0, 0.0))
    center, size = geometry.bounding_box(geometry.read_atoms(stranger))
    with pytest.raises(ValueError) as exc:
        geometry.check_box_occupancy(receptor, center, size, label="X_A_P_1")
    assert "not the same structure" in str(exc.value)


def test_a_pocket_of_its_own_receptor_is_accepted(tmp_path):
    """The other half: the guard must not reject real work."""
    receptor = _pdb(tmp_path / "r.pdb", [LEU] * 30)
    atoms = geometry.read_atoms(receptor)
    center, size = geometry.bounding_box(atoms[:60])
    assert geometry.check_box_occupancy(receptor, center, size) >= geometry.MIN_BOX_ATOMS


# --- 2. a ligand the engine silently never poses ---------------------------------------

def test_an_element_the_scoring_function_cannot_type_is_excluded_at_the_library():
    """Uni-Dock answers a boron ligand by writing no pose, with no error in the log. Six
    such drugs were absent from all 102 net1 replicates, and since a target is complete
    only when every ligand has poses, NO target could ever complete. The library states
    the limitation once instead of letting it surface as unexplained incompleteness."""
    prof = {"heavy_atoms": 28, "mw": 384.2, "rotatable_bonds": 11,
            "elements": ["B", "C", "N", "O"]}
    why = ligands.exclusion(prof)
    assert "cannot type" in why and "B" in why


def test_the_rule_is_the_element_set_not_a_list_of_known_drugs():
    """Pinned against the mechanism: any element outside the scoring function's types is
    excluded, not just the six boron drugs that happened to be in e-Drug3D."""
    assert "B" not in ligands.DOCKABLE_ELEMENTS
    for element in ("SI", "SE", "PT", "AS"):
        why = ligands.exclusion({"heavy_atoms": 10, "mw": 200.0, "rotatable_bonds": 1,
                                 "elements": ["C", element]})
        assert "cannot type" in why, element
    ordinary = ligands.exclusion({"heavy_atoms": 20, "mw": 300.0, "rotatable_bonds": 3,
                                  "elements": ["C", "H", "N", "O", "S", "CL", "F"]})
    assert ordinary == ""


# --- 3. an input PocketMatch cannot type -----------------------------------------------

def test_a_pocket_with_one_unmodelled_residue_never_reaches_pocketmatch(tmp_path):
    """Step3-PM_serial does not fail on an input it cannot handle, it SPINS: 100% CPU,
    ~1 MB RSS, no output, no error. One UNK residue in an otherwise healthy pocket is
    enough, and six such pockets stalled 11 of Cluster1's 66 jobs for 2.5 hours."""
    good = _pdb(tmp_path / "good.pdb", [LEU] * 8)
    one_unk = _pdb(tmp_path / "one_unk.pdb", [LEU] * 7 + [UNK])
    found = {p.name for p in similarity.untypable_pockets([good, one_unk])}
    assert found == {"one_unk.pdb"}


def test_a_job_that_will_not_finish_is_killed_rather_than_waited_on(tmp_path):
    """The general guard behind the specific one: PocketMatch spins rather than failing,
    so every job is bounded. Without this the stage waits for ever with no diagnostic."""
    assert similarity.PM_JOB_TIMEOUT > 0
    fake = tmp_path / "bin"
    fake.mkdir()
    for name in (similarity.CABBAGE, similarity.END_FILE):
        (fake / name).write_text("#!/bin/sh\nexit 0\n")
    (fake / similarity.PM_SERIAL).write_text("#!/bin/sh\nsleep 30\n")
    for name in similarity.PM_BINARIES:
        (fake / name).chmod(0o755)
    src = _pdb(tmp_path / "p.pdb", [LEU] * 4)
    with pytest.raises(similarity.PocketMatchError, match="did not finish"):
        similarity.run_pocketmatch([src], fake, tmp_path / "work", records=[src], timeout=1)


# --- the shared shape ------------------------------------------------------------------

def test_an_unmodelled_pocket_is_dropped_before_it_can_be_docked(tmp_path):
    """A pocket lined by traced main chain has no chemistry: docking it scores a backbone
    groove. Same root as the PocketMatch hang, one stage earlier."""
    pocket = _pdb(tmp_path / "p.pdb", [UNK] * 10)
    frac, bad, total = modelling.unmodelled_fraction(pocket)
    assert (frac, bad, total) == (1.0, 10, 10)
