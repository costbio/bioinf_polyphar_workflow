"""Detecting main chain that was traced but never sequenced."""

from polyphar import modelling


def _pdb(path, residues):
    """residues = [(resname, [atom names])]."""
    lines, i = [], 0
    for r, (name, atoms) in enumerate(residues, start=1):
        for atom in atoms:
            i += 1
            el = atom[0]
            lines.append(f"ATOM  {i:5d}  {atom:<3} {name} A{r:4d}    "
                         f"{r:8.3f}{i:8.3f}{r:8.3f}  1.00  0.00           {el}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n")
    return path


FULL = ("LEU", ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"])
STUB = ("LEU", ["N", "CA", "C", "O", "CB"])
UNK = ("UNK", ["N", "CA", "C", "O", "CB"])
GLY = ("GLY", ["N", "CA", "C", "O"])
ALA = ("ALA", ["N", "CA", "C", "O", "CB"])


def test_a_residue_with_its_side_chain_is_modelled():
    assert not modelling.is_unmodelled(*(FULL[0], FULL[1]))


def test_a_backbone_stub_is_unmodelled():
    assert modelling.is_unmodelled(*(STUB[0], STUB[1]))


def test_unk_is_unmodelled_whatever_atoms_it_carries():
    """Its identity is unknown, so it has no chemistry even with a full atom set."""
    assert modelling.is_unmodelled("UNK", ["N", "CA", "C", "O", "CB", "CG", "CD1"])


def test_glycine_and_alanine_are_not_unmodelled():
    """They have nothing beyond CB to model; calling them stubs would reject real sites."""
    assert not modelling.is_unmodelled(*(GLY[0], GLY[1]))
    assert not modelling.is_unmodelled(*(ALA[0], ALA[1]))


def test_unmodelled_fraction_counts_residues_not_atoms(tmp_path):
    p = _pdb(tmp_path / "m.pdb", [FULL, FULL, UNK, STUB])
    frac, bad, total = modelling.unmodelled_fraction(p)
    assert (bad, total) == (2, 4)
    assert frac == 0.5


def test_a_fully_modelled_file_is_zero(tmp_path):
    p = _pdb(tmp_path / "m.pdb", [FULL, GLY, ALA])
    assert modelling.unmodelled_fraction(p) == (0.0, 0, 3)


def test_an_all_unk_pocket_is_entirely_unmodelled(tmp_path):
    """The shape that hung PocketMatch: P42356/6bq1_A_res_P_10.pdb, 26 UNK residues of
    N/CA/C/O/CB."""
    p = _pdb(tmp_path / "pocket.pdb", [UNK] * 26)
    frac, bad, total = modelling.unmodelled_fraction(p)
    assert (frac, bad, total) == (1.0, 26, 26)


def test_an_empty_file_is_not_a_division_by_zero(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("END\n")
    assert modelling.unmodelled_fraction(p) == (0.0, 0, 0)


def test_residues_are_keyed_by_chain_number_and_insertion_code(tmp_path):
    """Two residues sharing a number but differing by insertion code are two residues."""
    p = tmp_path / "i.pdb"
    p.write_text(
        "ATOM      1  N   LEU A  10     1.000   1.000   1.000  1.00  0.00           N\n"
        "ATOM      2  N   LEU A  10A    2.000   2.000   2.000  1.00  0.00           N\n"
        "END\n")
    assert len(modelling.residue_atoms(p)) == 2


def test_a_truncated_record_does_not_raise(tmp_path):
    """Real files carry short ATOM lines (hand-edited, or a fixture); reading one must
    not take the stage down with an IndexError."""
    p = tmp_path / "short.pdb"
    p.write_text("ATOM      1  N   LEU\nATOM      2  CA  LEU A   1\nEND\n")
    frac, bad, total = modelling.unmodelled_fraction(p)
    assert total >= 1 and 0.0 <= frac <= 1.0
