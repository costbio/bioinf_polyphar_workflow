"""The box-occupancy guard: a pocket and a receptor from different structures must never
reach the docking engine."""

import pytest

from polyphar import geometry


def _pdb(tmp_path, name, coords):
    lines = []
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


def _cube(origin, step=1.0, n=4):
    ox, oy, oz = origin
    return [(ox + i * step, oy + j * step, oz + k * step)
            for i in range(n) for j in range(n) for k in range(n)]


def test_bounding_box_centre_and_size():
    center, size = geometry.bounding_box([(0, 0, 0), (2, 4, 6)])
    assert center == (1, 2, 3)
    assert size == (2, 4, 6)


def test_atoms_in_box_counts_only_what_is_inside():
    xyz = [(0, 0, 0), (10, 10, 10)]
    assert geometry.atoms_in_box(xyz, (0, 0, 0), (2, 2, 2)) == 1


def test_a_receptor_that_occupies_its_box_passes(tmp_path):
    receptor = _pdb(tmp_path, "r.pdb", _cube((0, 0, 0)))
    n = geometry.check_box_occupancy(receptor, (1.5, 1.5, 1.5), (4, 4, 4))
    assert n == 64


def test_a_box_in_empty_space_is_refused_and_says_how_far(tmp_path):
    """The real failure: the pocket came from one structure and the receptor from
    another, so the box sits off the protein and every ligand would dock in vacuum."""
    receptor = _pdb(tmp_path, "r.pdb", _cube((0, 0, 0)))
    with pytest.raises(ValueError) as exc:
        geometry.check_box_occupancy(receptor, (100.0, 0.0, 0.0), (10, 10, 10),
                                     label="2MJW_B_P_1")
    msg = str(exc.value)
    assert "2MJW_B_P_1" in msg and "0 receptor atom" in msg
    assert "not the same structure" in msg
    assert "98 A" in msg or "99 A" in msg       # 100 - the cube's own centre at 1.5


def test_a_box_holding_a_handful_of_atoms_is_refused(tmp_path):
    """Grazing the edge of the protein is still not a pocket."""
    receptor = _pdb(tmp_path, "r.pdb", _cube((0, 0, 0)))
    with pytest.raises(ValueError, match="minimum"):
        geometry.check_box_occupancy(receptor, (-2.0, -2.0, -2.0), (2.2, 2.2, 2.2))


def test_an_empty_receptor_file_is_refused(tmp_path):
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    with pytest.raises(ValueError, match="no ATOM/HETATM"):
        geometry.check_box_occupancy(empty, (0, 0, 0), (10, 10, 10))
