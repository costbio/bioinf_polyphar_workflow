"""Coordinate geometry shared by `select` and `dock`: the docking box, and the check
that a receptor actually occupies it.

A pocket and its receptor must be the same structure. When they are not, the box is
computed in one structure's frame and the receptor sits in another's, so the box can
land in empty space and every ligand is docked into vacuum -- scored, written and
indistinguishable from a real result. `select` pairs them and `dock` uses the pair, so
both call `check_box_occupancy` and both refuse the target.
"""

from __future__ import annotations

from pathlib import Path

# A box holding fewer receptor atoms than this is not a pocket. A real pocket is lined by
# its protein: the smallest occupancy seen across net1's correctly paired targets was 86
# atoms, while every mispaired one held exactly 0. The threshold is deliberately far below
# that: this guard exists to catch a box in empty space, not to judge pocket quality.
MIN_BOX_ATOMS = 10


def read_atoms(path) -> list[tuple[float, float, float]]:
    """ATOM/HETATM coordinates of a PDB or PDBQT file, in file order."""
    xyz = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return xyz


def bounding_box(xyz) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """(centre, size) of a coordinate list, unpadded."""
    if not xyz:
        raise ValueError("no atoms")
    lo = [min(c[k] for c in xyz) for k in range(3)]
    hi = [max(c[k] for c in xyz) for k in range(3)]
    return (tuple((hi[k] + lo[k]) / 2 for k in range(3)),
            tuple(hi[k] - lo[k] for k in range(3)))


def atoms_in_box(xyz, center, size) -> int:
    """How many of `xyz` fall inside the axis-aligned box."""
    lo = [center[k] - size[k] / 2 for k in range(3)]
    hi = [center[k] + size[k] / 2 for k in range(3)]
    return sum(1 for c in xyz
               if all(lo[k] <= c[k] <= hi[k] for k in range(3)))


def check_box_occupancy(receptor, center, size, *, label="", minimum=MIN_BOX_ATOMS) -> int:
    """Receptor atoms inside the box. Raises ValueError when the box is (near) empty.

    The message names both files, because the cause is almost always that they are
    different structures."""
    xyz = read_atoms(receptor)
    if not xyz:
        raise ValueError(f"{Path(receptor).name}: no ATOM/HETATM records")
    n = atoms_in_box(xyz, center, size)
    if n < minimum:
        cen, _ = bounding_box(xyz)
        gap = sum((center[k] - cen[k]) ** 2 for k in range(3)) ** 0.5
        raise ValueError(
            f"{label + ': ' if label else ''}the docking box holds {n} receptor atom(s) "
            f"(minimum {minimum}). The box centre is {gap:.0f} A from the centre of "
            f"{Path(receptor).name}, so the pocket and the receptor are not the same "
            f"structure and every ligand would be docked into empty space.")
    return n
