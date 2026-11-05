"""How much of a structure is actually modelled.

A deposited chain can carry stretches of main chain with no side chains and no identity:
`UNK` residues, or poly-ALA traces, where the experiment resolved a backbone path but not
the sequence. They are real coordinates, so nothing upstream rejects them, but they carry
no chemistry:

- PocketMatch types every atom by its residue, and `Step3-PM_serial` SPINS FOREVER on a
  pocket holding one (measured 2026-11-05: a single UNK residue is enough).
- Docking against them is meaningless: the side chains that would line the site are not
  there, so any affinity is an artefact of a backbone groove.

So a pocket lined by unmodelled residues is dropped where it is produced, and a chain that
is unmodelled throughout is rejected at extraction.
"""

from __future__ import annotations

# Atoms every residue has (plus CB, which GLY lacks and ALA ends at). A residue carrying
# nothing beyond these has no side chain modelled.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT", "CB"})
# GLY and ALA genuinely have no side chain past CB, so their absence is not a gap.
NO_SIDECHAIN_RESIDUES = frozenset({"GLY", "ALA"})
UNKNOWN_RESIDUE = "UNK"


def residue_atoms(path) -> dict:
    """{(chain, resseq, icode): (residue name, {atom names})} of a PDB file."""
    out: dict = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            # a record can be truncated (a fixture, a hand-edited file): pad rather than
            # raise IndexError on the chain or insertion-code column
            line = line.ljust(27)
            key = (line[21], line[22:26].strip(), line[26])
            name = line[17:20].strip().upper()
            if key not in out:
                out[key] = (name, set())
            out[key][1].add(line[12:16].strip().upper())
    return out


def is_unmodelled(name: str, atoms) -> bool:
    """True when this residue carries no modelled side chain.

    `UNK` always counts: its identity is unknown whatever atoms it has. Otherwise a
    residue counts when nothing beyond the backbone (and CB) is present, except GLY and
    ALA, which have nothing more to model."""
    if name == UNKNOWN_RESIDUE:
        return True
    if name in NO_SIDECHAIN_RESIDUES:
        return False
    return not (set(atoms) - BACKBONE_ATOMS)


def unmodelled_fraction(path) -> tuple[float, int, int]:
    """(fraction, unmodelled residues, residues) of a PDB file. An empty file is (0, 0, 0)."""
    res = residue_atoms(path)
    if not res:
        return 0.0, 0, 0
    bad = sum(1 for name, atoms in res.values() if is_unmodelled(name, atoms))
    return bad / len(res), bad, len(res)
