"""Merge-path hygiene: how a merge copies, and how it remembers what it dropped.

Both matter more now that `polyphar dock` merges a bundle that is still filling up.
Neither loses a draw; both quietly fill the discard archive with debris, which is worse
than it sounds — a discard is supposed to mean "a real sample was thrown away", and a log
full of fragments makes the real ones unreadable.
"""

import os
from pathlib import Path

import pytest

from polyphar.stages import dock


def test_a_copy_lands_whole_or_not_at_all(tmp_path):
    """A plain copyfile leaves a partial file when the merge is interrupted. The next
    merge then sees something standing in the unit's place and archives it as a discard."""
    src = tmp_path / "src.sdf"
    src.write_text("record\n$$$$\n")
    dst = tmp_path / "dst.sdf"
    dock._copy_atomic(src, dst)
    assert dst.read_text() == "record\n$$$$\n"
    assert not list(tmp_path.glob(".*tmp")), "a temporary file was left behind"


def test_a_failed_copy_leaves_nothing_in_place(tmp_path, monkeypatch):
    src = tmp_path / "src.sdf"
    src.write_text("record\n$$$$\n")
    dst = tmp_path / "dst.sdf"

    import shutil

    def boom(*_a, **_k):
        raise OSError("disk full")

    monkeypatch.setattr(shutil, "copyfile", boom)
    with pytest.raises(OSError):
        dock._copy_atomic(src, dst)
    assert not dst.exists(), "a partial file stands in the unit's place"
    assert not list(tmp_path.glob(".*tmp")), "a temporary file was left behind"


def test_a_copy_over_an_existing_file_replaces_it_wholly(tmp_path):
    src = tmp_path / "src.sdf"
    src.write_text("new\n$$$$\n")
    dst = tmp_path / "dst.sdf"
    dst.write_text("old and much longer than the new one\n$$$$\n")
    dock._copy_atomic(src, dst)
    assert dst.read_text() == "new\n$$$$\n"


# --- remembering what was dropped -----------------------------------------------------

def _unit(label="A", replicate=0, key="k0"):
    return dock.Unit(label=label, replicate=replicate, seed=42, key=key)


def test_logged_drops_reads_the_key_and_the_source_as_fields(tmp_path):
    """The old test was two substring scans over the whole log per dropped unit: O(log
    size) each, and loose enough that two units of one bundle could satisfy it between
    them without either having been logged."""
    dock.log_discard(tmp_path, _unit(key="key-aaa"), "dropped on merge", "/b/A/r0.sdf")
    got = dock._logged_drops(tmp_path)
    assert got == {("key-aaa", "/b/A/r0.sdf")}


def test_a_drop_is_matched_on_the_pair_not_on_either_half(tmp_path):
    """Same key from another bundle, or another key from the same file, is a different
    drop and must still be logged."""
    dock.log_discard(tmp_path, _unit(key="key-aaa"), "dropped", "/bundle1/A/r0.sdf")
    logged = dock._logged_drops(tmp_path)
    assert ("key-aaa", "/bundle2/A/r0.sdf") not in logged
    assert ("key-bbb", "/bundle1/A/r0.sdf") not in logged
    assert ("key-aaa", "/bundle1/A/r0.sdf") in logged


def test_the_header_is_not_read_as_a_drop(tmp_path):
    dock._ensure_discard_log(tmp_path)
    assert dock._logged_drops(tmp_path) == set()


def test_a_reason_containing_tabs_cannot_shift_the_columns(tmp_path):
    """log_discard cleans the free-text fields; if it did not, a tab in a reason would
    move the archived path and the dedup would silently stop matching."""
    dock.log_discard(tmp_path, _unit(key="key-aaa"), "dropped\ton\tmerge", "/b/A/r0.sdf")
    assert dock._logged_drops(tmp_path) == {("key-aaa", "/b/A/r0.sdf")}
