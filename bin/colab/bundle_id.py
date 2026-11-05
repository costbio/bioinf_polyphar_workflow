"""Colab exec: name the bundle on this VM, so a re-run of polyphar-colab can re-attach.

Prints BUNDLE_ID=<sha256 of run_key.json> (or BUNDLE_ID=none). The host compares it with
its own unmerged bundle before polling, so it never packs a different bundle's results.
"""
import hashlib
import pathlib

try:
    print("BUNDLE_ID=" + hashlib.sha256(
        pathlib.Path("/content/bundle/run_key.json").read_bytes()).hexdigest())
except OSError:
    print("BUNDLE_ID=none")
