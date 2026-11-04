"""tmfetch talks to UniProt's accessions endpoint, not /search with an OR-query.

Why this is pinned: on 2026-11-04 UniProt's /uniprotkb/search began returning HTTP 500
for "accession:A OR accession:B ..." queries of 20 terms or more, which aborted the
`select` stage on a live run. /uniprotkb/accessions is the documented endpoint for
fetching a known list and is not subject to that limit.
"""

import pytest

from polyphar import tmfetch


class _Resp:
    def __init__(self, text):
        self.text = text

    def raise_for_status(self):
        pass


class _Session:
    """Records every request and replies with a TSV covering the accessions asked for."""

    def __init__(self, tm=None):
        self.calls = []
        self.tm = tm or {}

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, dict(params or {})))
        accs = (params or {}).get("accessions", "").split(",")
        lines = ["Entry\tTransmembrane"]
        for a in accs:
            if a:
                lines.append(f"{a}\t{self.tm.get(a, '')}")
        return _Resp("\n".join(lines))


def test_fetch_uses_the_accessions_endpoint_with_a_comma_separated_list():
    s = _Session()
    tmfetch.fetch(["P10244", "P35222"], session=s)

    assert len(s.calls) == 1
    url, params = s.calls[0]
    assert url == "https://rest.uniprot.org/uniprotkb/accessions"
    # the accessions travel as a list, sorted and comma-joined
    assert params["accessions"] == "P10244,P35222"
    assert params["fields"] == "accession,ft_transmem"
    assert params["format"] == "tsv"
    # the OR-query form is what broke; it must not come back
    assert "query" not in params
    assert "OR" not in str(params)


def test_fetch_batches_at_most_batch_accessions_per_call():
    accs = [f"P{i:05d}" for i in range(tmfetch.BATCH + 5)]
    s = _Session()
    tmfetch.fetch(accs, session=s)

    assert len(s.calls) == 2
    first = s.calls[0][1]["accessions"].split(",")
    second = s.calls[1][1]["accessions"].split(",")
    assert len(first) == tmfetch.BATCH
    assert len(second) == 5
    # UniProt documents a 100-accession ceiling for this endpoint
    assert tmfetch.BATCH <= 100


def test_fetch_parses_transmem_segments_and_leaves_soluble_proteins_empty():
    s = _Session(tm={"P08183": "TRANSMEM 50..70; TRANSMEM 120..140", "P10244": ""})
    got = tmfetch.fetch(["P08183", "P10244"], session=s)

    assert got["P08183"] == [(50, 70), (120, 140)]
    assert got["P10244"] == []


def test_an_accession_uniprot_does_not_return_is_an_error_not_a_silent_no_tm():
    class _Partial(_Session):
        def get(self, url, params=None, timeout=None):
            self.calls.append((url, dict(params or {})))
            return _Resp("Entry\tTransmembrane\nP10244\t")

    with pytest.raises(RuntimeError, match="transmembrane status is unknown"):
        tmfetch.fetch(["P10244", "P08183"], session=_Partial())
