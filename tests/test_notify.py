import logging

import pytest

from polyphar import notify as notify_mod
from polyphar.notify import notify


class Resp:
    def __init__(self, body):
        self._body = body

    def json(self):
        return self._body


class Session:
    def __init__(self, body=None, exc=None):
        self.posts = []
        self.body = body if body is not None else {"status": 1, "request": "abc"}
        self.exc = exc

    def post(self, url, data=None, timeout=None):
        self.posts.append({"url": url, "data": data, "timeout": timeout})
        if self.exc:
            raise self.exc
        return Resp(self.body)


@pytest.fixture
def creds(monkeypatch):
    monkeypatch.setenv("PUSHOVER_APP_TOKEN", "tok")
    monkeypatch.setenv("PUSHOVER_USER_KEY", "usr")


def _stub(monkeypatch, session):
    monkeypatch.setattr(notify_mod, "http_session", lambda *a, **k: session)
    return session


def test_notify_without_credentials_returns_false(monkeypatch, caplog):
    monkeypatch.setattr(notify_mod, "_warned_no_credentials", False)
    monkeypatch.delenv("PUSHOVER_APP_TOKEN", raising=False)
    monkeypatch.delenv("PUSHOVER_USER_KEY", raising=False)

    def no_network(*a, **k):
        raise AssertionError("network must not be touched without credentials")

    monkeypatch.setattr(notify_mod, "http_session", no_network)
    with caplog.at_level(logging.WARNING, logger="polyphar.notify"):
        assert notify("t", "m") is False
    assert len(caplog.records) == 1


def test_missing_credentials_warn_once_per_process(monkeypatch, caplog):
    monkeypatch.setattr(notify_mod, "_warned_no_credentials", False)
    monkeypatch.delenv("PUSHOVER_APP_TOKEN", raising=False)
    monkeypatch.delenv("PUSHOVER_USER_KEY", raising=False)
    with caplog.at_level(logging.WARNING, logger="polyphar.notify"):
        for _ in range(3):
            assert notify("t", "m") is False
    assert len(caplog.records) == 1


def test_notify_sends_expected_request(monkeypatch, creds):
    s = _stub(monkeypatch, Session())
    monkeypatch.setenv("POLYPHAR_RUN_NAME", "kinases")
    assert notify("fix 25%", "hello", priority=1) is True
    post = s.posts[0]
    assert post["url"] == "https://api.pushover.net/1/messages.json"
    assert post["timeout"] and post["timeout"] > 0
    d = post["data"]
    assert d["token"] == "tok" and d["user"] == "usr"
    assert d["title"] == "[kinases] fix 25%"
    assert d["message"] == "hello" and d["priority"] == 1


def test_notify_prefix_uses_run_dir_name(monkeypatch, creds, tmp_path):
    s = _stub(monkeypatch, Session())
    monkeypatch.delenv("POLYPHAR_RUN_NAME", raising=False)
    run = tmp_path / "gpcr_run"
    run.mkdir()
    notify("t", "m", run_dir=run)
    assert s.posts[0]["data"]["title"] == "[gpcr_run] t"


def test_notify_truncates(monkeypatch, creds):
    s = _stub(monkeypatch, Session())
    notify("T" * 1000, "M" * 5000)
    d = s.posts[0]["data"]
    assert len(d["message"]) <= 1024
    assert len(d["title"]) <= 250


def test_notify_never_raises_on_network_error(monkeypatch, creds, caplog):
    _stub(monkeypatch, Session(exc=ConnectionError("down")))
    with caplog.at_level(logging.WARNING, logger="polyphar.notify"):
        assert notify("t", "m") is False
    assert len(caplog.records) == 1


def test_notify_rejected_by_api(monkeypatch, creds, caplog):
    _stub(monkeypatch, Session(body={"status": 0, "errors": ["user identifier is invalid"]}))
    with caplog.at_level(logging.WARNING, logger="polyphar.notify"):
        assert notify("t", "m") is False
    assert "invalid" in caplog.text


@pytest.mark.pushover
def test_notify_live():
    from polyphar.common import http_session
    import os
    r = http_session().post(
        notify_mod.PUSHOVER_URL,
        data={"token": os.environ["PUSHOVER_APP_TOKEN"], "user": os.environ["PUSHOVER_USER_KEY"],
              "title": "polyphar", "message": "polyphar test notification"},
        timeout=10)
    assert r.json()["status"] == 1
    assert notify("live test", "polyphar test notification") is True
