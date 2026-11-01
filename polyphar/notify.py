"""Pushover notifications. Optional: never raises, never fails a run.

Credentials come only from PUSHOVER_APP_TOKEN and PUSHOVER_USER_KEY.
Limits (https://pushover.net/api, checked 2026-10-31): message 1024 chars,
title 250 chars.
"""

from __future__ import annotations

import logging
import os
from datetime import datetime
from pathlib import Path

from polyphar.common import http_session

PUSHOVER_URL = "https://api.pushover.net/1/messages.json"
MESSAGE_LIMIT = 1024
TITLE_LIMIT = 250
TIMEOUT_S = 10

log = logging.getLogger("polyphar.notify")


def run_name(run_dir=None) -> str:
    """Name shown in the title prefix: $POLYPHAR_RUN_NAME (the container's
    /work has no useful name), else the run folder name, else the cwd name."""
    name = os.environ.get("POLYPHAR_RUN_NAME")
    if name:
        return name
    return Path(run_dir).resolve().name if run_dir else Path.cwd().name


NOTIFICATIONS_LOG = Path("logs") / "notifications.log"


def _record(run_dir, stamp: str, priority: int, sent: bool, title: str, message: str) -> None:
    """Append one line to <run_dir>/logs/notifications.log: timestamp, priority,
    sent / not sent, title, message (newlines and tabs escaped). Never raises."""
    if run_dir is None:
        return
    try:
        esc = lambda x: str(x).replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")
        path = Path(run_dir) / NOTIFICATIONS_LOG
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{stamp}\t{priority}\t{'sent' if sent else 'not sent'}\t"
                     f"{esc(title)}\t{esc(message)}\n")
    except Exception as exc:
        log.warning("could not record the notification: %s: %s", type(exc).__name__, exc)


def notify(title: str, message: str, priority: int = 0, run_dir=None) -> bool:
    """Send one push notification. True only if Pushover accepted it. With run_dir,
    every notification (sent or not, e.g. without credentials) is also recorded in
    <run_dir>/logs/notifications.log."""
    stamp = datetime.now().astimezone().isoformat(timespec="microseconds")
    sent = _send(title, message, priority, run_dir)
    _record(run_dir, stamp, priority, sent, title, message)
    return sent


def _send(title: str, message: str, priority: int, run_dir) -> bool:
    try:
        token = os.environ.get("PUSHOVER_APP_TOKEN")
        user = os.environ.get("PUSHOVER_USER_KEY")
        if not token or not user:
            log.warning("pushover credentials not set; notification skipped")
            return False
        data = {
            "token": token,
            "user": user,
            "title": f"[{run_name(run_dir)}] {title}"[:TITLE_LIMIT],
            "message": message[:MESSAGE_LIMIT],
        }
        if priority:
            data["priority"] = priority
        r = http_session(retries=2, backoff=1.0).post(PUSHOVER_URL, data=data,
                                                      timeout=TIMEOUT_S)
        body = r.json()
        if body.get("status") == 1:
            return True
        log.warning("pushover rejected the notification: %s", body.get("errors", body))
        return False
    except Exception as exc:  # network, JSON, anything: notifications are optional
        log.warning("pushover notification failed: %s: %s", type(exc).__name__, exc)
        return False
