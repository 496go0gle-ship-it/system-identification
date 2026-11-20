"""Record of a SysID run, written next to its DATA files (used by SysID_corrector_and_mover_GUI.py).

In the SysID data folder:
  sysid_settings.json   "settings": every option of the window (the latest run, used to restore it)
                        "runs": one entry per START or resume (start, end, status, reason, pairs measured,
                                and that run's settings); "running" left in the last entry means the program
                                did not reach the end
                        "excitations": one entry per measured +/- pair (kick, largest orbit change and where,
                                its signal-to-noise, the next kick)
  sysid_log.txt         every line printed while the measurement runs, with the time

Writing the record never stops the measurement: a failed write prints a warning and goes on.
"""

import json
import os
import subprocess
import sys
import threading
from datetime import datetime

import numpy as np

SETTINGS_FILE = "sysid_settings.json"
LOG_FILE = "sysid_log.txt"
RECORD_VERSION = 1


def now():
    return datetime.now().isoformat(timespec="seconds")


def _package_folder(module_name):
    import importlib
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        return None
    if getattr(module, "__file__", None):
        return os.path.dirname(os.path.abspath(module.__file__))
    paths = list(getattr(module, "__path__", []))  # namespace package: no __file__
    return os.path.abspath(paths[0]) if paths else None


def code_versions():
    """{repository: "<commit> (<branch>)" [+ local changes]}, or why git could not tell."""
    folders = {"system-identification": os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
               "response-matrix": _package_folder("response_matrix"),
               "machine-interfaces": _package_folder("machine_interfaces"),
               "machine-interfaces (ATF2 backend)": _package_folder("machine_interfaces_atf2")}
    out = {}
    for name, folder in folders.items():
        if not folder:
            continue
        try:
            def git(*args):
                return subprocess.run(["git", "-C", folder, *args], capture_output=True, text=True, timeout=5)
            head = git("rev-parse", "--short", "HEAD")
            if head.returncode != 0 or not head.stdout.strip():
                reason = (head.stderr.strip().splitlines() or ["no output"])[0]
                out[name] = f"unknown: git says '{reason}' ({folder})"
                continue
            branch = git("rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
            dirty = bool(git("status", "--porcelain", "--untracked-files=no").stdout.strip())
            out[name] = f"{head.stdout.strip()} ({branch})" + (" + local changes" if dirty else "")
        except Exception as exc:
            out[name] = f"unknown ({exc})"
    return out


def _clean(value):
    """JSON-friendly copy: numpy numbers to float, NaN to None."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def read_record(folder):
    path = os.path.join(folder, SETTINGS_FILE)
    if not os.path.isfile(path):
        return None
    with open(path) as file:
        return json.load(file)


class SysIDRecord:
    """sysid_settings.json of one data folder; safe to call from the worker thread."""

    def __init__(self, folder):
        self.folder = folder
        self.path = os.path.join(folder, SETTINGS_FILE)
        self._lock = threading.Lock()
        self.data = {"record_version": RECORD_VERSION, "settings": {}, "runs": [], "excitations": []}
        try:
            old = read_record(folder)
            if old:
                self.data.update({k: old.get(k, self.data[k]) for k in ("settings", "runs", "excitations")})
        except Exception as exc:
            print(f"WARNING: could not read {self.path} ({exc}); a new record is started")
        self.run_index = None

    def _write(self):
        try:
            os.makedirs(self.folder, exist_ok=True)
            tmp = self.path + ".tmp"
            with open(tmp, "w") as file:
                json.dump(_clean(self.data), file, indent=2)
            os.replace(tmp, self.path)  # never leave a half-written file
        except Exception as exc:
            print(f"WARNING: could not write {self.path} ({exc}); the measurement goes on")

    def begin(self, settings, resumed):
        with self._lock:
            self.data["settings"] = settings
            self.data["runs"].append({"started_at": now(), "status": "running", "resumed": bool(resumed),
                                      "pairs_measured": 0, "settings": settings})
            self.run_index = len(self.data["runs"]) - 1
            self._write()

    def add_excitation(self, entry):
        with self._lock:
            entry = {"time": now(), "run": self.run_index, **entry}
            self.data["excitations"].append(entry)
            if self.run_index is not None:
                self.data["runs"][self.run_index]["pairs_measured"] += 1
            self._write()

    def end(self, status, reason=""):
        with self._lock:
            if self.run_index is None:
                return
            run = self.data["runs"][self.run_index]
            run.update({"ended_at": now(), "status": status})
            if reason:
                run["reason"] = reason
            self._write()


def excitation_entry(actuator, axis, iteration, half_kick, unit, plus, minus, bpms, planes, nshots, next_kick):
    """Largest orbit change (half of plus - minus) over the planes, where, and its signal-to-noise.

    The noise of the half difference at a BPM is sqrt(std+^2 + std-^2) / (2 sqrt(shots)), from the
    shot-to-shot spread the two states hold.
    """
    op, om = plus.get_orbit(bpms), minus.get_orbit(bpms)
    best = {"max_orbit_change": 0.0, "at_bpm": None, "plane": None, "noise": None, "snr": None}
    best_snr = {"snr": 0.0, "at_bpm": None, "plane": None}
    shots = max(1, int(nshots or 1))
    for plane in planes:
        diff = (np.asarray(op[plane], float).ravel() - np.asarray(om[plane], float).ravel()) / 2.0
        noise = np.sqrt(np.asarray(op["std" + plane], float).ravel() ** 2
                        + np.asarray(om["std" + plane], float).ravel() ** 2) / (2.0 * np.sqrt(shots))
        if not np.any(np.isfinite(diff)):
            continue
        i = int(np.nanargmax(np.abs(diff)))
        if abs(diff[i]) > best["max_orbit_change"]:
            best = {"max_orbit_change": float(abs(diff[i])), "at_bpm": str(bpms[i]), "plane": plane,
                    "noise": float(noise[i]), "snr": float(abs(diff[i]) / noise[i]) if noise[i] > 0 else None}
        ratio = np.where(noise > 0, np.abs(diff) / np.where(noise > 0, noise, 1.0), np.nan)
        if np.any(np.isfinite(ratio)):
            j = int(np.nanargmax(ratio))
            if ratio[j] > best_snr["snr"]:
                best_snr = {"snr": float(ratio[j]), "at_bpm": str(bpms[j]), "plane": plane}
    return {"actuator": actuator, "axis": axis, "iteration": int(iteration), "half_kick": float(half_kick), "unit": unit,
            **best, "best_snr": best_snr, "shots": shots, "next_half_kick": float(next_kick)}


class LogTee:
    """Stands in for sys.stdout while a run is on: prints as before and appends each line, with the time,
    to the log file. Lines ending in '\\r' (the pause reminder) are logged once."""

    def __init__(self, path, stream):
        self.path = path
        self.stream = stream
        self._buffer = ""
        self._last = None
        self._lock = threading.Lock()

    def write(self, text):
        try:
            self.stream.write(text)
        except Exception:
            pass
        with self._lock:
            self._buffer += text
            lines = self._buffer.replace("\r", "\n").split("\n")
            self._buffer = lines.pop()
            for line in lines:
                key = line.rstrip(".").strip()
                if not key or (key == self._last and "PAUSED" in key):
                    continue
                self._last = key
                try:
                    with open(self.path, "a", encoding="utf-8") as file:
                        file.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {line}\n")
                except OSError:
                    pass
        return len(text)

    def flush(self):
        try:
            self.stream.flush()
        except Exception:
            pass

    def __getattr__(self, name):  # isatty, encoding, fileno, ... of the real stream
        return getattr(self.stream, name)


def install_log(folder):
    """Send everything printed from now on to folder/sysid_log.txt as well (replaces an earlier tee)."""
    stream = sys.stdout.stream if isinstance(sys.stdout, LogTee) else sys.stdout
    sys.stdout = LogTee(os.path.join(folder, LOG_FILE), stream)


def remove_log():
    if isinstance(sys.stdout, LogTee):
        if sys.stdout._buffer.strip():
            sys.stdout.write("\n")
        sys.stdout = sys.stdout.stream
