"""SysID application that excites correctors and quadrupole movers in one run.

Correctors are measured as in SysID_GUI.py and movers as in SysID_QM_GUI.py, into the
same folder and with the same file names:

- corrector:        DATA_<corrector>_<p|m><iteration>.pkl
- quadrupole mover: DATA_<quadrupole>_<x|y>_<p|m><iteration>.pkl

so ComputeResponseMatrix_GUI.py (correctors) and ComputeResponseMatrix_QM_GUI.py
(quadrupole movers) can both be run on the folder. The selection is saved as
actuators.txt (everything, used to resume), correctors.txt and quadrupole_movers.txt.

Unlike SysID_QM_GUI.py, a mover excitation whose DATA file already exists is not
measured again in later iterations, so a stopped run resumes like the corrector SysID.
"""

import os
import re
import sys
import time
from datetime import datetime
from enum import Enum

import numpy as np

try:
    from PyQt6.QtWidgets import QApplication, QDoubleSpinBox, QFileDialog, QLabel, QLineEdit, QMessageBox
except ImportError:
    from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QFileDialog, QLabel, QLineEdit, QMessageBox

from SysID_GUI import (
    MainWindow as CorrectorSysIDWindow,
    Worker as CorrectorWorker,
    finite_abs_max,
    update_amplitude,
)
from SysID_QM_GUI import QMWorker

# Quadrupoles of the ATF2 extraction line whose movers can be driven: the final-focus
# quadrupoles. Their mover read-backs (x, y, roll) are present in the machine data of
# 2026-11-19; the extraction-line quadrupoles QF1X ... QF21X have no mover.
ATF2_EXT_MOVERS = (
    "QM16FF", "QM15FF", "QM14FF", "QM13FF", "QM12FF", "QM11FF",
    "QD10BFF", "QD10AFF", "QF9BFF", "QF9AFF", "QD8FF", "QF7FF", "QD6FF",
    "QF5BFF", "QF5AFF", "QD4BFF", "QD4AFF", "QF3FF", "QD2BFF", "QD2AFF",
    "QF1FF", "QD0FF",
)
# Interfaces (get_name()) for which only the movers listed here are offered.
USABLE_MOVERS = {"ATF2_Ext": ATF2_EXT_MOVERS, "ATF2_Ext_RFT": ATF2_EXT_MOVERS}

DEFAULT_MOVER_STEP_UM = 100.0
DEFAULT_MOVER_RANGE_UM = 1000.0
MOVER_TOLERANCE_UM = 15.0  # same tolerance as apply_qmag_xyroll in the ATF2 interfaces


class CombinedMode(Enum):
    """Takes the place of ActuatorMode in folder names and saved settings."""
    KickerQM = "Correctors and quadrupole movers"


def read_selection_file(path):
    """Names in a selection file: one or more per line, separated by commas/spaces; '#' starts a comment."""
    names = []
    with open(path) as file:
        for line in file:
            names += [token for token in re.split(r"[,;\s]+", line.split("#", 1)[0]) if token]
    return names


def match_selection(requested, offered, movers, quadrupoles=()):
    """Match names from a selection file to the actuators offered in the list.

    Accepted spellings for a mover: QM15FF, qm15ff, QM15FF:x, QM15FF_y, and MQM15FF (the BPM on it).
    Returns (found, ignored): found in the order requested, ignored as {name: reason}.
    """
    offered = [str(name) for name in offered]
    by_key = {name.upper(): name for name in offered}
    mover_keys = {str(name).upper() for name in movers}
    quad_keys = {str(name).upper() for name in quadrupoles}
    found, ignored = [], {}
    for raw in requested:
        key = str(raw).strip().upper()
        candidates = [key]
        plane = re.match(r"^(.+?)[:_]([XY])$", key)
        if plane:
            candidates.append(plane.group(1))
        candidates += [candidate[1:] for candidate in list(candidates) if candidate.startswith("M")]
        match = next((by_key[c] for c in candidates if c in by_key and (c == key or c in mover_keys)), None)
        if match is not None:
            if match not in found:
                found.append(match)
        elif any(c in quad_keys for c in candidates):
            ignored[raw] = "quadrupole without a usable mover in this interface"
        else:
            ignored[raw] = "not in this window's list"
    return found, ignored


def _ignored_text(ignored):
    if not ignored:
        return ""
    return "Not selected:\n" + "\n".join(f"  {name}: {reason}" for name, reason in ignored.items())


def selection_report(filename, found, ignored, movers):
    """Text shown (and printed) after loading a selection file."""
    mover_set = set(movers)
    correctors = [name for name in found if name not in mover_set]
    chosen_movers = [name for name in found if name in mover_set]
    lines = [f"Loaded {os.path.basename(filename)}:",
             f"  correctors ({len(correctors)}): {', '.join(correctors) or 'none'}",
             f"  quadrupole movers ({len(chosen_movers)}): {', '.join(chosen_movers) or 'none'}"]
    if ignored:
        lines.append(_ignored_text(ignored))
        if any("mover" in reason for reason in ignored.values()) and movers:
            lines.append(f"  Movers offered here: {', '.join(movers)}")
    return "\n".join(lines)


def usable_movers(interface):
    """Quadrupole movers this interface can drive, in the interface's beamline order."""
    if not hasattr(interface, "apply_qmag_xyroll") or not hasattr(interface, "get_quadrupole_movers_names"):
        return []
    names = [str(name) for name in interface.get_quadrupole_movers_names()]
    allowed = USABLE_MOVERS.get(interface.get_name())
    if allowed is not None:
        names = [name for name in names if name in allowed]
    return names


class CombinedWorker(QMWorker):
    """Run the corrector and quadrupole-mover excitations of one SysID in beamline order."""

    def __init__(self, *args, movers=(), mover_initial_um=(DEFAULT_MOVER_STEP_UM, DEFAULT_MOVER_STEP_UM),
                 mover_max_range_um=(DEFAULT_MOVER_RANGE_UM, DEFAULT_MOVER_RANGE_UM), corrector_unit="", bpm_unit="mm", **kwargs):
        super().__init__(*args, **kwargs)
        self.corrector_unit = corrector_unit
        self.bpm_unit = bpm_unit
        self.actuator_mode = CombinedMode.KickerQM
        self.movers = set(movers)
        self.mover_max_x, self.mover_max_y = (float(value) for value in mover_max_range_um)
        self.hkicks = np.asarray(self.hkicks, dtype=float)
        self.vkicks = np.asarray(self.vkicks, dtype=float)
        for index, name in enumerate(self.correctors):
            if name in self.movers:
                self.hkicks[index], self.vkicks[index] = (float(value) for value in mover_initial_um)

    # --- bookkeeping -------------------------------------------------------------------
    def _checkpoint(self):
        """Wait while paused; return False once the run has been stopped."""
        if self.running and self.paused:
            self._await_user()
        return self.running

    def _excitations(self, name):
        return (("x", "+"), ("x", "-"), ("y", "+"), ("y", "-")) if name in self.movers else ((None, "+"), (None, "-"))

    def _file(self, name, axis, sign, iteration):
        if axis is None:
            return os.path.join(self.output_dir, f"DATA_{name}_{'p' if sign == '+' else 'm'}{iteration:04d}.pkl")
        return self._data_filename(name, axis, sign, iteration)

    def _step_done(self, name, axis, iteration):
        return all(os.path.isfile(self._file(name, axis, sign, iteration)) for sign in "+-")

    def _report(self, diff_x, err_x, diff_y, err_y, label):
        self.plot_data.emit(diff_x, err_x, diff_y, err_y, self.bpms, label)
        self.progress_value += 1
        self.progress.emit(int(self.progress_value / self._total_steps * 100))

    def _orbit_difference(self, plus, minus):
        orbit_plus, orbit_minus = plus.get_orbit(self.bpms), minus.get_orbit(self.bpms)
        diff_x = (orbit_plus["x"] - orbit_minus["x"]) / 2.0
        diff_y = (orbit_plus["y"] - orbit_minus["y"]) / 2.0
        nsamples = max(1, int(np.asarray(orbit_plus["stdx"]).size))
        err_x = np.sqrt(np.square(orbit_plus["stdx"]) + np.square(orbit_minus["stdx"])) / np.sqrt(nsamples)
        err_y = np.sqrt(np.square(orbit_plus["stdy"]) + np.square(orbit_minus["stdy"])) / np.sqrt(nsamples)
        return diff_x, err_x, diff_y, err_y

    def _pause_after_measurement(self, seconds):
        start = time.monotonic()
        while self.running and time.monotonic() - start < seconds:
            time.sleep(0.05)

    def _mover_readback(self, magnet, axis):
        """Mover position after a move [um]: encoder (xact/yact) if the interface has one, else the setpoint."""
        try:
            q = self.interface.get_quadrupoles([magnet])
        except Exception as exc:
            print(f"  {magnet}: could not read the mover back ({exc})")
            return None
        names = [str(name) for name in q.get("names", [])]
        if magnet not in names:
            return None
        index = names.index(magnet)
        for key in (f"{axis}act", f"{axis}des"):
            if key in q:
                value = float(np.asarray(q[key], dtype=float).ravel()[index])
                if np.isfinite(value):
                    return value
        return None

    # --- main loop ---------------------------------------------------------------------
    def _run_impl(self):
        self.running = True
        self.paused = False
        self.progress_value = 0
        pending = sum(
            not self._step_done(name, axis, iteration)
            for iteration in range(self.Niter)
            for name in self.correctors
            for axis in (("x", "y") if name in self.movers else (None,))
        )
        self._total_steps = max(pending, 1)
        movers = [name for name in self.correctors if name in self.movers]
        correctors = [name for name in self.correctors if name not in self.movers]
        print(f"SysID correctors + movers: {len(correctors)} correctors [{self.corrector_unit}] "
              f"({', '.join(correctors) or 'none'}), {len(movers)} quadrupole movers [um] ({', '.join(movers) or 'none'}), "
              f"{self.Niter} iteration(s), {pending} excitation(s) to measure")

        for iteration in range(self.Niter):
            if not self._checkpoint():
                break
            for index, name in enumerate(self.correctors):
                if not self._checkpoint():
                    break
                if name in self.movers:
                    self._measure_mover(index, name, iteration)
                else:
                    self._measure_corrector(index, name, iteration)
        self.running = False

    # --- correctors (same procedure as SysID_GUI.Worker) --------------------------------
    def _measure_corrector(self, index, corrector, iteration):
        interface = self.interface
        if corrector in self.hcorrs:
            kick, target, max_current, plane = float(self.hkicks[index]), float(self.max_osc_h), float(self.max_curr_h), "x"
        elif corrector in self.vcorrs:
            kick, target, max_current, plane = float(self.vkicks[index]), float(self.max_osc_v), float(self.max_curr_v), "y"
        else:
            print(f"Corrector {corrector} not in h/v correctors")
            return

        def clamp(value):
            return value if max_current == 0.0 else max(-max_current, min(value, max_current))

        def set_corrector(current):
            return interface.set_correctors(corrector, current) is not False

        initial = float(np.asarray(interface.get_correctors(corrector)["bdes"], dtype=float).ravel()[0]) + 0.0  # no "-0" in the log
        changed = measured = skipped = False
        states = {}
        for sign in "+-":
            filename = self._file(corrector, None, sign, iteration)
            if os.path.isfile(filename):
                states[sign] = self.state_class(filename=filename)
                continue
            if not self._checkpoint():
                break
            value = clamp(initial + (kick if sign == "+" else -kick))
            print(f"Corrector {corrector} '{sign}' excitation: {initial:+.5g} -> {value:+.5g} {self.corrector_unit} "
                  f"(step {value - initial:+.5g} {self.corrector_unit})")
            if not set_corrector(value):
                print(f"Skipping {corrector}: '{sign}' readback failed after retries")
                skipped = True
                break
            changed = True
            state = interface.get_state()
            state.save(filename=filename)
            states[sign] = state
            measured = True

        if (changed or skipped) and not set_corrector(initial):
            print(f"Stopping SysID: {corrector} did not return to its initial current")
            self.running = False
        if skipped or len(states) < 2:
            return

        diff_x, err_x, diff_y, err_y = self._orbit_difference(states["+"], states["-"])
        if measured:
            self._report(diff_x, err_x, diff_y, err_y, corrector)
        observed = finite_abs_max(diff_x if plane == "x" else diff_y)
        new_kick = 0.8 * update_amplitude(kick, observed, target, max_current) + 0.2 * kick
        print(f"  {corrector}: max|orbit change| in {plane} = {observed:.4g} {self.bpm_unit} for +/-{kick:.5g} {self.corrector_unit}; "
              f"next excitation {new_kick:.5g} {self.corrector_unit}")
        if plane == "x":
            self.hkicks[index] = new_kick
        else:
            self.vkicks[index] = new_kick
        self._save_kicks(self.hkicks, self.vkicks)
        if measured:
            self._pause_after_measurement(1.0)

    # --- quadrupole movers (same procedure as SysID_QM_GUI.QMWorker) --------------------
    def _measure_mover(self, index, magnet, iteration):
        interface = self.interface
        state0 = self._quadrupole_state(interface, magnet)
        if state0 is None:
            return
        x0, y0, roll0 = state0

        for axis in ("x", "y"):
            if not self._checkpoint():
                return
            amplitude = float(self.hkicks[index] if axis == "x" else self.vkicks[index])
            target = float(self.max_osc_h if axis == "x" else self.max_osc_v)
            max_range = self.mover_max_x if axis == "x" else self.mover_max_y
            if max_range > 0:
                amplitude = min(amplitude, max_range)

            states = {}
            moved = measured = False
            try:
                for sign in "+-":
                    filename = self._file(magnet, axis, sign, iteration)
                    if os.path.isfile(filename):
                        states[sign] = self.state_class(filename=filename)
                        continue
                    if not self._checkpoint():
                        break
                    displacement = amplitude if sign == "+" else -amplitude
                    start = x0 if axis == "x" else y0
                    print(f"Quadrupole mover {magnet} {axis} '{sign}' excitation: {axis} {start:+.1f} -> {start + displacement:+.1f} um "
                          f"(step {displacement:+.1f} um)")
                    moved = True
                    if axis == "x":
                        interface.apply_qmag_xyroll(magnet, x0 + displacement, y0, roll0)
                    else:
                        interface.apply_qmag_xyroll(magnet, x0, y0 + displacement, roll0)
                    readback = self._mover_readback(magnet, axis)
                    if readback is not None:
                        off = readback - (start + displacement)
                        print(f"  {magnet} {axis} readback {readback:+.1f} um" + (
                            f"  WARNING: {off:+.1f} um from the request" if abs(off) > MOVER_TOLERANCE_UM else ""))
                    state = interface.get_state()
                    state.save(filename=filename)
                    states[sign] = state
                    measured = True
            finally:
                if moved:
                    try:
                        interface.apply_qmag_xyroll(magnet, x0, y0, roll0)
                    except Exception as exc:
                        print(f"WARNING: failed to restore {magnet} mover state ({exc})")
            if len(states) < 2:
                return

            diff_x, err_x, diff_y, err_y = self._orbit_difference(states["+"], states["-"])
            if measured:
                self._report(diff_x, err_x, diff_y, err_y, f"{magnet}:{axis}")
            observed = max(finite_abs_max(diff_x), finite_abs_max(diff_y))
            used = amplitude
            amplitude = 0.8 * update_amplitude(amplitude, observed, target, max_range) + 0.2 * amplitude
            print(f"  {magnet}:{axis}: max|orbit change| x = {finite_abs_max(diff_x):.4g}, y = {finite_abs_max(diff_y):.4g} {self.bpm_unit} "
                  f"for +/-{used:.1f} um; next excitation {amplitude:.1f} um")
            if axis == "x":
                self.hkicks[index] = amplitude
            else:
                self.vkicks[index] = amplitude
            self._save_kicks(self.hkicks, self.vkicks)
            if measured:
                self._pause_after_measurement(0.2)


class MainWindow(CorrectorSysIDWindow):
    """SysID window offering correctors and quadrupole movers in one list."""

    ACTUATOR_MODE = CombinedMode.KickerQM

    def __init__(self, interface, dir_name):
        self.mover_names = usable_movers(interface)
        self._mover_settings = None
        super().__init__(interface, dir_name)
        self._add_mover_controls()
        self.setWindowTitle(f"{self._window_title_prefix()} - {self.interface.get_name()} [Idle]")
        self.save_correctors_button.setText("Save actuators")
        self.load_correctors_button.setText("Load actuators")
        self.clear_correctors_button.setText("Clear actuators")
        self.pattern_corrs_input.setPlaceholderText("e.g. ZH*, ZV*, QM1*FF, QD*FF")
        # Own save/load: an exception in a Qt slot aborts the whole application (PyQt6), and the
        # selection may be split between correctors.txt and quadrupole_movers.txt.
        for button, slot in ((self.save_correctors_button, self._save_actuators_clicked),
                             (self.load_correctors_button, self._load_actuators_clicked)):
            button.clicked.disconnect()
            button.clicked.connect(slot)

    def WORKER_CLASS(self, *args, **kwargs):  # called by MainWindow._launch_worker like a class
        initial_x, initial_y, range_x, range_y = self._mover_settings
        return CombinedWorker(*args, movers=self.mover_names, mover_initial_um=(initial_x, initial_y),
                              mover_max_range_um=(range_x, range_y), corrector_unit=self.corrs_unit,
                              bpm_unit=self.bpm_unit, **kwargs)

    def _window_title_prefix(self):
        return "SYSID Correctors + QM"

    def _add_mover_controls(self):
        grid = self.options_grid
        row = grid.rowCount()
        self.mover_initial_x_label = QLabel("Mover initial X [um]")
        self.mover_initial_x_input = QLineEdit(str(DEFAULT_MOVER_STEP_UM))
        self.mover_initial_y_label = QLabel("Mover initial Y [um]")
        self.mover_initial_y_input = QLineEdit(str(DEFAULT_MOVER_STEP_UM))
        self.mover_range_label = QLabel("Max mover range [um]")
        self.mover_range_x_label = QLabel("X:")
        self.mover_range_y_label = QLabel("Y:")
        self.mover_range_x_spinbox = QDoubleSpinBox()
        self.mover_range_y_spinbox = QDoubleSpinBox()
        for spinbox in (self.mover_range_x_spinbox, self.mover_range_y_spinbox):
            spinbox.setMaximum(1e6)
            spinbox.setSingleStep(10.0)
            spinbox.setValue(DEFAULT_MOVER_RANGE_UM)
        grid.addWidget(self.mover_initial_x_label, row, 0)
        grid.addWidget(self.mover_initial_x_input, row, 1)
        grid.addWidget(self.mover_range_label, row, 2)
        grid.addWidget(self.mover_range_x_label, row, 3)
        grid.addWidget(self.mover_range_x_spinbox, row, 4)
        grid.addWidget(self.mover_range_y_label, row, 5)
        grid.addWidget(self.mover_range_y_spinbox, row, 6)
        grid.addWidget(self.mover_initial_y_label, row + 1, 0)
        grid.addWidget(self.mover_initial_y_input, row + 1, 1)
        visible = bool(self.mover_names)
        for widget in (self.mover_initial_x_label, self.mover_initial_x_input, self.mover_initial_y_label,
                       self.mover_initial_y_input, self.mover_range_label, self.mover_range_x_label,
                       self.mover_range_y_label, self.mover_range_x_spinbox, self.mover_range_y_spinbox):
            widget.setVisible(visible)

    def _available_actuators(self):
        return super()._available_actuators() + list(self.mover_names)

    def _sort_actuators(self, names):
        names = [str(name) for name in names]
        movers = [name for name in self.mover_names if name in names]
        correctors = self._sort_elements([name for name in names if name not in self.mover_names], which="corrs")
        return correctors + movers

    def _refresh_actuator_labels(self):
        super()._refresh_actuator_labels()
        self.correctorsGroup.setTitle("Correctors and quadrupole movers")

    def _validate_start(self):
        try:
            initial_x = float(self.mover_initial_x_input.text())
            initial_y = float(self.mover_initial_y_input.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid mover setting", "Mover initial X and Y must be numbers (um).")
            return False
        if initial_x <= 0 or initial_y <= 0:
            QMessageBox.critical(self, "Invalid mover setting", "Mover initial X and Y must be positive (um).")
            return False
        self._mover_settings = (initial_x, initial_y, self.mover_range_x_spinbox.value(), self.mover_range_y_spinbox.value())
        selected = {item.text() for item in self.correctors_list.selectedItems()}
        if selected & set(self.mover_names) and not hasattr(self.interface, "apply_qmag_xyroll"):
            QMessageBox.critical(self, "Movers not available", "This interface does not expose quadrupole-mover controls.")
            return False
        return True

    def _actuator_selection_filename(self):
        return "actuators.txt"

    # --- saving and loading the selection ------------------------------------------------
    SPLIT_SELECTION_FILES = ("correctors.txt", "quadrupole_movers.txt")

    def _selection_folder(self):
        folder = os.path.expanduser(os.path.expandvars(self.working_directory_input.text().strip()))
        return folder if os.path.isdir(folder) else os.path.expanduser("~/CERN-Flight_Simulator-Data")

    @staticmethod
    def _read_names(path):
        with open(path) as file:
            return [line.strip() for line in file if line.strip()]

    def _save_actuators_clicked(self):
        try:
            names = [item.text() for item in self.correctors_list.selectedItems()] or self._available_actuators()
            names = self._sort_actuators(names)
            start = os.path.join(self._selection_folder(), self._actuator_selection_filename())
            filename, _ = QFileDialog.getSaveFileName(self, "Save actuators", start, "Text (*.txt)")
            if not filename:
                return
            os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
            with open(filename, "w") as file:
                file.writelines(f"{name}\n" for name in names)
            movers = sum(name in self.mover_names for name in names)
            QMessageBox.information(self, "Save actuators",
                                    f"Saved {len(names) - movers} correctors and {movers} quadrupole movers to\n{filename}")
        except Exception as exc:
            QMessageBox.critical(self, "Save actuators", f"Could not save the selection:\n{exc}\n\n"
                                 "Choose a folder you can write to (e.g. under ~/CERN-Flight_Simulator-Data).")

    def _load_actuators_clicked(self):
        try:
            start = os.path.join(self._selection_folder(), self._actuator_selection_filename())
            filename, _ = QFileDialog.getOpenFileName(self, "Load actuators", start, "Text (*.txt)")
            if not filename:
                return
            names = read_selection_file(filename)
            folder = os.path.dirname(filename)
            if os.path.basename(filename) in self.SPLIT_SELECTION_FILES:
                # correctors.txt and quadrupole_movers.txt of one SysID folder are loaded together
                for other in self.SPLIT_SELECTION_FILES:
                    path = os.path.join(folder, other)
                    if other != os.path.basename(filename) and os.path.isfile(path):
                        names += read_selection_file(path)
            names = list(dict.fromkeys(names))
            items = {self.correctors_list.item(i).text(): self.correctors_list.item(i) for i in range(self.correctors_list.count())}
            quadrupoles = getattr(self.interface, "quadrupoles", None) or getattr(self.interface, "Qmagnames", None) or []
            found, ignored = match_selection(names, items, self.mover_names, quadrupoles)
            print(selection_report(filename, found, ignored, self.mover_names))
            if not found:
                QMessageBox.warning(self, "Load actuators", f"No actuator of this window is listed in\n{filename}\n\n"
                                    + _ignored_text(ignored))
                return
            if self.pattern_corrs_input.text().strip():
                # A pattern would re-select the list the next time it is applied
                self.pattern_corrs_input.blockSignals(True)
                self.pattern_corrs_input.clear()
                self.pattern_corrs_input.blockSignals(False)
            self.correctors_list.clearSelection()
            for name in found:
                items[name].setSelected(True)
            QMessageBox.information(self, "Load actuators", selection_report(filename, found, ignored, self.mover_names))
        except Exception as exc:
            QMessageBox.critical(self, "Load actuators", f"Could not load the selection:\n{exc}")

    def _actuator_label(self):
        return "Actuator"

    def _save_names_if_missing(self, base_dir, filename, names):
        super()._save_names_if_missing(base_dir, filename, names)
        if filename == self._actuator_selection_filename():
            # Selection files of the corrector and QM response-matrix GUIs
            names = [str(name) for name in names]
            correctors = [name for name in names if name not in self.mover_names]
            movers = [name for name in names if name in self.mover_names]
            if correctors:
                super()._save_names_if_missing(base_dir, "correctors.txt", correctors)
            if movers:
                super()._save_names_if_missing(base_dir, "quadrupole_movers.txt", movers)

    def _device_position_on_bpm_axis(self, device_name, bpm_names):
        # Movers are named after their quadrupole (QM15FF); the sequence may list the BPM on it (MQM15FF).
        sequence = list(self.interface.get_sequence())
        for candidate in (device_name, f"M{device_name}"):
            if candidate in sequence:
                return super()._device_position_on_bpm_axis(candidate, bpm_names)
        return -0.5


def main():
    app = QApplication(sys.argv)
    from machine_interfaces.Backend import SelectInterface

    interface = SelectInterface.choose_acc_and_interface()
    if interface is None:
        return 1

    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = os.path.expanduser(f"~/CERN-Flight_Simulator-Data/{interface.get_name()}_{time_str}_KickerQM")
    window = MainWindow(interface, dir_name)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
