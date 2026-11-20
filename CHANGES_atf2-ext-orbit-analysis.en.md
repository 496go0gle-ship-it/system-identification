# Changes on branch `atf2-ext-orbit-analysis` (system-identification)

Base: `flight-simulator-framework/system-identification` `master` at `3e7f4ab`.
Japanese version: [CHANGES_atf2-ext-orbit-analysis.ja.md](CHANGES_atf2-ext-orbit-analysis.ja.md)

New file: `SysID_corrector_and_mover_GUI.py`. `SysID_GUI.py` and `SysID_QM_GUI.py` are unchanged.

## What it does
One SysID run excites correctors and quadrupole movers together, in the order of the list (correctors first, then movers in beamline order). Each actuator is measured with the same procedure as the existing applications:

| Actuator | Procedure (as in) | Data file |
|---|---|---|
| Corrector | +/- kick around the current setting, adaptive kick (`SysID_GUI.py`) | `DATA_<corrector>_<p\|m><iter>.pkl` |
| Quadrupole mover | +/- displacement in x, then in y, mover restored after each axis (`SysID_QM_GUI.py`) | `DATA_<quad>_<x\|y>_<p\|m><iter>.pkl` |

All files go to one folder, so `ComputeResponseMatrix_GUI.py` (correctors) and `ComputeResponseMatrix_QM_GUI.py` (quadrupole movers) can both be run on it. The selection is saved as `actuators.txt` (used to resume), `correctors.txt` and `quadrupole_movers.txt`.

## Window
- One list with every corrector and every usable mover; the pattern field accepts both (e.g. `ZH*, QM1*FF`).
- Corrector settings: the existing fields (initial kick, max strength).
- Mover settings: new fields **Mover initial X/Y [um]** (default 100) and **Max mover range X/Y [um]** (default 1000), the same defaults as SysID QM.
- The target orbit excursion, number of cycles, BPM samples and measurement modes (orbit / energy / intensity) apply to both.

## Usable movers
For `ATF2_Ext` and `ATF2_Ext_RFT` only the 22 final-focus quadrupoles are offered: `QM16FF`, `QM15FF`, `QM14FF`, `QM13FF`, `QM12FF`, `QM11FF`, `QD10BFF`, `QD10AFF`, `QF9BFF`, `QF9AFF`, `QD8FF`, `QF7FF`, `QD6FF`, `QF5BFF`, `QF5AFF`, `QD4BFF`, `QD4AFF`, `QF3FF`, `QD2BFF`, `QD2AFF`, `QF1FF`, `QD0FF`.
- These are the quadrupoles whose mover read-backs are present in the 2026-11-19 machine data; the extraction-line quadrupoles `QF1X` ... `QF21X` have none (SysID QM still lists them).
- The research setup (`orbit_correction_setup_1101`) uses 21 of them: every one except `QM16FF`, which is upstream of the last correctors `ZH1FF` / `ZV1FF`.
- Other interfaces get every mover their `get_quadrupole_movers_names()` returns.

## Differences from running SysID and SysID QM separately
- **Resume:** a mover excitation whose DATA file exists is not measured again, in every iteration. SysID QM measures again in iterations after the first.
- **Plot marker:** the position of a mover on the plot falls back to the BPM on the same quadrupole (`QM15FF` → `MQM15FF`), because the ATF2 Ext sequence lists BPM names.

## Verification (RF-Track simulation, `ATF2_Ext_RFT`)
- Worker run with `ZH1X`, `ZV1X`, `QM15FF`, `QD0FF`, 1 iteration: 12 DATA files, 6 plots (`ZH1X`, `ZV1X`, `QM15FF:x`, `QM15FF:y`, `QD0FF:x`, `QD0FF:y`), correctors and movers back at their initial values, a second run measures nothing again.
- From that folder: corrector response matrix (`ZH1X`, `ZV1X`) and quadrupole-mover response matrix (`QM15FF`, `QD0FF`, x and y) are both computed.
- The window opens with 51 actuators (29 correctors + 22 movers) and writes the three selection files.
- Not checked on the real machine or by pressing START in the window.

## Selection files and logs (added)
- **Selection files:** "Load actuators" reads names separated by newlines, spaces or commas, with `#` comments. A mover may be written `QM15FF`, `qm15ff`, `QM15FF:x` or `MQM15FF` (the BPM on the quadrupole). After loading, the selected correctors and movers and the names not selected (quadrupole without a mover, or not in the list) are shown and printed. A leftover name pattern, which would override the selection, is cleared. Same rules as the combined BBA.
- **Terminal log:** a header with the correctors and movers and their units; each excitation prints "start -> excited value (step)", correctors in their strength unit (A on the machine), movers in um; movers also print the position read back after the move (`WARNING` beyond 15 um); after each actuator, the largest orbit change and the next amplitude.
- Checked (RF-Track): 6 correctors + 5 movers measured; the mover x/y DATA files are written and the log above is printed.

## Mover read-back check (added)
- After each move, the setpoint read back and the encoder change are printed. A `WARNING` is printed when the setpoint is more than 15 um from the request, or the encoder moved by an amount more than 15 um different from the request. The encoder is compared by magnitude; a reversed one is reported as "encoder counts opposite to the setpoint" (on the machine: `QF1FF` and `QD0FF` in x).

## Settings, progress and log record, and restoring settings (added)
- **Why:** in the machine SysID of 2026-11-19, many correctors' ±0.01 A excitations were buried in BPM noise, and the ZH7X column was spoiled by an orbit jump during its measurement. Afterwards there was no way to tell which settings were used and how well each excitation worked. The record is the same idea as the combined BBA's.
- **Change** (new file `system_identification/run_record.py`):
  - `sysid_settings.json` in the data folder, written when START is pressed.
    - `settings`: every option of the window:
      - interface, folder, measuring mode
      - actuators, BPMs
      - initial kicks, target orbit excursion, max current
      - mover initial step and range, Niter, BPM samples, units
      - code commits (or git's reason when it cannot tell)
    - `runs`: one entry per START or resume, appended, never overwritten.
      - start and end time, pairs measured
      - how the run ended: finished, stopped (STOP, or a corrector did not return), error (with the reason)
      - still `running` means the program did not reach the end
    - `excitations`: one per measured ± pair, with the kick, the largest orbit change and its BPM, its signal/noise, and the next kick.
  - Each pair prints its `signal/noise`, with a WARNING below 3 ("this response is unreliable; increase the excitation or BPM samples"). The noise is √(σ₊² + σ₋²) / (2√shots) from the shot spread of the two states.
  - Every line printed during the run is appended, with the time, to `sysid_log.txt` (mode change and restore included).
  - Choosing a folder holding `sysid_settings.json` with the folder "…" button puts the last run's settings back in the window. It also shows how that run ended and which excitations were buried in noise.
  - A failure to write the record never stops the measurement (a warning is printed).
- **Checked (RF-Track, `ATF2_Ext_RFT`, 2 correctors and 2 movers):**
  - A finished run gave `finished` with 6 pairs and a 50-line log.
  - STOP after the first pair gave `stopped` (STOP pressed).
  - `set_correctors` raising during the measurement gave `error` with the reason.
  - A resume in the same folder with Niter 2 gave a second `runs` entry (resumed), skipped iteration 0 and measured the 6 pairs of iteration 1.
  - Restoring into a new window brought back the settings and selection.
  - Not yet checked on the machine.
- Noticed, unchanged: `SysID_GUI.py` divides its error bars by the square root of the number of BPMs instead of the number of shots (`nsamples = stdx.size`). This affects the plotted error bars only; the record's signal/noise uses the shot count.
