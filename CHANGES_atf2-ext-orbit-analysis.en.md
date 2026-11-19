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
