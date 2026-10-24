import os
import sys
from datetime import datetime
import numpy as np

try:
    from PyQt6.QtWidgets import QApplication
except ImportError:
    from PyQt5.QtWidgets import QApplication

from SysID_GUI import MainWindow as CorrectorSysIDWindow, Worker as CorrectorWorker, finite_abs_max, update_amplitude


def fit_response(C, O, c0):
    C, O = np.asarray(C, dtype=float), np.asarray(O, dtype=float)
    R, err = np.full((O.shape[1], C.shape[1]), np.nan), np.full((O.shape[1], C.shape[1]), np.nan)
    cols = np.isfinite(O).mean(axis=0) > 0.5
    rows = np.all(np.isfinite(O[:, cols]), axis=1) & np.all(np.isfinite(C), axis=1)
    A = np.column_stack([C[rows] - c0, np.ones(rows.sum())])
    if rows.sum() <= A.shape[1]:
        return None, None
    coef = np.linalg.lstsq(A, O[rows][:, cols], rcond=None)[0]
    sigma2 = np.sum((O[rows][:, cols] - A @ coef) ** 2, axis=0) / (rows.sum() - A.shape[1])
    R[cols] = coef[:-1].T
    err[cols] = np.sqrt(np.outer(sigma2, np.diag(np.linalg.pinv(A.T @ A))[:-1]))
    return R, err


class RandomWorker(CorrectorWorker):
    def _run_impl(self):
        self.running = True
        self.paused = False
        I = self.interface
        corr = I.get_correctors(self.correctors)
        names = [str(name) for name in corr['names']]
        n = len(names)
        c0 = np.asarray(corr['bdes'], dtype=float)
        h = np.array([name in self.hcorrs for name in names])
        plane = np.where(h, "x", "y")
        index = [self.correctors.index(name) for name in names]
        kicks = np.where(h, np.asarray(self.hkicks, dtype=float)[index], np.asarray(self.vkicks, dtype=float)[index])
        target = np.where(h, self.max_osc_h, self.max_osc_v)
        max_curr = np.where(h, self.max_curr_h, self.max_curr_v)
        nb, npairs = len(self.bpms), self.Niter * n
        rows = {"x": slice(0, nb), "y": slice(nb, 2 * nb)}
        excited = np.zeros(n, dtype=bool)
        C, O, R = [], [], None

        def excitation(k):
            if k < n or R is None:
                return np.where(np.arange(n) == k % n, kicks, 0.0)
            s = np.random.uniform(-1, 1, n) * excited
            d = kicks * s
            for p, r in rows.items():
                m = (plane == p) & excited
                if m.any():
                    predicted = finite_abs_max(np.nan_to_num(R[r][:, m]) @ d[m])
                    d[m] *= min(target[m][0] / max(predicted, 1e-12), 2 / np.abs(s[m]).max())
            return d

        def data_file(tag, k):
            return os.path.join(self.output_dir, f'DATA_RND_{tag}{k:04d}.pkl')

        measured = [all(os.path.isfile(data_file(tag, k)) for tag in "pm") for k in range(npairs)]
        total_steps = max(npairs - sum(measured), 1)
        for k in range(npairs):
            if not self.running: break
            if self.paused:      self._await_user()

            if measured[k]:
                states = [self.state_class(filename=data_file(tag, k)) for tag in "pm"]
            else:
                d = excitation(k)
                states = []
                for sign in (1, -1):
                    values = c0 + sign * d
                    if I.set_correctors(names, np.where(max_curr > 0, np.clip(values, -max_curr, max_curr), values)) is False:
                        print(f"Skipping random excitation {k}: readback failed after retries")
                        break
                    if not self.running: break
                    if self.paused:      self._await_user()
                    states.append(I.get_state())
                if len(states) < 2:
                    continue
                for tag, state in zip("pm", states):
                    state.save(filename=data_file(tag, k))
                self.progress_value = self.progress_value + 1
                self.progress.emit(int(self.progress_value / total_steps * 100))

            orbits = [state.get_orbit(self.bpms) for state in states]
            for state, orbit in zip(states, orbits):
                readbacks = state.get_correctors(names)
                bact = dict(zip(map(str, readbacks['names']), readbacks['bact']))
                C.append([bact[name] for name in names])
                O.append(np.r_[orbit['x'], orbit['y']])

            j = k % n
            if k < n:
                excited[j] = True
                Op, Om = orbits
                Diff = (O[-2] - O[-1]) / 2.0
                nsamples = max(1, int(np.asarray(Op['stdx']).size))
                Err = np.r_[np.hypot(Op['stdx'], Om['stdx']), np.hypot(Op['stdy'], Om['stdy'])] / np.sqrt(nsamples)
                kicks[j] = update_amplitude(kicks[j], finite_abs_max(Diff[rows[plane[j]]]), target[j], max_curr[j])
            if k >= n - 1:
                R, err = fit_response(C, O, c0)
            if k >= n:
                if R is None:
                    continue
                R0 = np.abs(np.nan_to_num(R))
                single = np.where(h, R0[:nb].max(axis=0), R0[nb:].max(axis=0)) * kicks
                kicks = np.where(excited, [update_amplitude(*args) for args in zip(kicks, single, target, max_curr)], kicks)
                Diff, Err = R[:, j] * kicks[j], err[:, j] * kicks[j]
            for jj, i in enumerate(index):
                (self.hkicks if h[jj] else self.vkicks)[i] = kicks[jj]
            with open(os.path.join(self.output_dir, 'kicks.txt'), 'w') as f:
                for i, c in enumerate(self.correctors):
                    f.write(f'{c} {self.hkicks[i]} {self.vkicks[i]}\n')
            if excited[j]:
                self.plot_data.emit(Diff[:nb], Err[:nb], Diff[nb:], Err[nb:], self.bpms, names[j])

        self.running = False


class MainWindow(CorrectorSysIDWindow):
    WORKER_CLASS = RandomWorker

    def _window_title_prefix(self):
        return "SYSID RND"

    def _mode_directory_name(self, mode):
        return "Random_" + super()._mode_directory_name(mode)

def main():
    app = QApplication(sys.argv)
    from machine_interfaces.Backend import SelectInterface
    interface = SelectInterface.choose_acc_and_interface()
    if interface is None:
        return 1

    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = os.path.expanduser(f"~/CERN-Flight_Simulator-Data/{interface.get_name()}_{time_str}_Random")
    window = MainWindow(interface=interface, dir_name=dir_name)
    window.show()
    return app.exec()

if __name__ == "__main__":
    sys.exit(main())
