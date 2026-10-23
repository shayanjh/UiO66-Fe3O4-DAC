"""Desorption metrics (Section 3.7 of the article and Table S11 of the SI).

The time series written by lammps/desorption/in.desorption_magnetic and in.desorption_uniform
(desorption_timeseries.txt) contains the step, the mean MOF temperature, the oxide temperature, and the number
of adsorbed CO2 molecules.

    released fraction  f(t) = (N0 - N(t)) / dN_ref
    t50, t90           times at which f reaches 0.5 and 0.9, where dN_ref is the
                       release reached in the uniform heating run
    t360               time at which the mean MOF temperature first reaches 360 K

Example
-------
    python desorption_metrics.py desorption_timeseries.txt 232.5
"""
import sys
import numpy as np

DT_FS = 0.5  # time step (fs)


def load_timeseries(path):
    """Read the fix ave/time output: step, T_MOF, T_Fe3O4, N_CO2."""
    a = np.loadtxt(path, comments="#")
    t_ps = a[:, 0] * DT_FS * 1e-3
    return t_ps, a[:, 1], a[:, 2], a[:, 3]


def first_crossing(t, y, level):
    """First time at which y reaches level (linear interpolation); NaN if never reached."""
    idx = np.where(y >= level)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return t[0]
    return t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1])


def desorption_metrics(t, T_mof, N, dN_ref, N0=None, n_avg=5):
    """Return released percentage, t50, t90, and t360 (ps)."""
    N0 = np.mean(N[:n_avg]) if N0 is None else N0
    f = (N0 - N) / dN_ref
    return {
        "released_pct": 100.0 * (N0 - np.mean(N[-n_avg:])) / N0,
        "t50_ps": first_crossing(t, f, 0.5),
        "t90_ps": first_crossing(t, f, 0.9),
        "t360_ps": first_crossing(t, T_mof, 360.0),
    }


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    t, T_mof, T_fe, N = load_timeseries(sys.argv[1])
    for k, v in desorption_metrics(t, T_mof, N, float(sys.argv[2])).items():
        print(f"{k:>14s} = {v:.1f}")
