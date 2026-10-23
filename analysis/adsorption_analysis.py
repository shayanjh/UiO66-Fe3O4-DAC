"""Analysis of computed adsorption isotherms (Sections 3.6 and 4.2 to 4.4 of the article).

Functions
---------
uptake_from_raw    mean uptake (mmol/g) from a raw_P=*.txt file written by the LAMMPS input scripts
langmuir_fit       Langmuir fit q = Nsat*b*P/(1+b*P) with 95 percent confidence intervals and R^2
vant_hoff          enthalpy of adsorption from the temperature dependence of b
gibbs_entropy      standard Gibbs free energy and entropy of adsorption (standard pressure P0)
qst_clausius       isosteric heat at constant loading from isotherms at several temperatures

Example
-------
    python adsorption_analysis.py ../data/computed/gcmc_dry_full_data.csv
"""
import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats

R = 8.314462618e-3  # kJ/(mol K)


def uptake_from_raw(path, cell_mass, column=3):
    """Mean uptake and its standard deviation (mmol/g) from a raw_P=<P>.txt file.

    path      : output of fix ave/time in lammps/adsorption/in.*
    cell_mass : molar mass of the simulated framework (g/mol), i.e., the sum of
                the masses of all framework atoms (of the UiO-66 atoms only for
                the composites, so that the uptake is per gram of UiO-66)
    column    : column of the number of molecules; the files contain
                TimeStep, step, P, N, ... so that N is column 3 for pure gases
                (for mixtures, use the column of the species of interest)
    """
    a = np.loadtxt(path, comments="#")
    a = np.atleast_2d(a)
    n = a[:, column]
    return 1e3 * n.mean() / cell_mass, 1e3 * n.std() / cell_mass


def langmuir(P, Nsat, b):
    return Nsat * b * P / (1.0 + b * P)


def langmuir_fit(P, q):
    """Return dict with Nsat, b (1/bar), their 95 percent confidence intervals, and R^2."""
    P, q = np.asarray(P, float), np.asarray(q, float)
    p0 = [1.2 * q.max(), 1.0 / np.median(P)]
    popt, pcov = curve_fit(langmuir, P, q, p0=p0, bounds=(0, np.inf), maxfev=20000)
    t = stats.t.ppf(0.975, max(len(P) - 2, 1))
    ci = t * np.sqrt(np.diag(pcov))
    res = q - langmuir(P, *popt)
    r2 = 1.0 - np.sum(res**2) / np.sum((q - q.mean())**2)
    return {"Nsat": popt[0], "Nsat_CI95": ci[0], "b": popt[1], "b_CI95": ci[1], "R2": r2, "n": len(P)}


def vant_hoff(T, b):
    """Slope of ln b versus 1/T gives -dH/R. Returns dH (kJ/mol), its 95 percent CI, and R^2."""
    x, y = 1.0 / np.asarray(T, float), np.log(np.asarray(b, float))
    fit = stats.linregress(x, y)
    t = stats.t.ppf(0.975, max(len(x) - 2, 1))
    return {"dH": -R * fit.slope, "dH_CI95": R * t * fit.stderr, "R2": fit.rvalue**2}


def gibbs_entropy(T, b, dH, P0=10.0):
    """dG = -RT ln(b P0) and dS = (dH - dG)/T; returns dG (kJ/mol) and dS (J/(mol K))."""
    T, b = np.asarray(T, float), np.asarray(b, float)
    dG = -R * T * np.log(b * P0)
    dS = 1e3 * (dH - dG) / T
    return dG, dS


def qst_clausius(isotherms, loadings):
    """Isosteric heat Qst = -R d(ln P)/d(1/T) at constant loading.

    isotherms : dict {T: (P, q)} with q increasing in P
    loadings  : loadings (same units as q) at which Qst is evaluated
    Pressures at each loading are obtained by linear interpolation of ln P versus q.
    """
    Ts = sorted(isotherms)
    out = []
    for qi in loadings:
        lnP = []
        for T in Ts:
            P, q = map(np.asarray, isotherms[T])
            order = np.argsort(q)
            lnP.append(np.interp(qi, q[order], np.log(P[order]), left=np.nan, right=np.nan))
        lnP = np.array(lnP)
        ok = np.isfinite(lnP)
        if ok.sum() < 2:
            out.append(np.nan)
            continue
        slope = stats.linregress(1.0 / np.array(Ts)[ok], lnP[ok]).slope
        out.append(-R * slope)
    return np.array(out)


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "../data/computed/gcmc_dry_full_data.csv"
    d = pd.read_csv(path)
    for (system, comp), g in d.groupby(["system", "component"]):
        if system != comp:
            continue  # pure components
        rows = []
        for T, gt in g.groupby("T"):
            f = langmuir_fit(gt["P_bar"], gt["uptake"])
            rows.append((T, f))
        T = np.array([r[0] for r in rows])
        b = np.array([r[1]["b"] for r in rows])
        vh = vant_hoff(T, b)
        dG, dS = gibbs_entropy(T, b, vh["dH"])
        print(f"\n{system} ({comp}): dH = {vh['dH']:.2f} +/- {vh['dH_CI95']:.2f} kJ/mol")
        for (Ti, f), g_, s_ in zip(rows, dG, dS):
            print(f"  T = {Ti:.0f} K  Nsat = {f['Nsat']:.3f} +/- {f['Nsat_CI95']:.3f}  "
                  f"b = {f['b']:.4f} +/- {f['b_CI95']:.4f} 1/bar  R2 = {f['R2']:.4f}  "
                  f"dG = {g_:.2f} kJ/mol  dS = {s_:.1f} J/(mol K)")
