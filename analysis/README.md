# Analysis

| File | Content |
|---|---|
| `fugacity_mbwr.ipynb` | fugacities of CO₂, N₂, and their mixture from the modified Benedict–Webb–Rubin equation of state, written in a form that can be passed to `fix gcmc` |
| `adsorption_analysis.py` | uptakes from the `raw_P=<P>.txt` files of the LAMMPS input scripts, Langmuir fits with 95 percent confidence intervals, van't Hoff enthalpy, Gibbs free energy and entropy (standard pressure 10 bar), and the isosteric heat from the Clausius–Clapeyron equation |
| `desorption_metrics.py` | released fraction, t50, t90, and t360 from the time series of a desorption run |

Running `python adsorption_analysis.py` reproduces the Langmuir and thermodynamic parameters of pure CO₂ and pure N₂ listed in `data/computed/gcmc_dry_summary.csv` and in Table S7 of the SI.
