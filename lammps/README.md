# LAMMPS input scripts

The input scripts follow the simulation protocol of Section 3 of the article and Sections S2 to S8 of the Supporting Information. The force field, the MD–GCMC logic, and the heating protocols are given in full. Everything that depends on the particular system is collected at the top of each script in a USER INPUTS block, as variables whose values are placeholders such as `{datafile}` or `{P_list}`. Each placeholder must be replaced, either in the script or on the command line, for example

```
lmp -var datafile my.data -var P_list "0.1 0.5 1.0" -var f_list "0.0987 0.4925 0.9832" ... -in in.gcmc_CO2
```

The molecule files `CO2.txt`, `N2.txt`, and `H2O.txt` are in `structures/` and must be copied to the folder in which LAMMPS is run.

## Contents

| File | Purpose |
|---|---|
| `adsorption/in.gcmc_CO2` | isotherm of pure CO₂ in UiO-66 |
| `adsorption/in.gcmc_N2` | isotherm of pure N₂ in UiO-66 |
| `adsorption/in.gcmc_CO2_N2` | coadsorption from a CO₂/N₂ mixture in UiO-66 |
| `adsorption/in.gcmc_CO2_H2O` | coadsorption of CO₂ and H₂O in UiO-66 at a given relative humidity |
| `adsorption/in.gcmc_composite_CO2` | isotherm of CO₂ in a Fe₃O₄/UiO-66/Fe₃O₄ composite |
| `composite/in.UiO66_relax` | relaxation of the UiO-66 cell and replication along z |
| `composite/in.Fe3O4_slab` | replication and relaxation of the Fe₃O₄ slab |
| `composite/build_composite.py` | cutting and H capping of the UiO-66 slab, and placement of the Fe₃O₄ slabs |
| `composite/in.composite_equilibration` | equilibration of the composite at 303 K |
| `desorption/in.desorption_magnetic` | magnetic desorption: heating localized in the Fe₃O₄ slabs |
| `desorption/in.desorption_uniform` | uniform desorption (control): heating of the whole solid |
| `cross_terms_Fe3O4_UiO66.lmp` | Fe₃O₄–framework Lennard-Jones cross terms |

The construction of the composite is described step by step in `composite/README.md`.

## Force field

* Framework: UFF with DDEC charges. Fe₃O₄: ClayFF-based model of Gürsoy et al. Guests: rigid TraPPE CO₂ and N₂, and rigid SPC/E water.
* Lennard-Jones and real-space Coulomb cutoff 12.5 Å; Lorentz–Berthelot mixing for unlike pairs.
* Fe₃O₄–framework cross terms (`cross_terms_Fe3O4_UiO66.lmp`): Lorentz–Berthelot combinations of the UFF parameters with the Fe (ε = 0.30 kcal/mol, σ = 1.853 Å) and O (ε = 0.20 kcal/mol, σ = 3.127 Å) parameters of Kerisit et al. (SI, Table S4). The file is included after `read_data` in every script that contains Fe₃O₄.
* Ewald summation with an accuracy of 10⁻⁴; for the composites, which are nonperiodic in z, the slab correction of Yeh and Berkowitz.

## Adsorption

* GCMC with insertion and deletion moves only; the fugacity of each species is passed directly as the reservoir pressure, with `fugacity_coeff 1`.
* `P_list` holds the pressures in bar and is used only to label the output files. `f_list` (or `fC_list` and `fN_list` for the mixture) holds the fugacities in atm at the same pressures and temperature, computed with `analysis/fugacity_mbwr.ipynb`.
* At every pressure, each script writes `raw_P=<P>.txt` with the step, the pressure, the number of adsorbed molecules N (and of each species for mixtures), N², the potential energy U, and U·N. `analysis/adsorption_analysis.py` converts these files into uptakes.
* The humid script starts from a framework equilibrated beforehand with water alone at the partial pressure RH × P_sat(T) (2 × 10⁶ MD steps; Section 3.4 of the article).
* In the composite, CO₂ insertions are restricted to the interior of the MOF slab.

## Desorption

Both scripts start from the CO₂-loaded composite at the end of the composite adsorption run at 1 bar, keep the exchange with a CO₂ reservoir at 1 bar, restrict insertions to the MOF, and use the current mean MOF temperature as the acceptance temperature. They write the time series of the MOF and oxide temperatures and of the number of CO₂ molecules (`desorption_timeseries.txt`) and the temperature profile T(z) in 1.5 Å bins (`Tz_profile.txt`), which are analyzed with `analysis/desorption_metrics.py`. The heating protocols are those of Table S10 of the SI.
