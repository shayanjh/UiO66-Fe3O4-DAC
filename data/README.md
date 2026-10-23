# Data

Uptakes are in mmol/g (per gram of UiO-66 for the composites), pressures in bar, temperatures in K, and energies in kJ/mol (entropies in J mol⁻¹ K⁻¹), unless stated otherwise.

## computed/

| File | Content |
|---|---|
| `gcmc_dry_full_data.csv` | computed uptakes of pure CO₂, pure N₂, and of CO₂ and N₂ in their mixture (CO₂ mole fraction 0.15) in UiO-66 at every pressure and temperature, with the Langmuir and thermodynamic quantities |
| `gcmc_dry_summary.csv` | Langmuir parameters (N_sat, b) with 95 percent confidence intervals, R², isosteric heat, and ΔG, ΔH, and ΔS for each system and temperature |
| `gcmc_CO2H2O_full_data.csv` | computed uptakes of CO₂ and H₂O in UiO-66 at relative humidities of 0 to 1 |
| `gcmc_CO2H2O_summary.csv` | Langmuir and thermodynamic parameters under humidity (fits of the total uptake) |
| `composite_adsorption_data.csv` | computed CO₂ isotherms of UiO-66 and of Composites 1, 2, and 3 (50, 20, and 15 wt % Fe₃O₄) at 303 K |
| `composite_fits.csv` | Langmuir fits of the computed and measured isotherms of UiO-66 and the composites |
| `desorption_summary.csv` | desorption metrics of Composite 3 (number of CO₂ molecules, amount released, t50, t90, and peak temperature gradient) |
| `desorption_table.csv` | Fe₃O₄-localized versus uniform heating of Composite 3 (SI, Table S11) |

Main column names: `P_bar` pressure; `T` temperature; `uptake`, `uptake_std` mean uptake and its standard deviation over the production run; `Nsat`, `b_bar` Langmuir saturation capacity and affinity constant (1/bar); `R2_langmuir` coefficient of determination of the fit; `qst_kJ` isosteric heat; `dG_kJ`, `dH_kJ`, `dS_J` Gibbs free energy, enthalpy, and entropy of adsorption (standard pressure 10 bar); `rh` relative humidity.

## measured/

| File | Content |
|---|---|
| `UiO-66-isotherm.xlsx` | CO₂ isotherms of UiO-66 at 30, 50, and 70 °C |
| `Fe15@UiO-66-isotherms.xlsx` | CO₂ isotherms of Fe₃O₄@UiO-66 (15 wt % Fe₃O₄) at 30, 50, and 70 °C |
| `breakthrough_temperature_SAR_raw.xlsx` | breakthrough curves for MISA and TSA, under humid conditions, for the empty bed, and for CO₂/N₂ mixtures; bed temperatures during MISA and TSA; and SAR heating curves |
