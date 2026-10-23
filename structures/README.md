# Structure files

| File | Content |
|---|---|
| `data.UiO-66.ddec` | UiO-66 cell (3,648 atoms) with UFF parameters and DDEC charges |
| `Fe3O4_min.data` | Fe₃O₄ cell (ClayFF-based model of Gürsoy et al.) |
| `Fe3O4_slab.data` | relaxed Fe₃O₄ slab |
| `Composite3_UiO66_Fe3O4.data.gz` | equilibrated Composite 3 (compressed; decompress with `gunzip Composite3_UiO66_Fe3O4.data.gz`); 15 wt % Fe₃O₄, 75,158 atoms, neutral; types 1 C_R, 2 H_, 3 O_R, 4 O_3_z, 5 Zr, 6 Fe(A), 7 Fe(B), 8 O, and 9 C, 10 O reserved for CO₂; can be used directly in `lammps/adsorption/in.gcmc_composite_CO2` |
| `CO2.txt` | molecule file for CO₂ (TraPPE, rigid) |
| `N2.txt` | molecule file for N₂ (TraPPE, rigid) |
| `H2O.txt` | molecule file for H₂O (SPC/E, rigid) |
