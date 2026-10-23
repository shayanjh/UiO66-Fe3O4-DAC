# Building the Fe₃O₄/UiO-66/Fe₃O₄ composite

The composite is a UiO-66 slab, H-capped at its two faces, placed between two Fe₃O₄ slabs, with the system periodic in x and y and nonperiodic in z (Section 3.2 of the article and Section S3 of the SI). The Fe₃O₄ weight fraction is set by the thickness of the UiO-66 slab. The equilibrated Composite 3 (15 wt % Fe₃O₄) is given in `structures/Composite3_UiO66_Fe3O4.data.gz` (decompress it with `gunzip`); it already contains the CO₂ atom types and can be used directly in `../adsorption/in.gcmc_composite_CO2`, so that steps 1 to 5 are needed only to build other composites.

## Steps

1. **Relax and replicate UiO-66.**
   ```
   lmp -var datafile data.UiO-66.ddec -var nz {nz} -var outfile {UiO66_bulk.data} -var seed {seed} -in in.UiO66_relax
   ```
   The cell is relaxed (minimization with box relaxation, then NVT at 300 K) and replicated `nz` times along z.

2. **Prepare the Fe₃O₄ slab.**
   ```
   lmp -var datafile Fe3O4_min.data -var nx {nx} -var ny {ny} -var nz {nz} -var outfile {Fe3O4_slab.data} -var seed {seed} -in in.Fe3O4_slab
   ```
   Choose `nx` and `ny` so that the lateral box lengths of the oxide slab are as close as possible to those of the UiO-66 cell, and `nz` for the thickness of the oxide slab.

3. **Cut, cap, and stack.**
   ```
   python build_composite.py --mof {UiO66_bulk.data} --oxide {Fe3O4_slab.data} \
       --zlo {zlo} --zhi {zhi} --gap {gap} --vacuum {vacuum} --out {composite.data}
   ```
   * The slab is cut by building units: a Zr₆O₄(OH)₄ node is kept if all of its atoms lie between `zlo` and `zhi`, and a BDC linker is kept only if it lies in the same window and both nodes to which it is coordinated are kept, so that only Zr–O(carboxylate) bonds are broken.
   * Every Zr that lost a carboxylate oxygen is capped with a hydroxyl group placed along the broken bond (Zr–O 2.15 Å, O–H 0.97 Å, Zr–O–H 110°). The bonds, angles, and dihedrals of the caps are assigned by analogy with the bulk framework, and the cap charges are adjusted so that the slab is neutral.
   * The two Fe₃O₄ slabs are placed symmetrically below and above the UiO-66 slab, at a distance `gap` from its surfaces, with `vacuum` Å of empty space at each end of the box. The Fe₃O₄ types are appended after the five framework types (6 Fe(A), 7 Fe(B), 8 O).
   * The script stops if the lateral box lengths of the oxide and the MOF differ by more than `--lateral_tol` (2 percent by default), and prints an audit: numbers of atoms and caps, net charge, shortest MOF–oxide distance, and the Fe₃O₄ weight fraction. Check the audit before using the composite.

4. **Check the weight fraction.** Adjust `zlo` and `zhi` (and, if needed, `nz` in step 1) until the Fe₃O₄ weight fraction printed in step 3 is the desired one (50, 20, and 15 wt % for Composites 1, 2, and 3 of the article).

5. **Equilibrate.**
   ```
   lmp -var datafile {composite.data} -var outfile {composite_equil.data} -var seed {seed} -in in.composite_equilibration
   ```
   Minimization, a gentle start, heating from 50 to 303 K, and equilibration at 303 K. The Fe₃O₄–framework cross terms are read from `../cross_terms_Fe3O4_UiO66.lmp`, and atom types 9 and 10 (with bond type 8 and angle type 15) are reserved for CO₂, so that the equilibrated composite can be used directly in the adsorption script.

The equilibrated composite is then used in `../adsorption/in.gcmc_composite_CO2`, and the CO₂-loaded composite at 1 bar in `../desorption/in.desorption_magnetic` and `../desorption/in.desorption_uniform`.
