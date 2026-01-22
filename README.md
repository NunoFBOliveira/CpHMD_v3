# CpHMD_v3
Git commit of the improved stochastic-CpHMD base code.

The original methodology can be found in :
[1] Baptista et al., J. Chem. Phys. 117, 4184 (2002) DOI: 10.1063/1.1497164

This new version of CpHMD enjoys several upgrades from the previous st-CpHMD used in the Machuqueiro Computational BiopHysics Lab in Lisbon.

Some of the most notable are:
   - Verified support for Recent GROMACS codes.
   - Verified support for GPU usage in the MD steps.
   - Support for multiple force-fields, including:
       * GROMOS: G54a7pH - Base parameterization of the methodology;
       * CHARMM: CHARMM36pH - New parameterization for titratable amino acids in CHARMM (https://doi.org/10.1021/acs.jpcb.2c04529)
       * AMBER: AMBER14SBpH - New parameterization for titratable amino acids in AMBER (https://doi.org/10.1021/acs.jctc.5c00415)
   - Support for Reduced Titration on sites with a predicted pKa far from the simulation pH.
   - Support for plumed metadynamic simulations.
   

Inside this project you can find:
   - The base code (CpHMD_v3.4.2);
   - The DelphiTools with all the needed parameters and scripts for the PB calculation;
   - petit 1.6.1 needed for the Monte Carlo step, developed in ITQB (https://www.itqb.unl.pt/labs/molecular-simulation/in-house-software);
   - The force-field (FFs) parameterized for st-CpHMD;
   - The tautomer files (Sts) for the titrating molecules. 
