# MASC


See [A&A, 686, A243 (2024)](https://doi.org/10.1051/0004-6361/202348967).

Mock AMR Sunyaev-Zeldovich Calculator (MASC) is a package written in Fortran 90 and OpenMP. It produces Sunyavev-Zeldovich (both kinetic and thermal) maps for MASCLET simulations processed with the halo finder [ASOHF](https://github.com/dvallesp/ASOHF).

In order to work, MASC requires the following folders inside the masc.f90 directory: `output_files` for saving the maps, `ASOHF_results` to read the dark matter haloes (clusters) detected by ASOHF and `simu_masclet` to read the MASCLET simulations outputs. An example `compile` file is also provided.
