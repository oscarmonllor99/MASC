# MASC
* Mock AMR Sunyaev-Zeldovich Calculator (MASC). Fortran and OMP based package to produce Sunyavev-Zeldovich (both kinetic and thermal) maps for MASCLET simulations processed with the halo finder ASOHF (see https://github.com/dvallesp/ASOHF).
* In order to work, MASC requires the following folders in the masc.f90 directory: "output_files" for saving the maps, "ASOHF_results" to read the dark matter haloes (clusters) detected by ASOHF and "simu_masclet" to read the MASCLET simulations outputs.
