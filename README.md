# CE-Espresso

This is a customized version of [Quantum-Espresso](https://www.quantum-espresso.org). It includes few changes and features that are either experimental or in conflict with the official QE code.
I try to keep the code base as much as possible aligned to and compatible with the official release.

### Disclaimer
No warrany, use at your own risk!


## Summary of features
### PW
* March 2017: **Support for restarting and computing band structure with metaGGAs (e.g. TB09, TPSS)**
* support for Rydberg pseudopotentials, having larger radial mesh and large cutoff
* always force positive charge density
* `conv_thr_init` is enforced also from the 2nd SCF iteration on
* complete list of atoms in `hubbard_l` and `tabd`
* Ewald parameter computed consistently in all routines
* possibility to control what to mix: density, magnetization, Hubbard occupations, PAW projectors
* possibility to apply an external potential, read from file
* computation of _Z2_ invariants according to PRB 83, 235401 (2011), working only with NCPPs

### atomic
* possibility to generate old UPFv1 pseudopotentials
* shift energy of individual angular momentum channels
* in AE calculations, output the spin density on the nucleus
* output the _l_-dependent pseudo-potentials

### PP
* plot the kinetic energy density $\tau$ (`plot_num = 21`)


## Future plans
* multi-orbital Hubbard U (will break input format)
* generalize _Z2_ to ultrasoft and PAWs
* implement projectability according to PRB 88, 165127 (2013)




