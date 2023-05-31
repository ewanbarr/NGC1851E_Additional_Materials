# NGC 1851 Software

This directory contains the software necessary to reproduce the results
presented in the paper. The contents of the directory are listed below.

**mass_mass_plot.c**
- Used to reproduce the Mass-Mass diagram presented in the supplementary material
- Has PGPLOT as dependency
- Can be complied with `gcc -o ../bin/mass_mass_plot -I/aux/miraculix/psrsoft/pgplot/ mass_mass_plot.c -L/aux/miraculix/packages/linux/pgplot -L/usr/X11R6/lib64 -lcpgplot -lpgplot -lX11 -lm -lgfortran`

**calc_probability_3rd_body.java**
- Used to estimate the omega-dot contributions from possible 3rd bodies in the NGC 1851E system
- Has a dependency on the org.spaceroots.mantissa mathematics package
- Can be run in an IDE or from the command line using `javac/java`

 **Timing_residuals_and_observation_dist.ipynb**
- Used to reproduce the figures for the post-fit timing residuals (as a function of epoch and orbital phase) and the orbital coverage of all the observations presented in the supplementary materials

**Evolution_of_periastron_advance.ipynb**
- Used to reproduce the figure for the evolution of the advance of periastron as a function of observations presented in the supplementary material

**PeriastronAdvance_LTPrecession.ipynb**
- Used to reproduce the figure for the maximum contribution to the rate of advance of periastron from Lense-Thirring precession presented in the supplementary material





