Corrections
===========
corrections.pdf -- list of corrections to published articles in scientific journals



Various python scripts implementations:


OL07.py
=======
Turbulence relative velocities

Reference:
    Ormel & Cuzzi (2007)

Calculation of particle's relative velocity in turbulent flows. The turbulent
power spectrum is assumed to be Kolmogorov.


OL18.py
=======
Pebble accretion efficiencies (small planets)

Reference:
    Ormel & Liu (2018)
    Liu & Ormel (2018)

Calculation of the pebble accretion probability for a single planet for
pebble-sized particles of Stokes number St<1. Expressions account for planet on
a Keplerian orbit and for a quite general turbulence model. Expressions have
been calibrated against three-body simulations.


HO23.py
=======
Pebble accretion probabilities for large pebbles

Reference:
    H. Huang & Ormel (2023)

Calculation of pebble accretion probability (epsilon) for a single planet in
the limit of aerodynamically large pebbles (St>1). The expressions apply only
for a planet on a circular orbit.


HO23i.py
========
Trapping of planets in first order (j+1:j) resonances

Reference:
    S. Huang & Ormel (2023)


def ta_crit (j, ...)
-----------
Calculation of the critical migration time of a planet, below which it crosses the resonance

