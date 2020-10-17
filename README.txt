EAMPA Program
#####################################################
Embedded Atom Method Potential Analyser

Written in Python and Fortran (uses F2PY).


What does it do?
#####################################################

Can be used to calculate:
· energy
· energy, forces 
· energy, forces, stress
of a collection of atoms

Calculates some other properties:
· a0
· v0
· e0
· b0
· 9 independent elastic constants (Orthorhombic)
· shear modulus
· young's modulus
· estimated melting temperature



What configuration files can be used?
#####################################################

There is a native format (see examples) similar to an xyz type file.
It will read Quantum Espresso PWscf output files.




What potentials?
#####################################################

EAM potentials.

Works with potentials in the following form:


Pair:

V(r)

Density - one or multiple distinct functions:

rho_1(r), rho_2(r), rho_3(r)

Embedding term - one or multiple terms:

F_1(rho_1) + F_2(rho_2) + F_3(rho_3)


The potentials can be tabulated or analytic.  It's easy to add new, custom functions, either in Python or the F2PY module.





RSS Calculation
#####################################################

EFS

Four weights: 
1. config weight - multiplier of overall rss
2. energy weight - multiplier of energy rss
3. force weight - multiplier of force rss
4. stress weight - multiplier of stress rss

How each are calculated
1. energy rss - rss between overall energy of calculated config vs known added up over all configs
2. force rss - rss for all forces, all 3 directions added up over all configurations






Still To Do/Fix
#####################################################
Only calculate BP if a potential exists
Only calculate BP rss if it's in the BP data file and has a potential
Surface energy






