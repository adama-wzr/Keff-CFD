# Keff-CFD
Effective thermal conductivity estimation in binary discrete 2D and 3D systems.

Current model ignores radiation and convection, only first order thermal conduction effects are considered. The CFD algorithm calculates the steady-state temperature profile of the given input structure, and from the temperature map the heat fluxes can be calculated and the effective property derived. More information in the documentation pdf's.

This code is under development. The 2D code is functional within the bounds defined in the documentation, there are versions of functional 3D code, albeait at the moment lagging behind the 2D code in terms of robustness and ease of use. The current code is being developped using both C and C++ libraries, but must be compiled with a C++ compiler.
