# Keff-CFD
Effective thermal conductivity estimation in binary discrete 2D and 3D systems.

Current model ignores radiation and convection, only first order thermal conduction effects are considered. The CFD algorithm calculates the steady-state temperature profile of the given input structure, and from the temperature map the heat fluxes can be calculated and the effective property derived. More information in the documentation pdf's.

This code is under development. The 2D code is functional within the bounds defined in the documentation, there are versions of functional 3D code, albeit due to constraints in computational power has not been as thoroughly tested as the 2D code. The current code is being developed using both C and C++ libraries, but must be compiled with a C++ compiler (g++ recommended).

If this code is useful in your research, please acknowledge that this GitHub repository was used and consider citing the following work:

https://www.sciencedirect.com/science/article/pii/S2666546823000824

Direct any questions via email to one of the authors or directly here on GitHub.
