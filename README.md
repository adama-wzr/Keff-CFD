# Keff-CFD
Effective thermal conductivity estimation in binary discrete 2D and 3D systems.

Current model ignores radiation and convection, only first order thermal conduction effects are considered. The NHT algorithm calculates the steady-state temperature profile of the given input structure, and from the temperature map the heat fluxes can be calculated and the effective property derived.

All information regarding compiling and running can be found in the documentation pdf's.

This code is under development. The 2D code is functional within the bounds defined in the documentation, there are versions of functional 3D code, albeit due to constraints in computational power has not been as thoroughly tested as the 2D code. The current code is being developed using both C and C++ libraries, but must be compiled with a C++ compiler (g++ recommended).

If this package is useful in your research, please acknowledge that this GitHub repository was used and consider citing the following work:

Andre Adam, Huazhen Fang, Xianglin Li,
Effective thermal conductivity estimation using a convolutional neural network and its application in topology optimization,
Energy and AI,
Volume 15,
2024,
100310,
ISSN 2666-5468,
https://doi.org/10.1016/j.egyai.2023.100310.

Direct any questions via email to one of the authors or directly in the discussion section here on GitHub.
