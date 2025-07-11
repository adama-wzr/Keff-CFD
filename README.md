# Keff-CFD

**Updated README!**

Stagnant effective thermal conductivity estimation in binary discrete 2D and 3D systems.

Current model ignores radiation and convection, only first order thermal conduction effects are considered. The NHT algorithm calculates the steady-state temperature profile of the given input structure, and from the temperature map the heat fluxes can be calculated and the effective property derived. Be sure to check out the documentation!


# Table of Contents

1. [Documentation](#documentation)
2. [How to Cite](#publications)
3. [Authors](#authors)
4. [Acknowledgements](#acknowledgements)
5. [Upcoming Changes](#upcoming-changes)

## Documentation

All information regarding compiling and running can be found in the documentation pdf's.

This code is under development. The 2D code is functional within the bounds defined in the documentation, there are versions of functional 3D code, albeit due to constraints in computational power has not been as thoroughly tested as the 2D code. The current code is being developed using both C and C++ libraries, but must be compiled with a C++ compiler (g++ recommended).

## Publications

There isn't a consolidated publication for this package. If this is useful, please consider citing one of the publications below:

- **Adam, A.**, Fang, H., & Li, X. (2024). Effective thermal conductivity estimation using a convolutional neural network and its application in topology optimization. Energy and AI, 15, 100310. https://doi.org/10.1016/j.egyai.2023.100310

**Dataset Publications:**

- **Adam, Andre**; Li, Xianglin; Fang, Huazhen (2023), “2D Binary Images and Effective Thermal Conductivity CFD Results”, Mendeley Data, V2, doi: 10.17632/454dsrmdyf.2

## Authors

- Main developer: Andre Adam (The University of Kansas)
    - [ResearchGate](https://www.researchgate.net/profile/Andre-Adam-2)
    - [GoogleScholar](https://scholar.google.com/citations?hl=en&user=aP_rDkMAAAAJ)
    - [GitHub](https://github.com/adama-wzr)
- Advisor: Dr. Xianglin Li (Washingtion University in St. Louis)
    - [Website](https://xianglinli.wixsite.com/mysite)
    - [GoogleScholar](https://scholar.google.com/citations?user=8y0Vd8cAAAAJ&hl=en)

Direct any questions via email to one of the authors or directly in the discussion section here on GitHub.

## Acknowledgements

This work wouldn't be possible without the computational time awarded as part of the following grants:

This work was only possible with the Extreme Science and Engineering Discovery Environment (XSEDE) grants for computational time in the Bridges2 and Expanse supercomputers, under the awards MAT210007 and MAT210014.

This work used Expanse(GPU) at SDSC through allocations MAT230071 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program, which is supported by National Science Foundation grants #2138259, #2138286, #2138307, #2137603, and #2138296.

## Upcoming Changes

- Hey, changes are coming! As soon as I have a plan, I will make sure to update this page!
