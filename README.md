This repository contains the computer code (in Matlab and Singular) for all mathematical analysis, and for generating all figures, in the publication  "Ultrasensitivity and bistability in covalent-modification cycles with positive autoregulation", in Proceedings of the Royal Society A.

Authors: Cailan Jeynes-Smith and Robyn P. Araujo, July 2021, Queensland University of Technology, Brisbane, Australia.

The following scripts were used to generate figures for the above-noted publication. Specifically:
- Fig 3 uses MichaelisMentenKineticModel.m
- Figs 8,10,15,16 uses PARdNumericalSimulations.m
- Figs 1,6,7,8,9,10,11, uses PARdGroebnerRootFinder.m
- Figs 6,14,15,17 uses PARiNumericalSimulations.m
- Figs 7,12,13 uses PARiGroebnerRootFinder.m

These five Matlab scripts generate the plots in the associated figures, and will call other scripts as necessary. In order to generate the figures in the published article, the associated function is run using the parameters identified in the caption of the relevant figure. 

The contents of PARdSingularCode.txt and PARiSingularCode.txt can be run in Singular (https://www.singular.uni-kl.de) to obtain the equations and coefficients used by PARdGroebnerRootFinder.m and PARiGroebnerRootFinder.m, respectively, and their dependent functions.

For any further questions, please contact c2.jeynessmith@qut.edu.au. 


