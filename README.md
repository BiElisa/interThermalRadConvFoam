# interThermalRadConvFoam
A new OpenFOAM solver, called interThermalRadConvFoam, to simulate free-surface viscous fluids with temperature changing due to radiative, convective, and conductive heat exchanges. The solver is accompained by two test cases, one is a benchmark and the other one is the simulation of a lava flow on a real topography.

-----------------------

Welcome to this project related to OpenFOAM development!
We are glad to present a new OpenFOAM solver, called interThermalRadConvFoam, 
to simulate free-surface viscous fluids with temperature changes due 
to radiative, convective, and conductive heat exchanges and the application 
of the solver to a benchmark test case and a lava flow simulation on a real 
topography. The solver is based on interFoam (available with OpenFOAM) and 
thus on the Volume of Fluid (VoF) technique used to describe the multiphase 
dynamics of two incompressible, viscous, and immiscible fluids (based on 
the Interface Capturing strategy). In our model, the two fluids are the 
fluid of interest with high viscosity and the surrounding atmosphere. 
Being interested in temperature effects, we added to the mass and momentum 
equations from interFoam, an equation for the energy that models the thermal 
exchanges between the fluid and the environment. Furthermore, a 
temperature-dependent viscoplastic model is used for the final application 
to lava flows. 

We present two test cases: 

(1) a benchmark from literature based on a laboratory experiment [1,2] 

(2) an application to a real lava flow by simulating the early hours of 
the 23 November 2014 eruption of the Fogo volcano at Cape Verde, West 
Africa [3,4]. 

For the simulation of the laboratory experiment, we also present simulations 
executed using a dynamic mesh with adaptive refinement.

## Installation
1. Execute from inside the folder "src":

     ./Allwmake

2. Execute from inside the folder "applications":

     ./Allwmake

* Alternatively, to install all at one time, execute from the main folder:

  ./Allwmake


## Run a test case
Move inside the folder of the test case of interest and execute:

./Allrun


## Compatibility: 
OpenFoam-v1912


## Who do I talk to? 
Elisa Biagioli (elisa.biagioli@manchester.ac.uk or biagioli.elisa88@gmail.com)


## References:

[1] B. Cordonnier, E. Lev, and F. Garel, “Benchmarking lava-flow models,” 
Geological Society, London, Special Publications, vol. 426, pp. 425–445, 2015.

[2] F. Garel, E. Kaminski, S. Tait, and A. Limare, “An experimental study 
of the surface thermal signature of hot subaerial isoviscous gravity currents: 
Implications for thermal monitoring of lava flows and domes,” Journal of 
Geophysical Research: Solid Earth, vol. 117, 2012.

[3] N. Richter, M. Favalli, E. de Zeeuw-van Dalfsen, A. Fornaciai, R. M. 
da Silva Fernandes, N. M. P erez, J. Levy, S. V. Silva, and T. R. Walter, 
“Lava flow hazard at Fogo Volcano, Cabo Verde, before and after the 2014–2015 
eruption,” Nat. Hazards Earth Syst. Sci., vol. 16, pp. 1925–1951, 2016.

[4] A. Cappello, G. Ganci, S. Calvari, N. M. P ́erez, P. A. Hern ́andez, S. V. 
Silva, and J. C. D. N. Cabral, “Lava flow hazard modeling during the 2014–2015 
Fogo eruption, Cape Verde,” J. Geophys. Res. Solid Earth, vol. 121, 
pp. 2290–2303, 2016.
