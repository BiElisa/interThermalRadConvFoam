Pico do Fogo erption 2014-2015 (Cape Verde) - real lava flow test case
======================

Simulation of the early hours of the effusive eruption occurred 
at Pico do Fogo (Cape Verde) 2014-2015 [1,2]. 
This test case uses the real topography of Pico do Fogo and accounts for 
two vents as eruptive sources. 
Lava is modelled as a fluid with a temperature dependent viscosity that 
exchanges heat flux with the ground through conduction and with the air 
through radiation and convection.

Executing:

./Allrun

----------------------
Execution SUMMARY:

- blockMesh generates a rectangular block.
- surfaceFeatureExtract reads the stl file that contains the topography
of Pico do Fogo volcano and extrapolates the surface features.
- snappyHexMesh intersects the topography with the original mesh and creates 
the actual computational domain. Moreover, the mesh is refined close 
to the topography.
- topoSet and createPatch create the vents.
- The solver adopted for this test case in interThermalRadConvFoam, a new
variation of the classic interFoam.

----------------------
The code execution can be performed also in parallel:

./Allrun-parallel

======================
References:

[1] N. Richter, M. Favalli, E. de Zeeuw-van Dalfsen, A. Fornaciai, R. M. 
da Silva Fernandes, N. M. P erez, J. Levy, S. V. Silva, and T. R. Walter, 
“Lava flow hazard at Fogo Volcano, Cabo Verde, before and after the 2014–2015 
eruption,” Nat. Hazards Earth Syst. Sci., vol. 16, pp. 1925–1951, 2016.

[2] A. Cappello, G. Ganci, S. Calvari, N. M. P ́erez, P. A. Hern ́andez, S. V. 
Silva, and J. C. D. N. Cabral, “Lava flow hazard modeling during the 2014–2015 
Fogo eruption, Cape Verde,” J. Geophys. Res. Solid Earth, vol. 121, 
pp. 2290–2303, 2016.
