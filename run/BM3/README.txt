Axisymmetric cooling and spreading (BM3 from [1])
======================

The first test simulates a hot viscous fluid injected into a horizontal 
plane, which starts spreading and cooling in an axisymmetric fashion.
For this test we consider a Newtonian fluid injected onto the plane from 
a hole located on the lower wall.
The fluid is modelled as a hot Newtonian fluid that exchanges heat flux with
the ground through conduction and with the air through radiation and convection.

This test refers to an analogue experiment reported in [2], denoted as C14, 
where a hot silicone oil (Rhodorsil® 47V 5000 or 47V 12500, dyed red) is 
injected, at a constant supply rate R = 2.2 * 10^-8 m^3/s, onto a horizontal 
plane of polystyrene from a point source of 2-4 mm of radius.

Executing:

./Allrun

----------------------
Execution SUMMARY:

- blockMesh generates a rectangular block.
- snappyHexMesh refines the mesh around the centre of the lower wall.
- topoSet and createPatch create the hole at the centre of the lower wall.
- The solver adopted for this test case in interThermalRadConvFoam, a new
variation of the classic interFoam.
- The presence of the dictionary dynamicMeshDict inside the "constant" folder
produce the dynamic mesh refinement during the runtime. The refinement is 
performed at the interface between the two fluids, namely where the value of 
alpha is between 0.01 and 0.09.

----------------------
The code execution can be performed also in parallel:

./Allrun-parallel

======================
References:

[1] B. Cordonnier, E. Lev, and F. Garel, “Benchmarking lava-flow models,” 
Geological Society, London, Special Publications, vol. 426, pp. 425–445, 2015.

[2] F. Garel, E. Kaminski, S. Tait, and A. Limare, “An experimental study 
of the surface thermal signature of hot subaerial isoviscous gravity currents: 
Implications for thermal monitoring of lava flows and domes,” Journal of 
Geophysical Research: Solid Earth, vol. 117, 2012.
