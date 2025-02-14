surfaceScalarField kappaf = mixture.kappaf();

// EB: add 
volScalarField NoInterface ( mixture.nearInterface() );
forAll ( mesh.C(), celli )
{
	if ( NoInterface[celli] > 0)  NoInterface[celli] = 0;
	else  
	    NoInterface[celli] = 1;
}
// End add 

fvScalarMatrix TEqn
(	
	fvm::ddt(rhoCp,T)
	+ fvm::div(rhoCpPhi,T) 
	- NoInterface * fvm::laplacian(kappaf,T) // EB: diffusion modified
	// EB: add Source Terms
	== 	
	mixture.calcSourceRadiation(U,T,RadiativeCoeff) // EB: Radiative term
  	+ mixture.calcSourceForcedConvection(U,T) // EB: ForcedConvction term
	// End add Source Terms 	
);

// EB: add Rad
mixture.calcT4mTenv4(T,T4mTenv4);
// End add Rad

TEqn.relax();
TEqn.solve(); 
 
