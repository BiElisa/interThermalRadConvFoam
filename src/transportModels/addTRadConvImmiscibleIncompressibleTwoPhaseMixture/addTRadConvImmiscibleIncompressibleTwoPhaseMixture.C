/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addTRadConvImmiscibleIncompressibleTwoPhaseMixture.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture::
addTRadConvImmiscibleIncompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    addTIncompressibleTwoPhaseMixture(U, phi),
    interfaceProperties(alpha1(), U, *this),
    
    // EB: add Rad, add Conv
    emissivity_
    (
        "emissivity", dimensionSet(0,0,0,0,0,0,0), lookup("emissivity")
    ) ,
    sigma_SB_
    (
        "sigma_SB", dimensionSet(1,0,-3,-4,0,0,0), lookup("sigma_SB")
    ) ,
    convectiveHeatTransferParameter_
    (
        "convectiveHeatTransferParameter", dimensionSet(1,0,-3,-1,0,0,0), lookup("convectiveHeatTransferParameter") // kg m s K mol A cd
    ) ,
    fractionalAreaExposed_
    (
        "fractionalAreaExposed", dimensionSet(0,0,0,0,0,0,0), lookup("fractionalAreaExposed")
    ) ,
    T_env_
    (
        "T_env", dimensionSet(0,0,0,1,0,0,0), lookup("T_env")
    )
    // end add Rad, add Conv
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture::read()
{
    return
        addTIncompressibleTwoPhaseMixture::read()
     && interfaceProperties::read();
}

// EB: add Rad
Foam::tmp<Foam::fvScalarMatrix> Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture::calcSourceRadiation(
const volVectorField& U,
volScalarField& T ,
volScalarField& RadiativeCoeff 
)
{
    Info << "calcSourceRadiation" << endl;
    // Initialize the radiative coefficient epsilon 
    // (not exactly as the equation (8) from CostaM. 2005)
    // as 1 or 0 for on interface or not
    volScalarField epsilonMatrix(nearInterface());
    // Initialize the field for the area of the interface surface
    // as 1 or 0 for on interface or not
    volScalarField Afs(epsilonMatrix);

    // Estimate cell surface area as V^(2/3), to be calculated with isosurfaces
    // later
    forAll(Afs, cellI) {
		if ( Afs[cellI] > 0 )
    	    Afs[cellI] = pow(U.mesh().V()[cellI],0.67);
    }
    forAll(epsilonMatrix.boundaryField(), patchI) {
	    forAll(epsilonMatrix.boundaryField()[patchI], faceI ) {
	        epsilonMatrix.boundaryFieldRef()[patchI][faceI] = scalar(0.0) ;
	    }
    }
    // Calculate the volume fraction rate term
    forAll(epsilonMatrix, cellI) {
	    if((T[cellI] > T_env_.value())) {   
	        epsilonMatrix[cellI] = emissivity_.value() * fractionalAreaExposed_.value()
	                 * sigma_SB_.value() * Afs[cellI] / U.mesh().V()[cellI] ;
	    }else{
	        epsilonMatrix[cellI] = scalar(0.0) ; 
	    }
    }
    
    dimensionedScalar dimCorr("dimCorr",dimMass/(pow4(dimTemperature)*pow3(dimTime)*dimLength),1);
    RadiativeCoeff = epsilonMatrix * dimCorr ;
/*    return(
        epsilonMatrix * dimCorr * pow4(T_env_) 
        - Foam::fvm::Sp( epsilonMatrix * dimCorr * pow3(T) , T)   
    );  
*/    return(
        epsilonMatrix * dimCorr * pow4(T_env_)
        - Foam::fvm::Sp( 4 * epsilonMatrix * dimCorr * pow3(T) , T)
        + epsilonMatrix * dimCorr * 3 * pow4(T) 
    ) ;    

}

void Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture::calcT4mTenv4(
	const volScalarField& T ,
	volScalarField& T4mTenv4 
	)
	{
		Info << "calcT4mTenv4" << endl;
		forAll(T,cellI){
			T4mTenv4[cellI] = pow4(T[cellI]) - pow4(T_env_.value());
		}
}

// End add Rad

// EB: add Conv
Foam::tmp<Foam::fvScalarMatrix> Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture::calcSourceForcedConvection(
	const volVectorField& U,
	volScalarField& T  
	)
	{
		Info << "calcSourceForcedConvection" << endl;
		// Initialize the forced convective coefficient W 
		// as 1 or 0 for on interface or not
		volScalarField WMatrix(nearInterface());
		// Initialize the field for the area of the interface surface
		// as 1 or 0 for on interface or not
		volScalarField Afs(WMatrix);

		// Estimate cell surface area as V^(2/3), to be calculated with isosurfaces
		// later
		forAll(Afs, cellI) {
			if ( Afs[cellI] > 0 )
				Afs[cellI] = pow(U.mesh().V()[cellI],0.67);
		}
		forAll(WMatrix.boundaryField(), patchI) {
			forAll(WMatrix.boundaryField()[patchI], faceI ) {
				WMatrix.boundaryFieldRef()[patchI][faceI] = scalar(0.0) ;
			}
		}
		// Calculate the volume fraction rate term
		forAll(WMatrix, cellI) {
			if((T[cellI] > T_env_.value())) {   
				WMatrix[cellI] = convectiveHeatTransferParameter_.value() * fractionalAreaExposed_.value()
						  * Afs[cellI] / U.mesh().V()[cellI] ;
			}else{
				WMatrix[cellI] = scalar(0.0) ; 
			}
		}
		
		dimensionedScalar dimCorr("dimCorr",dimMass/(dimTemperature*pow3(dimTime)*dimLength),1);
	 //   ConvectiveCoeff = WMatrix * dimCorr ; 
	 /*  activate this in case you create the volScalarField ConvectiveCoeff to visualize it on ParaFoam */
		return(
			WMatrix * dimCorr * T_env_
			- Foam::fvm::Sp( WMatrix * dimCorr , T)
		) ;    

}
// End add Conv


// ************************************************************************* //
