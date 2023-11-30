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

Class
    Foam::addTRadConvImmiscibleIncompressibleTwoPhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    addTRadConvImmiscibleIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef addTRadConvImmiscibleIncompressibleTwoPhaseMixture_H
#define addTRadConvImmiscibleIncompressibleTwoPhaseMixture_H

#include "addTIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "fvMatricesFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class addTRadConvImmiscibleIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class addTRadConvImmiscibleIncompressibleTwoPhaseMixture
:
    public addTIncompressibleTwoPhaseMixture,
    public interfaceProperties
{

public:
    // EB: add Rad, add Conv
    // Private data
        const dimensionedScalar emissivity_;
        const dimensionedScalar sigma_SB_ ;
        const dimensionedScalar convectiveHeatTransferParameter_; 
        const dimensionedScalar fractionalAreaExposed_ ;
        const dimensionedScalar T_env_;        
    // end add Rad, add Conv

    // Constructors

        //- Construct from components
        addTRadConvImmiscibleIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~addTRadConvImmiscibleIncompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            addTIncompressibleTwoPhaseMixture::correct();
            interfaceProperties::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
        
        // EB: add Radiation
        //- Calculates radiative source term, returned as volum fraction rate
        tmp<fvScalarMatrix> calcSourceRadiation(
				  const volVectorField& U,
				  volScalarField& T ,
				  volScalarField& RadiativeCoeff 
		);
		//- Calculates the difference between the fourth powers 
		// of temperature and environmental temperature
        void calcT4mTenv4(
				  const volScalarField& T ,
				  volScalarField& T4mTenv4 
		);
	
        // EB: add Forced Convection
        //- Calculates convective source term, returned as volume fraction rate
        tmp<fvScalarMatrix> calcSourceForcedConvection(
				  const volVectorField& U,
				  volScalarField& T 
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
