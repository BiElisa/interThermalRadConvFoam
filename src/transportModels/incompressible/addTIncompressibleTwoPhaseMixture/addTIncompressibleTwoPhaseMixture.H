/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::addTIncompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    addTIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef addTIncompressibleTwoPhaseMixture_H
#define addTIncompressibleTwoPhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class addTIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class addTIncompressibleTwoPhaseMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
        // EB : Add T
	   dimensionedScalar cp1_; 
	   dimensionedScalar cp2_; 
	   dimensionedScalar Pr1_; 
	   dimensionedScalar Pr2_;
	   // End Add T
	   
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    TypeName("addTIncompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        addTIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~addTIncompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return *nuModel1_;
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return *nuModel2_;
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };
        
        // EB : Add T
       	//- Return const-access to phase1 heat capacity 
        	const dimensionedScalar& cp1() const 
        	{ 
           	return cp1_; 
        	} 

        	//- Return const-access to phase2 heat capacity  
        	const dimensionedScalar& cp2() const 
        	{ 
            	return cp2_; 
        	}; 

        	//- Return const-access to phase1 Prandtl number 
        	const dimensionedScalar& Pr1() const 
        	{ 
			
            	return Pr1_; 
        	} 
        // End Add T
        
        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the dynamic laminar viscosity on patch
        tmp<scalarField> mu(const label patchI) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        // Add T
	    //- Return the face-interpolated conductivity 
	    tmp<surfaceScalarField> kappaf() const; 
	    // end Add T

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
