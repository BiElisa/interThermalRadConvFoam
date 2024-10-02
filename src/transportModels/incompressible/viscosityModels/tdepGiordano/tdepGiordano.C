/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "tdepGiordano.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tdepGiordano, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        tdepGiordano,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tdepGiordano::calcNu() const
{
    
    const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T");
    
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            dimensionedScalar("n", dimensionSet(0,2,-1,0,0,0,0), 1.0)*(Foam::pow
			(
			10.0, 
			A_ + B_/max
			(
				(T+C_),
				dimensionedScalar("minT", dimensionSet(0,0,0,1,0,0,0), 100.0)
			))
			)/ 2700.0
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tdepGiordano::tdepGiordano
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tdepGiordanoCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    A_(tdepGiordanoCoeffs_.lookup("A")),
    B_(tdepGiordanoCoeffs_.lookup("B")),
    C_(tdepGiordanoCoeffs_.lookup("C")),
    nuMin_(tdepGiordanoCoeffs_.lookup("nuMin")),
    nuMax_(tdepGiordanoCoeffs_.lookup("nuMax")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::tdepGiordano::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    tdepGiordanoCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    tdepGiordanoCoeffs_.lookup("A") >> A_;
    tdepGiordanoCoeffs_.lookup("B") >> B_;
    tdepGiordanoCoeffs_.lookup("C") >> C_;
    tdepGiordanoCoeffs_.lookup("nuMin") >> nuMin_;
    tdepGiordanoCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
