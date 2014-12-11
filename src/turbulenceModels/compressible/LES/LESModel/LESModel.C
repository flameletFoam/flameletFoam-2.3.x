/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#include "LESModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LESModel, 0);
defineRunTimeSelectionTable(LESModel, dictionary);
addToRunTimeSelectionTable(turbulenceModel, LESModel, turbulenceModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void LESModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

LESModel::LESModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    fluidThermo& thermoPhysicalModel,           
    const word& turbulenceModelName
)
:
    turbulenceModel(rho, U, phi, thermoPhysicalModel, turbulenceModelName),

    IOdictionary
    (
        IOobject
        (
            "LESProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    transportVarZ_(lookupOrDefault<Switch>("transportVarZ", false)),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),
    CvarZ_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CvarZ",
            coeffDict_,
            1.0
        )
    ),
    Cchi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cchi",
            coeffDict_,
            2.0
        )
    ),
    Sc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sc",
            coeffDict_,
            1.0
        )
    ),
    Sct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sct",
            coeffDict_,
            1.0
        )
    ),
    runTime_(U.time()),
    kMin_("kMin", sqr(dimVelocity), SMALL),
    delta_(LESdelta::New("delta", U.mesh(), *this)),
    varZ_(thermoPhysicalModel.varZ()),       
    Chi_(thermoPhysicalModel.Chi())
{  
    kMin_.readIfPresent(*this);

    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<LESModel> LESModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    fluidThermo& thermoPhysicalModel,    
    const word& turbulenceModelName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "LESProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("LESModel")
    );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LESModel::New"
            "("
                "const volScalarField&, "
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "const fluidThermo&, "
                "const word&"
            ")"
        )   << "Unknown LESModel type "
            << modelType << nl << nl
            << "Valid LESModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<LESModel>
    (
        cstrIter()(rho, U, phi, thermoPhysicalModel, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LESModel::correct(const tmp<volTensorField>&)
{
    turbulenceModel::correct();
    delta_().correct();
}


void LESModel::correctVarZ()
{
	if (transportVarZ_)
	{
		tmp<fvScalarMatrix> varZEqn
	    (
	        (
	        	  fvm::ddt(rho_, varZ_)
	            + fvm::div(phi_, varZ_)
	            - fvm::laplacian(DZEff(), varZ_)
	            - 2.0*DZEff()*magSqr(fvc::grad(this->thermo().Z()))
	            + 2.0*rho_*Chi_
	        )
	    );

	    varZEqn().relax();
	    varZEqn().boundaryManipulate(varZ_.boundaryField());

	    solve(varZEqn);
	    varZ_.correctBoundaryConditions();
	    bound(varZ_, 0.0);
	}
	else
	{
	    varZ_ = CvarZ_ * sqr(delta()) * magSqr(fvc::grad(this->thermo().Z()));
	    varZ_.correctBoundaryConditions();
	}

    Info<< "varZ min/max   = " << min(varZ_).value() << ", "
        << max(varZ_).value() << endl;
}

void LESModel::correctChi()
{
	if (transportVarZ_)
	{
	    Chi_ = mu()/(Sc_*rho_)*magSqr(fvc::grad(this->thermo().Z())) + Cchi_*DZSgs()*varZ_/(2.0*rho_*sqr(delta()));
	    Chi_.correctBoundaryConditions();
	}
	else
	{
	    Chi_ = Cchi_ * DZEff()/rho_ * magSqr(fvc::grad(this->thermo().Z()));
	    Chi_.correctBoundaryConditions();
	}

    Info<< "chi min/max   = " << min(Chi_).value() << ", "
        << max(Chi_).value() << endl;
}


void LESModel::correct()
{
    correct(fvc::grad(U_));
}


bool LESModel::read()
{
    // Bit of trickery : we are both IOdictionary ('RASProperties') and
    // an regIOobject (from the turbulenceModel). Problem is to distinguish
    // between the two - we only want to reread the IOdictionary.

    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            IOdictionary::type()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        kMin_.readIfPresent(*this);

        delta_().read(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
