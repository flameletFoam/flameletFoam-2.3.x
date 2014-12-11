Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::rhoCombustionModel> combustion
(
    combustionModels::rhoCombustionModel::New
    (
        mesh
    )
);

combustion->correct();

rhoReactionThermo& thermo = combustion->thermo(); 

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

// Mixture Fraction
volScalarField& Z = thermo.Z();

#include "compressibleCreatePhi.H"

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

const volScalarField& psi = thermo.psi();
