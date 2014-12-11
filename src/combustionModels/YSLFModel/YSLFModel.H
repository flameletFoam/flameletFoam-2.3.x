/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModels::YSLFModel

Description
    Base class for YSLFModel.

SourceFiles
    YSLFModel.C

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef YSLFModel_H
#define YSLFModel_H


#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "tableSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class YSLFModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class YSLFModel
:
    public CombThermoType
{

protected:

    // Protected data

        // tables
	    tableSolver solver_;

        //- Species mass fractions
        PtrList<volScalarField>& Y_;

        //- Enthalpy
        volScalarField& he_;    

        //- Mixture fraction
        const volScalarField& Z_;

        //- Mixture fraction variance
        const volScalarField& varZ_;

        //- scalar dissipation rate
        const volScalarField& Chi_;

        //- Upper Bounds for table interpolation
        List<List<int> > ubIF_, ubP_;

        //- Table Position
        List<scalarList> posIF_, posP_;

private:

    // Private Member

        //- use scalar dissipation for interpolation switch
        Switch useScalarDissipation_;

        //- use mixture fraction variance for interpolation switch
        Switch useMixtureFractionVariance_;

        //- Disallow copy construct
	    YSLFModel(const YSLFModel&);

        //- Disallow default bitwise assignment
        void operator=(const YSLFModel&);

public:

    // Constructors

        //- Runtime type information
        TypeName("YSLFModel");

        hashedWordList tables();

        //- Construct from components
        YSLFModel(const word& modelType, const fvMesh& mesh);

        //- Destructor
        virtual ~YSLFModel();

        //- Correct combustion rate
        virtual void correct();

        //- Correct density during PISO
        Switch correctDensity();

        //- Fuel consumption rate matrix  
        inline tmp<fvScalarMatrix> R(volScalarField& Y) const;         

        //- Enthalpy source term
        inline tmp<Foam::volScalarField> Sh() const;

        //- Heat release rate calculated from fuel consumption rate matrix
        inline tmp<volScalarField> dQ() const;

        // I-O
        //- Update properties
        virtual bool read();


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "YSLFModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
