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
    Foam::combustionModels::tableSolver

Description
    Abstract class for the interface between table look-up and combustion model.

SourceFiles
    tableSolver.C

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "flameletTable.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "fvMesh.H"

namespace Foam
{
namespace combustionModels
{

class tableSolver
:
	public IOdictionary
{

protected:

	//- List of the used tables
	wordList tableNames_;

	//- List of the table parameters
	wordList paramNames_;

	//- Chemistry tables of all tabulated properties
	PtrList<flameletTable> tables_;

	//- List of table parameter distributions
	List<scalarList> params_;

	//- Scalar Dissipation Rate at extinction
	scalar chiExt_;

public:

	//- Runtime type information
    TypeName("tableSolver");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        tableSolver,
        dictionary,
        (
            const fvMesh &mesh
        ),
        (mesh)
    );

    //- Constructor
    tableSolver(const fvMesh& mesh, const wordList& tableNames);

    //- Destructor
    virtual ~tableSolver();

    //- Determine upper bounds for table interpolation
    List<int> upperBounds(const scalarList& x) const;

    //- Determine table position
    scalarList position(const List<int>& ub, const scalarList& x) const;

    //- Interpolate in lookup-table with mixture fraction, mixture fraction variance and scalar dissipation rate
    scalar interpolate(const List<int>& ub , const scalarList& pos, const label& i) const;

    //- Return the maximum scalar dissipation rate i.e. extinction
    scalar maxChi();

    int sizeTableNames() const;
};

} // End Foam namespace
} // End combustionModels namespace

#endif /* TABLESOLVER_H_ */
