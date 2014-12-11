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
    Foam::combustionModels::flameletSolver

Description
    Class that holds the table entries and performs the interpolation.

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef FLAMELETTABLE_H_
#define FLAMELETTABLE_H_

#include "fvMesh.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "scalarList.H"

namespace Foam
{

class fvMesh;

class flameletTable
:
	public IOdictionary
{
private:

	word tableName_;

public:

	//- Runtime type information
    TypeName("flameletTable");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        flameletTable,
        dictionary,
        (
            const fvMesh &mesh,
            const word& tableName
        ),
        (mesh, tableName)
    );

    //- Constructor
    flameletTable(const fvMesh& mesh, const word& tableName);

    //- Construct and return a clone
    autoPtr<flameletTable> clone() const;

    //- Destructor
    virtual ~flameletTable();


    // Selectors
    //- Return a reference to the selected tableSolver
    static autoPtr<flameletTable> New(const fvMesh& mesh, const word& tableName);

    virtual scalar interpolate(const List<int>& ub, const scalarList& pos) const = 0;
};

} // End Foam namespace

#endif /* FLAMELETTABLE_H_ */
