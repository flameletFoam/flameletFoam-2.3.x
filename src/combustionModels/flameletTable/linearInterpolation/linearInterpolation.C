/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/


#include "linearInterpolation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<List<scalarList> > linearInterpolation::defaultList(0.0);

defineTypeNameAndDebug(linearInterpolation, 0);
addToRunTimeSelectionTable(flameletTable, linearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolation::linearInterpolation(const fvMesh& mesh, const word& tableName)
:
    flameletTable(mesh, tableName),
	tableValues_(this->lookupOrDefault<List<List<scalarList> > >(tableName, defaultList))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearInterpolation::~linearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

inline scalar linearInterpolation::interpolate(const List<int>& ub, const scalarList& pos) const
{
   // Perform trilinear interpolation
   // for chi
   scalar c00 = tableValues_[ub[0] -1][ub[1] -1][ub[2] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2] -1]*pos[0];
   scalar c10 = tableValues_[ub[0] -1][ub[1]][ub[2] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2] -1]*pos[0];
   scalar c01 = tableValues_[ub[0] -1][ub[1] -1][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2]]*pos[0];
   scalar c11 = tableValues_[ub[0] -1][ub[1]][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2]]*pos[0];
   
   // for Zeta
   scalar c0 = c00*(1-pos[1]) + c10*pos[1];
   scalar c1 = c01*(1-pos[1]) + c11*pos[1];

   // for Z
   return c0*(1-pos[2]) + c1*pos[2];
}

} // End Foam namespace
