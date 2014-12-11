/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Description
    translates cantera table Data to OpenFOAM

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "fvCFD.H"
#include "rhoCombustionModel.H"
//#include "rhoChemistryModel.H"  
#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "canteraReader.H"
#include <stdio.h>
#include <stdlib.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"

   IOdictionary tableDict
   (
       IOobject
       (
          "tableProperties",
          runTime.constant(),
          runTime,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
       )
   );

   Info<< "Creating thermodynamics model\n" << endl;

   autoPtr<combustionModels::rhoCombustionModel> combustion                  
   (
      combustionModels::rhoCombustionModel::New                              
      (
         mesh
      )
   );

   rhoReactionThermo& thermo = combustion->thermo();
   basicMultiComponentMixture& composition = thermo.composition();

   //create dummy tables
   hashedWordList dummytable(thermo.composition().species());
   List<scalar> dummy(0);

   for (int i=0; i<dummytable.size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
            dummytable[i]+"_table",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      word dummyName=dummytable[i]+"_table";
      dictionary.set(dummyName, dummy);
      OFstream output("constant/"+dummytable[i]+"_table");
      dictionary.writeHeader(output);
      dictionary.writeData(output);
   }
   //end creating dummy tables

   //read the cantera-data
   canteraReader canteraRead(tableDict, thermo, composition);

   Info<<"lines = "<<canteraRead.numberOfLines()<<endl;
   Info<<"columns = "<<canteraRead.numberOfColumns()<<endl;

   /*Write the tables, constant/tables/ must exist*/
   for (int i=0; i<canteraRead.getNames().size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
            canteraRead.getNames()[i],
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      OFstream output("constant/"+canteraRead.getNames()[i]+"_table");
      canteraRead.write(i, dictionary, output);
   }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "End\n" << endl;

return 0;
}

