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
    ANY WARRANTY; without even the implied warranty of MERCHANTABSILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    
Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/
%{
#undef yyFlexLexer
#define yyFlexLexer xxFlexLexer
#include "canteraReader.H"
#include <fstream>

// flex input buffer size
int Foam::canteraReader::yyBufSize = YY_BUF_SIZE;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! @cond dummy
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! @endcond dummy

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorIn("yyFlexLexer::yylex()")
        << "should not have called this function"
        << abort(Foam::FatalError);

    return 0;
}
//! \endcond


#define YY_DECL int Foam::canteraReader::lex()
%}

int	num_lines=0, num_columns=0, col_iter=0;
bool noSecondLine=1;
 /* ------------------------------------------------------------------------- *\
                      -----  define Names -----
 \* ------------------------------------------------------------------------- */

columnSeparator		  ","
alpha                 [_A-Za-z]
digit                 [0-9]

exponent_part         [eEdD][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))
word                  [A-Za-z](([A-Za-z0-9)*+-])|("("[^+]))*
molefraction          "mole fr."
massfraction          "mass fr."
sharpie               "#"


 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%x readSecondLine
%x readBody
%%

 /* ------------------------------------------------------------------------- *\
                      -----  read the first Line  -----
 \* ------------------------------------------------------------------------- */

{word} {
         word name(YYText());
         if (tablesToBeRead_.contains(name))
         {	
            tableSorted_.append(name);
         }
         else
         {
            tableSorted_.append("IGNORE");
         }
       }
		
{columnSeparator}   ++col_iter;
. {}					
					
\n  {BEGIN(readBody); num_columns=col_iter+1;col_iter=0;
   for (int k=0; k<tableSorted_.size(); k++)
   {
       if (tableSorted_[k] !="IGNORE")
       {
          int i(tablesToBeRead_[tableSorted_[k]]);
          columns_[i]=k;
       }
   }
 }

 /* ------------------------------------------------------------------------- *\
                      -----  read the second Line  -----
 \* ------------------------------------------------------------------------- */
<readBody>{sharpie}   noSecondLine=0; Info<<"sharpie found!"<<endl; BEGIN(readSecondLine);
<readSecondLine>{massfraction}	{
   if (tableSorted_[col_iter]!="IGNORE")
   {
      /*this is the place in data_, where the list of scalars is stored*/
      int i(tablesToBeRead_[tableSorted_[col_iter]]);
      columns_[i]=col_iter;
   }
}
<readSecondLine>{molefraction}	{
   if (tableSorted_[col_iter]!="IGNORE")
   {
       /*this is the place in data_, where the list of scalars is stored*/
       //int i(tablesToBeRead_[tableSorted_[col_iter]]);
       //columnsMolFract_[i]=col_iter;
   }
 }
<readSecondLine>{word} {}

<readSecondLine>{columnSeparator}	++col_iter;
<readSecondLine>.			{}

<readSecondLine>\n	{BEGIN(readBody); col_iter=0;}	
 /* ------------------------------------------------------------------------- *\
                      -----  read the body  -----
 \* ------------------------------------------------------------------------- */
<readBody>{columnSeparator}	++col_iter;

<readBody>\n  {
    ++num_lines;
    col_iter=0;
   }


<readBody>{floatNum}  {
   for (int i=0; i<=tablesToBeRead_.size()-1; i++)
   {
       if (columns_[i]==col_iter)
       {
          singleData_[i].append(atof(YYText()));
       }
   }
  }
<readBody>a  {}
<<EOF>>	  {BEGIN(INITIAL);	
  yyterminate();}

%%
   
