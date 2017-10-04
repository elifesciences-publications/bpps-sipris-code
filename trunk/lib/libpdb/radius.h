/******************************************************************************************
    Copyright (C) 1997-2014 Andrew F. Neuwald, Cold Spring Harbor Laboratory
    and the University of Maryland School of Medicine.

    Permission is hereby granted, free of charge, to any person obtaining a copy of 
    this software and associated documentation files (the "Software"), to deal in the 
    Software without restriction, including without limitation the rights to use, copy, 
    modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the 
    following conditions:

    The above copyright notice and this permission notice shall be included in all 
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
    PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
    OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
    OTHER DEALINGS IN THE SOFTWARE.

    For further information contact:
         Andrew F. Neuwald
         Institute for Genome Sciences and
         Department of Biochemistry & Molecular Biology
         University of Maryland School of Medicine
         801 West Baltimore St.
         BioPark II, Room 617
         Baltimore, MD 21201
         Tel: 410-706-6724; Fax: 410-706-1482; E-mail: aneuwald@som.umaryland.edu
 ******************************************************************************************/

#include <stdio.h>
#include <string.h>
#include "stdinc.h"

/******************************** private **********************************/
FILE  *stream_from_string (char *string);
char  *atom_type_name(int index);
int   find_type_definition (char *name);
int   allocate_atom_type(char *name);
void  read_atom_types (FILE *ff);
char  *GetAtomNotation(char *residuename, char *atomname);
void  InitializeAtomRadii();
float GetAtomVDWRadius(char *resname, char *atomname);
float GetAtomEffectiveRadius(char *resname, char *atomname);
/***************************************************************************/

									 
/*------------------------------------------------------------------------*/

/* 
   
   These are the size parameters currently used by 
   Mark Gerstein's protein analysis programs.

   General references for my program systems and parameters are:

   Y Harpaz, M Gerstein, & Chothia (1994). 
   "Volume Changes on Protein Folding," Structure 2: 641-649.

   M Gerstein & R Lynden-Bell (1993). 
   "Simulation of Water around a Model Protein Helix. 2. The
   Relative Contributions of Packing, Hydrophobicity, and Hydrogen-Bonding," 
   Journal of Physical Chemistry 97: 2991-2999.

*/

/* We use RICHARDS_RADII currently (not really) */

/* ATOM_TYPE VDW_RADIUS COVALENT_RADIUS CHARGE */

#define RICHARDS_RADII "\
DEFAULT 1.80  0.77\n\
C4	0.00  0.77\n\
C4H	2.00  0.77\n\
C4HH	2.00  0.77\n\
CPRO	2.00  0.77\n\
C4HHH	2.00  0.77\n\
C3	1.70  0.77\n\
C3H	1.85  0.77\n\
C3HH	1.85  0.77\n\
OWAT    1.40  0.77\n\
OW      1.40  0.77\n\
O1	1.40  0.66\n\
O2H	1.60  0.66\n\
O1O2H	1.50  0.66\n\
N4	0.00  0.70\n\
N4H	2.00  0.70\n\
N4HH	2.00  0.70\n\
N4HHH	2.00  0.70\n\
N3	1.50  0.70\n\
N3H	1.70  0.70\n\
N3HH	1.80  0.70\n\
S2	1.85  1.04\n\
S2H	2.00  1.04\n\
P4	0.00  1.10\n\
O1N3HH	1.60  0.68\n\
Z2ION	0.74  0.74\n\
FE	1.70  0.70\n\
"

#define CHC_RADII "\
C4      0.00    0.77\n\
C4H     1.87    0.77\n\
C4HH    1.87    0.77\n\
CPRO    1.87    0.77\n\
C4HHH   1.87    0.77\n\
C3      1.76    0.77\n\
C3H     1.76    0.77\n\
C3HH    0.0     0.77\n\
OWAT    1.40    0.77\n\
OW      1.40    0.77\n\
O1      1.40    0.66\n\
O2H     1.40    0.66\n\
O1O2H   1.40    0.66\n\
N4      0.00    0.70\n\
N4H     1.50    0.70\n\
N4HH    1.50    0.70\n\
N4HHH   1.50    0.70\n\
N3      1.5     0.70\n\
N3H     1.65    0.70\n\
N3HH    1.65    0.70\n\
S2      1.85    1.04\n\
S2H     1.85    1.04\n\
P4      0.00    1.10\n\
O1N3HH  1.60    0.68\n\
Z2ION   0.74    0.74\n\
FE      1.70    0.70\n\
DEFAULT	1.80    0.77\n\
"

/*------------------------------------------------------------------------*/
