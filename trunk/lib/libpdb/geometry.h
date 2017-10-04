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


#if !defined(_GEOMETRY_H_)
#define _GEOMETRY_H_
#include "stdinc.h"
#include <math.h>


//typedef struct vector_type{

typedef struct {
  float x;
  float y;
  float z;
} vector_type;

typedef vector_type *vec_typ; 


/******************************** Public ***********************************/
vec_typ GetFootInPlane(vec_typ A,  vec_typ  X0, vec_typ X1);
vec_typ GetVectorProjectionOntoPlane(vec_typ A, vec_typ X0, vec_typ X1);
vec_typ GetVectorProduct(vec_typ X1,  vec_typ  X2);
float   GetVectorLength(vec_typ V);
float   GetScalarProduct(vec_typ V, vec_typ U);
int     GetSignOfDihedralAngle(vec_typ A, vec_typ B, vec_typ C);
vec_typ MkVector(float x, float y, float z);
void	NilVector(vec_typ V);
/******************************** MACROS ***********************************/
#define GetVectorX(V)              ((V)->x)
#define GetVectorY(V)              ((V)->y)
#define GetVectorZ(V)              ((V)->z)
/***************************************************************************/

#endif

