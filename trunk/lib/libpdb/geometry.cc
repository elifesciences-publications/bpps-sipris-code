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

#include "geometry.h"


/*------------------------------------------------------------------------*/

vec_typ GetFootInPlane(vec_typ A, vec_typ X0, vec_typ X1)
{
//   Find the foot in the plane P from the point X1=(x1, y1, z1), where
//   P is perpendicular to the vector A=(a, b, c) and contains the point
//   X0=(x0, y0, z0).

  float a, b, c, d;
  float x0, y0, z0, x1, y1, z1, xh, yh, zh;
  vec_typ Xh;

  a  = GetVectorX(A);   b  = GetVectorY(A);   c  = GetVectorZ(A);
  x0 = GetVectorX(X0);  y0 = GetVectorY(X0);  z0 = GetVectorZ(X0);
  x1 = GetVectorX(X1);  y1 = GetVectorY(X1);  z1 = GetVectorZ(X1);

  /* P can be defined as
     ax + by + cz + d =0 
  */

  d = - (a*x0 + b*y0 + c*z0);

  xh = ( (b*b+c*c) * x1 - a * (b*y1 + c*z1 + d) ) / (a*a + b*b + c*c);  
  yh = ( (c*c+a*a) * y1 - b * (a*x1 + c*z1 + d) ) / (a*a + b*b + c*c);  
  zh = ( (a*a+b*b) * z1 - c * (a*x1 + b*y1 + d) ) / (a*a + b*b + c*c);

  Xh = MkVector(xh, yh, zh);

  return Xh;
}
  


/*------------------------------------------------------------------------*/

vec_typ GetVectorProjectionOntoPlane(vec_typ A, vec_typ X0, vec_typ X1)
{
//   Find the vector X0h = X_h - X0, which is the projection of the vector 
//   X01 = X1 - X0, where X0 = (x0, y0, z0) is a point in the plane P.  
//   X1 = (x1, y1, z1) sits outside P; P is perpendicular to the 
//   vector A=(a, b, c).

  float a, b, c, d;
  float x0, y0, z0, x1, y1, z1, xh, yh, zh;
  vec_typ X0h;

  a  = GetVectorX(A);   b  = GetVectorY(A);   c  = GetVectorZ(A);
  x0 = GetVectorX(X0);  y0 = GetVectorY(X0);  z0 = GetVectorZ(X0);
  x1 = GetVectorX(X1);  y1 = GetVectorY(X1);  z1 = GetVectorZ(X1);

  /* P can be defined as ax + by + cz + d =0 */

  float LenAxLenA=(a*a + b*b + c*c);
  d = - (a*x0 + b*y0 + c*z0);

  xh = ( (b*b+c*c) * x1 - a * (b*y1 + c*z1 + d) ) / LenAxLenA;
  yh = ( (c*c+a*a) * y1 - b * (a*x1 + c*z1 + d) ) / LenAxLenA;
  zh = ( (a*a+b*b) * z1 - c * (a*x1 + b*y1 + d) ) / LenAxLenA;

  X0h = MkVector(xh - x0, yh - y0, zh - z0);

  return X0h;
}

/*------------------------------------------------------------------------*/

vec_typ GetVectorProduct(vec_typ V1,  vec_typ V2)

{

//  V3 = (x3, y3, z3) = V1 x V2, 
//  where V1 = (x1, y1, z1) and V2 = (x2, y2, z2) 

  float x1, y1, z1, x2, y2, z2, x3, y3, z3;
  vec_typ V3;

  x1 = GetVectorX(V1);  y1 = GetVectorY(V1);  z1 = GetVectorZ(V1);
  x2 = GetVectorX(V2);  y2 = GetVectorY(V2);  z2 = GetVectorZ(V2);

  x3 = y1*z2 - y2*z1;
  y3 = z1*x2 - z2*x1;
  z3 = x1*y2 - x2*y1;    

  V3 = MkVector(x3, y3, z3);

  return V3;
}

/*------------------------------------------------------------------------*/

float GetVectorLength(vec_typ V)
{

// Return the length of the vector V

  float x = V->x, y = V->y, z = V->z;

  return x*x+y*y+z*z;
}
  

/*------------------------------------------------------------------------*/

float GetScalarProduct(vec_typ V, vec_typ U)
{

// Return the scalar product of V and U

  float product = V->x * U->x + V->y * U->y + V->z * U->z;

  return product;
}

/*------------------------------------------------------------------------*/

int GetSignOfDihedralAngle(vec_typ A, vec_typ B, vec_typ C)
{

// Determine the sign of (A x B) * C
//    +  if A x B is parallel to C.
//    -  if A x B is anti-parallel to C.

  vec_typ V;

  V = GetVectorProduct(A, B);

  if      ( GetScalarProduct(V, C) < 0) return -1;
  else if ( GetScalarProduct(V, C) > 0) return 1;

  return 0;
}

void	NilVector(vec_typ V)
{ delete V; }

/*------------------------------------------------------------------------*/
vec_typ  MkVector(float x, float y, float z)
{
  vec_typ V;
  V = new vector_type;
  V->x = x; V->y = y; V->z = z;
  return V;
}
/*------------------------------------------------------------------------*/
