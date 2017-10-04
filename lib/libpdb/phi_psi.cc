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

#include "phi_psi.h"

#define TO_DEG  180.0/M_PI

/*****************************************************************************

                'J' - Yet Another Molecular Mechanics Program

                           (c) 1995 Jan T Pedersen


This program is not in the public domain and may not be used or distributed
by anyone without the written permission of the author. The functionality
of the whole or parts of the program is not guaranteed.

******************************************************************************
 $Revision: 1.1 $
******************************************************************************/

/******************************************************************************
           v2    v4
           N     C
          / \   / \
         /   \ /   \
        C     Ca    N
       v1     v3    v5

  phi = M_TorsionDouble(v1,v2,v3,v4);
  psi = M_TorsionDouble(v2,v3,v4,v5);

Alpha helix  phi (-64 +/- 7), psi(-41 +/- 7)
3.10 helix     phi( -74 +/- 2), psi(-4  +/- 2)
pi helix     phi(-57  +/-  1),psi(-67  +/- 1)

sheet        phi(-120 +/- 10), psi(120 +/- 10)

gamma turn   phi(-70),psi(60), i+1 residue phi(70),psi(60)

******************************************************************************/
#define TTOL 0.00001

double M_TorsionDouble(double *v1, double *v2, double *v3, double *v4)
/******************************************************************************
* Calculate the torsion angle in radians
*
* Calculate the angle as the angle between the two planes
* spanned by the two subsets of vectors: AB & BC and BC & CD
*       A function which calculates the torsion angle of vector    
* BC,devised by the vectors AB BC and CD.                          
*       v1- v4 are double numbers (X,Y,Z),coordinates              
* for the 4 points in space defining the vectors AB BC CD.        
*
* Revised for 'J' - March 1995
*
* return a value between -PI and PI
******************************************************************************/
{
	double	p[2][3], n12[3];              /* n1 x n2        */
      	double	nbc[3];              /* unit BC vector */
      	double	n12bc[3];            /* n12 + nbc      */
      	double	n[2][3],v[4][3];

	double vec[3][3];
	double l[2],l12bc=0.0,lb=0.0;        /* square of length of n12bc */
	double coangle=0.0, Angle=0.0;
	int i=0,j=0;

	for(i=0;i<3;i++) {
		v[0][i] = v1[i]; v[1][i] = v2[i];
		v[2][i] = v3[i]; v[3][i] = v4[i];
	}

	for(i=0;i<3;i++){
	  for(j=0;j<3;j++) vec[i][j] = v[i+1][j] - v[i][j];
	}

/*-----------------------------------------------------------------*/
/*      Calculate the normals to the two planes n1 and n2          */
/*      this is given as the cross products :                      */
/*                                                                 */
/*              AB x BC                                            */
/*              ------- = n1                                       */
/*             |AB x BC|                                           */
/*                                                                 */
/*              BC x CD                                            */
/*              ------- = n2                                       */
/*             |BC x CD|                                           */
/*-----------------------------------------------------------------*/

	for(i=0;i<2;i++) {
	  p[i][0] = ((vec[0+i][1]*vec[1+i][2] - vec[0+i][2]*vec[1+i][1]));
	  p[i][1] = ((vec[0+i][2]*vec[1+i][0] - vec[0+i][0]*vec[1+i][2]));
	  p[i][2] = ((vec[0+i][0]*vec[1+i][1] - vec[0+i][1]*vec[1+i][0]));
	}

/*-----------------------------------------------------------------*/
/*      Calculate the length of the two normals and n              */
/*-----------------------------------------------------------------*/
	for(i=0;i<2;i++)
		l[i]=(double)sqrt((p[i][0]*p[i][0])+(p[i][1]*p[i][1])+(p[i][2]*p[i][2]));

/*-----------------------------------------------------------------*/
/*      calculate the two normals n1 and n2                        */
/*-----------------------------------------------------------------*/
	for(i=0;i<2;i++) for(j=0;j<3;j++) n[i][j] = p[i][j]/l[i];

/*-----------------------------------------------------------------*/
/*      calculate cos to torsion angle : coangle = n1.n2           */
/*-----------------------------------------------------------------*/
	for(i=0;i<3;i++) coangle += n[0][i]*n[1][i];

/*-----------------------------------------------------------------*/
/*      In order to check where we are on the unit circle we       */
/*      can do two things we either have to know what sine to the  */
/*      torsion angle is, or what site of the plane n1 and n2 are  */
/*      pointion too.                                              */
/*      Calulate the cross product between n1 and n2 and see if the*/
/*      direction is the same as the unit BC vector.               */
/*-----------------------------------------------------------------*/
	n12[0] = ((n[0][1]*n[1][2]) - (n[0][2]*n[1][1])); 
	n12[1] = ((n[0][2]*n[1][0]) - (n[0][0]*n[1][2]));
	n12[2] = ((n[0][0]*n[1][1]) - (n[0][1]*n[1][0]));

	for(i=0;i<3;i++) lb += (vec[1][i]*vec[1][i]);
	for(i=0;i<3;i++) { nbc[i] = vec[1][i] / ((double)sqrt(lb)); }

/*-----------------------------------------------------------------*/
/*      to check the direction calculate the lenght of the         */
/*      vector nbc + n12.                                          */
/*-----------------------------------------------------------------*/
	for(i=0;i<3;i++) {
	   n12bc[i] = n12[i] + nbc[i]; l12bc += (n12bc[i]*n12bc[i]); 
	}

/*-----------------------------------------------------------------*/
/*      If n12 and nbc have the same directions l12bc == 2         */
/*      If n12 and nbc dont have the same directions l12bc == 0    */
/*      1 is used as the criterium                                 */
/*-----------------------------------------------------------------*/
	if(fabs(coangle+1.0) < TTOL) Angle = M_PI;
	else if ((1.0-coangle) < TTOL) Angle = 0.0;
	else Angle = ((double)acos(coangle));
        
	// fprintf(stdout,"DBG >> M_TorsionDouble() Angle : %f\n",Angle);
	if(l12bc > 1) { return(TO_DEG*Angle); }

/* This part (next 5 lines) is borrowed from ~jan/progs/mc_side/source/Tools.c.
 * The routine it is borrowed from is torr().
 * It is to correct for chi angles so that a value between 0 to 2PI is returned.
 * if you don't want to use this, just use the define appropriately.
 */

#if 1
	return(TO_DEG*(-1.0)*Angle);
#else
#define ZERO_TO_2PI	0

#if defined(ZERO_TO_2PI)
	if(l12bc < 1) { return((2*((double)M_PI))-Angle); }
	return (Angle);
#else
	return((-1.0)*Angle);
#endif
#endif

}

