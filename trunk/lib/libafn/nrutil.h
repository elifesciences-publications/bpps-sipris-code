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

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include "stdinc.h"

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

#if 0
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static Int4 lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static Int4 lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
#endif

// #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// #if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
double *vector(Int4 nl, Int4 nh);
int *ivector(Int4 nl, Int4 nh);
unsigned char *cvector(Int4 nl, Int4 nh);
UInt4 *lvector(Int4 nl, Int4 nh);
double *dvector(Int4 nl, Int4 nh);
double **matrix(Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
double **dmatrix(Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
int **imatrix(Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
double **submatrix(double **a, Int4 oldrl, Int4 oldrh, Int4 oldcl, Int4 oldch,
	Int4 newrl, Int4 newcl);
double **convert_matrix(double *a, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
double ***f3tensor(Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch, Int4 ndl, Int4 ndh);
void free_vector(double *v, Int4 nl, Int4 nh);
void free_ivector(int *v, Int4 nl, Int4 nh);
void free_cvector(unsigned char *v, Int4 nl, Int4 nh);
void free_lvector(UInt4 *v, Int4 nl, Int4 nh);
void free_dvector(double *v, Int4 nl, Int4 nh);
void free_matrix(double **m, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
void free_dmatrix(double **m, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
void free_imatrix(int **m, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
void free_submatrix(double **b, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
void free_convert_matrix(double **b, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch);
void free_f3tensor(double ***t, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch,
	Int4 ndl, Int4 ndh);

#endif /* _NR_UTILS_H_ */
