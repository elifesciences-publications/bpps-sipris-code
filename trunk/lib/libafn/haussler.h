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

/* haussler.h These are the data structures used in the adaptive prior method */
#if !defined(HAUSSLER)
#define HAUSSLER
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "probability.h"
#define AlphLength 20

/* This structure stores parameters of the prior distributions */
typedef struct{
	Int4	AlphaChar;	/* Number of alphabet characters */
	Int4	L;		/* Number of prior distributions */
	double	*Mix;		/* Mixture coefficents for each prior */
	double	**Distr;	/* Prior distributions. L Dirchlet's 
				   over AlphaChar positions:
				   Distribution[L][AlphaChar+1] */
	Int4	*FullUpdate;	/* !=0 re-estimate all alpha
				   ==0 re-estimate alpha 0 only */
	Int4	*QUpdate;	/* !=0 update mixture coefficents
				   ==0 do not update coefficents */
	char	**StructID;	/* Structure Tag */
	char	**Comment;
	char	let2code[127];	/* DMProior codes for letters */
} DMPriors;

typedef DMPriors *dmp_typ;

/*************************** private *********************************/
dmp_typ	alloc_DMPriors( Int4 L, Int4 Alpha );
void	free_DMPriors (dmp_typ lib);
dmp_typ	read_DMPriors( FILE *fp );
double	*mixture_regularizer(double *freq, dmp_typ Lib);
double	LogAddLog(double x, double y);
double	logpajgy(double *n, dmp_typ Library, Int4 j, BooLean calc);
double  logpygaj(double *n, double *a);
double  log_beta1_dmp(double *n);
double  log_beta2_dmp(double *a, double *n);

/*************************** PUBLIC *********************************/
dmp_typ	MkDMPriors(BooLean full);
void	NilDMPriors(dmp_typ D);
void	CalcDMPriors(double *counts, double *pseudo, dmp_typ D, a_type A);
void	CalcDMPriors2(Int4 *counts, double *pseudo, dmp_typ D, a_type A);

#endif

