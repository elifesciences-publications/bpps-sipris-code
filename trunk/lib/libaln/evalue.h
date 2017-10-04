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

#if !defined (EVALUE)
#define	EVALUE
#include "probability.h"
#include "stdinc.h"
#include "alphabet.h"
#include "random.h"
/*************************** ADT EVALUE ***************************
	
	pattern E-value data type 
E=MkEvalue(d_max,d_min,k_max,pmin,pmax,nsing,mode,N,freq,C_min,A);

**********************************************************************/

/*************************** evalue type **************************/
typedef struct {
	Int4		k_max;		/* maximum length */
	Int4		d_max;		/* maximum depth */
	Int4		d_min;		/* minimum depth */
	Int4		nsingle;	/* minimum # singles */
	Int4		mode;		/* number of pairs */
	Int4		C_min;		/* minset */
	Int4		*N;		/* number of segs */
	double		pmin,pmax;	/* adjusted prob */
	double		*freq;		/* residue freq */
	double		*pfreq;		/* residue freq */
	a_type		A;		/* alphabet */
	double		*pval,*pval0;
	double		pmin0,pmax0,adj_pmax;
	Int4		min,max;	/* (Int4) (pmax/10) */
	double		K,logK;		/* total patterns */
	double		***NP;		/* number patterns */
} evalue_type;
typedef evalue_type *evl_typ;

/******************************* private *****************************/
void    eval_error(char *s);
double  numpat_eval(Int4 d, Int4 k, Int4 s, evl_typ E);
void	estprob_eval(Int4 ntest, evl_typ E);
Int4     sample_pttrn_type(Int4 *D, Int4 *K, Int4 *S, register evl_typ E);
double  sample_pttrn_freq(register Int4 d, register Int4 s, 
						register evl_typ E);

/******************************* Public ******************************/
/***************************** operations ****************************/
evl_typ	MkEvalue(Int4 d_max, Int4 d_min, Int4 k_max, double pmin, double pmax,
	Int4 nsing, Int4 mode, Int4 *N, double *freq, Int4 C_min,a_type A);
double	FastEvalue(evl_typ E, double prob);
double	EstEvalue(Int4 ntest, evl_typ E, double prob);
void    PutEvalue(FILE *fptr, evl_typ E);
void	NilEvalue(evl_typ E);

/**************************** macro operations **********************/
#define	ProbEvalue(E)	((E)->pval)
#define MultEvalue(E,p)	((E)->logK + (p))
#define	pmax0Evalue(E)	((E)->adj_pmax)
#define	pmin0Evalue(E)	((E)->pmin0)
#define	pmaxEvalue(E)	((E)->pmax)
#define	pminEvalue(E)	((E)->pmin)
#define	KEvalue(E)	((E)->K)
#endif

