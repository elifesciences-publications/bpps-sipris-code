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

/****************** betaprior.h - .***************

  beta prior distribution for Bayesean statistics:

	mean = m = a/(a+b)
	var = S^2 = a*b/((a+b+1)(a+b)^2)

	a = -m - m^2*(m-1)/s^2
	b = -1 +m + m*(m-1)^2/s^2

*************************************************************************/
#if !defined(BETAPRIOR)
#define BETAPRIOR
#include "afnio.h"
#include "stdinc.h"
#include "probability.h"

/********************************* PRIVATE ********************************/
typedef struct {
	double		A;		/* total pseudo trials */
	double		N;		/* total real trials */
	double		T;		/* Total trials */
	double  	alpha;		/* #successful pseudo trials */
	double		beta;		/* #failed pseudo trials */
	UInt4	success;	/* #successful real trials */
	UInt4	expect; 	/* expected number of success */
	double		p;		/* posterior probability of sucess */
	double		weight;		/* A/N = w/(1-w) - fractional weight */ 
	BooLean		calc;
} beta_prior_type;

typedef beta_prior_type *bp_type;

/******************************** private ********************************/
void    bprior_error(char *s);

/********************************* PUBLIC ********************************/
bp_type MkBPrior(Int4 expect, double weight, double N);
bp_type CopyBPrior(bp_type B);
void	ClearBPrior(bp_type B);
void	SetBPriorS(Int4 success, bp_type B);
void	SetBPriorN(double N, bp_type B);
void    PutBPrior(FILE *fptr, bp_type B);
double  BPriorMAP(register bp_type B);
double  LogLikeBPrior(bp_type B);
double  BPriorNullMAP(register bp_type B);
double  RatioBPrior(register bp_type B1, register bp_type B2);
/********************************* MACROS ********************************/
#define NilBPrior(B)	free((B))
#define SuccessBPrior(B) ((B)->success)
#define BPriorN(B)	((B)->N)
#define ExpectBPrior(B)	((B)->expect)
#define alphaBPrior(B)	((B)->alpha)
#define betaBPrior(B)	((B)->beta)
#define PostProbBPrior(B)  (((B)->calc)? ((B)->calc)=FALSE,\
			((B)->p=(B->alpha+(double)B->success)/(B->T-1)):\
			(B)->p)
#define AddBPrior(B)	(((B)->calc)=TRUE, (B)->success++)
#define RmBPrior(B)	(((B)->calc)=TRUE, (B)->success--)

#endif

