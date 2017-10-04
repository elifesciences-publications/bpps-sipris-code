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

/*   ncbimath.h
* ===========================================================================
* File Name:  ncbimath.h
* Author:  Gish, Kans, Ostell, Schuler
* Version Creation Date:   10/23/91
* $Revision: 6.0 $
* File Description:
*   	prototypes for portable math library
* Modifications:
* $Log: ncbimath.h,v $
* Revision 6.0  1997/08/25 18:16:37  madden
* Revision changed to 6.0 */

#ifndef _MYNCBIMATH_
#define _MYNCBIMATH_

#include "stdinc.h"
#include "my_ncbi.h"

#ifndef LN2
#define LN2 (0.693147180559945)
#endif
#ifndef LN10
#define LN10 (2.302585092994046)
#endif

/* log(x+1) for all x > -1 */
double GNlm_Log1p(double );

/* exp(x)-1 for all x */
double GNlm_Expm1 (double);

/* Factorial */
double GNlm_Factorial (Int4);

/* gamma(x) */
double GNlm_Gamma (double);

/* log(gamma(x)) */
double GNlm_LnGamma (double);

/* log(gamma(n)), integral n */
double GNlm_LnGammaInt (Int4);

/* digamma(x) 1st order derivative of log(gamma(x)) */
double GNlm_DiGamma (double);

/* trigamma(x) 2nd order derivative of log(gamma(x)) */
double GNlm_TriGamma (double);

/* Nth order derivative of log(gamma) */
double GNlm_PolyGamma (double x, Int4 order);

/* Change gamma coefficients */
void GNlm_GammaCoeffSet (double *coef, unsigned dimension);

/* Nth order derivative of ln(u) */
double GNlm_LogDerivative (Int4 order, double *u);


/* Combined Newton-Raphson and Bisection root solver */
double GNlm_NRBis (double y, double (*f) (double), 
	double (*df) (double), double p, double x, 
	double q, double tol);

/* Romberg numerical integrator */
double GNlm_RombergIntegrate (double (*f) (double, 
	void *), void *fargs, double p, double q, 
	double eps, Int4 epsit, Int4 itmin);

/* Greatest common divisor */
Int4 GNlm_Gcd (Int4, Int4);

/* Nearest integer */
Int4 GNlm_Nint (double);

/* Integral power of x */
double GNlm_Powi (double x, Int4 n);

/* Random no. seeder and generator */
void GNlm_RandomSeed (Int4 n);
Int4 GNlm_RandomNum(void);


#define Log1p	Nlm_Log1p
#define Expm1	Nlm_Expm1
#define Factorial	Nlm_Factorial
#define Gamma	Nlm_Gamma
#define LnGamma	Nlm_LnGamma
#define DiGamma	Nlm_DiGamma
#define TriGamma	Nlm_TriGamma
#define PolyGamma	Nlm_PolyGamma
#define GammaCoeffSet	Nlm_GammaCoeffSet
#define LogDerivative	Nlm_LogDerivative
#define NRBis	Nlm_NRBis
// #define RombergIntegrate	Nlm_RombergIntegrate
#define Gcd	Nlm_Gcd
#define Nint	Nlm_Nint
#define Powi	Nlm_Powi
#define RandomSeed	Nlm_RandomSeed
#define RandomNum	Nlm_RandomNum

#define LOGDERIV_ORDER_MAX	4
#define POLYGAMMA_ORDER_MAX	LOGDERIV_ORDER_MAX

#define NCBIMATH_PI	3.1415926535897932384626433832795
#define NCBIMATH_E	2.7182818284590452353602874713527
/* Euler's constant */
#define NCBIMATH_EULER 0.5772156649015328606065120900824
/* Catalan's constant */
#define NCBIMATH_CATALAN	0.9159655941772190150546035149324

/* sqrt(2) */
#define NCBIMATH_SQRT2	1.4142135623730950488016887242097
/* sqrt(3) */
#define NCBIMATH_SQRT3	1.7320508075688772935274463415059
/* sqrt(PI) */
#define NCBIMATH_SQRTPI 1.7724538509055160272981674833411
/* Natural log(2) */
#define NCBIMATH_LN2	0.69314718055994530941723212145818
/* Natural log(10) */
#define NCBIMATH_LN10	2.3025850929940456840179914546844
/* Natural log(PI) */
#define NCBIMATH_LNPI	1.1447298858494001741434273513531

#endif /* !_MYNCBIMATH_ */

