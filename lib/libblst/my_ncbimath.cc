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

/* ========================================================================
* File Name:  ncbimath.c
* Author:  Gish, Kans, Ostell, Schuler
* Version Creation Date:   10/23/91
* $Revision: 6.2 $
* File Description:
*   	portable math functions
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*/

// #define THIS_MODULE g_corelib
// #define THIS_FILE _this_file

#include "my_ncbimath.h"

#ifdef OS_MAC
#ifdef COMP_METRO
/*#include <fp.h>*/
#ifndef HUGE_VAL
#define HUGE_VAL __inf()
double_t __inf ( void );
#endif
#endif
#endif

// extern char * g_corelib;
// static char * _this_file = __FILE__;

/*
    GNlm_Expm1(x)
    Return values accurate to approx. 16 digits for the quantity exp(x)-1
    for all x.
*/
double GNlm_Expm1(register double x)
{
  register double absx;

  if ((absx = ABS(x)) > .33) return exp(x) - 1.;
  if (absx < 1.e-16) return x;
  return x * (1. + x *
       (0.5 + x * (1./6. + x *
             (1./24. + x * (1./120. + x *
                   (1./720. + x * (1./5040. + x *
                        (1./40320. + x * (1./362880. + x *
                            (1./3628800. + x * (1./39916800. + x *
                                (1./479001600. + x/6227020800.)
                                     ))
                                          ))
                                               ))
                                         )))));
}

// GNlm_Log1p(x) Return accurate values for the quantity log(x+1) for all x > -1.
double GNlm_Log1p(register double x)
{
	register Int4	i;
	register double sum, y;

	if (ABS(x) >= 0.2) return log(x+1.);
	/* Limit the loop to 500 terms. */
	for (i=0, sum=0., y=x; i<500 ; ) {
		sum += y/++i;
		if (ABS(y) < DBL_EPSILON) break;
		y *= x;
		sum -= y/++i;
		if (y < DBL_EPSILON) break;
		y *= x;
	} return sum;
}


/*
** Special thanks to Dr. John ``Gammas Galore'' Spouge for deriving the
** method for calculating the gamma coefficients and their use.
** (See the #ifdef-ed program included at the end of this file).
**/

/*
    For discussion of the Gamma function, see "Numerical Recipes in C",
    Press et al. (1988) pages 167-169.
*/
static double general_lngamma (double x,Int4 n);

static double _default_gamma_coef [] = {
	 4.694580336184385e+04,
	-1.560605207784446e+05,
	 2.065049568014106e+05,
	-1.388934775095388e+05,
	 5.031796415085709e+04,
	-9.601592329182778e+03,
	 8.785855930895250e+02,
	-3.155153906098611e+01,
	 2.908143421162229e-01,
	-2.319827630494973e-04,
	 1.251639670050933e-10
	 };
static double *gamma_coef = _default_gamma_coef;
static unsigned	gamma_dim = DIM(_default_gamma_coef);
static double xgamma_dim = DIM(_default_gamma_coef);

void Nlm_GammaCoeffSet(double *cof, unsigned dim) /* changes gamma coeffs */
{
	if (dim < 3 || dim > 100) /* sanity check */
		return;
	gamma_coef = cof;
	xgamma_dim = gamma_dim = dim;
}


static double general_lngamma(double x, Int4 order)
/* nth derivative of ln[gamma(x)] */
/* x is 10-digit accuracy achieved only for 1 <= x */
/* order is order of the derivative, 0..POLYGAMMA_ORDER_MAX */
{
	Int4		i;
	double		xx, tx;
	double		y[POLYGAMMA_ORDER_MAX+1];
	register double	tmp, value, *coef;

	xx = x - 1.;  /* normalize from gamma(x + 1) to xx! */
	tx = xx + xgamma_dim;
	for (i = 0; i <= order; ++i) {
		tmp = tx;
		/* sum the least significant terms first */
		coef = &gamma_coef[gamma_dim];
		if (i == 0) {
			value = *--coef / tmp;
			while (coef > gamma_coef)
				value += *--coef / --tmp;
		}
		else {
			value = *--coef / GNlm_Powi(tmp, i + 1);
			while (coef > gamma_coef)
				value += *--coef / GNlm_Powi(--tmp, i + 1);
			tmp = GNlm_Factorial(i);
			value *= (i%2 == 0 ? tmp : -tmp);
		}
		y[i] = value;
	}
	++y[0];
	value = GNlm_LogDerivative(order, y);
	tmp = tx + 0.5;
	switch (order) {
	case 0:
		value += ((NCBIMATH_LNPI+NCBIMATH_LN2) / 2.)
					+ (xx + 0.5) * log(tmp) - tmp;
		break;
	case 1:
		value += log(tmp) - xgamma_dim / tmp;
		break;
	case 2:
		value += (tmp + xgamma_dim) / (tmp * tmp);
		break;
	case 3:
		value -= (1. + 2.*xgamma_dim / tmp) / (tmp * tmp);
		break;
	case 4:
		value += 2. * (1. + 3.*xgamma_dim / tmp) / (tmp * tmp * tmp);
		break;
	default:
		tmp = GNlm_Factorial(order - 2) * GNlm_Powi(tmp, 1 - order)
				* (1. + (order - 1) * xgamma_dim / tmp);
		if (order % 2 == 0)
			value += tmp;
		else
			value -= tmp;
		break;
	}
	return value;
}


double	GNlm_PolyGamma(double x, Int4 order) /* ln(ABS[gamma(x)]) - 10 digits of accuracy */
	/* x is and derivatives */
	/* order is order of the derivative */
/* order = 0, 1, 2, ...  ln(gamma), digamma, trigamma, ... */
/* CAUTION:  the value of order is one less than the suggested "di" and
"tri" prefixes of digamma, trigamma, etc.  In other words, the value
of order is truly the order of the derivative.  */
{
	Int4		k;
	register double value, tmp;
	double 	y[POLYGAMMA_ORDER_MAX+1], sx;

	if (order < 0 || order > POLYGAMMA_ORDER_MAX) {
		print_error("PolyGamma: unsupported derivative order");
		return HUGE_VAL;
	}
	if (x >= 1.) return general_lngamma(x, order);
	if (x < 0.) {
		value = general_lngamma(1. - x, order);
		value = ((order - 1) % 2 == 0 ? value : -value);
		if (order == 0) {
			sx = sin(NCBIMATH_PI * x);
			sx = ABS(sx);
			if ( (x < -0.1 && (ceil(x) == x || sx < 2.*DBL_EPSILON))
					|| sx == 0.) {
				print_error("PolyGamma: log(0)");
				return HUGE_VAL;
			}
			value += NCBIMATH_LNPI - log(sx);
		} else {
			y[0] = sin(x *= NCBIMATH_PI);
			tmp = 1.;
			for (k = 1; k <= order; k++) {
				tmp *= NCBIMATH_PI;
				y[k] = tmp * sin(x += (NCBIMATH_PI/2.));
			}
			value -= GNlm_LogDerivative(order, y);
		}
	} else {
		value = general_lngamma(1. + x, order);
		if (order == 0) {
			if (x == 0.) {
				print_error("PolyGamma: log(0)");
				return HUGE_VAL;
			}
			value -= log(x);
		} else {
			tmp = GNlm_Factorial(order - 1) * GNlm_Powi(x,  -order);
			value += (order % 2 == 0 ? tmp : - tmp);
		}
	} return value;
}

double GNlm_LogDerivative(Int4 order, double *u) /* nth derivative of ln(u) */
	/* order is order of the derivative */
	/* u is values of u, u', u", etc. */
{
  Int4		i;
  double 	y[LOGDERIV_ORDER_MAX+1];
  register double value, tmp;

  if (order < 0 || order > LOGDERIV_ORDER_MAX) {
InvalidOrder:
    print_error("LogDerivative: unsupported derivative order");
    return HUGE_VAL;
  }
  if (order > 0 && u[0] == 0.) {
    print_error("LogDerivative: divide by 0"); return HUGE_VAL;
  }
  for (i = 1; i <= order; i++) y[i] = u[i] / u[0];
  switch (order) {
  case 0:
    if (u[0] > 0.) value = log(u[0]);
    else { print_error("LogDerivative: log(x<=0)"); return HUGE_VAL; }
    break;
  case 1: value = y[1]; break;
  case 2: value = y[2] - y[1] * y[1]; break;
  case 3: value = y[3] - 3. * y[2] * y[1] + 2. * y[1] * y[1] * y[1];
    break;
  case 4:
    value = y[4] - 4. * y[3] * y[1] - 3. * y[2] * y[2]
      + 12. * y[2] * (tmp = y[1] * y[1]);
    value -= 6. * tmp * tmp;
    break;
  default: goto InvalidOrder;
  }
  return value;
}


double GNlm_Gamma(double  x)	/* ABS[gamma(x)] - 10 digits of accuracy */
{ return exp(GNlm_PolyGamma(x, 0)); }

double GNlm_LnGamma(double x)	/* ln(ABS[gamma(x)]) - 10 digits of accuracy */
{ return GNlm_PolyGamma(x, 0); }

double GNlm_DiGamma(double x)	/* digamma, 1st order derivative of log(gamma) */
{ return GNlm_PolyGamma(x, 1); }

double GNlm_TriGamma(double x)	/* trigamma, 2nd order derivative of log(gamma) */
{ return GNlm_PolyGamma(x, 2); }

#define FACTORIAL_PRECOMPUTED	36

double GNlm_Factorial(Int4 n)
{
	static double precomputed[FACTORIAL_PRECOMPUTED]
		= { 1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.};
	static Int4	nlim = 10;
	register Int4	m;
	register double x;

	if (n >= 0) {
		if (n <= nlim) return precomputed[n];
		if (n < DIM(precomputed)) {
			for (x = precomputed[m = nlim]; m < n; ) {
				++m;
				precomputed[m] = (x *= m);
			}
			nlim = m;
			return x;
		} return exp(GNlm_LnGamma((double )(n+1)));
	} return 0.0; /* Undefined! */
}

/* GNlm_LnGammaInt(n) -- return log(Gamma(n)) for integral n */
double GNlm_LnGammaInt(Int4 n)
{
	static double precomputed[FACTORIAL_PRECOMPUTED];
	static Int4	nlim = 1; /* first two entries are 0 */
	register Int4 m;

	if (n >= 0) {
		if (n <= nlim) return precomputed[n];
		if (n < DIM(precomputed)) {
			for (m = nlim; m < n; ++m) {
				precomputed[m+1] = log(GNlm_Factorial(m));
			}
			return precomputed[nlim = m];
		}
	} return GNlm_LnGamma((double )n);
}

/*
	Combined Newton-Raphson and Bisection root-finder

	Original Function Name:  Inv_Xnrbis()
	Author:  Dr. John Spouge
	Location:  NCBI
	Received:  July 16, 1991
*/
#define F(x)  ((*f)(x)-y)
#define DF(x) ((*df)(x))
#define NRBIS_ITMAX 100

double GNlm_NRBis(double y, double (*f)(double ), 
	double (*df)(double ), double p, double x, 
	double q, double tol)		/* tolerance */
{
	double temp;	/* for swapping end-points if necessary */
	double dx;		/* present interval length */
	double dxold;	/* old interval length */
	double fx;		/* f(x)-y */
	double dfx;	/* Df(x) */
	Int4	j;                       /* loop index */
	double fp, fq;

	/* Preliminary checks for improper bracketing and end-point root. */
	if ((fp = F(p)) == 0.) return p;
	if ((fq = F(q)) == 0.) return q;
	if ((fp > 0. && fq > 0.) || (fp < 0. && fq < 0.)) {
		print_error("NRBis: root not bracketed"); return HUGE_VAL;
	}
	/* Swaps end-points if necessary to make F(p)<0.<F(q) */
	if (fp > 0.) { temp = p; p = q; q = temp; }
/* Set up the Bisection & Newton-Raphson iteration. */
	if ((x-p) * (x-q) > 0.) x = 0.5*(p+q);
	dxold = dx = p-q;
	for (j = 1; j <= NRBIS_ITMAX; ++j) {
		fx = F(x);
		if (fx == 0.) /* Check for termination. */ return x;
		if (fx < 0.) p = x; else q = x;
		dfx = DF(x);
/* Check: root out of bounds or bisection faster than Newton-Raphson? */
		if ((dfx*(x-p)-fx)*(dfx*(x-q)-fx) >= 0. || 2.*ABS(fx) > ABS(dfx*dx)) {
			dx = dxold;               /* Bisect */
			dxold = 0.5*(p-q);
			x = 0.5*(p+q);
			if (ABS(dxold) <= tol) return x;
		} else {
			dx = dxold;               /* Newton-Raphson */
			dxold = fx/dfx;
			x -= dxold;
			if (ABS(dxold) < tol && F(x-SIGN(dxold)*tol)*fx < 0.)
				return x;
		}
	}
	print_error("NRBis: iterations > NRBIS_ITMAX"); return HUGE_VAL;
}
#undef F /* clean-up */
#undef DF /* clean-up */



/*
	Romberg numerical integrator

	Author:  Dr. John Spouge, NCBI
	Received:  July 17, 1992
	Reference:
		Francis Scheid (1968)
		Schaum's Outline Series
		Numerical Analysis, p. 115
		McGraw-Hill Book Company, New York
*/
#define F(x)  ((*f)((x), fargs))
#define ROMBERG_ITMAX 20

double GNlm_RombergIntegrate(double (*f) (double ,void *), void *fargs,
	double p, double q, double eps, Int4 epsit, Int4 itmin)

{
	double romb[ROMBERG_ITMAX];	/* present list of Romberg values */
	double h;	/* mesh-size */
	Int4		i, j, k, npts;
	Int4	n;	/* 4^(error order in romb[i]) */
	Int4		epsit_cnt = 0, epsck;
	register double y;
	register double x;
	register double sum;

	/* itmin = min. no. of iterations to perform */
	itmin = MAX(1, itmin);
	itmin = MIN(itmin, ROMBERG_ITMAX-1);

	/* epsit = min. no. of consecutive iterations that must satisfy epsilon */
	epsit = MAX(epsit, 1); /* default = 1 */
	epsit = MIN(epsit, 3); /* if > 3, the problem needs more prior analysis */

	epsck = itmin - epsit;

	npts = 1;
	h = q - p;
	x = F(p);
	if (ABS(x) == HUGE_VAL) return x;
	y = F(q);
	if (ABS(y) == HUGE_VAL) return y;
	romb[0] = 0.5 * h * (x + y);	/* trapezoidal rule */
	for (i = 1; i < ROMBERG_ITMAX; ++i, npts *= 2, h *= 0.5) {
		sum = 0.;	/* sum of ordinates for */
		/* x = p+0.5*h, p+1.5*h, ..., q-0.5*h */
		for (k = 0, x = p+0.5*h; k < npts; k++, x += h) {
			y = F(x);
			if (ABS(y) == HUGE_VAL) return y;
			sum += y;
		}
		romb[i] = 0.5 * (romb[i-1] + h*sum); /* new trapezoidal estimate */

		/* Update Romberg array with new column */
		for (n = 4, j = i-1; j >= 0; n *= 4, --j)
			romb[j] = (n*romb[j+1] - romb[j]) / (n-1);

		if (i > epsck) {
			if (ABS(romb[1] - romb[0]) > eps * ABS(romb[0])) {
				epsit_cnt = 0;
				continue;
			}
			++epsit_cnt;
			if (i >= itmin && epsit_cnt >= epsit) return romb[0];
		}
	}
	print_error("RombergIntegrate: iterations > ROMBERG_ITMAX");
	return HUGE_VAL;
}

/*
	GNlm_Gcd(a, b)
	Return the greatest common divisor of a and b.
	Adapted 8-15-90 by WRG from code by S. Altschul.
*/
Int4 GNlm_Gcd(register Int4 a, register Int4 b)
{
	register Int4	c;

	b = ABS(b);
	if (b > a) c=a, a=b, b=c;
	while (b != 0) { c = a%b; a = b; b = c; }
	return a;
}

/* Round a floating point number to the nearest integer */
Int4 GNlm_Nint(register double x)	/* argument */
{ x += (x >= 0. ? 0.5 : -0.5); return (Int4)x; }

/* integer power function
Original submission by John Spouge, 6/25/90
Added to shared library by WRG */
double GNlm_Powi(double x, Int4 n)	/* power */
{
	double y;

	if (n == 0) return 1.;
	if (x == 0.) {
		if (n < 0) { print_error("Powi: divide by 0"); return HUGE_VAL; }
		return 0.;
	}
	if (n < 0) { x = 1./x; n = -n; }
	while (n > 1) {
		if (n&1) { y = x; goto Loop2; }
		n /= 2; x *= x;
	} return x;

Loop2:
	n /= 2;
	x *= x;
	while (n > 1) { if (n&1) y *= x; n /= 2; x *= x; }
	return y * x;
}

/*
	Additive random number generator

	Modelled after "Algorithm A" in
	Knuth, D. E. (1981). The art of computer programming, volume 2, page 27.

	7/26/90 WRG
*/

static Int4	state[33] = {
	(Int4)0xd53f1852,  (Int4)0xdfc78b83,  (Int4)0x4f256096,  (Int4)0xe643df7,
	(Int4)0x82c359bf,  (Int4)0xc7794dfa,  (Int4)0xd5e9ffaa,  (Int4)0x2c8cb64a,
	(Int4)0x2f07b334,  (Int4)0xad5a7eb5,  (Int4)0x96dc0cde,  (Int4)0x6fc24589,
	(Int4)0xa5853646,  (Int4)0xe71576e2,  (Int4)0xdae30df,  (Int4)0xb09ce711,
	(Int4)0x5e56ef87,  (Int4)0x4b4b0082,  (Int4)0x6f4f340e,  (Int4)0xc5bb17e8,
	(Int4)0xd788d765,  (Int4)0x67498087,  (Int4)0x9d7aba26,  (Int4)0x261351d4,
	(Int4)0x411ee7ea,  (Int4)0x393a263,  (Int4)0x2c5a5835,  (Int4)0xc115fcd8,
	(Int4)0x25e9132c,  (Int4)0xd0c6e906,  (Int4)0xc2bc5b2d,  (Int4)0x6c065c98,
	(Int4)0x6e37bd55 };

#define r_off	12

static Int4	*rJ = &state[r_off],
			*rK = &state[DIM(state)-1];

void GNlm_RandomSeed(Int4 x)
{
	register size_t i;

	state[0] = x;
	/* linear congruential initializer */
	for (i=1; i<DIM(state); ++i) {
		state[i] = 1103515245*state[i-1] + 12345;
	}
	rJ = &state[r_off];
	rK = &state[DIM(state)-1];
	for (i=0; i<10*DIM(state); ++i) (void) GNlm_RandomNum();
}


/* GNlm_RandomNum --  return value in the range 0 <= x <= 2**31 - 1 */
Int4 GNlm_RandomNum(void)
{
	register Int4	r;

	r = *rK;
	r += *rJ--;
	*rK-- = r;

	if (rK < state) rK = &state[DIM(state)-1];
	else if (rJ < state) rJ = &state[DIM(state)-1];
	return (r>>1)&0x7fffffff; /* discard the least-random bit */
}

