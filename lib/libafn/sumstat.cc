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

#include "sumstat.h"

static double    funct_f(register double, void *);
static double    funct_g(register double, void *);

static double RombergIntegrate(double (*f) (double, void *), void *fargs,
	double p, double q, double eps, Int4 epsit, Int4 itmin);

double SumPStd(Int4 r, double s)
/********************************************************************
 Estimate the Sum P-value by calculation or interpolation, as appropriate.
	Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
	r = number of segments
	s = total score (in nats), adjusted by -r*log(KN)
/********************************************************************/
{
	Int4		i, r1, r2;
	double	a;
	static float	tab2[] = { /* table for r == 2 */
		0.01669, 0.0249, 0.03683, 0.05390, 0.07794, 0.1111, 
		0.1559, 0.2146, 0.2890, 0.3794, 0.4836, 0.5965,  
		0.7092, 0.8114, 0.8931, 0.9490, 0.9806, 0.9944, 0.9989 };

	static float	tab3[] = { /* table for r == 3 */
		0.9806, 0.9944, 0.9989, 0.0001682,0.0002542,0.0003829,
		0.0005745,0.0008587,0.001278, 0.001893, 0.002789,
		0.004088, 0.005958, 0.008627, 0.01240, 0.01770, 
		0.02505, 0.03514, 0.04880, 0.06704, 0.09103, 0.1220,
		0.1612, 0.2097,0.2682, 0.3368, 0.4145, 0.4994, 0.5881,
		0.6765, 0.7596, 0.8326,0.8922, 0.9367, 0.9667, 0.9846,
		0.9939, 0.9980};

	static float	tab4[] = { /* table for r == 4 */
		2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,
		3.302e-06,4.990e-06,7.524e-06,1.132e-05,1.698e-05,2.541e-05,
		3.791e-05,5.641e-05,8.368e-05,0.0001237,0.0001823,0.0002677,
		0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
		0.003494, 0.004942, 0.006948, 0.009702, 0.01346, 0.01853,
		0.02532, 0.03431,0.04607, 0.06128, 0.08068, 0.1051, 0.1352,
		0.1719, 0.2157, 0.2669,0.3254, 0.3906, 0.4612, 0.5355, 
		0.6110, 0.6849, 0.7544, 0.8168,0.8699, 0.9127, 0.9451, 
		0.9679, 0.9827, 0.9915, 0.9963};

	static float *table[] = { tab2, tab3, tab4 };
	static short tabsize[] = { 18, 37, 54 };

	if (r == 1) return (s>8. ? exp(-s) : 1.0-exp(-exp(-s)));

	if (r <= 4) {
		if (r < 1) return 0.;
		r1 = r - 1;
		if (s >= r*r+r1) {
			a=exp(r1*log(s)-s)/r;
			for (i=2;i<r;++i) a/=i*i;
			return a;
		}
		if (s > -2*r) {
			/* interpolate */
			i = a = s+s+(4*r);
			a -= i;
			i = tabsize[r2 = r - 2] - i;
			return a*table[r2][i-1] + (1.-a)*table[r2][i];
		}
		return 1.;
	}

	return SumPCalc(r, s);
}

double SumPCalc(Int4 r, double s)
/*********************************************************************
 Evaluate the following double integral, where r = number of segments
 and s = the adjusted score in nats:

      (r-2)   oo   oo
  Prob(r,s) = r    -   - (r-2)
     ------------- | exp(-y) | x exp(-exp(x - y/r)) dx dy
     (r-1)! (r-2)! U   U
         s   0
/*********************************************************************/
{
	Int4		i, r1, itmin;
	double	t, d, a;
	double	mean, stddev, stddev4;
	double	xr, xr1, xr2, logr;
	double	args[6];

	xr = r;

	/* stddev in the limit of infinite r, but quite good for even small r */
	stddev = sqrt(xr+xr-0.5);
	stddev4 = 4.*stddev;
	xr1 = r1 = r - 1;

	/* mean in the limit of infinite r, but quite good for even small r */
	logr = log(xr);
	mean = xr * (1. - logr) - 0.5;
	if (s <= mean - stddev4) return 1.;

	if (s >= mean) { t = s + 6.*stddev; itmin = 1; }
	else { t = mean + 6.*stddev; itmin = 2; }

#define ARG_R args[0]
#define ARG_R2 args[1]
#define ARG_ADJ1 args[2]
#define ARG_ADJ2 args[3]
#define ARG_SDIVR args[4]
#define ARG_EPS args[5]

	ARG_R = xr;
	ARG_R2 = xr2 = r - 2;
	a=log(xr1);
	for (i=2;i<r1;++i) a+=2*log((double) i);
	ARG_ADJ1 = xr2*logr - a;
	ARG_EPS = sump_epsilon;

	do {
		d=RombergIntegrate(funct_g,args,s,t,sump_epsilon,0,itmin);
	} while (s < mean && d < 0.4 && itmin++ < 4);

	return (d < 1. ? d : 1.);
}

static double  funct_g(register double s, void *vp)
{
	register double *args = (double *) vp;
	double	mx;
	
	ARG_ADJ2 = ARG_ADJ1 - s;
	ARG_SDIVR = s / ARG_R;	/* = s / r */
	mx = (s > 0. ? ARG_SDIVR + 3. : 3.);
	return RombergIntegrate(funct_f, vp, 0., mx, ARG_EPS, 0, 1);
}

static double  funct_f(register double x, void *vp)
{
	register double *args = (double *) vp;
	register double	y;

	y = exp(x - ARG_SDIVR);
	if (ARG_R2 == 0.) return exp(ARG_ADJ2 - y);
	if (x == 0.) return 0.;
	return exp(ARG_R2*log(x) + ARG_ADJ2 - y);
}

/*********************************************************************
	Romberg numerical integrator

	Author: Dr. John Spouge, NCBI		Received: July 17, 1992
	Reference:
		Francis Scheid (1968)
		Schaum's Outline Series
		Numerical Analysis, p. 115
		McGraw-Hill Book Company, New York
*********************************************************************/

#define F(x) ((*f)((x), fargs))
#define ROMBERG_ITMAX 20
#define ABS(x)	((x) >= 0. ? (x) : -(x))

static double RombergIntegrate(double (*f) (double, void *), void *fargs,
	double p, double q, double eps, Int4 epsit, Int4 itmin)
{
	double	romb[ROMBERG_ITMAX];	/* present list of Romberg values */
	double	h;	/* mesh-size */
	int	i, j, k, npts;
	Int4	n;	/* 4^(error order in romb[i]) */
	int	epsit_cnt = 0, epsck;
	register double	y;
	register double	x;
	register double	sum;

	/* itmin = min. no. of iterations to perform */
	itmin = 1>itmin ? 1:itmin;
	itmin = itmin<ROMBERG_ITMAX-1 ? itmin : ROMBERG_ITMAX-1;

	/* epsit = min. no. of consecutive iterations that must satisfy epsilon */
	epsit = epsit>1 ? epsit:1; /* default = 1 */
	epsit = epsit<3 ? epsit:3; /* if > 3, the problem needs more prior analysis */

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
			if (ABS(y) == HUGE_VAL)
				return y;
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
			if (i >= itmin && epsit_cnt >= epsit)
				return romb[0];
		}
	}
	fprintf(stderr,"iterations > ROMBERG_ITMAX\n");
	return HUGE_VAL;
}

