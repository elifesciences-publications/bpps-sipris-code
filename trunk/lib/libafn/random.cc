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

#include "random.h"

Int4    random_integer(register Int4 max)
/** return a random 0 <= integer < max **/
{ return Random()%max; }

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

static Int4	*rJ = &state[r_off], *rK = &state[DIM(state)-1];

void	SampleMultinomial(UInt4 Number, unsigned short classes, 
	UInt4 *observed, double *freq)
{
        double		rand_no,total;
	unsigned short	b;

        for(total=0.0,b=1; b<=classes; b++){ total+=freq[b]; observed[b]=0; }
        for(UInt4 i=0; i < Number; i++){
           rand_no=total*SampleUniformProb();
           for(b=1; b<=classes; b++){
                if((rand_no-=freq[b]) <= 0.0){ observed[b]++; break; }
           } assert(rand_no <= 0.0);
        }
}

void sRandom(Int4 x)
{
	register Int4	i;

	state[0] = x;
	/* linear congruential initializer */
	for (i=1; i<DIM(state); ++i) {
		state[i] = 1103515245*state[i-1] + 12345;
	}

	rJ = &state[r_off];
	rK = &state[DIM(state)-1];

	for (i=0; i<10*DIM(state); ++i) (void) Random();
}

Int4	*RandomArray(Int4 N)
/** return a random array of integers from 1 to N **/
{
	dh_type	H;
	Int4	i,j,*L;

	H = dheap(N,3);
	NEW(L,N+2,Int4);
	for(i=1; i<=N; i++) insrtHeap(i,(keytyp) Random(),H);
	for(i=1;(j=delminHeap(H)) != 0; i++) L[i]=j;
	Nildheap(H);
	return L;
}

Int4 Random(void)
/** Random --  return value in the range 0 <= x <= 2**31 - 1 **/
{
	register Int4	r;

	r = *rK;
	r += *rJ--;
	*rK-- = r;

	if (rK < state) rK = &state[DIM(state)-1];
	else if (rJ < state) rJ = &state[DIM(state)-1];
	return (r>>1)&0x7fffffff; /* discard the least-random bit */
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(Int4 *idum)
{
        Int4 j;
        Int4 k;
        static Int4 iy=0;
        static Int4 iv[NTAB];
        double temp;

        if (*idum <= 0 || !iy) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ;
                        *idum=IA*(*idum-k*IQ)-IR*k;
                        if (*idum < 0) *idum += IM;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = *idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

#define PI 3.141592654

double bnldev(double pp, Int4 n, Int4 *idum)
/*********************************************************************
 Returns as a double number an integer value that is a random deviate
 drawn from a binomial distribution of n trials, each of probability
 pp, using ran1(idum) as a source of uniform deviates.
/*********************************************************************/
{
	Int4 j;
	static Int4 nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=lngamma(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-lngamma(em+1.0)
				-lngamma(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

double expdev(Int4 *idum)
{
	double dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

double	DirichletDev(Int4 *a, double *p, Int4 dim, Int4 *seed)
// generalization of betadev( )
{
	double	total=0.0;
	Int4	d;
	for(d=1; d <= dim; d++){
	   if(a[d] > 0){ p[d]=gamdev(a[d],seed); total+=p[d]; }
	   else p[d]=0.0;
	}
	for(d=1; d <= dim; d++){ p[d] = p[d]/total; }
	return total;
}

double	betadev(Int4 a, Int4 b, Int4 *seed)
/***************************************************************************
 See Knuth, vol. 2, p. 130 for a published routine to get around 
 nr copyright.  (Cheng ... CACM (1978?)...)
/***************************************************************************/
{
	double	X1;

	X1 = gamdev(a,seed);
	return X1/(X1+gamdev(b,seed));
}

double gamdev(Int4 ia, Int4 *idum)
{
	Int4 j;
	double am,e,s,v1,v2,x,y;

	if (ia < 1){
	    fprintf(stderr,"ia = %d\n",ia);
	    assert(ia > 0);
	    nrerror_afn("Error in routine gamdev");
	}
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran1(idum);
		x = -log(x);
	} else {
	   do {
		do {
		   do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
		   } while (v1*v1+v2*v2 > 1.0);
		   y=v2/v1;
		   am=ia-1;
		   s=sqrt(2.0*am+1.0);
		   x=s*y+am;
		} while (x <= 0.0);
		e=(1.0+y*y)*exp(am*log(x/am)-s*y);
	   } while (ran1(idum) > e);
	}
	return x;
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

double gasdev(Int4 *idum)
{
	static Int4 iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

#define PI 3.141592654

double poidev(double xm, Int4 *idum)
{
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-lngamma(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-lngamma(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
