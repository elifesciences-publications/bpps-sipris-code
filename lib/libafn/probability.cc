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

#include "probability.h"

#if 0
/*************** lngamma function adapted from Spouge *****************
** Special thanks to Dr. John ``Gammas Galore'' Spouge for deriving the
** method for calculating the gamma coefficients and their use.

 *  For discussion of the Gamma function, see "Numerical Recipes in C",
    Press et al. (1988) pages 167-169.  
 ***********************************************************************/

static double	lnpI_plus_ln2_over2 = ((MATHNCBI_LNPI+MATHNCBI_LN2) / 2.);
static double	lnpI_plus_ln2_over2X = (((MATHNCBI_LNPI+MATHNCBI_LN2) / 2.) - 10.5);

double	SpougeLnGamma(register double xx)
/* ln(ABSOLUTE[gamma(x)]) - 10 digits of accuracy */
/* ln[gamma(x)]; x is 10-digit accuracy achieved only for 1 <= x */
{
	register double	value,tx;

	if(xx >= 1.) {
		/* sum the least significant terms first */
		xx += 10.;  /** dimensions of gamma = 11 **/
		/* sum the least significant terms first */
		value = 1.251639670050933e-10/xx; xx -= 1.;
		value += -2.319827630494973e-04/xx; xx -= 1.;
		value += 2.908143421162229e-01/xx; xx -= 1.;
		value += -3.155153906098611e+01/xx; xx -= 1.;
		value += 8.785855930895250e+02/xx; xx -= 1.;
		value += -9.601592329182778e+03/xx; xx -= 1.;
		value += 5.031796415085709e+04/xx; xx -= 1.;
		value += -1.388934775095388e+05/xx; xx -= 1.;
		value += 2.065049568014106e+05/xx; xx -= 1.;
		value += -1.560605207784446e+05/xx; xx -= 1.;
		value = log(value + 1. + 4.694580336184385e+04/xx); 
		value += lnpI_plus_ln2_over2X + (xx - 0.5)*log(xx+10.5) - xx;
		return value;
	} else if(xx > 0.){
		tx = log(xx);
		xx += 11.;  /** dimensions of gamma = 11 (input xx + 1.) **/
		/* sum the least significant terms first */
		value = 1.251639670050933e-10/xx; xx -= 1.;
		value += -2.319827630494973e-04/xx; xx -= 1.;
		value += 2.908143421162229e-01/xx; xx -= 1.;
		value += -3.155153906098611e+01/xx; xx -= 1.;
		value += 8.785855930895250e+02/xx; xx -= 1.;
		value += -9.601592329182778e+03/xx; xx -= 1.;
		value += 5.031796415085709e+04/xx; xx -= 1.;
		value += -1.388934775095388e+05/xx; xx -= 1.;
		value += 2.065049568014106e+05/xx; xx -= 1.;
		value += -1.560605207784446e+05/xx; xx -= 1.;
		value = log(value + 1. + 4.694580336184385e+04/xx); 
		value += lnpI_plus_ln2_over2X + (xx - 0.5)*log(xx+10.5) - xx;
		return (value - tx);
	} else print_error("Input error in SpougeLnGamma( )");
}

/********************* end lngamma from Spouge **************************/
#endif

double	LogSumInvLogs(register double z, register double x)
/**************************************************************
  for z >= x the following is true

  exp[z + log(1+exp(x-z))] = exp(z)*(1+exp(x-z)) = exp(z) + exp(x)

 so this can compute log[exp(z) + exp(x)] without loss of precision.
 **************************************************************/
{
        if(z >= x) return (z + log(1+exp(x-z)));
        return (x + log(1+exp(z-x)));
}

double  WtdCumBnmlKLDivergence(double N, double n, double M, double m, double b1, double b2)
// compute the weighted Cumulative Binomial Kullback–Leibler divergence in nats.
// p = (n+1)/(N+2); q = (m+1)/(M+2);
// b1 and b2 are pseudocounts for n and not n,respectively.
{

        double fg,bg,P,p,q;

	q = (n+b1)/(N+b1+b2); p = (m+b1)/(M+b1+b2);
	if(p != q){
	 fg=CumBinomProb(n,N,q)*(LnCBP(n,N,q)-LnCBP(n,N,p)); // selection for
         fg += CumBinomProb(N-n,N,1-q)*(LnCBP(N-n,N,1-q)-LnCBP(N-n,N,1-p)); // selection against
	 bg=CumBinomProb(m,M,p)*(LnCBP(m,M,p)-LnCBP(m,M,q)); // selection for
         bg += CumBinomProb(M-m,M,1-p)*(LnCBP(M-m,M,1-p)-LnCBP(M-m,M,1-q)); // selection against
	 P = (M/(M+N))*fg + (N/(M+N))*bg;
	} else P=0.0; 
	if(q < p) return -P; else return P;
}

double  LnBinomialProb(double N, double p, double n)
{
	return (lngamma(N+1)-lngamma(n+1) - lngamma(N-n+1)
                                          + n*log(p) + (N-n)*log(1.0 - p));
}

double CumBinomProb(double k, double N, double x)
{
	if(k == 0.0) return 1.0;
	return exp(LnCBP(k, N, x)); 
}

double Log10CBP(double k, double N, double x)
{ return LnCBP(k, N, x)/NATURAL_LOG_FACTOR; }

double LnCBP(double a, double N, double x)
/***********************************************************************
  Calculation of Log base 10 of the cumulative binomial probability (CBP)
  using the continued fraction expression for the incomplete Beta function.
                    a       b                                              
	           x (1 - x)   / 1   d[1] d[2]     \ 
        I (a,b) = ----------- | ---- ---- ----  ... |
         x         a B(a,b)    \ 1+   1+   1+      / 

                         - (a + m)(a + b + m) x
	where d[2m+1] =  -----------------------
                          (a + 2m)(a + 2m + 1)

                           m (b - m) x
	where d[2m] =  -------------------    (see 26.5.8, p.944 in ref. 1)
                       (a + 2m -1)(a + 2m)

  (The nth convergent of this continued fraction is computed by making
  use of the theorems given in 3.10.1 (p.19) of ref. 1.)

  This is based on the following relationship: 

		I (a,N-a+1) = CBP(N,a,x).   (see 26.5.24, p.945 in ref. 1)
                 x
                                       a - 1
  Best results are obtained when x < --------- (ref. 1) 
                                     a + b - 2

                                  a + 1
	or alternately when x < ---------  (ref. 2).
                                a + b + 2

  Thus the symmetry relationship 

	I (a,b) = I (b,a)		(see 26.5.2, p.944 in ref. 1)
         x         x

  is used to obtain optimum performance.

  Parameters given in ref. 2 were used.
  
  References: 
   1. Abramowitz & Stegun, "Handbook of Mathematical Functions",
      1964 - Reprinted by Dover Press 1972; pp.944,.

   2. Press, Flannery, Teukolsky & Vetterling, "Numerical Recipes 
      in C", 1988, Cambridge University Press pp.178-180.

   normalization constant = ln[x^a(1-x)^b/aB(a,b)]  

	note: V = a + 2*m;
 ***********************************************************************/
{
	register double m,V,d;
	double C,a1,b,b1,b2,fn,old;
	BooLean	flag=TRUE;

  if(x > 0.0 && x < 1.0){
  	if(a==0.0) return 0.0;  // fix error: CBP = 1 --> lnCBP = 0; afn: 11/7/07
	if(x < (a+1.0)/(N+3.0)) b= N-a+1;
	else {flag=FALSE;b=a;a=N-a+1;x=1.0-x;}
	old=m=1.0; C=a+b; V=a+2.0;
	fn=1.0-C*x/(a+1.0);
	d=(b-1.0)*x/((V-1.0)*(V));  /** d[2] **/
	a1=1.0+d; b1=fn+d;
	d= -(a+1)*(C+1)*x/((V+1)*(V)); /** d[3] **/
	b2=b1+d*fn;
	fn=(a1+d)/b2;
	while(fabs(fn-old) >= (3.0e-7*fabs(fn))){
		if((m+=1.0) > 100.0) return ILLEGAL;
		V+=2;
		d=m*(b-m)*x/((V-1.0)*(V));  /** d[2m] **/
		a1=fn+d*(a1/b2);
		b1=1.0+d*(b1/b2);
		d= -(a+m)*(C+m)*x/((V+1.0)*(V)); /* d[2m+1] */
		b2=b1+d;
		old=fn; fn=(a1+d*fn)/b2; 
	}
	if(flag) return (lngamma(a+b)-lngamma(a)-lngamma(b)
			+a*log(x)+b*log(1.0-x) + log(fn/a));
	return (log(1.0-(fn/a)*
		exp(lngamma(a+b)-lngamma(a)-lngamma(b)+a*log(x)+b*log(1.0-x))));
  } else if(x < 0.0 || x > 1.0) return ILLEGAL;
  else if(x ==0.0){ if(a== 0.0) return 0.0; else return ILLEGAL; }
  else {if(a==N) return 0.0; else return ILLEGAL; }
}

double  HypGeomProb(Int4 N1,Int4 N2,Int4 n,Int4 x)
{
	double	p;
        p = (lnfact(N1)-lnfact(x)-lnfact(N1-x));
        p += (lnfact(N2)-lnfact(n-x)-lnfact(N2-n+x));
        p -= (lnfact(N2+N1)-lnfact(n)-lnfact(N2+N1-n));
        return exp(p);
}

double  CumHypGeomProb(Int4 N1,Int4 N2, Int4 n,Int4 x)
#if 0   //****************************************************
N total balls with N1 red balls and N2 = N-N1 black balls.
Choose n balls at random.  The probability that the
group so chosen will contain x or more red balls is
given by: p=CumHypGeomProb(N1,N2,n,x).
#endif  //****************************************************
{
        Int4	end;
        double	p,K;

	if(x == 0) return 1.0;
        end = MINIMUM(Int4,N1,n);
        K = (lnfact(N1)+lnfact(N2)-lnfact(N2+N1)+lnfact(n)+lnfact(N2+N1-n));
        for(p=0.0,x; x<=end; x++){
           p += exp(K-lnfact(x)-lnfact(N1-x)-lnfact(n-x)-lnfact(N2-n+x));
        }
        return p;
}

double  lnbico(register Int4 N, register Int4 k)
{ return lnfact(N)-lnfact(k)-lnfact(N-k); }

double  bico(Int4 N,Int4 k)
{ return floor(0.5+exp(lnfact(N)-lnfact(k)-lnfact(N-k))); }

double  factrl(register int n)
/* caution: assumes  n >= 0; 
   if (n < 0) print_error("illegal input to factrl( )"); ****/
{
        static int ntop=4;
        static double a[33]={1.0,1.0,2.0,6.0,24.0};
        register int j;

        if (n > 32) return exp(lngamma(n+1.0));
        while (ntop<n) {
                j=ntop++;
                a[ntop]=a[j]*ntop;
        }
        return a[n];
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

double	FisherExactTest(FILE *fp, Int4 red_in, Int4 black_in, Int4 red_out, Int4 black_out,
		double &onetail)
/*************** Fisher's Exact Test for a 2x2 contingency table *****
 *  observed:
 *
 *             red    black
 *      ----+-------+-------+-------
 *      out |   a   |   b   | a + b
 *      ----+-------+-------+-------
 *      in  |   c   |   d   | c + d
 *      ----+-------+-------+-------
 *          | a + c | b + d |a+b+c+d
 *
 *      see p. 81 of my thesis.
 *
 **********************************************************************/
{ 
	Int4	arg,s,time1;
	Int4	red,black,in,out;

	double	K,P,Q,end,p,OneTail,TwoTail,total;
	Int4	i,a,b,c,d;
	a = red_out; b = black_out;
	c = red_in; d = black_in;
	red = a+c; black = b+d; out = a+b; in = c+d;
	end = MINIMUM(double,a,d);
        K = lngamma(out+1.0) + lngamma(red+1.0)
                + lngamma(black+1.0) + lngamma(in+1.0) - lngamma(a+b+c+d+1.0);
        P = lngamma(a+1.0) + lngamma(b+1.0) + lngamma(c+1.0) + lngamma(d+1.0);
	OneTail=TwoTail=total=0.0;
        for(i=MAXIMUM(Int4,-b,-c); i <= end; i = i + 1.0) {
            Q = lngamma(a-i+1.0) + lngamma(b+i+1.0)
                        + lngamma(c+i+1.0) + lngamma(d-i+1.0);
	    p = exp(K-Q);
            if(fp) fprintf(fp," [%d : %d | %d : %d] (Q = %g)\n", 
			(Int4)(a-i), (Int4) (b+i),(Int4) (c+i), (Int4) (d-i),p);
            if(P <= Q){
		TwoTail += p;
		if((a-i) >= a) OneTail += p;
	    } total+=p;
        }
	if(fp){
	  fprintf(stderr,"Drew %d balls out of an urn with %d red and %d black balls\n",
			out,red,black);
	  fprintf(stderr,"Among these there were %d red balls.\n",red_out);
	  fprintf(stderr,"One tail prob = %g\n",OneTail);
	  fprintf(stderr,"Two tail prob = %g\n",TwoTail);
	  fprintf(stderr,"total = %.3f\n",total);
	} onetail=OneTail; 
	return TwoTail;
}

double	lnfact(Int4 n)
/* static variables are guaranteed to be initialized to zero */
{
	static double lnft[101];

	if (n <= 1) return 0.0;
	if (n <= 100) return lnft[n] ? lnft[n] : (lnft[n]=lngamma(n+1.0));
	else return lngamma(n+1.0);
}

/***** numerical recipes code - don't give away!!! **********/
double ErfC(double x)
{ return x < 0.0 ? ((1.0 + gammp(0.5,x*x))/2.0) : (gammq(0.5,x*x)/2.0); }

double gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror_afn("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
double gammq(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror_afn("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	Int4 n;
	double sum,del,ap;

	*gln=lngamma(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror_afn("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror_afn("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	void nrerror_afn(char error_text[]);
	Int4 i;
	double an,b,c,d,del,h;

	*gln=lngamma(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror_afn("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

#define ITMAX   100
#define EPS     3.0e-7

Int4	gser2(double *gamser,double a,double x,double *gln)
/* Returns the incomplete gamma function P(a,x) evaluated by its
series representation as gamser.  Also returns lnGamma(a) as gln. */
{
	Int4 n;
	double	sum, del,ap;

	*gln=lngamma(a);
	if(x <= 0.0) {
		if(x < 0.0) nrerror_afn("x less than 0 in routine GSER");
		*gamser=0.0;
		return TRUE;
	} else {
		ap=a;
		del=sum=1.0/a;
		for(n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if(fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return (TRUE);
			}
		}
		return (FALSE);
	}
}

Int4	gcf2(double *gammcf,double a,double x,double *gln)
/* Returns the incomplete gamma function Q(a,x) evaluated by its
continued fraction representation as gammcf.  Also returns Gamma(a) as
gln. */
{
	Int4 n;
	double	gold=0.0,g,fac=1.0,b1=1.0;	/* fac = the renormalization */
	double	b0=0.0,anf,ana,an,a1,a0=1.0;	/* factor for preventing */

	*gln=lngamma(a);
	a1=x;	/* We are here setting up the A's and B's of eqn(5.2.4)  */
	for(n=1;n<=ITMAX;n++) {   /* for evaluating the continued fraction */
		an=(double)n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;	/* One step of the recurrence (5.2.5)*/
		b0=(b1+b0*ana)*fac;	
		anf=an*fac;
		a1=x*a0+anf*a1;		/* Next step of the recurrence */
		b1=x*b0+anf*b1;		
		if(a1) {		/* Shall we renormalize? */
			fac=1.0/a1;	/* Yes. Set fac so it happens */
			g=b1*fac;	/* New value of answer */
			if(fabs((g-gold)/g) < EPS) {  /* Converge? If so exit*/
				*gammcf=exp(-x+a*log(x)-(*gln))*g; 
					/* put factors in front */
				return(TRUE);
			}
			gold=g;		/* if not, save value. */
		}
	}
	return (FALSE);
}

void	nrerror_afn(char *error_text)
/* Numerical recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

#include "nrutil.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(UInt4 n, double arr[])
{
	UInt4 i,ir=n,j,k,l=1;
	Int4	jstack=0,*istack;
	double	a,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				} arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) { SWAP(arr[l+1],arr[ir]) }
			if (arr[l] > arr[ir]) { SWAP(arr[l],arr[ir]) }
			if (arr[l+1] > arr[l]) { SWAP(arr[l+1],arr[l]) }
			i=l+1; j=ir; a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) { istack[jstack]=ir; istack[jstack-1]=i; ir=j-1; }
			else { istack[jstack]=j-1; istack[jstack-1]=l; l=i; }
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

#define EPS1 0.001
#define EPS2 1.0e-8

double probks(double alam)
{
	Int4 j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for (j=1;j<=100;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}
#undef EPS1
#undef EPS2
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

void kstwo(double data1[], UInt4 n1, double data2[], UInt4 n2,
	double *d, double *prob)
{
	UInt4 j1=1,j2=1;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	sort(n1,data1);
	sort(n2,data2);
	en1=n1;
	en2=n2;
	*d=0.0;
	while (j1 <= n1 && j2 <= n2) {
		if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
		if (d2 <= d1) fn2=j2++/en2;
		if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
	}
	en=sqrt(en1*en2/(en1+en2));
	*prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */


