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

#include "karlin.h"

/**************** Statistical Significance Parameter Subroutine ****************

	Version 1.0	February 2, 1990
	Version 2.0	March 18,   1993

	Program by:	Stephen Altschul

	Address:	National Center for Biotechnology Information
			National Library of Medicine
			National Institutes of Health
			Bethesda, MD  20894

	Internet:	altschul@ncbi.nlm.nih.gov

See:	Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
	Significance of Molecular Sequence Features by Using General Scoring
	Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

	Computes the parameters lambda and K for use in calculating the
	statistical significance of high-scoring segments or subalignments.

	The scoring scheme must be integer valued.  A positive score must be
	possible, but the expected (mean) score must be negative.

	A program that calls this routine must provide the value of the lowest
	possible score, the value of the greatest possible score, and a pointer
	to an array of probabilities for the occurence of all scores between
	these two extreme scores.  For example, if score -2 occurs with
	probability 0.7, score 0 occurs with probability 0.1, and score 3
	occurs with probability 0.2, then the subroutine must be called with
	low = -2, high = 3, and pr pointing to the array of values
	{ 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
	pointers to lambda and K; the subroutine will then calculate the values
	of these two parameters.  In this example, lambda=0.330 and K=0.154.

	The parameters lambda and K can be used as follows.  Suppose we are
	given a length N random sequence of independent letters.  Associated
	with each letter is a score, and the probabilities of the letters
	determine the probability for each score.  Let S be the aggregate score
	of the highest scoring contiguous segment of this sequence.  Then if N
	is sufficiently large (greater than 100), the following bound on the
	probability that S is greater than or equal to x applies:
	
		P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].
	
	In other words, the p-value for this segment can be written as
	1-exp[-KN*exp(-lambda*S)].

	This formula can be applied to pairwise sequence comparison by assigning
	scores to pairs of letters (e.g. amino acids), and by replacing N in the
	formula with N*M, where N and M are the lengths of the two sequences
	being compared.

	In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
	distinct segments all with score >= S is given by:

                               2             m-1           -y
		1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

	Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
	as the previous formula.

	Int4 low;		* Lowest score (must be negative)    *
	Int4 high;		* Highest score (must be positive)   *
	double *pr;		* Probabilities for various scores   *
	double *lambda;		* Pointer to parameter lambda        *
	double *K;		* Pointer to parmeter K              *
	double *H;		* Pointer to parmeter H              *
*******************************************************************************/

BooLean	karlin(Int4 low,Int4 high,double *pr,double *lambda,double *K,double *H)
{
	Int4 i,j,range,lo,hi,first,last;
	double up,New,sum,Sum,av,beta,oldsum,ratio,ftemp;
	double *p,*P,*ptrP,*ptr1,*ptr2;

	/* Check that scores and their associated probabilities are valid     */

	if (low>=0) {
	   fprintf(stderr,"Lowest score must be negative.\n");
	   return FALSE;
	}
	for (i=range=high-low;i> -low && !pr[i];--i);
	if (i<= -low) {
	   fprintf(stderr,"A positive score must be possible.\n");
	   return FALSE;
	}
	for (sum=i=0;i<=range;sum+=pr[i++]) if (pr[i]<0) {
	   fprintf(stderr,"Negative probabilities not allowed.\n");
	   return FALSE;
	}
	if (sum<0.99995 || sum>1.00005)
	   fprintf(stderr,"Probabilities sum to %.4f.  Normalizing.\n",sum);
	NEW(p,range+1,double);
	for (Sum=low,i=0;i<=range;++i) Sum+=i*(p[i]=pr[i]/sum);
	if (Sum>=0) {
	   fprintf(stderr,"Invalid (non-negative) expected score:  %.3f\n",Sum);
	   free(p);
	   return FALSE;
	}

	/* Calculate the parameter lambda */

	up=0.5;
	do {
		up*=2;
		ptr1=p;
		beta=exp(up);
		ftemp=exp(up*(low-1));
		for (sum=i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);
	}
	while (sum<1.0);
	for (*lambda=j=0;j<25;++j) {
		New=(*lambda+up)/2.0;
		beta=exp(New);
		ftemp=exp(New*(low-1));
		ptr1=p;
		for (sum=i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);
		if (sum>1.0) up=New;
		else *lambda=New;
	}

	/* Calculate the pamameter K */

	ptr1=p;
	ftemp=exp(*lambda*(low-1));
	for (av=0,i=low;i<=high;++i) av+= *ptr1++ *i*(ftemp*=beta);
        *H= *lambda*av;
	Sum=lo=hi=0;
	NEW(P,KARLINMAXIT*range+1,double);
	for (*P=sum=oldsum=j=1;j<=KARLINMAXIT && sum>0.001;Sum+=sum/=j++) {
		first=last=range;
		for (ptrP=P+(hi+=high)-(lo+=low);ptrP>=P;*ptrP-- =sum) {
			ptr1=ptrP-first;
			ptr2=p+first;
			for (sum=0,i=first;i<=last;++i) sum+= *ptr1-- * *ptr2++;
			if (first) --first;
			if (ptrP-P<=range) --last;
		}
		ftemp=exp(*lambda*(lo-1));
		for (sum=0,i=lo;i;++i) sum+= *++ptrP * (ftemp*=beta);
		for (;i<=hi;++i) sum+= *++ptrP;
		ratio=sum/oldsum;
		oldsum=sum;
	}
	for (;j<=200;Sum+=oldsum/j++) oldsum*=ratio;
	for (i=low;!p[i-low];++i);
	for (j= -i;i<high && j>1;) if (p[++i-low]) j=karlin_gcd(j,i);
	*K = (j*exp(-2*Sum))/(av*(1.0-exp(- *lambda*j)));
	free(p);
	free(P);
	return TRUE;		/* Parameters calculated successfully */
}

double ExpectedInformation(a_type A, double lambda, double *freq)
{
        Int4 i,j;
        double sum,tot,fij,eij,mij;

        sum = tot = 0.0;

        for (i=1; i<=20; i++){
           for (j=1; j<=20; j++) {
                mij = valAlphaR(i,j,A);
                fij = freq[i]*freq[j];
                tot += fij;
                eij = mij*fij*exp(lambda*mij);
                sum += eij;
           }
        }
        return(lambda*sum/tot);
}

Int4	karlin_gcd(Int4 a,Int4 b)
{
	Int4 c;

	if (b<0) b= -b;
	if (b>a) { c=a; a=b; b=c; }
	for (;b;b=c) { c=a%b; a=b; }
	return a;
}

