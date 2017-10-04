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

#include "binomial.h"

bn_type  MakeBinomial(Int4 N, Int4 k, char *name)
{
	bn_type	B;

	NEW(B,1,binomial_type);
	if((B->seed=(Int4)Random()) >= 0) B->seed *= -1;
	B->name = AllocString(name); B->N = N;
	SetBinomial(k, B);
	return B;
}

void    NilBinomial(bn_type B) { free(B->name); free(B); }

void	SetBinomial(Int4 x, bn_type B)
{	B->k = x; B->p = (double) B->k/(double)B->N; }

void	PutBinomial(FILE *fp, Int4 size, bn_type B)
{
	h_type  H;

	H = Histogram(B->name,0,B->N,1.0);
	while(size > 0){ 
		IncdHist(SampleBinomial(0,1000,B),H); size--; 
	}
	PutHist(fp,60,H); NilHist(H);
	fprintf(stderr,"mean = %d; total = %d; p = %g\n",
		B->k,B->N,B->p);
}

Int4	SampleBinomial(Int4 min, Int4 max, bn_type B)
{
	Int4	x,i=0;

	if(max < min) print_error("SampleBinomial( ) input error");
	do {
		if(i++ > 500){
			fprintf(stderr,"SampleBinomial( ) possible infinite loop\n");
			fprintf(stderr,"  min=%d; max=%d\n",min,max);
			if(min == max) return min;
			x = random_integer(1 + (max-min));
			return min + x;
		}
		x = (Int4) bnldev(B->p, B->N, &B->seed); 
	} while(x < min || x > max);
	return x;
}

