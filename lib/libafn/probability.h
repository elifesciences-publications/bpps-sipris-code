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

#if !defined(PROBABILITY)
#define PROBABILITY
#include "stdinc.h"

double  lnfact(Int4 n);
double  factrl(register int n);
double  CumHypGeomProb(Int4 N1,Int4 N2, Int4 n,Int4 x);
double  HypGeomProb(Int4 N1,Int4 N2,Int4 n,Int4 x);
double Log10CBP(double k, double N, double x);
double CumBinomProb(double k, double N, double x);
double LnCBP(double a, double N, double x);
double  lnbico(register Int4 N, register Int4 k);
double  bico(Int4 N,Int4 k);
double  LnBinomialProb(double N, double p, double n);
double  WtdCumBnmlKLDivergence(double N, double n, double M, double m, double b1, double b2);
double  FisherExactTest(FILE *fp, Int4 red_in, Int4 black_in, Int4 red_out, Int4 black_out,
                double &onetail);
#if 0
double  SpougeLnGamma(register double xx);
#endif

double probks(double alam);
void kstwo(double data1[], UInt4 n1, double data2[], UInt4 n2,
        double *d, double *prob);

/*********************************************************/
double  LogSumInvLogs(register double z, register double x);
double ErfC(double x);
Int4	gser2(double *gamser,double a,double x,double *gln);
Int4	gcf2(double *gammcf,double a,double x,double *gln);

void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammp(double a, double x);
double gammq(double a, double x);

void	nrerror_afn(char *error_text);
/*********************************************************/

#define NATURAL_LOG_FACTOR 2.3025850929940459
#define LnBeta(a,b)             (lngamma((double)(a))+lngamma((double)(b))\
                                        - lngamma((double)((a)+(b))))

#endif

