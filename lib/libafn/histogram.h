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

/* histogram.h - histogram data type */

#if !defined(HIST)
#define HIST
#include "stdinc.h"

/**********************************************************************

          min                                             max
         start   o----> <inc>            o---->o---->   end 
Int4 |_____|_____|_____|_____|___ ... ___|_____|_____|___:_|_____|
                                      
        0     1     2     3                           nbins  over
     (under)
                                                    round up end to max
                                               last bin includes overflow

 *************************** histogram type ***************************/
typedef struct {
	char	*id;		/* identifier string for histogram */
	Int4	*bin0;		/* array for histogram bin0ribution */
	Int4	*bin;		/* bins for histogram bin0ribution */
	Int4	nbins;		/* number of bin0inct values in histogram */
	double	maxval,minval;	/* maximum and minimum values */
	double	max;		/* maximum value */
	double	total;		/* total for mean */
	double	total_sq;	/* total X^2 for variance */
	double	min;		/* minimum value */
	double	inc;		/* bin0ance between values in histogram */
	Int4	n;		/* number of variables in bin0ribution */
} histogram_type;

typedef	histogram_type	*h_type;

/**************************** Public **********************************/
h_type  Histogram(const char *id,Int4 start,Int4 end,double inc);
Int4    IncdHistBin(double x,h_type H);
void	IncdHist(double x,h_type H);
void    IncdMHist(double x, Int4 number, h_type H);
void    PutHistX(FILE *fp, char *str, double inc, Int4 maxLength, Int4 *dist);
void	PutHist(FILE *fptr,Int4 line_leng,h_type H);
Int4	*RtnHist(Int4 line_leng,h_type H);
Int4    RtnHistArray(double **xval, double **yval, Int4 *Total, h_type H);
void    PutHistArray(FILE *fptr,h_type H);
void    NilHist(h_type H);
double  MeanHist(h_type H);
double  MedianHist(h_type H);
double  VarianceHist(h_type H);
Int4    TotalDataHist(h_type H);

#define indexHist(x,H)		((Int4)floor(((double)(x)-(H)->min)/(H)->inc))
#define IncfHist(x,H)		IncdHist((double)(x),H)
#define TotalHist(H)		((H)->total)
#define NumDataHist(H)		((H)->n)
#define MinValHist(H)		((H)->minval)
#define MinimumHist(H)		((H)->min)
#define MaximumHist(H)		((H)->max)
#define MaxValHist(H)		((H)->maxval)
#define NumBinsHist(H)		((H)->nbins+2)

#endif

