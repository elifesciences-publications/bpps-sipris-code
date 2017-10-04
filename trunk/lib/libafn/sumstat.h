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

/* sumstat.h - */
#if !defined (SUMSTAT)
#define SUMSTAT

/******************* High-scoring Segment P-value Subroutine *******************

	Version 1.0	October 15, 1993

	Program by:	Warren Gish, John Spouge & Stephen Altschul

See:	Karlin, S. & Altschul, S.F. (1993). "Applications and Statistics for
  Multiple High-scoring Segments in Molecular Sequences." Proc. Natl.
	Acad. Sci. USA 90:5873-5877.

	Internet:	altschul@nih.gov

	The subroutine SumPStd(r,s) calculates the probability that the sum
	of the normalized scores of the r highest-scoring segments is >= s.
	r is an integer >= 1, and s is a double precision real number.

	Any appropriate scoring system has an assocaited pair of parameters
	lambda & K [Karlin & Altschul (1990) PNAS 87:2264-2268]. The norm-
	alized score of a segment is lambda times its score minus ln(KN),
	where N is the size of the search space. For pairwise sequence
	comparison, N is the product of the lengths of the two sequences
	compared.

*****************************************************************************/
#include "stdinc.h"

double SumPCalc(Int4 r, double s);
double SumPStd(Int4 r, double s);
#if 0
double RombergIntegrate(double (*f)(register double, void *), void *fargs,
	double p, double q, double eps, Int4 epsit, Int4 itmin);
static double	 funct_f(register double, void *);
static double 	 funct_g(register double, void *);
#endif

#define	sump_epsilon	0.002

#endif

