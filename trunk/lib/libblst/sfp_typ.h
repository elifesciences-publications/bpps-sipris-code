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

// sfp_typ: Adapted from blastkar.h by Tom Madden
#ifndef __SFP_TYPE__
#define __SFP_TYPE__

#include "stdinc.h"
#include "my_ncbimath.h"

#define GBLAST_SCORE_MIN INT2_MIN
#define GBLAST_SCORE_MAX INT2_MAX

#define GBLAST_SCORE_1MIN (-10000)
#define GBLAST_SCORE_1MAX ( 1000)

#define GBLAST_SCORE_RANGE_MAX   (GBLAST_SCORE_1MAX - GBLAST_SCORE_1MIN)

/****************************************************************************
For more accuracy in the calculation of K, set K_SUMLIMIT to 0.00001.
For high speed in the calculation of K, use a K_SUMLIMIT of 0.001
Note:  statistical significance is often not greatly affected by the value
of K, so high accuracy is generally unwarranted.
*****************************************************************************/
// K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK()
#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.01

// LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd
#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5)

// LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)
#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17

// Initial guess for the value of Lambda in BlastKarlinLambdaNR
#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5

#define BLAST_KARLIN_K_ITER_MAX 100
#define GBLAST_SUMP_EPSILON_DEFAULT 0.002 /* accuracy for SumP calculations */

class sfp_typ {		// _blast_score_freq pointer = BLAST_ScoreFreq, *sfp_typ; 
public:
		sfp_typ(Int4,Int4);
		sfp_typ(Int4,Int4,Int4,Int4,unsigned char *,
			unsigned char,Int4 **,double *, Int4); // for psiblast matrix...
		~sfp_typ(){ Free( ); }
	double	KarlinLHtoK(double lambda, double H);
	double	KarlinLtoH(double lambda);
	double	KarlinLambdaBis( );
	double	KarlinLambdaNR( );
	Int4	ScoreMin( ){ return score_min; }
	Int4	ScoreMax( ){ return score_max; }
	short   ScoreFreqCalc(Int4 **, short , Int4 , double *,double *);
	Boolean	CheckScore(Int4 lo, Int4 hi);
	void	SetObs(Int4 omin, Int4 omax){ obs_min=omin; obs_max=omax; }
	
private:
	void	Free();
	Int4	score_min, score_max;
	Int4	obs_min, obs_max;
	double	score_avg;
	double	*sprob0, *sprob;
};

#endif /* !__BLASTKAR__ */

