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

#include "sfp_typ.h" // Taken from: blastkar.c by Tom Madden

Boolean	sfp_typ::CheckScore(Int4 lo, Int4 hi)
{ if(lo < score_min || hi > score_max) return TRUE; else return FALSE; }

static short BlastScoreChk(Int4 lo, Int4 hi)
{
   if (lo>=0 || hi<=0 || lo<GBLAST_SCORE_1MIN || hi>GBLAST_SCORE_1MAX) return 1;
   if (hi - lo > GBLAST_SCORE_RANGE_MAX) return 1;
   return 0;
}

sfp_typ::sfp_typ(Int4 scr_min, Int4 scr_max)
{
	Int4 	range;

	assert(BlastScoreChk(scr_min, scr_max) == 0);
	score_min = scr_min; score_max = scr_max;
	range = score_max - score_min + 1;
	sprob = (double *) GMemNew(sizeof(double) * range);
	assert(sprob != NULL);
	sprob0 = sprob; sprob -= score_min;
	obs_min = obs_max = 0; score_avg = 0.0;
}

sfp_typ::sfp_typ(Int4 scr_min, Int4 scr_max, Int4 dim1, Int4 dim2,
	unsigned char *query,unsigned char Xchar,Int4 **matrix,
	double *standardProb, Int4 effectiveLength)
// constructor for psiblast matrices...
// Xchar = character for low-complexity columns.
{
   Int4 c;  /*index on characters*/
   Int4 p;  /*index on positions*/
   Int4 s;  /*index on scores */
   // Int4 numberOfScores; // number of distinct scores (NOT USED)
   Int4		score;  	// one score in the matrix
   double	increment;	// Increment in probability due to one score
   
   assert(BlastScoreChk(scr_min, scr_max) == 0);
   score_min = scr_min; score_max = scr_max;
   Int4 range = score_max - score_min + 1;
   sprob = (double *) GMemNew(sizeof(double) * range);
   assert(sprob != NULL);
   sprob0 = sprob; sprob -= score_min;
   score_avg = 0.0;
   obs_min = scr_min; obs_max = scr_max; 
   for(score=score_min; score<=score_max; score++) sprob[score] = 0.0;
   // numberOfScores = (score_max) - (score_min) + 1; // never used...
   for (p = 0; p < dim1; p++){
     if (Xchar != query[p]){
       for (c = 0; c < dim2; c++) {
         // Increment the weight for the score in position [p][c].
         score = matrix[p][c];
#if 1	// AFN: the following is my attempt to fix purify error!!...
	if(score > GBLAST_SCORE_MIN){
           increment = (standardProb[c]/ effectiveLength);
           sprob[score]+= increment;
	}
#endif	// ANF: END 

#if 0	// ORIGINAL: 
// std::cerr << score; std::cerr << std::endl;
         increment = (standardProb[c]/ effectiveLength);
	 if(score <= GBLAST_SCORE_MIN) // AFN DEBUG...
	   fprintf(stderr,
		"p = %d(%d); c = %d; score = %d; sprob[score] = %g; increment=%g\n",
			p,query[p],c,score,sprob[score],increment);
         sprob[score]+= increment;
#endif	// ORIGINAL: END 
       }
     }
   }
   score_avg = 0.0;
   for(s = score_min; s <= score_max; s++) score_avg += s*sprob[s];
}


short	sfp_typ::ScoreFreqCalc(Int4 **matrix, short alphabet_end, 
	Int4 loscore, double *rfp1,double *rfp2)
{
        Int4    score;
        double  score_sum;
        short	index1, index2;

	for(score=score_min; score<=score_max; score++) sprob[score] = 0.0;
        for (index1=0; index1<alphabet_end; index1++) {
             for (index2=0; index2<alphabet_end; index2++) {
                score = matrix[index1][index2];
                if (score >= loscore) sprob[score] += rfp1[index1] * rfp2[index2];
             }
        } score_sum = 0.;
        obs_min = obs_max = GBLAST_SCORE_MIN;
        for (score = score_min; score <= score_max; score++) {
                if (sprob[score] > 0.) {
                        score_sum += sprob[score];
                        obs_max = score;
                        if (obs_min == GBLAST_SCORE_MIN) obs_min = score;
                }
        }
        score_avg = 0.0;
        if (score_sum > 0.0001 || score_sum < -0.0001) {
                for (score = obs_min; score <= obs_max; score++) {
                        sprob[score] /= score_sum;
                        score_avg += score * sprob[score];
                }
        }
	return 0;
}

void	sfp_typ::Free() { GMemFree(sprob0); }

//   Everything below here was (more or less) copied from the old
//   karlin.c and could work separately from the stuff above.
/**************** Statistical Significance Parameter Subroutine ****************

    Version 1.0     February 2, 1990
    Version 1.2     July 6,     1990

    Program by:     Stephen Altschul

    Address:        National Center for Biotechnology Information
                    National Library of Medicine
                    National Institutes of Health
                    Bethesda, MD  20894

    Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
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

*******************************************************************************/

#define DIMOFP0	(iter*range + 1)
#define DIMOFP0_MAX (BLAST_KARLIN_K_ITER_MAX*GBLAST_SCORE_RANGE_MAX+1)

double	sfp_typ::KarlinLHtoK(double lambda, double H)
{
#ifndef BLAST_KARLIN_STACKP
	double	*P0 = (double *) NULL;
#else
	double	P0 [DIMOFP0_MAX];
#endif
	Int4 	low;	/* Lowest score (must be negative) */
	Int4 	high;	/* Highest score (must be positive) */
	double	K;			/* local copy of K */
	double	ratio;
	int		i, j;
	Int4 	range, lo, hi, first, last;
	register double	sum;
	double	Sum, av, oldsum, oldsum2;
	int		iter;
	double	sumlimit;
	double	*p, *ptrP, *ptr1, *ptr2, *ptr1e;
	double	etolami, etolam;

	if (lambda <= 0. || H <= 0.) { return -1.; }
	if (score_avg >= 0.0) { return -1.; }
	low = obs_min;
	high = obs_max;
	range = high - low;
	av = H/lambda;
	etolam = exp((double)lambda);
	if (low == -1 || high == 1) {
		if (high == 1) K = av;
		else K = (score_avg * score_avg) / av;
		return K * (1.0 - 1./etolam);
	}
	sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
	iter = BLAST_KARLIN_K_ITER_MAX;
	if (DIMOFP0 > DIMOFP0_MAX) { return -1.; }
#ifndef BLAST_KARLIN_STACKP
	P0 = (double *)GMemNew(DIMOFP0 * sizeof(*P0));
	if (P0 == NULL) return -1.;
#else
	Nlm_MemSet((CharPtr)P0, 0, DIMOFP0*sizeof(P0[0]));
	// this is turned off...
#endif
	Sum = 0.;
	lo = hi = 0;
	p = &sprob[low];
	P0[0] = sum = oldsum = oldsum2 = 1.;
    for (j = 0; j < iter && sum > sumlimit; Sum += sum /= ++j) {
        first = last = range;
		lo += low;
		hi += high;
        for (ptrP = P0+(hi-lo); ptrP >= P0; *ptrP-- =sum) {
            ptr1 = ptrP - first;
            ptr1e = ptrP - last;
            ptr2 = p + first;
            for (sum = 0.; ptr1 >= ptr1e; ) sum += *ptr1--  *  *ptr2++;
            if (first) --first;
            if (ptrP - P0 <= range) --last;
        }
		etolami = GNlm_Powi((double)etolam, lo - 1);
        for (sum = 0., i = lo; i != 0; ++i)
		{ etolami *= etolam; sum += *++ptrP * etolami; }
        for (; i <= hi; ++i) sum += *++ptrP;
	oldsum2 = oldsum; oldsum = sum;
    }
	/* Terms of geometric progression added for correction */
	ratio = oldsum / oldsum2;
	if (ratio >= (1.0 - sumlimit*0.001)) { K = -1.; goto CleanUp; }
	sumlimit *= 0.01;
	while (sum > sumlimit) { oldsum *= ratio; Sum += sum = oldsum / ++j; }

/* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
Karlin&Altschul (1990) */
    	for (i = 1, j = -low; i <= range && j > 1; ++i)
        	if (p[i]) j = GNlm_Gcd(j, i);

	if (j*etolam > 0.05) {
		etolami = GNlm_Powi((double)etolam, -j);
    		K = j*exp((double)-2.0*Sum) / (av*(1.0 - etolami));
	} else
	    K = -j*exp((double)-2.0*Sum) / (av*GNlm_Expm1(-j*(double)lambda));

CleanUp:
#ifndef BLAST_KARLIN_K_STACKP
	if (P0 != NULL) GMemFree(P0);
#endif
	return K;
}

double	sfp_typ::KarlinLambdaBis( )
// BlastKarlinLambdaBis: Calculate Lambda using the bisection method (slow).
{
	double	lambda, up, newval;
	Int4 	i, low, high;
	int		j;
	register double	sum, x0, x1;

	if (score_avg >= 0.) { return -1.; }
	low = obs_min; high = obs_max;
	if (BlastScoreChk(low, high) != 0) return -1.;
	up = BLAST_KARLIN_LAMBDA0_DEFAULT;
	for (lambda=0.; ; ) {
		up *= 2;
		x0 = exp((double)up);
		x1 = GNlm_Powi((double)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		} else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(up * i);
		} if (sum >= 1.0) break;
		lambda = up;
	}
	for (j=0; j<BLAST_KARLIN_LAMBDA_ITER_DEFAULT; ++j) {
		newval = (lambda + up) / 2.;
		x0 = exp((double)newval);
		x1 = GNlm_Powi((double)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		} else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(newval * i);
		}
		if (sum > 1.0) up = newval;
		else lambda = newval;
	}
	return (lambda + up) / 2.;
}

/******************* Fast Lambda Calculation Subroutine ************************
	Version 1.0	May 16, 1991
	Program by:	Stephen Altschul

	Uses Newton-Raphson method (fast) to solve for Lambda, given an initial
	guess (lambda0) obtained perhaps by the bisection method.
*******************************************************************************/
double	sfp_typ::KarlinLambdaNR( )
{
	Int4 	low;		// Lowest score (must be negative) 
	Int4 	high;		// Highest score (must be positive)
	int	j;
	Int4 	i;
	double	lambda0, sum, slope, temp, x0, x1, amt;

	low = obs_min; high = obs_max;
	if (score_avg >= 0.) {	/* Expected score must be negative */
		return -1.0;
	}
	if (BlastScoreChk(low, high) != 0) return -1.;
	lambda0 = BLAST_KARLIN_LAMBDA0_DEFAULT;
	/* Calculate lambda */
	for (j=0; j<20; ++j) { /* limit of 20 should never be close-approached */
		sum = -1.0;
		slope = 0.0;
		if (lambda0 < 0.01) break;
		x0 = exp((double)lambda0);
		x1 = GNlm_Powi((double)x0, low - 1);
		if (x1 == 0.) break;
		for (i=low; i<=high; ++i) {
			sum += (temp = sprob[i] * (x1 *= x0));
			slope += temp * i;
		} lambda0 -= (amt = sum/slope);
		if (ABS(amt/lambda0) < BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT) {
			/*
			Does it appear that we may be on the verge of converging
			to the ever-present, zero-valued solution?
			*/
			if (lambda0 > BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT)
				return lambda0;
			break;
		}
	}
	return KarlinLambdaBis( );
}

double	sfp_typ::KarlinLtoH(double lambda)
// BlastKarlinLtoH Calculate H, the relative entropy of the p's and q's.
{
	Int4 	score;
	double	av, etolam, etolami;

	if (lambda < 0.) { return -1.; }
	if (BlastScoreChk(obs_min, obs_max) != 0) return -1.;
	etolam = exp((double)lambda);
	etolami = GNlm_Powi((double)etolam, obs_min - 1);
	if (etolami > 0.) {
	    av = 0.0;
	    for (score=obs_min; score<=obs_max; score++)
   			av += sprob[score] * score * (etolami *= etolam);
	} else {
	    av = 0.0;
	    for (score=obs_min; score<=obs_max; score++)
   			av += sprob[score] * score * exp(lambda * score);
	}
    	return lambda * av;
}

