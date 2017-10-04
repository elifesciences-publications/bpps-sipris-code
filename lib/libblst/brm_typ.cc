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

#include "brm_typ.h"  // adapted from posit.c by Alejandro Schaffer

// Allocates and fills the brm_typ *. 
brm_typ::brm_typ(Int4 alpha_size, Int4 qlength, unsigned char *q,double *stdProb,
	Int4 **mtx, Int4 **private_mtx, kbp_typ *kbpstd, kbp_typ *kbppsi,
	kbp_typ *kbpgapstd, kbp_typ *kbpgappsi,double lambda, double K)
{
	alphabet_size = alpha_size;
	query_length = qlength; query = q;
	standardProb = stdProb;
	matrix = mtx; private_matrix = private_mtx;
	kbp_std = kbpstd; kbp_psi = kbppsi;
	kbp_gap_std = kbpgapstd; kbp_gap_psi = kbpgappsi;
	lambda_ideal = lambda; K_ideal = K;
}

// sfp_typ *brm_typ::ComputeProbs(Boolean position_dependent)
sfp_typ *brm_typ::ComputeProbs( )
// Compute probabilities for each score in posMatrix, also sets minScore and maxScore
// AFN: Can set query_length to alphabet_size to rescale regular matrix
{
   Int4		c;  		// index on characters
   Int4		p;  		// index on positions
   Int4		s;  		// index on scores 
   Int4		score_min, score_max;
   Int4		score;  // one score in the matrix*/
   double	increment;  /*Increment in probability due to one score*/
   Int4		effectiveLength;

   score_min = 0; score_max = 0; effectiveLength = 0;
   for (p = 0; p < query_length; p++) if (Xchar != query[p]) effectiveLength++;
   for (p = 0; p < query_length; p++){
     if (Xchar != query[p])
       for (c = 0; c < alphabet_size; c++) {
// std::cerr << matrix[p][c]; std::cerr << std::endl; 
	 if (matrix[p][c] <= GBLAST_SCORE_MIN || matrix[p][c] >= GBLAST_SCORE_MAX) continue;
	 if (matrix[p][c] < (score_min)) (score_min) = matrix[p][c];
	 if (matrix[p][c] > (score_max)) (score_max) = matrix[p][c];
       }
   }
   sfp_typ *sfp=new sfp_typ(score_min,score_max,query_length,
		alphabet_size,query,Xchar,matrix,standardProb,effectiveLength);
   return(sfp);
}

// void	brm_typ::updateLambdaK(Boolean position_dependent) // always true.
void	brm_typ::updateLambdaK( )
{
  sfp_typ *sfp= ComputeProbs( );
  GBlastKarlinBlkCalc(kbp_psi[0], sfp);
  kbp_gap_psi[0]->K = (kbp_psi[0]->K)*(kbp_gap_std[0]->K)/K_ideal;
  kbp_gap_psi[0]->logK = log(kbp_gap_psi[0]->K);
  delete sfp;
}

double	brm_typ::Scale( )
// former Boolean position_dependent always True.
{
   Int4		a,c;				// loop indices
   Int4		index;				// loop index for binary search
   Boolean	too_high=TRUE, done, first_time;
   double	factor, factor_low=1.0, factor_high=1.0;
   double	lambda, new_lambda; 		// Karlin-Altschul parameter

   // Bracket the values.
   lambda=lambda_ideal; done=FALSE; first_time=TRUE; factor=1.0;
   while (done != TRUE) {
   	for(c = 0; c < query_length; c++) {
       	    for(a = 0; a < alphabet_size; a++) {
		if(private_matrix[c][a]==GBLAST_SCORE_MIN) matrix[c][a]=GBLAST_SCORE_MIN;
		else matrix[c][a]=(factor*private_matrix[c][a])/POSIT_SCALE_FACTOR;
	    }
        }
        updateLambdaK( );
	new_lambda = kbp_psi[0]->Lambda;
	if(new_lambda > lambda) {
		if (first_time) {
			factor_high = 1.0 + POSIT_PERCENT;
			factor = factor_high; factor_low = 1.0;
			too_high = TRUE; first_time = FALSE;
		} else {
			if (too_high == FALSE) break;
			factor_high += (factor_high-1.0);
			factor = factor_high;
		}
	} else {
		if (first_time) {
			factor_high = 1.0;
			factor_low = 1.0 - POSIT_PERCENT;
			factor = factor_low;
			too_high = FALSE; first_time = FALSE;
		} else {
			if (too_high == TRUE) break;
			factor_low += (factor_low-1.0);
			factor = factor_low;
		}
	}
   } // binary search for ten times. 
   for (index=0; index<POSIT_NUM_ITERATIONS; index++) {
        factor = 0.5*(factor_high+factor_low);
   	for(c = 0; c < query_length; c++) {
       	    for(a = 0; a < alphabet_size; a++) {
		if(private_matrix[c][a]==GBLAST_SCORE_MIN) matrix[c][a]=GBLAST_SCORE_MIN;
		else matrix[c][a] = (factor*private_matrix[c][a])/POSIT_SCALE_FACTOR;
	    }
   	}
        updateLambdaK( );
	new_lambda = kbp_psi[0]->Lambda;
	if(new_lambda > lambda) factor_low=factor; else factor_high=factor; 
    }
   for(a = 0; a < alphabet_size; a++) matrix[query_length][a] = GBLAST_SCORE_MIN;
   return factor;
}

