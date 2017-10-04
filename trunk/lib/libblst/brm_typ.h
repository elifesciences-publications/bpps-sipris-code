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

// Adapted from posit.h by Alejandro Schaffer
#ifndef __BLST_RESCALE_MTX__
#define __BLST_RESCALE_MTX__

#include "my_blastkar.h"

/*************************** FROM BLASTOOL.C *******************************/
#define POSIT_PERCENT 0.05
#define POSIT_NUM_ITERATIONS 10

#define POSIT_SCALE_FACTOR 1000
#define Xchar   0    // AFN character for low-complexity columns
// #define Xchar   21    /*character for low-complexity columns*/

class brm_typ {	// _blast_matrix_rescale: Structure used for matrix rescaling.
public:
		brm_typ(Int4,Int4,unsigned char *,double *,Int4 **,Int4 **,
			kbp_typ *,kbp_typ *,kbp_typ *,kbp_typ *,double,double);
	double	Scale( );
	Int4	**Matrix( ){ return matrix; }
	// ~brm_typ( );  // not needed for now. == BlastMatrixRescaleDestruct( )
private:
	void		updateLambdaK( );
	sfp_typ		*ComputeProbs( );
	Int4		alphabet_size,
			query_length;	// length of query.
	unsigned char	*query;
	double		*standardProb;
	Int4		**matrix;
	Int4		**private_matrix;
	kbp_typ 	*kbp_std, 
			*kbp_psi, 
			*kbp_gap_std, 
			*kbp_gap_psi;
	double		lambda_ideal,
                	K_ideal;
}; 		// BlastMatrixRescale, *BlastMatrixRescalePtr;

#endif /* __BLST_RESCALE_MTX__ */

