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

/*****************************************************************************
File name: blastkar.h
Author: Tom Madden
Contents: definitions and prototypes used by blastkar.c to calculate GBLAST
	statistics.
******************************************************************************/
/* $Revision: 6.12 $ * $Log: blastkar.h,v $ * * */

#ifndef __MYGBLASTKAR__
#define __MYGBLASTKAR__

#include "stdinc.h"
#include "my_ncbi.h"
#include "my_ncbimath.h"
#include "afnio.h"
#include "sfp_typ.h"
#include "alphabet.h"

// Defines for the matrix 'preferences' (as specified by S. Altschul). 
#define GBLAST_MATRIX_NOMINAL 0
#define GBLAST_MATRIX_PREFERRED 1
#define GBLAST_MATRIX_BEST 2

/*************************************************************************
   Structure to the Karlin-Blk parameters.
   This structure was (more or less) copied from the old karlin.h.
**************************************************************************/

typedef struct {
	double	Lambda;		// Lambda value used in statistics.
	double	K, logK;	// K value used in statistics.
	double	H; 		// H value used in statistics.
	// "real" values are ones actually found, may be replaced by above.
	double	Lambda_real, K_real, logK_real, H_real;
	double	paramC;	/* for use in seed. */
} GBLAST_KarlinBlk, *kbp_typ;

/********************************************************************
	Structures relating to scoring or the GBLAST_ScoreBlk
SCORE_MIN is (-2**31 + 1)/2 because it has been observed more than once that
a compiler did not properly calculate the quantity (-2**31)/2.  The list
of compilers which failed this operation have at least at some time included:
NeXT and a version of AIX/370's MetaWare High C R2.1r.
For this reason, SCORE_MIN is not simply defined to be INT4_MIN/2.
********************************************************************/

#define GBLAST_MATRIX_SIZE 32

typedef struct _gblast_scoreblk {
	a_type	AB;			// AFN alphabet type
        double	searchsp_eff;		// AFN new....
	Int2 	alphabet_size;  	// size of alphabet.
	CharPtr name;			// name of matrix. 
	Int4	*matrix[GBLAST_MATRIX_SIZE];
	Int4	Int4_matrix[GBLAST_MATRIX_SIZE*GBLAST_MATRIX_SIZE];
	Int4	*maxscore; 		// Max. score for each letter 
	Int4	loscore, hiscore; 	// Min. & max. substitution scores
	sfp_typ	**sfp;			// score frequencies. 
		// kbp & kbp_gap are ptrs that should be set to kbp_std, kbp_psi, etc. 
	kbp_typ *kbp; 			// Karlin-Altschul parameters. 
	kbp_typ *kbp_gap; 		// K-A parameters for gapped alignments. 
		// Below are the Karlin-Altschul parameters for non-position based ('std')
		// and position based ('psi') searches. 
	kbp_typ *kbp_std,*kbp_psi,*kbp_gap_std,*kbp_gap_psi;
	kbp_typ kbp_ideal;  		// Ideal values 
					//(for query with average database composition). 
} GBLAST_ScoreBlk, *sbp_typ;

//===================  AFN defined routines: =====================

Int4	WordXDropoffSBP(double word_extend_dropoff_in_bits, sbp_typ sbp);
Int4	GapTriggerSBP(double gap_trigger_bits,sbp_typ sbp);
Int4    GapXDropoffSBP(double x_parameter_in_bits, sbp_typ sbp);
double	GappedBitScoreSBP(Int4 score, sbp_typ sbp);
void	PutKarlinAltschulSBP(FILE *fp,sbp_typ sbp);
double  GappedScoreToEvalueSBP(Int4 S, sbp_typ sbp);
Int4	GappedEvalueToScoreSBP(double Eval, sbp_typ sbp);
void	SetPsiStatsSBP(sbp_typ sbp);
void	SetStdStatsSBP(sbp_typ sbp);

#define	SMatrixSBP(sbp)	((sbp)->matrix)

/************************* DECLARATIONS *****************************/

sbp_typ GBLAST_ScoreBlkDestruct(sbp_typ sbp);
Int2	GBlastScoreBlkFill(sbp_typ sbp, unsigned char *string, Int4 length);
Int2	GBlastScoreBlkMaxScoreSet(sbp_typ sbp);
Int4	*GBlastResCompNew(sbp_typ sbp);
Int4	*GBlastResCompDestruct(Int4 *rcp);
Int2	GBlastResCompStr(sbp_typ sbp, Int4 *rcp, unsigned char *str, Int4 length);

// Produces a Karlin Block, and parameters, with standard protein frequencies.
Int2	GBlastKarlinBlkStandardCalc(sbp_typ sbp,Int2 context_start,Int2 context_end);
kbp_typ GBlastKarlinBlkStandardCalcEx(sbp_typ sbp);

// Functions taken from the OLD karlin.c 
kbp_typ GBlastKarlinBlkCreate(void);
kbp_typ GBlastKarlinBlkDestruct(kbp_typ);
Int2	GBlastKarlinBlkCalc(kbp_typ kbp, sfp_typ *sfp);
Int2	GBlastKarlinBlkGappedCalc(kbp_typ kbp,Int4 gap_open,Int4 gap_extend);
Int4	GBlastKarlinEtoS(double  E, kbp_typ kbp, double  qlen, double  dblen);
Int4	GBlastKarlinEtoS_simple(double  E, kbp_typ kbp, double searchsp); 
double	GBlastKarlinPtoE(double p);
double	GBlastKarlinEtoP(double x);
double	GBlastKarlinStoP(Int4 S, kbp_typ kbp, double  qlen, double  dblen);
double	GBlastKarlinStoP_simple(Int4 S, kbp_typ kbp, double  searchsp);
double	GBlastKarlinStoE_simple(Int4 S, kbp_typ kbp, double  searchsp);
double	GBlastKarlinStoLen(kbp_typ kbp, Int4 S);
Int2	GBlastCutoffs(Int4 *S, double *E, kbp_typ kbp, double qlen, double dblen, 
		Nlm_Boolean dodecay);
Int2	GBlastCutoffs_simple(Int4 *S, double *E, kbp_typ kbp, double search_sp, 
		Nlm_Boolean dodecay);

// SumP function. Called by GBlastSmallGapSumE and GBlastLargeGapSumE.
double GBlastSumP(Int4 r, double s);

// Functions to calculate SumE(for large and small gaps). 
double	GBlastSmallGapSumE(kbp_typ kbp, Int4 gap, double gap_prob,
	  double gap_decay_rate, Int2 num, Int4 sum, double xsum, 
	  Int4 query_length, Int4 subject_length, Boolean min_length_one);
double	GBlastLargeGapSumE(kbp_typ kbp, double gap_prob, 
	  double gap_decay_rate, Int2 num, Int4 sum,  double xsum, 
	  Int4 query_length, Int4 subject_length, Boolean old_stats);
CharPtr GBlastRepresentativeResidues(Int2 length); // Used for random sequences. 

sbp_typ GBLAST_ScoreBlkNew2(a_type alphabet);
sbp_typ GBLAST_ScoreBlkNew(a_type AB, unsigned char *query, Int4 query_length,
        Int4 gap_open, Int4 gap_extend, Int4 dblen, Int4 dbseq_num);

Int2	GBlastResFreqNormalize(sbp_typ sbp,double *rfp,double norm);
double	*GBlastResFreqNew(sbp_typ sbp);
double	*GBlastResFreqDestruct(double *rfp);
Int2	GBlastResFreqString(sbp_typ sbp,double *rfp,unsigned char *string,Int4 length);
Int2	GBlastResFreqStdComp(sbp_typ sbp, double *rfp);
Int2	GBlastResFreqResComp(sbp_typ sbp, double *rfp, Int4 *rcp);
Int2	GBlastResFreqClr(sbp_typ sbp, double *rfp);

// Functions used to convert between Stephen's pseudo scores and E or p-values.
Int2	GConvertPtoPseudoS(double p, double n);
Int2	GConvertEtoPseudoS(double E, double searchsp);
double	GConvertPseudoStoE(Int2 s, double n);

/****************************************************************************
Obtains arrays of the allowed opening and extension penalties for gapped GBLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  The pref_flags field is
used for display purposes, with the defines: 
  GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_PREFERRED, and GBLAST_MATRIX_BEST.

Any of these fields that are not required should be set to NULL.  
The Int2 return value is the length of the arrays.
******************************************************************************/
Int2 GBlastKarlinGetMatrixValues(Int4Ptr *open, Int4Ptr *extension, 
	double **lambda, double **K, double **H, Int4Ptr *pref_flags);


#endif /* !__MYGBLASTKAR__ */

