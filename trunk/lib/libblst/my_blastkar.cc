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

#include "my_blastkar.h"

/*****************************************************************************
Taken from: blastkar.c by Tom Madden
Contents: Functions to calculate GBLAST probabilities etc.

	- allocate and deallocate structures used by GBLAST to calculate
	probabilities etc.
	- calculate residue frequencies for query and "average" database.
	- read in matrix.
        - calculate sum-p from a collection of HSP's, for both the case
	of a "small" gap and a "large" gap, when give a total score and the
	number of HSP's.
	- calculate expect values for p-values.
	- calculate pseuod-scores from p-values.
******************************************************************************/
/* $Revision: 6.29 $ = $Log: blastkar.c,v $ */

static Int2 GBlastScoreFreqCalc(sbp_typ sbp, sfp_typ *sfp, double  *rfp1, double  *rfp2);

/* OSF1 apparently doesn't like this. */
#if defined(HUGE_VAL) && !defined(OS_UNIX_OSF1)
#define GBLASTKAR_HUGE_VAL HUGE_VAL
#else
#define GBLASTKAR_HUGE_VAL 1.e30
#endif

static double GBlastSumPCalc(int r, double s); // Calculates sump for GBlastSumPStd 

typedef double array_of_5[5];		// Used in GBlastKarlinBlkGappedCalc.
// Used to temporarily store matrix values for retrieval. 

typedef struct _matrix_info {
	CharPtr		name;			// name of matrix (e.g., BLOSUM90). 
	array_of_5 	*values;		// The values (below).
	Int4		*prefs;			// Preferences for display. 
	Int4		max_number_values;	// number of values (e.g., BLOSUM90_VALUES_MAX). 
} MatrixInfo, *mip_typ;

/**************************************************************************************
		NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!
	The arrays below list all the matrices allowed for gapped GBLAST.
	If new ones are added, remember to edit the function GBlastLoadMatrixValues!!!
	(i.e., MatrixInfoNew( );)
***************************************************************************************/
	
#define BLOSUM62_20_VALUES_MAX 8
static double blosum62_20_values[BLOSUM62_20_VALUES_MAX][5] = {
	{(double) INT2_MAX, (double) INT2_MAX, 0.03391, 0.125, 0.45},
	{100.0, 10.0, 0.0293, 0.054, 0.26},
	{97.0, 10.0, 0.0286, 0.046, 0.24},
	{95.0, 10.0, 0.0280, 0.041, 0.22},
	{93.0, 10.0, 0.0273, 0.035, 0.21},
	{90.0, 10.0, 0.0266, 0.031, 0.19},
	{88.0, 10.0, 0.0261, 0.029, 0.17},
	{87.0, 10.0, 0.0256, 0.026, 0.16}
}; 

static Int4 blosum62_20_prefs[BLOSUM62_20_VALUES_MAX] = {
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL
};

#define BLOSUM62_VALUES_MAX 7

#if 1	// WARNING: change posComputePseudoFreqs if change this...
// gap-open, gap-extend, lambda, K, H.
static double blosum62_values[BLOSUM62_VALUES_MAX][5] = {
	{(double) INT2_MAX, (double) INT2_MAX,     0.3176,     0.134,    0.40},
	{  9.0,  2.0,     0.285,     0.075,    0.27},
	{  8.0,  2.0,     0.265,     0.046,    0.22},
	{  7.0,  2.0,     0.243,     0.032,    0.17},
	{ 12.0,  1.0,     0.281,     0.057,    0.27},
	{ 11.0,  1.0,     0.270,     0.047,    0.23},
	{ 10.0,  1.0,     0.250,     0.033,    0.17},
}; 
#else	// Use adjusted values...
static double blosum62_values[BLOSUM62_VALUES_MAX][5] = {
	{(double) INT2_MAX, (double) INT2_MAX,     0.3176,     0.134,    0.40},
	{  9.0,  2.0,     0.279,     0.058,    0.19},
	{  8.0,  2.0,     0.264,     0.045,    0.15},
	{  7.0,  2.0,     0.239,     0.027,    0.10},
	{ 12.0,  1.0,     0.283,     0.059,    0.19},
	{ 11.0,  1.0,     0.267,     0.041,    0.14},	// best...
	{ 10.0,  1.0,     0.243,     0.024,    0.10},
}; 
#endif
				
static Int4 blosum62_prefs[BLOSUM62_VALUES_MAX] = {
GBLAST_MATRIX_NOMINAL, 
GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_PREFERRED, 
GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_PREFERRED, 
GBLAST_MATRIX_BEST,
GBLAST_MATRIX_PREFERRED
};

#define BLOSUM50_VALUES_MAX 13
static double blosum50_values[BLOSUM50_VALUES_MAX][5] = {
	{(double) INT2_MAX, (double) INT2_MAX,     0.232,     0.11,      0.34},
	{12.0,  3.0,     0.206,     0.055,      0.23},
	{11.0,  3.0,     0.198,     0.046,      0.20},
	{10.0,  3.0,     0.189,     0.038,      0.17},
	{ 9.0,  3.0,     0.177,     0.030,      0.14},
	{15.0,  2.0,     0.211,     0.062,      0.25},
	{14.0,  2.0,     0.205,     0.053,      0.22},
	{13.0,  2.0,    0.197,     0.043,      0.19},
	{12.0,  2.0,     0.183,     0.028,      0.15},
	{18.0,  1.0,     0.208,     0.055,      0.23},
	{17.0,  1.0,     0.200,     0.042,      0.20},
	{16.0,  1.0,     0.189,     0.030,      0.17},
	{15.0,  1.0,     0.175,     0.020,      0.13},
};

static Int4 blosum50_prefs[BLOSUM50_VALUES_MAX] = {
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_NOMINAL,
GBLAST_MATRIX_NOMINAL };

#define BLOSUM45_VALUES_MAX 13
static double blosum45_values[BLOSUM45_VALUES_MAX][5] = {
	{(double) INT2_MAX, (double) INT2_MAX,     0.2291,     0.092,      0.25},
	{13.0,  3.0,     0.209,     0.057,      0.19},
	{12.0,  3.0,     0.203,     0.049,      0.17},
	{11.0,  3.0,     0.193,     0.037,      0.15},
	{10.0,  3.0,     0.182,     0.029,      0.12},
	{15.0,  2.0,     0.206,     0.049,      0.18},
	{14.0,  2.0, 	 0.199,     0.040,      0.16},
	{13.0,  2.0,	 0.190,     0.032,      0.14},
	{12.0,  2.0,     0.177,     0.023,      0.11},
	{19.0,  1.0,     0.209,     0.049,      0.19},
	{18.0,  1.0,     0.202,     0.041,      0.17},
	{17.0,  1.0,     0.195,     0.034,      0.14},
	{16.0,  1.0,     0.183,     0.024,      0.12}
};

static Int4 blosum45_prefs[BLOSUM45_VALUES_MAX] = {
GBLAST_MATRIX_NOMINAL, GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_BEST, GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED, GBLAST_MATRIX_PREFERRED,
GBLAST_MATRIX_PREFERRED };

sbp_typ GBLAST_ScoreBlkNew2(a_type AB)
// Allocates memory for the sbp_typ. 
{
	sbp_typ sbp;

	assert(nAlpha(AB) < GBLAST_MATRIX_SIZE);
	NEW(sbp,1,GBLAST_ScoreBlk);
	sbp->AB=AB;	// my alphabet structure.
	if(sbp != NULL) {
		sbp->alphabet_size = NAlpha(AB); // includes gap ('-') character...
		for(Int4 index=0; index < GBLAST_MATRIX_SIZE-1; index++) {
		  sbp->matrix[index] = sbp->Int4_matrix + index*GBLAST_MATRIX_SIZE;
		} NEW(sbp->maxscore,GBLAST_MATRIX_SIZE,Int4);
		sbp->sfp = new sfp_typ* [1];
		sbp->kbp_std = (kbp_typ*) GMemNew(sizeof(kbp_typ));
		sbp->kbp_gap_std = (kbp_typ *) GMemNew(sizeof(kbp_typ));
		sbp->kbp_psi = (kbp_typ*) GMemNew(sizeof(kbp_typ));
		sbp->kbp_gap_psi = (kbp_typ *) GMemNew(sizeof(kbp_typ));
	} // Fill in the GBLAST_ScoreBlk structure.  
	Int4	**matrix = sbp->matrix;	
	for(Int4 index1 =0; index1 <= nAlpha(AB); index1++){
	   for(Int4 index2 = 0; index2 <= nAlpha(AB); index2++){
		matrix[index1][index2] = (Int4) valAlphaR(index1,index2,AB);
	   } // Initialize an sbp_typ object.
	} GBlastScoreBlkMaxScoreSet(sbp);
	return sbp;
}

//================= AFN DEFINED ROUTINES ====================
static const double sbp_typ_ln2 = 0.693147180559945309;

Int4    WordXDropoffSBP(double word_x_dropoff_in_bits, sbp_typ sbp)
{ return (Int4) ceil(word_x_dropoff_in_bits*sbp_typ_ln2/sbp->kbp[0]->Lambda); }

Int4    GapTriggerSBP(double gap_trigger_bits,sbp_typ sbp)
{ return (Int4) floor((gap_trigger_bits*sbp_typ_ln2+sbp->kbp[0]->logK)
                			/ sbp->kbp[0]->Lambda); }

Int4	GapXDropoffSBP(double x_parameter_in_bits, sbp_typ sbp)
{ return (Int4) ((x_parameter_in_bits*sbp_typ_ln2)/sbp->kbp_gap[0]->Lambda);
} // use ceil??

double	GappedBitScoreSBP(Int4 score, sbp_typ sbp)
{ return ((score*sbp->kbp_gap[0]->Lambda) - sbp->kbp_gap[0]->logK)/sbp_typ_ln2; }

void    PutKarlinAltschulSBP(FILE *fp,sbp_typ sbp)
{
	fprintf(fp,"\neffective search space: %-15.0f\n", sbp->searchsp_eff);
	fprintf(fp,"ungapped: Lambda =%.3f; K = %.3f; H = %.3f\n",
           sbp->kbp[0]->Lambda, sbp->kbp[0]->K, sbp->kbp[0]->H);
  	fprintf(fp,"gapped:   Lambda =%.3f; K = %.3f; H = %.3f\n",
           sbp->kbp_gap[0]->Lambda, sbp->kbp_gap[0]->K, sbp->kbp_gap[0]->H);
}

double  GappedScoreToEvalueSBP(Int4 S, sbp_typ sbp)
{ return GBlastKarlinStoE_simple(S, sbp->kbp_gap[0], sbp->searchsp_eff); }

Int4	GappedEvalueToScoreSBP(double Eval, sbp_typ sbp)
{ return GBlastKarlinEtoS_simple(Eval,sbp->kbp_gap[0],sbp->searchsp_eff); }

void    SetPsiStatsSBP(sbp_typ sbp)
{ sbp->kbp = sbp->kbp_psi; sbp->kbp_gap = sbp->kbp_gap_psi; }

void    SetStdStatsSBP(sbp_typ sbp)
{ sbp->kbp = sbp->kbp_std; sbp->kbp_gap = sbp->kbp_gap_std; }

//===========================================================

Int2 GBlastScoreBlkMaxScoreSet(sbp_typ sbp)
{
	Int4 	score, maxscore;

	sbp->loscore = GBLAST_SCORE_1MAX;
        sbp->hiscore = GBLAST_SCORE_1MIN;
	Int4	**matrix = sbp->matrix;
	for(short index1=0; index1<=nAlpha(sbp->AB); index1++) {
	  maxscore=GBLAST_SCORE_MIN;
	  for (short index2=0; index2<=nAlpha(sbp->AB); index2++) {
		score = matrix[index1][index2];
		if (score <= GBLAST_SCORE_MIN || score >= GBLAST_SCORE_MAX)
				continue;
		if (score > maxscore) { maxscore = score; }
		if (sbp->loscore > score) sbp->loscore = score;
		if (sbp->hiscore < score) sbp->hiscore = score;
	  } sbp->maxscore[index1] = maxscore;
	}
	// If the lo/hi-scores are GBLAST_SCORE_MIN/GBLAST_SCORE_MAX, 
	// (i.e., for gaps), then use other scores. */
	if (sbp->loscore < GBLAST_SCORE_1MIN) sbp->loscore = GBLAST_SCORE_1MIN;
	if (sbp->hiscore > GBLAST_SCORE_1MAX) sbp->hiscore = GBLAST_SCORE_1MAX;
	return 0;
}

sbp_typ GBLAST_ScoreBlkDestruct(sbp_typ sbp)
{
	if (sbp == NULL) return NULL;
	if(sbp->sfp) delete sbp->sfp[0];
	if(sbp->kbp_std) sbp->kbp_std[0] = 
		GBlastKarlinBlkDestruct(sbp->kbp_std[0]);
	if(sbp->kbp_gap_std) sbp->kbp_gap_std[0]=
		GBlastKarlinBlkDestruct(sbp->kbp_gap_std[0]);
	if(sbp->kbp_psi) sbp->kbp_psi[0] = 
		GBlastKarlinBlkDestruct(sbp->kbp_psi[0]);
	if(sbp->kbp_gap_psi) sbp->kbp_gap_psi[0] = 
		GBlastKarlinBlkDestruct(sbp->kbp_gap_psi[0]);
	if(sbp->kbp_ideal) sbp->kbp_ideal = GBlastKarlinBlkDestruct(sbp->kbp_ideal);
	delete [] sbp->sfp;
	sbp->kbp_std = (kbp_typ *) GMemFree(sbp->kbp_std);
	sbp->kbp_psi = (kbp_typ *) GMemFree(sbp->kbp_psi);
	sbp->kbp_gap_std = (kbp_typ *) GMemFree(sbp->kbp_gap_std);
	sbp->kbp_gap_psi = (kbp_typ *) GMemFree(sbp->kbp_gap_psi);
	sbp->maxscore = (Int4 *) GMemFree(sbp->maxscore);
	sbp->name = (Nlm_CharPtr) GMemFree(sbp->name);
	sbp = (sbp_typ) GMemFree(sbp);
	return sbp;
}

Int2	GBlastScoreBlkFill(sbp_typ sbp, unsigned char *query, Int4 query_length)
// Calculate the Karlin parameters.  This function should be called once
// for each context, or frame translated.
// The rfp and stdrfp are calculated for each context, this should be fixed. 
{
	double  *rfp, *stdrfp;
	Int2 retval=0;
	sfp_typ	*sfp;

	rfp = GBlastResFreqNew(sbp);
	stdrfp = GBlastResFreqNew(sbp);
	GBlastResFreqStdComp(sbp, stdrfp);
	GBlastResFreqString(sbp, rfp, query, query_length);
	sfp = new sfp_typ(sbp->loscore, sbp->hiscore);
	sbp->sfp[0] = sfp;

	GBlastScoreFreqCalc(sbp, sfp, rfp, stdrfp);
	sbp->kbp_std[0] = GBlastKarlinBlkCreate();
	retval = GBlastKarlinBlkCalc(sbp->kbp_std[0],sbp->sfp[0]);
	if (retval) {
		rfp = GBlastResFreqDestruct(rfp);
		stdrfp = GBlastResFreqDestruct(stdrfp);
		return retval;
	}
	sbp->kbp_psi[0] = GBlastKarlinBlkCreate();
	retval = GBlastKarlinBlkCalc(sbp->kbp_psi[0], sbp->sfp[0]);
	rfp = GBlastResFreqDestruct(rfp);
	stdrfp = GBlastResFreqDestruct(stdrfp);
	return retval;
}

kbp_typ GBlastKarlinBlkStandardCalcEx(sbp_typ sbp)
//	Calculates the standard Karlin parameters.  This is used
//	if the query is translated and the calculated (real) Karlin
//	parameters are bad, as they're calculated for non-coding regions.
{
	kbp_typ kbp_ideal;
	double  *stdrfp;

	stdrfp = GBlastResFreqNew(sbp);
	GBlastResFreqStdComp(sbp, stdrfp);
	sfp_typ *sfp = new sfp_typ(sbp->loscore, sbp->hiscore);
	GBlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
	kbp_ideal = GBlastKarlinBlkCreate();
	GBlastKarlinBlkCalc(kbp_ideal,sfp);
	stdrfp = GBlastResFreqDestruct(stdrfp);
	delete sfp;
	return kbp_ideal;
}

Int2 GBlastKarlinBlkStandardCalc(sbp_typ sbp, Int2 context_start, 
	Int2 context_end)
{
	kbp_typ kbp_ideal, kbp;
	Int2 index;

	kbp_ideal = GBlastKarlinBlkStandardCalcEx(sbp);
/* Replace the calculated values with ideal ones for blastx, tblastx. */
	for (index=context_start; index<=context_end; index++){
		kbp = sbp->kbp[index];	
		if (kbp->Lambda >= kbp_ideal->Lambda) {
			kbp->Lambda = kbp_ideal->Lambda;
			kbp->K = kbp_ideal->K;
			kbp->logK = kbp_ideal->logK;
			kbp->H = kbp_ideal->H;
		}
	} kbp_ideal = GBlastKarlinBlkDestruct(kbp_ideal);
	return 0;
}

kbp_typ GBlastKarlinBlkCreate(void) // Creates the Karlin Block. 
{ kbp_typ kbp = (kbp_typ) GMemNew(sizeof(GBLAST_KarlinBlk)); return kbp; }

kbp_typ GBlastKarlinBlkDestruct(kbp_typ kbp) // Deallocates the Karlin Block. 
{ kbp = (kbp_typ) GMemFree(kbp); return kbp; }

Int4	*GBlastResCompNew(sbp_typ sbp)
// 	Allocated the rcp for a given alphabet.  Only the
// 	alphabets ncbistdaa and ncbi4na should be used by GBLAST.
{
	Int4 *rcp;
	rcp = (Int4 *) GMemNew(GBLAST_MATRIX_SIZE*sizeof(Int4));
	assert(rcp != NULL); return rcp;
}

Int4	*GBlastResCompDestruct(Int4 *rcp)
{
	if (rcp == NULL) return NULL;
	if (rcp != NULL) rcp = (Int4 *) GMemFree(rcp);
	return rcp;
}

Int2	GBlastResCompStr(sbp_typ sbp, Int4 *rcp, unsigned char *str, Int4 length)
// Store the composition of a (query) string.  
{
	unsigned char	*lp,*lpmax;
	Int2 index;

	if(sbp == NULL || rcp == NULL || str == NULL) return 1;
	for(index=0; index<(sbp->alphabet_size); index++) rcp[index] = 0;
	for(lp = str, lpmax = lp+length; lp < lpmax; lp++) { ++rcp[*lp]; }
	// Don't count ambig. residues.
	rcp[0]=0;	// assumes that AlphaCode( ) == 0 is ambig. residue!!!
	// AFN: NEED TO FIX THIS! Want something like rcp[UndefAlpha(AB)]=0;
	return 0;
}

static Int2 GBlastScoreChk(Int4 lo, Int4 hi)
{
   if (lo>=0 || hi<=0 || lo<GBLAST_SCORE_1MIN || hi>GBLAST_SCORE_1MAX) return 1;
   if (hi - lo > GBLAST_SCORE_RANGE_MAX) return 1;
   return 0;
}

static Int2 GBlastScoreFreqCalc(sbp_typ sbp, sfp_typ *sfp, double *rfp1, double *rfp2)
{
	if (sbp == NULL || sfp == NULL) return 1;
	if(sfp->CheckScore(sbp->loscore,sbp->hiscore)) return 1;
	sfp->ScoreFreqCalc(sbp->matrix, sbp->alphabet_size,sbp->loscore, rfp1,rfp2);
	return 0;
}

typedef struct {
	Char		ch;
	double	p;
} GBLAST_LetterProb;

#define STD_AMINO_ACID_FREQS Robinson_prob

#if STD_AMINO_ACID_FREQS == Dayhoff_prob
/*  M. O. Dayhoff amino acid background frequencies   */
static GBLAST_LetterProb	Dayhoff_prob[] = {
	{ 'A', 87.13 }, { 'C', 33.47 }, { 'D', 46.87 }, { 'E', 49.53 },
	{ 'F', 39.77 }, { 'G', 88.61 }, { 'H', 33.62 }, { 'I', 36.89 },
	{ 'K', 80.48 }, { 'L', 85.36 }, { 'M', 14.75 }, { 'N', 40.43 },
	{ 'P', 50.68 }, { 'Q', 38.26 }, { 'R', 40.90 }, { 'S', 69.58 },
	{ 'T', 58.54 }, { 'V', 64.72 }, { 'W', 10.49 }, { 'Y', 29.92 } };
#endif

#if STD_AMINO_ACID_FREQS == Altschul_prob
/* Stephen Altschul amino acid background frequencies */
static GBLAST_LetterProb Altschul_prob[] = {
	{ 'A', 81.00 }, { 'C', 15.00 }, { 'D', 54.00 }, { 'E', 61.00 },
	{ 'F', 40.00 }, { 'G', 68.00 }, { 'H', 22.00 }, { 'I', 57.00 },
	{ 'K', 56.00 }, { 'L', 93.00 }, { 'M', 25.00 }, { 'N', 45.00 },
	{ 'P', 49.00 }, { 'Q', 39.00 }, { 'R', 57.00 }, { 'S', 68.00 },
	{ 'T', 58.00 }, { 'V', 67.00 }, { 'W', 13.00 }, { 'Y', 32.00 } };
#endif

#if STD_AMINO_ACID_FREQS == Robinson_prob
/* amino acid background frequencies from Robinson and Robinson */
static GBLAST_LetterProb Robinson_prob[] = {
	{ 'A', 78.00 }, { 'C', 19.00 }, { 'D', 54.00 }, { 'E', 63.00 },
	{ 'F', 39.00 }, { 'G', 74.00 }, { 'H', 22.00 }, { 'I', 52.00 },
	{ 'K', 57.00 }, { 'L', 90.00 }, { 'M', 22.00 }, { 'N', 45.00 },
	{ 'P', 52.00 }, { 'Q', 43.00 }, { 'R', 51.00 }, { 'S', 71.00 },
	{ 'T', 59.00 }, { 'V', 64.00 }, { 'W', 13.00 }, { 'Y', 32.00 } };
#endif

double	*GBlastResFreqNew(sbp_typ sbp)
// Allocates the rfp_typ and fill in the frequencies in the probabilities.
{
	if (sbp == NULL) { return NULL; }
	double	*rfp = (double *) GMemNew(sizeof(double) * sbp->alphabet_size);
	assert(rfp != NULL);
	return rfp;
}

/* Normalize the frequencies to "norm".  */
Int2 GBlastResFreqNormalize(sbp_typ sbp, double *rfp, double norm)
{
	Int2	alphabet_stop, index;
	double	sum = 0., p;

	if (norm == 0.) return 1;
	alphabet_stop = sbp->alphabet_size;
	for (index=0; index<alphabet_stop; index++) {
		p = rfp[index];
		if (p < 0.) return 1;
		sum += p;
	}
	if (sum <= 0.) return 0;
	for (index=0; index<alphabet_stop; index++) {
		rfp[index] /= sum; rfp[index] *= norm;
// std::cerr << rfp[index]; std::cerr << std::endl;
	}
	return 0;
}

Int2	GBlastResFreqStdComp(sbp_typ sbp, double *rfp)
{
	Int2	index, retval;
	Int4	cd;

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) {
		cd =AlphaCode(STD_AMINO_ACID_FREQS[index].ch,sbp->AB);
		rfp[cd] = STD_AMINO_ACID_FREQS[index].p;
// std::cerr << rfp[cd]; std::cerr << std::endl;
	}
	GBlastResFreqNormalize(sbp, rfp, 1.0);
	return 0;
}

CharPtr  GBlastRepresentativeResidues(Int2 length)
{
	CharPtr buffer, ptr;
	Int2 index, total;
	Int4 number;
	total=0;

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) {
		total += (Int2) STD_AMINO_ACID_FREQS[index].p;
	}
	buffer = (CharPtr) GMemNew((length+1)*sizeof(Char));
	ptr = buffer;
	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) {
		number = GNlm_Nint((STD_AMINO_ACID_FREQS[index].p)
			*((double) length)/((double) total));
		while (number > 0) {
			*ptr = STD_AMINO_ACID_FREQS[index].ch;
			ptr++; number--;
		}
	} return buffer;
}

Int2	GBlastResFreqString(sbp_typ sbp, double *rfp, unsigned char *string, Int4 length)
{
	Int4 *rcp;
	rcp = GBlastResCompNew(sbp);
	GBlastResCompStr(sbp, rcp, string, length);
	GBlastResFreqResComp(sbp, rfp, rcp);
	rcp = GBlastResCompDestruct(rcp);
	return 0;
}

double  *GBlastResFreqDestruct(double  *rfp)
{
	if (rfp == NULL) return NULL;
	rfp = (double *) GMemFree(rfp);
	return rfp;
}

Int2	GBlastResFreqResComp(sbp_typ sbp, double *rfp, Int4 *rcp)
// Calculate the residue frequencies associated with the provided ResComp 
{
	Int2	alphabet_max, index;
	double	sum = 0.;

	if (rfp == NULL || rcp == NULL) return 1;
	alphabet_max = sbp->alphabet_size;
	for (index=0; index<alphabet_max; index++)
		sum += rcp[index];
	if (sum == 0.) { GBlastResFreqClr(sbp, rfp); return 0; }
	for (index=0; index<alphabet_max; index++)
		rfp[index] = rcp[index] / sum;
	return 0;
}

Int2	GBlastResFreqClr(sbp_typ sbp, double *rfp)
{
	Int2	alphabet_max, index;
 
	alphabet_max = sbp->alphabet_size;
	for (index=0; index<alphabet_max; index++)
                rfp[index] = 0.0;
        return 0;
}


/* Deallocates mip_typ */
static mip_typ MatrixInfoDestruct(mip_typ matrix_info)
{
	if (matrix_info == NULL) return NULL;
	GMemFree(matrix_info->name);
	return (mip_typ) GMemFree(matrix_info);
}

/* Makes New mip_typ */
static mip_typ MatrixInfoNew(CharPtr name, array_of_5 *values, 
	Int4Ptr prefs, Int4 max_number)
{
	mip_typ matrix_info;
	matrix_info = (mip_typ) GMemNew(sizeof(MatrixInfo));
	matrix_info->name = AllocString(name);
	matrix_info->values = values;
	matrix_info->prefs = prefs;
	matrix_info->max_number_values = max_number;
	return matrix_info;
}

Int2 GBlastKarlinGetMatrixValues(Int4Ptr *open, Int4Ptr *extension, 
	double **lambda, double **K, double **H, Int4Ptr *pref_flags)
/***************************************************************************
Obtains arrays of the allowed opening and extension penalties for gapped GBLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
****************************************************************************/
{
	array_of_5	*values;
	Int4		index, max_number_values=0;
	Int4Ptr		open_array=NULL, extension_array=NULL, pref_flags_array=NULL, prefs;
	double		*lambda_array=NULL,*K_array=NULL,*H_array=NULL;
	mip_typ	matrix_info;

	matrix_info=MatrixInfoNew("BLOSUM62",blosum62_values,blosum62_prefs,
				BLOSUM62_VALUES_MAX);
	values = matrix_info->values;
	max_number_values = matrix_info->max_number_values;
	prefs = matrix_info->prefs;

	if(open) *open = open_array = 
		(Nlm_Int4Ptr) GMemNew(max_number_values*sizeof(Int4));
	if(extension) *extension = extension_array = 
		(Nlm_Int4Ptr) GMemNew(max_number_values*sizeof(Int4));
	if(lambda) *lambda = lambda_array = (double *) GMemNew(max_number_values*sizeof(double ));
	if(K) *K = K_array = (double *) GMemNew(max_number_values*sizeof(double ));
	if(H) *H = H_array = (double *) GMemNew(max_number_values*sizeof(double ));
	if(pref_flags) *pref_flags = pref_flags_array = 
		(Nlm_Int4Ptr) GMemNew(max_number_values*sizeof(Int4));
	for (index=0; index<max_number_values; index++) {
		if (open) open_array[index] = values[index][0];
		if (extension) extension_array[index] = values[index][1];
		if (lambda) lambda_array[index] = values[index][2];
		if (K) K_array[index] = values[index][3];
		if (H) H_array[index] = values[index][4];
		if (pref_flags) pref_flags_array[index] = prefs[index];
	}
	MatrixInfoDestruct(matrix_info);
	return max_number_values;
}
	
Int2 GBlastKarlinBlkGappedCalc(kbp_typ kbp,Int4 gap_open,Int4 gap_extend)
/**************************************************************************
	Supplies lambda, H, and K values, as calcualted by Stephen Altschul 
	in Methods in Enzy. (vol 266, page 474).
	if kbp is NULL, then a validation is perfomed.
	// ValNodePtr *error_return as input parameter eliminated.
***************************************************************************/
{
	Boolean		found_values;
	array_of_5	*values;
	Int4		index, max_number_values=0;
	mip_typ		matrix_info;

	values = (array_of_5 *) NULL;
	found_values = FALSE;
	// AFN: ELIMINATED COMPLEXITY BY ALLOWING ONLY BLOSUM62 MATRIX for NOW.
	matrix_info=MatrixInfoNew("BLOSUM62",blosum62_values,blosum62_prefs,
					BLOSUM62_VALUES_MAX);
	values = matrix_info->values;
	max_number_values = matrix_info->max_number_values;
	for (index=0; index < max_number_values; index++) {
		if (GNlm_Nint(values[index][0]) == gap_open &&
			GNlm_Nint(values[index][1]) == gap_extend) {
			if(kbp) {
				kbp->Lambda_real = kbp->Lambda = values[index][2];
				kbp->K_real = kbp->K = values[index][3];
				kbp->logK_real = kbp->logK = log(kbp->K);
				kbp->H_real = kbp->H = values[index][4];
			}
			found_values = TRUE; break;
		}
	}
	MatrixInfoDestruct(matrix_info);
	if(found_values==TRUE) return 0;
	fprintf(stderr,
	  "Input gap opening and extension values of %d and %d not supported",
		(Int4) gap_open, (Int4) gap_extend);
	print_error("GBlastKarlinBlkGappedCalc error");
	return 1;
}

//      Everything below here was (more or less) copied from the old 
//      karlin.c and could work separately from the stuff above. 
Int2	GBlastKarlinBlkCalc(kbp_typ kbp, sfp_typ *sfp)
{
	if (kbp == NULL || sfp == NULL) return 1;
	/* Calculate the parameter Lambda */
	kbp->Lambda_real = kbp->Lambda = sfp->KarlinLambdaNR( );
	if (kbp->Lambda < 0.) goto ErrExit;
	/* Calculate H */
	kbp->H_real = kbp->H = sfp->KarlinLtoH(kbp->Lambda);
	if (kbp->H < 0.) goto ErrExit;
	/* Calculate K and log(K) */
	kbp->K_real = kbp->K = sfp->KarlinLHtoK(kbp->Lambda, kbp->H);
	if (kbp->K < 0.) goto ErrExit;
	kbp->logK_real = kbp->logK = log(kbp->K);
	/* Normal return */
	return 0;
ErrExit:
	kbp->Lambda = kbp->H = kbp->K
		= kbp->Lambda_real = kbp->H_real = kbp->K_real = -1.;
#ifdef GBLASTKAR_HUGE_VAL
	kbp->logK_real = kbp->logK = GBLASTKAR_HUGE_VAL;
#else
	kbp->logK_real = kbp->logK = 1.e30;
#endif
	return 1;
}

static double GBlastGapDecayInverse(double pvalue, unsigned nsegs, 
	double decayrate)
{
	if (decayrate <= 0. || decayrate >= 1. || nsegs == 0) return pvalue;
	return pvalue * (1. - decayrate) * GNlm_Powi(decayrate, nsegs - 1);
}

static double GBlastGapDecay(double pvalue, unsigned nsegs, 
	double decayrate)
{
	if (decayrate <= 0. || decayrate >= 1. || nsegs == 0) return pvalue;
	return pvalue / ((1. - decayrate) * GNlm_Powi(decayrate, nsegs - 1));
}

//	GBlastCutoffs
// 	Calculate the cutoff score, S, and the highest expected score.
//	WRG (later modified by TLM).
// 	S = cutoff score; E = expected no. of HSPs scoring at or above S
//	qlen = length of query sequence; 
//	dblen = length of database or database sequence
//	dodecay == TRUE ==> use gapdecay feature.
Int2 GBlastCutoffs(Int4 *S, double *E, kbp_typ kbp, double qlen,
	double dblen, Nlm_Boolean dodecay) 
{ return GBlastCutoffs_simple(S, E, kbp, qlen*dblen, dodecay); }

Int2 GBlastCutoffs_simple(Int4 *S, double *E, kbp_typ kbp, double searchsp,
	Nlm_Boolean dodecay)
{
	Int4	s = *S, es;
	double	e = *E, esave;
	Boolean	s_changed = FALSE;

	if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.) return 1;
	/* Calculate a cutoff score, S, from the Expected
	(or desired) number of reported HSPs, E.  */
	es = 1; esave = e;
	if (e > 0.) {
		if (dodecay) e = GBlastGapDecayInverse(e, 1, 0.5);
		es = GBlastKarlinEtoS_simple(e, kbp, searchsp);
	}
	/* Pick the larger cutoff score between the user's choice
	and that calculated from the value of E.  */
	if (es > s) { s_changed = TRUE; *S = s = es; }

	/* Re-calculate E from the cutoff score, if E going in was too high */
	if (esave <= 0. || !s_changed) {
		e = GBlastKarlinStoE_simple(s, kbp, searchsp);
		if (dodecay) e = GBlastGapDecay(e, 1, 0.5);
		*E = e;
	}
	return 0;
}

// GBlastKarlinEtoS() -- given an Expect value, return the associated cutoff score
// 	Error return value is GBLAST_SCORE_MIN
Int4 GBlastKarlinEtoS(double E, kbp_typ	kbp, double qlen, double dblen)
{ return GBlastKarlinEtoS_simple(E, kbp, qlen*dblen); }


Int4	GBlastKarlinEtoS_simple(double	E, kbp_typ kbp, double	searchsp)
{

	double	Lambda, K, H; /* parameters for Karlin statistics */
	Int4	S;
	Lambda = kbp->Lambda;
	K = kbp->K; H = kbp->H;
	if (Lambda < 0. || K < 0. || H < 0.) { return GBLAST_SCORE_MIN; }
	S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
	return S;
}

// 	GBlastKarlinStoP:  Calculate the probability (as opposed to expectation)
//	of achieving a particular score.  On error, return value is -1.
double GBlastKarlinStoP(Int4 S, kbp_typ kbp,
		double	qlen,	/* length of query sequence */
		double	dblen)	/* length of database */
{ return GBlastKarlinStoP_simple(S, kbp, qlen*dblen); }

double GBlastKarlinStoP_simple(Int4 S, kbp_typ kbp, double searchsp)
{
	double x = GBlastKarlinStoE_simple(S, kbp, searchsp);
	if (x == -1.) return x;
	return GBlastKarlinEtoP(x);
}

double GBlastKarlinStoE_simple(Int4 S, kbp_typ kbp, double searchsp)
{
	double Lambda = kbp->Lambda; double K = kbp->K; 
	double H = kbp->H;
	if(Lambda < 0. || K < 0. || H < 0.) { return -1.; }
	return searchsp * exp((double)(-Lambda * S) + kbp->logK);
}

double GBlastKarlinPtoE(double p)
// GBlastKarlinPtoE -- convert a P-value to an Expect value
{
        if (p < 0. || p > 1.0) { return INT4_MIN; }
	if (p == 1) return INT4_MAX;
        return -GNlm_Log1p(-p);
}

// GBlastKarlinEtoP -- convert an Expect value to a P-value.
double GBlastKarlinEtoP(double x) { return -GNlm_Expm1((double)-x); }

// Given a score, return the length expected for an HSP of that score
double GBlastKarlinStoLen(kbp_typ kbp, Int4 S){ return kbp->Lambda*S/kbp->H; }

static double	tab2[] = { /* table for r == 2 */
0.01669,  0.0249,   0.03683,  0.05390,  0.07794,  0.1111,   0.1559,   0.2146,   
0.2890,   0.3794,   0.4836,   0.5965,   0.7092,   0.8114,   0.8931,   0.9490,   
0.9806,   0.9944,   0.9989
		};

static double	tab3[] = { /* table for r == 3 */
0.9806,   0.9944,   0.9989,   0.0001682,0.0002542,0.0003829,0.0005745,0.0008587,
0.001278, 0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240,  0.01770,  
0.02505,  0.03514,  0.04880,  0.06704,  0.09103,  0.1220,   0.1612,   0.2097,   
0.2682,   0.3368,   0.4145,   0.4994,   0.5881,   0.6765,   0.7596,   0.8326,   
0.8922,   0.9367,   0.9667,   0.9846,   0.9939,   0.9980
		};

static double	tab4[] = { /* table for r == 4 */
2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,3.302e-06,4.990e-06,
7.524e-06,1.132e-05,1.698e-05,2.541e-05,3.791e-05,5.641e-05,8.368e-05,0.0001237,
0.0001823,0.0002677,0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
0.003494, 0.004942, 0.006948, 0.009702, 0.01346,  0.01853,  0.02532,  0.03431,
0.04607,  0.06128,  0.08068,  0.1051,   0.1352,   0.1719,   0.2157,   0.2669,
0.3254,   0.3906,   0.4612,   0.5355,   0.6110,   0.6849,   0.7544,   0.8168,
0.8699,   0.9127,   0.9451,   0.9679,   0.9827,   0.9915,   0.9963
		};

static double *table[] = { tab2, tab3, tab4 };
static short tabsize[] = { DIM(tab2)-1, DIM(tab3)-1, DIM(tab4)-1 };

static double f(double,void *);
static double g(double,void *);

/*
    Estimate the Sum P-value by calculation or interpolation, as appropriate.
	Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
	r = number of segments
	s = total score (in nats), adjusted by -r*log(KN)
*/
double GBlastSumP(Int4 r, double s)
{
	Int4		i, r1, r2;
	double	a;

	if (r == 1) return -GNlm_Expm1(-exp(-s));
	if (r <= 4) {
		if (r < 1) return 0.;
		r1 = r - 1;
		if (s >= r*r+r1) {
			a = GNlm_LnGammaInt(r+1);
			return r * exp(r1*log(s)-s-a-a);
		}
		if (s > -2*r) {
			/* interpolate */
			i = (Int4) (a = s+s+(4*r));
			a -= i;
			i = tabsize[r2 = r - 2] - i;
			return a*table[r2][i-1] + (1.-a)*table[r2][i];
		} return 1.;
	} return GBlastSumPCalc(r, s);
}

/***************************************************************************
    GBlastSumPCalc

    Evaluate the following double integral, where r = number of segments
    and s = the adjusted score in nats:

                    (r-2)         oo           oo
     Prob(r,s) =   r              -            -   (r-2)
                 -------------   |   exp(-y)  |   x   exp(-exp(x - y/r)) dx dy
                 (r-1)! (r-2)!  U            U
                                s            0
 ***************************************************************************/
static double GBlastSumPCalc(int r, double s)
{
	int		r1, itmin;
	double	t, d, epsilon;
	double	est_mean, mean, stddev, stddev4;
	double	xr, xr1, xr2, logr;
	double	args[6];

	epsilon = GBLAST_SUMP_EPSILON_DEFAULT; /* accuracy for SumP calcs. */
	if (r == 1) {
		if (s > 8.) return exp(-s);
		return -GNlm_Expm1(-exp(-s));
	}
	if (r < 1) return 0.;
	xr = r;
	if (r < 8) { if (s <= -2.3*xr) return 1.; }
	else if (r < 15) { if (s <= -2.5*xr) return 1.; }
	else if (r < 27) { if (s <= -3.0*xr) return 1.; }
	else if (r < 51) { if (s <= -3.4*xr) return 1.; }
	else if (r < 101) { if (s <= -4.0*xr) return 1.; }

	/* stddev in the limit of infinite r, but quite good for even small r */
	stddev = sqrt(xr);
	stddev4 = 4.*stddev;
	xr1 = r1 = r - 1;

	if (r > 100) {
		/* Calculate lower bound on the mean using inequality log(r) <= r */
		est_mean = -xr * xr1;
		if (s <= est_mean - stddev4)
			return 1.;
	}

	/* mean is rather close to the mode, and the mean is readily calculated */
	/* mean in the limit of infinite r, but quite good for even small r */
	logr = log(xr);
	mean = xr * (1. - logr) - 0.5;
	if (s <= mean - stddev4) return 1.;
	if (s >= mean) { t = s + 6.*stddev; itmin = 1; }
	else { t = mean + 6.*stddev; itmin = 2; }

#define ARG_R args[0]
#define ARG_R2 args[1]
#define ARG_ADJ1 args[2]
#define ARG_ADJ2 args[3]
#define ARG_SDIVR args[4]
#define ARG_EPS args[5]

	ARG_R = xr;
	ARG_R2 = xr2 = r - 2;
	ARG_ADJ1 = xr2*logr - GNlm_LnGammaInt(r1) - GNlm_LnGammaInt(r);
	ARG_EPS = epsilon;

	do {
		d = GNlm_RombergIntegrate(g, args, s, t, epsilon, 0, itmin);
#ifdef GBLASTKAR_HUGE_VAL
		if (d == GBLASTKAR_HUGE_VAL) return d;
#endif
	} while (s < mean && d < 0.4 && itmin++ < 4);
	return (d < 1. ? d : 1.);
}

static double g(double	s, void *vp)
{
	register double *args = (double *) vp;
	double	mx;
	
	ARG_ADJ2 = ARG_ADJ1 - s;
	ARG_SDIVR = s / ARG_R;	/* = s / r */
	mx = (s > 0. ? ARG_SDIVR + 3. : 3.);
	return GNlm_RombergIntegrate(f, vp, 0., mx, ARG_EPS, 0, 1);
}

static double f(double	x, void *vp)
{
	register double *args = (double *)vp;
	register double	y;

	y = exp(x - ARG_SDIVR);
#ifdef GBLASTKAR_HUGE_VAL
	if (y == GBLASTKAR_HUGE_VAL) return 0.;
#endif
	if (ARG_R2 == 0.) return exp(ARG_ADJ2 - y);
	if (x == 0.) return 0.;
	return exp(ARG_R2*log(x) + ARG_ADJ2 - y);
}

//	Calculates the p-value for alignments with "small" gaps (typically
//	under fifty residues/basepairs) following ideas of Stephen Altschul's.
//	"gap" gives the size of this gap, "gap_prob" is the probability
//	of this model of gapping being correct (it's thought that gap_prob
//	will generally be 0.5).  "num" is the number of HSP's involved, "sum" 
//	is the "raw" sum-score of these HSP's. "subject_len" is the (effective)
//	length of the database sequence, "query_length" is the (effective) 
//	length of the query sequence.  min_length_one specifies whether one
//	or 1/K will be used as the minimum expected length.
double GBlastSmallGapSumE(kbp_typ kbp, Int4 gap, double gap_prob, 
	double gap_decay_rate, Int2 num, Int4 sum, double score_prime, 
	Int4 query_length, Int4 subject_length, Boolean min_length_one)
{

	double sum_p, sum_e;
		
	score_prime -= kbp->logK + log((double)subject_length
		*(double)query_length) + (num-1)*(kbp->logK + 2*log((double)gap));
	score_prime -= log(GNlm_Factorial((int) num)); 
	sum_p = GBlastSumP(num, score_prime);
	sum_e = GBlastKarlinPtoE(sum_p);
	sum_e = sum_e/((1.0-gap_decay_rate)*GNlm_Powi(gap_decay_rate, (num-1)));
	if (num > 1) {
		if (gap_prob == 0.0) sum_e = INT4_MAX;
		else sum_e = sum_e/gap_prob;
	}
	return sum_e;
}

//	Calculates the p-values for alignments with "large" gaps (i.e., 
//	infinite) followings an idea of Stephen Altschul's.
double GBlastLargeGapSumE(kbp_typ kbp, double gap_prob, 
	double gap_decay_rate, Int2 num, Int4 sum, double score_prime, 
	Int4 query_length, Int4 subject_length, Boolean old_stats)
{

	double	sum_p, sum_e;
	/* The next two variables are for compatability with Warren's code. */
	double	lcl_subject_length, lcl_query_length;

        lcl_query_length = (double) query_length;
        lcl_subject_length = (double) subject_length;

	score_prime -= num*(kbp->logK + log(lcl_subject_length*lcl_query_length)) 
			- log(GNlm_Factorial((int) num)); 
	sum_p = GBlastSumP(num, score_prime);
	sum_e = GBlastKarlinPtoE(sum_p);
	sum_e = sum_e/((1.0-gap_decay_rate)*GNlm_Powi(gap_decay_rate, (num-1)));
	if (num > 1) {
		if (gap_prob == 1.0) sum_e = INT4_MAX;
		else sum_e = sum_e/(1.0 - gap_prob);
	}
	return sum_e;
}

/********************************************************************
*
*	The following function, from Stephen Altschul, calculates a 
*	pseudoscore from a p-vlaue and n, the product of the database
*	and query lengths.
*	double	p;		 p-value	
*	double	n;		 search space 
*********************************************************************/
/* Move the following constant into blast.h??, or only the last one. */
#define		PSCALE 20.0
#define 	PSEUDO_SCORE_MAX 32767
#define 	SYBASE_MIN 1.0e-300

Int2 GConvertPtoPseudoS(double p, double n)
{
	Int2	s;
	double	E;

/* If p is 1.0, then E is very large and E/n is about one. */
	if (p > 0.99) return 0.5;
/* If p is very small, the highest score should be returned. */
	else if (p < SYBASE_MIN) return PSEUDO_SCORE_MAX;
/* E = -ln(1-p); the following calculates it. */
        E = -GNlm_Log1p(-p);
	s= GConvertEtoPseudoS (E, n);
	return s;
}

/*******************************************************************
*
*	This function calculates a pseudo-score from an E value.
*	The searchsp is the product of the query and database
*	lengths.  As the E value is directly related to the search
*	space, this effectively scales the database size out of
*	the calculation of the pseudo-score.
*******************************************************************/
Int2 GConvertEtoPseudoS(double E, double searchsp)
{
	Int2	s;

/* If E is very small, the highest score should be returned. */
	if (E < SYBASE_MIN) return PSEUDO_SCORE_MAX;
	if (E > searchsp) return 0;
/* 0.5 is added to make sure this is rounded up, not down. */
	s= 0.5-PSCALE*log(E/searchsp);
	return s;
}
 
/*
Given a pseudoscore, a subroutine for calculating an E-value for a comparison
of size n (the product of the sequence length) is:
*/
double GConvertPseudoStoE(Int2 s, double n) { return n*exp(-s/PSCALE); }

sbp_typ	GBLAST_ScoreBlkNew(a_type AB, unsigned char *query, Int4 query_length,
	Int4 gap_open, Int4 gap_extend, Int4 dblen, Int4 dbseq_num)
// dblen=TotalSeqSet(seqset);	dbseq_num=NSeqsSeqSet(seqset);
// query_length = LenSeq(queryE); query = SeqPtr(queryE) + 1;
{
        Int4            dblen_eff;      // effective length of the database.
        Int4            length_adjustment=0; // amount removed from query & db sequences.
        Int4            effective_query_length;  // adjusted length of query.

	sbp_typ sbp=GBLAST_ScoreBlkNew2(AB); // info on scoring.
	// following from BLASTSetUpSearchInternalByLoc( ) in blast.c
	// assumes matrix have been allocated above....
        if(GBlastScoreBlkFill(sbp,query,query_length)) 
		print_error("Karlin parameter calculation failed");
        sbp->kbp_gap_std[0] = GBlastKarlinBlkCreate();
        assert(!GBlastKarlinBlkGappedCalc(sbp->kbp_gap_std[0],gap_open,gap_extend));
        sbp->kbp_gap_psi[0] = GBlastKarlinBlkCreate();
        assert(!GBlastKarlinBlkGappedCalc(sbp->kbp_gap_psi[0],gap_open,gap_extend));
        sbp->kbp_gap = sbp->kbp_gap_std; sbp->kbp = sbp->kbp_std;
	//******************* WHAT IS THIS FOR??? *** NEEDED????
	Int4	*open, *extend,index;
	double	*lambda, *K, *H;
        Int4 array_size=GBlastKarlinGetMatrixValues(&open,&extend, 
				&lambda,&K,&H,NULL);
        if (array_size > 0) {
                for(index=0; index<array_size; index++) {
                      if(open[index] == INT2_MAX && extend[index] == INT2_MAX) {
                         sbp->kbp_ideal = GBlastKarlinBlkCreate();
                         sbp->kbp_ideal->Lambda = lambda[index];
                         sbp->kbp_ideal->K = K[index];
                         sbp->kbp_ideal->H = H[index];
                      }
                }
                GMemFree(open); GMemFree(extend);
                GMemFree(lambda); GMemFree(K); GMemFree(H);
        }
        if (sbp->kbp_ideal == NULL)
                sbp->kbp_ideal = GBlastKarlinBlkStandardCalcEx(sbp);
	//**************************** COMPUTE EFFECTIVE SearchSpace *****************
	Int4	full_query_length;
        Int4    length,last_length_adjustment,min_query_length;

	full_query_length=query_length;
	length = query_length;
	
        min_query_length = 1/(sbp->kbp_gap_std[0]->K); // for gapped_calculation
        last_length_adjustment = 0;
        for(index=0; index<5; index++) {
                length_adjustment = ((sbp->kbp[0]->logK)
                  +log((double)(length-last_length_adjustment)*(double)(
		  MAX(1,(dblen)-(dbseq_num*last_length_adjustment))
		  )))/(sbp->kbp[0]->H);
                if (length_adjustment >= length-min_query_length) {
                        length_adjustment=length-min_query_length; break;
                }
                if(ABS(last_length_adjustment-length_adjustment) <= 1) break;
                last_length_adjustment = length_adjustment;
        }
        length_adjustment=MAX(length_adjustment, 0);
        dblen_eff=MAX(1, dblen - dbseq_num*length_adjustment);
        effective_query_length=MAX(length - length_adjustment, min_query_length);
        effective_query_length= effective_query_length;
        sbp->searchsp_eff=((double) dblen_eff)*((double) effective_query_length);
	return sbp;
}

