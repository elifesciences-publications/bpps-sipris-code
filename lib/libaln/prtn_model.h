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

/* prtn_model.h - codes and constants for scan model. */
#if !defined (_PRTN_MODEL)
#define _PRTN_MODEL
#include "stdinc.h"
#include "afnio.h"
#include "scanheap.h"
#include "probability.h"
#include "histogram.h"
#include "sequence.h"
#include "smatrix.h"
#include "wmodel.h"
#include "pseg.h"
#include "guide.h"
#include "residues.h"
// #include "spouge.h"
#include "sma.h"
#include "dheap.h"
#include "fmodel.h"

typedef struct {
	a_type  A;
	BooLean	weights;	/** weight input sequences **/
	BooLean	segmask;	/** check for compositional bias? **/
	BooLean	nonglobular;	/** check for potential nonglobular regions? **/
	/****  SPOUGE GAPS ****/
	char	mode;		/** mode for spouge functions **/
	BooLean gapfunct;	/** use a gap function for search **/
	Int4	**observedGap;	/** gap lengths **/
	Int4	**scoreGap;	/** log-odds gap scores **/
	Int4	maxLength;	/** maximum sequence length **/
	Int4	totGaps;	/** total number of gaps **/
	/****  SPOUGE GAPS ****/
	// char	*snfile;
	Int4	N;		/** number of motif models **/
	Int4	maxrpts;	/** maximum number of repeats **/
	double	*freq;
	float	minmap;		// minimum block map for inclusion in model.
	wm_type  *M;		/** weighted model **/
	wm_type  *M2;		// alternative model (ignoring '^') **/
	smx_typ smx,*sM,*fsM;	// various models and matrices 
	Int4	*w;
	Int4	max_gap_score;
	Int4	tot_len;
	// js_type	Spg;
	char	method;		// method used for model: g,c,r
	/** TEST ***/
	Int4	*X[2];
	double	*temp;
	// OVERLAPS 
	Int4	tot_lenOL; // total compressed length of models
	Int4	totalOL;   // total overlaps for models
	char	*overlap;	// overlaps between models
	sma_typ	MA;
} prtn_model_type;
typedef prtn_model_type *ptm_typ;

/********************************* private ********************************/
BooLean ReadPrtnModelSMA(FILE *fptr, double pseudo,double *freq, double minprob,
	ptm_typ F);
void    AllocatePrtnModel(Int4 *,sma_typ,Int4, ptm_typ);

float   **edit_info_protein_model(float **info0, ptm_typ PM);
/********************************* PUBLIC ********************************/
ptm_typ *MakePrtnModels(char *snfile, a_type A, Int4 maxrpts, char method,
        char mode, Int4 maxLength, BooLean weight,float minmap,
        double pseudo, double minprob, double *freq);
ptm_typ MakePrtnModel(char *snfile, a_type A, Int4 maxrpts, char method,
	char mode, Int4 maxLength, BooLean weight, float minmap,double pseudo,
	double minprob, double *freq);
ptm_typ MakePrtnModel(FILE *fptr, a_type A, Int4 maxrpts, char method,
        char mode, Int4 maxLength, BooLean weight,float minmap,
        double pseudo, double minprob, double *freq);
BooLean CheckPrtnModel(Int4 *Score, e_type E,UInt8 total, double maxEval, ptm_typ F);
Int4    ComparePrtnModel(e_type E,UInt8 total,double maxEval,Int4 *p,float *pv,
        double *pvalue,ptm_typ F);
Int4    ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue,
        Int4 open, Int4 extend,ptm_typ F);
Int4    ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue,
        Int4 open, Int4 extend,ptm_typ F,char mode);
Int4    ComparePrtnModelRpts(e_type E,UInt8 total,double maxEval,Int4 *p,float *pv,
        double *pvalue,double rptsEval,ptm_typ F);
Int4    ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue,
	ptm_typ F);
wm_type MergePrtnModels(e_type E, ptm_typ PM);
void    PutSmxPrtnModel(Int4 t, ptm_typ PM);
void    NilPrtnModel(ptm_typ F);
double  SWHistPrtnModel(FILE *fp, Int4 a, Int4 b, e_type E, ptm_typ PM);
double	LocalHistPrtnModel(FILE *fp,e_type E, ptm_typ PM);
double	HistPrtnModel(FILE *fp,e_type E, ptm_typ PM);
Int4	PutSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM);
double  PutSWAlnPrtnModelRpts(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM);
Int4    PlotCRSPrtnModelSW(FILE *fptr, Int4 open,Int4 extend, double info_cutoff,
                        e_type E, fm_type *FM, ptm_typ PM,char *InputColors);
Int4    PlotCRSPrtnModelSW(FILE *fptr, Int4 open,Int4 extend, double info_cutoff, 
		                        e_type E, fm_type *FM, ptm_typ PM);
Int4    PutCRSSeqAlnSMatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2,
	  UInt4 offset, Int4 J, Int4 nmod, smx_typ *M,char color, Int4 m_use,
	  	Int4 rpt, fm_type FM);
BooLean PutSeqAlnPrtnModel(FILE *fp, char *operation, Int4 start, e_type E, ptm_typ PM);
double  PutSampledSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM);
double  PutFullSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM);
double  ShortGapAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM);
Int4	ScoreSWAlnPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);
Int4    LocalScoreSWAlnPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);

Int4    GapScorePrtnModel(Int4 block, Int4 gaplength, ptm_typ PM);
void    PutGapScoresPrtnModel(FILE *fp, ptm_typ PM);
char    *GapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);
char    *SampleGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);
char    *RptsGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);
char    *RptsSampleGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);
void    CompareScoresPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM);

float   **ExcessInfoPrtnModel(char *string, Int4 purge_cutoff, ptm_typ PM);
float   **InfoPrtnModel(Int4 purge_cutoff, ptm_typ PM);
double	*ObservedPrtnModel(Int4 m,Int4 i,ptm_typ PM);
double  *ObservedFilterPrtnModel(Int4 m,Int4 i,ptm_typ PM);

double  EvaluePrtnModel(Int4 score, UInt4 len, ptm_typ PM);

/********************************* MACROS ********************************/
#define NumSeqsPrtnModel(P)    		nseqSMA((P)->MA)
#define NumBlksPrtnModel(P)    		ntypSMA((P)->MA)
#define	BackgroundFreqPrtnModel(P)	((P)->freq)
#define SetMethodPrtnModel(x,F)		((F)->method = (char)(x))
#define NoMaskPrtnModel(F)		((F)->segmask=FALSE)
#define MaskNonGlobularPrtnModel(F)	((F)->nonglobular=TRUE)
#define NumModelsPrtnModel(F)    	((F)->N)
#define PrtnModelLengths(F)    		((F)->w)
#define PrtnModelA(F)    		((F)->A)
#define BlkMinMapPrtnModel(F)    	((F)->minmap)
#define ModePrtnModel(F)    		((F)->mode)
#define TotLenPrtnModel(F)    		((F)->tot_len)
#define NmaxPrtnModel(F)		((F)->maxrpts*(F)->N)
#define MaxRptsPrtnModel(F)		((F)->maxrpts)
#define maxLengthPrtnModel(F)		((F)->maxLength)
#define GapScoresPrtnModel(F)		((F)->scoreGap)
#define SMatricesPrtnModel(F)		((F)->sM)
#define SMatrixPrtnModel(m,F)		(((m) <= (F)->N && (m) > 0)?\
					GetSMatrixWModel((F)->M[(m)]): NULL)
#define WModelPrtnModel(m,F)		(((m) <= (F)->N && (m) > 0)? (F)->M[(m)]: NULL)
// don't give this out!! #define smaPrtnModel(F)	((F)->MA)
#define smaPrtnModel(F)			((F)->MA) // why not...
#define	DescriptionsmaPrtnModel(F)	DescriptionSMA((F)->MA)

#endif

