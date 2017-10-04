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

#if !defined (SMATRIX)
#define	SMATRIX
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "blosum62.h"
#include "probability.h"
#include "histogram.h"
#include "sequence.h"
#include "idp_typ.h"
/*************************** ADT PROFILE ***************************
	
			scoring matrix data type

**********************************************************************/

/*************************** SMatrix Type **************************/
typedef struct {
	Int4	*cmax,*max;		/* maximum */
	Int4	*cmin,*min;		/* minimum */
	Int4	K;			/* length */
	Int4	**score;		/* scoring matrix */
	double	mean,var,sd;
	double	nsd;			/* minimum # std. dev. above mean */
	double	*f,*f0,*fx;
	double	*freq;			/* freq[r] */
	a_type	A;			/* alphabet */
	Int4	nlet;
	Int4	neginf;			/* lowest permissible score */
	BooLean	changed,calc_stats,calc_prob;	/* update data */
	/** gap penalties ***/
	Int4	*oi,*ei,*od,*ed;	// position specific gap penalties
	/** posMatrix fields ***/
	Int4	**posMatrix;		/** psi-blast Matrix **/
	Int4	posMtrxLen;		/** psi-blast Matrix length **/
	Int4	posMtrxOffset;		/** psi-blast Matrix offset **/
	unsigned char	*qseq;		/** query sequence for Matrix **/
} smatrix_type;
typedef smatrix_type *smx_typ;

/******************************* private *****************************/
void	smatrix_prob(register Int4 k, register double *f, smx_typ M);
Int4	min_max_smatrix(smx_typ M);
void	stats_smatrix(smx_typ M);
double  smatrix_prob_fast(register Int4 k, register double *f, smx_typ M);
void    smatrix_error(char *s);
Int4    local_score_smatrix(register unsigned char *seq, register Int4 j,
        register Int4 **scr);
Int4    score_smatrix(register unsigned char *seq, register Int4 j, 
        register Int4 **scr);
Int4    gap_func_aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *pos, Int4 **gapscore,
        Int4 (*score_smtrx)(unsigned char *, Int4 , smx_typ));
Int4    aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
        Int4 *pos, Int4 (*score_smtrx)(register unsigned char *, register Int4,
        register Int4 **));
Int4    fast_aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
        Int4 (*score_smtrx)(register unsigned char *, register Int4, 
	register Int4 **));

char    *gap_aln_trace_smatrixSW(Int4 a,Int4 b,Int4 n2,unsigned char *seq2,
        Int4 nmod,smx_typ *M,Int4 *J,Int4 *alnscore,Int4 *start);
char	*gapped_aln_seq_smatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *J, Int4 *alnscore);
char	*sample_gapped_aln_seq_smatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *J, Int4 *alnscore);
Int4    put_seqaln_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
        UInt4 offset, Int4 J, Int4 nmod, smx_typ *M,char mode);
Int4    put_seqaln_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
        UInt4 offset, Int4 J, Int4 nmod, smx_typ *M);
Int4    FastLocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *X[2]);
Int4    FastAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *X[2]);
Int4    gap_func_aln_seq_overlap_smatrix_Cterm(Int4 len,unsigned char *seq,
        Int4 nmod, smx_typ *M,Int4 *pos,Int4 **gapscore, char *overlap,
        Int4 (*score_smtrx)(unsigned char *,Int4,smx_typ));
char    *gapped_aln_trace_smatrixSW(Int4 a,Int4 b,Int4 n2,unsigned char *seq2,
        Int4 nmod,smx_typ *M,Int4 **gapscore,Int4 *J,Int4 *alnscore,Int4 *start);
/**************************** NEW: Spouge **********************/
Int4    extend_core_right_smatrix(register unsigned char *seq, register Int4 end,
        register Int4 **scr);
Int4    extend_core_left_smatrix(register unsigned char *seq, register Int4 j,
        register Int4 **scr);
Int4    local_core_smatrix(register unsigned char *seq, register Int4 len,
        register Int4 **scr);
/******************************* Public ******************************/
Int4    put_cmaseq_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2,
        UInt4 offset, Int4 J, Int4 nmod, smx_typ *M);
/***************************** operations ****************************/
smx_typ ReverseSMatrix(smx_typ smx);
void    PutSMatrix(FILE *fptr, smx_typ M);
smx_typ MkSMatrixN(Int4 N, Int4 K, double *freq, a_type A);
smx_typ MkSMatrix(double nsd, Int4 K, double *freq, a_type A);
void    NilSMatrix(smx_typ M);
Int4	MaxScoreSMatrix(smx_typ M);
Int4    MinScoreSMatrix(smx_typ M);
Int4    MaxSegSMatrix(unsigned char *maxseg, smx_typ M);
Int4    MinSegSMatrix(unsigned char *minseg, smx_typ M);
void	SetSMatrix(Int4 r, Int4 row , Int4 score, smx_typ M);
Int4	ScoreSMatrix(register unsigned char *seq, Int4 start, smx_typ M);
Int4	SplitScoreSMatrix(unsigned char *seq, Int4 n, Int4 *start, Int4 *leng, 
	smx_typ M);
double  SMatrixProb(Int4 score, smx_typ M);
double  SMatrixProbFast(Int4 score, smx_typ M);
Int4    AlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 *pos);
/*** NEW ****/
double  ExpScoreSMatrix(smx_typ M);
double  ExpScoreSMatrix(Int4 i,smx_typ M);
double  ExpScoreSMatrix(double **freq, smx_typ M);
/**************************** NEW: Spouge **********************/
Int4    LocalScoreSMatrix(register unsigned char *seq, register Int4 start, 
	register smx_typ M);
Int4    LocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 *pos);
Int4    GapFuncAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
        Int4 *pos, Int4 **gapscore);
Int4    LocalGapFuncAlnSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos, Int4 **gapscore);
/**************************** NEW: Spouge **********************/
/**************************** NEW: OVERLAPS **********************/
Int4    FastAlnSeqOverlapSMatrix(Int4 len, unsigned char *seq,
        Int4 nmod, char *overlap, smx_typ *M, Int4 *X[2]);
Int4    AlnSeqOverlapSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
        char *overlap, smx_typ *M, Int4 *pos);
Int4    GapFuncAlnSeqOverlapSMatrix(Int4 len,unsigned char *seq,
	Int4 nmod,smx_typ *M,Int4 *pos,Int4 **gapscore,char *overlap);
// private 
static Int4 fast_aln_seq_overlap_smatrix(Int4 len, unsigned char *seq,
        Int4 nmod, smx_typ *M, char *overlap,
        Int4 (*score_smtrx)(register unsigned char *, register Int4,
        register Int4 **));
static Int4 aln_seq_overlap_smatrix(Int4 len, unsigned char *seq,
        Int4 nmod, smx_typ *M, Int4 *pos, char *overlap,
        Int4 (*score_smtrx)(register unsigned char *, register Int4,
        register Int4 **));
Int4    SubScoreSMatrix(unsigned char *seq, Int4 start, Int4 smx_start,
        Int4 smx_end, smx_typ M);
Int4    OptInsertScoreSMatrix(unsigned char *seq, Int4 start, Int4 overlap,
        smx_typ M1, smx_typ M2);
Int4    BestAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *pos);
Int4    BestLocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *pos);

double	HistSeqSMatrix(FILE *fp, Int4 len, unsigned char *seq, smx_typ M);
double	LocalHistSeqSMatrix(FILE *fp, Int4 len, unsigned char *seq, smx_typ M);
/**************************** NEW: Gapped sequence **********************/
Int4    AlnSeqSMatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, Int4 nmod,
        smx_typ *M, Int4 *mpos, Int4 cutoff);
Int4    GapXDropSMatrix(Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        smx_typ M, Int4 *mpos, Int4 maxgap);
Int4    GappedAlnSeqSMatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore);
Int4    ScoreGappedAlnSeqSmatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore, BooLean local);
Int4    PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2,
        unsigned char *seq2, Int4 nmod, smx_typ *M, Int4 **gapscore);
Int4    PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore,Int4 offset);
Int4    PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore,Int4 offset,char mode);
void    PutGappedSeqAlnSMatrix(FILE *fp, char *operation, Int4 offset, Int4 lenseq,
        unsigned char *seq, Int4 nmod, smx_typ *M);
Int4    PutSampledSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2,
        UInt4 offset, Int4 nmod, smx_typ *M, Int4 **gapscore);
Int4    PutFullSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2,
        unsigned char *seq2, Int4 nmod, smx_typ *M, Int4 **gapscore);
char    *GapOperationsSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore);
char    *GapAlnTraceSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start);
char    *SampleGapOperationsSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore);
Int4    GapFuncAlnSeqOverlapSMatrixCterm(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *pos, Int4 **gapscore,char *overlap);
char    *GapAlnTraceSMatrix2(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
        Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start,Int4 *score);

// from ALEX's code.
char    *ComplexAlnSMatrix(idp_typ *idp, Int4 seq_len, unsigned char *seq,
	Int4 Rpts, smx_typ *M, Int4 *Score, Int4 *Start, Int4 *Oper_len);
/**************************** NEW: PSI-BLAST **********************/
Int4    **SMatrix2GPSI(smx_typ M, Int4 neg_inf);
smx_typ GPSI2SMatrix(Int4 length, int **posMatrix, a_type A);
char    *SampleOperationSeqAlnSMatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
        UInt4 offset, Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *J,
        Int4 *alnscore);
Int4    **BlastCircularSMatrix(Int4 seqlen, smx_typ M);
Int4    CircularSMatrixPos(Int4 seqlen, Int4 mtxpos, smx_typ M);
unsigned char   *CircularSMatrixQuery(smx_typ M);
Int4    SMatrixQueryLen(smx_typ M);
void    PutCircularSMatrix(FILE *fptr, smx_typ M);

Int4    PutSeqAlnSMatrixSW2(FILE *fp, char *operation, Int4 n2, unsigned char *seq2,
        UInt4 offset, smx_typ M, Int4 *start);
/**************************** NEW: libalex **********************/
// e_type ConsensusSMatrix(smx_typ M);
e_type ConsensusSMatrix(smx_typ *M, Int4 nbr, Int4 rpts);
/**************************** macro operations **********************/
#define LenSMatrix(M)		((M)->K)
#define SMatrixA(M)		((M)->A)
#define meanSMatrix(M)	(((M)->changed)? stats_smatrix(M),(M)->mean:\
				(M)->mean)
#define sdSMatrix(M)	(((M)->changed)? stats_smatrix(M),(M)->sd:(M)->sd)
#define ValSMatrix(pos,r,M)	((M)->score[(pos)][(r)])
#define ValuesSMatrix(M)	((M)->score)
#define NegInfSMatrix(M)	((M)->neginf)

#endif

