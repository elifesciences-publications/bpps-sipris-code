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

/******************************** cmsa.h - *******************************
   		Colinear Multiple Sequence Alignment Data Type:

  
*************************************************************************/
#if !defined(_CCMSA_)
#define _CCMSA_
#include "afnio.h"
#include "fmodel.h"
#include "sites.h"
#include "dheap.h"
#include "mheap.h"
#include "histogram.h"
#include "probability.h"
#include "residues.h"
#include "random.h"
#include "colwt.h"
#include "sma.h"
#include "mdl.h"
#include "guide.h"
#include "gss_typ.h"
#include "spouge.h"
#include "dom_typ.h"
#include "sset.h"
#include "set_typ.h"
#include "dsets.h"
#include "prtn_model.h"

/******************** colinear multiple alignment type ******************/
typedef struct colinearmaln {
	double  	**lngamma;      // array for lngamma[n][r].
        st_type         sites;          // working multiple alignment.
        sti_typ         best;           // archived best alignment.
	mdl_typ		*mdl;		// protein statistical model.
        BooLean         *null,**bestnull;	// column configurations
// temporary arrays:
        Int4            *pos;           // temporary array
        Int4            *maxlen;	// 
// subfamily alignments... 
	Int4		num_subaln;	// number of subfamily alignments=NumSeqCMSA.
	struct colinearmaln **next;	// link to subfamilies (NULL -> leaf node)
//	struct psiblst_aln  **sub;	// could also use psiblast aligns?
//	Added information about full sequences for subseqs.
	ss_type		FullSeq;	// full size sequences;
	unsigned short	*FullRpts;	// number of repeats in each full seq.
	Int4		*SubToFull;	// map subseq to fullseq.
	Int4		*FirstRpt;	// map fullseq to first subseq.
//	Added information about other domains withing sequences.
	dom_typ		*Domains;
#if 0
	// Want to also add grease, secondary struct, etc...
	idp_typ		*idp;		// insert & deletion penalties...
#endif
	Int4		Level;		// level in alignment hierarchy
#if 1	
	double		alpha;
	Int4		A0,B0;	// BPPS tuning parameters. (afn: 12/2/09).
	char		set_mode;
#endif
} colinearmaln_type;
typedef colinearmaln_type	*cma_typ;

/******************************** private ********************************/
cma_typ MkMSA(st_type S);
#define	MAX_FLANK_CMSA	5

Int4    slide_left_dfs_cmsa(Int4 t, double *new_map, cma_typ cma, Int4 depth);
Int4    slide_right_dfs_cmsa(Int4 t, double *new_map, cma_typ cma, Int4 depth);

Int4    blks4cross_cmsa(Int4 numcross, Int4 *aln1, Int4 *aln2,
        Int4 *cross1, Int4 *cross2);
void    calc_lngamma_cmsa(cma_typ C);
Int4    cmsa_fix_alignment(Int4 lenM, Int4 N,char **alignment, BooLean *null);
BooLean comp_map_cmsa(double nmap, double omap,double temperature);
cma_typ delete_blk_cmsa(cma_typ ma1, Int4 tdel);
void    dfs_config_cmsa(Int4 *aln1, Int4 *aln2, Int4 **config, Int4 depth,
        Int4 *nc, Int4 *rc, Int4 len, Int4 nrc);
Int4	**gap_score_cmsa(cma_typ cma);
cma_typ make_cmsa(cma_typ L);
BooLean move_column_msa(cma_typ G, Int4 lemon, Int4 t);
double  num_pseudo_msa(double N);
double  RecombinantMapCMSA(cma_typ ma1, cma_typ ma2, Int4 *config);
double  single_rel_map_cmsa(cma_typ L, Int4 t);
double  transfer_weight(double wd, double cd, double nwr, double owr,
        double cr);
Int4    sites_from_operations_cmsa(char *operation, Int4 start, Int4 *pos,
        Int4 *indels,gsq_typ *gsq);
Int4    create_gsq_cmsa(cma_typ cma);

void    fill_models_cmsa(cma_typ L);
void    delete_models_cmsa(cma_typ cma);
void	create_models_cmsa(cma_typ cma,BooLean **null);

/********************************* PUBLIC ********************************/
//--------------- make cmsa -----------------
cma_typ MakeCMSA(st_type S, BooLean **null);
cma_typ MakeCMSA(e_type *ListE,Int4 N,char **operation,Int4 *start,cma_typ cma);
cma_typ MakeCMSA(e_type *ListE,Int4 N,char **operation,Int4 *start,Int4 nblks,
        Int4 *lengths,Int4 open,Int4 extend,double pernats,Int4 leftflank,
        Int4 rightflank,char *name,a_type A);
cma_typ MakeCMSA(e_type *ListE, Int4 N, char **operation, Int4 *start, Int4 nblks,
        Int4 *lengths, Int4 gapo, Int4 gapx, double pernats, Int4 leftflank, 
        Int4 rightflank, char *name, a_type A, e_type *FullE, unsigned short *FullR);
cma_typ OneSeqToCMSA(e_type Seq, char *name, a_type AB);
cma_typ MakeColinearMSA(st_type S, BooLean **null);
cma_typ EmptyCMSA(Int4 N, Int4 *lengths, ss_type data,Int4 open,
        Int4 extend, double pernats, Int4 leftflank, Int4 rightflank);
cma_typ CopyCMSA(cma_typ msa);

//----------------- init cmsa ------------------
void	InitCMSA(cma_typ G);
void	InitMAPCMSA(cma_typ G);

//----------------- destroy cmsa ------------------
void    NilCMSA(cma_typ L);
void    TotalNilCMSA(cma_typ L);
BooLean **NullCMSA(cma_typ ma);
BooLean	NullSiteCMSA(Int4 t,Int4 s,cma_typ cma);

//============== debugging ===================
BooLean SamePerSeqCMSA(cma_typ L);
void    PutTableOfCMSA(FILE *fptr, cma_typ cma);	// with no more than 100 seqs in cma.
//================ Information routines ===============

//------------------ probabilities ---------------
double  ProbFModelCMSA(Int4 t, Int4 n, Int4 s, cma_typ L);
double  GetProbCMSA(Int4 t, Int4 n, cma_typ L);
double  GetTotalProbCMSA(Int4 n, cma_typ L);
double  GetGappedProbCMSA(UInt4 t, UInt4 n, cma_typ cma);
double  GetTotalGappedProbCMSA(Int4 n, cma_typ cma);

//------------------- indels ------------------
double	ExpectedGapLengthCMSA(cma_typ cma, Int4 blk);
double  GetGappedSqLPR_CMSA(Int4 s, cma_typ cma);
BooLean	IsDeletedCMSA(UInt4 n, UInt4 r, cma_typ cma);
BooLean IsDeletedCMSA(UInt4 blk, UInt4 n, UInt4 r, cma_typ cma);
unsigned short InsertionCMSA(UInt4 blk,UInt4 n,UInt4 i, cma_typ cma);
unsigned short InsertionCMSA(UInt4 n,UInt4 i,cma_typ cma);
Int4    NumInsertsCMSA(Int4 blk, Int4 site, cma_typ cma, Int4 &ins_res);
Int4    ResBetweenCMSA(register Int4 t, register cma_typ cma);
Int4    ResBetweenCMSA(register Int4 t, register Int4 n, register cma_typ cma);
double   FractDeletionsCMSA(Int4 blk, Int4 pos, cma_typ cma);
Int4     NumDeletionsCMSA(Int4 blk, Int4 pos, cma_typ cma);

//------------------ blocks and columns -----------
Int4    ConfigCMSA(Int4 *ncols, cma_typ L);
Int4    ConfigCMSA2(Int4 *width, cma_typ L);
void	PutConfigCMSA(FILE *fp, cma_typ cma);
Int4    NumColumnsCMSA(cma_typ msa);
Int4    TotalLenCMSA(cma_typ msa);
float   *RelEntropyCMA(Int4 blk, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, BooLean *skip, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, Int4 ***observed, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, BooLean *skip, Int4 ***observed, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, double ***observed, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, BooLean *skip, double ***observed, cma_typ cma);
double  **ColResFreqsCMSA(Int4 t, set_typ Set, Int4 ***observed, cma_typ cma);

//----------------- sequence features -----------------
Int4    ResidueCMSA(register Int4 t, register Int4 n, register Int4 s, cma_typ cma);
Int4    FastDangerousResCMSA(register Int4 t, register Int4 n, register Int4 s, cma_typ cma);
Int4    SetResidueCMSA(Int4 t, Int4 n, Int4 s, unsigned char r, cma_typ cma);

//----------------- sequence positions -----------------
BooLean	SeqIsInCMSA(Int4 sq,cma_typ cmsa);
void    PutAlnSeqsCMSA(FILE *fp,cma_typ cma);
Int4    PosSiteCMSA(Int4 t, Int4 n, Int4 *pos, cma_typ cmsa);
Int4    RealToFakeCMSA(Int4 sq, Int4 s, cma_typ cma);
Int4    FakeToRealCMA(Int4 sq,Int4 s, cma_typ cma);
Int4    TruePosCMSA(Int4 sq, Int4 s, cma_typ cma);
Int4    TruePosCMSA(Int4 sq, Int4 blk, Int4 s, cma_typ cma);

//----------------- consensus sequences -----------------
char    *ConsensusSeqCMSA(Int4 t, cma_typ cma);
char    *ConsensusSeqCMSA(Int4 t, set_typ Set, cma_typ cma);
char    *ConsensusSeqCMSA(cma_typ cma);
e_type	MkConsensusCMSA(set_typ Set,cma_typ cma);
e_type	MkConsensusCMSA(cma_typ cma);
cma_typ AddConsensusCMSA(cma_typ cma);
Int4    ConsensusScoreCMSA(cma_typ cma, e_type E, Int4 gapo, Int4 gapx, Int4 sq);
cma_typ RtnConSqAsCMSA(char *name,sst_typ *xsst,e_type keyE,a_type AB);

//----------------- sequences-to-sequence scores  -----------------
Int4    PseudoAlnScoreCMSA(Int4 sq1, Int4 sq2, cma_typ cma);
double  PseudoAlnScoreCMSA(unsigned char *qseq, Int4 sq, Int4 *Rank, cma_typ cma);
Int4    PseudoAlnScoreSqToCMSA(e_type E, Int4 sq2, cma_typ cma);
Int4    PseudoAlnScoreTwoCMSA(Int4 sq1, cma_typ cma1, Int4 sq2, cma_typ cma2);

//-------------------- misc. routines -------------------
void    PutMinimalSeqCMA(FILE *fp, Int4 sq, cma_typ cma);
cma_typ MinimizeFirstSeqCMSA(cma_typ cma);

Int4	SubToFullCMSA(Int4 sq, cma_typ cma);
smx_typ SMXforSeqCMSA(Int4 t,Int4 sq,cma_typ cma);
e_type  MaxWordAndSubSeqCMSA(Int4 t, Int4 qsq, Int4 ssq, Int4 *Word_score,
                Int4 *Q_start, Int4 *S_start, Int4 *Ssq_start, cma_typ cma);
Int4    RepeatsInfoCMSA(Int4 *start,Int4 *end,Int4 *gap0,Int4 *gap1,Int4 sq, 
		cma_typ cma);

double *HenikoffWeightsCMSA(cma_typ cma);
double *HenikoffWeightsCMSA(cma_typ cma,Int4 ***RawCnts);

//---------------- residue diversity ------------------------
double  ResidueDiversityCMSA(FILE *fp,cma_typ cma);
double  ResidueDiversityCMSA(FILE *fp,BooLean *skip,cma_typ cma);
double  ResidueDiversityCMSA(register BooLean *skip, register cma_typ cma);
double  ResidueDiversityCMSA(register Int4 *sqid, register cma_typ cma);
double  ResidueDiversityNotCMSA2(set_typ set, cma_typ cma);
double  ResidueDiversityCMSA2(set_typ set, cma_typ cma);
double  ResidueDiversityCMSA(set_typ set, cma_typ cma); // WARNING: fix

unsigned char   *GetAlnResInSiteCMSA(Int4 t, Int4 sq, cma_typ cma);

//============== input and output routines cmsa_io.cc ==============
void    WriteCMSA(char *filename, cma_typ cma);
void    PutCMSA(FILE *fptr, cma_typ cma);

//--------------- put merged routines ---------------------
void    PutMergedCMSA(FILE *fp,unsigned short nsets,cma_typ *cma);
void    PutMergedCMSA(FILE *fp,unsigned short nsets,set_typ *set, cma_typ *cma,sst_typ *sstP);
void    PutMergedCMSA(FILE *fp,char *name, unsigned short nsets,set_typ *set, cma_typ *cma,
		sst_typ *sstP);
void    PutMergedCMSA(FILE *fp,char *name, unsigned short nsets,set_typ *set, cma_typ *cma,
		sst_typ *sstP,char Labeling);

//--------------- put good routines ---------------------
void    WriteGoodCMSA(char *filename,double cutoff, cma_typ cma);
void    PutGoodCMSA(FILE *fp,double cutoff, cma_typ cma);
void    PutBestCMSA(FILE *fp,Int4 Num, BooLean KeepFirst, cma_typ cma);
void    PutGoodRptsCMSA(FILE *fp,Int4 min_rpt,Int4 min_spacing,cma_typ cma);

//--------------- put select routines ---------------------
Int4	PutSelectCMSA(FILE *fp,BooLean *skip, cma_typ cma);
Int4	PutSelectCMSA0(FILE *fp,BooLean *skip, cma_typ cma);
Int4	PutSelectCMSA0(FILE *fp,BooLean *skip, BooLean put_csq, cma_typ cma);
Int4	PutSelectOneCMSA(FILE *fp,BooLean *skip, cma_typ cma);
Int4	PutSelectOneCMSA(FILE *fp,BooLean *skip, Int4 *sortedlist, cma_typ cma);
void    PutSelectCMSA(FILE *fp,FILE *fp2, const char *name1,const char *name2,
        set_typ Set, cma_typ cma);
void    PutSelectCMSA(FILE *fp,FILE *fp2, set_typ Set, cma_typ cma);

void    PutInSetCMSA(FILE *fp,set_typ set, cma_typ cma);
void    PutInSetCMSA(FILE *fp,set_typ set, sst_typ *sstP,cma_typ cma);
void    LabelPutInSetCMSA(FILE *fp,set_typ set, cma_typ cma);
void    UnLabelPutInSetCMSA(FILE *fp,set_typ set, cma_typ cma);
void    PutSeqsInSetCMSA(FILE *fp,set_typ set, cma_typ cma);

Int4    *SortByQueryOneCMSA(cma_typ cma);

//------------------- read in one or more cmsa -----------------
cma_typ	ReadCMSA(FILE *fptr,a_type A);
cma_typ	ReadCMSA(FILE *fptr,char **column_status,a_type A);
cma_typ	ReadCMSA2(char *filename,a_type A);
cma_typ *MultiReadCMSA(FILE *fp,Int4 *Number,char ***status, a_type A);
cma_typ *MultiReadCMSA(FILE *fp,Int4 *Number,a_type A);

void    PutConsensusCMSA(FILE *fp, cma_typ cma);
void    PutCsqWithCMSA(FILE *fp,cma_typ cma);

//------------------- output as fasta -----------------
void    PutFastaAlnCMSA(FILE *fp,cma_typ cma);
Int4	PutFastaCMSA(FILE *fp,cma_typ cma, BooLean add_consensus=FALSE);
void    PutStockholmCMSA(FILE *fp,cma_typ cma,set_typ set=0);

Int4	PutFullSeqCMSA(FILE *fp,BooLean *skip, cma_typ cma);
void    PutDiffSetsCMSA(char *filename, Int4 set_id,set_typ SubSet,set_typ SuperSet,
                cma_typ cma);	// not used...?

// misplaced?
cma_typ MakeConsensusCMSA(cma_typ cma);
double  AveRelEntropyCMSA(BooLean *skip, cma_typ cma);
double  AveRelEntropyCMSA(set_typ set, cma_typ cma);
double  AvePercentIdentityCMSA(register set_typ set, register unsigned char *qsq,
        	register cma_typ cma);
cma_typ MkMainFileCMSA(cma_typ cma, Int4 num_random,cma_typ &rcma);
//************** cmsa_io.cc **************

//*********************** cmsa.cc ********************

//*********************** cmsa_put.cc ********************
void    PutGoodAlnCMSA(char *name, cma_typ L,double cutoff, gd_type G); // not used...?
void    PutRepSetCMSA(FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma);
void    PutRepSetCMSA(FILE *fp_err,FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma);
BooLean	*RtnRepSetCMSA(FILE *fp_err,Int4 percent_ident,Int4 *Nset,cma_typ cma);
void    PutBlockSpacingCMSA(FILE *fp, Int4 block, cma_typ cma);
//*********************** cmsa_put.cc ********************

// From CHAIN analysis routines (these use set_typ.h)
cma_typ GetInSetCMSA(set_typ set, cma_typ cma);
cma_typ GetBestCsqCMSA(cma_typ mcma);
cma_typ GetBestCsqAsCMSA(set_typ set, cma_typ cma);
e_type	GetBestConSqCMSA(cma_typ mcma);
e_type  GetSeqAsCsqCMSA(Int4 i, cma_typ cma);
e_type  GetSeqAsCsqCMSA(set_typ set, cma_typ cma);
void    LabelSeqsCMSA(cma_typ cma);
void    UnLabelSeqsCMSA(cma_typ cma);

// Domain information for FullSeqs
void	AddDomainsCMSA(dom_typ *domains, cma_typ cma);
void    CopyDomainsCMSA(cma_typ cma2, cma_typ cma);
void	RmDomainsCMSA(cma_typ cma);

// CONVERSION FROM SMA (cmsa_sma.c == sma_typ )
cma_typ SMA2CMSA(char *fafile, sma_typ MA);
cma_typ SubSMA2CMSA(char *fafile, BooLean *remove, sma_typ MA);
cma_typ SubSMAtoCMSA(ss_type data, BooLean *remove, sma_typ MA);
cma_typ SMAtoCMSA(ss_type data, sma_typ MA);

//============ cmsa alignment routines ========
Int4    *AlignSeqCMSA(FILE *fptr,Int4 *Score, e_type  E,cma_typ cmsa);

//============ cmsa_operations.c:  stochastic & algebraic operations) ========

set_typ RtnFastRepSetCMSA(FILE *fp_err, Int4 percent_ident,set_typ InSet,cma_typ cma);

//-------------------- adding a aligned block to a sequence ----------------------
void    AddSiteCMSA(Int4 t,Int4 n,Int4 s, cma_typ L);
void    RmSiteCMSA(Int4 t,Int4 n,Int4 s, cma_typ L);
void    VacateSitesCMSA(Int4 n,cma_typ L);

//-------------------- realigning a sequence ----------------------
void	ReplaceCMSA(Int4 s, gsq_typ *gsq, cma_typ ma);
gsq_typ *SwapGsqCMSA(Int4 s, gsq_typ *gsq, cma_typ cma, BooLean swap=TRUE);
BooLean ReAlignGSqCMSA(Int4 s,char *operation, Int4 start, cma_typ *oldcma);

//-------------------- adding, removing or inserting columns ----------------------
Int4    AddColumnMSA(Int4 t, Int4 pos, cma_typ L);
Int4    RmColumnMSA(Int4 t, Int4 lemon, cma_typ L);
BooLean InsertColCMSA(Int4 t, BooLean right, cma_typ cma);
cma_typ InsertColumnsCMSA(cma_typ cma, Int4 Blk, Int4 start_ins, Int4 length);
cma_typ InsertColumnsCMA(Int4 start, Int4 N, cma_typ cma);
cma_typ RmOverHangsCMSA(cma_typ cma);

//------------------ columns to insertions --------------------
Int4    IronOutOperation(char *operation);
char    *AddInsertToOperationArray(Int4 start, Int4 end, char *operation); // goes with ColumnsToInsertCMSA()
BooLean ColumnsToInsertCMSA(cma_typ cma,Int4 start_ins, Int4 end_ins);  // single block only?
cma_typ ColumnsToInsertCMSA2(cma_typ cma,Int4 start_ins, Int4 end_ins);
cma_typ ConvertColsToInsertsCMSA(cma_typ cma, Int4 Blk, Int4 start_ins, Int4 end_ins); // newer routine.

//--------------------- block operations -----------------------
cma_typ RmBlkCMSA(Int4 t, cma_typ L);
BooLean TrimCMSA(Int4 blk,unsigned short RmLeft,unsigned short RmRight,cma_typ cmsa);
cma_typ TrimCMSA(float info_cut, Int4 *TrimLimit, Int4 *RmLeft, Int4 *RmRight, 
	cma_typ cmsa);

//--------------------- fuse two blocks -----------------------
st_type FuseElementsSites3(Int4 x, Int4 maxlen, st_type S,a_type A); // older; delete?
cma_typ SimpleFuseBlksCMSA(Int4 x, Int4 maxlen, cma_typ L); // older; delete?
cma_typ FuseBlksCMSA(Int4 x, Int4 maxlen, cma_typ L); // older; pushes blocks together...
// cma_typ  FuseBlocksCMSA(Int4 Blk, cma_typ cma); // newer, more elegant routine..

//---------------------cma_typ conversions -----------------------
cma_typ RemoveOverhangsCMSA(cma_typ in_cma, BooLean AddX2Ends=FALSE);
// cma_typ OneBlockCMSA(cma_typ cma);	// convert from multi- to single block.
cma_typ RmWrinklesCMSA(cma_typ cma);	// remove deletions next to insertions.
void    ExtendFakeToRealCMSA(cma_typ cma);
void    ExtendFakeToRealCMSA(Int4 sq,cma_typ cma);
// void    IronOutCMSA(cma_typ cma);

//================== randomize ===============
cma_typ ShuffleSeqCMSA(double fraction, cma_typ cma);
cma_typ SimSeqToCMSA(e_type *EList, Int4 N, a_type A);
cma_typ AddRandomCMSA(cma_typ cma, Int4 num_random);
cma_typ RandomCMSA(Int4 ntyps, Int4 *sitelen, gss_typ& gss);
void	PutRandomCMA(FILE *fp, double *freq, Int4 Length, Int4 NumSeqs, a_type AB);
e_type *SimulatedSeqsCMSA(cma_typ cma, Int4 nSimSeq, Int4 rpts, Int4 *Gap_Len);
void    ShuffleColumnsCMA(cma_typ cma);

//********************* Additional OPERATIONS ******************************* 

#if 0	// moved back to twkcma: 1/26/08 (afn)
Int4    PutClusterOfCMSA(char *name, Int4 percent_ident,Int4 min_size,
        BooLean IncludeFirst,cma_typ cma);
#endif
double  JunLiuHMM_PenaltyCMA(FILE *fp, cma_typ cma);

//============= sampling routines (cmsa_sample.c) =====================

//------------- sampling columns --------------------
BooLean SampAddColCMSA(Int4 t, cma_typ L);
BooLean SampAddColTempCMSA(Int4 t, double temperature, cma_typ L);
BooLean SampRmColCMSA(Int4 t, cma_typ L);
BooLean SampRmColTempCMSA(Int4 t, double temperature, cma_typ L);
BooLean SampSlideColRtCMSA(Int4 t, double *oldMap, double *newMap,
	cma_typ L, double temperature);
BooLean SampSlideColLtCMSA(Int4 t, double *oldMap, double *newMap,
	cma_typ L, double temperature);
BooLean MoveMultiColsMSA(cma_typ G, Int4 t, Int4 num);
BooLean MoveColumnCMSA(cma_typ G, Int4 t);
BooLean TransferColumnCMSA(cma_typ G);
BooLean ShiftCMSA(cma_typ G, Int4 t);
Int4    SampNewColMSA(Int4 t, cma_typ L);
BooLean dfsSlideColRtCMSA(Int4 t, double *oldMap, double *newMap,
        cma_typ L, double temperature);
BooLean dfsSlideColLtCMSA(Int4 t, double *oldMap, double *newMap,
        cma_typ L, double temperature);
BooLean SampleEndColCMSA(Int4 t, double temperature, cma_typ cma);

//----------------- sampling blocks ---------------------
cma_typ DeleteBlkCMSA(cma_typ L);
cma_typ DeleteBlkCMSA(Int4 b,cma_typ L);
cma_typ AddBlkCMSA(Int4 x, Int4 lenx, cma_typ L);
cma_typ SplitBlkCMSA(Int4 , Int4 , cma_typ );
cma_typ SplitBlkCMSA(Int4 x, Int4 left_leng, Int4 minlen, cma_typ L);

//----------------- sampling gaps ---------------------
BooLean SampleGapsGibbsCMSA(Int4 s,cma_typ *cma, double Temperature);
BooLean ClusterSampleGapsGibbsCMSA(Int4 *sqset,cma_typ *oldcma,double Temperature); // not used?

//------------------- sampling sequence alignments ---------------
double  ReAlignBestCMSA(cma_typ cmsa);
char    *GapAlignSeqCMSA(FILE *ftpr,Int4 a, Int4 b, Int4 *Score, e_type  E,
        				Int4 **gapscore,cma_typ cmsa);
cma_typ	RefineCMSA(cma_typ msa);
void    RefineCMSA2(char method, cma_typ msa);
double  GappedSimAnnealCMSA(cma_typ *M,double StartTemp,double EndTemp,double inc);

//======== recombination routines cmsa_recombine.c =======================
Int4    CrossPtsCMSA(cma_typ ma1, cma_typ ma2,Int4 *aln1,Int4 *aln2);
cma_typ RecombineCMSA(cma_typ ma1, cma_typ ma2);
cma_typ MkRecombinantCMSA(cma_typ ma1, cma_typ ma2, Int4 *config);
cma_typ IntersectionCMSA(cma_typ ma1, cma_typ ma2);
Int4    SeqCrossPtsCMSA(Int4 n,st_type S1,st_type S2,Int4 *aln1,Int4 *aln2);

//============== gapped recombination routines =======================
cma_typ GRecombineCMSA(cma_typ ma1, cma_typ ma2);
double  GRecombinantMapCMSA(cma_typ ma1, cma_typ ma2, Int4 *config);
cma_typ RecombineGapsCMSA(cma_typ cma1, cma_typ cma2, Int4 *config);

//=========  LOG-PROBABILITY RATIO ROUTINES: cmsa_map.c ============
double  DirichletRelMap(cma_typ L);
double  RelMapCMSA(cma_typ cma);
double  UnGappedRelMapCMSA(cma_typ cma);
double  RelMapCMSA2(cma_typ cma);
double  RelMapMinusCMSA(Int4 tdel, cma_typ L);
void	SaveBestCMSA(cma_typ G);
Int4    *CntsFieldCMSA(cma_typ L, Int4 t, Int4 *len);
double  FieldRelMapCMSA(cma_typ L, Int4 t);
void    SetPseudoToMapCMSA(cma_typ cmsa);
double  PutRelMapCMSA(FILE *fp, cma_typ L);
double  RelMapMinusSeqCMSA(Int4 sq, cma_typ cma);

//================ OUTPUT ROUTINES: cmsa_put.c ====================
void    PrintCMSA(FILE *fptr, cma_typ G);
void    PutAlnCMSA(FILE *fptr, cma_typ cma);
void    PutAlnCMSA(char *name, cma_typ cma,gd_type G,Int4 Rpts=1);
void    PutSeqAlnCMSA(FILE *fptr, cma_typ L);
void    PutSitesCMSA(FILE *fptr,Int4 t,double **site_prob, cma_typ L);
void    PrintAlnCMSA(FILE *fptr, char *name, cma_typ L, gd_type G);
void    WriteMtfCMSA(char *name, cma_typ msa, gd_type G);
void    PutModelsCMSA(FILE *fp, cma_typ cma);
void    WritePhylipCMSA(char *name, cma_typ L);
void    WriteConsensusSeqCMSA(char *name, cma_typ L);
void    PutInterBlockCMSA(FILE *fptr, Int4 t1, Int4 t2, cma_typ L);
void    PutStuffedSeqSetCMSA(FILE *fptr, Int4 mingap, cma_typ L);
void    PutGappedAlnCMSA(FILE *fp,cma_typ cma,gd_type G,Int4 Rpts=1);
void    PutMSA(FILE *fptr, cma_typ G);
void    PutAsDomTypCMSA(FILE *fp,char *name,cma_typ cma);

void    PutRptSpacingCMSA(FILE *fp, cma_typ cma);
void    PutRptSpacingsCMSA(FILE *fp, cma_typ cma);
void    PutRptPartionCMSA(FILE *fp, cma_typ cma);

//=============== fullSeq information for SubSeq repeats. =================
void	AddFullCountsCMSA(ss_type FullSeq, unsigned short  *FullRpts, cma_typ cma);
void    CopyFullCountsCMSA(cma_typ cma2, cma_typ cma);
void	RmFullCountsCMSA(cma_typ cma);

//--------------------- repeats --------------------------
Int4	FullRptsCMSA(Int4 sq, cma_typ cma);
unsigned short *FullRptsCMSA(cma_typ cma);
Int4    *FirstRptFullCMSA(cma_typ cma);
Int4    RptNumCMSA(Int4 sq, cma_typ cma);

//=======================  rich text format routines: cmsa_rtf.cc =====================
void    CMA2RTF(FILE *fptr, double cbp_cut, double purge_cut,
	double infoLO, double infoHI, ss_type key_seq, char PageSetUp,
	cma_typ cma,BooLean LabeledOnly);
BooLean *purge_sma_cma(double cutoff,sma_typ MA,cma_typ cma);
void    CMA2RTF(FILE *fptr,double cbp_cut,double purge_cut,
	double infoLO,double infoHI,ss_type key_seq,char PageSetUp,
	double ***freq,sma_typ MA,cma_typ cma,BooLean LowerOnly=FALSE);
BooLean *ConservedCMA(Int4 *observed,double cutoff,double freq,cma_typ cma);
void    PutSchematicCMSA(FILE *fp, char FillColor, unsigned char ShapeType,
        char PageSetUp, BooLean group, cma_typ cma, ss_type keydata=0);

void	ReNameCMSA(char *newname, cma_typ cma);

//===================================================================
BooLean *FindCloseCMA(double cutoff,e_type E,Int4 *start,cma_typ cma);
BooLean *RepSetCMA(double cutoff,cma_typ cma);
e_type  **GetRepSetCMSA(FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma);
char    *NewConservedCMA(Int4 *observed,double cutoff,double *freq,a_type A);

void	SetBPPS_CMA(double alpha,Int4 A0, Int4 B0, char set_mode, cma_typ cma);
double	GetBPPS_CMA(Int4 *A0, Int4 *B0, char *set_mode, cma_typ cma);

//======================= junk ==============================
double DirichletRelMap(char method, cma_typ L);
double DirichletRelMap(cma_typ L);

/****************************** MACROS ******************************/
#define GetPosSitesCMSA(sq,c)   GetPosSites((sq),SitesCMSA(c))
#define	LevelCMSA(L)		((L)->Level)
#define	SetLevelCMSA(x,L)	((L)->Level=(x))
#define	FullSeqCMSA(L)		((L)->FullSeq)
#define	DomainsCMSA(L)		((L)->Domains)
#define	gssCMSA(L)		SitesGSS((L)->sites)
#define	gsqCMSA(s,L)		(SitesGSS((L)->sites)->GetGSQ((s)))
#define	gsqForceCMSA(s,L)	(SitesGSS((L)->sites)->ForcedGetGSQ((s)))
#define	mdlCMSA(L)		((L)->mdl)
#define	SetPenaltyCMSA(o,x,L)	(gssCMSA(L)->SetIndelPenalty((o),(x)))
#define	SetPerNatsCMSA(pn,L)	(gssCMSA(L)->SetPerNats((pn)))
#define	PerNatsCMSA(L)		(gssCMSA(L)->PerNats())
#define nullCMSA(L)             ((L)->null)
#define DataCMSA(L)		SitesSeqSet((L)->sites)
#define TrueDataCMSA(L)		SitesTrueSqSet((L)->sites)
#define SitesCMSA(L)		((L)->sites)
#define AlphabetCMSA(L)		SeqSetA(TrueDataCMSA(L))
#define NameCMSA(L)		NameSeqSet(TrueDataCMSA(L))
#define NameSeqSetCMSA(L)	NameSeqSet(DataCMSA((L)))
#define CountsCMSA(L)		CountsSites((L)->sites)
#define nBlksCMSA(L)		nTypeSites(SitesCMSA((L)))
#define NumSeqsCMSA(L)		NSeqsSeqSet(TrueDataCMSA((L)))
#define tFreqCMSA(L)		tFreqSeqSet(TrueDataCMSA((L)))
#define LenSeqCMSA(n,L)		SqLenSeqSet((n),DataCMSA((L)))
#define XSeqPtrCMSA(n,L)	XSeqSeqSet((n),DataCMSA((L)))
#define SeqPtrCMSA(n,L)		SeqSeqSet((n),DataCMSA((L)))
#define ModelsCMSA(L)           ((L)->mdl->Models())
#define	ModelCMSA(t,L)		((L)->mdl->Model((t)))
#define	LengthCMSA(t,L)		SiteLen((t),(L)->sites)
#define	ColResCntsCMSA(t,i,L)	GetSiteFreq((L)->sites,(t),(i)-1)
#define	LengthsCMSA(L)		SiteLengths((L)->sites)
#define	ContigBlkCMSA(t,L)	ContigFModel(ModelCMSA((t),(L)))
#define MaxSeqCMSA(L)		MaxSeqSeqSet(TrueDataCMSA((L)))
#define TrueSeqCMSA(n,L)	SeqSetE(n,TrueDataCMSA(L))
#define FakeSeqCMSA(n,L)	SeqSetE(n,DataCMSA(L))
#define PutSqDefLineCMSA(f,n,L)	PutSeqInfo2((f),SeqSetE((n),TrueDataCMSA(L)))
#define MaxTrueSeqCMSA(L)	MaxSeqSeqSet(TrueDataCMSA((L)))
#define MinTrueSeqCMSA(L)	MinSeqSeqSet(TrueDataCMSA((L)))
#define MinSeqCMSA(L)		MinSeqSeqSet(DataCMSA((L)))
#define PutSeqSetCMSA(fp,L)	PutSeqSetEs(fp,DataCMSA(L))
#define PutIDsCMSA(fp,L)	PutSeqSetPIDs(fp,DataCMSA(L))
#define PSEUDO_CMSA		0.1
#define AlignedCMSA(L)		((L)->best != NULL)
#define	RenameCMSA(name,L)	RenameSeqSet(name,TrueDataCMSA(L))

#endif
