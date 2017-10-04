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

/*   my_ncbi.h
* ===========================================================================
* File Name:  ncbimath.h
* Author:  Gish, Kans, Ostell, Schuler
* Version Creation Date:   10/23/91
* $Revision: 6.0 $
* File Description:
*   	prototypes for portable math library
* Modifications:
* $Log: ncbimath.h,v $
* Revision 6.0  1997/08/25 18:16:37  madden
* Revision changed to 6.0 */

#ifndef _MY_NCBI_
#define _MY_NCBI_

#include "stdinc.h"
#include "afnio.h"
#include "sequence.h"
#include "dsets.h"
#include "dheap.h"
#include "histogram.h"

#ifndef Char
typedef char            Nlm_Char, *Nlm_CharPtr;
#define Char            Nlm_Char
#define CharPtr         Nlm_CharPtr
#endif

#ifndef Boolean
typedef unsigned char   Nlm_Boolean, *Nlm_BoolPtr;
#define Boolean         Nlm_Boolean
#define BoolPtr         Nlm_BoolPtr
#endif

#ifndef Uint1
typedef unsigned char   Nlm_Uint1, *Nlm_Uint1Ptr;
#define Uint1           Nlm_Uint1
#define Uint1Ptr        Nlm_Uint1Ptr
#define UINT1_MAX       UCHAR_MAX
#endif

#ifndef Int2
typedef short           Nlm_Int2, *Nlm_Int2Ptr;
#define Int2            Nlm_Int2
#define Int2Ptr         Nlm_Int2Ptr
#define INT2_MIN        SHRT_MIN
#define INT2_MAX        SHRT_MAX
#endif

typedef int		Nlm_Int4, *Nlm_Int4Ptr;

#ifndef Int4
#define Int4            Nlm_Int4
#endif

#ifndef Int4Ptr
#define Int4Ptr         Nlm_Int4Ptr
#endif


/*----------------------------------------------------------------------*/
/*      Misc Common Macros                                              */
/*----------------------------------------------------------------------*/

#ifndef TRUE
#define TRUE	1
#endif

#ifndef FALSE
#define FALSE	0
#endif

#ifndef MIN
#define MIN(a,b)        ((a)>(b)?(b):(a))
#endif

#ifndef MAX
#define MAX(a,b)        ((a)>=(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a)  ((a)>=0?(a):-(a))
#endif

#ifndef SIGN
#define SIGN(a) ((a)>0?1:((a)<0?-1:0))
#endif

#ifndef DIM
#define DIM(A) (sizeof(A)/sizeof((A)[0]))
#endif

void    *GMemNew(size_t size);
void	*GMemFree(void *ptr);

typedef struct gdenseg {
    Int2		dim, numseg;
    UInt4	subject_id; 	// AFN addition
    UInt4	query_id; 	// AFN addition
    e_type		sE;		// AFN addition (sequence pointer)
    e_type		qE;		// AFN addition (sequence pointer)
    Int4		*starts;        // dimension is dim * numseg.
    Int4		*lens;          // dimension is numseg.
} GDenseSeg, *dsp_typ;

dsp_typ GDenseSegNew (void);
void    GDenseSegPut(FILE *fp, dsp_typ gdsp);
dsp_typ GDenseSegFree(dsp_typ dsp);
unsigned char *GetSequenceWithGDenseSeg(dsp_typ dsp,Boolean query,Int4Ptr start, 
		Int4Ptr length);
BooLean DeleteGDenseSeg(BooLean first, dsp_typ dsp);
Int4    ComputeQueryScoreDSP(FILE *fp,e_type qE,Int4 **mx);
Int4    ComputeScoreDSP(FILE *fp,dsp_typ dsp,Int4 **mx);
Int4    ComputeScoreDSP(dsp_typ dsp,Int4 **mx);
#if 1	// Routines for HMM-BLAST
Int4    MaxSegScoreDSP(dsp_typ dsp, Int4 **mtx, e_type sE, Int4 *sp_nmbr);
Int4    MaxWordScoreDSP(dsp_typ dsp,Int4 sp_nmbr,Int4 word_len,Int4 **mtx,
		Int4 mtx_len,e_type sE);
#endif

typedef struct gseqalign {
    double	evalue,bit_score; 	// AFN addition
    Int4	score;	 		// AFN addition
    char	label;			// AFN addition
    Int2	dim;
    dsp_typ segs;   
    struct gseqalign *next;
} sap_type, *sap_typ; // = AFN SeqAlign, *SeqAlignPtr;

sap_typ MakeGSeqAlign(Int2 numseg, UInt4 query_id, 
        UInt4 subject_id, e_type qE, e_type sE, Int4 *starts,
        Int4 *lens);
BooLean NumberOfGapsGSAP(h_type HG, sap_typ sap);
Int4    DriverNumberOfGapsGSAP(FILE *fp,sap_typ head);
BooLean StartEndScoreGSAP(sap_typ sap, Int4 *Start, Int4 *End, Int4 **mx,
	a_type A);
BooLean DeleteOverlapGSAP(sap_typ sap1,sap_typ sap2, Int4 **mx,a_type A);
// sap_typ RmOverlapsGSAP(sap_typ head,Int4 **mx,a_type A);
// void    RemoveOverlapsGSAP(sap_typ head,Int4 **mx,a_type A);
BooLean FixAlnRunOverGSAP0(FILE *fp, sap_typ sap, Int4 **mx, Int4 minsubscore,
        a_type A);
Int4    FixAlnRunOverGSAP(sap_typ head,Int4 **mx,a_type A,Int4 minfix);
BooLean FixSplitsGSeqAlign(sap_typ sap); // AFN FIX;
BooLean IsSplitGSeqAlign(register sap_typ sap);
void    MkGlobalSeqAlign(sap_typ sap,a_type AB);
void    MakeGlobalSeqAlign(sap_typ head,a_type AB);

void    GSeqAlignPut(FILE *fp, sap_typ gsap);
double  PercentIdentGSAP(Int4 *net_len, sap_typ sap);
Int4	IdentitiesGSAP(Int4 *alnlen, Int4 *Ident, Int4 *Inserts,sap_typ sap);
sap_typ SortBySeqIDGSAP(sap_typ sap);
BooLean IsQueryGSAP(sap_typ sap);
Int4	PutSeqsListGSAP(FILE *fp, sap_typ head, a_type AB);
void    PutSameSeqsListGSAP(FILE *fp, sap_typ head, Int4 MinLen, double fract_IDs, a_type AB);
void    PutGSeqAlign(FILE *fp, sap_typ gsap, Int4 width, a_type A);
void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB);
void    PutOneSeqAlign(FILE *fp, sap_typ sap, Int4 width, a_type AB);
sap_typ MinEvalSAP(sap_typ head, double *MinEval);
void    PutMultiGSeqAlign(FILE *fp, sap_typ sap, Int4 width, a_type AB);
Int4    *PutScwrlSeqGSAP(FILE *fp,sap_typ sap, e_type keyE, a_type A);

Int4    *PutSCGenSeqGSAP(FILE *fp,sap_typ sap, e_type keyE, a_type A);
Int4	DelimitGSeqAlignList(sap_typ head,Int4 **SS,Int4 **SP,Int4 **ES,Int4 **EP);
void    DelimitGSeqAlign(sap_typ sap, Int4 *SS,Int4 *SP,Int4 *ES,Int4 *EP);
Int4    QueryStartGSeqAlnList(sap_typ sap);
sap_typ GSeqAlignNew(void);
sap_typ GSeqAlignFree(sap_typ anp);
void    FreeGSeqAlignList(sap_typ sap);
sap_typ ToGSeqAlign(Int4 numopers, char *operation, e_type qE, e_type sE,
        Int4 start1, Int4 start2);
void    PutRasMolGSAP(FILE *fp, sap_typ sap, Int4 N_Colors, const char **Colors);
void    PutSubSeqsListGSAP(FILE *fp, sap_typ head, Int4 left, Int4 right,
	a_type AB);
void    PutGSubSeq(FILE *fp, sap_typ sap, Int4 left, Int4 right,
	a_type AB);
Int4    InfoGSAP(Int4 *QS, Int4 *QE, Int4 *SS, Int4 *SE, sap_typ sap);
// Routines to convert SeqAlign to CMA output file...
void    SeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank, Int4 rightflank,
        a_type AB);
void    PutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank,
	a_type AB);
Int4    NumSeqsListGSAP(sap_typ head);
Int4    NumHSPsListGSAP(sap_typ head);
Int4    LengthListGSAP(sap_typ head);
BooLean OverlapGSAP(sap_typ sap1, sap_typ sap2);
void    PutMergedGSAP(FILE *fp, sap_typ head, Int4 left, Int4 right,
        a_type AB);
sap_typ AlignToGSAP(e_type qE,e_type sE,char *operation,Int4 start_qE,Int4 start_sE);
sap_typ ConcatenateGSAP(sap_typ head, sap_typ sap);
sap_typ PathHMMToGSAP(e_type hmmE, e_type sE, char *hmm_path, double *Evalue,
        double *bit_score,Int4 score,Int4 *scores);

Int4    GetStartEndGSAP(sap_typ sap, UInt4 *Start, UInt4 *End);
void    PutDeleteBestHSPs(FILE *fp, sap_typ head, a_type AB);
void    GetSeqAlignTable(unsigned short *TabP, unsigned char *TabR, sap_typ sap,
        a_type AB);
void    SeqAlignToTable(FILE *fp,sap_typ head, a_type AB);

sap_typ TrimSAP(sap_typ sap, Int4 trimleft, Int4 trimright);
sap_typ MakeSelfAlignedGSAP(e_type E);

// DEFINES:
#define EvalueGSAP(sap)		((sap)->evalue)
#define ScoreGSAP(sap)		((sap)->score)
#define SetLabelGSAP(c,sap)	((sap)->label=(c))
#define LabelGSAP(sap)		((sap)->label)
#define NextGSAP(sap)		((sap)->next)
#define QuerySeqGSAP(sap)       ((sap)->segs->qE)
#define SubjectSeqGSAP(sap)     ((sap)->segs->sE)


#endif /* !_MY_NCBI_ */

