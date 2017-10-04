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

#if !defined (_CMA_GMB_)
#define _CMA_GMB_

#include "sma.h"
#include "cmsa.h"
#include "residues.h"
#include "gibbs.h"
#include "prtn_model.h"
#include "gss_typ.h"
#include "seqset.h"
#include "gsq_typ.h"
#include "histogram.h"
#include "random.h"
#include "gmb_typ.h"
#include "gsm_typ.h"
#include "dheap.h"

//====================== tweakcma.cc ====================
// Int4    IronOutOperation(char *operation);
// char    *AddInsertToOperationArray2(Int4 start, Int4 end, char *operation);

//====================== gsq_typ_main.cc ====================
void    PutShrinkFlanksCMSA(FILE *fp, Int4 left, Int4 right, cma_typ cma);
set_typ IndelBreakPntsCMSA(Int4 Blk, Int4 MinHalfBlk,cma_typ cma,Int4 cut=0);
void	PutIndelsCMSA(FILE *fp, cma_typ cma);
void	PutIndelBrkPtsCMSA(FILE *rfp,cma_typ cma);
// cma_typ RmOverHangsCMSA(cma_typ cma); // moved to cmsa_operations.cc

cma_typ SplitUpBlkCMSA(Int4 Blk, set_typ Set,cma_typ cma);
cma_typ ExtendBlkCMSA(cma_typ cma, Int4 Blk, Int4 AddLeft, Int4 AddRight);
cma_typ TrimDownCMSA(cma_typ cma, char mode='R', Int4 limit=4);
cma_typ ShiftBlkCMSA(Int4 Blk,Int4 shift,cma_typ cma);
cma_typ RmBlockCMSA(cma_typ cma, Int4 Blk);
cma_typ TrimBlkCMSA(cma_typ cma, Int4 Blk, Int4 RmLeft, Int4 RmRight, Int4 limit);
cma_typ SplitUpBlocksCMSA(set_typ Set,cma_typ cma);
cma_typ TransferColumnCMSA(Int4 Blk,BooLean right,Int4 numCols,cma_typ cma);

cma_typ OneBlockCMSA(cma_typ cma);
// cma_typ RmWrinklesCMSA(cma_typ cma);
void    ExtendFakeToRealCMSA(cma_typ cma);
void    ExtendFakeToRealCMSA(Int4 sq, cma_typ cma);
cma_typ FuseBlocksCMSA(Int4 Blk, cma_typ cma);
void    IronOutCMSA(cma_typ cma);

// gsq_typ *gsq_typ::RmOverHangs(Int4 nblk,Int4 *sites, Int4 *len,Int4 start,Int4 end);
// char    *gsq_typ::RmFlankingInOperation(char *operation);

cma_typ ContractColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj,double minFrq2Del=0.33);
cma_typ ExtendColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj,double minFrq2Del=0.33);
cma_typ ModifyColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj,
					double frq_cutoff=0.5,char mode='A');
cma_typ InsertColumnsCMSA(cma_typ cma, Int4 Blk, Int4 start_ins, Int4 length);
cma_typ MoveColumnsCMSA(cma_typ cma, char *operation);
cma_typ RemoveTheseColumnsCMSA(char *Operation, cma_typ CMA);
Int4    MoveAlignmentCMSA(set_typ Used, cma_typ From, cma_typ To);

void    PutSubGroupFootPrints(FILE *fp, set_typ Set, cma_typ cma);
void    PutSubGroupCsqScores(FILE *fp, set_typ Set, cma_typ cma);

//======================== gsm_plus.cc ===============================
cma_typ run_gismo_plus(int argc, char *argv[],a_type AB);
cma_typ run_gambit(FILE *rfp,char *name, Int4 iter, ssx_typ *isst,Int4 StageStart=0, Int4 Stages=4);
cma_typ *GetCoreCMSA(cma_typ cma1, cma_typ cma2,double &ave_frq,double freq_cut=0.34,FILE *fp=0);
cma_typ MvColumnsCMSA(ssx_typ *ssx, set_typ &SetSq);
// cma_typ hieraln_gismo_plus(int argc, char *argv[],a_type AB, ss_type data=0);
//======================== gsm_plus.cc ===============================

//======================== gmb_smpl.cc ===============================
set_typ *RtnTargetSizeSetsCMSA(Int4 &numSets,Int4 &percent_ident,cma_typ CMA,double MaxFraction=0.10);
set_typ *RtnFastClustersCMSA(Int4 &numSets, Int4 percent_ident,set_typ InSet,cma_typ cma);
Int4    NumAlnSeqsCMSA(cma_typ cma);
// void    PutAlnSeqsCMSA(FILE *fp,cma_typ cma);
cma_typ CreateRandomCMA(Int4 nblk, Int4 *len,ss_type data,double per_nats);

//========================== gmb_subaln.cc ========================
cma_typ SampleSubAlignsCMSA(char *name,cma_typ cma,char dms_mode=' ');

//======================== hieraln.cc ===============================
// set_typ RtnFastRepSetCMSA(FILE *fp_err, Int4 percent_ident,set_typ InSet,cma_typ cma);
set_typ RtnRepSubSetCMSA(FILE *fp_err, Int4 percent_ident,set_typ Set,cma_typ cma);
cma_typ InSetMkCMSA(set_typ set, cma_typ cma);
void    ExtendAlnToEndsCMSA(cma_typ cma);
set_typ AutoPurgeCMSA(Int4 s, set_typ SubGrp,char *name,cma_typ cma);
double  ResidueDiversityCMSA3(set_typ Set, cma_typ cma);
cma_typ MkSubCMSA(set_typ Set, BooLean  AddCsq, cma_typ cma);
e_type  MKCsqSubCMSA(set_typ Set,cma_typ cma);
Int4    PutWithCsqCMSA(FILE *fp,e_type CsqE, cma_typ cma);
void    PutCsqScoresCMSA(FILE *fp,set_typ set,cma_typ cma);
cma_typ	AddThisCsqCMSA(e_type CsqE, cma_typ cma);

char    *FixTmplOperation(char *operation, e_type CsqE, Int4 &Start, Int4 &TrimNt,
                Int4 &TrimCt);
double	ComputSeqWtsCMSA(cma_typ CMA,set_typ Set=0);

//========================== will send to cmsa.cc ========================
Int4	NumInsertsCMSA(Int4 blk, Int4 site, cma_typ cma, Int4 &ins_res);
void    ColumnQualityCMSA(ssx_typ *ssx);
double  **BILD_ScoresCMSA(ssx_typ *ssx);
cma_typ RmWorstColsCMSA(Int4 mincol,ssx_typ *ssx);
Int4    RmGappyColumnsCMSA(double cutoff, cma_typ &cma);
Int4    NumberDeletionsCMSA(Int4 Column, cma_typ cma);
Int4    ModeLengthSeqSet(ss_type data);
set_typ *IdenticallyAlignedCMSA(cma_typ cmaA, cma_typ cmaB);
Int4    *CommonColsCMSA(cma_typ cmaA, cma_typ cmaB, double *&DD,double *&Del,double frq_cut);
#if 0
	{
	   Int4	*ColI2J,i,lenA=LengthCMSA(1,cmaA); 
	   set_typ *setA=CommonColsCMSA(cmaA,cmaB, DD,Del,frq_cut,ColI2J);
	   for(i=1; i <= lenA; i++){ NilSet(setA[i]);} free(setA);
	}
set_typ	*CommonColsCMSA(cma_typ cmaA, cma_typ cmaB, double *&DD,double *&Del,double frq_cut, Int4 *&ColI2J);
#endif
set_typ **ParticipatingColPairsCMSA(cma_typ cmaA, cma_typ cmaB, double frq_cut);
cma_typ SortByCMA(cma_typ cma, cma_typ scma);

//========================== junk routines ======================
// Int4    PutNewCMSA(FILE *fp,set_typ Set, BooLean put_csq, cma_typ cma);
cma_typ LengthenBlksCMSA(cma_typ cma);

// cma_typ gmb_typ::Optimize();

#endif

