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

#if !defined (_GPSI2CMA_)
#define _GPSI2CMA_

#include "gpsi_typ.h"
#include "alphabet.h"
#include "residues.h"
#include "smatrix.h"
#include "cmsa.h"

/********************************* private ********************************/

/********************************* PUBLIC ********************************/
Int4    NumSeqsListGSAP(sap_typ head, Int4 max_sq_id, Int4 *min_qs, Int4 *max_qe);
Int4    NumSeqsListGSAP(sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
		Int4 *min_qs, Int4 *max_qe);
Int4    NumHSPsListGSAP(sap_typ head, Int4 max_sq_id, Int4 *min_qs, Int4 *max_qe);
Int4    NumHSPsListGSAP(sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
		Int4 *min_qs, Int4 *max_qe);
Int4    InDelsSeqAlign(sap_typ sap, Int4 *ins, Int4 *del,Int4 min_qs, Int4 max_qe);
void    tPutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank,
        	Int4 min_qs, Int4 max_qe, a_type AB);
void    xPutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank,
        	Int4 min_qs, Int4 max_qe, Int4 SeqID, a_type AB);

void    tSeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank, Int4 rightflank,
		a_type AB);
void    tSeqAlignToCMA(FILE *fp, sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
	Int4 leftflank, Int4 rightflank, a_type AB);
void    xSeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank, Int4 rightflank,
        a_type AB);
void    xSeqAlignToCMA(FILE *fp, sap_typ head, Int4 max_sq_id, Int4 leftflank,
	Int4 rightflank, a_type AB);
void    xSeqAlignToCMA(FILE *fp, sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
	Int4 leftflank, Int4 rightflank, a_type AB);
void    QueryFirstSeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank,
        Int4 rightflank, a_type AB);
void    QueryFirstSeqAlignToCMA(FILE *fp, sap_typ head, Int4 max_sq_id,
	Int4 leftflank, Int4 rightflank, a_type AB);
sap_typ RtnSAP_PsiBLAST(e_type qE, ss_type data, int argc, char *argv[],
        const char *USAGE, Int4 maxrounds, a_type A,Int4 minfix,BooLean use_gibbs,
			cma_typ mtfcma,BooLean use_checkin);
cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
        const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
        Int4 minfix,BooLean UseAllHSPs,Int4 MaxSqID,cma_typ mtfcma,
        BooLean use_checkin);
cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
		const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
		Int4 minfix,BooLean UseAllHSPs,Int4 MaxSqID,cma_typ mtfcma);
cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
        const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
        Int4 minfix,BooLean UseAllHSPs,cma_typ mtfcma);

cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
		const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
		Int4 minfix,BooLean UseAllHSPs);
#if 0	// These aren't used...
cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
        const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
        Int4 minfix,BooLean UseAllHSPs,Int4 MaxSqID);
cma_typ PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[],
        Int4 maxrounds, a_type A);
#endif


/********************************* MACROS ********************************/

#endif

