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

/* goscan.h - codes and constants for ordered scan program. */
#if !defined (_EDITCMA_)
#define _EDITCMA_
#include "cmsa.h"
#include "brh_typ.h" // ncbi-type blast routines.
#include "dsets.h"
#include "wdigraph.h"
#include "smooth.h"
#include "goscan.h"
#include "cluster.h"
#include "swaln.h"
#include "gpsi_typ.h"
#include "ema_typ.h"
#include "HMM_typ.h"
#include "hma_typ.h"

/********************************* private ********************************/
ema_typ CMAtoEMA(cma_typ cma);
cma_typ EMAtoCMA(ema_typ ema, char *pssm_arg, cma_typ cma);

/********************************* PUBLIC ********************************/
cma_typ AddRelatedCMSA(double cluster_cut, double add_cut, char *pssm_arg,
	cma_typ cma);

char    SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature, 
	UInt4 min_rpt,cma_typ *cma);
char    SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature,
        double *OldMap, UInt4 min_rpt,cma_typ *cma);
char    SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature,
        double prob_accept, double *OldMap, UInt4 min_rpt,
	cma_typ *cma,psm_typ *psm);
char    SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature,
        double prob_accept, double *OldMap, UInt4 min_rpt,
	hma_typ *hma);

char    SampleGSeqGibbsCMSA(Int4 s,BooLean *skipseq, cma_typ *oldcma,
        double *OldMap, double Temperature, double prob_accept,
	UInt4 min_rpt,psm_typ *psm);

/********************************* MACROS ********************************/

#endif

