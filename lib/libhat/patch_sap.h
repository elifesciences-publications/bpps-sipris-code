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

#if !defined (_PATCH_SAP_)
#define _PATCH_SAP_
#include "my_ncbi.h"
#include "alphabet.h"
#include "wdigraph.h"
#include "sph_typ.h"

#if 0   //************************ Problem SAPs ****************************

1. Repeats in Subject sequence.

2. Repeats in Query sequence.

3. Non-colinear regions of homology (e.g., circularly permuted).

4. Multiple artifacts (combinations of the above.

//**************************** Ideal situation: ****************************
--> Non-repetitive colinear SAPs

    ========================================================= Query
     qs1          qe1
     ---------------
     ss1          se1

#endif  //******************************************************************

Int4    *LongestPathGSAP(Int4 num_saps, sap_typ *SAP,UInt4 *QS,
                UInt4 *QE,UInt4 *SS,UInt4 *SE);
Int4    RawInfoGSAP(Int4 *QS, Int4 *QE, Int4 *SS, Int4 *SE, sap_typ sap);
Int4    *RunningScoreSAP_L(sap_typ sap,Int4 **mx,a_type A);
Int4    *RunningScoreSAP_R(sap_typ sap,Int4 **mx,a_type A);
Int4    QueryToSubjectSiteSAP(sap_typ sap,Int4 Qsite);
BooLean TrimOverlappingSAPs(sap_typ tail_sap, sap_typ head_sap,Int4 crosspoint);
Int4    *EvolvingRmOverlapsGSAP(Int4 overlap_cutoff, sap_typ *SAP,Int4 *path,
                                        Int4 num,Int4 **mx,a_type A);
char    *OperationFromGSAP(sap_typ sap,UInt4 *QS, UInt4 *SS,
                UInt4 *QE,UInt4 *SE);
sap_typ PatchSAPs(Int4 overlap_cutoff, Int4 **mtx,sap_typ raw_sap,Int4 num_seqs,
                        e_type queryE, a_type AB);
sap_typ PatchSAPs(Int4 overlap_cutoff, Int4 **mtx,sap_typ raw_sap,e_type queryE, a_type AB);

double  FractionOverlapSAP(sap_typ sap1, sap_typ sap2);
BooLean AtLeastOneLabeledSAPS(sap_typ sap);
void    RemoveRejectsSAP(sap_typ *HEAD, Int4 num_sap,a_type AB);
Int4    LabelBestSAP(sap_typ *sap,UInt4 start,UInt4 N,a_type AB);
Int4    FindBestOverlappingSAPs(sap_typ *head, Int4 num_sap_lists,Int4 num_seqs,a_type AB);
Int4    FindBestOverlappingSAPs(sap_typ *head, Int4 num_sap_lists,a_type AB);
Int4    LabelBestSAPs(sap_typ *HEAD, Int4 num_sap, Int4 num_seqs, BooLean get_all_hits,
		a_type AB, UInt4 *NUM_EACH_CMA);

#define MAX_ALLOWED_GHSPs       1000


#endif

