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

#if !defined(_CLUSTER_)
#define _CLUSTER_

#include "seqset.h"
#include "afnio.h"
#include "gblast.h"
#include "gpsi_typ.h"
#include "swaln.h"
#include "karlin.h"
#include "residues.h"
#include "dsets.h"
#include "pairaln.h"

/******************************** private **********************************/

/******************************** PUBLIC **********************************/
e_type  **ClusterSeqs(e_type *ListE, char method, double cutoff, Int4 *counts,
        a_type A, Int4 *Nset,Int4 T);
e_type  **ClusterSeqs(e_type *ListE, char method, double cutoff, Int4 *counts,
        a_type A, Int4 *Nset);
int     PutRepSeq(FILE *fptr, char c, BooLean pval, e_type *EList, 
		BooLean UseLabeled, a_type A);

Int4    RepSetCluster(FILE *fp,char *DBS_NAME, Int4 sizedbs, BooLean xnu,
        char mode, Int4 cutoff, BooLean UseLabeled, a_type  A, Int4 T);
Int4    RepSetCluster(FILE *fp,char *DBS_NAME, Int4 sizedbs, BooLean xnu,
        char mode, Int4 cutoff, a_type  A);
Int4    RepSetCluster(FILE *fp,char *DBS_NAME, Int4 sizedbs, BooLean xnu,
        char mode, Int4 cutoff, BooLean UseLabeled, a_type A);
ds_type ClusterGPSI(char method, double cutoff, ss_type data);
e_type  **ClusterGPSI(char method,double cutoff,ss_type data,Int4 *Nset);
e_type  **ClusterGPSI(char method,double cutoff,ss_type data,a_type A,
	Int4 *Nset);
e_type  **ClusterGPSI(char method, double cutoff, ss_type data,a_type A,
	Int4 *Nset,Int4 Threshold);
Int4    *ClusterIntoSets(double cutoff, ss_type data, Int4 *Nset);

#endif
