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

#if !defined(_JACKKNIFE_)
#define _JACKKNIFE_

#include "cluster.h"
#include "seqset.h"
#include "sma.h"
#include "purge.h"
#include "oscan.h"
#include "goscan.h"

/****************************** Jackknife **********************************
  * use both PSI-BLAST and probe to sort out real from bogus hits.

 ***************************************************************************/

/******************************** private **********************************/

/******************************** PUBLIC **********************************/
Int4    Jackknife(ss_type data, ss_type msadata, sma_typ MA, sma_typ *MA2,
        ss_type *data2);
Int4    Jackknife(char method, double Cluster_cutoff, double cutoff, Int4 dbs_size,
		ss_type data, ss_type msadata, sma_typ MA, sma_typ *MA2, 
		ss_type *data2,Int4 MaxBadBlks);
e_type	*MergeSeqSet(Int4 *N, ss_type data1, ss_type data2);
BooLean *IntersectSeqs(e_type *ListE, ss_type data, Int4 *Nnew);
BooLean *DifferSeqs(e_type *ListE, ss_type data, Int4 *Nnew);

#endif

