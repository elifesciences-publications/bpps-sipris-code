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

/* msaheap.h - sequence min/max heap. */
#if !defined (_HSPHEAP_)
#define _HSPHEAP_
#include "gpxdrop.h"
#include "my_blastkar.h"
#include "histogram.h"
#include "residues.h"
#include "mheap.h"
#include "random.h"

// The next two structures are the final output produced by GBLAST.
// Formatters should then convert the data into GSeqAligns or the GBLAST ASN.1 spec.
typedef struct _gblast_results_hsp {
        Int4    number;         // number of HSP's used to calculate the p-value.
        Int4    score;          // score of this HSP.
        double  e_value,        // expect value of this set of HSP's.
                p_value,        // p-value of this set of HSP's.
                bit_score;      // above score * lambda/ln2
        Int4    query_start,   // Start of the query HSP.
                query_length,   // Length of the query HSP.
                query_end,	// end of the query HSP.
                subject_start, // Start of the subject HSP.
                subject_length, // Length of the subject HSP.
                subject_end,	// end of the subject HSP.
                hspset_cnt;     // which set of HSP's?
        Int4    query_gapped_start,   // Starting points (on original HSP)
                subject_gapped_start; // for a gapped extension with X dropoff.
        gxeb_typ gap_info;      // ALL gapped alignment is here
} hsp_type, *hsp_typ; // GBLASTResultHsp, *GBLASTResultHspPtr;

hsp_typ MakeHSPTyp(gab_typ gap_align, sbp_typ sbp);
hsp_typ MakeHSPTyp(gab_typ gap_align, double e_value,double bit_score);
BooLean WithinHSPTyp(hsp_typ hsp1, hsp_typ hsp2);
double  OverlapHSPTyp(hsp_typ hsp1,hsp_typ hsp2);
void    PutHSPTyp(FILE *fp, hsp_typ hsp);
void    NilHSPTyp(hsp_typ hsp);

/*********************************************************************/

/*************************** msaheap type **************************/
typedef struct {
	mh_type		mH;		// minmax heap for hsps
	hsp_typ		*hsp_array;		
	Int4		size;		/** heap size **/
} hspheap_type;

typedef hspheap_type *hhp_typ;
/*********************************************************************/

/******************************* private *******************************/
BooLean WithinHSP(Int4 query_start, Int4 subject_start,
        Int4 query_end, Int4 subject_end, hsp_typ hsp);
/******************************* PUBLIC *******************************/
hhp_typ	MkHSPHeap(Int4 hpsz);
void    NilHSPHeap(hhp_typ H);
Int4    InsertHSPHeap(hsp_typ M, hhp_typ H);
hsp_typ SeeHSPHeap(Int4 item, hhp_typ H);
BooLean WithinHSPHeap(Int4 qs, Int4 ss, Int4 qe, Int4 se, hhp_typ H);
Int4	PurgeHSPHeap(hhp_typ H);
hsp_typ	DelMinHSPHeap(double *key, hhp_typ H);
hsp_typ	DelMaxHSPHeap(double *key, hhp_typ H);

/*********************************************************************/
#define nHSPHeap(H)		(ItemsInMheap((H)->mH))
#define	FullHSPHeap(H)		FullMheap((H)->mH)
#define EmptyHSPHeap(H)		EmptyMheap((H)->mH)
#define	keyHSPHeap(i,H)		(-keyMheap((i),(H)->mH))
// #define MinItemMSAheap(H)	MinItemMheap((H)->mH) 
#define SizeHSPHeap(H)		SizeMheap((H)->mH) 

#endif

