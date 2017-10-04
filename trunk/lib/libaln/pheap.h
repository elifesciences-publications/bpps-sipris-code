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

#if !defined (PHEAP)
#define	PHEAP
#include "afnio.h"
#include "histogram.h"
#include "mheap.h"
#include "probability.h"
#include "pattern.h"
#include "segment.h"
#include "stdinc.h"
#include "alphabet.h"
/*************************** ADT GHEAP ***************************
	
	pattern min-max heap

**********************************************************************/

/*************************** GHEAP type **************************/
typedef struct {
	/***** GENERATOR *****/
	mh_type		mheap;		/* min-maxheap for motifs */
	ptn_typ		*pattern;	/* pattern array */
	Int4		*C;		/* cardinality for G */
	s_type		**S;		/* lists of matching segments */
	Int4		d_max;		/* maximum depth */
	Int4		d_min;		/* minimum depth */
	Int4		k_max;		/* maximum length */
	a_type		A;		/* alphabet */
} pheap_type;
typedef pheap_type *ph_type;

/******************************* private *****************************/
BooLean	QuickPurgePheap(ptn_typ Q, Int4 C, ph_type H);

/******************************* Public ******************************/
/***************************** operations ****************************/
ph_type MkPheap(Int4 d_max, Int4 d_min, Int4 k_max, Int4 hpsz, a_type A);
ptn_typ  DelMinPheap(s_type **S, ph_type H);
ptn_typ  DelMaxPheap(s_type **S, ph_type H);
Int4	PurgePheap(ph_type H);
Int4	InsertPheap(ptn_typ G, ph_type H, s_type *S, Int4 C,double prob);
Int4	PutPheap(FILE *fptr,ph_type H);
ptn_typ *CombinePheap(ptn_typ **motif, s_type ***slist, Int4 *ngroup, 
        Int4 *npat, Int4 N, ph_type H);
Int4     Combine1PHeap(ptn_typ *GA, ptn_typ *SGA, s_type **SGL, Int4 ngrps, 
	Int4 N, s_type **List, Int4 TotCard, ph_type H);
Int4     CombineOSPHeap(ptn_typ *GA, ptn_typ *SGA, s_type **SGL, Int4 ngrps,
        Int4 N, s_type **List, Int4 TotCard);
ph_type	NilPheap(ph_type H);

/**************************** macro operations **********************/
#define ItemsInPheap(H)		ItemsInMheap((H)->mheap)
#define MinKeyPheap(H)		MinKeyMheap((H)->mheap)
#define MaxKeyPheap(H)		MaxKeyMheap((H)->mheap)
#define MinCardPheap(H)		((H)->C[MinItemMheap((H)->mheap)])
#define EmptyPheap(H)		EmptyMheap((H)->mheap)
#define	PheapA(H)		((H)->A)
#endif

