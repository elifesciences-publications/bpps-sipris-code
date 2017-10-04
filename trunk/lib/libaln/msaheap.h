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
#if !defined (MSAHEAP)
#define MSAHEAP
#include "cmsa.h"
#include "histogram.h"
#include "residues.h"
#include "mheap.h"
#include "random.h"

/********************* CMSA Heap Type ***********************

/*********************************************************************/

/*************************** msaheap type **************************/
typedef struct {
	mh_type		mH;		/** minmax heap for msas **/
	cma_typ 	*msa;		/** msas**/
	Int4		size;		/** heap size **/
} msaheap_type;

typedef msaheap_type *mah_typ;
/*********************************************************************/

/******************************* private *******************************/

/******************************* PUBLIC *******************************/
mah_typ	MkMSAHeap(Int4 hpsz);
void    NilMSAHeap(mah_typ H);
Int4    InsertMSAHeap(cma_typ M, double key, mah_typ H);
cma_typ RmMSAHeap(Int4 item, mah_typ H);
cma_typ SeeMSAHeap(Int4 item, mah_typ H);
Int4	RandMSAHeap(mah_typ H);
cma_typ	DelMinMSAHeap(double *map, mah_typ H);
cma_typ	DelMaxMSAHeap(double *map, mah_typ H);
cma_typ	DelBestMSAHeap(double *map, mah_typ H);
double  PurgeMSAHeap(mah_typ H);
double	ConvergedMSAHeap(mah_typ H);
BooLean KeyInMSAHeap(double key, mah_typ H);

/*********************************************************************/
#define nMSAHeap(H)		(ItemsInMheap((H)->mH))
#define	FullMSAHeap(H)		FullMheap((H)->mH)
#define EmptyMSAHeap(H)		EmptyMheap((H)->mH)
#define	keyMSAHeap(i,H)		(-keyMheap((i),(H)->mH))
#define BestItemMSAheap(H)	MinItemMheap((H)->mH) 
#define SizeMSAHeap(H)		SizeMheap((H)->mH) 



#endif

