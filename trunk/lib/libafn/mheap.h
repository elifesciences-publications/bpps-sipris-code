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

#if !defined (MMHEAP)
#define	MMHEAP
#include "dheap.h"
#include "stdinc.h"
/*************************** ADT MMHEAP ***************************
	
	min-max heap

**********************************************************************/

/*************************** MMHEAP type **************************/
typedef struct {
	dh_type		heap;		/* minheap */
	dh_type		maxheap;	/* maxheap */
	Int4		hpsz;		/* heap size */
	Int4		nfree;		/* #items available */
	Int4		*avail;		/* list of available item no */
} mheap_type;

typedef mheap_type *mh_type;

/******************************* private *****************************/

/******************************* Public ******************************/
/***************************** operations ****************************/
mh_type Mheap(Int4 hpsz, Int4 d);
Int4	DelMinMheap(mh_type H);
Int4	DelMaxMheap(mh_type H);
Int4	RmMheap(Int4 i, mh_type H);
Int4	InsertMheap(keytyp key, mh_type H);
mh_type	NilMheap(mh_type H);

/**************************** macro operations **********************/
#define ItemsInMheap(H)		ItemsInHeap((H)->heap)
#define MinKeyMheap(H)		minkeyHeap((H)->heap)
#define MinItemMheap(H)		minItemHeap((H)->heap)
#define MaxKeyMheap(H)		(-minkeyHeap((H)->maxheap))
#define EmptyMheap(H)		emptyHeap((H)->heap)
#define SizeMheap(H)		((H)->hpsz)
#define FullMheap(H)		fullHeap((H)->heap)
#define	keyMheap(i,H)		keyHeap((i),(H)->heap)
#define MemMheap(i,H)		memHeap(i,(H)->heap)
#endif

