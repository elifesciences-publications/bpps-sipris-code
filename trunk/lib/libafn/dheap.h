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

#if !defined(DHEAP)
#define DHEAP
/********* Header file for d-heap data structure. *********
 Maintains a subset of items in {1,...,m}, where item has a key.
 *********************************************************/
#include "stdinc.h"

typedef double keytyp;

typedef	struct {
	Int4	N;			/* max number of items in heap */
	Int4	n;			/* number of items in heap */
	Int4	d;			/* base of heap */
	Int4	*h;			/* {h[1],...,h[n]} is set of items */
	Int4	*pos;			/* pos[i] gives position of i in h */
	keytyp	*kvec;			/* kvec[i] is key of item i */
} dheap_type;

typedef dheap_type *dh_type;

/************************ private operations ***************************/
/* parent of item, leftmost and rightmost children */
#define pHeap(x,H)          (((x)+((H)->d-2))/(H)->d)
#define leftHeap(x,H)       ((H)->d*((x)-1)+2)
#define rightHeap(x,H)      ((H)->d*(x)+1)
#define	MAX_KEY		    DBL_MAX

Int4	minchildHeap(Int4 i,dh_type H);	/* returm smallest child of item */
void	siftupHeap(Int4 i ,Int4 x,dh_type H);
				/* move item up to restore heap order */
void	siftdownHeap(Int4 i,Int4 x,dh_type H);	
				/* move item down to restore heap order */
void    dheap_error(char *s);

/************************ public operations ***************************/
dh_type	dheap(Int4 N,Int4 D);
void	Nildheap(dh_type H);
void	insrtHeap(Int4 i,keytyp k, dh_type H);	
					/* insert item with specified key */
Int4	rmHeap(Int4 i,dh_type H);	/* remove item from heap */
Int4	delminHeap(dh_type H);		/* delete and return smallest item */
void	chkeyHeap(Int4 i,keytyp k, dh_type H);		
					/* change the key of an item */
void	PutHeap(FILE *fptr, dh_type H);		/* print the heap */

/************************* macros definitions **************************/
#define	minHeap(H)	((H)->n==0?NULL:(H)->h[1])
#define	ItemsInHeap(H)	((H)->n)
#define minkeyHeap(H)	((H)->kvec[(H)->h[1]])
#define minItemHeap(H)	((H)->h[1])
#define keyHeap(i,H)	((H)->kvec[(i)])
#define memHeap(i,H)	((H)->pos[(i)]!=NULL?TRUE:FALSE)
#define	emptyHeap(H)	((H)->n==0)
#define	fullHeap(H)	((H)->n >= (H)->N)
#endif



