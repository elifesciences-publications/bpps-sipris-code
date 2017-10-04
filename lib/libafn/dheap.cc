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

#include "dheap.h"

dh_type	dheap(Int4 N,Int4 D)
/* Initialize a heap to store items in {1,...,N}.*/
{
	dh_type	H; Int4 i;

	if(N >= INT4_MAX) print_error("dheap size over the limit; exiting");
	NEW(H,1,dheap_type);
	H->N = N; H->d = D; H->n = 0;
	NEW(H->h,N+1,Int4); NEW(H->pos,N+1,Int4); NEW(H->kvec,N+1,keytyp);
	for (i=1; i<= N; i++) H->pos[i] = 0;
	return H;
}

void	Nildheap(dh_type H)
{ free(H->h); free(H->pos); free(H->kvec); free(H); }

void	insrtHeap(Int4 i,keytyp k, dh_type H)
/* insert item i with specified key */
{
	assert(i > 0);
	assert(isfinite(k));
	// if(i<1) dheap_error("fatal! item to insert is < 1");
	if(H->pos[i]!=0) rmHeap(i,H);
	H->kvec[i] = k; H->n++; siftupHeap(i,H->n,H);
}

Int4	rmHeap(Int4 i,dh_type H)
/* Remove item i from heap. */
{
	Int4 j;
	if(H->pos[i]==0) return 0;
	j = H->h[H->n--];
	if (i != j && H->kvec[j] <= H->kvec[i]) siftupHeap(j,H->pos[i],H);
	else if (i != j && H->kvec[j]>H->kvec[i]) siftdownHeap(j,H->pos[i],H);
	H->pos[i]=0;
	return i;
}

Int4	delminHeap(dh_type H)
/* delete and return item with smallest key */
{
	Int4 i;
	if (H->n == 0) return 0;
	i = H->h[1];
	rmHeap(H->h[1],H);
	return i;
}


void	siftupHeap(Int4 i ,Int4 x,dh_type H)
/* Shift i up from position x to restore heap order.*/
{
	Int4 px = pHeap(x,H);
	while (x > 1 && H->kvec[H->h[px]] > H->kvec[i]) {
		H->h[x] = H->h[px]; H->pos[H->h[x]] = x;
		x = px; px = pHeap(x,H);
	}
	H->h[x] = i; H->pos[i] = x;
}

void	siftdownHeap(Int4 i,Int4 x,dh_type H)
/* Shift i down from position x to restore heap order.*/
{
	Int4 cx = minchildHeap(x,H);
	while (cx != 0 && H->kvec[H->h[cx]] < H->kvec[i]) {
		H->h[x] = H->h[cx]; H->pos[H->h[x]] = x;
		x = cx; cx = minchildHeap(x,H);
	}
	H->h[x] = i; H->pos[i] = x;
}

Int4	minchildHeap(Int4 x,dh_type H)
/* Return the position of the child of the item at position x
   having minimum key. */
{
	Int4 y, minc;
	if ((minc = leftHeap(x,H)) > H->n) return 0;
	for (y = minc + 1; y <= rightHeap(x,H) && y <= H->n; y++) {
		if (H->kvec[H->h[y]] < H->kvec[H->h[minc]]) minc = y;
	}
	return minc;
}

void	chkeyHeap(Int4 i,keytyp k, dh_type H)
/* Change the key of i and restore heap order.*/
{
	keytyp ki;
	if(H->pos[i]==0) return;
	ki = H->kvec[i]; H->kvec[i] = k;
	     if (k < ki) siftupHeap(i,H->pos[i],H);
	else if (k > ki) siftdownHeap(i,H->pos[i],H);
}

void	PutHeap(FILE *fptr,dh_type H)
/* Print the contents of the heap. */
{
	Int4 x;
	fprintf(fptr,"   h:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %2d",H->h[x]);
	fprintf(fptr,"\nkvec:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %3f",H->kvec[H->h[x]]);
	fprintf(fptr,"\n pos:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %2d",H->pos[H->h[x]]);
	fprintf(fptr,"\n");
}

void    dheap_error(char *s){fprintf(stderr,"dheap: %s\n",s);exit(1);}

