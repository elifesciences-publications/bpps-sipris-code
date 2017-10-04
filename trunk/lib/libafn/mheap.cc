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

#include "mheap.h"

mh_type	Mheap(Int4 hpsz, Int4 d)
{
	mh_type H;
	Int4	i;

	NEW(H,1,mheap_type); 
	NEW(H->avail,hpsz+1,Int4);
	H->heap=dheap(hpsz,d);
	H->maxheap=dheap(hpsz,d);
	H->hpsz = hpsz;
	H->nfree = hpsz;
	for(i=1; i<=hpsz; i++) H->avail[i] = i;
	return H;
}

mh_type	NilMheap(mh_type H)
{
	if(H==NULL) return (mh_type) NULL;
	Nildheap(H->heap);
	Nildheap(H->maxheap);
	free(H->avail);
	free(H); 
	return (mh_type) NULL;
}

Int4	DelMinMheap(mh_type H)
{
	Int4 i;
	
	if((i=delminHeap(H->heap))!=0){
		rmHeap(i,H->maxheap); 
		H->nfree++;
		H->avail[H->nfree] = i;
		return i;
	} else return 0;
}

Int4	DelMaxMheap(mh_type H)
{
	Int4 i;
	
	if((i=delminHeap(H->maxheap))!=0){
		rmHeap(i,H->heap); 
		H->nfree++;
		H->avail[H->nfree] = i;
		return i;
	} else return 0;
}

Int4	RmMheap(Int4 i, mh_type H)
{
	if(rmHeap(i,H->heap) == 0) return 0; 
	rmHeap(i,H->maxheap); 
	H->nfree++;
	H->avail[H->nfree] = i;
	return i;
}

Int4	InsertMheap(keytyp key, mh_type H)
/* NULL == not inserted; -1 == inserted but none deleted  */
{	
	Int4 i;

	if(H->nfree > 0){
		i = H->avail[H->nfree];
		H->nfree--;
	} else if(minkeyHeap(H->maxheap) < -key) {
		i=delminHeap(H->maxheap);
		rmHeap(i,H->heap);
	} else return 0;
	insrtHeap(i,key,H->heap);
	insrtHeap(i,-key,H->maxheap);
	return i;
}

