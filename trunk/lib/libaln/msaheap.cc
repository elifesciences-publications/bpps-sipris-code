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

#include "msaheap.h"

mah_typ MkMSAHeap(Int4 hpsz)
{
	mah_typ	H;

	NEW(H,1,msaheap_type);
	NEW(H->msa,hpsz+2,cma_typ);
	H->mH=Mheap(hpsz, 3);
	H->size = hpsz;
	return H;
}

double	ConvergedMSAHeap(mah_typ H)
{
	Int4	i;
	double	v;

	for(i=1; i<=SizeMheap(H->mH); i++){
	    if(MemMheap(i,H->mH)){
		assert(H->msa[i] != 0);
		fprintf(stderr,"key=%g ",-keyMheap(i,H->mH));
		fprintf(stderr,"(%d blocks; %d columns)\n",
			nBlksCMSA(H->msa[i]), NumColumnsCMSA(H->msa[i]));
	    }
	}
	if(!FullMheap(H->mH)) { v = DBL_MAX; }
	else {
	    v = MinKeyMheap(H->mH) - MaxKeyMheap(H->mH);
	    v = fabs(v);
	    // fprintf(stderr,"difference = %g\n",v);
	} return v;
}

BooLean	KeyInMSAHeap(double key, mah_typ H)
{
	Int4	i;
	double	k;

	for(i=1; i<=SizeMheap(H->mH); i++){
            if(MemMheap(i,H->mH)){
		k = keyMSAHeap(i,H); if(k==key) return TRUE;
	    }
        }
	return FALSE;
}

double	PurgeMSAHeap(mah_typ H)
/*** remove all but one of those items with identical keys ***/
{
	Int4	n,i,j,N=SizeMheap(H->mH);
	double	*map,mp;
	cma_typ ma,*msa;
	gss_typ *gssX;

	if(nMSAHeap(H) < 2) return nMSAHeap(H);
	NEW(msa, N+2, cma_typ); NEW(map, N+2, double);
	for(n=0; (ma=DelMinMSAHeap(&mp, H)) != NULL; ){
		n++; map[n] = mp; msa[n] = ma;
	}
	for(i=1; i < n; i++){
	  if(msa[i] != NULL){
	    for(j=i+1; j <= n; j++){
		if(msa[j] != NULL && map[i] == map[j]){
			gssX=gssCMSA(msa[j]); 
			if(gssX) gssX->~gss_typ();
			NilCMSA(msa[j]); msa[j]=NULL;
		}
	    }
	  }
	}
	for(i=1; i <= n; i++){
	   if(msa[i] != NULL) {
		if(InsertMSAHeap(msa[i], map[i], H) == NULL){
			gssX=gssCMSA(msa[i]); 
			if(gssX) gssX->~gss_typ();
			NilCMSA(msa[i]); msa[i]=NULL;
		}
	   }
	}
	free(msa); free(map);
	return nMSAHeap(H);
}

void    NilMSAHeap(mah_typ H)
{
	cma_typ	M;
	double	key;
	gss_typ *gssX;

	while((M=DelMinMSAHeap(&key,H)) != NULL){
		gssX=gssCMSA(M); if(gssX) gssX->~gss_typ(); NilCMSA(M);
	} NilMheap(H->mH); free(H->msa); free(H);
}

Int4    InsertMSAHeap(cma_typ M, double map, mah_typ H)
{
        Int4    item;
	gss_typ *gssX;

        item=InsertMheap(-map, H->mH);
        if(item==NULL) return NULL;
        if(H->msa[item]!= NULL){
	    gssX=gssCMSA(H->msa[item]); 
	    if(gssX) gssX->~gss_typ(); 
	    NilCMSA(H->msa[item]);
	} H->msa[item]=M; 
        return item;
}

cma_typ RmMSAHeap(Int4 item, mah_typ H)
{
	cma_typ M=NULL;

	if(RmMheap(item, H->mH) != NULL){
		M = H->msa[item]; H->msa[item] = NULL;
	} 
	return M;
}

cma_typ SeeMSAHeap(Int4 item, mah_typ H) { return H->msa[item]; }

Int4    RandMSAHeap(mah_typ H)
{
	Int4	item;

	do {
		item = random_integer(H->size) + 1;
	} while(H->msa[item]== NULL);
	return item;
}

cma_typ DelBestMSAHeap(double *map, mah_typ H){ return DelMinMSAHeap(map,H); }

cma_typ DelMinMSAHeap(double *map, mah_typ H)
{
        Int4    item;
	cma_typ M;

        *map= -MinKeyMheap(H->mH);
        item = DelMinMheap(H->mH);
	M = H->msa[item]; H->msa[item]=NULL;
        return M; 
}

cma_typ DelMaxMSAHeap(double *map, mah_typ H)
{
        Int4    item;
	cma_typ M;

        *map= -MaxKeyMheap(H->mH);
        item = DelMaxMheap(H->mH);
	M = H->msa[item]; H->msa[item]=NULL;
        return M; 
}

