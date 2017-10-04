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

#include "oscan2msa.h"

cma_typ ScanHeap2CMSA(snh_typ sH, Int4 cutoff, char *name, a_type A, 
	BooLean fragment)
/** destroy scan heap and return cmsa; cutoff = purge cutoff score **/
/** NULL is returned if none of the sequences qualify **/
{
	BooLean *use,**null;
	e_type	*E,*E2;
	sni_typ	*IList;
	Int4	N,ntyps,*lens,i,j,s,t,site;
	st_type	S;
	ss_type	data;
	cma_typ	msa;
	BooLean swiss,*rm;
	
	/** get alignment information from scanheap **/
	ntyps = nblksScanHeap(sH);
	NEW(lens,ntyps +2,Int4);
	NEWP(null,ntyps +2,BooLean);
	for(i=1; i <= ntyps; i++) lens[i]=BlkLenScanHeap(i,sH);
	IList = RtnNilScanHeap(sH); 

	/** get sequences and sites from scanheap **/
	N=0; while(IList[N+1] != NULL) N++; 
	if(N < 1) {
		for(s=1; s <= N; s++) NilScanInfo(IList[s]);
		free(IList); return NULL;
	}
	NEW(E,N+2,e_type); NEW(E2,N+2,e_type);
	for(s=1; s <= N; s++){ E[s] = SeqScanInfo(IList[s]); } 

	if(cutoff < 0){ 
	        NEW(use,N+3, BooLean);
		for(s=1; s <= N; s++) use[s] = TRUE;
	} else if(cutoff == 0){
	  NEW(use,N+3, BooLean);
          NEW(rm,N+2,BooLean);
          for(i=1; i<= N; i++){
            if(!rm[i]){
             swiss = SwissSeq(E[i]);
             if(swiss) use[i] = TRUE;
             for(j=i+1; j<= N; j++){
                if(!rm[j]){
                  if(IdentSeqs(E[i],E[j])){ 
                    if(!swiss) { 
                        if((swiss=SwissSeq(E[j]))) use[j] =TRUE;
                    }
                    rm[j]=TRUE;
                  }
                }
             }
             if(!swiss) use[i]=TRUE;
             rm[i]=TRUE;
            }
          }
	  free(rm);
#if 0
	     NEW(use,N+3, BooLean);
	     for(s=1; s<= N; s++){
		e1 = E[s]; use[s] = TRUE;
                for(j=s+1; j<= N; j++){
                    if(IdentSeqs(e1,E[j])){ use[s]=FALSE; break; }
		}
             }
#endif
	} else use = RmHomologsList(cutoff, 'B', 1, N, FALSE, 11, N, E, A);
	for(i=0, s=1; s <= N; s++) if(use[s]) { i++; E2[i] = E[s]; }
	free(E); 
	if(i < 2) {
		for(s=1; s <= N; s++) NilScanInfo(IList[s]);
		free(IList); return NULL;
	}

	/** create sites **/
	data = Array2SeqSet(E2, i, name, A); /** E2 converted to SeqSet **/
	S = MkSites(ntyps, lens, data,10,5,5.0,9,9);
	for(i=0, s=1; s <= N; s++){
	    if(use[s]) { 
		i++; 
		for(t=1; t <= ntyps; t++){
		     site = SiteScanInfo(t,IList[s]);
		     AddSite(t, i, site, S);
		}
	    	NilScanInfoRtnE(IList[s]); /** this will not delete sequence. **/
	    } else NilScanInfo(IList[s]); /** this will delete sequence. **/
	}
	free(IList); free(lens); free(use);
	if(fragment) print_error("ScanHeap2CMSA( ) fragment flag not implemented");

	/** create cmsa **/
	msa = MakeCMSA(S, null);
	if(fragment) for(i=1; i <= ntyps; i++){ free(null[i]); } 
	free(null);
	return msa;
}

