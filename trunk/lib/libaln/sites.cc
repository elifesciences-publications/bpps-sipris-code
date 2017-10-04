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

#include "sites.h"

st_type	MakeSites(Int4 ntyps, Int4 *len_elem, gss_typ &gss0)
// makes a copy of gapped sequence seq constructor 
{
    st_type	S;

    if(ntyps > MAX_NO_TYPESITES){
         // assert(!(ntyps > MAX_NO_TYPESITES));
	 sites_error("Too many types.");
    }
    S = new sites_type[1]; // replaces NEW(S,1,sites_type);
    S->gss=gss0;	// copy constructor called.
    mk_sites(ntyps, len_elem, S);
    return S;
}

st_type CreateNewSites(Int4 ntyps, Int4 *len_elem, st_type S)
// Create a "copy" of S using the same gapped sequences. 
// WARNING: NO SITE ARE ADDED!
{ return MakeSites(ntyps,len_elem, S->gss); }

st_type	MkSites(Int4 ntyps, Int4 *len_elem, ss_type data,Int4 open,
	Int4 extend, double pernats, Int4 leftflank, Int4 rightflank)
// creates a gapped sequence set from ss_type sequence set.
{
    st_type	S;

    if(ntyps > MAX_NO_TYPESITES) sites_error("Too many types.");
    S = new sites_type[1]; // replaces NEW(S,1,sites_type);
    // S->gss is initialized by constructor function.
    S->gss.initialize(data,open,extend,pernats,leftflank,rightflank);
    // initialize to full dataset.
    mk_sites(ntyps, len_elem, S);
    return S;
}

void	mk_sites(Int4 ntyps, Int4 *len_elem, st_type S)
/*********************************************************************
 Initialize a site structure containing no sites. 
							 nsites: A B C
  seq[1]    = .................................................* 0 0 0
  seq[2]    = ...........................................*	 0 0 0
     :
  seq[nseq] = ...............................................*   0 0 0
**********************************************************************/
{
    Int4	n,m,t,max;
    ss_type	data=SitesSeqSet(S);

    S->ntyp = ntyps;
    S->nseq = NSeqsSeqSet(data);
    S->len_seq = LengthsSeqSet(data);
    NEWP(S->type, S->nseq+2, char);
    NEW(S->pos, S->nseq+2, ol_type);
    for(max=0,n = 1; n <= S->nseq; n++){
	max = MAXIMUM(Int4,max,S->len_seq[n]);
	NEW(S->type[n],(S->len_seq[n]+2),char);
	S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	S->pos[n] = Olist(S->len_seq[n]+1);
    }
    NEW(S->tmp,max+2,Int4);
    NEW(S->len_elem, ntyps+2, Int4);
    NEW(S->totsites, ntyps+2, Int4);
    NEWP(S->nsites,ntyps+2,Int4);
    NEWPP(S->site_pos, ntyps+2,unsigned short);
    for(t = 1; t <= ntyps; t++) {
	S->totsites[t] = 0;
        S->len_elem[t] = len_elem[t];
        NEW(S->nsites[t], S->nseq+2,Int4);
	NEWP(S->site_pos[t], NSeqsSeqSet(data)+2,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		S->nsites[t][n] = 0;
                m = ((Int4) SqLenSeqSet(n,data)/(Int4)len_elem[t]) +3;
                NEW(S->site_pos[t][n],m+2,unsigned short);
	}
	if(t==1) S->maxinc = S->len_elem[t] - 1;
	else S->maxinc = MINIMUM(Int4,S->maxinc,(S->len_elem[t]-1));
    }
    S->maxinc = MAXIMUM(Int4,1,S->maxinc);
}

st_type	StartSites(Int4 ntyps, Int4 *sitelen, gss_typ& gss)
// get starting sites from argument strings.
{
   st_type 	S;
   FILE	   	*fptr;
   Int4    	l,arg,s,t,i,n,len;
   dh_type	H;
	
   if(ntyps > MAX_NO_TYPESITES) sites_error("too many element types");
   for(t=1; t<=ntyps; t++){
	if(sitelen[t] < 2 || sitelen[t] > MAX_LENG_SITES){
   		for(t=1; t<=ntyps; t++){
		  fprintf(stderr,"sitelen[%d] = %d\n",t,sitelen[t]);
		} print_error("element length out of range\n"); 
	}
   }
   S=MakeSites(ntyps, sitelen, gss); 
   ss_type	data = SitesSeqSet(S);
   n = MaxSeqSeqSet(data);
   H=dheap(n+2,3);
   for(n = 1; n <= NSeqsSeqSet(data); n++) {
	    len = S->len_seq[n];
	    for(t = 1; t <= ntyps; t++){
		len -= S->len_elem[t];
		insrtHeap(t,(keytyp)Random(),H);
	    }
            for(s=t; len > 0 ; s++,len--) insrtHeap(s,(keytyp)Random(),H);
            for(t=s=1; (i=delminHeap(H)) != NULL; ){
		// WARNING: need to generalize beyond one colinear! 
	    	if(i<=S->ntyp){ AddSite(t,n,s,S); s+=S->len_elem[t]; t++; }
	    	else s++;
	    }
    }
    Nildheap(H);
    return S;
}

BooLean ColinearSites(st_type S)
{
    Int4 t,n,s;

    /** 1. Check for colinearity. **/
    if(S->ntyp == 0) return TRUE;
    for(n=1; n<=S->nseq; n++) if(S->nsites[1][n]!=1) return FALSE;
    for(t=2; t<=S->ntyp; t++) {
          for(n=1; n<=S->nseq; n++) {
            if(S->nsites[t][n] != 1) return FALSE;
            if(S->site_pos[t-1][n][1] > S->site_pos[t][n][1]) return FALSE;
          }
    }
    return TRUE;
}

st_type	FuseElementsSites(Int4 x, Int4 maxlen, st_type S)
/*****************************************************************
 if element x and x+1 are <= maxlen then fuse them into one 
 element; create and return the resultant object S2.
 *****************************************************************/
{
    Int4        t,n,s,N,T,T2,*len_elem,len,len2,newlen;
    st_type	S2;

    // 0. Test for input errors
    T = nTypeSites(S); T2=T-1; 
    if(x < 1 || x >= T) return NULL;
    len = SiteLen(x,S); len2 = SiteLen(x+1,S);
    if(len > maxlen || len2 > maxlen) return NULL;

    // 1. Create new sites object 
    NEW(len_elem, T2+3, Int4);
    for(t=1; t <= T2; t++) {
        if(t < x) len_elem[t]=SiteLen(t,S);
        else if(t == x) { 
		len_elem[t]=len+len2; 
	} else len_elem[t]=SiteLen(t+1,S);
    }
    S2=MakeSites(T2,len_elem, S->gss);
    free(len_elem);

    // 2. Add sites to new object 
    N=NSeqsSites(S2);
    for(n=1; n<=N; n++){
	for(t=1; t <= T2; t++) {
          if(t < x) AddSite(t,n,SitePos(t,n,1,S),S2);
	  else if(t == x){
	    switch(random_integer(2)) {
		case 0: // slide 2nd element to the left
		  AddSite(t,n,SitePos(t,n,1,S),S2); break; 
		case 1: // slide 1st element to the right 
		  AddSite(t,n,SitePos(t+1,n,1,S)-len,S2); break;
		default: print_error("FuseElementsSites( ) - this should not happen!");
		break;
	    }
	  } else AddSite(t,n,SitePos(t+1,n,1,S),S2);
	}
    }
    return S2;
}

st_type	FuseElementsSites2(Int4 x, Int4 maxlen, st_type S)
/*****************************************************************
 NEW: 5/1/04 - Andy Neuwald
 if element x and x+1 are <= maxlen then fuse them into one 
 element; create and return the resultant object S2.
 WARNING: NOT TESTED.
 *****************************************************************/
{
    Int4        t,n,s,s2,N,T,T2,*len_elem,len,len2,newlen;
    st_type	S2;

    print_error("FuseElementsSites2 not yet implemented");
    // 0. Test for input errors
    T = nTypeSites(S); T2=T-1; 
    if(x < 1 || x >= T) return NULL;
    len = SiteLen(x,S); len2 = SiteLen(x+1,S);
    if(len > maxlen || len2 > maxlen) return NULL;

    // 1. Create new sites object 
    NEW(len_elem, T2+3, Int4);
    for(t=1; t <= T2; t++) {
        if(t < x) len_elem[t]=SiteLen(t,S);
        else if(t == x) len_elem[t]=len+len2; 
	else len_elem[t]=SiteLen(t+1,S);
    } S2=MakeSites(T2,len_elem, S->gss);
    free(len_elem);

#if 1
    // 2. Conver region between fused blocks into an insert
    N=NSeqsSites(S2);
    for(n=1; n<=N; n++){
      s=EndSitePos(x,n,1,S);
      s2=SitePos(x+1,n,1,S);
#if 0
      char *operation=gss->Operation(n);
#endif
      // if((s2-s) > 0) S2->gss.InsertGap(n,s,(s2-s));
      S->gss.Put(stdout,n);
    }
#endif

    // 3. Add sites to new object 
    N=NSeqsSites(S2);
    for(n=1; n<=N; n++){
	for(t=1; t <= T2; t++) {
          if(t < x) AddSite(t,n,SitePos(t,n,1,S),S2);
	  else if(t == x) AddSite(t,n,SitePos(t,n,1,S),S2); 
	  else AddSite(t,n,SitePos(t+1,n,1,S),S2);
	}
    } return S2;
}

void	InsertGapSites(Int4 n, Int4 pos, unsigned short gap, st_type S)
// Insert a gap at pos within sequence n; retain all site positions.
// WARNING: Vacates all sites in n prior to inserting a gap.
{
	Int4	t,m,max1,max2,len2;

	assert(n >= 1 && n <= S->nseq); 
        VacateSites(n,S);
	max1=MaxSeqSeqSet(SitesSeqSet(S));
	S->gss.InsertGap(n,pos,gap);
	len2=LenSeq(S->gss.FakeSeq(n));
	max2=MaxSeqSeqSet(SitesSeqSet(S));
	if(max2 > max1){ free(S->tmp); NEW(S->tmp,max2+2,Int4); }
	S->len_seq[n]=len2;
	free(S->type[n]); NEW(S->type[n],(S->len_seq[n]+2),char);
	NilOlist(S->pos[n]); 
	S->pos[n] = Olist(S->len_seq[n]+1);
        S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	for(t=1; t<=nTypeSites(S); t++){
            free(S->site_pos[t][n]);
            m = ((Int4) len2/(Int4)S->len_elem[t]) + 3;
            NEW(S->site_pos[t][n],m+2,unsigned short);
	}
}

void	ReplaceSeqSites(Int4 n, gsq_typ *gsq, st_type S)
// Replaces nth sequence in S with E; all sites in nth sequence are removed.  
// Returns replaced sequence if successful; otherwise NULL is returned.
{
	Int4	m,max1,max2,len1,len2;
	e_type	newE;

	assert(n >= 1 && n <= S->nseq); 
        VacateSites(n,S);
	len1=LenSeq(S->gss.FakeSeq(n)); 
	max1=MaxSeqSeqSet(SitesSeqSet(S));
#if 0	// old routine...
	newE=gsq->FakeSeq(); 
	len2=LenSeq(newE);
	S->gss.Replace(n,gsq);	// replaces and destroys old Sequence.
#endif
#if 1	// new routine to allow domain sampling.
	if(gsq == NULL){  // then remove gapped seq.
	   S->gss.RmFake(n);	// restores TrueSeq
	   newE=S->gss.TrueSeq(n); len2=LenSeq(newE);
	} else {
	   newE=gsq->FakeSeq(); len2=LenSeq(newE);
	   S->gss.Replace(n,gsq);	// replaces and destroys old Sequence.
	}
#endif
	max2=MaxSeqSeqSet(SitesSeqSet(S));
	if(max2 > max1){ free(S->tmp); NEW(S->tmp,max2+2,Int4); }
	S->len_seq[n] = LenSeq(newE);
	if(len1 < len2){
	  free(S->type[n]); NEW(S->type[n],(S->len_seq[n]+2),char);
	} 
	NilOlist(S->pos[n]); 
	S->pos[n] = Olist(S->len_seq[n]+1);
        S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	for(Int4 t=1; t<=nTypeSites(S); t++){
            free(S->site_pos[t][n]);
            m = ((Int4) LenSeq(newE)/(Int4)S->len_elem[t]) + 3;
            NEW(S->site_pos[t][n],m+2,unsigned short);
	}
}

gsq_typ *SwapSeqSites(Int4 n, gsq_typ *gsq, st_type S)
// Replaces nth sequence in S with E; all sites in nth sequence are removed.  
// Returns replaced sequence if successful; otherwise NULL is returned.
{
	Int4	m,max1,max2,len1,len2;
	e_type	newE;
	gsq_typ	*rtn_gsq=0;

	assert(n >= 1 && n <= S->nseq); 
        VacateSites(n,S);
	len1=LenSeq(S->gss.FakeSeq(n)); 
	max1=MaxSeqSeqSet(SitesSeqSet(S));
	assert(gsq != NULL);
	newE=gsq->FakeSeq(); len2=LenSeq(newE);
	rtn_gsq=S->gss.Swap(n,gsq);	// replaces and returns old Sequence.
	max2=MaxSeqSeqSet(SitesSeqSet(S));
	if(max2 > max1){ free(S->tmp); NEW(S->tmp,max2+2,Int4); }
	S->len_seq[n] = LenSeq(newE);
	if(len1 < len2){ free(S->type[n]); NEW(S->type[n],(S->len_seq[n]+2),char); } 
	NilOlist(S->pos[n]); 
	S->pos[n] = Olist(S->len_seq[n]+1);
        S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	for(Int4 t=1; t<=nTypeSites(S); t++){
            free(S->site_pos[t][n]);
            m = ((Int4) LenSeq(newE)/(Int4)S->len_elem[t]) + 3;
            NEW(S->site_pos[t][n],m+2,unsigned short);
	} return rtn_gsq;
}

st_type SplitElementSites(Int4 x, Int4 length1, Int4 minlen, st_type S)
/*****************************************************************
 if element x is >= minlen + length1 then split it into two blocks; 
 create and return the resultant object S2.
 *****************************************************************/
{
    Int4        t,n,s,N,T,T2,*len_elem,len,leftleng,rightleng;
    st_type     S2;
    double      r;

    // 0. Test for input errors
    T = nTypeSites(S); T2=T+1;
    if(x < 1 || x > T) return NULL;
    len = SiteLen(x,S) - length1;
    // if(len < minlen || len < 3 || length1 < 3) return NULL;
    // modified minimum length to accommodate curated_srch routines...
    if(len < minlen || len < 1 || length1 < 1) return NULL;
    leftleng = length1; rightleng = len;

    // 1. Create new sites object  (same as below; eventually merge...)
    NEW(len_elem, T2+3, Int4);
    for(t=1; t <= T2; t++) {
        if(t < x) len_elem[t]=SiteLen(t,S);
        else if(t == x) len_elem[t]=leftleng;
        else if(t == (x+1)) len_elem[t]=rightleng;
        else len_elem[t]=SiteLen(t-1,S);
    }
    S2=MakeSites(T2,len_elem, S->gss);
    N=NSeqsSites(S2);
    free(len_elem);

    // 2. Add sites to new object
    for(n=1; n<=N; n++){
        for(t=1; t <= T2; t++) {
          if(t <= x) AddSite(t,n,SitePos(t,n,1,S),S2);
          else if(t == (x+1)) AddSite(t,n, leftleng+SitePos(x,n,1,S),S2);
          else AddSite(t,n,SitePos(t-1,n,1,S),S2);
        }
    }
    return S2;
}

st_type	SplitElementSites(Int4 x, Int4 minlen, st_type S)
/*****************************************************************
 if element x is >= minlen then evenly split it into two elements; create
 and return the resultant object S2.
 *****************************************************************/
{
    Int4        t,n,s,N,T,T2,*len_elem,len,leftleng,rightleng;
    st_type	S2;
    double      r;

    // 0. Test for input errors
    T = nTypeSites(S); T2=T+1; 
    if(x < 1 || x > T) return NULL;
    len = SiteLen(x,S);
    if(len < minlen || len < 6) return NULL;
    leftleng = rightleng = len/2;
    if(len%2 == 1) if(random_integer(2)==0) leftleng++; else rightleng++; 

    // 1. Create new sites object 
    NEW(len_elem, T2+3, Int4);
    for(t=1; t <= T2; t++) {
        if(t < x) len_elem[t]=SiteLen(t,S);
        else if(t == x) len_elem[t]=leftleng;
        else if(t == (x+1)) len_elem[t]=rightleng;
        else len_elem[t]=SiteLen(t-1,S);
    }
    S2=MakeSites(T2,len_elem, S->gss);
    N=NSeqsSites(S2);
    free(len_elem);

    // 2. Add sites to new object 
    for(n=1; n<=N; n++){
	for(t=1; t <= T2; t++) {
          if(t <= x) AddSite(t,n,SitePos(t,n,1,S),S2);
          else if(t == (x+1)) AddSite(t,n, leftleng+SitePos(x,n,1,S),S2);
          else AddSite(t,n,SitePos(t-1,n,1,S),S2);
	}
    }
    return S2;
}

st_type	AddElementSites(Int4 lenx, Int4 x, Int4 min_free, 
	Int4 min_squeeze, st_type S)
/*****************************************************************
 If there is room add another element between the xth and x+1st element.
 requires that each sequence have at least min_free space after adding
 the element and that no more than min_squeeze sequences have to 
 shift adjacent sites to perform this operation.
 If the operation fails of if sites are not colinear then zero is returned.
 *****************************************************************/
{
    Int4        j,t2,t,n,s,newsite,end;
    Int4        N,T,T2,total,tight,lastsite,nextsite,space;
    Int4        *len_elem,*start,*nsites,*lastgap,*gap,*nextgap;
    st_type	S2;
    ss_type     data;
    double      r;

    // 1. Check for colinearity. 
    if(!ColinearSites(S)) return NULL;
    T = nTypeSites(S); T2=T+1;
    total=lenx + min_free;
    for(t=1; t<= T; t++){ total += SiteLen(t,S); }
    data=SitesSeqSet(S); N = NSeqsSeqSet(data);
    for(n=1; n<=N; n++) {
	if(SqLenSeqSet(n,data) < total) return NULL;
    }

    // 2. Check whether enough room for new sites. 
    NEW(start,N+3,Int4); NEW(nsites,N+3,Int4);
    NEW(lastgap,N+3,Int4); NEW(gap,N+3,Int4); NEW(nextgap,N+3,Int4); 
    if(x >= T) x=T;
    for(tight=0, n=1; n<=N; n++) {
	// gap determination 
        if(x < 1) { 
		start[n]=1; 
		gap[n]=SitePos(1,n,1,S)-1; 
	} else if(x < T){
          start[n] = SitePos(x,n,1,S) + SiteLen(x,S);
          gap[n] = SitePos(x+1,n,1,S) - start[n];
        } else {
          start[n] = SitePos(T,n,1,S) + SiteLen(T,S);
          gap[n] = SeqLenSites(n,S) - start[n] + 1;
        }
        if(gap[n] < lenx) { 
	   Int4 room; tight++;
	   if(T > 0){
	     // determine lastgap (between x-1 & x)
	     if(x > 1){
	       s = SitePos(x-1,n,1,S) + SiteLen(x-1,S);
	       lastgap[n] = SitePos(x,n,1,S) - s;
	     } else if(x == 1){
	       lastgap[n] = SitePos(1,n,1,S) - 1;
	     } // else lastgap[n] = 0;
	     // determine nextgap (between x+1 & x+2)
	     if(x == (T - 1)){
	       s = SitePos(x+1,n,1,S) + SiteLen(x+1,S);
	       nextgap[n] = (SeqLenSites(n,S) + 1) - s;
	     } else if(x < (T-1)){
	       s = SitePos(x+1,n,1,S) + SiteLen(x+1,S);
	       nextgap[n] = SitePos(x+2,n,1,S) - s;
	     } // else nextgap[n] = 0;
	   }
	   room = lastgap[n] + gap[n] + nextgap[n];
	   if(tight > min_squeeze || lenx > room){
		free(gap); free(lastgap); free(nextgap);
		free(start); free(nsites); return NULL;
	   } else { nsites[n]=gap[n]-lenx+1; }
        } else { nsites[n]=gap[n]-lenx+1; }
    }

    // 3. Create new sites object 
    NEW(len_elem, T+3, Int4);
    for(t=1; t <= T+1; t++) {
        if(t<=x) len_elem[t]=SiteLen(t,S);
        else if(t==(x+1)) len_elem[t]=lenx;
        else len_elem[t]=SiteLen(t-1,S);
    }
    S2=MakeSites(T2,len_elem, S->gss);
    free(len_elem); N=NSeqsSites(S2);

    // 4. First add sites to easy fitting sequences 
    for(n=1; n<=N; n++){
      if(nsites[n] > 0){	// = 'easy' fitting sequences
	for(t=1; t <= T2; t++) {
          if(t<=x) AddSite(t,n,SitePos(t,n,1,S),S2);
	  else if(t==(x+1)){
           do{
                  r = (double) Random()/(double) RANDOM_MAX;
                  newsite=(Int4)(r*(double)nsites[n]);
           } while(newsite>=nsites[n]);
           newsite+=start[n]; AddSite(t,n,newsite,S2);
	  } else AddSite(t,n,SitePos(t-1,n,1,S),S2);
	}
      }
    }

    // 5. if some sequeces are tight then squeeze these in too.
    if(tight > 0){
      // first determing addition room to left & right of new site
      for(n = 1; n <= N; n++) {
        if(nsites[n] <= 0){	// tight fitting sequences
          for(t=1; t < x; t++) AddSite(t,n,SitePos(t,n,1,S),S2);
	  for(t=x+3; t <= T2; t++) AddSite(t,n,SitePos(t-1,n,1,S),S2);
	  // need to move site x to the left and/or site x+2 to the right
	  if(x > 0) lastsite=SitePos(x,n,1,S); else lastsite=0;
	  newsite = start[n]; 
	  if(x < T) nextsite=SitePos(x+1,n,1,S); else nextsite=0;
	  for(Int4 g=gap[n]; g < lenx; ) {
	    switch(random_integer(2)) {
	      case 0: // move next site one space to the right 
		if(nextgap[n] > 0){ g++; nextsite++; nextgap[n]--; break; }
	      case 1: // move previous site one space to the left
		if(lastgap[n] > 0) { 
			g++; newsite--; 
			lastsite--; lastgap[n]--; break; 
		}
	      default:
		if(nextgap[n] > 0) { g++; nextsite++; nextgap[n]--; break; }
		else print_error("AddElementSites( ): this should not happen");
		break;
	    }
	  }
	  // add these sites if present
	  AddSite(x+1,n,newsite,S2);
          if(lastsite > 0) AddSite(x,n,lastsite,S2); 
          if(nextsite > 0) AddSite(x+2,n,nextsite,S2); 
        } // end tight fitting sequences
      }
    }

    free(start); free(nsites); 
    free(gap); free(lastgap); free(nextgap);
    return S2;
}

BooLean	AddRandomCLSites(Int4 n, st_type S)
// Add random colinear sites to sequence n in S.
// 
{
    Int4	i,j,t,s,len;
    dh_type	H;
	
    if(n < 1 || n > S->nseq) return FALSE;
    VacateSites(n,S);
    i = MaxSeqSeqSet(SitesSeqSet(S));
    H=dheap(i+2,3);
    len = S->len_seq[n];
    for(t=1; t <= S->ntyp; t++){
	    len -= S->len_elem[t]; insrtHeap(t,(keytyp)Random(),H);
    }
    if(len <= 0) sites_error("AddRandomCLSites( ) input error");
    for(s=t; len > 0 ; s++,len--) insrtHeap(s,(keytyp)Random(),H);
    for(t=s=1; (i=delminHeap(H)) != NULL; ){
	    if(i<=S->ntyp){ AddSite(t,n,s,S); s+=S->len_elem[t]; t++; }
	    else s++;
    }
    Nildheap(H);
    return TRUE;
}

BooLean	ShuffleCLSites(st_type S)
// shuffle colinear sites retaining one site in each sequence.
// WARNING: must have one and only one of each type in each sequence.
{
    Int4	i,t,n,s,len;
    dh_type	H;
	
    n = MaxSeqSeqSet(SitesSeqSet(S));
    H=dheap(n+2,3);
    for(n = 1; n <= S->nseq; n++) {
	len = S->len_seq[n];
	for(t=1; t <= S->ntyp; t++){
	    if(nSites(t,n,S)!=1) sites_error("ShuffleCLSites: input error");
	    s=S->site_pos[t][n][1];
	    VacateSite(t, n, s, S);
	    len -= S->len_elem[t];
	    insrtHeap(t,(keytyp)Random(),H);
	}
        for(s=t; len > 0 ; s++,len--) insrtHeap(s,(keytyp)Random(),H);
        for(t=s=1; (i=delminHeap(H)) != NULL; ){
	    if(i<=S->ntyp){ AddSite(t,n,s,S); s+=S->len_elem[t]; t++; }
	    else s++;
	}
    }
    Nildheap(H);
#if 0
    PutTypeSites(stderr, S); 
    for(t=1; t <= S->ntyp; t++){ PutSites(stderr, t, S, NULL,NULL); }
#endif
    return TRUE;
}

/****************************** Archive Sites **************************/
sti_typ	ArchiveSites(st_type S)
{
	sti_typ	X;
	Int4	i,t,n,m;

    X = new sites_info_type[1]; // replaces NEW(X,1,sites_info_type);
    X->ntyp = S->ntyp;
    X->gss = S->gss; // Copy constructor called.
    NEW(X->len_elem, X->ntyp+2, Int4);
    NEWP(X->nsites,X->ntyp+2,Int4);
    NEWPP(X->site_pos, X->ntyp+2,unsigned short);
    ss_type data=X->gss.FakeSqSet();  // Note that this was initialized above.
    for(t = 1; t <= X->ntyp; t++) {
        X->len_elem[t] = S->len_elem[t];
        NEW(X->nsites[t], S->nseq+2,Int4);
	NEWP(X->site_pos[t], NSeqsSeqSet(data)+2,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		X->nsites[t][n] = S->nsites[t][n];
                m = ((Int4) SqLenSeqSet(n,data)/(Int4)S->len_elem[t]) +1;
                NEW(X->site_pos[t][n],m+2,unsigned short);
		for(i = 1; i <= S->nsites[t][n]; i++) {
			X->site_pos[t][n][i] = S->site_pos[t][n][i];
		}
	}
    }
    return X;
}

sti_typ	CopyArchiveSites(sti_typ S)
{
    sti_typ	X;
    Int4	i,t,n,m;

    X = new sites_info_type[1]; // replaces NEW(X,1,sites_info_type);
    X->ntyp = S->ntyp;
    X->gss = S->gss; // Copy constructor called.
    NEW(X->len_elem, X->ntyp+2, Int4);
    NEWP(X->nsites,X->ntyp+2,Int4);
    NEWPP(X->site_pos, X->ntyp+2,unsigned short);
    ss_type data=X->gss.FakeSqSet();  // Note that this was initialized above.
    for(t = 1; t <= X->ntyp; t++) {
        X->len_elem[t] = S->len_elem[t];
        NEW(X->nsites[t], NSeqsSeqSet(data)+2,Int4);
	NEWP(X->site_pos[t], NSeqsSeqSet(data)+2,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		X->nsites[t][n] = S->nsites[t][n];
                m = ((Int4) SqLenSeqSet(n,data)/(Int4)X->len_elem[t]) +1;
                NEW(X->site_pos[t][n],m+2,unsigned short);
		for(i = 1; i <= X->nsites[t][n]; i++) {
			X->site_pos[t][n][i] = S->site_pos[t][n][i];
		}
	}
    }
    return X;
}

st_type	ExtractSites(sti_typ X)
{
	st_type	S;
	Int4	t,n,i;

    S = MakeSites(X->ntyp, X->len_elem, X->gss);
    ss_type data=SitesSeqSet(S);
    for(t = 1; t <= X->ntyp; t++) {
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		for(i = 1; i <= X->nsites[t][n]; i++) {
			AddSite(t, n, X->site_pos[t][n][i], S);
		}
	}
    }
    return S;
}

void	NilArchiveSites(sti_typ X)
{
    Int4   t,n;

    ss_type data=X->gss.FakeSqSet();
    for(t = 1; t <= X->ntyp; t++) {
	for(n = 1; n <= NSeqsSeqSet(data); n++) free(X->site_pos[t][n]);
        free(X->nsites[t]); free(X->site_pos[t]);
    }
    free(X->len_elem); free(X->nsites); free(X->site_pos);
    delete []X; // replaces free(X);
}

/**************************** end Archive Sites ************************/
void	InitSites(st_type S)
/* inititalize S to contain no sites; i.e., vacate sites. */
{
    Int4		n,t,i;

    for(n = 1; n <= S->nseq; n++) {
	ClearOlist(S->pos[n]);
	for(i = 0; i<= S->len_seq[n]; i++) S->type[n][i] = VACANT;
    }
    for(t = 1; t <= S->ntyp; t++) {
	for(n = 1; n <= S->nseq; n++) {
	   S->nsites[t][n] = 0;
	   S->site_pos[t][n][1] = NULL;
	}
   }
}

st_type CopySites(st_type S)
// Create and return an exact copy of S.
{
    Int4	t,n,s;
    st_type	S2;

    S2 = MakeSites(S->ntyp,S->len_elem, S->gss); 
    for(t = 1; t <= S->ntyp; t++) {
	for(n = 1; n <= S->nseq; n++) {
	    for(s = 1; s <= S->nsites[t][n]; s++){
		  AddSite(t, n, S->site_pos[t][n][s], S2);
	    }
	}
    }
    return S2;
}

void	NilSites(st_type S)
/* destroy sites structure S. */
{
    Int4	t,n;

   free(S->len_elem); free(S->len_seq); free(S->tmp); free(S->totsites); 
   for(n = 1; n <= S->nseq; n++) { NilOlist(S->pos[n]); free(S->type[n]); }
   free(S->type); free(S->pos);
   for(t = 1; t <= S->ntyp; t++) {
   	for(n = 1; n <= S->nseq; n++) free(S->site_pos[t][n]);
	free(S->site_pos[t]); free(S->nsites[t]);
   }
   free(S->site_pos); free(S->nsites);
   delete []S; // replaces free(S);
}

Int4    *GetPosSites(Int4 sq,st_type S)
{
        Int4    b,*rtn,nblk=nTypeSites(S);
        assert(sq > 0 & sq <= NSeqsSites(S));
        NEW(rtn,nblk+3,Int4);
        for(b=1; b <= nblk; b++) rtn[b]=SitePos(b,sq,1,S);
        return rtn;
}

Int4	GetEdgeBlocksSite(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
// Get the edges that are blocked in S.
// Should run ExtendFakeToRealCMSA(cma): GK{()---  ->  {(GK)---- 
// ...prior to calling this to avoid problems.
// Should fix this routine instead eventually. 
{
        Int4    n,b,p,len,num_blk=0;

        len = SiteLen(t,S);
        for(n=1; n <=NSeqsSites(S); n++){
	   if(nSites(t,n,S)==1){
                p=SitePos(t,n,1,S);
                if(right) p += len; else p--;
                if(p > SeqLenSites(n,S) || p < 1 || !OpenPos(n,p,S)){
			blocked[n]='T'; num_blk++;
                } else blocked[n]='F';
		if(right) pos[n]=p-1; else pos[n]=p;
           } else blocked[n]='M'; // element is missing from seq.
        } return num_blk;
}
// WARNING: BOTH GetEdgeBlocksSite() AND AlwaysGetSiteFreq() probably
// SET POS TOO FAR OVER WHEN RIGHT == TRUE.

Int4	*AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
/*********************************************************************
   Note: From calling environment:
		 case 1: -inf...-1 is left of site.
		 case 2:  0...+inf is within or beyond site
 *********************************************************************/
{
        Int4    *site_freq,n,b,p,len;
        ss_type P=SitesSeqSet(S);
        e_type  E;

        len = SiteLen(t,S);
        MEW(site_freq, nAlpha(SeqSetA(P))+2,Int4);
        for(b=0;b<= nAlpha(SeqSetA(P)); b++) site_freq[b] = 0;
        for(n=1; n <=NSeqsSeqSet(P); n++){
           E=SeqSetE(n,P);
	   if(nSites(t,n,S)==1){
                p = SitePos(t,n,1,S);
                if(right) p += len; else p--;
                if(p > SqLenSeqSet(n,P) || p < 1 || !OpenPos(n,p,S)){
			blocked[n]='T'; site_freq[0]++; // 'X' residue.
                } else {
			blocked[n]='F';
			b = ResSeq(p,E); site_freq[b]++; 
		}
		if(right) pos[n]=p-1; else pos[n]=p;
           } else blocked[n]='M'; // element is missing from seq.
        }
        return site_freq;
}

Int4	*GetSiteFreq(st_type S, set_typ Set, Int4 t, Int4 d)
/*********************************************************************
   Note: From calling environment:
		 case 1: -inf...-1 is left of site.
		 case 2:  0...+inf is within or beyond site
   returns the residue frequencies for a site d residues right (or left
   if d negative) of sites of type t in S.  If position is blocked or 
   off ends of a sequence NULL is returned. 
   d = 0..len-1 means column is within the site.
   Option use_blocked = TRUE ignores blocked sites.
 *********************************************************************/
{
	BooLean	no_blocked=TRUE;
	Int4	*site_freq,n,k,b,p,len=SiteLen(t,S);
	ss_type	P=SitesSeqSet(S);
	e_type	E;

	MEW(site_freq, nAlpha(SeqSetA(P))+2,Int4);
	for(b=0;b<= nAlpha(SeqSetA(P)); b++) site_freq[b] = 0;
	for(n=1; n <=NSeqsSeqSet(P); n++){
	   if(!MemberSet(n,Set)) continue;
	   E = SeqSetE(n,P);
	   for(k=1; k<=nSites(t,n,S); k++){
		p = SitePos(t,n,k,S) + d;
		if(p > SqLenSeqSet(n,P) || p < 1){ free(site_freq); return (Int4*) NULL; }
		if(no_blocked && (d >= len || d < 0) && !OpenPos(n,p,S)){
			free(site_freq); return (Int4*) NULL;
		} b = ResSeq(p,E); site_freq[b]++;
	   }
	} return site_freq;
}

Int4	*GetSiteFreq(st_type S, BooLean *skip, Int4 t, Int4 d)
/*********************************************************************
   Note: From calling environment:
		 case 1: -inf...-1 is left of site.
		 case 2:  0...+inf is within or beyond site
   returns the residue frequencies for a site d residues right (or left
   if d negative) of sites of type t in S.  If position is blocked or 
   off ends of a sequence NULL is returned. 
   d = 0..len-1 means column is within the site.
   Option use_blocked = TRUE ignores blocked sites.
 *********************************************************************/
{
	BooLean	no_blocked=TRUE;
	Int4	*site_freq,n,k,b,p,len;
	ss_type	P=SitesSeqSet(S);
	e_type	E;

	len = SiteLen(t,S);
	MEW(site_freq, nAlpha(SeqSetA(P))+2,Int4);
	for(b=0;b<= nAlpha(SeqSetA(P)); b++) site_freq[b] = 0;
	for(n=1; n <=NSeqsSeqSet(P); n++){
	   if(skip[n]) continue;
	   E = SeqSetE(n,P);
	   for(k=1; k<=nSites(t,n,S); k++){
		p = SitePos(t,n,k,S) + d;
		if(p > SqLenSeqSet(n,P) || p < 1){ free(site_freq); return (Int4*) NULL; }
		if(no_blocked && (d >= len || d < 0) && !OpenPos(n,p,S)){
			free(site_freq); return (Int4*) NULL;
		} b = ResSeq(p,E); site_freq[b]++;
	   }
	} return site_freq;
}

Int4	*GetSiteFreq(st_type S, Int4 t, Int4 d)
/*********************************************************************
   Note: From calling environment:
		 case 1: -inf...-1 is left of site.
		 case 2:  0...+inf is within or beyond site
 *********************************************************************/
{ return get_site_freq(S, t, d, TRUE); }

Int4	*get_site_freq(st_type S, Int4 t, Int4 d, BooLean no_blocked)
/*********************************************************************
   returns the residue frequencies for a site d residues right (or left
   if d negative) of sites of type t in S.  If position is blocked or 
   off ends of a sequence NULL is returned. 
   d = 0..len-1 means column is within the site.
   Option use_blocked = TRUE ignores blocked sites.
 *********************************************************************/
{
	Int4	*site_freq,n,k,b,p,len;
	ss_type	P=SitesSeqSet(S);
	e_type	E;

	len = SiteLen(t,S);
	MEW(site_freq, nAlpha(SeqSetA(P))+2,Int4);
	for(b=0;b<= nAlpha(SeqSetA(P)); b++) site_freq[b] = 0;
	for(n=1; n <=NSeqsSeqSet(P); n++){
	   E = SeqSetE(n,P);
	   for(k=1; k<=nSites(t,n,S); k++){
		p = SitePos(t,n,k,S) + d;
		if(p > SqLenSeqSet(n,P) || p < 1){
			free(site_freq); return (Int4*) NULL;
		}
		if(no_blocked && (d >= len || d < 0) && !OpenPos(n,p,S)){
			free(site_freq); return (Int4*) NULL;
		}
		b = ResSeq(p,E); site_freq[b]++;
	   }
	}
	return site_freq;
}

void	ShiftSitesM(st_type S, Int4 t, Int4 d)
/* shift all type t sites d spaces to the left(?). */
{
	Int4	i;
	BooLean	left;
	// char	c; 

	if(d < 0){ /** c = '-'; /**/ d *= -1; left=FALSE; }
	else if(d > 0){ /** c = '+';/**/ left = TRUE; }
	else return;
	// fprintf(stderr,"[");
	for(i=1; i <=d; i++) {
		// fprintf(stderr,"%c",c);
		ShiftSites(S, t, left);
	}
	// fprintf(stderr,"]"); 
}

void	ShiftSites(st_type S, Int4 t, BooLean left)
/*********************** shift left *******************
 Shift all type t sites one position to the left or right.
 	WARNING: assumes that new positions are available.
 SHIFT LEFT:

     site 		  type[n][site] = VACANT ('o')
      |			  type[n][site + 1] = t; ('A')
      A   a   a   a   ...   a   o       	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-1| w |
      o   A   a   a   ...   a   a       	
          |		  type[n][site + w] = BLOCKED (t)
        site + 1
	
 SHIFT RIGHT:

         site 		  type[n][site] = BLOCKED (t)
          |		  type[n][site - 1] = t; ('A')
      o   A   a   a   ...   a   a       	<- sequence n
    |-1 |+0 |+1 |+2 | ... |w-2|w-1|
      A   a   a   a   ...   a   o       	
      |			  type[n][site + w - 1] = VACANT ('o')
     site - 1

WARNING: This operation (apparently) won't work if the model 
  lengths are < 2.
*******************************************************/
{
	Int4	n,k,site;

    if(left){		/* free 1; 1..w-1 = 2..w; w = new */
#if 0
	fprintf(stderr,"\nmodel %d ********\n", t);
	PutTypeSites(stderr, S); PutSites(stderr, t, S, NULL,NULL);
#endif
	for(n=1; n<= S->nseq; n++){
	 if(S->nsites[t][n] > 0){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
	   	site=S->tmp[k];
		if(S->type[n][site]==t){
		   RmOlist(site,S->pos[n]);	   /* move pattern left */
		   InsertOlist(site +1, S->pos[n]);
		}
	   }
	   bubble_sort_sites(S->site_pos[t][n]);
	   for(k=S->nsites[t][n]; k > 0; k--){
	   	site = S->site_pos[t][n][k];
		S->type[n][site] = VACANT; 
		S->type[n][site+1] = t;	 	/* move over +1 */
		if(S->type[n][site+S->len_elem[t]] != VACANT){
			PutTypeSites(stderr, S);
			PutSites(stderr, t, S, NULL,NULL);
			fprintf(stderr,"n = %d; site = %d; t = %d; len = %d\n",
				n,site,t,S->len_elem[t]);
			sites_error("shift left operation is blocked.");
		}
		S->type[n][site+S->len_elem[t]] = BLOCKED(t);
		S->site_pos[t][n][k]++; /** move all pos. ahead one **/
	   }
	 }
	}
#if 0
	PutTypeSites(stderr, S); PutSites(stderr, t, S, NULL,NULL);
#endif
    } else {		/* free w; w..2 = w-1..1; 1 = new */
	for(n=1; n<= S->nseq; n++){
	 if(S->nsites[t][n] > 0){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t){
		   RmOlist(site,S->pos[n]);	   /* move pattern left */
		   InsertOlist(site-1, S->pos[n]);
		}
	   }
	   bubble_sort_sites(S->site_pos[t][n]);
	   for(k=1; k <= S->nsites[t][n]; k++){
	   	site = S->site_pos[t][n][k];
		if(S->type[n][site-1] != VACANT)
			sites_error("shift right operation is blocked.");
		S->type[n][site-1] = t;	   /* move over -1 */
		S->type[n][site] = BLOCKED(t);
		S->type[n][site+S->len_elem[t]-1] = VACANT;
		S->site_pos[t][n][k]--;  /** move all pos. back one **/
	   }
	 }
	}
    }
}

void    GrowSites(Int4 t, st_type S)
/*********************** grow right *************************
 Lengthen all type t sites one position to the right.
	 WARNING: assumes that new positions are available.
 GROW RIGHT:
     site 
      |
      A   a   a   a   ...   a   o    	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-1|	w |	 len_elem[t]++;
      A   a   a   a   ...   a   a       	
           		 type[n][site + w] = BLOCKED(t)
***************************************************************/
{
	Int4	w,n,k,site;

	w = S->len_elem[t]++;		/* w = old length */
	for(n=1; n<= S->nseq; n++){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t){
		   if(S->type[n][site + w] != VACANT){
			fprintf(stderr,"site=%d; t=%d; n=%d\n",site,t,n);
			PutTypeSites(stderr, S);
			PutSites(stderr, t, S, NULL,NULL);
			sites_error("grow operation is blocked.");
		   }
		   S->type[n][site + w] = BLOCKED(t);
		}
	   }
	}
	for(S->maxinc=w+1,t=1; t<=S->ntyp; t++) 
		S->maxinc = MINIMUM(Int4, S->maxinc,(S->len_elem[t]-1));
	S->maxinc = MAXIMUM(Int4,1,S->maxinc);	// afn: 2-21-08
}

void    ShrinkSites(Int4 t, st_type S)
/*********************** shift left *******************
 Shortens all type t sites one position on the right.
 SHRINK RIGHT:
     site 
      |
      A   a   a   a   ...   a   a      	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|		 len_elem[t]--;
      A   a   a   a   ...   a   o           	
           		 type[n][site + w - 1] = VACANT 

*******************************************************/
{
	Int4	w,n,k,site;

	w = S->len_elem[t]--;		/* w = old length */
	if(w<2) sites_error("cannot shrink element to length < 2");
	for(n=1; n<= S->nseq; n++){
           GetOlist(S->tmp,S->pos[n]);
	   for(k=1;(site=S->tmp[k]) != 0; k++){
		if(S->type[n][site]==t) S->type[n][site + w-1] = VACANT;
	   }
	}
	S->maxinc = MINIMUM(Int4, S->maxinc,(S->len_elem[t]-1));
	S->maxinc = MAXIMUM(Int4,1,S->maxinc);	// afn: 2-21-08
}

void	VacateSites(Int4 n, st_type S)
// remove all sites from nth sequence.
{
   Int4		i,t,s,k;

   for(t=1; t<=nTypeSites(S); t++){
	k = nSites(t,n,S);
	for(i=1;i<=k;i++){ s=S->site_pos[t][n][i]; VacateSite(t,n,s,S);}
   }
}

void	VacateSite(Int4 t, Int4 n, Int4 site, st_type S)
/*********************** vacate site *******************
 Remove site in sequence n of type t by vacating all positions. 
     site 
      |
      A   a   a   a   ...   a   a      	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      o   o   o   o   ...   o   o           	
        for p = site ... site + w - 1 -> type[n][p] = VACANT 
*******************************************************/
{
	Int4	i,end;
	unsigned short  *site_pos=S->site_pos[t][n];

	if(S->type[n][site] != t){
		PutTypeSites(stderr, S);
		PutSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,"ELEMENT %c; seq %d; site %d\n",
			'A' +t -1, n,site);
		sites_error("attempt to remove site where none exists.");
	}
	RmOlist(site,S->pos[n]);
	end = site + S->len_elem[t] -1; 
	for(i=site; i <= end ; i++) S->type[n][i] = VACANT; 
/****************************************************************
            i= 1   2   3   4    (nsites = 4)       
  remove(35, [12, 35, 78, 104, NULL]) 

	-> [12, 104, 78, 104, NULL] -> [12, 104, 78, NULL ] 

  remove(35, [35, NULL]) 

	-> [35, NULL] -> [NULL]

 ****************************************************************/
	for(i=1; i <= S->nsites[t][n]; i++){
		if(site_pos[i] == site){
			site_pos[i] = site_pos[S->nsites[t][n]];
			site_pos[S->nsites[t][n]] = NULL;
			S->nsites[t][n]--;
			S->totsites[t]--;
			return ;
		}
	}
	sites_error("VacateSite( ) - this should not happen");
}

//============ afn: 7/29/2014 ================
Int4	MkRoomForSite(Int4 blk, Int4 sq, Int4 site, st_type S)
// If new site is blocked on ends, then lengthen FakeSeq with gaps '-' to make room.
// check for ends only!!!
#if 0
		           Ccccc
          sq: .....AaaaBbbb...--
#endif
{
	Int4	end,gap=0;
        assert(sq > 0 & sq <= NSeqsSites(S));
	end = site + S->len_elem[blk] - 1;	// site...end = region for new block.
        if(end > S->len_seq[sq]){	// need to lengthen sequence...
	    Int4  k,b,nblk,*rtn;
	    nblk=nTypeSites(S);
            NEW(rtn,nblk+3,Int4);
            for(b=1; b <= nblk; b++){
		k = nSites(b,sq,S); assert(k==0 || k==1);
		if(k == 1) rtn[b]=SitePos(b,sq,1,S); else rtn[b]=-1;
	    }
            VacateSites(sq,S);
	    gap= end - S->len_seq[sq];
	    gss_typ *gss=SitesGSS(S);

	    // gss->FakeSqSet();	// makes the fake seq set if currently null.
            gsq_typ *gsq=gss->GetGSQ(sq);
            // gsq->Put(stderr,AB);
	    e_type E = 0;
	    if(gsq) E=gsq->FakeSeq(); else E=gss->TrueSeq(sq);
	    Int4 p=LenSeq(E);

	    Int4 m,max1,max2,len2;
	    max1=MaxSeqSeqSet(SitesSeqSet(S));
	    // fprintf(stderr,"sq=%d; p=%d; len=%d; gap=%d\n", sq,p,LenSeq(E),gap);
	    gss->InsertGap(sq,p,gap);

	    // Need to reset Site parameters for sq.
	    E=S->gss.FakeSeq(sq);
	    // PutSeq(stderr,E,AlphabetSites(S));
	    len2=LenSeq(E);
            max2=MaxSeqSeqSet(SitesSeqSet(S));
            if(max2 > max1){ free(S->tmp); NEW(S->tmp,max2+2,Int4); }
            S->len_seq[sq]=len2;
            free(S->type[sq]); NEW(S->type[sq],(S->len_seq[sq]+2),char);
            NilOlist(S->pos[sq]);
            S->pos[sq] = Olist(S->len_seq[sq]+1);
            S->type[sq][(S->len_seq[sq]+1)] = ENDTYPESITE;
            for(b=1; b<=nTypeSites(S); b++){
              free(S->site_pos[b][sq]);
              m = ((Int4) len2/(Int4)S->len_elem[b]) + 3;
              NEW(S->site_pos[b][sq],m+2,unsigned short);
            }

            // gsq->Put(stderr,AB);
            for(b=1; b <= nblk; b++){
	   	if(rtn[b] > 0) AddSite(b,sq,rtn[b],S);
    	    } free(rtn);
	} 
#if 0	// later add gaps for other issues if they arise...
	for(p=site; p < end ; p++){
		else if(S->type[sq][p]) gap++; 
	{
        if(end > S->len_seq[sq]){	
                assert(end <= S->len_seq[sq]);   
        }
        if(S->type[sq][end]){	// blocked on end...
		 return TRUE;
	}
        return FALSE;
#endif
	return gap;
}

BooLean	IsBlockedSite(Int4 blk, Int4 sq, Int4 site, st_type S)
// Return TRUE if new site is blocked else return FALSE.
{
	register Int4	p,end;
	end = site + S->len_elem[blk] - 1;	// site...end = region for new block.
        if(end > S->len_seq[sq]) return TRUE;
	for(p=site; p < end ; p+=S->maxinc) if(S->type[sq][p]) return TRUE;
        if(S->type[sq][end]) return TRUE;
        else return FALSE;
}

//============ afn: 7/29/2014 ================

void	AddSite(Int4 blk, Int4 sq, Int4 site, st_type S)
/*************************** add site ********************************
 Add type t site in sequence n.
     site 
      |
      o   o   o   o   ...   o   o           	
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      A   a   a   a   ...   a   a      	<- sequence n
	type[n][site] = t  ('A')
        for p = site + 1 ... site + w - 1 -> type[n][p] = BLOCKED (t) 
***********************************************************************/
{
	Int4	p,end;

// fprintf(stderr,"DEBUG: OccupiedSite()\n");
	if(IsBlockedSite(blk,sq,site,S)) MkRoomForSite(blk,sq,site,S);
// assert(site>0);  // for debugging omission of sequences from cmsa...adding to 0 when abscent!
	if(OccupiedSite(blk, sq, site, S)){
		PutTypeSites(stderr, S);
		PutSites(stderr, blk, S, NULL,NULL);
		fprintf(stderr,
			"ELEMENT %c(length=%d); seq %d; site %d\n",
			'A' + blk -1, S->len_elem[blk], sq,site);
		S->gss.Put(stderr,sq);
		// S->gss.PutFA(stderr);
		fprintf(stderr,
			"sequence length = %d; numseqs = %d\n", 
			SeqLenSites(sq,S),S->nseq);
		assert(!OccupiedSite(blk, sq, site, S));
		sites_error("attempt to add site where one exists.");
	}
// fprintf(stderr,"DEBUG sites 1 (S->pos[%d]=%d)\n",n,S->pos[n]);
	InsertOlist(site, S->pos[sq]);
	S->type[sq][site] = blk;
	end = site + S->len_elem[blk] - 1; 
// fprintf(stderr,"DEBUG sites 1 (%d..%d)\n",site+1,end);
	for(p=site+1; p <= end ; p++) S->type[sq][p] = BLOCKED(blk); 
// fprintf(stderr,"DEBUG sites 2\n");
	S->nsites[blk][sq]++;
	S->totsites[blk]++;
	S->site_pos[blk][sq][S->nsites[blk][sq]] = (unsigned short) site;
}

BooLean OccupiedSite(register Int4 t, register Int4 n, register Int4 site, 
	register st_type S)
/* determine if a site is blocked. */
{
	register Int4	p,end;

	end = site + S->len_elem[t] - 1; 
// fprintf(stderr,"DEBUG: OccupiedSite() (%d..%d; %d)\n",site,end,S->maxinc);
	for(p=site; p<end ; p+=S->maxinc){if(S->type[n][p])return TRUE;}
	if(end > S->len_seq[n]){
		ss_type	data = SitesSeqSet(S);
		PutSeqSetE(stderr,n,data);
		fprintf(stderr,"block=%d; seq=%d; site=%d\n",t,n,site);
		assert(end <= S->len_seq[n]);	// afn:2/19/08 - Valgrind complaining...
	}
	if(S->type[n][end]) return TRUE;
	return FALSE;
}

Int4	PosTSites(Int4 t, Int4 n, Int4 *pos, st_type S)
/* modifies array pos to contain the positions of type t sites in seq. n */
{
	Int4	s,site,i;

	GetOlist(S->tmp, S->pos[n]); 
	for(i=0,s=1; (site=S->tmp[s]) != 0; s++) {
		if(S->type[n][site]==t){ i++; pos[i] = S->tmp[s]; }
	} pos[i+1] = 0;
	return i;
}

Int4	PosSites(Int4 n, Int4 *pos, st_type S)
/* modifies array pos to contain the positions of sites in sequence n */
{
	Int4	s;
	GetOlist(S->tmp, S->pos[n]); 
	for(s=1; S->tmp[s] != 0; s++) pos[s] = S->tmp[s];
	pos[s] = 0;
	return (s-1);
}

Int4	GapBetweenSites(Int4 n, Int4 t, st_type S)
// return the gap length between sites
{
	assert(t >= 0 && t <= S->ntyp);
	if(t==0) return (SitePos(1,n,1,S) - 1);
	else if(t==S->ntyp) return (SeqLenSites(n,S) - EndSitePos(t,n,1,S));
	else return (SitePos(t+1,n,1,S) - EndSitePos(t,n,1,S) - 1);
}

/********************************** output ******************************/
void	PutSites(FILE *fptr,Int4 t,st_type S,double **site_prob, BooLean *off)
{
	ss_type	P=SitesSeqSet(S);
	Int4	i,n,s,e,end,N,start,length;
	BooLean	are_sites = FALSE;
	e_type	E;
	char	r,c;
        Int4	flank=10; 
   
	fprintf(fptr,"\n\n");
	length = SiteLen(t,S);
	for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
	   bubble_sort_sites(S->site_pos[t][n]);
	   E=SeqSetE(n,P);
	   if(nSites(t,n,S) > 0) N++; 
	   for(s=1; s<= nSites(t,n,S); s++){
		are_sites = TRUE;
	   	fprintf(fptr,"%2d-%-2d ",n,s);
		start= S->site_pos[t][n][s]; 
		fprintf(fptr,"%4d  ",start);
		e = start + length - 1;
		end = e + flank;
		for(i=start-flank; i <= end; i++){
		   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
		   else{
			r = ResSeq(i,E);
			if(OpenPos(n,i,S)) c = AlphaCharLow(r,SeqSetA(P));
			else c = AlphaChar(r,SeqSetA(P));
			if(i == e) fprintf(fptr,"%c ", c);
			else if(i == start) fprintf(fptr," %c", c);
			else fprintf(fptr,"%c", c);
		   }
		}
		fprintf(fptr," %4d",e);
		if(site_prob != NULL) {
			fprintf(fptr," (%0.2f)",site_prob[n][start]);
			// PutSeqID(fptr,E); /*** TEMP ***/
		} 
		fprintf(fptr,"\n");
	   }
	}
	if(are_sites){
	   if(off != NULL){
	   	fprintf(fptr,"sites:%17c",' ');
	   	for(s=1; s <= S->len_elem[t]; s++){
			if(off[s])fprintf(fptr," "); 
			else fprintf(fptr,"*"); 
	  	 }
	   }
	   fprintf(fptr,"\n%*s", 23, "");
	   for(i=1; i<= SiteLen(t,S); i++) {
		if(i%5==0) fprintf(fptr,"%5d",i);
	   }
	}
	fprintf(fptr,"\n\t(%d sites in %d sequences)\n\n",S->totsites[t],N);
}

void    PutTypeSites(FILE *fptr, st_type S)
/*  	e.g.    3: D(24)-A(35)-C(45)-[59]  */
{
	Int4 i,n,t,tot;

	fprintf(fptr,"\n");
	for(n=1; n<=S->nseq; n++){
	   for(tot=0,t=1; t<=S->ntyp; t++) tot += S->nsites[t][n];
	   if(tot > 0){
	      fprintf(fptr,"%3d: ",n);
	      for(i=1; i<=S->len_seq[n]; i++){
		if(S->type[n][i] > 0){
		   t = S->type[n][i];
	   	   if(tot == 0) fprintf(fptr,"-"); tot = 0;
		   if(S->ntyp <= 26){
			fprintf(fptr,"%c(%d)",('A'+t-1),i);
		   } else { fprintf(fptr,"%d(%d)",t,i); }
		}
	      }
	      fprintf(fptr,"-[%d]\n",S->len_seq[n]);
	   }
	}
	fprintf(fptr,"\n");
}

void	PutScanSites(FILE *fptr, Int4 t, st_type S, BooLean *off)
{ PutScanSitesProb(fptr, t, S, NULL, off, 0.0); }

void	PutScanSitesProb(FILE *fptr, Int4 t, st_type S, double **prob, BooLean *off,
	double cutoff)
/** Create a Scan file (*.sn) eliminating poor matches from the aligned block **/
{
        ss_type P=SitesSeqSet(S);
        Int4    i,n,s,e,end,N,start,length;
        e_type  E;
        char    r,c;

        if(off != NULL){
                for(s=1; s <= S->len_elem[t]; s++){
                        if(off[s]==TRUE)fprintf(fptr,".");
                        else if(off[s]==FALSE)fprintf(fptr,"*");
                        else fprintf(fptr,"^");
                 }
        } else for(s=1; s <= S->len_elem[t]; s++) fprintf(fptr,"*");
        fprintf(fptr,"\n");
        length = SiteLen(t,S);
        for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
           bubble_sort_sites(S->site_pos[t][n]);
           E=SeqSetE(n,P);
           if(nSites(t,n,S) > 0) N++;
           for(s=1; s<= nSites(t,n,S); s++){
             start= S->site_pos[t][n][s];
             if(prob == NULL || prob[n][start] >= cutoff) {
                e = start + length - 1;
                end = e;
                for(i=start; i <= end; i++){
                   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        c = AlphaChar(r,SeqSetA(P));
                        fprintf(fptr,"%c", c);
                   }
                }
                fprintf(fptr,"\n");
             }
           }
        }
        fprintf(fptr,"\n");
}

/********************************** private ******************************/
void	print_sites(unsigned short *L)
{
	Int4	i;
	L++;
	fprintf(stderr,"\n");
	for(i=0; L[i]!=NULL; i++) fprintf(stderr," %d",L[i]);
	fprintf(stderr," (null)\n");
}

Int4     bubble_sort_sites(unsigned short *L)
{
        unsigned short i,j,t,n;
	Int4 temp=0;

	L++;
        for(n=i=0; L[i]!=NULL; i++) n++;
        for(i=1; i < n; i++){
                for(j=i;j > 0 && (L[j-1] > L[j]);j--){
                        t=L[j]; L[j]=L[j-1]; L[j-1]=t;
                }
        }
	return temp;
}

void	sites_error(const char *s){fprintf(stderr,"sites error: %s\n",s);exit(1);}

void    PutPhylipSites(FILE *fptr,st_type S, BooLean **off)
/** PutSites to create a PHYLIP alignment input file **/
/** WARNING: this assumes that sites are colinear (checked) **/
/** WARNING: this assumes that off has the correct dimensions (not checked) **/
{
    ss_type	P=SitesSeqSet(S);
    Int4	i,j,n,s,t,r,end,N,start,length,ntyp;
    e_type	E;
    a_type	A=SeqSetA(P);
    char	c,*str;

    if(!ColinearSites(S)) sites_error("PutPhylipSites( ) input error");
    ntyp = nTypeSites(S);
    for(length=0, t=1; t <= ntyp; t++) {
	 for(i=1; i<=SiteLen(t,S); i++) {
		if(off == NULL || !off[t][i]) length++;
	}
    }
    N = NSeqsSeqSet(P);
    fprintf(fptr,"  %d  %d\n",N,length);
    for(t=1; t <= ntyp; t++){
	length = SiteLen(t,S);
	for(n=1; n<= N; n++){
	   E=SeqSetE(n,P);
	   if(t==1){
		for(str = SeqKey(E), j=0; j < 10; j++){
		   if(str[j]==0 || isspace(str[j])){
			do { fprintf(fptr," "); j++; } while(j < 10); break;
		   } else fprintf(fptr,"%c",str[j]);
		}
	   }
	   bubble_sort_sites(S->site_pos[t][n]);
	   if(nSites(t,n,S) > 1) sites_error("PutPhylipSites( ) input error");
	   for(s=1; s<= nSites(t,n,S); s++){
		start= S->site_pos[t][n][s];
		end = start + length - 1;
		for(j=1,i=start; i <= end; i++,j++){
			if(i < 1 || i > (Int4) LenSeq(E)) 
				sites_error("PutPhylipSites( ) error");
			else if(off== NULL || !off[t][j]){
			   r = ResSeq(i,E); 
			   fprintf(fptr,"%c", AlphaChar(r,A));
			}
		}
		fprintf(fptr,"\n");
	    }
	}
	fprintf(fptr,"\n");
    }
    fprintf(fptr,"\n");
}

