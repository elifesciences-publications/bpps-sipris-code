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

/* purge.c - homology purge program */
#include "purge.h"

Int4	PurgeFiles(char *infile, char *outfile, Int4 cutoff, Int4 inc, 
	Int4 minseq, Int4 maxseq, a_type A)
{
    Int4	N,s,m,n,c,min;
    ss_type	data;
    FILE	*fp;
    e_type	*List;
    BooLean	*Use;

    /*** data = SeqSet(infile,A); /******/
    data = MkXSeqSet(infile,A);
    N = NSeqsSeqSet(data);
    NEW(List,N+3,e_type);
    for(s=1; s <= N; s++){ List[s] = SeqSetE(s,data); }
    if(N > 100){
	/*******************************************************
	If more than 100 sequences then (to speed up purge) first remove 
	at cutoff=500 with T = 18; then purge at cutoff and write to file.
	/*******************************************************/
	c = MAXIMUM(Int4,cutoff,500);
	Use=RmHomologsList(c,'B',1,N, FALSE, 18, N, List, A);
        for(n=0,s=1; s <= N; s++){ 
	    if(Use[s]){ n++; List[n] = SeqSetE(s,data); }
	}
	free(Use);
    } else { n = N; }
    for(c=cutoff; TRUE; ){
	   Use=RmHomologsList(c,'B',1,n, FALSE, 11, n, List, A);
           for(m=0,s=1; s <= n; s++) if(Use[s]) m++;
	   c-=inc;
	   if(m <= maxseq || c < 20) break; else free(Use);
    }
    c+=inc;
    /******* NEW *******/
    if(m < minseq && c==cutoff){  /** then try to bring up to minseq **/
       for(c=cutoff+inc; TRUE; ){
	   Use=RmHomologsList(c,'B',1,n, FALSE, 11, n, List, A);
           for(m=0,s=1; s <= n; s++) if(Use[s]) m++;
	   c+=inc;
	   if(m >= minseq || c > 300) break; else free(Use);
       }
    }
    /******* NEW *******/
    /******* Write to output file *******/
    fp = open_file(outfile,"","w");
    min = MaxSeqSeqSet(data);
    for(s=1; s <= n; s++){ 
	    if(Use[s]){ 
		PutSeq(fp,List[s],A); 
		min = MINIMUM(Int4,min,LenSeq(List[s]));
	    }
    }
    fclose(fp); free(Use); free(List);
    NilSeqSet(data); 
    return m;
}

BooLean	*RmHomologsList(Int4 cutoff, char method, Int4 minimum, Int4 maximum,
	BooLean query, Int4 T, Int4 N, e_type *List, a_type A)
/**************************************************************************
  Purge from a list of sequences.
/**************************************************************************/
{
	Int4	maxleng,i,k,entry,*lengs;
	Int4	item,item1,*hits,num_seqs=0;
	b_type	*edges;
	dh_type	H;
	gb_typ	B;
	keytyp	key;
	e_type	E,E1,E2;
	BooLean	*itemList;
	BooLean	related,report=TRUE;

	if(method == 'B') { report= FALSE; method = 'b'; }
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); NEW(lengs,N+2,Int4); NEW(hits,N+2,Int4);
	for(maxleng=0,entry =1 ;entry <= N;entry++) {
		E=List[entry]; 
		lengs[entry]=LenSeq(E); edges[entry]=Block(N+1);
		maxleng = MAXIMUM(Int4,maxleng,LenSeq(E));
	}
	maxleng *= 2;
	entry = 1;
	if(cutoff < 99999) {
	  for(item =1 ;item <= N;item++) {
	    E1=List[item];
	    if(method == 'b') B=MakeGBlast(T, E1, A);
	    for(item1 =item+1 ;item1 <= N;item1++) {
	        E2=List[item1];
		if(method == 'h'){	/** fasta heuristic method **/
		     related = RelateSeqFastp(E1,E2,A,cutoff);
		} else if(method == 'b'){ /** gblast heuristic method **/
			related = FastMatcherGBlast(E2, B, cutoff);
/*****
			if(item==3 && item1==40) MatcherGBlast(stderr,E2, B);
/*****/
		} else related = RelatedSeqs(cutoff, E1, E2, A);
		if(related){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    }
	    if(method == 'b') NilGBlast(B);
	    if(report) fprintf(stderr,"\r%d",entry++);
	    else entry++;
	  }
	}
	fprintf(stderr,"\n%d items compared; cutoff %d\n", entry-1,cutoff); 
 if(!query){
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		if(report && hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key); 
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RmHomologs( )");
	    else if(minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=1;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
  } else {	/** don't purge #1 from set; i.e., don't put in heap. **/
	for(entry=2; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		if(report && hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key); 
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H) || minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=2;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
	insrtHeap(1,-1e99,H);
  }
	num_seqs = ItemsInHeap(H);
	if(query) num_seqs++;	/** i.e., first not in heap **/
	NEW(itemList,N+3,BooLean);
	if(num_seqs >= minimum && num_seqs <= maximum){
	  for(k=0,entry=1; entry <= N; entry++){
	    E=List[entry]; 
	    if(memHeap(entry,H)){ itemList[entry]=TRUE; k++; }
	    else itemList[entry] = FALSE;
	  }
	}
	fprintf(stderr,"%d sequences left after purging\n",k);
	for(entry=1; entry <= N; entry++) { NilBlock(edges[entry]); }
	free(lengs); free(hits); free(edges); Nildheap(H);
	return itemList;
}

Int4	RmHomologs(Int4 cutoff, char method, Int4 minimum, Int4 maximum,
	BooLean query, ss_type P, char *outfile, Int4 T)
{
	BooLean	*use;
	Int4	i,N,k;
	e_type  *List;
	FILE    *fptr;

	N = NSeqsSeqSet(P);
	NEW(List, N+2, e_type);
	for(i=1; i<=N; i++) List[i] = SeqSetE(i,P);
	use = RmHomologsList(cutoff, method, minimum, maximum, query, 
		T, N, List, SeqSetA(P));
	fptr = open_file(outfile,"","w");
	for(i=1,k=0; i<= N; i++){
	    if(use[i]) { PutSeqSetE(fptr,i,P); k++; }
	}
	fclose(fptr);
	free(use); free(List);
	return k;
}

Int4	RmHomologs2(Int4 cutoff, char method, Int4 minimum, Int4 maximum,
	BooLean query, ss_type P, char *outfile, Int4 T)
/***	Use block bit array for edges to save space.  Return the number of 
	sequences in the resultant set. ***/
{
	Int4	maxleng,i,k,entry,*lengs;
	Int4	item,item1,*hits,N,num_seqs=0;
	b_type	*edges;
	FILE	*fptr;
	dh_type	H;
	gb_typ	B;
	keytyp	key;
	e_type	E,E1,E2;
	BooLean	related,report=TRUE;

	if(method == 'B') { report= FALSE; method = 'b'; }
	N = NSeqsSeqSet(P);
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); NEW(lengs,N+2,Int4); NEW(hits,N+2,Int4);
	for(maxleng=0,entry =1 ;entry <= N;entry++) {
		E=SeqSetE(entry,P);
		lengs[entry]=LenSeq(E);
		edges[entry] = Block(N+1);
		maxleng = MAXIMUM(Int4,maxleng,LenSeq(E));
	}
	maxleng *= 2;
	entry = 1;
	if(cutoff < 99999) {
	  for(item =1 ;item <= N;item++) {
	    E1=SeqSetE(item,P);
	    if(method == 'b') B=MakeGBlast(T, E1, SeqSetA(P));
	    for(item1 =item+1 ;item1 <= N;item1++) {
	        E2=SeqSetE(item1,P);
		if(method == 'h'){	/** fasta heuristic method **/
		     related = RelateSeqFastp(E1,E2,SeqSetA(P),cutoff);
		} else if(method == 'b'){ /** gblast heuristic method **/
			related = FastMatcherGBlast(E2, B, cutoff);
/*****
			if(item==3 && item1==40) MatcherGBlast(stderr,E2, B);
/*****/
		} else related = RelatedSeqs(cutoff, E1, E2, SeqSetA(P));
		if(related){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    }
	    if(method == 'b') NilGBlast(B);
	    if(report) fprintf(stderr,"\r%d",entry++);
	    else entry++;
	  }
	}
	fprintf(stderr,"\n%d items compared; cutoff %d\n", entry-1,cutoff); 
 if(!query){
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		if(report && hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key); 
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RmHomologs( )");
	    else if(minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=1;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
  } else {	/** don't purge #1 from set; i.e., don't put in heap. **/
	for(entry=2; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		if(report && hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key); 
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H) || minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=2;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
	insrtHeap(1,-1e99,H);
  }
	num_seqs = ItemsInHeap(H);
	if(query) num_seqs++;	/** i.e., first not in heap **/
	if(num_seqs  < N){
	 if(num_seqs >= minimum && num_seqs <= maximum){
	  fptr = open_file(outfile,"","w");
	  for(entry=1,k=0; entry <= N; entry++){
	    E=SeqSetE(entry,P); 
	    if(memHeap(entry,H)) { PutSeqSetE(fptr,entry,P); k++; }
	  }
	  fclose(fptr);
	  fprintf(stdout,"%d sequences put into file: '%s'\n",k,outfile);
	 }
	}
	for(entry=1; entry <= N; entry++) { NilBlock(edges[entry]); }
	free(lengs); free(hits); free(edges); Nildheap(H);
	return num_seqs;
}

e_type	*PurgeScanInfoRpts(Int4 cutoff,Int4 Nt_flank,Int4 Ct_flank,snh_typ sH,a_type A)
/********************************************************************************
  Purge a set of N sequences using the motif information to locate and compare
  similar regions.

  1. Cluster sequences at some cutoff level.  (cutoff is in % identity).
  2. Pick a representative sequence from each set.
	a. Purge cluster of sequences: take the one sequence that 
	   has the highest aggregate score with all other sequences?

/********************************************************************************/
{
        Int4    i,j,k,b,s1,s2,leng,num_seqs,*sum,max;
	Int4	m,M,n,s,score,set1,set2,ident;
        e_type  E1,E2,*E,*tmpE,*rtnE=0;
        unsigned char    *seq1,*seq2;
	char	**R=AlphaR(A);
        sni_typ I1,I2,*tmpI;
	sni_typ *rtn;
	h_type	HG;
	ds_type	sets;
	BooLean	verbose=TRUE;
	Int4	N,nblk,*blklen;

	print_error("PurgeScanInfoRpts( ) needs work...");
	sni_typ *I = DomainsScanHeap(&N, Nt_flank, Ct_flank, sH);
        nblk = nblksScanHeap(sH);
        blklen = LengthsScanHeap(sH);
	for(leng=0,b=1 ;b<=nblk; b++) { leng+=blklen[b]; }
	NEW(rtn,N+1,sni_typ);
	NEW(tmpI,N+1,sni_typ);
	NEW(sum,N+1,Int4);
	NEW(E,N+1,e_type);
	NEW(tmpE,N+1,e_type);
        for(i=1 ;i<= N;i++){
           E1 = SeqScanInfo(I[i]);
           s1 = SiteScanInfo(1,I[i]);
	   s1 = MAXIMUM(Int4, s1 - Nt_flank, 1);
           s2 = SiteScanInfo(nblk,I[i]) + blklen[nblk] - 1;
	   s2 = MINIMUM(Int4, s2 + 10, LenSeq(E1));
	   E[i] = MkSubSeq(s1, s2, E1);
	}
	HG = Histogram("pairwise scores",0,500,5.0);
	/*** First cluster sequences into related sets ***/
	sets = DSets(N);
        for(i=1 ;i<= N;i++){
	   set1 = findDSets(i,sets);	
	   // fprintf(stderr,"score[%d][%d] = %d\n",i,j,score); 
           I1 = I[i]; E1 = SeqScanInfo(I1);
           for(j=i+1; j<= N;j++) {
                I2 = I[j]; E2 = SeqScanInfo(I2);
                for(ident=0,b=1; b <= nblk; b++){
                   s1 = SiteScanInfo(b,I1);
                   s2 = SiteScanInfo(b,I2);
           	   seq1 = XSeqPtr(E1) + s1; 
		   seq2 = XSeqPtr(E2) + s2;
		   for(k=0; k < blklen[b]; k++){
			if(seq1[k] == seq2[k]) ident++;
		   }
                }
		score = (Int4) ((double) (100.0 *ident)/(double)leng);
		IncdHist(score,HG);
		if(score >= cutoff){
                   set2 = findDSets(j,sets);
                   if(set1 != set2){ set1 = linkDSets(set1,set2,sets); }
		}
           }
        }
	PutHist(stderr,60,HG);
	NilHist(HG);
	/*** PutDSets(stderr,sets); /*****/
	/*** Then pick the most representative sequence from that set ***/
	for(num_seqs=s=0,i=1; i <= N; i++){
	  for(m=0,n=1; n <= N; n++){
	    E2=SeqScanInfo(I[n]);
	    s1 = findDSets(n,sets);	
	    if(i==s1){
		if(m==0){ s++; if(verbose) printf("set %3d:\n",s); }
		m++; tmpI[m] = I[n]; tmpE[m] = E[n];
		if(verbose) { printf("%3d: ",n); PutSeqInfo(stdout,E2); }
	    }	
	  }
	  if(m > 0) {  /*** set not empty ***/
	   if(verbose) printf("\n");
	   M = m; 
           for(m=0; m<= M;m++) sum[m] =0;
	   /******* Select a representative sequence *********/
           for(m=1; m<= M;m++) {
             I1 = tmpI[m]; E1 = SeqScanInfo(I1);
             for(j=m+1; j<= M;j++) {
                I2 = tmpI[j]; E2 = SeqScanInfo(I2);
                for(ident=0,b=1; b <= nblk; b++){
           	   seq1 = XSeqPtr(E1) + SiteScanInfo(b,I1); 
		   seq2 = XSeqPtr(E2) + SiteScanInfo(b,I2);
		   for(k=0; k < blklen[b]; k++){
			if(seq1[k] == seq2[k]) ident++;
		   }
                } sum[m] += ident; sum[j] += ident;
	      }
           }
           for(max=-1,m=1; m<= M;m++) {
		if(max < sum[m]) { max = sum[m]; j=m; }
	   }
           for(m=1; m<= M;m++) { if(m != j) NilSeq(tmpE[j]); }
	   num_seqs++; rtn[num_seqs] = tmpI[j];
	/**************************************************/
/******* Output the sequence *********/
           I1 = tmpI[j]; E1 = SeqScanInfo(I1);
	   PutSeq(stderr,tmpE[j],A);
/******* Output the sequence *********/
	  }
	}
        for(i=1 ;i<= N;i++) NilSeq(E[i]); free(E);
	NilDSets(sets); free(tmpI); free(sum);
	return rtnE;
}

e_type *PurgeScanInfo(Int4 cutoff, Int4 Nt_flank, Int4 Ct_flank, snh_typ sH, a_type A)
/********************************************************************************
  Purge a set of N sequences using the motif information to locate and compare
  similar regions.

  1. Cluster sequences at some cutoff level (based on substitution scores).
  2. Pick a representative sequence from each set.
	a. Purge cluster of sequences: take the one sequence that 
	   has the highest aggregate score with all other sequences?

/********************************************************************************/
{
        Int4    i,j,k,b,s1,leng,num_seqs,*sum,max;
	Int4	m,M,n,s,score,set1,set2;
        e_type  E1,E2,*E,*rtnE;
        unsigned char    *seq1,*seq2;
        sni_typ I1,I2,*tmpI,*I;
	gb_typ  B;
	h_type	HG;
	ds_type	sets;
	BooLean	verbose=TRUE;
	Int4	nblk,*blklen, N;

	I = DomainsScanHeap(&N, Nt_flank, Ct_flank, sH);
	NEW(E,N+2,e_type);
	for(i=1 ;i<= N;i++){ E[i] = SeqScanInfo(I[i]); }
	nblk = nblksScanHeap(sH);
	blklen = LengthsScanHeap(sH);
	for(leng=0,b=1 ;b<=nblk; b++) { leng+=blklen[b]; }
	NEW(sum,N+1,Int4);
	NEW(tmpI,N+1,sni_typ);
	NEW(rtnE,N+1,e_type);
	HG = Histogram("pairwise scores",0,500,5.0); /*****/
	/*** First cluster sequences into related sets ***/
	sets = DSets(N);
        for(i=1 ;i<= N;i++){
	   set1 = findDSets(i,sets);	
           E1 = E[i];
	   B = MakeGBlast(11, E1, A);
           for(j=i+1; j<= N;j++) {
           	E2 = E[j];
		score = MatcherGBlast(NULL,E2,B);
		// related = FastMatcherGBlast(E2, B, cutoff);
		// score = FastAlnSeqSW(12, 4, E1, E2, A);
		IncdHist(score,HG);
		if(score >= cutoff){
                   set2 = findDSets(j,sets);
                   if(set1 != set2){ set1 = linkDSets(set1,set2,sets); }
		}
           }
        }
	PutHist(stderr,60,HG);
	NilHist(HG);
	/*** PutDSets(stderr,sets); /*****/
	/*** Then pick the most representative sequence from that set ***/
	for(num_seqs=s=0,i=1; i <= N; i++){
	  for(m=0,n=1; n <= N; n++){
	    s1 = findDSets(n,sets);	
	    if(i==s1){
		if(m==0) s++; 
		m++; tmpI[m] = I[n];
	    }	
	  }
	  if(m > 0) {  /*** if set not empty then select a representative ***/
	   M = m; 
           for(m=0; m<= M;m++) sum[m] =0;
           for(m=1; m<= M;m++) {
             I1 = tmpI[m]; E1 = SeqScanInfo(I1);
             for(j=m+1; j<= M;j++) {
                I2 = tmpI[j]; E2 = SeqScanInfo(I2);
		/** score = FastAlnSeqSW(12, 4, E1, E2, A); /***/
                for(score=0,b=1; b <= nblk; b++){
           	   seq1 = XSeqPtr(E1) + SiteScanInfo(b,I1); 
		   seq2 = XSeqPtr(E2) + SiteScanInfo(b,I2);
		   for(k=0; k < blklen[b]; k++){
			if(seq1[k] == seq2[k]) score++;
		   }
                }
		sum[m] += score; sum[j] += score;
	      }
           }
           for(max=-1,m=1; m<= M;m++) { if(max < sum[m]) { max = sum[m]; j=m; } }
           // for(m=1; m<= M;m++) { if(m != j) NilScanInfo(tmpI[m]); }
	   num_seqs++; 
	   rtnE[num_seqs] = SeqScanInfo(tmpI[j]); 
	   // rtnE[num_seqs] = NilScanInfoRtnE(tmpI[j]); 
	  }
	}
	free(E); free(sum); free(tmpI); free(I); NilDSets(sets); 
	return rtnE;
}

Int4	RmIdent(e_type *List)
{ 
	Int4	i,j,N,N2;
	ss_type	P = NULL;
	e_type	E1,E2,*List2;

	for(N2=N=0; List[N+1] != NULL; ) N++;
	NEW(List2,N+2,e_type);
        for(i=1; i<= N; i++){
	   if(List[i] != NULL){
             E1 = List[i];
	     N2++; List2[N2] = E1;
             for(j=i+1; j<= N; j++){
		if(List[j] != NULL){
		   E2 = List[j];
		   if(IdentSeqs(E1,E2)){ NilSeq(E2); List[j] = NULL; }
		}
             }
	     List[i] = NULL;
	   }
        }
	for(i=1; i<= N2; i++){ List[i] = List2[i]; }
	free(List2);
	return N2;
}


