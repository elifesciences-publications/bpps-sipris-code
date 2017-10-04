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

#include "tblast.h"

tb_typ  MakeTBlast(char *name, char *dbs, Int4 maxhit, double q_value, 
	a_type A)
{
	tb_typ	T;
	Int4	i,*tmp,total;

	NEW(T,1,tblast_type);
	T->counts=NULL;

	T->number=GetFastaInfo(dbs, 0, &T->counts, NULL, A);
	NEW(T->L,T->number+3,BooLean);
	T->maxhit=maxhit; T->A=A; 
	T->dbs=AllocString(dbs); T->name=AllocString(name); 
	T->maxdepth = 0; T->qval = q_value; T->HG = NULL;

	/** compute dbsfreqs for karlin-altschul parameters **/
	if(T->counts==NULL){
	   NEW(T->counts, nAlpha(A) + 2, Int4);
	   for(i=0; i<= nAlpha(A); i++) T->counts[i] = tmp[i]; 
	} 
	for(total=0, i=0; i<= nAlpha(A); i++) total += T->counts[i];
	T->total = total;
	T->avelen = (double) total/(double) T->number;

	NEW(T->dbsfreq,nAlpha(A)+2,double);
        for(i=0;i<=nAlpha(A);i++)
		T->dbsfreq[i]=(double)T->counts[i]/(double)total;
	/** for(i=0;i<=nAlpha(A);i++) T->dbsfreq[i]=blastfreq[i]; /****/
	return T;
}

void    NilTBlast(tb_typ T)
{
	free(T->L); free(T->dbs); free(T->counts); free(T->dbsfreq);
	free(T->name);
	if(T->HG != NULL) NilHist(T->HG);
	free(T);
}

Int4	RunMultiTBlast(e_type *Q, double expect, BooLean lowH, Int4 depth, tb_typ T)
/** run multiple query sequences using tblast ***/
{
	e_type	E,*List;
	FILE	*fp;
	Int4	i,s,item,hits,len,n;
	double	key;
	char	*id;
	a_type	A250 = MkAlpha(AMINO_ACIDS,SSEARCH_PAM250);

	time_t	time1=time(NULL);
	T->maxdepth = depth; T->SH = MakeSqHeap(T->maxhit); 
	for(hits=0, s = 1; Q[s] != NULL; s++){
	   E = CopySeq(Q[s]); len = LenSeq(Q[s]);
	   /** printf("\n\t%d sequences (cutoff=%f)\n\n",
		T->number,expect/(float) T->number); /********/
	   hits += run_tblast(E, E, expect, lowH, 0,T); NilSeq(E);
	}
	fprintf(stderr,"\n\t%d hits\n\n",hits);

	NEW(List,hits+3,e_type);

	/*********** get hits ***********/
	fp = open_file(T->name,".blst","w");
	for(n=0; (E=DelMinSqHeap(&key,&item,T->SH))!=NULL; ){
		n++; List[n] = E;
		id = SeqKey(E);
		for(i=0; i<=55; i++){
		   if(id[i] == 0){ while(i++ <= 55) fprintf(fp,"."); break; }
		   else fprintf(fp,"%c",id[i]);
		} fprintf(fp,"   %g\n",key);
	}
	fclose(fp);
	/******* Save sequences *******
	/******************************************************************/
	fp = open_file(T->name,".seq","w");
	for(hits=0,i=1; i<=n; i++){
		E = List[i]; len = LenSeq(E);
		hits++; PutSeq(fp,E,T->A);
		NilSeq(E); 
	}
	fclose(fp); 
	free(List); NilSqHeap(T->SH); T->SH=NULL;
	NilAlpha(A250); 

	double runtime=difftime(time(NULL),time1);
	fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	return hits;
}

Int4	RunTBlast(e_type Q, double expect, BooLean lowH, Int4 depth, tb_typ T)
{
	e_type	E,*List;
	FILE	*fp,*fp2;
	Int4	i,item,hits,len;
	double	key,ratio;
	char	*id;
	Int4	*numaln,n,minlen;
	Int4	a=8,b=4,start,end;
	a_type	A250 = MkAlpha(AMINO_ACIDS,SSEARCH_PAM250);
	/************ NEW ***************
	h_type	HG,lenHG;
	/************ NEW ***************/

	time_t	time1=time(NULL);
	T->maxdepth = depth; T->SH = MakeSqHeap(T->maxhit); 
	E = CopySeq(Q); len = LenSeq(Q);
	/** printf("\n\t%d sequences (cutoff=%f)\n\n",
		T->number,expect/(float) T->number); /********/
	hits = run_tblast(E, E, expect, lowH, 0,T); NilSeq(E);
	fprintf(stderr,"\n\t%d hits\n\n",hits);

	NEW(numaln,len+3,Int4); NEW(List,hits+3,e_type);
	/************ NEW ***************
	HG = Histogram("SW aligned regions",0,len,10);
	lenHG = Histogram("sequence lengths",0,len,5);
	/************ NEW ***************/

	/*********** look at where hits align ***********/
	fp = open_file(T->name,".blst","w");
	for(n=0; (E=DelMinSqHeap(&key,&item,T->SH))!=NULL; ){

		SubSeqSW(a, b, Q, E, A250,&start,&end);
		for(i=start; i<=end; i++) numaln[i]++; 
		n++; List[n] = E;

		/************ NEW ***************
		for(i=start; i<=end; i++) IncdHist(i,HG); 
		IncdHist(LenSeq(E),lenHG);
		/************ NEW ***************/

		id = SeqKey(E);
		for(i=0; i<=55; i++){
		   if(id[i] == 0){ while(i++ <= 55) fprintf(fp,"."); break; }
		   else fprintf(fp,"%c",id[i]);
		} fprintf(fp,"   %g\n",key);
	}
	fclose(fp);
	/******************** Delineate minimum domain. *******************
	 Minimum domain is that region that region that aligns with at least
	 90% of the sequences.  
	/******************************************************************/
	for(start=end=0,i=1; i<= len; i++){
		ratio = (double) numaln[i]/(double)n;
		if(ratio >= 0.90){ if(start==0) start=i; end = i; } 
	}
	minlen = end - start + 1;
	fprintf(stderr,"\n\tmin domain: %d-%d; %d aa\n\n",start,end,minlen);
	/******* Save sequences at least as Int4 as minimum domain. *******
	 Assume that sequences shorter than minimum domain are fragments.
	/******************************************************************/
	fp = open_file(T->name,".seq","w");
	fp2 = open_file(T->name,".fsq","w");
	for(hits=0,i=1; i<=n; i++){
		E = List[i]; len = LenSeq(E);
		if(len >= minlen){ hits++; PutSeq(fp,E,T->A); }
		else {
			PutSeq(fp2,E,T->A);
			fprintf(stderr,"sequence %d: too short (%d)\n",i,len);
		}
		NilSeq(E); 
	}
	fclose(fp); fclose(fp2);
	free(numaln); free(List); NilSqHeap(T->SH); T->SH=NULL;
	NilAlpha(A250); 

	/************ NEW ***************
	PutHist(stderr,60,HG); NilHist(HG); 
	PutHist(stderr,60,lenHG); NilHist(lenHG); 
	/************ NEW ***************/

	double runtime=difftime(time(NULL),time1);
	fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	return hits;
}

Int4	run_tblast(e_type Q, e_type fullQ, double expect,BooLean lowH,
	Int4 depth,tb_typ T)
{ 
	Int4	v,i,j,hits,len,lenQ,flenQ,cut;
	Int4	s,flank=20,low,high,cutoff=200,numseq;
	double	*pr,*freq;
	double	pval,N,p,key,adjust;
	char	*id,str[200];
	unsigned char	*seq;
	FILE	*fptr=NULL,*tfp=NULL;
        e_type  E,fullE,*FullQ,*ListE,*ListQ;
	sh_typ	SH;
	a_type	A=T->A,A250;
	gb_typ  B,fB;
	Int4	overlap,item;
	Int4	fstart,fend,flen,fscore;
	Int4	start,end,slen,score;
	double	percent;
	BooLean	flag,*Use;
	Int4	a=8,b=4;	

	A250 = MkAlpha(AMINO_ACIDS,SSEARCH_PAM250);

	if(depth < T->maxdepth){
		NEW(FullQ,T->maxhit+2,e_type); SH = MakeSqHeap(T->maxhit); 
	} else SH = NULL;

	/***** A. Compute karlin-altschul parameters *****/
        low=lowAlphaR(A); high=highAlphaR(A); len=high - low + 1;
        NEW(pr,len+1,double);
        NEW(freq,nAlpha(A)+2,double);	/** query sequence frequency **/
	FreqResSeq(Q, freq, A);
        for(i=0; i<=nAlpha(A); i++){
           for(j=0; j<=nAlpha(A); j++){
                v=(Int4)valAlphaR(i,j,A)-low; pr[v]+=T->dbsfreq[i]*freq[j];
           }
        }
        if(!karlin(low,high,pr,&T->lambda,&T->K,&T->H)){print_error("fatal");}
	free(pr); free(freq);
	/***** end compute karlin-altschul parameters *****/

	/** B. Compute threshold score: from Stephen Altschul ******/
	/** T->T = (Int4) ceil((4.676 + 1.3*log(T->H))/T->lambda - 0.5); /***/
	T->T = (Int4) floor((4.676 + 1.3*log(T->H))/T->lambda - 0.5);
        fprintf(stderr,"K=%g; lambda=%g; H=%g; T=%d; database %d residues\n", 
		T->K,T->lambda,T->H,T->T,T->total);
	/****** end Compute threshold score ******/

	/***** C. Open database file. ******/
	if(depth == 0){
          tfp = open_file(T->name,".tmp","w");
          if((fptr = fopen(T->dbs,"r")) == NULL) {
             fprintf(stderr,"Could not open file \"%s\"\n",T->dbs);
             print_error("File does not exist!");
          }
	} else fptr = open_file(T->name,".tmp","r");
	if(setvbuf(fptr,NULL,_IOFBF, BUFSIZ*200))
                print_error("tblast() buffer error");
	/***** end Open database file. ******/

	/***** D. Prepare query sequence for search. *****/
	if(lowH) { 
	   /** Note: seqset.c does ProcessSeqPSeg(17, 2.2,2.5,100,E,A); **/
	   ProcessSeqPSeg(12,2.2,2.5,100,Q,A); 
	   // ProcessSeqPSeg(45,3.4,3.75,100,Q,A); /** mask coiled-coils: **/
	   if(depth > 0){
		ProcessSeqPSeg(12,2.2,2.5,100,fullQ,A); 
		// ProcessSeqPSeg(45,3.4,3.75,100,fullQ,A); /** coiled-coils: **/
	   }
	}
	/***** end Prepare query sequence for search. *****/

	/***** E. Prepare finite state machines, etc. for BLAST search. *****/
	lenQ = LenSeq(Q);
	cut=(Int4)((log((double)(lenQ*T->avelen))+log(T->K/0.5))/T->lambda);
	/***/ fprintf(stderr,"cutoff = %d\n",cut); /****/
	B = MkGBlast(50, cut, T->T, Q, A);

        if(depth > 0){	/** then make full length FSM to check hits **/
	  flenQ = LenSeq(fullQ);
	  cut=(Int4)((log((double)(flenQ*T->avelen))+log(T->K/0.5))/T->lambda);
	  /***/ fprintf(stderr,"fullcutoff = %d\n",cut); /****/
	  fB = MkGBlast(50, cut, T->T, fullQ, A);
        } else { fB=NULL; flenQ=0; }

	NEW(id,MAX_ID_LENG_GB,char); NEW(seq,MAX_SEQ_LENG_GB,unsigned char);
	/***** end Prepare for BLAST search. *****/

	/***** F. Do BLAST search. *****/
        for(hits=0,i=1; i<=T->number;i++) {
	    len = ReadSeqStr(fptr,id,seq,MAX_SEQ_LENG_GB,A);
	    if(FALSE && i%1000==0) fprintf(stderr,"\r%d",i);
	    if(T->L[i] < 1){
	      score = MatcherGBlastStr(len, seq, B, AlphaR(A));
	      N = len*lenQ;
	      pval=SumStatGBlastStr(lenQ,len,T->lambda,T->H,T->K,B);
	      adjust = (double)T->total/(double)len;
			// old: adjust = (double)T->number; 
	      if((pval*adjust) <= expect){
		 T->L[i]=1;
             	 /*** fprintf(stdout,"\rpval=%g; i=%d\n",pval,i); /***/
		 /** give low complexity junk only one chance **/
 		 fullE = MkSeq(id, len, seq); EqSeqI(0,fullE);
		 /** don't want high scores due to low complexity junk **/
		 ProcessSeqPSeg(12,2.2,2.5,100,fullE,A);
		 score = MatcherGBlastStr(len, XSeqPtr(fullE), B, AlphaR(A));
		 p = SumStatGBlastStr(lenQ,len,T->lambda,T->H,T->K,B); 
		 /** Allow at most a 5-fold increase with SEG'ed sequence. **/
		 if((p*adjust) <= (5.0*expect)){ /** allow 5-fold seg increase **/
		    if(depth > 0){   /** then check p-value for full query ***/
			score = MatcherGBlastStr(len, seq, fB, AlphaR(A));
			p = SumStatGBlastStr(flenQ,len,T->lambda,T->H,T->K,fB); 
			if((p*adjust) <= (5.0*expect)) flag = TRUE;
		    	else flag = FALSE;
		    } else flag = TRUE;
		 } else flag = FALSE;
		 if(flag){
		   hits++;
		   /******* don't use seg'ed sequence ********/
	           score = MatcherGBlastStr(len, seq, B, AlphaR(A));
		   /******* ^don't use seg'ed sequence ********/
	  	   pval *= adjust; 
#if 0		   //************************************************ 
		   fprintf(stderr,"\n%d: pval=%g; p=%g; raw p = %g\n",
			i,pval,p*adjust,pval/adjust);
	      	   e = -expm1(-T->K*N*exp(-T->lambda*(double)score));
	  	   e *= adjust;
		   fprintf(stderr," (S=%d; E=%1.2g: %d aa; SumP = %g) ", 
			score,e,len,pval);
		   for(k=0; id[k]!=0 && !isspace(id[k]); k++){
                	fprintf(stderr,"%c", id[k]);
           	   } fprintf(stderr,"\n");
#endif		   //************************************************ 
		   if(depth < T->maxdepth && score <= cutoff){
		     AlnSeqSW(a, b, fullE, Q, A250);
		     score=SubSeqSW(a, b, fullE, Q, A250,&start,&end);
		     slen = end - start +1;
		     if(depth > 0){	/** if down at least one level **/
		        /***** TRY FULL SEQUENCE IF TRANSITIVE TEST ***/
			// AlnSeqSW(a, b, fullE, fullQ, A250); 
			fscore=SubSeqSW(a,b,fullE,fullQ,A250,&fstart,&fend);
			flen=fend-fstart+1;
			/** A. See if at least 50% overlap with subQ **/
			if(start > fstart){
			   if(end > fend) overlap=fend-start+1;
			   else overlap=slen;
			} else {
			   if(end > fend) overlap=flen;
			   else overlap=end-fstart+1; 
			}
			percent = (double) overlap/(double)slen;
			if(percent > 0.50){	/** then use fullQ **/
			   E=MkSeq(SeqKey(fullE),(fend-fstart+1),(seq+fstart-1));
			} else {		/** use subQ **/
			   E=MkSeq(SeqKey(fullE),(end-start+1),(seq+start-1));
			}
		        /***** TRY FULL SEQUENCE IF TRANSITIVE TEST ***/
		     } else E=MkSeq(SeqKey(fullE),(end-start+1),(seq+start-1));
		     // i = MINIMUM(Int4,n1,max_i+flank);
		     // j = MAXIMUM(Int4, 1,i-flank+1);
		     // PutSeq(stderr,E,A); 
		     score = LenSeq(E); 
		     if(score >= 50) { 	/** don't use short sequences **/
		     	/** use high scores first **/
		     	item=InsertSqHeap(E,-score,SH); 
		        if(item==NULL) NilSeq(E);
			else {
			   if(FullQ[item] != NULL) NilSeq(FullQ[item]);
			   FullQ[item] = CopySeq(fullE); 
			}
		     } else { NilSeq(E); }
		     if(InsertSqHeap(fullE,pval,T->SH)==NULL) NilSeq(fullE); 
		   } else {
		     if(InsertSqHeap(fullE,pval,T->SH)==NULL) NilSeq(fullE);
		   }
		 } else {	/** low complexity junk **/
             		/*** fprintf(stdout,"pval=%g; p=%g\n",pval,p);
			PutSeq(stdout,fullE,A); /****/
			NilSeq(fullE);
		 }
	         if(depth == 0){ fprintf(tfp,">null \nX\n"); }
	      } else if(depth == 0){	/** first level of search **/
		   if(pval > T->qval){ 
			fprintf(tfp,">null \nX\n");
			T->L[i]=1; 
		   } else { /** save in a file **/
 		 	E = MkSeq(id, len, seq);
		        PutSeq(tfp,E,A);
		        NilSeq(E);
		   }
	      }
	    }
	}
	if(depth == 0) {
	  fclose(tfp);
          tfp = open_file(T->name,".hits","w");
	  fprintf(tfp,"\n\t%d out of %d sequences found (depth=%d)\n",
		hits,T->number,depth); fclose(tfp);
	}
	NilAlpha(A250); free(id); free(seq);
	fprintf(stdout,"\n\t%d out of %d sequences found (depth=%d)\n",
		hits,T->number,depth);
	fprintf(stdout,"query: "); PutSeqID(stdout,Q); fprintf(stdout,"\n");
	fclose(fptr); NilGBlast(B);
	if(fB!=NULL) NilGBlast(fB);
	/***** G. Do transitive BLAST search at lower levels. *****/
	if(depth < T->maxdepth){
	  numseq=nSqHeap(SH);
	  if(numseq > 0){
	    /*************************************************************
	    Start with weakest significant hits and look for transitive
	    hits.  Purge list to reduce redundancy.
	    /*************************************************************/
	    NEW(ListE,numseq+3,e_type); NEW(ListQ,numseq+3,e_type);
	    for(s=1; (E=DelMinSqHeap(&key,&item,SH))!=NULL; ){ 
		ListE[s] = E;
		ListQ[s] = FullQ[item]; 
		s++; 
	    }
	    NilSqHeap(SH);
	    Use=RmHomologsList(cutoff,'b',1,T->maxhit,
					FALSE,T->T,numseq,ListQ,A);
	    fprintf(stderr,"%d in ListQ\n",numseq);
	    for(s=1; s <= numseq; s++){
	      E = ListE[s];
	      if(Use[s]) hits += run_tblast(E,ListQ[s],expect,lowH,depth+1,T);
	      NilSeq(E); NilSeq(ListQ[s]);
	    }
	    free(Use);
	    free(ListQ); free(ListE); 
	  }
	  free(FullQ); 
	}
	if(depth == 0) {
		strcpy(str,T->name); strcat(str,".tmp"); remove(str);
	}
	return hits;
}

