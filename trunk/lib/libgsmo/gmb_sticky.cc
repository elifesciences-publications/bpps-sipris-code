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

#include "gmb_typ.h"

//************************ Single Sequence Sampling Routines ***********************
BooLean	gmb_typ::SampleTogether(set_typ Set,double Temp,double &MAP, dh_type dH)
// sample insertions and deletions in sequence sq.
{
	cma_typ CMA=SSX->RtnCMA(); 
	Int4	sq,i,start,score,*newpos,**oldpos,NN=NumSeqsCMSA(CMA);
	char	*operation;
	gsq_typ	*ogsq,*gsq;
	e_type	sbjE=0;

	MAP=SSX->GapMap(Temp);
	//======== 1. Remove sq from alignment, sample scoring matrix and realign. =====
	NEWP(oldpos,NN+5,Int4);
        for(sq=1; sq <= NN; sq++){
	    if(!MemberSet(sq,Set)) continue;
	    oldpos[sq]=SSX->RemoveFromAlign(sq); assert(oldpos[sq] != 0);
	} gss_typ *gss=gssCMSA(CMA); SSX->InitNDL(Temp);  // initialize model parameters?

	// update model after each addition....
	//========== 2. Create a gapped sequence using 'operation'. ===========
	while(!emptyHeap(dH)){
	   sq=delminHeap(dH); assert(sq <= NN);
	   // if(!MemberSet(sq,Set)) continue;
	   assert(MemberSet(sq,Set)); DeleteSet(sq,Set);
	   sbjE = gss->TrueSeq(sq); 	
	   operation=SSX->GapAlnTrace(sbjE,Temp,start,score);  // Recalculates HMM parameters...
	   Int4 trace_length=strlen(operation);
	   gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(CMA)+3,Int4);  
           gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,sbjE,newpos);
	   // fprintf(stderr,"%d: %s\n",sq,operation);
	   //========== 3. If qsq is unchanged then retain old alignment and return. =========
	   BooLean the_same=FALSE;
	   if(gss->Identical(sq,*gsq)){
	        the_same=TRUE;
           	for(Int4 b=1; b<= nBlksCMSA(CMA); b++)
              	   { if(oldpos[sq][b] != newpos[b]){ the_same=FALSE; break; } }
	   }
           if(the_same){ SSX->AddToAlign(sq,oldpos[sq]); delete []gsq; }
	   else { ogsq=SwapGsqCMSA(sq,gsq,CMA); SSX->AddToAlign(sq,newpos); delete []ogsq; }
	   free(oldpos[sq]); free(newpos); free(operation); 
	} free(oldpos); assert(CardSet(Set) == 0);
	SSX->InitNDL(Temp);	// initialize model parameters?
        MAP=SSX->GapMap(Temp); 
	return TRUE;
}

/********************************* private ********************************/
/*****************************************************************
  Return the likelihood of site at pos in seq being in model M
  WARNING: variables s and likelihood are assumed to have private
        information:
        i.e., s = M->start - M->end.
        i.e., likelihood = M->likelihood + M->start.
/*****************************************************************/
double  log_likelihood_fmodel(register unsigned char *seq, register Int4 s,
        register double **likelihood)
{
        register double L=0.0;
	FILE	*efp=0; 
	Int4	S;

        for(S=s; s >= 0; s--){ if(likelihood[s]!=NULL) L +=log10(likelihood[s][(seq[s])]); }
        if(!isfinite(L)){
        	for(L=0.0,s=S; s >= 0; s--){
                  if(likelihood[s]!=NULL) L += log10(likelihood[s][(seq[s])]);
		  if(efp) fprintf(stderr,"likelihood[%d][%d] = %g; sum=%g.\n",
			s,seq[s],likelihood[s][(seq[s])],L);
		} assert(isfinite(L));
	} return L;
}

#define LogLikelihoodFModel(s,p,M) ( ((M)->update) ? update_fmodel(M), \
        log_likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)):\
        log_likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)))

double	GetProbCMSA_X(Int4 t, Int4 n, cma_typ L)
{
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    fm_type		*model=ModelsCMSA(L);
    Int4		s,end;
    double		prob;
    unsigned char	*seq;
    mdl_typ		*mdl=mdlCMSA(L);

    if(n < 1 || n > NSeqsSeqSet(data) || t < 1 || t > nTypeSites(sites))
			print_error("GetProbCMSA( ) input error");
    end = SqLenSeqSet(n,data) + 1;
    end += 1 - LenFModel(model[t]);
    seq = SeqSeqSet(n,data);
    s = SitePos(t,n,1,sites);
    mdl->Remove(t,seq,s);
    prob = (double) LogLikelihoodFModel(seq, s, model[t]);
    mdl->Add(t,seq,s);
#if 0
    if(prob > 0.0) prob=log10(prob);
    else {
	fprintf(stderr,"Seq %d: prob = %g!\n",n,prob);
	// PutSeq(stderr,TrueSeqCMSA(n,L),AlphabetCMSA(L));
	prob=(double) INT4_MIN;
	// assert(prob > 0.0);
    }
    if(!isfinite(prob)){
	mdl->Remove(t,seq,s);
	PutFModel(stderr,model[t]);
	PutSeq(stderr,TrueSeqCMSA(n,L),AlphabetCMSA(L));
	assert(isfinite(prob));
    }
#else
    assert(isfinite(prob));
#endif
    return prob;
}

double  GetTotalProbCMSA_X(Int4 n, cma_typ cma)
{
	assert(n > 0 && n <= NumSeqsCMSA(cma));
        assert(SeqIsInCMSA(n,cma));
	double	prob=0.0;
	for(Int4 b=1; b<=nBlksCMSA(cma); b++){ prob += GetProbCMSA_X(b,n,cma); }
	return prob;
}

double	gmb_typ::RandomTogether(Int4 MinSticky, double MaxFrctn, double temp)
// sample with gaps using simulated annealing...
{
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
        assert(nBlksCMSA(CMA) == 1);
	Int4    blk=1;
	Int4	stop,n,i,j,r,s,t,sq,NN=NumSeqsCMSA(CMA);
        double  d,D,best,map;
	FILE	*efp=0;
	Int4	nsets=10;	// create 10 partitions...
	nsets=(Int4) ceil(1.0/MaxFrctn);
	assert(MaxFrctn >= 0.001 && MaxFrctn <= 0.50);	// 2 to 1000 paritions allowed...
	
	SetPseudoToMapCMSA(CMA); 
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	best=SSX->GapMap(temp); 
	// fprintf(stderr,"best map = %.2f; Temp = %.1f\n",best,temp);

   set_typ  *set;  NEW(set,nsets+5,set_typ); 
   for(s=1; s <= nsets; s++){ set[s]=MakeSet(NN+5); ClearSet(set[s]); }
   for(sq=1; sq <= NumSeqsCMSA(CMA); sq++) {
#if 1
	if(SetEx && MemberSet(sq,SetEx)) continue;
	else {
	   s=random_integer(nsets)+1; AddSet(sq,set[s]);	// 0 < s <= nsets.
	}
#else
	s=random_integer(nsets)+1; AddSet(sq,set[s]);	// 0 < s <= nsets.
#endif
   }
   for(n=1; n <= nsets; n++){

	dh_type dH=dheap(NN+5,4);
        for(sq=1; sq <= NumSeqsCMSA(CMA); sq++) {
	    if(!SeqIsInCMSA(sq,CMA)) continue;
	    if(!MemberSet(sq,set[n])) continue;
	    insrtHeap(sq,(keytyp) Random(),dH); 
	} 
	dh_type dH2=dheap(NN+5,4);
	for(i=1; !emptyHeap(dH); i++){
		keytyp key=minkeyHeap(dH); sq=delminHeap(dH);
		// fprintf(stderr,"%d: sq %d score = %.3f\n",i,sq,key);
		insrtHeap(sq,-key,dH2);	  // put in random order...
	} Nildheap(dH); dH=dH2; // best=0.0;
	if(efp) fprintf(stderr,"Sampling %d seqs together.\n",stop);
	if(!emptyHeap(dH) && SampleTogether(set[n],temp,map,dH)){
            if(map > best){ 
		if(efp) fprintf(stderr,"map changes from %.1f to %.1f (%.1f K)\n",best,map,temp);
         	best=map; SaveBestCMSA(SSX->RtnCMA());
            }
	} // fprintf(stderr,"final temperature = %.1f K\n %.1f\n", temp,best);
	NilSet(set[n]); Nildheap(dH);
   } free(set);
	double map0=SSX->GapMap(temp);
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	map=SSX->GapMap(temp);
	if(efp) fprintf(stderr,"map1 = %.1f; restored map = %.1f; best = %.1f\n",map0,map,best);
	return best;
}

double	gmb_typ::WorstTogether(Int4 MinSticky, double MaxFrctn, double temp)
// Find the most poorly aligned sequences and resample them ...
{
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
        assert(nBlksCMSA(CMA) == 1);
	Int4    blk=1;
	Int4	stop,n,i,j,r,s,t,sq,NN=NumSeqsCMSA(CMA);
	FILE	*efp=0;

	stop=(Int4) floor(MaxFrctn*(double)NumAlnSeqsCMSA(CMA));
        if(stop < MinSticky) stop=MinSticky; 

	SetPseudoToMapCMSA(CMA); 
	SSX->RestoreBest(); // reset to optimum CMA configuration.
        double  d,D,lpr,best,lpr0; lpr0=best=SSX->GapMap(temp); 
	// fprintf(stderr,"best lpr = %.2f; Temp = %.1f\n",best,temp);

	// pull out the worst based on scores; sample back in from better to worse.
	dh_type dH=dheap(NN+5,4);
        for(sq=1; sq <= NumSeqsCMSA(CMA); sq++) {
	    if(!SeqIsInCMSA(sq,CMA)) continue;
#if 0	// This is working fine for tweakcma...
	   d=GetGappedProbCMSA(1,sq,CMA); if(!isfinite(d)) continue;
#elif 1	// This takes nearly twice as long but avoids overflow errors...
	   d=GetTotalProbCMSA_X(sq,CMA); assert(isfinite(d));
#else
	   d=GetGappedProbCMSA(1,sq,CMA); 
	   if(!isfinite(d)){
		if(isinf(d) == 1) d=DBL_MAX;
		else if(isinf(d) == -1) d=(double) LONG_MIN;
		else assert(isfinite(d));
	   } insrtHeap(sq,(keytyp)d,dH);
#endif
	} 
	dh_type dH2=dheap(NN+5,4);
	set_typ	Set=MakeSet(NN+5);  ClearSet(Set);
	for(i=1; !emptyHeap(dH); i++){
		keytyp key=minkeyHeap(dH); sq=delminHeap(dH);
		// fprintf(stderr,"%d: sq %d score = %.3f\n",i,sq,key);
		if(i <= stop){ AddSet(sq,Set); insrtHeap(sq,-key,dH2); } else break;
	} Nildheap(dH); dH=dH2; // best=0.0;
	if(efp) fprintf(stderr,"'W': Sampling %d seqs together.\n",stop);
	if(SampleTogether(Set,temp,lpr,dH)){
            if(lpr > best){ 
		if(efp) fprintf(stderr,"lpr changes from %.1f to %.1f (%.1f K)\n",best,lpr,temp);
         	best=lpr; SaveBestCMSA(SSX->RtnCMA());
            }
	} // fprintf(stderr,"final temperature = %.1f K\n %.1f\n", temp,best);
	NilSet(Set); Nildheap(dH); 
	SSX->RestoreBest(); lpr=SSX->GapMap(temp);
	if(efp) fprintf(stderr,"lpr1 = %.1f; restored lpr = %.1f; best = %.1f\n",lpr0,lpr,best);
	return best;
}

set_typ	gmb_typ::GambitTogether(Int4 II, Int4 MinSticky, double MaxFrctn, double temp,
		set_typ LstSet, Int4 min_ins)
// sample with gaps using simulated annealing...
{
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
        assert(nBlksCMSA(CMA) == 1);
	Int4    blk=1,n_ins,n_del,ins_res;
	Int4	n,i,j,r,s,t,x,NN=NumSeqsCMSA(CMA),sq;
        double  d,D,best,map;
	set_typ	Set=0;
	dh_type dH=0;
	FILE 	*efp=0; // efp=stderr;

	//===================  1. Setting up... ============================
	SetPseudoToMapCMSA(CMA); 
	SSX->RestoreBest(); // reset to optimum CMA configuration; retains this as best...
	best=SSX->GapMap(temp);
	// fprintf(stderr,"best map = %.2f; Temp = %.1f\n",best,temp);
        assert(II > 0 && II <= LengthCMSA(blk,CMA));

    if(min_ins > 0){	
	//===================  2a. Look for insertions... ============================
	Set=MakeSet(NN+5);  ClearSet(Set); dH=dheap(NN+5,4);
        for(n_ins=ins_res=0,sq=1; sq <= NumSeqsCMSA(CMA); sq++){
		if(!SeqIsInCMSA(sq,CMA)) continue;
                n = InsertionCMSA(blk,sq,II,CMA);
                if(n >= min_ins){
		   AddSet(sq,Set); n_ins++; ins_res += n; 
		   insrtHeap(sq,(keytyp)n,dH);	// sample seqs with shortest inserts first.
		}
        } d=(double)n_ins/(double)NumAlnSeqsCMSA(CMA);
        if(n_ins < MinSticky || d > MaxFrctn){ NilSet(Set); Nildheap(dH); return 0; }
	if(efp) fprintf(stderr,"%d: n_ins = %d; FractIns = %.3f\n",II,n_ins,d);
     } else {		
	//===================  2a. Look for deletions... ============================
	// n_del=NumDeletionsCMSA(blk,II,CMA); d=(double)n_del/(double)NN;
	n_del=NumDeletionsCMSA(blk,II,CMA); d=(double)n_del/(double)NumAlnSeqsCMSA(CMA);
        if(n_del < MinSticky || d > MaxFrctn) return 0; 
	Set=MakeSet(NN+5);  ClearSet(Set); dH=dheap(NN+5,4);
        for(n=0,sq=1; sq <= NN; sq++){
	   if(!SeqIsInCMSA(sq,CMA)) continue;
	   if(IsDeletedCMSA(blk,sq,II,CMA)){ n++; AddSet(sq,Set); insrtHeap(sq,(keytyp) Random(),dH); }
	}
	if(LstSet){
	    d = (double) CardInterSet(Set,LstSet)/(double) n;
	    if(d > 0.75){ NilSet(Set); Nildheap(dH); return 0; }
        } if(efp) fprintf(stderr,"%d: n_del = %d; FractDel = %.3f\n",II,n_del,d);
#if 0
FILE *wfp=open_file("junk_debug",".cma","w"); PutInSetCMSA(wfp,Set,CMA); fclose(wfp);
wfp=open_file("junk_debug2",".cma","w"); PutAlnSeqsCMSA(wfp,CMA); fclose(wfp);
#endif
     }
	if(SampleTogether(Set,temp,map,dH)){
            if(map > best){ 
		if(efp) fprintf(stderr,"%d: map changes from %.1f to %.1f (%.1f K)\n",II,best,map,temp);
         	best=map; SaveBestCMSA(SSX->RtnCMA());
            }
	} // fprintf(stderr,"final temperature = %.1f K\n %.1f\n", temp,best);
	Nildheap(dH); // NilSet(Set); 
	double map0=SSX->GapMap(temp);
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	map=SSX->GapMap(temp);
	if(efp) fprintf(stderr,"map1 = %.1f; restored map = %.1f; best = %.1f\n",map0,map,best);
	return Set;
}

#if 0	// adapt this to finding the number of deletions within each sequence.
Int4	PutUniqueMergedMinColCMSA(char *filename, char *matstr, double MinCol, a_type AB)
{
	Int4 Number,file;
	FILE *fp=OpenFileToRead(filename);
	cma_typ cma,tcma,*IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);

	// -m option
	assert(nBlksCMSA(IN_CMA[1]) == 1);
	if(Number > 1){
	   fp=tmpfile();
	   PutMergedCMSA(fp,Number,IN_CMA); 
	   for(file=1; file <= Number; file++) TotalNilCMSA(IN_CMA[file]); free(IN_CMA);
	   rewind(fp); cma=ReadCMSA(fp,AB); fclose(fp);
	} else { cma=IN_CMA[1]; free(IN_CMA); }

	// -mincol=0.75 option
	h_type HG=Histogram("fraction aligned",0,1,0.025);
	Int4	i,j,k,s,J,I,n,Len,na;
	Int4    N = NumSeqsCMSA(cma);
        BooLean *skip; NEW(skip,N+3,BooLean);
        for(J=1; J <= N; J++){ skip[J]=TRUE; }
	Len=LengthCMSA(1,cma);
        for(n=0,J=1; J <= N; J++){
		// e_type E=FakeSeqCMSA(J,cma);
		for(na=0,s=1 ; s <= Len;s++){
			Int4 r=ResidueCMSA(1,J,s,cma);
			if(r != UndefAlpha(A)) na++; 
		}
		double fr=(double)na/(double)Len;
		if(fr >= MinCol){ skip[J]=FALSE; n++; }
		IncdHist(fr, HG);
	} 
	if(n > 0){
		fp=tmpfile(); PutSelectCMSA(fp,skip,cma); 
		TotalNilCMSA(cma); rewind(fp); cma=ReadCMSA(fp,AB); fclose(fp);
		PutHist(stdout,60,HG); NilHist(HG); fflush(stdout); free(skip);
	} else {
		PutHist(stdout,60,HG); NilHist(HG); fflush(stdout);
		free(skip); TotalNilCMSA(cma);
		return 0;
	}

	// -U option...
	ss_type data = TrueDataCMSA(cma);
	N=NSeqsSeqSet(data); NEW(skip,N+3,BooLean); 
        for(i=1;i < N; i++) {
	   if(skip[i]) continue;
	   e_type  qE=SeqSetE(i,data);
	   if(i % 1000 == 0) fprintf(stderr,"\r%.1f",100.0*((double)i/(double)N));
       	   for(j=i+1;j <= N; j++) {
		if(skip[j]) continue;
		if(IdentSeqs(qE,SeqSetE(j,data))) skip[j]=TRUE;
	   }
	}
	fp = open_file(filename,matstr,"w");
	PutSelectCMSA(fp,skip,cma); free(skip); fclose(fp);
	TotalNilCMSA(cma);
	return n;
}
#endif

BooLean	gmb_typ::SamplePurged(Int4 percent_id, double MaxFrctn,double Temp, double &strt_map, 
	double &end_map, Int4 stage)
// This doesn't appear to help; not sure why...
{
        cma_typ bcma=SSX->RtnCMA(),cma=CopyCMSA(bcma);
        Int4    x,numSets,sq,i,j,k,N=NumSeqsCMSA(cma),*oldpos;
        double  d,map,best,bst_lpr;
	dh_type	dH=0,pdH=0;
	FILE 	*efp=0; // efp=stderr;

        SaveBestCMSA(cma); assert(MaxFrctn <= 0.95);

	//=========== 1. Find related sequence sets... ==============
        set_typ *SetSq=RtnTargetSizeSetsCMSA(numSets,percent_id,cma,MaxFrctn);
	// ^ this resets percent_id so that largest set <= MaxFractn of seqs.
	x=CardSet(SetSq[1]);
	d=((double) x/(double)N);	// maximum fraction of seqs sampled together...
	if(efp) fprintf(stderr," max size=%d/%d=%.3f; %d%c identity.\n",x,N,d,percent_id,'%');

	//======== 2. Select one sequence from each non-singleton set. =============
	set_typ Used=MakeSet(SetN(SetSq[1])),Sample=MakeSet(SetN(SetSq[1]));
	ClearSet(Used); ClearSet(Sample); pdH = dheap(N+4,3);	// purged sequence heap.
        for(i=1; i <= numSets; i++){
            x=CardSet(SetSq[i]);  assert(x > 0);
	    if(efp && x >= 2) fprintf(efp,"%d ",x);
	    if(x < 2) UnionSet(Used,SetSq[i]); // Used = Used U SetSq[i].
	    else {	// pick one sequence at random to keep in the alignment
		k=random_integer(x) + 1;	// 1 <= k <= x.
                for(j=0,sq=1; sq <= N; sq++){
                   if(MemberSet(sq,SetSq[i])){
		     j++;   // d=GetProbCMSA(1,sq,cma); if(!isfinite(d)) d=0;
		     if(j == k){	// use only one sequence from each set.
			AddSet(sq,Used); AddSet(sq,Sample); DeleteSet(sq,SetSq[i]); 
		     } else {
		        d=SampleUniformProb(); assert(d >= 0 && d <= 1.0);
		        d += (double) x;	// add back the largest sets last...
		        insrtHeap(sq,(keytyp)d,pdH);
		     }
		   }
            	} assert(j == x && k > 0 && k <= x);
	    }
        } if(efp) fprintf(efp," (%d sets)\n",numSets); 
	x=CardSet(Used) + ItemsInHeap(pdH); assert(x == N);

	//========= 3. create a new gmb_typ.
	Int4	pio,pdo,eie,ede;
	double	pwt;
	char	dms_md;

	SSX->GetParameters(pio,pdo,eie,ede,pwt,dms_md); // this->Free(); NilCMSA(cma); 
	gmb_typ *gmb = new gmb_typ(pio,pdo,eie,ede,cma,dms_md,pwt,0); 
	gmb->RestoreBest(); bst_lpr=gmb->RtnMap(); 
	ssx_typ	*ssx=gmb->RtnSSX();

	//========= 3. Store rep. seqs on a heap; remove other seqs from alignment. ==========
	// The 'Sample' set contains the representative sequences for the larger sets only.
	for(dH=dheap(N+4,3),sq=1; sq <= N; sq++){
	     if(MemberSet(sq,Used)){
		if(stage == 0) insrtHeap(sq,((keytyp)Random()),dH);
		else if(MemberSet(sq,Sample)) insrtHeap(sq,((keytyp)Random()),dH);
		assert(!memHeap(sq,pdH));
	     } else { 
		assert(memHeap(sq,pdH));
		oldpos=ssx->RemoveFromAlign(sq); free(oldpos); 
	     }
	} SaveBestCMSA(cma);

	//=============== 4. Sample remaining sequences one-at-a-time. ==============
	best=gmb->RtnMap();
	for(j=1; (sq=delminHeap(dH)) != 0; j++){
	    assert(MemberSet(sq,Used));
	    if(stage != 0) assert(MemberSet(sq,Sample));
	    if(gmb->SampleSingle(sq,Temp,map)){
            	if(map > best){ best=map; SaveBestCMSA(cma); }
	    }
	} Nildheap(dH); gmb->RestoreBest(); // reset to optimum CMA configuration.

	//========== 5. Sample the removed sequences back in one-at-a-time. ===============
	Int4	start=0,score,trace_length,*newpos;
	gss_typ *gss=gssCMSA(cma);
	for(j=1; (sq=delminHeap(pdH)) != 0; j++){
	    assert(!MemberSet(sq,Used));
	    assert(!SeqIsInCMSA(sq,cma));	// i.e., sq has no aligned sites.
	    e_type  sbjE=gss->TrueSeq(sq); ssx->InitNDL(Temp);     
	    char    *operation=ssx->GapAlnTrace(sbjE,Temp,start,score);
	    trace_length=strlen(operation);
	    gsq_typ *ogsq,*gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(cma)+3,Int4);
	    gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation,trace_length,start,sbjE,newpos);
	    free(operation);
	    ogsq=SwapGsqCMSA(sq,gsq,cma); ssx->AddToAlign(sq,newpos); free(newpos);
	    if(ogsq) delete []ogsq; 
	} Nildheap(pdH); NilSet(Used); NilSet(Sample);
	SaveBestCMSA(cma);	// save this configuration to restore best below...

	//========= 6. Sample purged sequences some more ===========
	// ...this is like 'S' mode, but with repseq left aligned.
#if 1	
	best=this->RtnMap();
        for(i=1; i <= numSets; i++){
	    best=ssx->GapMap(0);
	    dH = dheap(N+4,3);
            for(sq=1; sq <= N; sq++){
                if(MemberSet(sq,SetSq[i])){ insrtHeap(sq,((keytyp)Random()),dH); }
	    } gmb->SampleTogether(SetSq[i],Temp,map,dH); Nildheap(dH);
	    if(map > best){ best=map; SaveBestCMSA(cma); }
	    NilSet(SetSq[i]); 
	} free(SetSq); 
#else
        for(i=1; i <= numSets; i++) NilSet(SetSq[i]); free(SetSq); 
#endif

	//========= 7. If alignment has improved then use it. ===========
	ssx->RestoreBest(); // reset to optimum CMA configuration.
        map=gmb->RtnMap(); 
	delete gmb; strt_map=bst_lpr; end_map=map;
	if(map > bst_lpr){
		SSX->GetParameters(pio,pdo,eie,ede,pwt,dms_md);
		this->Free(); NilCMSA(bcma); 
		this->init(pio,pdo,eie,ede,cma,dms_md,0); this->SetPriorWt(pwt);
		SSX->RestoreBest(); return TRUE; 
	} else { NilCMSA(cma); return FALSE; } // keep original; copy destroyed.
}

double	gmb_typ::SimilarTogether(Int4 MinSticky, double MaxFrctn,double Temp)
{
	cma_typ	cma=SSX->RtnCMA(); 
        Int4    i,x,numSets,percent_id=40,N=NumSeqsCMSA(cma);
        double  d,map,best=this->RtnMap(); // MaxFrctn = 0.15;
	FILE 	*efp=0; // efp=stderr;

        set_typ *SetSq=RtnTargetSizeSetsCMSA(numSets,percent_id,cma,MaxFrctn);
	// ^ this resets percent_id so that largest set <= MaxFractn of seqs.
	x=CardSet(SetSq[1]);
	d=((double) x/(double)N);	// maximum fraction of seqs sampled together...
	if(efp) fprintf(stderr," max size = %d/%d (%.3f; %d%c)\n",x,N,d,percent_id,'%');
	if(x > MinSticky && d <= MaxFrctn){	// RtnTargetSizeSetsCMSA() ensures this.
          for(i=1; i <= numSets; i++){
	    x=CardSet(SetSq[i]); if(x <= MinSticky) break;
            dh_type dH=dheap(NumSeqsCMSA(cma)+5,4);
            for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
                 if(MemberSet(sq,SetSq[i]))
                        { gsq_typ *gsq=gsqCMSA(sq,cma); insrtHeap(sq,(keytyp)gsq->NumDel(),dH); }
            } this->SampleTogether(SetSq[i],Temp,map,dH); Nildheap(dH); 
	    if(map > best){ best=map; SaveBestCMSA(SSX->RtnCMA()); }
	    if(efp) fprintf(stderr,"%d ",x);
          } if(efp) fprintf(stderr,"\n");
	  SSX->RestoreBest(); // reset to optimum CMA configuration.
	} for(i=1; i <= numSets; i++) NilSet(SetSq[i]); free(SetSq); 
	return SSX->GapMap(0);
}

cma_typ gmb_typ::SampleStickyTogether(FILE *fp,Int4 MinSticky, char m[20], double f[20],
                                double T, char mode)
// perform various operations as specified by char array m and frequence array f. 
{
      char  c;
      cma_typ rcma=0;
      Int4  i,num=0;
      for(i=1; isupper(m[i]); i++){ num++; assert(i < 18 && f[i] > 0.0 && f[i] < 0.50); }
      if(mode=='R'){        // call in random order.
         dh_type    dH=dheap(15,4);
         for(i=1; i<= num; i++) insrtHeap(i,(keytyp) Random(),dH);
         while((i=delminHeap(dH)) != 0){
            rcma=SampleStickyTogether(fp,MinSticky,f[i],0,T,m[i]);
         } Nildheap(dH);
      } else {              // 'O' == call in given order.
         for(i=1; i < num; i++){
            SampleStickyTogether(fp,MinSticky,f[i],0,T,m[i]);
         } return SampleStickyTogether(fp,MinSticky,f[i],0,T,m[i]);
      } return rcma;
}

cma_typ	gmb_typ::SampleStickyTogether(FILE *fp, Int4 MinSticky, double MaxFrctn,
						char *name,double Temperature, char mode)
{
   double	strt_map,map;
   cma_typ	CMA=SSX->RtnCMA(); 
   assert(nBlksCMSA(CMA) == 1);
   Int4		ins_len=0,i,n,dt,timeR=time(NULL),nDelSampled=0;
   SSX->InitNDL(Temperature);	// initialize model parameters?
   if(fp) strt_map=this->RtnMap();
   set_typ set=0,rset=0; 
	
   switch(mode){
     case 'c': { // Sample sequences clusters together ...
	this->SetMaxIter(1); Int4 similarity=1;
	this->Sample(0,'S',similarity,Temperature,Temperature,50,0);
	map=SSX->GapMap(Temperature); } break;
     case 'O': { // Sample sequences one by one...
	this->SetMaxIter(1); Int4 similarity=0;
	this->Sample(0,'S',similarity,Temperature,Temperature,50,0);
	map=SSX->GapMap(Temperature); } break;
     case 'C': { // Sample co-conserved sequences together ...
	map=this->ConservedTogether(MaxFrctn,Temperature); } break;
     case 'S': { // Sample similar sequences together ...
	map=this->SimilarTogether(MinSticky,MaxFrctn,Temperature); } break;
     case 'R': { // Randomly sample sequences together ...
	map=this->RandomTogether(MinSticky,MaxFrctn,Temperature); } break;
     case 'W':	// raw scores against the other sequences...
	map=this->WorstTogether(MinSticky, MaxFrctn, Temperature); break;
     case 'D': 	// deletions between conserved blocks.
      {
       dh_type dH=dheap(TotalLenCMSA(CMA)+5,4);
       for(i=1; i <= LengthCMSA(1,CMA); i++){	// deletions...
	   insrtHeap(i,(keytyp)Random(),dH);	// columns picked in random order...
       }
       for(set=0,n=0; (i=delminHeap(dH)) != 0; ){	// deletions...
	rset=this->GambitTogether(i,MinSticky,MaxFrctn,Temperature,set,0);
#if 0	// attempt to trim the search...
	if(set){ NilSet(set); set=0; }	
	if(rset){ set=rset; n++; }
#else	// ignore sets....
	if(rset){ n++; NilSet(rset); }
#endif
       } if(set) NilSet(set); Nildheap(dH);
	nDelSampled=n; // fprintf(fp,"   (%d/%d columns sampled)\n",n,LengthCMSA(1,CMA));
      } break;
     case 'P': {	// Sample a purged sequence set..
	Int4 percent_id=40; // percent_id=30; // WARNING: this can change CMA!!!
	if(this->SamplePurged(percent_id,MaxFrctn,Temperature,strt_map,map)) CMA=SSX->RtnCMA();
      } break;
     case 'I': ins_len=1;
     case 'V': if(ins_len==0) ins_len=5;
     case 'X': if(ins_len==0) ins_len=10;
     case 'F': if(ins_len==0) ins_len=15;
     case 'N': if(ins_len==0) ins_len=19;
     case 'T': if(ins_len==0) ins_len=30;
      for(i=1; i <= LengthCMSA(1,CMA); i++){
	rset=this->GambitTogether(i,MinSticky,MaxFrctn,Temperature,0,ins_len); map = SSX->GapMap(0);
	if(rset) NilSet(rset);
      } break;
     default:  fprintf(stderr,"mode='%c'\n",mode);
	print_error("SampleStickyTogether() input error"); break;
   } 
   SSX->RestoreBest(); // reset to optimum CMA configuration.
   if(fp){  
	double end_map=this->RtnMap();
	if(mode == 'D'){
	  fprintf(fp,"   StickySampling ('%c' %.3f)(%d/%d cols): LLR = %.2f (%.2f);",
		mode,MaxFrctn,nDelSampled,LengthCMSA(1,CMA),end_map,end_map-strt_map);
	} else {
	  fprintf(fp,"   StickySampling ('%c' %.3f)(%d cols): LLR = %.2f (%.2f);",
		mode,MaxFrctn,LengthCMSA(1,CMA),end_map,end_map-strt_map);
	} dt=time(NULL)-timeR; double WtSq=this->RtnWtNumSeqs();
	if(dt < 60) fprintf(fp," %.1f wt.sq. (%d secs)(%.1f K)\n", WtSq,dt,Temperature);
	else fprintf(fp," %.2f WtSq (%0.2f mins)(%.1f K)\n", WtSq,(float)(dt)/60.0,Temperature);
   } CMA=SSX->RtnCMA(); return CMA;
}

double	gmb_typ::ConservedTogether(double MaxFrctn, double Temp)
// Find sequence sets that are sample together candidates...
{
	SetPseudoToMapCMSA(SSX->RtnCMA()); 
	SSX->RestoreBest(); // reset to optimum CMA configuration; retains this as best...
	cma_typ CMA=SSX->RtnCMA();
	assert(nBlksCMSA(CMA)==1);
	double	lpr,best=SSX->GapMap(0);
	Int4	*obs,sq,i,j,J,r,rr,x,NN=NumSeqsCMSA(CMA);
	double  target,bild,d,dd,frq,*tfrq=tFreqCMSA(CMA);
	double TotalWtSq=SSX->TotalWtSeq();
	rst_typ *rst = new rst_typ('G',AB);
	FILE *efp=0; // efp=stderr;
	dh_type  dH,dHI;

	//============ 1. Sort columns by # of deletions... ===============
	dH=dheap(TotalLenCMSA(CMA)+5,4);
	h_type HG=Histogram("BILD scores",-500,500,0.10);
	for(i=1; i <= LengthCMSA(1,CMA); i++){
	    dd=NumDeletionsCMSA(1,i,CMA); 
	    insrtHeap(i,(keytyp)dd,dH); // want to use non-deleted columns...
	    d=SSX->BildScore(1,i)/TotalWtSq; IncdHist(d,HG);
	} if(efp) PutHist(efp,60,HG); 
	double stdev=sqrt(VarianceHist(HG)),mean=MeanHist(HG); NilHist(HG);
	double low=mean - stdev,high=mean + stdev;

	//============ 2. Save columns with fewest deletions... ===============
	dHI=dheap(TotalLenCMSA(CMA)+5,4);
	Int4 Stop=TotalLenCMSA(CMA)/2; assert(Stop > 0);
	for(J=1; !emptyHeap(dH); J++){
		if(J > Stop) break;
		d = minkeyHeap(dH); assert((i=delminHeap(dH)) != 0);
		// float   *RelEntropyCMA(1, CMA);
		dd=SSX->BildScore(1,i)/TotalWtSq; 
		// if(dd >= mean && dd <= high) insrtHeap(i,(keytyp)-dd,dHI); // want to remove best first...
		if(dd >= mean) insrtHeap(i,(keytyp)-dd,dHI); // remove best first...
	} Nildheap(dH);

	// Sort by most conserved columns...
	dH=dheap(TotalLenCMSA(CMA)+5,4);
	Stop=10;	
	for(J=1; !emptyHeap(dHI); J++){
	   if(J > Stop) break;
	   bild = -minkeyHeap(dHI); assert((i=delminHeap(dHI)) != 0);
	   // insrtHeap(i,(keytyp)-bild,dH);	// best columns first...
	   insrtHeap(i,(keytyp)Random(),dH);	// columns picked in random order...
	} Nildheap(dHI);  dHI=dH; dH=0;
	
	// sample sequences with a common conserved residue at specific positions...
	sst_typ fsst=0; fsst=SsetLet(AlphaCode('L',AB));
	fsst=UnionSset(fsst,SsetLet(AlphaCode('I',AB)));
	fsst=UnionSset(fsst,SsetLet(AlphaCode('V',AB)));
	if(efp){ PutSST(stderr,fsst,AB); fprintf(efp," are forbidden residues.\n"); }
	for(J=1; !emptyHeap(dHI); J++){
	   bild = -minkeyHeap(dHI); assert((i=delminHeap(dHI)) != 0);
	   obs=ColResCntsCMSA(1,i,CMA);
	   if(efp) fprintf(efp,"%d: bild=%.3f\n",i,bild);
	   dH=dheap(nAlpha(AB)+3,4);
	   for(r=1; r <= nAlpha(AB); r++){
		if(MemSset(r,fsst)) continue;
		frq=(double) obs[r]/(double) NN; 
		if(frq==0) continue; else d=frq/tfrq[r];
		if(d > 1.2){ insrtHeap(r,(keytyp)-d,dH); }
	   }
	   if(emptyHeap(dH)){ free(obs); Nildheap(dH);  continue; } 
#if 0
	   rst->IsLegalSet(sst);
	   Int4 r1 = ResidueCMSA(1,sq,Position0[j],IN_CMA[i]);
	   if(MemSset(r1,Residues[j])){ Match[i][j]++; hits++; }
#endif
	   set_typ	*set; NEW(set,nAlpha(AB)+3,set_typ);
	   sst_typ xsst=0,sst=0;	// MemSset(i,s) sst=SsetLet(r); sst=UnionSset(s,r);
	   for( ; !emptyHeap(dH); ){
		d=-minkeyHeap(dH); r=delminHeap(dH); 	
		dd=(double) obs[r]/(double) NN;
		if(dd <= MaxFrctn){
		   xsst=SsetLet(r); sst=UnionSset(sst,xsst);
		   set[r]=MakeSet(NN+5); ClearSet(set[r]);
	    	   if(efp) fprintf(efp,"  %c=%.3f (%d) bg: %.3f; LR = %.3f.\n",
			AlphaChar(r,AB),dd,obs[r],tfrq[r],d);
		}
	   } if(efp) fprintf(efp,"\n"); Nildheap(dH); free(obs);

	   // Add relevant sequences to each set
	   for(sq=1; sq <= NN; sq++){
		r=ResidueCMSA(1,sq,i,CMA);
		if(MemSset(r,sst)){ AddSet(sq,set[r]); }
	   } 

	   // Sample sequences at this site.
	   for(r=1; r <= nAlpha(AB); r++){
		if(set[r]){
	   	   dH=dheap(NN+10,4);
	   	   for(sq=1; sq <= NN; sq++){
		      if(!MemberSet(sq,set[r])) continue;
		      // d=GetGappedProbCMSA(1,sq,CMA); if(!isfinite(d)) d=0;
#if 0	// "Seq sq: prob = 0!" happening...
		      d=GetProbCMSA(1,sq,CMA); if(!isfinite(d)) d=0;
		      insrtHeap(sq,(keytyp)d,dH);	// sample the worst first; the best are sticky.
#else
	              insrtHeap(sq,((keytyp)Random()),dH); 
#endif
		   }
	           this->SampleTogether(set[r],Temp,lpr,dH);
		   // if(efp) fprintf(efp,"lpr = %.3f\n",lpr);
           	   if(lpr > best){ 
			  if(efp) fprintf(efp,"LLR changes from %.1f to %.1f (%.1f K)\n",
					best,lpr,Temp);
       			  best=lpr; SaveBestCMSA(CMA);
           	   } else { SSX->RestoreBest(); }
		   Nildheap(dH);
		   NilSet(set[r]);
		}
	   } free(set); 
	} Nildheap(dHI);
	delete rst;
	SSX->RestoreBest(); // reset to best CMA configuration; retains stored best...
	return SSX->GapMap(0);
}

BooLean	gmb_typ::SampleByLayers(Int4 MinSize,double Temp, double &strt_map,double &end_map, Int4 stage)
// Sample most dissimilar first.
{
        cma_typ bcma=SSX->RtnCMA(),cma=CopyCMSA(bcma);
        Int4    u,x,numSets,sq,i,j,k,n,o,N=NumSeqsCMSA(cma),*oldpos,percent_id;
        double  d,map,best,bst_lpr;
	Int4	end_pid=30;
	dh_type	dH=0,pdH=0;
        SaveBestCMSA(cma);
	assert(MinSize >= 40);
	FILE	*efp=0; // efp=stderr;

	//=========== 1. Find related sequence sets... ==============
        set_typ *SetSq=0; NEW(SetSq, 20, set_typ);	// need 14 sets.
        Int4 *pid; NEW(pid, 20, Int4);	// need 14 sets.
	SetSq[0]=MakeSet(N+5); FillSet(SetSq[0]); pid[0]=100;
	for(i=N+1; i < SetN(SetSq[0]); i++) DeleteSet(i,SetSq[0]); DeleteSet(0,SetSq[0]);
	for(i=0,percent_id=95; percent_id >= end_pid; percent_id -=5){
		i++; SetSq[i]=RtnFastRepSetCMSA(0, percent_id,SetSq[i-1],cma);
		pid[i]=percent_id;

	} numSets=i; 
	for(i=numSets; i > 0; i--){
	   x=CardSet(SetSq[i]);
	   if(x < MinSize){
		UnionSet(SetSq[i-1],SetSq[i]);	// combine set i and i-1.
		NilSet(SetSq[i]);  SetSq[i]=0; numSets--;
	  } else break;
	}
	if(numSets <= 1){
           for(i=0; i <= numSets;  i++){ NilSet(SetSq[i]); } free(SetSq);  free(pid);
	   NilCMSA(cma); return FALSE; 
	}
	//=========== 2. Sort these by % identity cutoff... ==============
	set_typ SetU=0,*SqSet=0; NEW(SqSet, numSets+5, set_typ);
	SqSet[0]=MakeSet(SetN(SetSq[1])); FillSet(SqSet[0]);
	for(i=N+1; i < SetN(SqSet[0]); i++) DeleteSet(i,SqSet[0]); DeleteSet(0,SqSet[0]);
        for(j=0,i=numSets; i >= 0;  i--){
	    j++; x=CardSet(SetSq[i]);  assert(x > 0);
	    if(i==numSets){ SetU=CopySet(SetSq[i]); SqSet[j]=CopySet(SetSq[i]);}
	    else {
		SqSet[j]=MakeSet(SetN(SetU));
		IntersectNotSet(SetSq[i],SetU,SqSet[j]); //  SqSet[j] = SetSq[i] & not SetU.
		UnionSet(SetU,SetSq[i]);
	    }
	    n=CardSet(SqSet[j]); o=CardInterSetINotJ(SetU,SqSet[j]); u=CardSet(SetU);
	    if(efp) fprintf(stderr,"%d: %d%c ident; %d seqs (%d new; %d old; union = %d).\n",
				i,pid[i],'%',x,n,o,u);
	}
        for(i=0; i <= numSets;  i++){ NilSet(SetSq[i]); } free(SetSq); 
	numSets=j;

	//========= 3. create a new gmb_typ.
	Int4	pio,pdo,eie,ede;
	double	pwt;
	char	dms_md;
	SSX->GetParameters(pio,pdo,eie,ede,pwt,dms_md);
	// this->Free(); NilCMSA(cma); 
	gmb_typ *gmb = new gmb_typ(pio,pdo,eie,ede,cma,dms_md,pwt,0); // gmb->SetPriorWt(pwt);
	gmb->RestoreBest(); bst_lpr=gmb->RtnMap(); // MaxFrctn = 0.15;
	ssx_typ	*ssx=gmb->RtnSSX();

#if 1	// Get Keys for heap...Add best matches to pre-aligned sequences first...
     	double  *Key; NEW(Key,N+5,double); ClearSet(SetU);
     	cma_typ *scma; NEW(scma, numSets+5, cma_typ);
     	for(i=1; i <= numSets;  i++){
	   UnionSet(SetU,SqSet[i]);
	   scma[i]=MkSubCMSA(SetU,TRUE,cma); 
	   // scma[i]=MkSubCMSA(SqSet[i],TRUE,cma); 
	}
        gss_typ *gss=gssCMSA(cma);
     	for(j=0,i=1; i <= numSets; j++,i++){
           for(sq=1; sq <= N; sq++){
	     if(!MemberSet(sq,SqSet[i])) continue;
	     if(j==0) Key[sq]=(double)Random();
	     else {
#if 0		// Don't use this: best in first; these are more likely to be stuck???
		Key[sq] = -(double)PseudoAlnScoreTwoCMSA(1,scma[j],sq,cma);
		Key[sq] -= (double) (gss->NumDel(sq)*lowAlphaR(AB));
#elif 0		// not as good as worst first?
	        Key[sq]=(double)Random();
#else		//  worst in first; more consistent...
		Key[sq] = (double)PseudoAlnScoreTwoCMSA(1,scma[j],sq,cma);
		Key[sq] += (double) (gss->NumDel(sq)*lowAlphaR(AB));
#endif
	     }
	   }
     	}
     	for(i=1; i <= numSets;  i++){ TotalNilCMSA(scma[i]); } free(scma);
#endif
	NilSet(SetU);

     for(i=1; i <= numSets;  i++){
	//========= 3. Store rep. seqs on a heap; remove other seqs from alignment.
	set_typ Used=SqSet[i];
	dH = dheap(N+4,3);
        for(sq=1; sq <= N; sq++){
	     // if(MemberSet(sq,Used)){ insrtHeap(sq,((keytyp)Random()),dH); }
	     if(MemberSet(sq,Used)){ insrtHeap(sq,(keytyp)Key[sq],dH); }
	     else if(i == 1){ oldpos=ssx->RemoveFromAlign(sq); free(oldpos); }
	} 
	if(i > 1){	// 5. Sample the new sequences back in one-at-a-time.
	  Int4	start=0,score,trace_length,*newpos;
	  gss_typ *gss=gssCMSA(cma);
	  for(j=1; (sq=delminHeap(dH)) != 0; j++){
		assert(MemberSet(sq,Used));
		assert(!SeqIsInCMSA(sq,cma));	// i.e., sq has no aligned sites.
		e_type  sbjE=gss->TrueSeq(sq);
		ssx->InitNDL(Temp);     
		char    *operation=ssx->GapAlnTrace(sbjE,Temp,start,score);
		trace_length=strlen(operation);
#if 0
	if(PdbSeq(sbjE)){
	   put_seqaln_smatrixSW(stderr,operation,LenSeq(sbjE)-start+1,
                SeqPtr(sbjE)+start-1,OffSetSeq(sbjE)+start-1,trace_length,nBlksCMSA(cma),ssx->RtnSMX());
	   PutSeqID(stderr,sbjE);
	   fprintf(stderr,"%d. operation=%s\n",sq,operation);
	}
#endif
		gsq_typ *ogsq,*gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(cma)+3,Int4);
		gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation,trace_length,start,sbjE,newpos);
		free(operation);
		ogsq=SwapGsqCMSA(sq,gsq,cma); ssx->AddToAlign(sq,newpos); free(newpos);
	        if(ogsq) delete []ogsq; 
	  } 
          for(sq=1; sq <= N; sq++){
	     if(MemberSet(sq,Used)){ insrtHeap(sq,((keytyp)Random()),dH); } 
	     // if(MemberSet(sq,Used)){ insrtHeap(sq,(keytyp)Key[sq],dH); }
	  }
	} SaveBestCMSA(cma); // wipes away previous best...

	//=============== 4. Sample added sequences one-at-a-time.
	best=gmb->RtnMap();
	double T=Temp;
#if 0
	for(Int4 I=1; I <= 2; T=0.0,I++){
	  for(j=1; (sq=delminHeap(dH)) != 0; j++){
	    assert(MemberSet(sq,Used)); ssx->InitNDL(Temp);     
	    // if(gmb->SampleSingle(sq,Temp,map))
	    if(gmb->SampleSingle(sq,T,map)){
            	if(map > best){ best=map; SaveBestCMSA(cma); }
	    }
	  }
	  if(i > 1) break;
	  else if(I == 1){	// sample twice over initial set.
            for(sq=1; sq <= N; sq++){
	      if(MemberSet(sq,Used)){ insrtHeap(sq,((keytyp)Random()),dH); } 
	    }
	  }
	}
#else
	for(j=1; (sq=delminHeap(dH)) != 0; j++){
	    assert(MemberSet(sq,Used)); ssx->InitNDL(Temp);     
	    // if(gmb->SampleSingle(sq,Temp,map))
	    if(gmb->SampleSingle(sq,T,map)){
            	if(map > best){ best=map; SaveBestCMSA(cma); }
	    }
	}
#endif
	Nildheap(dH); gmb->RestoreBest(); // reset to optimum CMA configuration.
	if(0) fprintf(stderr,"%d: %d new seqs added.\n",i,CardSet(SqSet[i]));
	NilSet(SqSet[i]); 
     } NilSet(SqSet[0]); free(SqSet); free(pid); free(Key);
    
	ssx->RestoreBest(); // reset to optimum CMA configuration.
        map=gmb->RtnMap(); // MaxFrctn = 0.15;
	delete gmb; strt_map=bst_lpr; end_map=map;
	if(map > bst_lpr){
		SSX->GetParameters(pio,pdo,eie,ede,pwt,dms_md);
		this->Free(); NilCMSA(bcma); 
		this->init(pio,pdo,eie,ede,cma,dms_md,0); this->SetPriorWt(pwt);
		SSX->RestoreBest(); return TRUE; 
	} else { NilCMSA(cma); return FALSE; } // Destroy copy & revert to previous cma.
}

