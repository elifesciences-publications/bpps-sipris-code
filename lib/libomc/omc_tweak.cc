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

#include "omc_typ.h"

void	omc_typ::PartitionNode(set_typ *set, Int4 pI, Int4 I, e_type bstE, e_type csqE)
// Partition the set[pI] into two sets, the other being set[I] based on sequence bstE.
{
	Int4	sq,n,x,score,sc,top_sc,top_sq,max,wst_sc,wst_sq;

	for(sc=0,score=0,x=1; x <= LenSeq(bstE); x++){
	   if(ResSeq(x,bstE)==0) EqSeq(x,AlphaCode('A',AB),bstE);  // need for non-map_typ routine.
	   // assert(ResSeq(x,bstE) != 0);
	   score+=valAlphaR(ResSeq(x,bstE),ResSeq(x,bstE),AB); // self-score of best vs best.
	   sc+=valAlphaR(ResSeq(x,bstE),ResSeq(x,csqE),AB); // score of best vs consensus.
	} 
#if 0	// debug
	double	PerCent=100.0*(double)sc/(double)score;
	fprintf(stderr,"Score = %d vs %d (%.1f%c)\n",sc,score,PerCent,'%');
#endif
     	for(n=0,top_sq=0,top_sc=INT4_MIN,sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
           if(MemberSet(sq,set[pI])){        
		score=PseudoAlnScoreSqToCMSA(bstE,sq,TrueMainCMA);
		if(score > top_sc){ top_sc=score; top_sq=sq; }
		sc=PseudoAlnScoreSqToCMSA(csqE,sq,TrueMainCMA);
		if(score > sc){ DeleteSet(sq,set[pI]); AddSet(sq,set[I]); } 
	   }
	} 
	if(CardSet(set[I]) < 1){ DeleteSet(top_sq,set[pI]); AddSet(top_sq,set[I]); } 
// fprintf(stderr,"Card set[%d] = %d; Card parent set [%d] = %d\n",I,CardSet(set[I]),pI,CardSet(set[pI]));
}

e_type	omc_typ::FindSeedSeq(Info *info, Int4 pI, Int4 I,sst_typ *osst,FILE *efp)
//========= Get seed sequences based on pattern matches =====
{
	e_type	bstE;
	Int4	n,m,top_sc,top_sq,sq,k,x,r,len=mcs->RtnLengthMainCMSA();
	double	sd,mean,median;

	ClearSet(info->set[I]);

	h_type HG=Histogram("pattern matches",-10,len,1.0);
	for(top_sc=top_sq=0,sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
          if(!MemberSet(sq,info->set[pI])) continue;
          for(n=0,m=0,k=1; k <= len; k++){
                if(osst[k]==0) continue;
                if(MemSset(ResidueCMSA(1,sq,k,TrueMainCMA),osst[k])) n++; else m++;
          } // Put sequences with more matches than mismatches into child set.
          IncdHist(n,HG);
          if(n > top_sc){ top_sc=n; top_sq=sq; }
	}  if(efp) PutHist(efp,60,HG); 
	median=MedianHist(HG); mean=MeanHist(HG); sd=sqrt(VarianceHist(HG)); NilHist(HG);
	if(efp) fprintf(efp,"mean = %.2f; + 1 sd = %.2f; + 2 std = %.2f; median = %.2f\n",
			mean,mean+sd,mean+2*sd,median);
	if(efp) fprintf(efp,"top_sq = %d; top_sc = %d\n",top_sq,top_sc);

	if(sd > 8.0) x=top_sc-2; else if(sd > 3.0) x=top_sc-1; else x=top_sc;
	for(sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
                   if(!MemberSet(sq,info->set[pI])) continue;
                   for(n=0,k=1; k <= len; k++){
                     if(osst[k]==0) continue;
                     if(MemSset(ResidueCMSA(1,sq,k,TrueMainCMA),osst[k])) n++;
                   } if(n >= x){ AddSet(sq,info->set[I]); }
	} bstE=GetSeqAsCsqCMSA(info->set[I],TrueMainCMA); 
	if(efp) fprintf(efp,"  %d seqs selected\n",CardSet(info->set[I]));
	ClearSet(info->set[I]);
        char Ala=AlphaCode('A',AB);
        for(k=1; k <= len; k++){
            if(osst[k]==0){ // then if 'X' change to 'A'
               if(ResSeq(k,bstE) == 0){ EqSeq(k,Ala,bstE); EqXSeq(k,Ala,bstE); }
            } else if(!MemSset(ResSeq(k,bstE),osst[k])){
               for(r=1; r <= nAlpha(AB); r++){
                  if(MemSset(r,osst[k])){ EqSeq(k,r,bstE); EqXSeq(k,r,bstE); break; }
               }
            }
        } return bstE;
}

set_typ	omc_typ::MkSubTreeSet(Info *info,Int4 pI)
{
	assert(pI > 0 && pI < Hpt->NumSets());
	if(pI != 1 && Hpt->TypeOfSet(pI) != '?') return 0;
	set_typ setR=CopySet(info->set[pI]); ClearSet(setR);
	set_typ setH=Hpt->MkSubTreeSet(pI);
	for(Int4 i=1; i < Hpt->NumSets(); i++){
	   if(MemberSet(i,setH)){ UnionSet(setR,info->set[i]); }
	} NilSet(setH); return setR;
}

void	omc_typ::PatternBasedPartition(set_typ *set, sst_typ *sst,Int4 pI, Int4 I, double alpha)
// Partition the sequences based on the likelihood that they would have matched the sequences by chance...
{
	Int4	r,sq,j,i,miss,hits,del,**observed;
	Int4	N=NumSeqsCMSA(TrueMainCMA),len=LengthCMSA(1,TrueMainCMA);
	double	*p,*q,P,Q,d,D,**freq=ColResFreqsCMSA(1,set[pI],&observed,TrueMainCMA);
	double	a=alpha,m1a=1.0-alpha; // alpha prior parameter...freq of matches to pattern.

	h_type HG=Histogram("pattern matches",-100,100,5.0);
	NEW(p,len+3,double); NEW(q,len+3,double);
        for(j=1; j <= len; j++){
	    if(sst[j]==0){ free(freq[j]); free(observed[j]); continue; }
	    // fprintf(stderr,"%d: ",j); PutSST(stderr,sst[j],AB); // fprintf(stderr,"\n");
	    PutSST(stderr,sst[j],AB); fprintf(stderr,"%d ",j);
	    for(D=freq[j][0],d=0.0,r=1; r<= nAlpha(AB); r++){ 
		D+= freq[j][r];
		// fprintf(stderr," %c %.3f\n",AlphaChar(r,AB),freq[j][r]);
		if(MemSset(r,sst[j])) d+=freq[j][r]; 
	    } free(freq[j]); free(observed[j]);
	    // fprintf(stderr," total = %.3f (%.3f)\n",D,d);
	    // use alpha parameter here...
	    // D = ((m1a*d)+a);	// (1 - alpha)*d + alpha;
	    p[j]=log(a/d); q[j]=log(m1a/(1.0-d)); 
	} free(freq); free(observed);
	// D=-log((double)CardSet(set[pI]));
	// D = (double)CardSet(set[pI]);
        for(sq=1; sq <= N; sq++){
              if(MemberSet(sq,set[pI])){        
                for(P=Q=0.0,miss=del=0,j=1; j <= len; j++){
		   if(sst[j]==0) continue;
                   r = ResidueCMSA(1,sq,j,TrueMainCMA);
                   if(MemSset(r,sst[j])){ P += p[j]; } else {
			P+=q[j];
			if(r != 0){ miss++; } else { del++; }
		   }
                } IncdHist(P,HG);
                // if((miss + del) <= MaxMisMatches){ DeleteSet(sq,set[pI]); AddSet(sq,set[I]); }
                if(P > 0.0){ DeleteSet(sq,set[pI]); AddSet(sq,set[I]); }
              }
	} free(p); free(q);
	fprintf(stderr,"\n"); PutHist(stderr,60,HG); NilHist(HG); 
}

Int4	omc_typ::AdjustSetIDs()
// change set ids so that the max_id == NumBPPS; return old ID and new ID.
{
     Int4    last_id=0,max_id,id,i,j,oID,nID,N;
     set_typ idSet=MakeSet(MaxNumNodesPlus+2); ClearSet(idSet);
// mcs->PutHyperPartition(stderr);
     dh_type dH=dheap(MaxNumNodesPlus+2,4);
// Hpt->Put(stderr); if(stable) PutSet(stderr,stable);
     set_typ nStable=0;
     if(stable){ nStable=CopySet(stable); ClearSet(nStable); }
     // 1. Use a Help to sort rows in Hpt by their set ids.
     for(max_id=0,i=1; i <= Hpt->NumBPPS(); i++){
	id=Hpt->ItoSetID(i); AddSet(id,idSet);	// Store this in the idSet
	insrtHeap(i,(keytyp)id,dH);	
	if(id > max_id){ max_id=id; } // Find set with maximum Set ID.
     }
     // stk_typ stack(max_id+2); 	// last item in is the first item out.
     stk_typ *stack = new stk_typ(max_id+2); 	// last item in is the first item out.
     for(j=1; (i=delminHeap(dH)) != NULL; j++){ 
	id=Hpt->ItoSetID(i); assert(stack->Push(id)); 
     } Nildheap(dH); stack->Put(stderr); 
// PutSet(stderr,idSet); fprintf(stderr,"max_id = %d\n",max_id);
     for(N=0,id=1; id <= max_id; id++){	// find missing ids.
	last_id=id;
	if(MemberSet(id,idSet)){
		if(stable && MemberSet(id,stable)){ AddSet(id,nStable); } continue; 
	} else N++;
	nID=id; oID=stack->Pop(); 	// returns 0 if stack is empty.
	if(stable && MemberSet(oID,stable)){ AddSet(nID,nStable); }
	sprintf(str,"Set%d",nID);  assert(Hpt == mcs->GetHpt());
	i=Hpt->SetIDtoI(oID); Hpt->ReNameSet(i,str); mcs->RenameDisplayCMA(i,str);
	max_id=stack->Top();	// returns the item on the top of the Stack w/o removing it.
fprintf(stderr,"%d: Changed Set%d to Set%d (%s); max = %d.\n",i,oID,nID,str,max_id);
	stack->Put(stderr);
	AddSet(nID,idSet); DeleteSet(oID,idSet);
     } NilSet(idSet); // max_id=stack->Pop(); 
// mcs->PutSubLPRs(stderr);
// if(N > 0) mcs->PutHyperPartition(stderr);
     if(stable){ NilSet(stable); stable=nStable; mcs->SetStringency(MinimumLLR,MinimumSetSize,stable); }
     delete stack;
     assert(last_id == Hpt->NumBPPS()); 
// Hpt->Put(stderr); if(stable) PutSet(stderr,stable);
// if(N>0) assert(5 == 3);
}

