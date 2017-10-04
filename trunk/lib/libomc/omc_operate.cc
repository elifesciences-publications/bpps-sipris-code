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
#include "blosum62.h"

//*********************************** AddLeaf() ********************************************
Int4	omc_typ::AddLeaves(double Temp,char Action)
{
	Int4	i,j,x,id,n,card,NumAdded=0,*tries,iter,*Ntry;
        double d,D,*map = mcs->RtnSubLPR();      // map[0] == total LPR.
	dh_type dH=dheap(MaxNumNodes+2,4);
	BooLean	RootOnly=FALSE;
	// if(Action == 'e'){ RootOnly=TRUE; Action ='a'; } // ResurrectNodes calls this...
	if(Action == 'e'){ RootOnly=TRUE; Action ='A'; }
#if 1
	else if(Action == 'z'){ RootOnly=TRUE; Action ='a'; }
	else if(Action == 'Z'){ RootOnly=TRUE; Action ='A'; }
#endif
	assert(Action=='A' || Action == 'a');
	NEW(tries,MaxNumNodesPlus+2,Int4);
	NEW(Ntry,MaxNumNodesPlus+2,Int4);
        for(i=1; i < Hpt->NumSets(); i++){  // NumBPPS == NumSets - one Reject set.
	    if(Hpt->NodeDepth(i) >= this->MaxDepthHierarchy) continue;
	    id=Hpt->ItoSetID(i); 
	    if(stable && MemberSet(id,stable)) continue;
	    card=mcs->SetCard(i);
	    // fprintf(stderr,"%d: id=%d; card = %d\n",i,id,card);
            if(card >= MinimumSplitSize){
		tries[id]=card/MinimumSplitSize; assert(tries[id] > 0);
		tries[id]= MINIMUM(Int4,tries[id],3);	// 3 strikes and you're out.
		insrtHeap(id,-(keytyp)card,dH); 
	    } if(RootOnly) break;
        }
	h_type HG=Histogram("the number of retries by id",0,MaxNumNodes,1);
	for(NumAdded=0,iter=1; (id=delminHeap(dH)) != NULL; iter++){
	    n=Hpt->SetIDtoI(id); 
	    if(n == 0) continue;  // this node has been deleted...!!!
            // fprintf(stderr,"======= Sampling a leaf attached to node %d (Set%d). =======\n",n,id);
            Int4 x,I,ID;
            for(ID=1; ID <= Hpt->NumSets(); ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID:
	    assert(ID <= Hpt->NumSets());
	    mcs_typ *rtn_mcs =0;	// rtn_mcs == attached leaf node ('Set<ID>') to node pID.
	    Ntry[id]++;
	    rtn_mcs=this->Operate(Action,id,ID,Temp); 
	    IncdHist(0.001+(double)id,HG); 
	    if(rtn_mcs == 0){ tries[id]=0; fprintf(stderr," ------Set%d' aborted. ------\n",ID); }
	    else {
	      hpt_typ *hpt=rtn_mcs->GetHpt(); i=hpt->SetIDtoI(ID); CalcLLR(rtn_mcs);
              double *Map=rtn_mcs->RtnSubLPR( ),D=mcs->RtnBestLPR(); 
	      if(testHG){ IncdHist(AddLeafLLR-Map[i],testHG); PutHist(stderr,60,testHG); }
	      fprintf(stderr,
		"NewLLR=%.2f; subLLR=%.2f (%.2f min); OldLLR=%.2f; node %d ('Set%d %c'); %d seqs (%d min) %.1f npws.\n",
                  Map[0],Map[i],MinimumLLR,D,i,ID,hpt->TypeOfSet(i),rtn_mcs->SetCard(i),MinimumSetSize,
		rtn_mcs->RtnNatsPerWtSeq(i));
	      if((Map[0] > D && Map[i] > MinimumLLR && (hpt->TypeOfSet(i) == '?'
				|| rtn_mcs->SetCard(i) >= MinimumSetSize))
				&& rtn_mcs->OkayNatsPerWtSq(i)){
		if(triesHG){ IncdHist(Ntry[id],triesHG); PutHist(stderr,60,triesHG); } Ntry[id]=0;
		rtn_mcs->PutHyperPartition(stderr);  
		fprintf(stderr,
		   "!!!!!!!!!!!!!!!!!! Sampling node %d ('Set%d') sucessful. !!!!!!!!!!!!!!!!!!!!\n",i,ID);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); this->RestoreBest();

	        if(goodHG && AddLeafLLR > 0.0) IncdHist(AddLeafLLR,goodHG); 
		NumAdded++;  if(NumAdded % 8 == 0) mcs->Sample(2,2,2,2,mcs->TargetMod);
	        i=Hpt->SetIDtoI(id); card=mcs->SetCard(i);
	        if(card >= MinimumSplitSize){   // Then insert the  parent node.
			tries[id]=card/MinimumSplitSize; assert(tries[id] > 0);
			tries[id]= MINIMUM(Int4,tries[id],3);	// no more than 3 failures in a row.
			insrtHeap(id,-(keytyp)card,dH); 
		}
	        if(!RootOnly && !fullHeap(dH)){	// insert child into heap; should never be full!
	      	  i=Hpt->SetIDtoI(ID); card=mcs->SetCard(i);
	      	  if(card >= MinimumSplitSize){
			tries[ID]=card/MinimumSplitSize; assert(tries[ID] > 0);
			tries[ID]= MINIMUM(Int4,tries[ID],3);	// no more than 3 failures in a row.
			insrtHeap(ID,-(keytyp)card,dH); 
		  }
	        } else if(!RootOnly) fprintf(stderr,"WARNING: leaf heap is full\n");
	      } else {
		DeleteMCS(rtn_mcs); tries[id]--;  
	        if(Action != 'a' && tries[id] > 0){
	          i=Hpt->SetIDtoI(id); card=mcs->SetCard(i); 
		  if(card >= MinimumSplitSize) insrtHeap(id,-(keytyp)card,dH); 
	        } fprintf(stderr," ------ #!#!# --- Node 'Set%d' failed. --- #!#!# ------\n",ID); 
	        if(badHG && AddLeafLLR > 0.0) IncdHist(AddLeafLLR,badHG); 
	      }
	    }
	} free(tries); Nildheap(dH); PutHist(stderr,60,HG); NilHist(HG); free(Ntry);
	return NumAdded;
}

// mcs_typ *omc_typ::AddLeaf(Int4 pID, Int4 ID,double T){ return this->Operate('a',pID,ID,T); }
mcs_typ *omc_typ::AddLeaf(Int4 pID, Int4 ID,double T){ return this->Operate('A',pID,ID,T); }
//*********************************** end AddLeaf() ********************************************

//*********************************** DeleteNodes() ********************************************
Int4	omc_typ::DeleteNodes(double Temp, double minLLR, Int4 minCard)
// Delete nodes with LLRs < minLLR
{
	Int4	id,i,j,worst_i,Ndel=0,Worst=-1,card;
	double	d,d_min=DBL_MAX,worst=DBL_MAX,worst_npws=DBL_MAX,npws;
	mcs->RestoreBest();
	for(i=2; i <= Hpt->NumBPPS(); i++){	// 1. Eliminate empty leaf nodes...
	    if(Hpt->TypeOfSet(i) != '!') continue;
	    if(mcs->SetCard(i) == 0){
		id=Hpt->ItoSetID(i); mcs->SaveBest=TRUE; 
	        if(stable && MemberSet(id,stable)) DeleteSet(id,stable);
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); Ndel++;
fprintf(stderr,"Removing empty leaf node %d ('Set%d').\n",i,id);
	    }
	} 
#if 1	// Remove leaf nodes too far down the hierarchy
	for(i=2; i <= Hpt->NumBPPS(); i++){	
	    if(Hpt->TypeOfSet(i) != '!') continue;  // Skip non-leaf nodes...
	    if(Hpt->NodeDepth(i) > this->MaxDepthHierarchy){
		id=Hpt->ItoSetID(i); mcs->SaveBest=TRUE; 
	        if(stable && MemberSet(id,stable)) DeleteSet(id,stable);
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); Ndel++;
fprintf(stderr,"Removing too-deep leaf node %d ('Set%d').\n",i,id);
	    }
	} 
#endif
#if 1
	// Eliminate nodes with too few pattern positions.
	for(i=2; i <= Hpt->NumBPPS(); i++){	// 1. Eliminate empty leaf nodes...
	    // if(Hpt->TypeOfSet(i) != '!') continue;
	    if(mcs->NumColumns(i) < mcs->MinNumColumns(i)){
		id=Hpt->ItoSetID(i); mcs->SaveBest=TRUE; 
	        if(stable && MemberSet(id,stable)) DeleteSet(id,stable);
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); Ndel++;
fprintf(stderr,"Removing node %d ('Set%d') due to too few columns.\n",i,id);
	    }
	} 
#endif
#if 1	// new addition...
	do {			// 2. Eliminate nodes with poor NatsPerWtSeq.
	   //!!!!!!!!!!!!! Eventually only use MinNatsPerWtSq as determinant: delete minLLR !!!!!!!!!!!!!
	   double  *Map=mcs->RtnSubLPR( );
	   for(worst_npws=DBL_MAX,worst_i=0, i=2; i <= Hpt->NumBPPS(); i++){
	      if(stable && MemberSet(Hpt->ItoSetID(i),stable)) continue;
	      npws=mcs->RtnNatsPerWtSeq(i);
	      if(!mcs->OkayNatsPerWtSq(i) && npws < worst_npws){ worst_i=i; worst_npws=npws; } 
	   }
	   if(worst_npws < mcs->MinNatsPerWtSq()){
		id=Hpt->ItoSetID(worst_i); mcs->SaveBest=TRUE; 
		mcs->SampleColumns(); mcs->RestoreBest();
// mcs->PutHyperPartition(stderr); // exit(1);
fprintf(stderr,"Removing node %d ('Set%d') due to low nats per wt_seq (%.2f).\n",worst_i,id,worst_npws);
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
		Ndel++; // mcs->PutHyperPartition(stderr);
	   }
	} while(worst_npws < mcs->MinNatsPerWtSq());
#endif
	do {			// 2. Eliminate nodes with poor LLRs.
	   //!!!!!!!!!!!!! Eventually only use MinNatsPerWtSq as determinant: delete minLLR !!!!!!!!!!!!!
	   double  *Map=mcs->RtnSubLPR( );
	   for(worst=DBL_MAX, worst_i=0, i=2; i <= Hpt->NumBPPS(); i++){
	      if(stable && MemberSet(Hpt->ItoSetID(i),stable)) continue;
	      if(Map[i] < worst){ worst = Map[i]; worst_i=i; } 
	   }
	   if(worst < minLLR){
		id=Hpt->ItoSetID(worst_i); mcs->SaveBest=TRUE; 
		mcs->SampleColumns(); mcs->RestoreBest();
// mcs->PutHyperPartition(stderr); // exit(1);
fprintf(stderr,"Removing node %d ('Set%d') due to low LLR (%.2f).\n",worst_i,id,worst);
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
		Ndel++; // mcs->PutHyperPartition(stderr);
	   }
	} while(worst < minLLR);
	do {		// 3. Eliminate nodes with small sets.
	   for(Worst=minCard, worst_i=0, i=2; i <= Hpt->NumBPPS(); i++){
              if(Hpt->TypeOfSet(i) != '!') continue;
	      if(stable && MemberSet(Hpt->ItoSetID(i),stable)) continue;
	      card=mcs->SetCard(i); 
	      if(card < Worst){ Worst = card; worst_i=i; }
	   }
	   if(Worst < minCard){
		id=Hpt->ItoSetID(worst_i); mcs->SaveBest=TRUE; 
		// mcs->SampleColumns(); 
		mcs->RestoreBest();
fprintf(stderr,"Removing node %d ('Set%d') due to low #seqs (%d).\n",worst_i,id,Worst);
// mcs->PutHyperPartition(stderr); // exit(1);
#if 0	// DEBUG
		Hpt->Put(stderr);
		fprintf(stderr,"Hpt->ItoSetID(%d)=%d.\n",id,Hpt->SetIDtoI(id));
		assert(Hpt->SetIDtoI(id) > 0);
#endif
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
		Ndel++;
		// mcs->PutHyperPartition(stderr);
	   }
	} while(Worst < minCard);
	if(Ndel > 0){
		// need to resample when an internal node changes to a leaf node; check for this!!!
		mcs_typ	*rtn_mcs = this->RtnCopy( ); assert(rtn_mcs != 0);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
		// this->AdjustSetIDs(); 
		// this->Sample(Temp,2,2,2,2); 
		mcs->ResetTemperature(Temp); mcs->Sample(2,2,2,2); d=mcs->RestoreBest();  CalcLLR(mcs);
		// mcs->PutHyperPartition(stderr); 
		return Ndel; 
	} else return 0;
}

mcs_typ	*omc_typ::DeleteNode(Int4 ID,double T)
{
	if(stable && MemberSet(ID,stable)) DeleteSet(ID,stable);
	return this->Operate('d',0,ID,T); 
}
// Create a mcs_typ with node ('Set<ID>') deleted and its seqs merged into 'parent'.
// then attach the nodes children to the parent.
// Need to make sure that, if a leaf node, parent has other children; otherwise change parent to a leaf node.

//********************************* end DeleteNodes() ******************************************

mcs_typ *omc_typ::InsertNode(Int4 pID, Int4 ID, hsi_typ *hsi,BooLean sample)
// Create a mcs_typ with a node ('Set<ID>') inserted between the row 'parent' and that row's children.
// then move those children whose set ids are in the BG set back up to the parent.
{
	char mode='I'; 
	if(sample) mode='i';
	mcs_typ *rtn_mcs=this->Operate(mode,pID,ID,300,hsi); hpt_typ *hpt=rtn_mcs->GetHpt();
	fprintf(stderr,"\nNode %d ('Set%d') inserted between node %d ('Set%d') and it's children\n",
		hpt->SetIDtoI(ID),ID,hpt->SetIDtoI(pID),pID);
	return rtn_mcs;
}

mcs_typ	*omc_typ::RtnCopy(char mode) { mcs->RestoreBest(); return Operate(mode,0,0); }
// mode == -2 --> use _[A-E].out file name; mode == -3 --> use _bst.out filename.

mcs_typ *omc_typ::Optimize( )
{ mcs_typ *xmcs=this->Operate('O',0,0); xmcs->SaveBest=TRUE; xmcs->StoreBest(); return xmcs; }

e_type	omc_typ::SubGrpCsq(set_typ set)
{
	Int4	i,j,x,s,sq,Sq,**cnts,d,lenM;
	st_type	S=SitesCMSA(TrueMainCMA);
	a_type	AB=AlphabetCMSA(TrueMainCMA);
	BooLean	*skip;
	e_type	csqE=0;

    	assert(nBlksCMSA(TrueMainCMA) == 1);
	NEW(skip,NumSeqsCMSA(TrueMainCMA) + 3, BooLean);
	for(i=1; i <= NumSeqsCMSA(TrueMainCMA); i++) if(!MemberSet(i,set)) skip[i]=TRUE;

	double **freq=ColResFreqsCMSA(1,skip,&cnts,TrueMainCMA);

        Int4 N,N0; N0=NumSeqsCMSA(TrueMainCMA);
        for(N=0,sq=1; sq <= N0; sq++) if(!skip[sq]) N++; 
	assert(N == CardSet(set));
        lenM=LengthCMSA(1,TrueMainCMA);
	char 	r,R,*seq; NEW(seq,lenM +3,char);
#if 0	// Working on this code; not done yet...not sure if I will need it.
	double  best,total,*WtCnt,**SqWt=swt->Weights();
        for(i=1; i <= lenM; i++){
	   NEW(WtCnt,nAlpha(AB) +3, double);
	   for(sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
	   	if(skip[sq]) continue;
	   	s=SitePos(1,sq,1,SitesCMSA(TrueMainCMA)); // start of fake sequence; first blk first repeat.
	   	unsigned char *fsq = SeqSeqSet(sq,DataCMSA(TrueMainCMA));
		r = fsq[s+i-1];
		WtCnt[r] += SqWt[i][sq];
	   } 
	   for(total=0.0,r=0; r <= nAlpha(AB); r++){ total += WtCnt[r]; }
	   double d,bst;
	   for(bst=0,R=0,r=0; r <= nAlpha(AB); r++){
		d=WtCnts[r]/total;
		if(bst < d){ bst=d; R=r; }
	   } if(R == 0) seq[i]='A'; else seq[i]=AlphaChar(R,AB);
	   free(WtCnts);
	} seq[i]=seq[0]=0; csqE=StringToSeq(seq+1,"Consensus seq", 1, AB); 
#else
        for(i=1; i <= lenM; i++){
	   Int4 bst;
	   for(bst=0,R=0,r=0; r <= nAlpha(AB); r++){
		if(bst < cnts[i][r]){ bst=cnts[i][r]; R=r; }
	   } if(R == 0) seq[i]='A'; else seq[i]=AlphaChar(R,AB);
	} seq[i]=seq[0]=0; csqE=StringToSeq(seq+1,"Consensus seq", 1, AB); 
#endif
	// PutSeq(stderr,csqE,AB); 
        for(i=1; i <= lenM; i++){ free(freq[i]); free(cnts[i]); }
	free(freq); free(cnts); free(skip); free(seq);
	return csqE;
}

Int4    IdentAlignSeq(e_type sq1, e_type sq2, a_type AB)
{
        register Int4 i,score=0,len=LenSeq(sq1);
        assert(LenSeq(sq2) ==len);
        unsigned char *seq1 = XSeqPtr(sq1),*seq2=XSeqPtr(sq2);
        char **R=AlphaR(AB);
        for(i=1; i <= len; i++){ score += R[seq1[i]][seq2[i]]; }
        return score;
}

void	omc_typ::Info::Init(Int4 N)
{
	sma=0; set=0; hpt=0; T=0; xsst=0; pttrn=0;
	NEW(sma,N + 3,cma_typ); NEW(set,N + 3,set_typ); NEWP(pttrn,N+3,char);
}

void	omc_typ::Info::Free()
{
	if(this->sma) free(this->sma); this->sma=0; // individual cma files deleted with mcs_typ.
	if(this->set) free(this->set); this->set=0; // individual sets deleted with mcs_typ.
	if(this->pttrn) free(pttrn);
	if(xsst) free(this->xsst); 
	if(this->hpt) this->hpt=0;		// hpt is freed with mcs_typ.
}

mcs_typ *omc_typ::Operate(char action, Int4 pID,Int4 ID, double Temp,hsi_typ *found)
// main routine to perform operations on mcs_typ mcs.
{
	Int4	*P=0,id,sq,i,j,pI,I,J,r,R,m,n,p,x,cID=0;
	double	d,D;
	Info	*info = new Info; info->Init(MaxNumNodesPlus); 
	hpt_typ	*hpt=0;
	// char    **pttrn; NEWP(pttrn,MaxNumNodesPlus+3,char);
	BooLean	Abort=FALSE;

	if(action != 'S' && action != 's' && action != 'u'){ 
		// mcs->PutHyperPartition(stderr); // Hpt->Put(stderr); 
		assert(Hpt->IsTree(P)); 
	} else P=0;
	// if(action == 'u') Hpt->PutSorted(stderr);
	set_typ *set=mcs->CopyOfBestTreeSets(); // calls mcs->CopyOfSeqSets(); if Best not restored.
	// sst_typ **xsst=mcs->RtnCopyOfSSTs();
	info->xsst=mcs->RtnCopyOfSSTs();
	// if(action == 'u') Hpt->PutSorted(stderr);
#if 0	// leave failed nodes in as they can help reach an optimum.
	for(i=1; i < Hpt->NumSets(); i++){
		if(Hpt->TypeOfSet(i) == '!'){
		   if(!((action == 'D' || action == 'd') && Hpt->SetIDtoI(ID) == i)){ 
			if(CardSet(set[i]) <= 0){  // should be eliminated before calling...
				fprintf(stderr,"action = %c; i = %d; ID = %d\n",action,i,ID);
				mcs->PutHyperPartition(stderr);
			} assert(CardSet(set[i]) > 0); 	
		   }	// shouldn't happen unless this is an empty leaf node to be deleted!
		}	
	} 
#endif
	switch (action) {
	   case 'a': {	//******** Use Cont. Tab to add a leaf node to node p == Set<pID>. ************
	   	// Retrieved from /local/projects/aneuwald/JunkDir/omc_operate1_31_13.cc
fprintf(stderr,"**************** action == 'a' ****************\n");
		e_type	bstE,wstE,csqE;
		cma_typ	cma=0;
		if(Hpt->NumBPPS() >= MaxNumNodes) print_error("Maximum # nodes limit reached");
		assert(ID <= MaxNumNodes);
		p=Hpt->SetIDtoI(pID); assert(Hpt->SetIDtoI(ID)==0); assert(ID <= MaxNumNodes); info->T=Temp;
		hpt=Hpt->AddLeaf(p,ID); info->hpt=hpt; 
		r=hpt->SetIDtoI(ID); assert(r != 0);       // r is the new row/column just added.
		fprintf(stderr,"\n==== Attaching leaf node %d ('Set%d') to node %d ('Set%d') =====\n\n",
				r,ID,p,pID);
		for(i=1; i <= Hpt->NumBPPS(); i++){
		   id=Hpt->ItoSetID(i); // assert(id < ID); 
		   sprintf(str,"Set%d",id); j=hpt->SetIDtoI(id); 
		   info->sma[j]=mcs->RtnCsqAsCMSA(i,str,info->xsst[i]); // redefines Csq-incompatible xsst[i].
		   info->set[j]=set[i];
		   info->pttrn[j]=SST2ArgStrAlpha(info->xsst[i],mcs->RtnLengthMainCMSA(),AB);
		   // fprintf(stderr,"%d: %s\n",j,info->pttrn[j]);
		} info->set[hpt->NumSets()] = set[i]; // Random set.
		I = hpt->SetIDtoI(ID); pI = hpt->SetIDtoI(pID); // set the new node's pattern below...
     		info->set[I] = MakeSet(SetN(info->set[1])); ClearSet(info->set[I]);

		//================ find pattern and matching seq. set for new leaf node ===============
		cma=GetInSetCMSA(info->set[pI],TrueMainCMA); // MSA of parent set.
		csqE=MkConsensusCMSA(cma); TotalNilCMSA(cma); // csq of parent set.
		// Use map_typ to find initial pattern & partition.
		Int4	score,sc,top_sc,top_sq,max;
		sst_typ *osst,**sst_rtn=0;    // look at sstList in mad_typ::CliqueClusters();
		Int4    k,NumClust,MinClique=4;	// dont' go below 4.
        	// double  MinKeyFrq=0.4;      // minimum frequency of pattern matches required in cma...
        	double  MinKeyFrq=0.3;      // minimum frequency of pattern matches required in cma...
        	double  MaxGapFrq=0.5;      // maximum number of gap residues (x) allowed.
		char    sets_mode='R',**sst_str;	// 'R' = root
		double  pcut=0.0001;		// for clustering into cliques.
		double	Fisher_pcut=0.001;	// for continency tables.
		double	alpha=0.98;
		AddLeafLLR=0;
	// START WITH MinKeyFrq=0.7 and go down to MinKeyFrq=0.1??? (largest first)!!!???
#if 0
		// switch (random_integer(3)){	// random int from 0..2
		switch (0){	// random int from 0..2
		  case 0: {
			MinKeyFrq=0.5; MaxGapFrq=0.5; pcut=0.01;Fisher_pcut=0.01; } break;
		  case 1: {
			MinKeyFrq=0.2; MaxGapFrq=0.5; pcut=0.0001; Fisher_pcut=0.001; } break;
		  case 2: {
			MinKeyFrq=0.2; MaxGapFrq=0.5; pcut=0.005;Fisher_pcut=0.05; } break;
		  default:
			MinKeyFrq=MP->MinKeyFrq;
			MaxGapFrq=MP->MaxGapFrq; pcut=0.0001;Fisher_pcut=0.001; break;
		}
#endif
                mad_typ *mad=new mad_typ(swt,info->set[pI],TrueMainCMA,MinKeyFrq,MaxGapFrq,sets_mode);
					// MP->MinKeyFrq,MP->MaxGapFrq,MP->sets_mode);
		// mad->Verbose( );
		Fisher_pcut = Fisher_pcut/mad->SearchSpace(); 
		sst_str = mad->CliqueStrings(MinClique,&NumClust,pcut,Fisher_pcut,sst_rtn);
		if(NumClust <= 0){ Abort=TRUE; delete mad; assert(1 || sst_rtn==0); NilSeq(csqE); break; }

		//================ Found one or more patterns at this point ==================
		char	TypNode='L'; // 'm' == Misc mode; 'l' == leaf nodes.
		set_typ setFG=0,setBG=0,SubTreeP=0,setP=0;
		Int4	card=0;
		FILE	*efp=0;

		SubTreeP=MkSubTreeSet(info,pI); // SubTreeP is the parent subtree.
		if(SubTreeP==0) SubTreeP=CopySet(info->set[pI]); // SubTreeP == 0 --> pI is a leaf node...
		setP=CopySet(SubTreeP); setFG=MakeSet(SetN(SubTreeP)); setBG=MakeSet(SetN(SubTreeP));  
		set_typ *tmpSet; NEW(tmpSet,MaxNumNodesPlus+2,set_typ);

		for(osst=0,J=0,D=0.0,j=1; j <= NumClust; j++){
#if 1 //********************** DEBUG **************************
                  fprintf(stderr,"%d: %s\n%d: ",j,sst_str[j],j);
                  // fprintf(stderr,"%d: ",j);
                  for(k=1; k <= mcs->RtnLengthMainCMSA(); k++){
                    if(sst_rtn[j][k]){ PutSST(stderr,sst_rtn[j][k],AB); fprintf(stderr,"%d ",k); }
                  } fprintf(stderr,"\n");
#endif //********************** DEBUG **************************
                  for(n=0,k=1; k <= mcs->RtnLengthMainCMSA(); k++) if(sst_rtn[j][k]) n++;
		  // MaxMisMatches = MAXIMUM(Int4,0,n-(1 + (n/2))); 
		  ClearSet(setFG); CopySet(setP,info->set[pI]); // setP = info->set[pI].
		  tmpSet[I]=setFG; tmpSet[pI]=setP;
		  PatternBasedPartition(tmpSet,sst_rtn[j],pI,I,alpha);
		  IntersectNotSet(SubTreeP,tmpSet[I],setBG); // setBG = SubTreeP && not set[I].
		  sst_typ *zsst=this->GetOptPttrnLPR(efp,tmpSet[I],setBG,d,n,TypNode); 
#if 1 //********************** DEBUG **************************
                  fprintf(stderr,"%d*: ",j);
                  for(k=1; k <= mcs->RtnLengthMainCMSA(); k++){
                    if(zsst[k]){ PutSST(stderr,zsst[k],AB); fprintf(stderr,"%d ",k); }
                  } fprintf(stderr,"\n");
#endif //********************** DEBUG **************************
		  card=CardSet(tmpSet[I]);
		  fprintf(stderr,"FG = %d; BG = %d; LLR = %.2f\n",card,CardSet(setBG),d);
		  if(d > D && card > 1){ D=d; J=j; if(osst) free(osst); osst=zsst; } else free(zsst);
		} // free(sst_str[j]); done by mad_typ.
		AddLeafLLR=D;
		for(j=1; j <= NumClust; j++){ free(sst_rtn[j]); }
		NilSet(setP); NilSet(setFG); NilSet(setBG); NilSet(SubTreeP);  free(tmpSet);
		if(D <= 0.0){ Abort=TRUE; delete mad; NilSeq(csqE); break; }

		PatternBasedPartition(info->set,osst,pI,I,alpha);
		if(CardSet(info->set[I]) == 0){ Abort=TRUE; delete mad; NilSeq(csqE); break; }
		cma=GetInSetCMSA(info->set[I],TrueMainCMA); // MSA of parent set.
		bstE=MkConsensusCMSA(cma); TotalNilCMSA(cma); // csq of parent set.
		// PartitionNode(info->set,pI,I,bstE,csqE);
		for(k=1; k <= mcs->RtnLengthMainCMSA(); k++){
		    for(r=1; r <= nAlpha(AB); r++){
			if(MemSset(r,osst[k])){ EqSeq(k,r,bstE); EqXSeq(k,r,bstE); break; }
		    } if(ResSeq(k,bstE) == 0) EqSeq(k,AlphaCode('A',AB),bstE); // change 'X' to 'A'
		} 
		sprintf(str,"Set%d",ID); cID=ID; // mcs->PutHyperPartition(stderr);
		info->sma[I]=RtnConSqAsCMSA(str,osst,bstE,AB); 
		info->pttrn[I]=SST2ArgStrAlpha(osst,mcs->RtnLengthMainCMSA(),AB); 
#if 0	// DEBUG...
		PutSeq(stderr,bstE,AB); fprintf(stderr,"%s\n",info->pttrn[I]);

setBG=0;SubTreeP=MkSubTreeSet(info,pI); // SubTreeP is the parent subtree.
if(SubTreeP){ setBG=CopySet(SubTreeP); } // SubTreeP == 0 --> pI is a leaf node...
else { setBG=CopySet(info->set[pI]); SubTreeP=info->set[pI]; }
IntersectNotSet(SubTreeP, info->set[I], setBG); // SetBG = SubTreeP && not set[I]

		sst_typ *zsst=this->GetOptPttrnLPR(efp,info->set[I],info->set[pI],d,n,TypNode); 
#if 1 //********************** DEBUG **************************
                  fprintf(stderr,"%d*: ",j);
                  for(k=1; k <= mcs->RtnLengthMainCMSA(); k++){
                    if(zsst[k]){ PutSST(stderr,zsst[k],AB); fprintf(stderr,"%d ",k); }
                  } fprintf(stderr,"\n");
#endif //********************** DEBUG **************************
		free(zsst);
Int4 scC=IdentAlignSeq(csqE,csqE,AB),scE=IdentAlignSeq(bstE,csqE,AB);
Int4 new_card=CardSet(info->set[I]); 
fprintf(stderr,"FG = %d; BG = %d; LLR = %.2f; score = %d vs %d (%.3f)\n",
	new_card,CardSet(setBG),d,scE,scC,(double)scE/(double)scC);
NilSet(setBG);
#if 0	// prior check to speed up...doesn't seem to help...
		if(new_card < (MinimumSetSize - 15)) Abort=TRUE; 
		if(d < (MinimumLLR*0.5)) Abort=TRUE; 
#endif
#if 0
AlnSeqSW(11, 1, bstE,csqE,AB);
PutCMSA(stderr,info->sma[I]); fprintf(stderr,"pttrn = %s\n",info->pttrn[I]);
#endif
#endif
		delete mad; free(sst_rtn); NilSeq(bstE); NilSeq(csqE); free(osst);
fprintf(stderr,"**************** end action == 'a' ****************\n");
	    } break;
	//*******************************************************************************
	   case 'A': {	//********** Add a leaf node to node p == Set<pID>. **************
fprintf(stderr,"**************** action == 'A' ****************\n");
		info->T=Temp;
		if(!AddNewLeaf(info, set, pID,ID)) Abort=TRUE; else cID=ID;
fprintf(stderr,"**************** end action == 'A' ****************\n");
	//*******************************************************************************
	    } break;
	   case 'r': 	// copy of best seq in each set as .
	   case 'b': {	// copy of best seq in each set as .
		assert(Hpt->NumSets() <=3);
		hpt=Hpt->Copy();  info->hpt=hpt;  info->T=0; 
mcs->PutHyperPartition(stderr);
		id=Hpt->ItoSetID(1); sprintf(str,"Set%d",id);
		info->sma[1]=mcs->RtnCsqAsCMSA(1,str,info->xsst[1]); 
		info->set[1]=set[1];

		id=Hpt->ItoSetID(2); sprintf(str,"Set%d",id);
		if(action=='b') info->sma[2]=mcs->RtnBstAsCMSA(2,str,info->xsst[2]);  // may redefine xsst[i].
		else {	// pick a seed sequence at random...
			Int4 r,rx,s,card=CardSet(set[2]); assert(card>0);
			r=random_integer(card)+1;
			for(rx=0,s=1; s <= NumSeqsCMSA(TrueMainCMA); s++){
			    if(MemberSet(s,set[2])) rx++;
			    if(rx == r) break;
			} assert(s <= NumSeqsCMSA(TrueMainCMA));
			info->sma[2]=mcs->RtnSeqAsCMSA(s,str,info->xsst[2]);
#if 0	// debug...
			for(Int4 rs=1; rs <= LengthCMSA(1,TrueMainCMA); rs++)
			   if(info->xsst[2][rs]) PutSST2Alpha(stderr,info->xsst[2][rs],AB);
#endif
		}
		info->set[2]=set[2];
		for(Int4 x=1; x <= LengthCMSA(1,TrueMainCMA); x++){
		   if(info->xsst[1][x] == info->xsst[2][x]){ info->xsst[2][x]=0; } // xsst[1][x]=0;
#if 0
PutSST2Alpha(stderr,info->xsst[1][x],AB);
fprintf(stderr,": xsst[1][%d]=%d\n",x,info->xsst[1][x]);
PutSST2Alpha(stderr,info->xsst[2][x],AB);
fprintf(stderr,": xsst[2][%d]=%d\n",x,info->xsst[2][x]);
#endif
		} info->set[hpt->NumSets()] = set[3]; // info->set=set; set=0;

		info->pttrn[1]=SST2ArgStrAlpha(info->xsst[1],mcs->RtnLengthMainCMSA(),AB);
		info->pttrn[2]=SST2ArgStrAlpha(info->xsst[2],mcs->RtnLengthMainCMSA(),AB);
#if 0
fprintf(stderr,"\nAlphaCode('C',AB) =%d; sst = %d\n",AlphaCode('C',AB),SsetLet(AlphaCode('C',AB)));
fprintf(stderr,"\nAlphaCode('X',AB) =%d; sst = %d\n",AlphaCode('X',AB),SsetLet(AlphaCode('X',AB)));
#endif

if(info->pttrn[1]) fprintf(stderr,"\npttrn[1]=%s\n",info->pttrn[1]);
if(info->pttrn[2]) fprintf(stderr,"\npttrn[2]=%s\n",info->pttrn[2]);
	    } break;
	   case 'B':	// copy for best_mcs.
	   case 'c':	// copy of mcs_typ but with sequences partitioned anew.
	   case 'C': {	// Exact copy of mcs_typ; ID indicates mode!!!
		hpt=Hpt->Copy();  info->hpt=hpt;  info->T=0; 
		for(j=1; j < hpt->NumSets(); j++){
		   id=Hpt->ItoSetID(j); sprintf(str,"Set%d",id);
		   // info->sma[j]=mcs->RtnCsqAsCMSA(j,str,info->xsst[j]);  // ensures xsst[j] compatible w/ seq.
	           // info->sma[j]=mcs->RtnCsqSstAsCMSA(j,str,info->xsst[j]);  // ensures seq. compatible w/ xsst[i]
		   info->sma[j]=mcs->RtnCsqAsCMSA(j,str);  // won't redefine xsst[i].
		   info->pttrn[j]=SST2ArgStrAlpha(info->xsst[j],mcs->RtnLengthMainCMSA(),AB);
		   info->set[j]=set[j];
		} info->set[hpt->NumSets()] = set[j]; // info->set=set; set=0;
	    } break;
	   case 'd':	// simple column sampling only afterwards...
	   case 'D': {	//*********** Delete node with ID. *****************
		R=Hpt->SetIDtoI(ID); assert(R != 0); info->T=Temp;      // R is the new row/column to be deleted.
		p=P[R]; assert(p != 0); pID=Hpt->ItoSetID(p); assert(pID != 0);
		fprintf(stderr,"  Deleting row %d ('Set%d') with parent row %d ('Set%d')\n\n",R,ID,p,pID);
		hpt=Hpt->Delete(R); info->hpt=hpt;
	        for(i=1; i < Hpt->NumSets(); i++){
	          id=Hpt->ItoSetID(i);
	          if(i > R) j=i-1; else if(i < R) j=i;
	          else { // i == R...
	                if(info->pttrn[i]) free(info->pttrn[i]); pI=Hpt->SetIDtoI(pID);
	                UnionSet(set[pI],set[i]); // Add deleted node's set to its parent's set.
	                NilSet(set[i]); continue;
	          } info->set[j]=set[i]; sprintf(str,"Set%d",id); 
	          info->sma[j]=mcs->RtnCsqAsCMSA(i,str,info->xsst[i]);    // may redefine xsst[i].
	          info->pttrn[j]=SST2ArgStrAlpha(info->xsst[i],mcs->RtnLengthMainCMSA(),AB);
     		} info->set[i-1]=set[i]; pI=hpt->SetIDtoI(pID); n=hpt->NumDescendants(pI);
	        if(n == 0 && hpt->TypeOfSet(pI) == '?') 
		    { hpt->ChangeTypeOfSet(pI,'!'); assert(hpt->SetIDtoI(pID) == pI); }
		// WARNING: internal --> leaf: -A=40:10 --> -A=48:2 -Ri=0.4 --> -Ri=0.003.
		fprintf(stderr,"node Set%d deleted from Set%d (%d)\n",ID,pID,hpt->SetIDtoI(pID));
	    } break;
	   case 'F': {	// Fuse nodes. 
#if 0
		// Move nodes down and then delete node; or perhaps call these instead.
#endif
			print_error("FuseNodes() not yet implemented");
		} break;
	   case 'i':
		cID=ID;
	   case 'I': {	// Insert internal node. 
		// Create a mcs_typ with a node ('Set<ID>') inserted between the row 'pID' and that 
		// row's children.  Then move those children whose set ids are in the BG set back up 
		// to the parent.  
		set_typ BG=found->RtnSetBG(),FG=found->RtnSetFG();;
		assert(Hpt->SetIDtoI(ID) == 0);
		p=Hpt->SetIDtoI(pID); hpt=Hpt->Insert(p,ID); // PutSet(stderr,FG); PutSet(stderr,BG);
		R=hpt->SetIDtoI(ID); assert(R != 0); info->T=Temp;      // R == new row/column just added.
		// fprintf(stderr,"Inserted row %d ('Set%d') after parent row %d ('Set%d')\n",R,ID,p,pID);
		for(i=2; i <= hpt->NumBPPS(); i++){  // move children of Set<ID> in BG set back to grandparent
		   if(i == R) continue;
	           if(hpt->Cell(i,R) == '+'){
		     id=hpt->ItoSetID(i); I=Hpt->SetIDtoI(id); 
	             if(I > 0 && Hpt->ItoSetID(P[I]) == pID && MemberSet(id,BG)) {
	                hpt->MoveNodeUp(i); // hpt->Put(stderr); // debug...
	             }
	           }
		}
		hpt_typ *tmpHpt=hpt->Sort(); delete hpt; hpt=tmpHpt; info->hpt=hpt;
		// hpt->Put(stderr);
		sst_typ	*fsst=found->GetSST();
		for(i=1; i <= Hpt->NumBPPS(); i++){
	          id=Hpt->ItoSetID(i); 
	          sprintf(str,"Set%d",id);       // modify sma by adding internal node consensus sequence here.
	          j=hpt->SetIDtoI(id);
		  if(fsst==0){
	             info->sma[j]=mcs->RtnCsqAsCMSA(i,str,info->xsst[i]);  // redefines incompatible xsst[i]
		  } else {
	             info->sma[j]=mcs->RtnCsqSstAsCMSA(i,str,info->xsst[i]);  // makes seq. compatible w/ xsst[i]
		  }
	          info->set[j]=set[i];
		  assert(SetN(set[1]) == SetN(set[i]));
#if 0		  // original code...modified for InsertNodes...
		  if(!MemberSet(id,FG)) info->pttrn[j]=SST2ArgStrAlpha(info->xsst[i],mcs->RtnLengthMainCMSA(),AB);
#else
		  info->pttrn[j]=SST2ArgStrAlpha(info->xsst[i],mcs->RtnLengthMainCMSA(),AB);
	          // fprintf(stderr,"%d ('Set%d'): %s\n",i,id,info->pttrn[j]);
#endif
		} info->set[hpt->NumSets()] = set[Hpt->NumSets()]; // Random set.
		sprintf(str,"Set%d",ID); j=hpt->SetIDtoI(ID); assert(info->set[j] == 0);
#if 0
		info->sma[j]= mcs->RtnCsqAsCMSA(p,str);   // make the new node csq the same as it's parent.
#else
		if(fsst==0){ // use for ResurrectRejects();
			info->sma[j]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA(str,info->sma[j]);
			e_type csqE=MkConsensusCMSA(TrueMainCMA);
			fsst=SST_FromSeq(csqE); NilSeq(csqE);
		} else {
		        info->sma[j]= mcs->RtnCsqSstAsCMSA(p,str,fsst);   
// e_type csqE = TrueSeqCMSA(1,info->sma[j]); PutSeq(stderr,csqE,AB);
			info->pttrn[j]=SST2ArgStrAlpha(fsst,mcs->RtnLengthMainCMSA(),AB);
		}
		
#endif
		info->set[j] = MakeSet(SetN(set[1])); ClearSet(info->set[j]); // leave new set empty..
		assert(SetN(set[1]) == SetN(info->set[j])); 
		// hpt->Put(stderr);
	    } break;
	   case 'O': {	// Optimize the display set.
		FILE	*fp=0; info->T=100;
		hpt=Hpt->Copy();  info->hpt=hpt; 
	        for(i=1; i < hpt->NumSets(); i++){
		  Int4	hits=0;
		  cma_typ setcma=0,tmp_sma[4],cma=TrueMainCMA; 

		  id=hpt->ItoSetID(i); // assert(id <= hpt->NumBPPS());
	          info->pttrn[i]=SST2ArgStrAlpha(info->xsst[i],mcs->RtnLengthMainCMSA(),AB);

		  sprintf(str,"Set%d",id); tmp_sma[1]=mcs->RtnCsqAsCMSA(i,str);
		  if(CardSet(set[i]) > 0){
		    fp=tmpfile(); PutInSetCMSA(fp,set[i],cma); rewind(fp);
		    setcma=ReadCMSA(fp,AB); fclose(fp); fp=tmpfile();
		    hits=PutBestRepsCMSA(fp,10,TRUE,setcma,TRUE,tmp_sma[1]);
		  } else { fp=0; hits=0; }
         	  if(hits > 0){
                    rewind(fp); tmp_sma[2]=ReadCMSA(fp,AB); fclose(fp);
                    fp=tmpfile(); PutMergedCMSA(fp,2,tmp_sma); rewind(fp); info->sma[i] = ReadCMSA(fp,AB);
                    TotalNilCMSA(tmp_sma[1]); TotalNilCMSA(tmp_sma[2]);
                    sprintf(str,"Set%d",id);  // modify sma by adding internal node consensus sequence here.
                    RenameCMSA(str,info->sma[i]); fclose(fp);
         	  } else { if(fp) fclose(fp); info->sma[i]=tmp_sma[1]; }
		  if(setcma) TotalNilCMSA(setcma);

		  info->set[i]=set[i];

		} info->set[hpt->NumSets()] = set[i]; // info->set=set; set=0;
	    } break;
	   case 'R': {	// Randomize sequence assignments; 1/pID of the seq. assignments are randomized.
		assert(pID > 0); // This operations uses pID for fraction of assignment to randomize.
		hpt=Hpt->Copy();  info->hpt=hpt;  info->T=0; 
		for(j=1; j < hpt->NumSets(); j++){
		   id=Hpt->ItoSetID(j); sprintf(str,"Set%d",id);
		   info->sma[j]=mcs->RtnCsqAsCMSA(j,str,info->xsst[j]);  // may redefine xsst[i].
		   info->pttrn[j]=SST2ArgStrAlpha(info->xsst[j],mcs->RtnLengthMainCMSA(),AB);
		   info->set[j]=set[j];
		} info->set[hpt->NumSets()] = set[j]; // info->set=set; set=0;
		set_typ SetR=MakeSet(SetN(set[1])); ClearSet(SetR);
		for(i=1; i <= NumSeqsCMSA(TrueMainCMA); i++){
		   if(random_integer(pID) != 0) continue;
		   for(j=1; j <= hpt->NumSets(); j++){
			if(MemberSet(i,set[j])){ DeleteSet(i,set[j]); break; }
		   } x=random_integer(hpt->NumSets()) + 1; AddSet(i,set[x]);
		}
	    } break;
	   case 'u':
	   case 's':
	   case 'S': {	// Sort Hpt based on set IDs; checks that Hpt is a tree too.
		hpt=Hpt->Sort();  info->hpt=hpt;  info->T=Temp; 
		for(j=1; j < hpt->NumSets(); j++){
		   id=hpt->ItoSetID(j); assert(id > 0); I=Hpt->SetIDtoI(id);
		   // id=Hpt->ItoSetID(j); 
		   sprintf(str,"Set%d",id);
		   info->sma[j]=mcs->RtnCsqAsCMSA(I,str,info->xsst[I]);
		   info->pttrn[j]=SST2ArgStrAlpha(info->xsst[I],mcs->RtnLengthMainCMSA(),AB);
		   info->set[j]=set[I];
		} info->set[hpt->NumSets()] = set[Hpt->NumSets()]; // info->set=set; set=0;
	    } break;
	   case 'U': {	// move node 'ID' one level up the hiearchy.
	    } break;
	   default: 
		// assert(isalpha(action));
		fprintf(stderr,"action = '%c'(%d)\n",action,action);
		print_error("Fatal: omc_typ::Operate( ) input 'action' error"); 
		break;
	}
	hpt=info->hpt;
	if(hpt){ // if Aborted then hpt=0;
	  // assert(hpt == info->hpt);
     	  for(i=1; i <= hpt->NumBPPS(); i++){
            if(info->pttrn[i]){
		char *tmp[3]; tmp[0]=info->pttrn[i]; hpt->ReSetArgv(i,1,tmp); free(info->pttrn[i]); 
		id=hpt->ItoSetID(i); sprintf(str,"Set%d",id); hpt->ReNameGroup(i,str);
	    }
	  }
     	} 
 // CompareInput(stderr,info->set,info->hpt,info->sma);  // debug..
	mcs_typ *rtn_mcs=0;
	if(!Abort) rtn_mcs=CreateMCS(info->hpt,info->set,info->sma,info->T,action,pID,cID);
	else if(info->hpt){ 
	  if(info->sma){
		for(i=1; i < info->hpt->NumSets(); i++){
		   if(info->sma[i]) TotalNilCMSA(info->sma[i]); info->sma[i] = 0;
		} 
	  } if(info->set){
		for(i=1; i <= info->hpt->NumSets(); i++){
		   if(info->set[i]) NilSet(info->set[i]); info->set[i] = 0;
		}
	  } delete info->hpt; 
	}
	for(i=1; i <= Hpt->NumBPPS(); i++){ if(info->xsst[i]) free(info->xsst[i]); } 
	info->Free(); delete info; // rtn_mcs->PutHyperPartition(stderr); 
	if(set) free(set); 
	if(P) free(P);
	return rtn_mcs;
}

