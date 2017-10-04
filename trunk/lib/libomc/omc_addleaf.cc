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


Int4	omc_typ::AddFocusedLeaf(double Temp)
{
// cerr << "debug 1\n";
        Int4	x,I,ID,i,j,id,n,card,NumAdded=0;
        double d,D,*map = mcs->RtnSubLPR();      // map[0] == total LPR.
	assert(FocusNode > 0 && FocusNode < Hpt->NumSets());
	// assert(FocusSeq > 0 && FocusSeq <= NumSeqs);
	if(Hpt->NodeDepth(FocusNode) < this->MaxDepthHierarchy){
	      id=Hpt->ItoSetID(FocusNode); 
	      if(stable && MemberSet(id,stable)) print_error("AddFocusedLeaf() input error");
	      card=mcs->SetCard(FocusNode);
	      fprintf(stderr,"%d: id=%d; card = %d\n",i,id,card);
              if(card < MinimumSplitSize) print_error("AddFocusedLeaf() input error");
        } else return 0;
	{
	    n=Hpt->SetIDtoI(id); 
	    if(n == 0) return 0;  // this node has been deleted...!!!
            // fprintf(stderr,"======= Sampling a leaf attached to node %d (Set%d). =======\n",n,id);
            for(ID=1; ID <= Hpt->NumSets(); ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID:
	    assert(ID <= Hpt->NumSets());
	    mcs_typ *rtn_mcs=this->AddFocusedLeaf(id,ID,Temp); 
	    if(rtn_mcs == 0){ fprintf(stderr," ------Set%d' aborted. ------\n",ID); }
	    else {
	      hpt_typ *hpt=rtn_mcs->GetHpt(); i=hpt->SetIDtoI(ID); CalcLLR(rtn_mcs);
              double *Map=rtn_mcs->RtnSubLPR( ),D=mcs->RtnBestLPR(); 
	      fprintf(stderr,
		"LLR=%.2f --> %.2f; subLLR=%.2f (%.2f min); node %d('Set%d %c'); %d seqs (%d min) %.1f npws\n",
                  D,Map[0],Map[i],MinimumLLR,i,ID,hpt->TypeOfSet(i),
		  rtn_mcs->SetCard(i),MinimumSetSize, rtn_mcs->RtnNatsPerWtSeq(i));
	      if((Map[0] > D && Map[i] > MinimumLLR && (hpt->TypeOfSet(i) == '?'
				|| rtn_mcs->SetCard(i) >= MinimumSetSize))
				&& rtn_mcs->OkayNatsPerWtSq(i)){
		// rtn_mcs->PutHyperPartition(stderr);  
		fprintf(stderr,
		   "!!!!!!!!!!!!!!!!!! Sampling node %d ('Set%d') sucessful. !!!!!!!!!!!!!!!!!!!!\n",i,ID);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); this->RestoreBest();
		NumAdded++;  // if(NumAdded % 8 == 0) mcs->Sample(2,2,2,2,mcs->TargetMod);
	        i=Hpt->SetIDtoI(id); card=mcs->SetCard(i);
#if 0	// This may not be changing anything...
		this->ReSetUpFocusedSrch(i); // FocusNode=i; FocusSeq=0; // to speed up sampling...
#endif
	      } else DeleteMCS(rtn_mcs);  
	    }
// cerr << "debug 4\n";
	} return NumAdded;
}

BooLean omc_typ::ReSetUpFocusedSrch(Int4 NewI)
// set stable set so that only sequences assigned to Rejected, FocusNode and NewI node get reassigned.
// sample between parent and child nodes only; stable should not include parent and new child.
{
        Int4 i,id;
        if(FocusNode == 0) return FALSE;
        else {
            assert(Hpt);
            if(FocusNode < 1 || FocusNode >= Hpt->NumSets()) print_error("ReSetUpFocusedSrch() error");
            if(FocusNode <= 1  && FocusSeq == 0) print_error("ReSetUpFocusedSrch() error");
            set_typ subtree=Hpt->MkSubTreeSet(FocusNode);
            for(i=1; i < Hpt->NumSets(); i++){
                if(i == FocusNode) continue;
                else if(i == NewI) continue;
                else if(MemberSet(i,subtree)){
		   id=Hpt->ItoSetID(i); AddSet(id,stable);
                } 
            } // Hpt->Put(stderr); PutSet(stderr,subtree); PutSet(stderr,stable);
            NilSet(subtree);
        } // exit(1);
        return TRUE;
}

mcs_typ *omc_typ::AddFocusedLeaf(Int4 pID, Int4 ID,double T)
{
	double Temp=T;
	Int4	*P=0,id,i,j,I,r,R,m,n,p,x,cID=0;
	double	d,D;
	Info	*info = new Info; info->Init(MaxNumNodesPlus); 
	hpt_typ	*hpt=0;
	BooLean	Abort=FALSE;

	assert(Hpt->IsTree(P)); 
	set_typ *set=mcs->CopyOfBestTreeSets(); // calls mcs->CopyOfSeqSets(); if Best not restored.
	info->xsst=mcs->RtnCopyOfSSTs(); info->T=Temp;
	if(!AddNewLeaf(info, set, pID,ID)) Abort=TRUE; else cID=ID;
	//*******************************************************************************
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
	mcs_typ *rtn_mcs=0;
	if(!Abort) rtn_mcs=CreateMCS(info->hpt,info->set,info->sma,info->T,'A',pID,cID);
	else if(info->hpt){ 
	  if(info->sma){
		for(i=1; i < info->hpt->NumSets(); i++)
		   { if(info->sma[i]) TotalNilCMSA(info->sma[i]); info->sma[i] = 0; } 
	  } if(info->set){
		for(i=1; i <= info->hpt->NumSets(); i++)
		   { if(info->set[i]) NilSet(info->set[i]); info->set[i] = 0; }
	  } delete info->hpt; 
	}
	for(i=1; i <= Hpt->NumBPPS(); i++){ if(info->xsst[i]) free(info->xsst[i]); } 
	info->Free(); delete info; // rtn_mcs->PutHyperPartition(stderr); 
	if(set) free(set); if(P) free(P);
	return rtn_mcs;
}

BooLean	omc_typ::GrowFocusedLeaf(FILE *fp,double Temp, Int4 iter)
{
	Int4 i=Hpt->NumSets();
	if(FocusNode == 0 || FocusSeq == 0){ print_error("GrowFocusedLeaf() input error 1"); } 
#if 1	// move only FocusSeq to root node.
	if(FocusNode == 1 && mcs->ResurrectRejectSeq(FocusSeq) == FALSE){  
	   // true if sq is in root or rejected .
	   print_error("GrowFocusedLeaf() input error 2"); // will only 
	}
#else	// move all rejects to root node.
	if(FocusNode == 1 && mcs->ResurrectRejectSeq(0) == FALSE){  // true if sq is in root or rejected .
	   print_error("GrowFocusedLeaf() input error 2"); // will only 
	}
#endif
        if((status=this->AddFocusedLeaf(Temp)) > 0){
                if(fp) fprintf(fp,"GrowFocusedLeaf() added %d nodes.\n",status);
#if 1
                Sample(Temp,2,2,2,2,mcs,0.03);
#else
		
#endif
                this->RestoreBest(); // mcs->PutHyperPartition();
                // return TrimNodes(fp,300,iter);
        	// return this->CalcLLR(mcs);
                // return this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
		return TRUE;
        } else return FALSE; // return this->CalcLLR(mcs);
}

BooLean	omc_typ::AddNewLeaf(Info *info, set_typ *set, Int4 pID,Int4 ID)
// main routine to perform operations on mcs_typ mcs.
{
	Int4	id,i,j,k,pI,I,J,r,R,m,n,p,x,sq,NumClust,len=mcs->RtnLengthMainCMSA();
	Int4    NumCols=mcs->RtnDefaultMaxCol();
	hpt_typ	*hpt=0;
	BooLean	Success=TRUE;
	double d,D,dd,bstD=-999999999999.9;
	e_type	tmpE,TmpE,bstE,wstE,csqE,CsqE,bstCsqE;
	sst_typ *osst=0;
	mad_typ *mad=0;

	// mcs->PutHyperPartition(stderr); // Hpt->Put(stderr); 
	//=============== Check input values =======================
	if(Hpt->NumBPPS() >= MaxNumNodes) print_error("Maximum # nodes limit reached");
	assert(ID <= MaxNumNodes);
	p=Hpt->SetIDtoI(pID); assert(Hpt->SetIDtoI(ID)==0); assert(ID <= MaxNumNodes); 
	if(CardSet(set[p]) == 0){ return FALSE; }
	hpt=Hpt->AddLeaf(p,ID); info->hpt=hpt; 
	r=hpt->SetIDtoI(ID); assert(r != 0);       // r is the new row/column just added.

	fprintf(stderr,"\n==== Attaching leaf node %d ('Set%d') to node %d ('Set%d') =====\n\n",
			r,ID,p,pID);

	//====================================================================================
	//=================== Create a new hpt by adding the leaf to pID node =================
	//====================================================================================
	for(i=1; i <= Hpt->NumBPPS(); i++){
	   id=Hpt->ItoSetID(i); j=hpt->SetIDtoI(id); 	// j = new node index for old node i.
	   sprintf(str,"Set%d",id); 
	   info->sma[j]=mcs->RtnCsqAsCMSA(i,str,info->xsst[i]); // redefines Csq-incompatible xsst[i].
	   // if(i==1) info->sma[j]=MakeConsensusCMSA(TrueMainCMSA);
	   // else info->sma[j]=mcs->RtnBstAsCMSA(i,str,info->xsst[i]); // redefines incompatible xsst[i].
	   info->set[j]=set[i];
	   info->pttrn[j]=SST2ArgStrAlpha(info->xsst[i],len,AB);
	   // fprintf(stderr,"%d: %s\n",j,info->pttrn[j]);
	} info->set[hpt->NumSets()] = set[i]; // Random set.
	I = hpt->SetIDtoI(ID); pI = hpt->SetIDtoI(pID); // set the new node's pattern below...
     		info->set[I] = MakeSet(SetN(info->set[1])); ClearSet(info->set[I]);

	// info->PrintSetSizes(stderr);
	// cma_typ cma=GetInSetCMSA(info->set[pI],TrueMainCMA); // MSA of parent set.
	//************* ^ Creating this cma_typ uses up too much memory...!!! *****************

	//====================================================================================
	//=================== Load sequences to be sampled onto a heap =======================
	//====================================================================================
	AddLeafLLR=0;	// this is a class variable 
	sst_typ *bsst=0,*zsst;
	set_typ setBG=0,setP=MkSubTreeSet(info,pI); // setP is the parent subtree.
	if(setP){ setBG=CopySet(setP); } // setP == 0 --> pI is a leaf node...
	h_type HG=Histogram("preliminary LLR increase in nats",0,1000,40);
	dh_type dH=dheap(NumSeqsCMSA(TrueMainCMA)+2,4);
#if 0
	for(sq=1; sq <= NumSeqsCMSA(TrueMainCMA) ; sq++){
		if(MemberSet(sq,info->set[pI])) insrtHeap(sq,(keytyp)Random(),dH); 
	}
#else
	if(FocusNode > 0 && FocusSeq > 0){
	     if(!MemberSet(FocusSeq,info->set[pI])){
		fprintf(stderr,"FocusNode = %d; FocusSeq = %d\n",FocusNode,FocusSeq);
		print_error("Focus sequence not assigned to focus node");
	     } insrtHeap(FocusSeq,(keytyp)Random(),dH);
	} else for(sq=1; sq <= NumSeqsCMSA(TrueMainCMA) ; sq++){
		if(MemberSet(sq,info->set[pI])) insrtHeap(sq,(keytyp)Random(),dH); 
	}
#endif
        mad=new mad_typ(swt,info->set[pI],TrueMainCMA,MP->MinKeyFrq,MP->MaxGapFrq,MP->sets_mode);
	// NumClust=MINIMUM(Int4,50,ItemsInHeap(dH)); // Randomly select up to 50 sequences as seeds.
	NumClust=MinimumSetSize; 	// MinimumSetSize >= ItemsInHeap(dH)/tries; 
	NumClust=MINIMUM(Int4,NumClust,ItemsInHeap(dH)); // Make sure 10 isn't > heapsize.
	// NumClust=ItemsInHeap(dH)/3;
#if 1
	// Int4 MinSize=MinimumSetSize-2;
	Int4 MinSize=(Int4) ceil((double)MinimumSetSize/10.0);
	MinSize = MAXIMUM(Int4,MinSize,5);
	Int4 MinCol=mcs->RtnDefaultMinCol();
#else
	Int4 MinSize=(Int4) ceil((double)MinimumSetSize/2.0);
	MinSize = MAXIMUM(Int4,MinSize,5);
#endif
	//===================================================================================
	//======================== Sample sequences & partitions  ===========================
	//===================================================================================
	FILE *efp=0; // efp=stderr;
	Int4 IbstN=0,pIbstN=0,ItmpN=0,pItmpN=0,tmpC=0,TmpC=0,bstCol=0;
	for(D=d=0.0,bstE=0,bstCsqE=0,j=1; j <= NumClust; j++){
	  sq=delminHeap(dH); 
	  tmpE=GetSeqAsCsqCMSA(sq,TrueMainCMA); 
#if 0	// debug...
	  osst=SST_FromSeq(tmpE); // same as sequence.
	  if(j > 0){ n=PrintPttrn(stderr,osst); n=MAXIMUM(Int4,n,30); }
	  else { n=30; fprintf(stderr," (using consensus as pattern).\n"); }
	  free(osst);
#endif
	  csqE=CopySeq(mad->CnsSeq());
	  // if(setP){ csqE=SubGrpCsq(setP); } else csqE=CopySeq(mad->CnsSeq());
	  osst=0; CsqE=0; D=0.0; d=0.0; bstCol=0;
	  char TypNode='L'; // 'L' == mode for leaf nodes.
	  double pRi=0.25;
	  e_type useE=0;	// or set to tmpE to constrain pattern.
	  Int4 MaxIter=2;	// WARNING: don't set to 3 unless change the routine below.
	  // if(FocusSeq > 0) MaxIter=3;
	  for(Int4 Iter=1; Iter <= MaxIter; Iter++){
	    PartitionNode(info->set,pI,I,tmpE,csqE); 
	    // Allow for greater contamination early on.... look at bpps_typ::Initialize().
	    if(FocusSeq > 0){	
		if(Iter==1){ TypNode='R'; pRi=0.25; useE=0; }
		else if(Iter==2){ TypNode='I'; pRi=0.10; useE=0; }
		else { TypNode='m'; pRi=0.05; useE=0; }
	    } else { 
		if(Iter==1){ TypNode='I'; pRi=0.10; useE=0; }
	        else if(Iter==2){ TypNode='m'; pRi=0.05; useE=0; }
	    	else { TypNode='l'; pRi=0.01; useE=0; }
	    }
	    if(setP){
	      IntersectNotSet(setP, info->set[I], setBG); // SetBG = setP && not set[I]
	      zsst=this->GetOptPttrnLPR(efp,info->set[I],setBG,d,NumCols,TypNode,useE,pRi); 
	      if(FocusSeq > 0){
		fprintf(stderr,"FG = %d; BG = %d; LLR = %.2f\n",CardSet(info->set[I]),CardSet(setBG),d);
	      }
	    } else zsst=this->GetOptPttrnLPR(efp,info->set[I],info->set[pI],d,NumCols,TypNode,useE,pRi); 
	    for(tmpC=0,k=1; k <= len; k++) if(zsst[k]) tmpC++;
	    IncdHist(d,HG);
#if 1	// DEBUG...
	    if(tmpC > 0){
		if(FocusSeq > 0){
		   for(Int4 r=1; r <= LenSeq(tmpE); r++){
			if(zsst[r]){ PutSST(stderr,zsst[r],AB); fprintf(stderr,"%d ",r); }
		   } fprintf(stderr," (%d pttrn pos.)\n",tmpC);
		}
	    } else d = -9999999.0; // ignore these...
	    // MinCol=3;
#endif
	    if(Iter==2) break;

 	    osst=zsst; D=d; TmpE=tmpE; CsqE=csqE; TmpC=tmpC;
	    ItmpN=CardSet(info->set[I]); pItmpN=CardSet(info->set[pI]); 
#if 0
	    if(FocusSeq > 0){	// skip this for now...
	        info->Move(pI,I);  // move subnode seqs I back to parent node.
	        e_type *rtnE=FindDivergentSeqs(info,pI,I,zsst);	// 
	        csqE=rtnE[0]; tmpE=rtnE[1]; free(rtnE);
		continue;
	    }
#endif
#if 0	    // Find seed sequence based on matches to the pattern from first iteration.
	    if(setP){ csqE=SubGrpCsq(setP); } else csqE=CopySeq(mad->CnsSeq());
	    info->Move(pI,I);  // move subnode seqs I back to parent node.
	    tmpE=FindSeedSeq(info,pI,I,zsst);	// 
#else 	    // FindDivergentSeqs( ) routine below...Right now this appears to work best.
	    info->Move(pI,I);  // move subnode seqs I back to parent node.
	    // e_type *rtnE=FindDivergentSeqs(info,pI,I,zsst,stderr);	// 
	    e_type *rtnE=FindDivergentSeqs(info,pI,I,zsst);	// 
	    csqE=rtnE[0]; tmpE=rtnE[1]; free(rtnE);
	    // NilSeq(csqE); csqE=SubGrpCsq(info->set[pI]);
#endif
	  }  // D > d --> Iter==1 (found seed) is better than Iter==2 (SubGrp Csq).
	  Int4 bstIter=0; dd=-99999999.0;
	  if(D > d && ItmpN >= MinSize && TmpC >= MinCol){	// only save if leaf size >= MinSize...
	    bstIter=1; NilSeq(tmpE); tmpE=TmpE; free(zsst); dd=D; NilSeq(csqE); csqE=CsqE; 
	    bstCol=TmpC;
	  } else {	// then use iter 2 is better; free iter 1 data.
	    ItmpN=CardSet(info->set[I]); 
	    if(ItmpN >= MinSize && tmpC >= MinCol){
	       bstIter=2; free(osst); osst=zsst; NilSeq(TmpE); NilSeq(CsqE); dd=d; 
	       pItmpN=CardSet(info->set[pI]); bstCol=tmpC;
	    } else { free(zsst); NilSeq(TmpE); NilSeq(CsqE); }
	  }
#if 0
	  if(bstIter > 0) fprintf(stderr,"  ----- %.2f nats; %d cols; Child=%d >= %d; Parent=%d (%d) -----\n",
                                dd,bstCol,ItmpN,MinSize,pItmpN,bstIter);
#endif
	  // if(dd >= MinimumPreLLR && dd > bstD && ItmpN >= MinSize)
	  if(bstIter > 0 && dd > bstD && ItmpN >= MinSize)
	  {
		// if(CsqE) AlnSeqSW(12, 4,csqE,CsqE,AB);
		if(bsst!=0) free(bsst); if(bstCsqE) NilSeq(bstCsqE); if(bstE) NilSeq(bstE);
	        IbstN=ItmpN; pIbstN=pItmpN;
		bstCsqE=csqE; bstD=dd; bsst=osst; bstE=tmpE;
	        fprintf(stderr,"  ===== %.2f nats; %d cols; Child=%d >= %d; Parent=%d (%d) =====\n",
                                bstD,bstCol,IbstN,MinSize,pIbstN,bstIter);
	  } else { free(osst); NilSeq(tmpE); NilSeq(csqE); }
	  info->Move(pI,I);		// move subnode seqs I back to parent node.
	} if(setP) NilSet(setP); if(setBG) NilSet(setBG);

	//===================================================================================
	//=========================== Recompute the best found ==============================
	//===================================================================================
	PutHist(stderr,60,HG); 
	if(bstE==0){ Success=FALSE; }
	else {
	  // info->Move(pI,I);		// move subnode seqs I back to parent node.
	  AddLeafLLR=bstD;
	  PartitionNode(info->set,pI,I,bstE,bstCsqE);  // reset this using bstE;
	  sprintf(str,"Set%d",ID); // mcs->PutHyperPartition(stderr);
	  info->sma[I]=RtnConSqAsCMSA(str,bsst,bstE,AB); 
	  info->pttrn[I]=SST2ArgStrAlpha(bsst,len,AB); 
#if 0	// DEBUG
{
 Int4	k,len=LengthCMSA(1,TrueMainCMA),r;
 e_type	xE=FakeSeqCMSA(1,info->sma[I]);
 for(k=1; k <= len; k++){
       if(bsst[k]==0) continue;
       r=ResSeq(k,xE);
       if(!MemSset(r,bsst[k])){
            fprintf(stderr,"%d.",I); PutSST(stderr,bsst[k],AB);
            fprintf(stderr,"%d != '%c'\n",k,AlphaChar(r,AB));
            assert(MemSset(r,bsst[k]));
       }
 } fprintf(stderr,"%d: %s\n",I,info->pttrn[I]);
}
#endif
	  fprintf(stderr,"  !!!!! %.2f nats; Child=%d >= %d; Parent=%d !!!!!\n",
                         bstD,CardSet(info->set[I]),MinSize,CardSet(info->set[pI]));
	}
	
	//===================================================================================
	//============================= Free memory and return ==============================
	//===================================================================================
	delete mad; if(bsst) free(bsst); Nildheap(dH); NilHist(HG);
	// info->PrintSetSizes(stderr);
#if 0
AlnSeqSW(11, 1, bstE,csqE,AB);
PutCMSA(stderr,info->sma[I]); fprintf(stderr,"pttrn = %s\n",info->pttrn[I]);
exit(1);
#endif
        if(bstE) NilSeq(bstE); if(bstCsqE) NilSeq(bstCsqE);
    	return Success;
}

e_type	*omc_typ::FindDivergentSeqs(Info *info, Int4 pI, Int4 I,sst_typ *osst,FILE *efp)
//========= Obtain divergent sequences based on pattern matches =====
{
	e_type	tmpE,*rtnE; NEW(rtnE,5,e_type);
	Int4	n,m,top_sc,top_sq,sq,k,x,r,len=mcs->RtnLengthMainCMSA();
	Int4	i,N=NumSeqsCMSA(TrueMainCMA),one4th,onehalf;
	double	sd,mean,median;
	dh_type	dH[2]; 

	dH[0]=dheap(N,3); dH[1]=dheap(N,3);

	h_type HG=Histogram("pattern matches",-10,len,1.0);
	for(top_sc=top_sq=0,sq=1; sq <= N; sq++){
          if(!MemberSet(sq,info->set[pI])) continue;
          for(n=0,m=0,k=1; k <= len; k++){
                if(osst[k]==0) continue;
                if(MemSset(ResidueCMSA(1,sq,k,TrueMainCMA),osst[k])) n++; else m++;
          } // Put sequences with more matches than mismatches into child set.
          IncdHist(n,HG);
	  insrtHeap(sq,(keytyp)n,dH[0]); insrtHeap(sq,-(keytyp)n,dH[1]);
          if(n > top_sc){ top_sc=n; top_sq=sq; }
	}  if(efp) PutHist(efp,60,HG); 
	median=MedianHist(HG); mean=MeanHist(HG); sd=sqrt(VarianceHist(HG)); NilHist(HG);
	if(efp) fprintf(efp,"mean = %.2f; + 3 sd = %.2f; + 4 std = %.2f; median = %.2f\n",
			mean,mean+3*sd,mean+4*sd,median);
	if(efp) fprintf(efp,"top_sq = %d; top_sc = %d\n",top_sq,top_sc);

	n=ItemsInHeap(dH[0]);
	one4th=(Int4) ceil((double)n/4.0);
	one4th=MINIMUM(Int4,one4th,3);
	onehalf=(Int4) ceil((double)n/2.0);
	for(i=0; i <=1; i++){		// i == 1 --> top; i == 0 = bottom.
	   ClearSet(info->set[I]);
	   for(n=0;(sq=delminHeap(dH[i])) != 0; ){
	     n++; AddSet(sq,info->set[I]);
	     if(i==1 && n >= one4th) break;
	     else if(i==0 && n >= onehalf) break;
	   } 
	    // rtnE[i]=GetSeqAsCsqCMSA(info->set[I],TrueMainCMA); 
	    rtnE[i]=SubGrpCsq(info->set[I]);
	   if(efp) fprintf(efp,"  %d seqs selected; one4th = %d\n",CardSet(info->set[I]),one4th);
           char Ala=AlphaCode('A',AB);
           for(k=1; k <= len; k++){
	     tmpE=rtnE[i];
             if(ResSeq(k,tmpE) == 0){ EqSeq(k,Ala,tmpE); EqXSeq(k,Ala,tmpE); }
           } Nildheap(dH[i]);
	} return rtnE;
}

