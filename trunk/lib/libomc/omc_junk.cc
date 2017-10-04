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

double  omc_typ::RaiseMultiBranches(FILE *fp,double Temp, Int4 iter)
// RaiseMultiBranches(stderr,temp,iter); No+=RecordChange("RaiseMultiBranches");
{
        PrintOperationHeader(fp,"Moving multiple nodes Up",iter); // mcs->PutHyperPartition(stderr);
        if((status=MultiMoveNodesUp(Temp)) > 0){
            this->RestoreBest(); mcs->PutHyperPartition(); return TrimNodes(fp,Temp,iter);
            // this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
        } return this->CalcLLR(mcs);
}

Int4	omc_typ::MultiMoveNodesUp(double Temp, set_typ pre_skip, double minLLR)
// move groups of negative child nodes up one level...
// effects background to move one at a time (especially if empty parent node between child and grandpa.
{
	Int4	pID,cID,id,i,j,n,p,r,c,worst_p,NumMoves=0,*Parent;
	double	d,worst,Worst;
	set_typ skip = MakeSet(MaxNumNodesPlus + MaxNumNodes); ClearSet(skip);
	set_typ WrstKids=0,BadKids=MakeSet(MaxNumNodesPlus + MaxNumNodes);
    do {	// Find the worst parent...
	mcs->RestoreBest(); mcs->PutHyperPartition();
	for(worst_p=0,Worst=minLLR,p=2; p <= Hpt->NumBPPS(); p++){
	   if(Hpt->TypeOfSet(p) != '?') continue;
	   pID=Hpt->ItoSetID(p); if(MemberSet(pID,skip)) continue;
	   if(stable && MemberSet(pID,stable)) continue;
	   assert(Hpt->IsScrambledTree(Parent));
	   ClearSet(BadKids);
	   // for(worst=0.0,c=p+1; c < Hpt->NumSets(); c++)	// find first parent node...in tree format.
	   for(worst=0.0,c=1; c < Hpt->NumSets(); c++)	// find first parent node...in tree format.
	   {
		if(Parent[c] != p) continue;
		if(mcs->RtnContribLLR(c,p,d)){
		   if(d < minLLR){
			cID=Hpt->ItoSetID(c); worst+=d; 
			AddSet(cID,BadKids);
			fprintf(stderr,"If node %d ('Set%d') moved up to parent of node %d ('Set%d'):",
				c,Hpt->ItoSetID(c),p,pID);
			fprintf(stderr," %.2f nats.\n",-d);
		   }
		}
	   }
	   if(worst < Worst){
		if(WrstKids) NilSet(WrstKids);
		worst_p=p; Worst=worst; WrstKids=CopySet(BadKids);
	   } free(Parent);
        } 
	if(worst_p > 0){
	   pID=Hpt->ItoSetID(worst_p);
	   fprintf(stderr,"/////// Worst = %d ('Set%d'): %.2f ///////\n",worst_p,pID,Worst);
	   PutSet(stderr,WrstKids);
mcs->PutMapContributions(stderr);
mcs->PutHyperPartition(stderr);  
	   mcs_typ *rtn_mcs = this->MultiMoveUp(WrstKids,0);
	   if(rtn_mcs){
		NumMoves++; 

		double *Map=rtn_mcs->RtnSubLPR( ),D=mcs->RtnBestLPR();
		hpt_typ *hpt=rtn_mcs->GetHpt(); i=hpt->SetIDtoI(pID);
		if(i == 0){	// pID was deleted...
		  fprintf(stderr,
		  "NewLLR=%.2f; subLLR=%.2f (%.2f min); OldLLR = %.2f; node %d ('Set%d %c').\n",
                	Map[0],Map[i],MinimumLLR,D,i,pID,hpt->TypeOfSet(worst_p));
		} else {
		  fprintf(stderr,
		  "NewLLR=%.2f; subLLR=%.2f (%.2f min); OldLLR = %.2f; node %d ('Set%d %c'); size: %d (%d min).\n",
                	Map[0],Map[i],MinimumLLR,D,i,pID,hpt->TypeOfSet(worst_p),
				rtn_mcs->SetCard(i),MinimumSetSize);
		}
		// rtn_mcs->PutHyperPartition(stderr);  // new hpt
            	fprintf(stderr,
		"!!!!!!!!!!!!!!!!!! Sampling parent node %d ('Set%d') sucessful. !!!!!!!!!!!!!!!!!!!!\n",i,pID);

		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); mcs->SaveBest=TRUE; 
		mcs->PutMapContributions(stderr);
		// id=Hpt->ItoSetID(worst_p); AddSet(id,skip); 
	   } else { id=Hpt->ItoSetID(worst_p);  AddSet(id,skip); }  // don't try this node again..
	}
    } while(worst_p > 0); NilSet(skip); return NumMoves;
}

mcs_typ	*omc_typ::MultiMoveUp(set_typ BadKids,double MinDeltaLLR)
{
	Int4	gpID,pID,P,GP,node,id,cID,t,i,j,x,parent,gp,*Parent,NumChild_P=0;
	double	lpr0,lpr1,DeltaLLR;
	set_typ	subtree=0; 
    	mcs->IsTreeHpt=TRUE;	// Need to let mcs object know hpt is a tree.

	subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
	lpr0=mcs->RtnBestLPR();
	assert(stable == 0);
	for(gpID=P=GP=0,node=2; node < Hpt->NumSets(); node++){
	  cID=Hpt->ItoSetID(node);
	  if(!MemberSet(cID,BadKids)) continue;
	  ReSetSubTree(subtree,node); // Hpt->Put(stderr,FALSE);
	  assert(Hpt->IsScrambledTree(Parent));
	  parent=Parent[node]; gp = Parent[parent]; assert(parent > 0); 
	  if(stable && MemberSet(Hpt->ItoSetID(gp),stable)) continue;
	  if(gp==0){ free(Parent); NilSet(subtree); return 0; }
	  if(P==0){
	  	for(NumChild_P=0,x=1; x <= Hpt->NumBPPS(); x++){
			if(Parent[x] == parent) NumChild_P++; 
		} P = parent;  GP=gp; gpID=Hpt->ItoSetID(gp); pID=Hpt->ItoSetID(P);
	  } else { assert(P==parent && GP==gp); }
	  NumChild_P--; free(Parent);
	  fprintf(stderr,"LLR = %.3f (gp=%d p=%d node = %d)\n",lpr0,gp,parent,node);
	  fprintf(stderr," ... moving %d up from %d to %d\n",node,parent,gp);
	  if(mcs->MoveUp(gp,parent,node,subtree)) ReSetSubTree(subtree,node);
	  // ^ this should not change the index of the nodes...
	  else print_error("MoveUp( ) operation failed");
	} // check for two internal nodes...
	assert(GP > 0 && P > 0);
	if(NumChild_P == 1){	// then move up  last child to GP and remove parent node
	    assert(Hpt->IsScrambledTree(Parent));
	    for(node=2; node < Hpt->NumSets(); node++){
		cID=Hpt->ItoSetID(node);
		if(Parent[node]==P && !MemberSet(cID,BadKids)){
	  	   ReSetSubTree(subtree,node); // Hpt->Put(stderr,FALSE);
	  	   fprintf(stderr,"(Only one child left...)\n");
	  	   fprintf(stderr,"LLR = %.3f (gp=%d p=%d node = %d)\n",lpr0,GP,P,node);
	  	   fprintf(stderr," ... moving %d up from %d to %d\n",node,P,GP);
	  	   if(mcs->MoveUp(GP,P,node,subtree)) ReSetSubTree(subtree,node);
	  	   else print_error("MoveUp( ) operation failed"); break;
		}
	    } NumChild_P--; free(Parent);
	} 
	NilSet(subtree); 
#if 0	// Sort the Hpt (a new operation without optimization) and trim vacant nodes...
#endif
	// Hpt->Put(stderr,FALSE); 
	lpr1=mcs->CalcTotalLPR( );
	// if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
	fprintf(stderr,"Initial LLR = %.2f (diff = %.3f)\n",lpr1,lpr1-lpr0);
	if(lpr1 <= 0 || (lpr1/lpr0) < 0.95){
		this->RevertToBest();  
		fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		return 0;
	} 
//!!!!!!!!!!! add a new operation to run CreateMCS() on parent subtree only!!!!!
	// mcs_typ *rtn_mcs=this->Operate('s',0,pID);
	mcs_typ *rtn_mcs=this->Operate('s',gpID,0);
	// mcs_typ *rtn_mcs=this->Operate('S',0,0);
	// lpr1=mcs->CalcTotalLPR(0,FALSE);
	if(NumChild_P == 0){ // then if node lacks minimum #nodes or LLR too low, delete it.
	    hpt_typ *hpt=rtn_mcs->GetHpt();
	    P=hpt->SetIDtoI(pID); // pID; id=Hpt->ItoSetID(P); 
	    Int4 card=rtn_mcs->SetCard(P);
	    double d=rtn_mcs->RtnSubLPR(P);
	    fprintf(stderr,"Delete node %d (Set%d)? (subLPR: %.3f; %d seqs)\n",P,pID,d,card);
	    if(d < MinimumLLR || card < MinimumSetSize){
		mcs_typ *tmp=mcs; mcs=rtn_mcs; Hpt=rtn_mcs->GetHpt( );
		mcs_typ *xmcs = this->Operate('D',0,pID,100);
		// mcs_typ *xmcs = this->DeleteNode(pID,100);
		DeleteMCS(mcs); mcs=tmp; Hpt=mcs->GetHpt(); rtn_mcs=xmcs; 
	    }
	}
	lpr1=rtn_mcs->CalcTotalLPR(0,FALSE);
	fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,lpr1-lpr0);
	DeltaLLR = lpr1 - lpr0;
	if(DeltaLLR < MinDeltaLLR || mcs->NumFailedNodes() > 0){
		this->RevertToBest(); 
		fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		DeleteMCS(rtn_mcs); return 0;
	} 
	// if(NumChild_P==1) Hpt->ChangeTypeOfSet(parent,'!'); // parent is now a leaf node.
	return rtn_mcs;
}


void    omc_typ::StoreRestoreBest(mcs_typ *xmcs)
{
        if(!xmcs->IsBestRestored( )){
             double d=xmcs->CalcTotalLPR();
             BooLean  md=xmcs->NoFailureMode;
             xmcs->NoFailureMode=FALSE;
             if(mcs->RtnBestLPR() < d){ mcs->StoreBest(); }
             mcs->RestoreBest(); mcs->NoFailureMode=md;
        }
}

e_type	omc_typ::SetToBest(set_typ set, Int4 n, mcs_typ *xmcs)
{
	hpt_typ *hpt=xmcs->GetHpt(); assert(n > 0 && n < hpt->NumSets());
        a_type  AB=AlphabetCMSA(TrueMainCMA);
        cma_typ cma=GetInSetCMSA(set,TrueMainCMA);
        FILE    *fp=tmpfile(); PutBestCMSA(fp,1,FALSE,cma); rewind(fp);
#if 1
	cma_typ cma2=ReadCMSA(fp,AB); fclose(fp);
        e_type  csqE=MkConsensusCMSA(cma2); 
	for(Int4 s=1; s <= LenSeq(csqE); s++){ 
	  if(ResSeq(s,csqE)==0) EqSeq(s,AlphaCode('A',AB),csqE);
	} xmcs->SetTheCSQ(0,n,csqE); 
#if 1
	cma_typ cmaX=MakeConsensusCMSA(cma),cmaY=MakeConsensusCMSA(cma2); 
	fprintf(stderr,"LenSeq = %d; length cmaX = %d; length cmaY = %d\n",
			LenSeq(csqE),LengthCMSA(1,cmaX),LengthCMSA(1,cmaY));
	Int4 score=PseudoAlnScoreSqToCMSA(csqE,1,cmaX);
	Int4 sc=PseudoAlnScoreSqToCMSA(csqE,1,cmaY);
	fprintf(stderr,"%d ('Set%d'): %d vs %d\n",n,hpt->ItoSetID(n),score,sc);
	TotalNilCMSA(cmaX); TotalNilCMSA(cmaY);
#else
FILE *ofp=open_file("junk",".chn","w"); 
	PutConsensusCMSA(ofp,cma2); PutCMSA(ofp,TrueMainCMA); PutConsensusCMSA(ofp,cma2); 
	// PutConsensusCMSA(ofp,cma2); PutCMSA(ofp,cma); PutConsensusCMSA(ofp,cma2); 
	exit(1);
#endif
	// PutSeq(stderr,csqE,AB);
	TotalNilCMSA(cma); TotalNilCMSA(cma2); // NilSeq(csqE); 
#else
        TotalNilCMSA(cma); cma=ReadCMSA(fp,AB); fclose(fp);
        e_type  csqE=MkConsensusCMSA(cma); TotalNilCMSA(cma);
	for(Int4 s=1; s <= LenSeq(csqE); s++){ 
	  if(ResSeq(s,csqE)==0) EqSeq(s,AlphaCode('A',AB),csqE);
	} xmcs->SetTheCSQ(0,n,csqE); // NilSeq(csqE);
#endif
	return csqE;
}

void	omc_typ::SortHpt( )
// Create a mcs_typ with Hpt sorted (apply after MoveUpAndDown())
{
     Int4    id,I,i,j,r;
     mcs_typ *rtn_mcs;  
     char    *tmp[3];

     assert(mcs->IsTreeHpt); mcs->RestoreBest();

     set_typ *set=mcs->CopyOfBestTreeSets(); 
     set_typ *in_set; NEW(in_set, MaxNumNodesPlus +2, set_typ);
     cma_typ *in_sma; NEW(in_sma, MaxNumNodesPlus +2, cma_typ);
     char **pttrn= mcs->RtnCopyOfPttrns( );
     hpt_typ *hpt=Hpt->Sort(); assert(Hpt->NumBPPS() == hpt->NumBPPS());
     for(i=1; i <= Hpt->NumBPPS(); i++){
         id=Hpt->ItoSetID(i); // assert(id <= Hpt->NumBPPS());
         sprintf(str,"Set%d",id); I=hpt->SetIDtoI(id);
	 in_sma[I]=mcs->RtnCsqAsCMSA(i,str);
	 in_set[I]=set[i];
	 // fprintf(stderr,"%d ('Set%d'): %s\n",i,id,pttrn[id]);
	 if(pttrn[i]){ tmp[0]=pttrn[i]; hpt->ReSetArgv(I,1,tmp); free(pttrn[i]); }
     } free(pttrn);

     Int4 NumSets=mcs->RtnNumElmntSetCMA( ); assert(NumSets == hpt->NumSets());
     in_set[NumSets]=set[NumSets];   // pass in set[NumSets] == RandomSet; will be recognized.
    
     rtn_mcs=CreateMCS(hpt,in_set,in_sma,100,'S',0);
     rtn_mcs->PutHyperPartition(stderr); // exit(1);
     DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
}

Int4	omc_typ::MoveUpAndDown(Int4 target,double MinDeltaLLR)
{
	Int4	s,t,i,j,x;
	double	lpr0,lpr1,DeltaLLR;
    	mcs->IsTreeHpt=TRUE;	// Need to let mcs object know hpt is a tree.
	assert(target > 0);
	mcs->SampledColumn=0;
	dh_type dH=dheap(Hpt->NumBPPS()+2,4);
	set_typ	*subtree=0; NEW(subtree,Hpt->NumBPPS()+2,set_typ);
	for(s=1; s <= Hpt->NumBPPS(); s++){ 
		subtree[s]=MakeSet(Hpt->NumBPPS()+2); // get set of sampled node subtree.
		ReSetSubTree(subtree[s],s);
		insrtHeap(s,(keytyp) Random(),dH);
	} Int4 NumOperations=0; Hpt->Put(stderr,FALSE);
	for(j=1; (s=delminHeap(dH)) != NULL; j++)
	// for(s=1; s <= Hpt->NumBPPS(); s++)
	{
	     Int4 *Parent,NumChild_P=0,NumChild_GP=0;
	     BooLean flag=FALSE;
#if 0
	     if(CardSet(subtree[target])==1 && Hpt->TypeOfSet(target) == '?'){
		Hpt->ChangeTypeOfSet(target,'!'); assert(Hpt->IsScrambledTree(Parent));
		Hpt->ChangeTypeOfSet(target,'?');
	     } else assert(Hpt->IsScrambledTree(Parent));
#endif
	     if(NumOperations > 0){ 
	        // fprintf(stderr,"***************** Sorted Hpt ***********************\n");
		// Hpt->PutSorted(stderr); // Hpt->PutAsTree(stderr); 
	        fprintf(stderr,"***************** node %d **************************\n",s);
		NumOperations=0; Hpt->Put(stderr,FALSE); 
	        fprintf(stderr,"****************************************************\n");
	     }
	     assert(Hpt->IsScrambledTree(Parent));
	     Int4 p = target, gp = Parent[target];
	     if(gp==0){ free(Parent); continue; }
	     for(NumChild_P=0,NumChild_GP=0,x=1; x <= Hpt->NumBPPS(); x++){
		if(Parent[x] == p) NumChild_P++; else if(Parent[x] == gp) NumChild_GP++;
	     }
	     BooLean NotMoved=TRUE;
	     if(Parent[s] == p && NumChild_P > 2){  // Move s up from parent to grandparent...
		// disallow internal nodes with less than 2 children!
		Int4 TargetSize=CardSet(subtree[p])-CardSet(subtree[s]);
		if(TargetSize < 1){	// don't make p childless!
		  lpr0=mcs->RtnBestLPR();
		  fprintf(stderr,"LLR = %.3f (gp=%d p=%d s = %d)\n",lpr0,gp,p,s);
		  fprintf(stderr," ... moving %d up from %d to %d\n",s,p,gp);
		  if(mcs->MoveUp(gp,p,s,subtree[s])) ReSetSubTree(subtree[s],s);
		  else print_error("MoveUp( ) operation failed");
		  NumOperations++; Hpt->Put(stderr,FALSE); 
		  // lpr1=mcs->CalcTotalLPR( );
		  lpr1=CalcLLR(mcs);
		  if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
		  fprintf(stderr,"Initial LLR = %.3f (diff = %.3f)\n",lpr1,lpr1-lpr0);
		  if(lpr1 <= 0 || (lpr1/lpr0) < 0.95){
			this->RevertToBest(); 
			fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		  } else {
		    // this->Sample(200.0,1,1,2,1); 
		    // this->Sample(200.0,2,2,3,2); 
		    // this->Sample(200.0,2,2,2,2); 
		    this->Sample(100.0,2,2,2,2); 
		    // mcs->NoFailureMode=FALSE; lpr1=mcs->CalcTotalLPR( ); mcs->NoFailureMode=TRUE;
		    // lpr1=mcs->CalcNoFailLPR( ); 
		    lpr1=CalcLLR(mcs); 
		    fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,lpr1-lpr0);
		    DeltaLLR = lpr1 - lpr0;
		    if(DeltaLLR < MinDeltaLLR || mcs->NumFailedNodes() > 0){
			this->RevertToBest();
			fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		    } else {
			if(TargetSize==2) Hpt->ChangeTypeOfSet(target,'!'); // now a leaf node.
// mcs->PutSetRelations(stderr);
			Parent[s] = gp; NumChild_GP++; NumChild_P--;
			NotMoved=FALSE;
			// mcs->PutHyperPartition(); // continue; 
		    }
	          }
		}
	     }
	     if(NotMoved && Parent[s] == gp && s != target && NumChild_GP > 2){ // Move Down...  
		// then try moving s from target's parent to target node.
		lpr0=mcs->RtnBestLPR();
		fprintf(stderr,"LLR = %.3f (gp=%d p=%d s = %d)\n",lpr0,gp,p,s);
		fprintf(stderr," ... moving %d down from %d to %d\n",s,gp,p);
		PutSet(stderr,subtree[s]);
		if(mcs->MoveDown(gp,p,s,subtree[s])) ReSetSubTree(subtree[s],s);
		else print_error("MoveDown( ) operation failed");
		NumOperations++; Hpt->Put(stderr,FALSE); 
		// mcs->NoFailureMode=FALSE; lpr1=mcs->CalcTotalLPR( ); mcs->NoFailureMode=TRUE;
		// lpr1=mcs->CalcTotalLPR( ); 
		lpr1=CalcLLR(mcs);
		if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
		fprintf(stderr,"Initial LLR = %.3f (diff = %.3f; %.1f%c)\n",
						lpr1,lpr1-lpr0,100.0*(lpr1/lpr0),'%');
		if(lpr1 <= 0 || (lpr1/lpr0) < 0.95){ 
			this->RevertToBest(); 
			fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		} else {
		   // this->Sample(200.0,1,1,2,1); 
		   // this->Sample(200.0,2,2,2,2); 
		   this->Sample(100.0,2,2,2,2); 
		   // lpr1=mcs->CalcTotalLPR( ); 
		   lpr1=mcs->CalcNoFailLPR( ); 
		   fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,lpr1-lpr0);
		   DeltaLLR = lpr1 - lpr0;
		   if(DeltaLLR < MinDeltaLLR || mcs->NumFailedNodes() > 0){ // Failed to improve....
			this->RevertToBest();
			fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		   } else {
fprintf(stderr,"node %d MovedUp&Down() operation improvement!\n",target);
			// mcs->PutHyperPartition();
		   }
// mcs->PutSetRelations(stderr);
		}
	     } free(Parent);
	} Nildheap(dH); 
	Int4 rtn=CardSet(subtree[target]);
	for(s=1; s <= Hpt->NumBPPS(); s++){ if(subtree[s]) NilSet(subtree[s]); } free(subtree);
	this->SortHpt( );
	return rtn;
}

hsi_typ	*omc_typ::ArrangeAsTree0(hsi_typ **hsi,Int4 num,set_typ *Set)
{
	double tLpr,oLpr,min_oLpr,tWtCnt,oWtCnt,ratio,d;
	FILE *fp=stderr;
	Int4	s,t;
	hsi_typ *root=0,*Y,*Z; 
	for(s=1; s <= num; s++){ if(hsi[s] == 0) continue; else { root=hsi[s]; break; } }
	for(t=s+1; t < num; t++){
	   if(hsi[t] == 0) continue;
	   else Y=hsi[t]; // Y->Put(stderr); Y->Put(stderr,Set,num);
	   Z=root->PlaceIntoTree(Y,Set,num,20); 
	   if(Z == 0){ hsi[t]=0; delete Y; } 
	   else if(root != Z){ root = Z; }
	} if(root) root->PutTree(stderr,Set,Hpt->NumBPPS());
	return root;
}

Int4    omc_typ::ResurrectRejects0(double Temp)
// If there are MinimumSplitSize  or more nodes rejected, try to resurrect them by creating a new node.
/***************************************************************
  (x)X <---  Random ---> X(x)                   
     |                   |           |
     |                   O           O(new root; empty)
     |    insert node    | 
     O                   * x (resurrected)
   / | \                /|\
  o  o  o              o o o
 ***************************************************************/
{
        if(mcs->SetCard(Hpt->NumSets()) < MinimumSplitSize) return 0;

        Int4    x,p,id,i,j,n,R,r,c,I,pI,pID,ID,xID,*Parent;

        // 1. Insert a node between the root and ALL of its current children.
        for(ID=1; ID <= MaxNumNodes; ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID:
        pID = Hpt->ItoSetID(1);
#if 1
	if(best_mcs != 0) RestoreFinalBest( );
	mcs_typ *xmcs=RtnCopy(); 
	mcs->MoveSeqs(Hpt->NumSets(),1);
	mcs->CalcTotalLPR(); mcs->SaveBest=TRUE; mcs->StoreBest();
        mcs->PutHyperPartition(stderr);
	// SetStringency(0);
	AddLeaves(300,'A');
        mcs->PutHyperPartition(stderr);
	SetStringency(stringency);
	TrimNodes(stderr,300,1);
        mcs->PutHyperPartition(stderr);
#else
	ClearSet(ChildSetBG); ClearSet(ChildSetFG); FillSet(ChildSetFG);
        hsi_typ *found= new hsi_typ(ChildSetFG,ChildSetBG,0,0,0,0,lpr);
        mcs_typ *rtn_mcs=InsertNode(pID,ID,found);
	hpt_typ *hpt=rtn_mcs->GetHpt();
        rtn_mcs->PutHyperPartition(stderr);

	// 2. Move the root node sequences to the inserted node.
	rtn_mcs->MoveSeqs(1,hpt->SetIDtoI(ID));
        rtn_mcs->PutHyperPartition(stderr);

	// 3. Move the root node sequences to the inserted node.
	rtn_mcs->MoveSeqs(hpt->NumSets(),1);
        rtn_mcs->PutHyperPartition(stderr);

        // 4. Add a leaf to the root node; sample over this...
	mcs_typ *xmcs=mcs; mcs=rtn_mcs; Hpt=mcs->GetHpt();
        rtn_mcs->PutHyperPartition(stderr);
        for(xID=1; xID <= MaxNumNodes; xID++){ if(Hpt->SetIDtoI(xID) == 0) break; } // find an unused ID:
	rtn_mcs=this->Operate('A',1,xID); hpt=rtn_mcs->GetHpt();
	rtn_mcs->MoveSeqs(1,hpt->SetIDtoI(xID));
        rtn_mcs->PutHyperPartition(stderr);
	rtn_mcs->Sample(2,2,2,2);
        rtn_mcs->PutHyperPartition(stderr);
#endif

}

hsi_typ	*omc_typ::Check4Bottom(Int4 depth, Int4 maxdepth)
{
        double	Lpr,d,T_neg=0.0,T_pos=0.0;
	Int4	i,s,c,cID,kID,N_neg=0,N_pos=0,NumCols=mcs->RtnDefaultMaxCol();
	sst_typ	*xsst=0;
	hsi_typ	*hsi=0;
	char	typ_node;
	
	if(depth >= maxdepth || (CardSet(ChildSetBoth) == CardSet(ChildSetU) 
		&& CardInterSet(ChildSetBoth,ChildSetU) == CardSet(ChildSetU))){
           // xsst=this->GetOptPttrnLPR(stderr,SetFG,SetBG,Lpr,NumCols,'m'); 
	   unsigned char *csq=0; 
           xsst=this->GetOptPttrnLPR(stderr,SetFG,SetBG,Lpr,NumCols,csq,'m'); 
	   hsi= new hsi_typ(ChildSetFG,ChildSetBG,xsst,Lpr,CardSet(SetFG),CardSet(SetBG),lpr);
	   if(depth != maxdepth){	// showing output only!!
		e_type pE= MkSeq("consensus", LengthCMSA(1,TrueMainCMA),csq); 
		PutSeq(stderr,pE,AB);
	        IncdHist(Lpr,dfsHG);
		N_neg=0; N_pos=0; T_neg=0.0; T_pos=0.0;
		fprintf(stderr,"+.Set+: (*) %.2f (%d:%d)\n",Lpr,CardSet(SetFG),CardSet(SetBG));
		set_typ ChildFG=MakeSet(SetN(ChildSetFG)),ChildBG=CopySet(ChildSetBG);
		ClearSet(ChildFG); 
	        for(c=1; c <= Hpt->NumBPPS(); c++){ // Calculate sub-LLRs for each set.
		  kID=Hpt->ItoSetID(c);
		  if(MemberSet(kID,ChildSetFG)){
			e_type cE=mcs->RtnQryBPPS(c);
			double D,DD,dd,Dn,Dp,dp,dn;
			Int4   B,C,Bp,Bn,Cp,Cn;
			d=CalcSetvsPttrnLPR(0,Set[c],SetBG,xsst,'m',pE);
		        C=CardSet(Set[c]); B=CardSet(SetBG);

			fprintf(stderr,"%d.Set%d(%d): (+) ",c,kID,C);
			if(d < 0.0) fprintf(stderr,"*%.2f*",d); else fprintf(stderr,"%.2f",d); 
			if(d < 0.0){ N_neg++; T_neg+=d; } else { N_pos++; T_pos+=d; }
			fprintf(stderr," (%d:%d) ",C,B);

			UnionSet3(SetBG,TmpSet,TmpBG);	// Add Parent Set (TmpSet) to BG.
			Dp=CalcSetvsPttrnLPR(0,Set[c],TmpBG,xsst,'m',pE);
           		free(this->GetOptPttrnLPR(stderr,SetFG,TmpBG,dp,NumCols,'m',pE)); 
			Bp=CardSet(TmpBG); Cp=CardSet(SetFG);

			if(Dp < 0.0){ fprintf(stderr," (++) *%.2f* ",Dp); }
			else { fprintf(stderr," (++) %.2f ",Dp); }
			fprintf(stderr,"(%d:%d) ",C,Bp);

			UnionSet3(TmpBG,Set[c],TmpBG);	// Put Set[c] in BG...
			IntersectNotSet(SetFG,Set[c],TmpFG);  // & remove from FG.
			Bn=CardSet(TmpBG); Cn=CardSet(TmpFG);
			Dn=CalcSetvsPttrnLPR(0,TmpFG,TmpBG,xsst,'m',pE);
           		free(this->GetOptPttrnLPR(stderr,TmpFG,TmpBG,dn,NumCols,'m',pE)); 

			sst_typ *csst=mcs->RtnCopyOfSST(c);
			if(Hpt->TypeOfSet(c) == '?') typ_node='m'; else typ_node='l';
			D=CalcSetvsPttrnLPR(0,Set[c],TmpFG,csst,typ_node,cE); free(csst);

#if 0
			DD=Dp-Dn; 
			if(DD < 0.0) fprintf(stderr,"(+-) *%.3f* ",DD); 
			else fprintf(stderr,"(-) %.3f ",DD); 
			fprintf(stderr," (%d:%d) (*) %.2f\n",Cn,Bn,dd);
#endif

			dd=dp-dn;
			if(dd < 0.0) fprintf(stderr,"(+-) *%.3f*",dd);
			else fprintf(stderr,"(+-) %.3f",dd); 
			fprintf(stderr," (%d:%d) (*) %.2f\n",Cp,Bp,D);
			if(dd > 0.0) AddSet(kID,ChildFG); else AddSet(kID,ChildBG);

			NilSeq(cE);
		  }
		} PrintPttrn(stderr,xsst); NilSeq(pE);
		// PutSet(stderr,ChildSetFG); PutSet(stderr,ChildSetBG); 
	        fprintf(stderr,"       --- LPR=%.3f (%.2f)---\n",Lpr,T_neg);
		if(CardSet(ChildSetBG) > 1){
		  N_neg=0; N_pos=0; T_neg=0.0; T_pos=0.0;
	          for(c=1; c <= Hpt->NumBPPS(); c++){ // Calculate sub-LLRs for each set.
		    kID=Hpt->ItoSetID(c);
		    if(MemberSet(kID,ChildSetBG)){
			IntersectNotSet(SetBG,Set[c],TmpBG); // TmpBG = SetBG intersect not Set[c] 
			d=CalcSetvsPttrnLPR(0,Set[c],TmpBG,xsst,'m');
			if(d < 0.0) continue;
			fprintf(stderr,"--- background node %d ('Set%d') (%d seqs): ",
				c,kID,CardSet(Set[c]));
			if(d < 0.0){ fprintf(stderr,"subLPR=*%.3f* ---\n",d); N_neg++; }
			else { fprintf(stderr,"subLPR=%.3f ---\n",d); N_pos++; }
			if(d < 0.0){ N_neg++; T_neg+=d; } else { N_pos++; T_pos+=d; }
		    }
		  }
		} fprintf(stderr,"\n");
		// Optimize final hsi based on above analysis...
#if 1
hsi=OptimizeInsertion(ChildFG,ChildBG,hsi);
NilSet(ChildFG); NilSet(ChildBG);
#else
		if(CardInterSet(ChildFG,ChildSetFG) != CardSet(ChildSetFG)){
		   ClearSet(TmpFG); ClearSet(TmpBG);
	           for(c=1; c <= Hpt->NumBPPS(); c++){ // Calculate sub-LLRs for each set.
		  	kID=Hpt->ItoSetID(c);
			if(MemberSet(kID,ChildFG)){ UnionSet(TmpFG,Set[c]); }
			if(MemberSet(kID,ChildBG)){ UnionSet(TmpBG,Set[c]); }
		   } free(csq); csq=0; 
       		   sst_typ *tsst=this->GetOptPttrnLPR(stderr,TmpFG,TmpBG,Lpr,NumCols,csq,'m'); 
		   // pE= MkSeq("consensus", LengthCMSA(1,TrueMainCMA),csq); 
		   hsi->Store(ChildFG,ChildBG,tsst,Lpr,CardSet(TmpFG),CardSet(TmpBG));
		} NilSet(ChildFG); NilSet(ChildBG);
#endif
	   }
	   free(csq); xsst=0; 

	} return hsi;
}

hsi_typ	*omc_typ::OptimizeInsertion(set_typ ChildFG, set_typ ChildBG, hsi_typ *hsi)
{
        double	Lpr,d;
	Int4	id,i,s,r,c,ID,kID,pID,parent=0,*P=0,NumCols=mcs->RtnDefaultMaxCol();
	sst_typ	*xsst=0;
	char	*pttrn;
	unsigned char *csq=0; 
	mcs_typ	*rtn_mcs,*xmcs=0;
	BooLean	changed;
	
     do {
	xmcs=this->RtnCopy(); 
	parent=0; changed=FALSE;
	assert(Hpt->IsTree(P));
	// Optimize final hsi based on above analysis...
	if(CardInterSet(ChildFG,ChildSetFG) != CardSet(ChildSetFG)){
		ClearSet(TmpFG); ClearSet(TmpBG);
	        for(c=1; c <= Hpt->NumBPPS(); c++){ // Calculate sub-LLRs for each set.
		  	kID=Hpt->ItoSetID(c);
			if(MemberSet(kID,ChildFG)){ if(parent==0) parent=P[c]; UnionSet(TmpFG,Set[c]); }
			if(MemberSet(kID,ChildBG)){ UnionSet(TmpBG,Set[c]); }
		} 
       		sst_typ *tsst=this->GetOptPttrnLPR(stderr,TmpFG,TmpBG,Lpr,NumCols,csq,'m'); 
		fprintf(stderr,"       --- LPR=%.3f ---\n",Lpr);
		// pE= MkSeq("consensus", LengthCMSA(1,TrueMainCMA),csq); 
		hsi->Store(ChildFG,ChildBG,tsst,Lpr,CardSet(TmpFG),CardSet(TmpBG));
		CopySet(ChildSetFG,ChildFG); CopySet(ChildSetBG,ChildBG); 
		free(csq);
	} free(P);
	if(parent==0){
		PutSet(stderr,ChildFG); PutSet(stderr,ChildBG);
     		if(CardSet(ChildFG) < 2){ delete hsi; return 0; }
     		if(CardSet(ChildBG) < 1){ delete hsi; return 0; }
	}
	assert(parent); pID=Hpt->ItoSetID(parent); assert(pID > 0 && pID < MaxNumNodes); 
	for(ID=1; ID <= MaxNumNodes; ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID
	assert(ID > 0 && ID <= MaxNumNodes);
	rtn_mcs=InsertNode(pID,ID,hsi); // pID is the Set identifier, not the index!
        hpt_typ *hpt=rtn_mcs->GetHpt();
	parent=hpt->SetIDtoI(pID);
	Int4 node=hpt->SetIDtoI(ID);
	d=rtn_mcs->SampleColumns( );
	// for(c=1; c <= hpt->NumBPPS(); c++) // Calculate sub-LLRs for each set.
	c=parent; c=node;
	PutSet(stderr,ChildFG); PutSet(stderr,ChildBG);
	char **HP=hpt->RtnHyperPartition();
	for(r=1; r < hpt->NumSets(); r++){ // Calculate sub-LLRs for each set.
	   if(r==parent) continue;
	   if(rtn_mcs->RtnContribLLR(r,c,d)){
	   	id=hpt->ItoSetID(r);
		if(HP[r][c] == '+'){
		    fprintf(stderr,"%d.Set%d: (+) %.2f\n",r,id,d);
		    if(d < 0.0){ DeleteSet(id,ChildFG); AddSet(id,ChildBG); changed=TRUE; }
		}
		if(HP[r][c] == '-'){
		    fprintf(stderr,"%d.Set%d: (-) %.2f\n",r,id,d);
		    if(d < 0.0){ DeleteSet(id,ChildBG); AddSet(id,ChildFG); changed=TRUE; }
		} assert(HP[r][c] != 'o');
	   }
	}
        rtn_mcs->PutMapContributions(stderr,lpr);
        // PutSeq(stderr,rtn_mcs->RtnKeySeq(hpt->SetIDtoI(ID)),AB);
        rtn_mcs->PutHyperPartition(stderr);
        // pttrn=rtn_mcs->RtnCopyOfPttrn(hpt->SetIDtoI(ID));
        // fprintf(stderr,"%s\n",pttrn); free(pttrn);
	if(changed){ 
	   DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
	} else {
	   DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
	   // DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
	}
     } while(changed);
     if(CardSet(ChildFG) < 2){ delete hsi; hsi=0; }
     else if(CardSet(ChildBG) < 1){ delete hsi; hsi=0; }
     return hsi;
}

