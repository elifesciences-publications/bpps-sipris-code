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

Int4    omc_typ::FuseNodes(double minLLR, double Temp)
// Merge two sibling nodes into a single node.
{
        Int4    Sibs,id,i,n,worst_i,worst_s,NumFused=0,*Parent,OldNumSets,NumCols=mcs->RtnDefaultMaxCol();
        double  dd,d,D,DD,worst,*subLLR;
	char	NdTyp;
	sst_typ	*xsst;
        set_typ skip = MakeSet(MaxNumNodesPlus + MaxNumNodes); ClearSet(skip);
	mcs->SaveBest=TRUE;
    do {
        mcs->RestoreBest(); 
	// mcs->PutHyperPartition(); 
	Hpt->IsScrambledTree(Parent);
	OldNumSets=Hpt->NumSets();
	assert(Set==0); NEW(Set,OldNumSets + 3,set_typ);
        for(i=1; i < OldNumSets; i++){
            for(Sibs=0,n=2; n < Hpt->NumSets(); n++){     // find sibling nodes...in tree format.
		if(Parent[i] == Parent[n]) Sibs++;	// count number of sibling nodes.
	    } if(Sibs == 0) AddSet(Hpt->ItoSetID(i),skip);
	    Set[i]=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(i,Set[i]); 
	} NEW(subLLR,OldNumSets + 3, double);
        for(i=2; i < OldNumSets; i++){
		if(MemberSet(Hpt->ItoSetID(i),skip)) continue;
	        if(Hpt->TypeOfSet(i) == '!') NdTyp='L'; else NdTyp='M';
		IntersectNotSet(Set[Parent[i]],Set[i],TmpBG); // tmpBG = SetParent & ! Set[i].
		if(CardSet(TmpBG) == 0){ AddSet(Hpt->ItoSetID(i),skip); }
		else {
	          xsst=this->GetOptPttrnLPR(0,Set[i],TmpBG,d,NumCols,NdTyp); free(xsst);  // Current LLR.
		  subLLR[i]=d;
		}
	}
	assert(preHG == 0); preHG=Histogram("preliminary change in LLR",0,2000,100);
        for(worst_i=0,worst=DBL_MAX,i=2; i < Hpt->NumSets(); i++){
           id=Hpt->ItoSetID(i); if(MemberSet(id,skip)) continue;
	   if(stable && MemberSet(id,stable)) continue;
	   if(stable && MemberSet(Hpt->ItoSetID(Parent[i]),stable)) continue;
           for(n=i+1; n < Hpt->NumSets(); n++){     // find sibling nodes...in tree format.
		if(Parent[i] != Parent[n]) continue;	// Sibling nodes only.

		UnionSet3(Set[i],Set[n],TmpFG); // TmpFG = Set[i] U Set[n].
	   	IntersectNotSet(Set[Parent[i]],TmpFG,TmpBG); // tmpBG = SetParent & !TmpFG.
	   	if(Hpt->TypeOfSet(i) == '!' && Hpt->TypeOfSet(n) == '!') NdTyp='L'; else NdTyp='M';
		xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,d,NumCols,NdTyp); free(xsst); // Proposed LLR.
		D = subLLR[i] + subLLR[n];
		if(1 || d > D){	// Proposed is better than current.
		  fprintf(stderr,"%d + %d: (%d + %d) = %.2f nats",
				i,n,CardSet(Set[i]),CardSet(Set[n]),D);
		  fprintf(stderr,"; fused --> (%d vs %d) = %.2f nats; diff = %.2f.\n",
				CardSet(TmpFG),CardSet(TmpBG),d,d-D);
		} DD = d - D;
		if((DD < 0.0)) continue;
                else { IncdHist(DD,preHG); if(-DD < worst){ worst=-DD; worst_i=i; worst_s=n; } }
           }
        }  free(Parent); PutHist(stderr,60,preHG); NilHist(preHG); preHG=0;
	for(i=1; i < OldNumSets; i++) NilSet(Set[i]); free(Set); Set=0; free(subLLR);
	if(worst_i > 0){
           id=Hpt->ItoSetID(worst_i);
           fprintf(stderr,"/////// Worst = %d ('Set%d'): %.2f ///////\n",worst_i,id,worst);
           mcs_typ *rtn_mcs = this->Fuse(worst_i,worst_s,minLLR);
           if(rtn_mcs){
	        rtn_mcs->SampleColumns();
		rtn_mcs->PutHyperPartition(stderr); // exit (1);
                NumFused++; DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); mcs->SaveBest=TRUE;
           } else { id=Hpt->ItoSetID(worst_i);  AddSet(id,skip); }  // don't try this node again..
        }
    } while(worst_i > 0); NilSet(skip); 
    return NumFused;
}

mcs_typ	*omc_typ::Fuse(Int4 node, Int4 sibling, double MinDeltaLLR)
// Fuse the two nodes into one node...
{
	fprintf(stderr,"Attempting to fuse nodes %d and %d\n",node,sibling);
        Int4    id,t,i,j,x,n,parent,gp,*Parent,NumChild_P=0;
        double  lpr0,lpr1,DeltaLLR;
	mcs_typ *rtn_mcs=0,*tmp_mcs=0,*old_mcs;
	// this->RevertToBest();
        set_typ subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
        mcs->IsTreeHpt=TRUE;    // Need to let mcs object know hpt is a tree.
        assert(node > 1 && node < Hpt->NumSets()); 
	assert(sibling > 1 && sibling < Hpt->NumSets()); 
	assert(Hpt->IsScrambledTree(Parent));

        gp=Parent[node]; 
        if(gp==0 || Parent[sibling] != gp){ free(Parent); return 0; }
	if(stable && MemberSet(Hpt->ItoSetID(gp),stable)){ free(Parent); return 0; }
        free(Parent);
	mcs->SaveBest=TRUE; lpr0=mcs->CalcTotalLPR( );
        // lpr0=mcs->RtnBestLPR();
	lpr0=CalcLLR(mcs);
        fprintf(stderr,"LLR = %.3f (gp=%d p=%d node = %d)\n",lpr0,gp,sibling,node);
        fprintf(stderr," ... moving node %d down from node %d to node %d\n",node,gp,sibling);

	tmp_mcs=this->Operate('C',0,0);   // store the current best state.
	old_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	id=Hpt->ItoSetID(node);
        ReSetSubTree(subtree,node); // Hpt->Put(stderr,FALSE);
        if(!mcs->MoveDown(gp,sibling,node,subtree)) print_error("MoveDown( ) operation failed");
        NilSet(subtree); 
        tmp_mcs=this->Operate('u',0,0);   // unscramble mcs without optimizing...
	DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt( ); // this->RevertToBest();  
	mcs->SampleColumns();
	mcs->PutHyperPartition(stderr); 
	rtn_mcs=this->Operate('d',0,id);
	rtn_mcs->SaveBest=TRUE;
	rtn_mcs->SampleColumns(); rtn_mcs->PutHyperPartition(stderr); 
        lpr1=rtn_mcs->CalcTotalLPR( );
        // if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
        fprintf(stderr,"Initial LLR = %.2f; old LLR = %.2f (diff = %.3f)\n",lpr1,lpr0,lpr1-lpr0);
 	DeleteMCS(mcs); mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
        if(lpr1 <= 0 || lpr0 > 0 && (lpr1/lpr0) < 0.95){
                fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
                return 0;
        }
	n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
        lpr1=rtn_mcs->CalcTotalLPR(0,FALSE); DeltaLLR = lpr1 - lpr0;
        fprintf(stderr,"Presampling LLR = %.3f (diff = %.3f)\n",lpr1,DeltaLLR);
        // if(DeltaLLR > MinDeltaLLR && n == 0){ return rtn_mcs; }
        if(DeltaLLR > MinDeltaLLR){ return rtn_mcs; }
        // mcs_typ *rtn_mcs=this->Operate('S',0,0);
	Sample(300.0,2,2,2,2,rtn_mcs,0.03);
        lpr1=rtn_mcs->CalcTotalLPR();
        DeltaLLR = lpr1 - lpr0;
        fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,DeltaLLR);
	n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
        if(DeltaLLR > MinDeltaLLR && n == 0){ return rtn_mcs; }
        fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
        DeleteMCS(rtn_mcs); return 0;
}

