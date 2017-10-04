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

double	omc_typ::OnTheFlyMvUp(Int4 Mv, Int4 *Parent)
// SetFG & SetBG are initialized from calling environment and correspond to node Mv's current FG & BG.
// Don't modify SetFG & SetBG!!!
// TmpFG & TmpBG are used for computing New FG/BG.
{
        Int4    Kids,p,gp,i,j,c,sibling,NumCols=mcs->RtnDefaultMaxCol();
	char	NdTyp;
        double  DD,D,dd,d,deltaLLR=0.0;
	sst_typ	*xsst;

	assert(Mv > 1 && Mv < Hpt->NumSets());
	p = Parent[Mv]; assert(p > 1);
	gp = Parent[p]; assert(gp > 0 && gp < Hpt->NumSets());

	// Calculate node p node's current contrast alignment LLR.
	mcs->CopySubTreeSeqSet(p,TmpFG); mcs->CopyBkGrndSeqSet(p,TmpBG);
	xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,DD,NumCols,'M'); free(xsst); // Current LLR.
	sprintf(str,"parent %d --> %d: (%d vs %d; %.2f nats)",p,gp,CardSet(TmpFG),CardSet(TmpBG),DD);

	// Calculate node p's proposed contrast alignment LLR.
	if(Hpt->TypeOfSet(Mv) == '!') NdTyp='L'; else NdTyp='M';
	IntersectNotSet(TmpFG,SetFG); // TmpFG=TmpFG & !SetFG = p's proposed BG 
	UnionSet(TmpBG,SetFG); // TmpBG = TmpBG U SetFG.
	for(Kids=0,i=2; i < Hpt->NumSets(); i++){     // find Mv's sibling nodes...in tree format.
                if(Parent[i] == p) Kids++;      // count number of sibling nodes.
        } if(Kids <= 1) NdTyp='L'; else NdTyp='M';
	xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,dd,NumCols,NdTyp); free(xsst); // Proposed LLR.
	D = dd-DD;
	if(D > 0.0) fprintf(stderr,"%s vs %d --> %d (%d vs %d; %.2f nats) change = %.2f.\n",
				str,p,gp,CardSet(TmpFG),CardSet(TmpBG),dd,D);
	deltaLLR += D;
	// compute influence on 'gp's other children (if any) as well...
	for(c=2; c < Hpt->NumSets(); c++){     
	   if(c == Mv || Parent[c] != p) continue;  // Examine other children of proposed parent. 
	   if(Hpt->TypeOfSet(c) == '!') NdTyp='L'; else NdTyp='M';
	   // c's current LLR.
	   mcs->CopySubTreeSeqSet(c,TmpFG); mcs->CopyBkGrndSeqSet(c,TmpBG);
	   xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,D,NumCols,NdTyp); free(xsst); // Current LLR.
	   sprintf(str," sibling %d --> %d: (%d vs %d; %.2f nats)",c,p,CardSet(TmpFG),CardSet(TmpBG),D);
	   // c's proposed LLR.
	   IntersectNotSet(TmpBG,SetFG); // TmpBG=TmpBG & !SetFG = c's proposed BG 
	   xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,d,NumCols,NdTyp); free(xsst); // Current LLR.
	   if(d-D > 0.0) fprintf(stderr,"%s vs %d --> %d - %d (%d vs %d; %.2f nats) change = %.2f.\n",
				str,c,p,Mv,CardSet(TmpFG),CardSet(TmpBG),d,d-D);
	   deltaLLR += d-D;
	} return deltaLLR;
}

Int4	omc_typ::MoveNodesUp(double Temp, double minLLR)
{
	Int4	id,i,j,p,gp,worst_i,worst_n,NumMoves=0,*Parent,NumCols=mcs->RtnDefaultMaxCol();
	double	D,d,dd,DD,worst;
	char	NdTyp;
	sst_typ	*xsst;
	set_typ skip = MakeSet(MaxNumNodesPlus + MaxNumNodes); ClearSet(skip);
    do {
	mcs->RestoreBest(); // mcs->PutHyperPartition(); 
	Hpt->IsScrambledTree(Parent);
	for(worst_i=0,worst=DBL_MAX,i=2; i < Hpt->NumSets(); i++){
	   id=Hpt->ItoSetID(i); if(MemberSet(id,skip)) continue;
	   if(stable && MemberSet(Hpt->ItoSetID(i),stable)) continue;
	   if(stable && MemberSet(Hpt->ItoSetID(Parent[i]),stable)) continue;
	   p=Parent[i];
	   if(p == 1) continue;  // node i attached at the root.
	   else gp=Parent[p];	// n == grandparent node.

	   // 1. compute node i's current contrast alignment LLR.
           if(Hpt->TypeOfSet(i) == '!') NdTyp='L'; else NdTyp='M';
	   mcs->CopySubTreeSeqSet(i,SetFG); mcs->CopyBkGrndSeqSet(i,SetBG);
           xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,D,NumCols,NdTyp); free(xsst);  // Current LLR.
           fprintf(stderr,"current %d --> %d: (%d vs %d) = %.2f nats.\n",i,p,CardSet(SetFG),CardSet(SetBG),D);

	   // 1. compute node i's proposed contrast alignment LLR.
           mcs->CopySubTreeSeqSet(gp,TmpBG);        // TmpBG == gp's current FG
	   IntersectNotSet(TmpBG,SetFG);	// TmpBG=TmpBG & !SetFG = i's proposed BG
           xsst=this->GetOptPttrnLPR(0,SetFG,TmpBG,d,NumCols,NdTyp); free(xsst); // Proposed LLR.
           fprintf(stderr,"  proposed --> %d (vs %d) = %.2f nats.\n",gp,CardSet(TmpBG),d);
           // if(D > d) continue; else D = d-D;
	   dd = this->OnTheFlyMvUp(i,Parent);
	   DD = (d - D) + dd;
	   if((DD < 0.0)) continue;
	   fprintf(stderr," ---> total change = (%.2f - %.2f) + %.2f = %.2f.\n",d,D,dd,DD);
#if 0	// if node is a leaf and an only child compute fusing it into parent node.
           if(Hpt->TypeOfSet(i) == '!' && Hpt->NumDescendants(p) == 1) {
             mcs->CopySubTreeSeqSet(p,TmpFG);        // TmpBG == gp's current FG
	     mcs->CopyBkGrndSeqSet(p,TmpBG);
             xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,d,NumCols,'L'); free(xsst); // Proposed LLR.
             fprintf(stderr," Fuse %d + %d: (%d vs %d) = %.2f nats.\n",i,p,CardSet(TmpFG),CardSet(TmpBG),d);
	   }	// It seems that this can't improve on this.
#endif
	   if(-DD < worst){ worst=-DD; worst_i=i; worst_n=gp; } 
        } free(Parent);
	if(worst_i > 0){
	   id=Hpt->ItoSetID(worst_i);
	   fprintf(stderr,"/////// Up: Worst = %d ('Set%d'): %.2f increase ///////\n",worst_i,id,-worst);
	   mcs_typ *rtn_mcs = this->MoveUp(worst_i,0.0);
	   if(rtn_mcs){
		NumMoves++; 
		double *Map=rtn_mcs->RtnSubLPR( ),D=mcs->RtnBestLPR();
		hpt_typ *hpt=rtn_mcs->GetHpt(); i=hpt->SetIDtoI(id);
		fprintf(stderr,
		"NewLLR=%.2f; subLLR=%.2f (%.2f min); OldLLR=%.2f; node %d ('Set%d %c'); %d seqs (%d min) %.1f npws.\n",
                	Map[0],Map[i],MinimumLLR,D,i,id,hpt->TypeOfSet(worst_i),
				rtn_mcs->SetCard(i),MinimumSetSize,rtn_mcs->RtnNatsPerWtSeq(i));
		// rtn_mcs->PutHyperPartition(stderr);  // new hpt
            	fprintf(stderr,
		"!!!!!!!!!!!!!!!!!! Sampling node %d ('Set%d') sucessful. !!!!!!!!!!!!!!!!!!!!\n",i,id);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); mcs->SaveBest=TRUE; 
	   } else { id=Hpt->ItoSetID(worst_i);  AddSet(id,skip); }  // don't try this node again..
	}
    } while(worst_i > 0); NilSet(skip); return NumMoves;
}

mcs_typ	*omc_typ::MoveUp(Int4 node, double MinDeltaLLR)
{
	Int4	id,t,i,j,n,x,parent,gp,*Parent;
	double	d,lpr0,lpr1,DeltaLLR;
	mcs_typ *rtn_mcs=0,*old_mcs=0,*tmp_mcs=0;

    	mcs->IsTreeHpt=TRUE;	// Need to let mcs object know hpt is a tree.
	assert(node > 1 && node < Hpt->NumSets());
	assert(Hpt->IsScrambledTree(Parent));
	parent=Parent[node]; gp = Parent[parent]; assert(parent > 0); 
	if(gp==0){ free(Parent); return 0; }
	if(stable && MemberSet(Hpt->ItoSetID(gp),stable)){ free(Parent); return 0; }
	free(Parent); 
	lpr0=mcs->RtnBestLPR();
	fprintf(stderr,"LLR = %.3f (gp=%d p=%d node = %d)\n",lpr0,gp,parent,node);
	fprintf(stderr," ... moving %d up from %d to %d\n",node,parent,gp);

	tmp_mcs=this->RtnCopy(); // MoveUp will change mcs!!! Need to save a copy!
        old_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );

	set_typ	subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
	ReSetSubTree(subtree,node); // Hpt->Put(stderr,FALSE);
	if(!mcs->MoveUp(gp,parent,node,subtree)) print_error("MoveUp( ) operation failed");
	// Caution: the BestMCS is not changed; make sure Best is not restored prior to unscrambling!!!
	NilSet(subtree); 
	rtn_mcs=this->Operate('u',0,0);	// unscramble hpt...
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
	mcs->SampleColumns();
// mcs->PutHyperPartition(stderr);
	// Hpt->Put(stderr,FALSE); 
	lpr1=mcs->CalcTotalLPR( );
	// if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
	d = lpr1-lpr0;
	fprintf(stderr,"Initial LLR = %.2f --> %.2f (diff = %.3f)\n",lpr0,lpr1,d);
	// if(lpr1 <= 0 || (lpr0 > 0.0 && (lpr1/lpr0) < 0.95))
	if(lpr1 <= 0 || d < 0.0)
	{
		// this->RevertToBest();  
		fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
		DeleteMCS(mcs); mcs=old_mcs; Hpt=mcs->GetHpt(); return 0;
	} 
// mcs->PutHyperPartition(stderr);
	rtn_mcs=mcs; mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
	rtn_mcs->SaveBest=TRUE;
        n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
        lpr1=rtn_mcs->CalcTotalLPR(0,FALSE);
        DeltaLLR = lpr1 - lpr0;
        fprintf(stderr,"Presampling LLR = %.3f (diff = %.3f; %d failed nodes)\n",
		lpr1,DeltaLLR,n);
        if(DeltaLLR > MinDeltaLLR && n == 0){ return rtn_mcs; }
        // if(DeltaLLR > MinDeltaLLR){ return rtn_mcs; }

	// Sample(300.0,2,2,2,2,rtn_mcs,0.03);
	rtn_mcs->Sample(2,2,2,2);	// do once only!
	lpr1=rtn_mcs->CalcTotalLPR(0,FALSE);
	fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,lpr1-lpr0);
	DeltaLLR = lpr1 - lpr0;
	n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
	if(DeltaLLR >= MinDeltaLLR && n == 0) return rtn_mcs;
	// this->RevertToBest();
	fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
	DeleteMCS(rtn_mcs);  return 0;
}

mcs_typ	*omc_typ::FlattenHiearchy()
// Flatten the hieerarchy by moving all nodes up to the root; delete empty nodes.
// Use to test MoveDown() operation.
{
	Int4	id,t,i,j,n,x,p,gp,*Parent,num_changed;
	double	lpr0,lpr1,DeltaLLR;
	mcs_typ *rtn_mcs=0,*old_mcs=0,*tmp_mcs=0;
	set_typ	subtree=0;

	tmp_mcs=this->RtnCopy(); // MoveUp will change mcs!!! Need to save a copy!
        old_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	for(j=2; j < Hpt->NumSets(); j++){
            if(mcs->SetCard(j) == 0){
                tmp_mcs=DeleteNode(Hpt->ItoSetID(j)); assert(tmp_mcs);
                DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt( ); j=1;
            }
	}
	do {
	    assert(Hpt->IsTree(Parent));
	    for(num_changed=0,i=2; i < Hpt->NumSets(); i++){
		p=Parent[i]; assert(Hpt->TypeOfSet(p) == '?');
		if(p == 1) continue;
		gp = Parent[p]; assert(Hpt->TypeOfSet(gp) == '?');
		if(gp == 1){
		   fprintf(stderr,"Moving node %d (\"%s\") up from %d (\"%s\") to %d (\"%s\")\n",
			i, Hpt->SetName(i),p,Hpt->SetName(p),gp,Hpt->SetName(gp));
		   subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
		   ReSetSubTree(subtree,i); // Hpt->Put(stderr,FALSE);
		   if(!mcs->MoveUp(gp,p,i,subtree)) print_error("MoveUp( ) operation failed");
		   // Caution: the BestMCS is not changed;
		   // make sure Best is not restored prior to unscrambling!!!
		   NilSet(subtree); 
		   // fprintf(stderr,"NumSets = %d; NumBPPS = %d\n",Hpt->NumSets(),Hpt->NumBPPS());
		   // Hpt->Put(stderr);
		   // Hpt->PutSorted(stderr);
		   rtn_mcs=this->Operate('u',0,0);	// unscramble hpt...
		   DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
		   // Hpt->PutSorted(stderr);
		   // mcs->SampleColumns(); 
		   num_changed++;
		   break;
		}
	    } free(Parent);
	} while(num_changed > 0);
	mcs->SampleColumns(); 
#if 0
	for(Int4 j=2; j < Hpt->NumSets(); j++){
            if(mcs->SetCard(j) == 0){
                tmp_mcs=DeleteNode(Hpt->ItoSetID(j)); assert(tmp_mcs);
                DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt( );
            }
        }
#endif
	rtn_mcs=mcs; mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
	return rtn_mcs;
}


