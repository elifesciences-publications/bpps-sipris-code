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

double	omc_typ::OnTheFlyMvDown(Int4 Mv, Int4 To, Int4 *Parent)
// SetFG & SetBG are initialized from calling environment and correspond to node Mv's current FG & BG.
// Don't modify SetFG & SetBG!!!
// TmpFG & TmpBG are used for computing New FG/BG.
{
        Int4    i,j,c,sibling,NumCols=mcs->RtnDefaultMaxCol();
	char	NdTyp;
        double  DD,D,dd,d,deltaLLR=0.0;
	sst_typ	*xsst;

	assert(To > 1 && To < Hpt->NumSets());
	assert(Mv > 1 && Mv < Hpt->NumSets());
	assert(To != Mv);
	assert(Parent[Mv] == Parent[To]);

	// Calculate node 'To's current contrast alignment LLR.
	if(Hpt->TypeOfSet(To) == '!') NdTyp='L'; else NdTyp='M';
	mcs->CopySubTreeSeqSet(To,TmpFG); mcs->CopyBkGrndSeqSet(To,TmpBG);
	xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,DD,NumCols,NdTyp); free(xsst); // Current LLR.
	sprintf(str,"%d --> %d: %d (%d vs %d; %.2f nats)",Mv,Parent[Mv],To,CardSet(TmpFG),CardSet(TmpBG),DD);

	// Calculate node 'To's proposed contrast alignment LLR.
	UnionSet(TmpFG,SetFG); // TmpFG = TmpFG U SetFG.
	IntersectNotSet(TmpBG,SetFG); // TmpBG=TmpBG & !SetFG = To's proposed BG (Set = Mv's current FG).
	xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,dd,NumCols,'M'); free(xsst); // Proposed LLR.
	deltaLLR = dd-DD;
	if(deltaLLR > 0.0) fprintf(stderr,"%s --> %d: %d (%d vs %d; %.2f nats) change = %.2f.\n",
				str,To,To,CardSet(TmpFG),CardSet(TmpBG),dd,deltaLLR);
	// compute influence on 'To's other children (if any) as well...
	for(c=2; c < Hpt->NumSets(); c++){     
	   if(Parent[c] != To) continue;	// Examine current children of proposed parent node. 
	   if(Hpt->TypeOfSet(c) == '!') NdTyp='L'; else NdTyp='M';
	   // c's current LLR.
	   mcs->CopySubTreeSeqSet(c,TmpFG); mcs->CopyBkGrndSeqSet(c,TmpBG);
	   xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,D,NumCols,NdTyp); free(xsst); // Current LLR.

	   sprintf(str,"    + %d (%d vs %d; %.2f nats)", c,CardSet(TmpFG),CardSet(TmpBG),D);

	   // c's proposed LLR.
	   UnionSet(TmpBG,SetFG); // TmpBG = TmpBG U SetFG. --> add moved node's subtree to BG.
	   xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,d,NumCols,NdTyp); free(xsst); // Current LLR.
	   deltaLLR += d-D;

	   if(deltaLLR > 0) 
	     fprintf(stderr,"%s --> (%d vs %d; %.2f nats) change = %.2f; sum = %.2f.\n",
				str,CardSet(TmpFG),CardSet(TmpBG),d,d-D,deltaLLR);
	}
	return deltaLLR;
}

Int4    omc_typ::MoveNodesDown(double Temp, double minLLR)
// by default: minLLR = 0.0. 
{
        Int4    Sibs,id,i,n,worst_i,worst_s,NumMoves=0,*Parent,NumCols=mcs->RtnDefaultMaxCol();
        double  dd,d,D,DD,worst;
	char	NdTyp;
	sst_typ	*xsst;
        set_typ skip = MakeSet(MaxNumNodesPlus + MaxNumNodes); ClearSet(skip);
	mcs->SaveBest=TRUE;
	char	tmpstr[200],TmpStr[200];
    do {
	assert(preHG == 0); preHG=Histogram("preliminary change in LLR",0,2000,100);
        mcs->RestoreBest(); // mcs->PutHyperPartition(); 
	Hpt->IsScrambledTree(Parent);
        for(worst_i=0,worst=DBL_MAX,i=2; i < Hpt->NumSets(); i++){
           id=Hpt->ItoSetID(i); if(MemberSet(id,skip)) continue;
	   if(stable && MemberSet(id,stable)) continue;
	   if(stable && MemberSet(Hpt->ItoSetID(Parent[i]),stable)) continue;
           for(Sibs=0,n=2; n < Hpt->NumSets(); n++){     // find sibling nodes...in tree format.
		if(Parent[i] == Parent[n]) Sibs++;	// count number of sibling nodes.
	   } if(Sibs == 0) continue;
if(MvNode && i!=MvNode) continue;
	   mcs->CopySubTreeSeqSet(i,SetFG);	// SetFG == current foreground for node i.
	   mcs->CopyBkGrndSeqSet(i,SetBG);
	   if(Hpt->TypeOfSet(i) == '!') NdTyp='L'; else NdTyp='M';
	   xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,D,NumCols,NdTyp); free(xsst);  // Current LLR.

	   sprintf(TmpStr,"old: %d --> %d: (%d vs %d) %.2f nats; ",
				i,Parent[i],CardSet(SetFG),CardSet(SetBG),D);

           for(n=2; n < Hpt->NumSets(); n++){     // find sibling nodes...in tree format.
		if(n == i) continue;
		if(Parent[i] != Parent[n]) continue;	// Sibling nodes only.
if(ToNode && n != ToNode) continue;
		mcs->CopySubTreeSeqSet(n,TmpBG);	// TmpSet == proposed BG == n's current FG
		xsst=this->GetOptPttrnLPR(0,SetFG,TmpBG,d,NumCols,NdTyp); free(xsst); // Proposed LLR.

		sprintf(tmpstr,"new:  --> %d: %d (%d vs %d) %.2f nats: change = %.2f ",
			n,i,CardSet(SetFG),CardSet(TmpBG),d,d-D);

		dd = this->OnTheFlyMvDown(i,n,Parent);
		DD = (d - D) + dd; 
// fprintf(stderr,"Move i=%d to n=%d; dd=%lf; DD = %lf\n",i,n,dd,DD);
if(MvNode) DD=1.0; else
		if((DD < 0.0)) continue;
#if 1
		fprintf(stderr,"%s %s - %.2f = %.2f nats.\n",TmpStr,tmpstr,dd,DD);
#endif
                { IncdHist(DD,preHG); if(-DD < worst){ worst=-DD; worst_i=i; worst_s=n; } }
           }
        }  free(Parent); PutHist(stderr,60,preHG); NilHist(preHG); preHG=0;
	if(worst_i > 0){
           id=Hpt->ItoSetID(worst_i);
           fprintf(stderr,"/////// Down: Worst = %d ('Set%d'): %.2f ///////\n",worst_i,id,-worst);
           mcs_typ *rtn_mcs = this->MoveDown(worst_i,worst_s,minLLR);
           if(rtn_mcs){
	        rtn_mcs->SampleColumns();
		rtn_mcs->PutHyperPartition(stderr); // exit (1);
                NumMoves++; DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); mcs->SaveBest=TRUE;
           } else { id=Hpt->ItoSetID(worst_i);  AddSet(id,skip); }  // don't try this node again..
        }
if(MvNode) return 1;
    } while(worst_i > 0); NilSet(skip); return NumMoves;
}

mcs_typ *omc_typ::MoveDown(Int4 node, Int4 sibling, double MinDeltaLLR)
// Move 'node' down from it's parent to 'sibling'.
{
        Int4    id,t,i,j,x,n,gp,*Parent,NumChild_P=0;
        double  lpr0,lpr1,DeltaLLR;
	mcs_typ *rtn_mcs=0,*old_mcs=0,*tmp_mcs=0;

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
        // lpr0=mcs->RtnBestLPR();
        lpr0=mcs->CalcTotalLPR();
        fprintf(stderr,"LLR = %.3f (gp=%d p=%d node = %d)\n",lpr0,gp,sibling,node);
        fprintf(stderr," ... moving node %d down from node %d to node %d\n",node,gp,sibling);

	tmp_mcs=this->Operate('C',0,0);   // store the current best state.
	old_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	
        ReSetSubTree(subtree,node); // Hpt->Put(stderr,FALSE);
        if(!mcs->MoveDown(gp,sibling,node,subtree)) print_error("omc_typ::MoveDown( ) failed");
        NilSet(subtree); 
        rtn_mcs=this->Operate('u',0,0);   // unscramble mcs without optimizing...
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
	// Hpt->Put(stderr); Hpt=mcs->GetHpt( );  // make sure Hpt corresponds to mcs after MoveDown...
	mcs->SampleColumns();
	mcs->PutHyperPartition(stderr); // exit (1);
        lpr1=mcs->CalcTotalLPR( );
        // if(lpr1 <= 0 || (lpr1/lpr0) < 0.95) lpr1=mcs->SampleColumns( );
        fprintf(stderr,"Initial LLR = %.2f; old LLR = %.2f (diff = %.3f)\n",lpr1,lpr0,lpr1-lpr0);
if(MvNode) lpr0=1.0; else
        if(lpr1 <= 0 || lpr0 > 0 && (lpr1/lpr0) < 0.95){
		DeleteMCS(mcs); mcs=old_mcs; Hpt=mcs->GetHpt( ); // this->RevertToBest();  
                fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
                return 0;
        }

#if 0	// Delete this code after testing.
        rtn_mcs=this->Operate('u',0,0);   // unscramble mcs without optimizing...
	DeleteMCS(mcs); mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
	rtn_mcs->SampleColumns();
	rtn_mcs->PutHyperPartition(stderr); // exit (1);
#endif
	rtn_mcs=mcs; mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
	rtn_mcs->SaveBest=TRUE;
	n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
	// n=rtn_mcs->NumRawFailedNodes(0,0);	// for testing with FlattenHiearchy() routine.
        lpr1=rtn_mcs->CalcTotalLPR(0,FALSE);
        DeltaLLR = lpr1 - lpr0;
        fprintf(stderr,"Presampling LLR = %.3f (diff = %.3f; cutoff = %.3f; %d failed)\n",
		lpr1,DeltaLLR,MinDeltaLLR,n);
#if 1
        if(DeltaLLR > MinDeltaLLR && n == 0){ return rtn_mcs; }
#else
        if(DeltaLLR > MinDeltaLLR){ return rtn_mcs; }
#endif
        // mcs_typ *rtn_mcs=this->Operate('S',0,0);
	// Sample(300.0,2,2,2,2,rtn_mcs,0.03);
	rtn_mcs->Sample(2,2,2,2);
        lpr1=rtn_mcs->CalcTotalLPR();
        DeltaLLR = lpr1 - lpr0;
        fprintf(stderr,"Postsampling LLR = %.3f (diff = %.3f)\n",lpr1,DeltaLLR);
	n=rtn_mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
        if(DeltaLLR >= MinDeltaLLR && n == 0){ return rtn_mcs; }
        fprintf(stderr," ... reverting to best (LPR: %.3f back to %.3f)\n",lpr1,lpr0);
        DeleteMCS(rtn_mcs); return 0;
}

#if 0
mcs_typ	*omc_typ::DeepenHiearchy()
// Deepen the hieerarchy by moving all nodes down to leaves; delete empty nodes.
// Use to test MoveUp() operation.
{
	Int4	id,t,i,j,n,x,p,gp,*Parent,num_changed;
	double	lpr0,lpr1,DeltaLLR;
	mcs_typ *rtn_mcs=0,*old_mcs=0,*tmp_mcs=0;
	set_typ	subtree=0;

	tmp_mcs=this->RtnCopy(); // MoveUp will change mcs!!! Need to save a copy!
        old_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	do {
	    assert(Hpt->IsTree(Parent));
	    for(num_changed=0,i=2; i < Hpt->NumSets(); i++){
		p=Parent[i];
		if(p == 1) continue;
		gp = Parent[p];
		if(gp == 1){
		   subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
		   ReSetSubTree(subtree,i); // Hpt->Put(stderr,FALSE);
		   if(!mcs->MoveUp(gp,p,i,subtree)) print_error("MoveUp( ) operation failed");
		   NilSet(subtree); 
		   rtn_mcs=this->Operate('u',0,0);	// unscramble hpt...
		   DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
		   mcs->SampleColumns(); num_changed++;
		   break;
		}
	    } free(Parent);
	} while(num_changed > 0);
	for(Int4 j=2; j < Hpt->NumSets(); j++){
            if(mcs->SetCard(j) == 0){
                tmp_mcs=DeleteNode(Hpt->ItoSetID(j)); assert(tmp_mcs);
                DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt( );
            }
        } rtn_mcs=mcs; mcs=old_mcs; Hpt=mcs->GetHpt( );  // restore original mcs_typ
	return rtn_mcs;
}
#endif



