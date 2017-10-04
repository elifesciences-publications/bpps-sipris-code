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
#include <iostream>

void	omc_typ::PrintOperationHeader(FILE *fp,char *msg, Int4 iter)
{
      if(fp){
      	 double D; if(mcs->RtnBestLPR() > Best_LPR) D=mcs->RtnBestLPR(); else D=Best_LPR;
	 fprintf(fp,"################ %s (%d: LLR = %.2f (%d/%d failed)) ###############\n",
			msg,iter,D,mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR),Hpt->NumBPPS());
      }
}

double	omc_typ::TrimNodes(FILE *fp,double Temp, Int4 iter)
{
      if(Hpt->NumSets() <= 3) return this->CalcLLR(mcs);	// don't delete down to one node!! 
      Int4 n; PrintOperationHeader(fp,"Deleting nodes",iter); t_stat=0;
      do {
         if((n=this->DeleteNodes(Temp)) > 0){
            if(fp) fprintf(fp,"Initial sampling: %d rejected nodes deleted.\n",n);
            this->AdjustSetIDs(); this->RestoreBest(); 
	    // mcs->PutHyperPartition(); // mcs->PutHyperPartition(stderr); 
	 } t_stat += n;
      } while(n > 0);
      return this->CalcLLR(mcs); 
}

double	omc_typ::MergeSiblings(FILE *fp,Int4 iter,double Temp)
{
	if(Hpt->NumSets() <= 3) return this->CalcLLR(mcs);	// don't fuse down to one node!! 
	Int4 n; PrintOperationHeader(fp,"Merging sibling nodes",iter); t_stat=0;
	if((status=this->FuseNodes( )) > 0){  
		if(fp) fprintf(fp,"FuseNodes() reduced tree by %d nodes.\n",status);
		Sample(Temp,2,2,2,2,mcs,0.03);
		this->RestoreBest(); // mcs->PutHyperPartition(); 
		return TrimNodes(fp,300,iter); 
		// return this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} return this->CalcLLR(mcs);
}

double	omc_typ::AddInternalNodes(FILE *fp,double Temp, Int4 iter)
{
	if(Hpt->NumSets() <= 2) return this->CalcLLR(mcs);
     	PrintOperationHeader(fp,"Inserting nodes",iter);
	if((status=this->InsertInternalNodes(fp,Temp)) > 0){  
		this->RestoreBest(); // mcs->PutHyperPartition(); 
		return TrimNodes(fp,300,iter); 
		// return this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} else return this->CalcLLR(mcs);
}

double	omc_typ::GrowLeaves(FILE *fp,double Temp, Int4 iter, char Action)
{
#if 0
	static Int4 num_calls=0;
	PrintOperationHeader(fp,"Sampling leaves",iter);
	if(num_calls==0 && (status=this->AddLeaves(Temp,'a')) > 0){  
		if(fp) fprintf(fp,"GrowLeaves('a') added %d nodes.\n",status);
	} num_calls++;
#endif
	assert(Action == 'A' || Action == 'a' || Action == 'Z' || Action == 'z');
	if((status=this->AddLeaves(Temp,Action)) > 0){  
		if(fp) fprintf(fp,"GrowLeaves('%c') added %d nodes.\n",Action,status);
		Sample(Temp,2,2,2,2,mcs,0.03);
		this->RestoreBest(); // mcs->PutHyperPartition(); 
		return TrimNodes(fp,300,iter); 
		// return this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} else return this->CalcLLR(mcs);
}

double	omc_typ::RaiseBranches(FILE *fp,double Temp, Int4 iter)
{
	if(Hpt->NumSets() <= 2) return this->CalcLLR(mcs);
     	PrintOperationHeader(fp,"Moving nodes Up",iter); // mcs->PutHyperPartition(stderr);
	if((status=MoveNodesUp(Temp)) > 0){
	    this->RestoreBest(); // mcs->PutHyperPartition(); 
	    return TrimNodes(fp,Temp,iter); 
	    // this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} return this->CalcLLR(mcs);
}

double	omc_typ::LowerBranches(FILE *fp,double Temp, Int4 iter)
{
	if(Hpt->NumSets() <= 2) return this->CalcLLR(mcs);
     	PrintOperationHeader(fp,"Moving nodes Down",iter); status=0;
	// mcs->PutHyperPartition(stderr);
	if(Hpt->NumSets() > 2 && (status=MoveNodesDown(Temp)) > 0){
	    this->RestoreBest(); return TrimNodes(fp,Temp,iter);
	    // this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} return this->CalcLLR(mcs);
}

Int4    omc_typ::ResurrectRejects(FILE *fp,double Temp,Int4 iter)
// If there are MinimumSplitSize or more rejected seqs, then attempt to resurrect them.
{
        // if(mcs->SetCard(Hpt->NumSets()) < MinimumSplitSize) return 0;
        if(mcs->NumRejected() < MinimumSplitSize){ status=t_stat=0; return 0; }

	double D,d;
	if(best_mcs != 0) RestoreFinalBest( );  // sets best_mcs=0; Best_LPR=0.0;
	// mcs->RestoreBest();
	mcs_typ *xmcs=RtnCopy(); 	// save a copy of best_mcs.
	mcs->MoveSeqs(Hpt->NumSets(),1);	// move rejected sequences to root node.
	mcs->SampleColumns();  
	mcs->CalcTotalLPR(); mcs->SaveBest=TRUE; mcs->StoreBest();
        PrintOperationHeader(fp,"Resurrecting Rejects",iter);
        xmcs->PutHyperPartition(stderr); mcs->PutHyperPartition(stderr);
        if((status=this->AddLeaves(Temp,'e')) > 0){
                if(fp) fprintf(fp,"Resurrecting rejects added %d nodes.\n",status);
                // Sample(Temp,2,2,2,2,mcs,0.03);
                this->RestoreBest(); // mcs->PutHyperPartition(); 
		D=TrimNodes(fp,300,iter);
                // this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
        } else { D=this->CalcLLR(mcs); }
        d=CalcLLR(xmcs);  // old mcs...
	Int4 n=mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
        fprintf(stderr,"NewLLR=%.2f; OldLLR = %.2f.\n",D,d);
        if(D > d && n == 0){	// new mcs is better than old xmcs!!!
            fprintf(stderr,"!!!!!!!!!!!!!!!!!! Resurrection sucessful. !!!!!!!!!!!!!!!!!!!!\n");
            mcs->PutHyperPartition(stderr); DeleteMCS(xmcs); 
        } else { DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt(); this->CalcLLR(mcs); }
	// this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
        return 1;
}

void    omc_typ::Sample(double T,Int4 p1,Int4 p2,Int4 p3,Int4 p4,mcs_typ *&xmcs,
							double min_gain, set_typ st)
// p1=IterStart,p2=IterEvolve,p3=NumRounds,p4=ColSampleStart
// Keep iterating as long as the fractional gain in the LLR is >= min_gain.
{
    if(xmcs==0){ 
	mcs->ResetTemperature(T); mcs->Sample(p1,p2,p3,p4); mcs->RestoreBest(); CalcLLR(mcs);
    } else {
      if(st) PutSet(stderr,st);
      double R,D,d=xmcs->CalcTotalLPR( );
      xmcs->ResetTemperature(T); xmcs->SampledSet=st;
      do {
         xmcs->Sample(p1,p2,p3,p4); 
	 // CheckSSTvsCSQ(xmcs);	// DEBUG...
#if 0
	 xmcs->PutHyperPartition(stderr);
	 xmcs->PutPttrnVsConsSeq(stderr,"omc->Sample() debug 0");
#endif
	 D=CalcLLR(xmcs); // calls xmcs->RestoreBest();
	 if(xmcs == mcs){
		Int4 n=xmcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
		if(n > 0){ D=TrimNodes(stderr,T,1); xmcs=mcs; } // xmcs was deleted...
	 }
	 R=(D-d)/d;
         fprintf(stderr," ////// LLR=%.1f (last=%.1f + %.1f%c; best=%.1f) //////\n",
                                        D,d,100*R,'%',Best_LPR); d=D;
      } while(R >= min_gain); xmcs->SampledSet=0;
    }
}

Int4    omc_typ::DeleteWorstNode(double Temp)
{
        double worst=DBL_MAX,*Map=mcs->RtnSubLPR( );
        for(Int4 i=2; i <= Hpt->NumBPPS(); i++){ if(Map[i] < worst){ worst = Map[i]; } }
        return DeleteNodes(300.0,worst+0.1);
}

Int4    omc_typ::PruneTree(FILE *fp,double Temp, Int4 iter,Int4 Retain)
{
     	PrintOperationHeader(fp,"Pruning Tree",iter); 
	while(Hpt->NumBPPS() > Retain) DeleteWorstNode(Temp);
	RePartitionSeqs(fp,iter); // RandomizeSeqAssign(fp,50,iter); 
	// p1=IterStart,p2=IterEvolve,p3=NumRounds,p4=ColSampleStart
	// this->Sample(300,2,2,3,2,mcs,0.05);  
	mcs->SampleColumns( ); RecordChange("\nPruneTree");
}

double	omc_typ::RandomizeSeqAssign(FILE *fp, Int4 F, Int4 iter)
// Randomly assigns one out of F sequences to a set.
{
     	PrintOperationHeader(fp,"Randomizing",iter);
	// mcs->PutHyperPartition(stderr);
	mcs_typ *rtn_mcs=Randomize(F); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
	return this->CalcLLR(mcs);
}

double	omc_typ::RePartitionSeqs(FILE *fp, Int4 iter)
// Repartitions sequences based on the consensus for each subgroup.
{
	mcs_typ *rtn_mcs;
     	PrintOperationHeader(fp,"RePartitioning",iter);
	// mcs->PutHyperPartition(stderr);
	// mcs_typ *rtn_mcs=RePartition( ); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
	if(iter == 0) rtn_mcs=RePartition( ); 
	else if(iter == 1) rtn_mcs=BstRePartition( ); 
	else rtn_mcs=RandSqRePartition( ); 
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
	return this->CalcLLR(mcs);
}

void	omc_typ::MixItUp(char mode, Int4 iter)
{
     switch(mode){
	case 'B': case 'b': 
	  { RePartitionSeqs(stderr,iter); RandomizeSeqAssign(stderr,10,iter); } break;
	case 'R': case 'r': { RandomizeSeqAssign(stderr,10,iter); } break;
	case 'P': case 'p': { RePartitionSeqs(stderr,iter); break; } break;
	default: break;
     }
     switch(mode){
	case 'B': case 'R': case 'P': mcs->SampleColumns( ); break;
	default: break;
     }
}

Int4	omc_typ::RunIter(Int4 iter,Int4 jter,double Temp, double temp)
{
	Int4	n,No=0;	// number of operations...
	mcs->RestoreBest();  
	if(iter==1){ TrimNodes(stderr,Temp,iter); }
	else if(jter > 1){ this->RestoreVeryBest(); }
	status=t_stat=0; log_time=time(NULL);  // Setting for RecordChange();

	GrowLeaves(stderr,temp, iter); n=RecordChange("\nGrowLeaves");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
#if 1
	MergeSiblings(stderr,iter,temp); n=RecordChange("FuseSiblings");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
#endif
	ResurrectRejects(stderr,temp,iter); n=RecordChange("ResurrectRejects");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
	RaiseBranches(stderr,temp, iter); n=RecordChange("RaiseBranches");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
	LowerBranches(stderr,temp, iter); n=RecordChange("LowerBranches");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
	// if(No <= 0) RePartitionSeqs(stderr,iter); 
	AddInternalNodes(stderr,temp,iter); n=RecordChange("AddInternalNodes");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
	GrowLeaves(stderr,temp, iter,'a'); n=RecordChange("AddLeaves");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
#if 1
	MergeSiblings(stderr,iter,temp); n=RecordChange("FuseSiblings");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
#endif
	RaiseBranches(stderr,temp, iter); n=RecordChange("RaiseBranches");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }
	LowerBranches(stderr,temp, iter); n=RecordChange("LowerBranches");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint();}
	RmInternalNodes(stderr,temp, iter); n=RecordChange("RmInternalNodes");
        if(n > 0){ this->RestoreVeryBest(); No+=n; this->PutToggleCheckPoint(); }

        fprintf(stderr,"################### Sampling partitions/patterns ###############\n");
	this->Sample(Temp,2,2,3,2,mcs,0.01); 
        this->RestoreVeryBest(); 
	mcs->PutHyperPartition(stderr); // this->PutToggleCheckPoint();
	return No;
}

double	omc_typ::RmInternalNodes(FILE *fp,double Temp, Int4 iter)
{
     	PrintOperationHeader(fp,"Removing internal nodes",iter);
	if((status=this->DeleteInternalNodes(fp,Temp)) > 0){  
		this->RestoreBest(); // mcs->PutHyperPartition(); 
		return TrimNodes(fp,300,iter); 
		// return this->RestoreVeryBest(); // last_lpr=this->RevertToBest();
	} else return this->CalcLLR(mcs);
}

Int4	omc_typ::DeleteInternalNodes(FILE *fp,double Temp)
// Try deleting internal nodes to see whether this improves the LLR.
{
	Int4	i,id,n=0;
	double	D,d;

	for(i=2; i <= Hpt->NumBPPS(); i++){
		if(Hpt->TypeOfSet(i) == '!') continue;
		id=Hpt->ItoSetID(i); 
		if(stable && MemberSet(id,stable)) continue;
		mcs->SaveBest=TRUE; D=mcs->CalcTotalLPR( );
		mcs_typ *rtn_mcs = this->DeleteNode(id,Temp);
		d=rtn_mcs->CalcTotalLPR( );
		if(d > D){
		   fprintf(stderr,"Deleting internal node %d ('Set%d') improves LLR from %.2f to %.2f\n",
				i,id,D,d); n++; i=1; // start over...
		   rtn_mcs->PutHyperPartition(stderr);
		   DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); this->CalcLLR(mcs);
		} else {
		   fprintf(stderr,"Failed: reverting back to originial hierarchy.\n");
		   DeleteMCS(rtn_mcs);
		}
	} return n;
}

double	omc_typ::SproutLeaves(FILE *fp,double Temp, Int4 iter)
{
	Int4 MinimumSetSize_bak=MinimumSetSize,MinimumSplitSize_bak=MinimumSplitSize;
	MinimumSetSize = (Int4) ceil((double)SizeTrueMainCMA * 0.05);
	MinimumSetSize = MAXIMUM(Int4,MinimumSetSize,1000);
	MinimumSplitSize = 2* MinimumSetSize;
	fprintf(stderr,"MinimumSetSize = %d; MinimumSplitSize = %d; SizeTrueMainCMA = %d\n",
		MinimumSetSize,MinimumSplitSize,SizeTrueMainCMA);
	if(MinimumSplitSize < SizeTrueMainCMA){
	   if(Hpt->NumSets() == 2){
		mcs->PutHyperPartition(stderr);
		AddLeaves(300.0, 'z');	// Zero decendants 
		mcs->PutHyperPartition(stderr);
	   }
	   if(Hpt->NumSets() == 2){
		AddLeaves(300.0, 'Z');	// Zero decendants 
		mcs->PutHyperPartition(stderr);
	   }
	} MinimumSetSize=MinimumSetSize_bak; MinimumSplitSize=MinimumSplitSize_bak;
	// exit(1);
        return this->CalcLLR(mcs); 
}

Int4	omc_typ::Run(FILE *fpmma, FILE *fphpt, FILE *ptrnfp)
{
	Int4	i,j,n,N,s,iter=0,result=0,No;
	double	best_lpr,last_lpr,d,D;

	assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
	mcs->NoFailureMode=TRUE;  // if node configuration subLPR <= 0.0, then reject it.
	mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
#if 1
	if(NumForced > 0){ 	// Don't move this below PrintMode!!
		status=this->InsertInternalNodes(stderr,300);
		mcs->PutHyperPartition(stderr);
		sprintf(str,"%s_C",infile); this->PutCheckPoint(str,FALSE);
		mcs->Put(TRUE);
		// mcs->Sample(2,2,2,2); // IterStart IterEvolve NumRounds ColSampleStart.
		mcs->PutHyperPartition(stderr);
		this->RestoreBest(); 
		mcs->PutHyperPartition(stderr);
		this->CalcLLR(mcs); this->RestoreFinalBest(); 
		// mcs->PutHyperPartition(stderr);
		// this->Put( );	// outputs checkpoint files by default.
		return result;
	}
	if(DelNode){ 		// Don't move this below PrintMode!!
		Int4 id=Hpt->ItoSetID(DelNode); 
		// if(stable && MemberSet(id,stable)) continue;
		mcs_typ *rtn_mcs = this->DeleteNode(id,300);
		DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
		sprintf(str,"%s_B",infile); this->PutCheckPoint(str,FALSE);
		mcs->PutHyperPartition(stderr);
		mcs->Put(TRUE);
		return 1;	
	}
	if(MvNode){ return TestSubRoutine('L');	}	// Don't move this below PrintMode!!
#endif
	if(Evolve) mcs->DoEvolve(); else mcs->DoNotEvolve();
	if(PrintMode && PrintMode != 'T') { PrintOptions(fpmma,fphpt,ptrnfp); return 1; }
	if(XC_line){ free(this->FindCrossConserved(XC_line)); return 1; }
	if(test_mode) return TestSubRoutine(test_mode);	// tests individual subroutines (omc_debug.cc).
        fprintf(stderr,"################### Sampling partitions/patterns ###############\n");
	if(flatten_hpt) this->FlattenHierarchy( );
	else if(random_start > 0) this->RandomStart(random_start);
	else {
#if 1
	  if(Hpt->NumSets() > 2){ TrimNodes(stderr,300.0,iter); }
	  // mcs->DoEvolve(); mcs->SampleColumns( );
	  // mcs->Sample(2,2,2,2); // IterStart IterEvolve NumRounds ColSampleStart.
	  SproutLeaves(stderr,300,iter);
#elif 0
	  mcs->PutHyperPartition(stderr);
	  AddLeaves(300.0, 'z');	// Zero decendants 
	  mcs->PutHyperPartition(stderr);
	  AddLeaves(300.0, 'Z');	// Zero decendants 
	  mcs->PutHyperPartition(stderr);
	  exit(1);
#endif
        }
	if(InitMode == 'z' && Hpt->NumSets() <= 3){
	   SetStringency(10); iter=0;
	   mcs_typ *rtn_mcs=0,*old_mcs;
	   do {	// get a starting hpt without failed nodes...
	  	if(logfp){ mcs->PutHyperPartition(logfp); fflush(logfp); }
		// mcs->Sample(2,2,3,3); // IterStart IterEvolve NumRounds ColSampleStart.
		mcs->Sample(2,2,2,2); // IterStart IterEvolve NumRounds ColSampleStart.
	  	if(logfp){ mcs->PutHyperPartition(logfp); fflush(logfp); }
                d=mcs->CalcTotalLPR(); 
		if(mcs->NumFailedNodes() == 0) break;
		// if(mcs->NumberFailedNodes() == 0) break;
		if(Hpt->NumSets() == 3){
		   old_mcs=mcs;
		   mcs->PutHyperPartition(stderr);
		   D= (double) mcs->NumRejected()/(double) NumSeqsCMSA(TrueMainCMA);
		   if(D > 0.25){
        		PrintOperationHeader(stderr,"Resurrecting Rejects",iter);
			mcs->MoveSeqs(1,2);	// move root sequences to leaf node.
			mcs->MoveSeqs(Hpt->NumSets(),1); // move rejects to root node.
#if 0
		   	mcs->PutHyperPartition(stderr);
			mcs->LoadUpColumns(); mcs->DoEvolve(); mcs->SampleColumns( );
		   	mcs->PutHyperPartition(stderr);
			mcs->Sample(2,2,3,3); // IterStart IterEvolve NumRounds ColSampleStart.
#endif
		   }
		   if(mcs->NumFailedNodes() == 0) break;
		   // if(mcs->NumberFailedNodes() == 0) break;
		   mcs->PutHyperPartition(stderr);
		   rtn_mcs=this->DeleteNode(2); mcs=rtn_mcs; Hpt=mcs->GetHpt();
		   mcs->PutHyperPartition(stderr);
Hpt->Put(stderr);
		   // rtn_mcs=this->AddLeaf(1,2); 
		   rtn_mcs=this->Operate('a',1,2); 
		   if(rtn_mcs){
			DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
		   } else { print_error("failed to find a significant pattern/partition"); }
		   mcs->PutHyperPartition(stderr);
	        }
// exit(1);
		mcs->Sample(4,2,4,3); // IterStart IterEvolve NumRounds ColSampleStart.
		if(mcs->NumFailedNodes() == 0) break;
		// if(mcs->NumberFailedNodes() == 0) break;
		if(Hpt->NumSets() == 3){
		   rtn_mcs=this->DeleteNode(2); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
		   mcs->PutHyperPartition(stderr);
		   // rtn_mcs=this->AddLeaf(1,2); 
		   rtn_mcs=this->Operate('a',1,2); 
		   if(rtn_mcs){
			DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); 
		   } else { print_error("failed to find a significant pattern/partition"); }
		   mcs->PutHyperPartition(stderr);
#if 1		// fix core dump with contaminated sequence sets.
		} else if(Hpt->NumSets() == 2 && iter > 0){
			double d1=(double)mcs->SetCard(1);
			double dR=(double)mcs->NumRejected();
			double d= d1/(d1+dR);
			if(d > 0.25){
			     FILE *fp=open_file(infile,"_clean.cma","w");
			     set_typ *xset=mcs->CopyOfSeqSets();
			     PutInSetCMSA(fp,xset[1], TrueMainCMA); fclose(fp);
			     fprintf(stderr,"omcBPPS input alignment is grossly contaminated.\n");
			     fprintf(stderr,"Saved decontaminated alignment as %s_clean.cma.\n",infile);
			     print_error("FATAL: try rerunning on cleaned up alignment.");
			} else print_error("omcBPPS failed to find a hierarchy.");
#endif
		}
	  	if(logfp){ mcs->PutHyperPartition(logfp); fflush(logfp); }
		// mcs->SampledColumn=0;
		mcs->LoadUpColumns();
		// MixItUp('p', iter); // == RePartitionSeqs(stderr,iter); mcs->SampleColumns( );
		MixItUp('p',iter); // == RePartitionSeqs(stderr,iter); mcs->SampleColumns( );
		// mcs->LoadUpColumns();
	  	if(logfp){ mcs->PutHyperPartition(logfp); fflush(logfp); }
		mcs->DoEvolve(); iter++;
	    } while(iter==1 || (time(NULL)-time1) < 120);
	} else {
	   // if(stable) mcs->Sample(2,2,3,3); // IterStart IterEvolve NumRounds ColSampleStart.
	   mcs->SampleColumns( ); mcs->DoEvolve();  
	   // mcs->StoreBest(); mcs->RestoreBest();
	   // mcs->Sample(2,3,2,1); mcs->DoEvolve();  // mcs->Sample(2,2,2,2);
	   // if(Hpt->NumSets() <= 3) mcs->Sample(2,2,2,2);	// let csq evolve first.
	}
	SetStringency(stringency);
	mcs->CalcTotalLPR(); mcs->NoFailureMode=FALSE; mcs->StoreBest(); mcs->RestoreBest(); 
	mcs->NoFailureMode=TRUE;
	double Temp=300.0,temp=300.0; 
	best_lpr=this->CalcLLR(mcs);
#if 1	// dealing with lots of rejects...
	if(PrintMode == 'T') { PrintOptions( ); return 1; }
	D= (double) mcs->NumRejected()/(double) NumSeqsCMSA(TrueMainCMA);
	if(D > 0.25) ResurrectRejects(stderr,300,1);
#endif
	Int4	MaxIter=1,Remove=2,Retain,Last;
	for(i=1; i <= MaxIter; i++){
	  switch(i){
	    case 1: break;
	    case 2: case 3:
	        this->RestoreVeryBest(); 
		TrimNodes(stderr,300.0,iter); // need to start all over again...
		RePartitionSeqs(stderr,0); mcs->SampleColumns( );
		break;
	    default: break;
	  } 
	  do {
	    last_lpr=best_lpr;
	    iter++; No=RunIter(iter,0,Temp,temp);  // don't restore very best
	    this->RestoreVeryBest(); best_lpr=this->CalcLLR(mcs);
	    // Temp = MAXIMUM(double,Temp-25,150);
	    d=100.0*((best_lpr/last_lpr) - 1.0);
	    if(d <= 0.0) break; // worse or exactly the same ratio.
	    else if(verbose) { sprintf(str,"%s_i",infile); PutCheckPoint(str); }
	    else { sprintf(str,"%s",infile); PutCheckPoint(str); }
	    // else { sprintf(str,"%s_chk",infile); PutCheckPoint(str); }
	    // else { this->PutToggleCheckPoint(); }
	    // else { sprintf(str,"%s_%d_chk",infile,iter); PutCheckPoint(str); }
	  } while(No > 0  || d >= 0.2);	// require increase in LLR by 1/5th percent.
	  // if(i == 1 && Hpt->NumSets() <= 10) // only 10 nodes generated; try harder...
	// 	{ if((time(NULL)-time1) < 120) i=0; }
	  if(i < MaxIter){
		if(logfp) fprintf(logfp,"\n** %d. Reinitializing the current hierachy **\n", i);
	  } 
	}
	this->Sample(150,2,2,2,2,mcs,0.05); 
	this->Sample(50,2,2,3,2,mcs,0.01); 
	fprintf(stderr,"############### Creating optimum hierarchy ###############\n");
	this->RestoreFinalBest(); 
	// mcs->PutHyperPartition(); 
	mcs->PutHyperPartition(stderr); 
	this->PutAsBestHpt(mcs);

        FILE *fp=0; if(verbose) fp=open_file(infile,"_sq.lpr","w");
        FILE *hfp=0; if(verbose) hfp=open_file(infile,"_sq.hst","w");
        FILE *sfp=open_file(infile,"_usr.sma","w");
        mcs->PutSeqCntrb(fp,hfp,sfp);
        if(fp) fclose(fp); if(hfp) fclose(hfp); if(sfp) fclose(sfp); 

	mcs_typ *rtn_mcs=this->Optimize(); 
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt(); this->CalcLLR(mcs); 
	this->Put( );	// outputs checkpoint files by default.
	this->PrintTime(stderr);
	return result;
}

void	omc_typ::CompareBPPS_BILD( )
//****************** Compare BPPS with BILD scores. *******************
{
	hpt_typ *Hpt=mcs->GetHpt();
	char	dms_mode='F',wt_factor=100;	// 
h_type	HG1=Histogram("BILD scores",-50,100,1);
h_type	HG2=Histogram("BPPS scores",-50,100,1);
	Int4	pernats=1000;
	dms_typ	*dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	che_typ **che=mcs->RtnChe( );	// need to create this function.
	UInt4	*WtCntsFB; NEW(WtCntsFB,nAlpha(AB)+3,UInt4);
        mcs->PutHyperPartition(stdout); 
#if 0	// compare times...
	Int4 Iters=2000000,time1=time(NULL);
    for(Int4 X=1; X <= Iters; X++)
	for(Int4 n=1; n < Hpt->NumSets(); n++){ bpps_typ *bpps=mcs->RtnBPPS(n); }
	fprintf(stderr, "\tBPPS time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

	time1=time(NULL);
exit(1);
    for(Int4 X=1; X <= Iters; X++)
	for(Int4 n=1; n < Hpt->NumSets(); n++){
	    sst_typ *sst=mcs->RtnCopyOfSST(n);
	    double d,Dd,FG,BG,DD,dD;
	    UInt4  **WtCntsFG=che[n]->GetResWtsFG();
	    UInt4  **WtCntsBG=che[n]->GetResWtsBG();
	    // fprintf(stdout,"\n=========== %d.%s ==========\n",n,Hpt->SetName(n));
	    for(Int4 i=1; i <= mcs->RtnLengthMainCMSA(); i++){
		if(sst[i] != 0) FG=dms->bild(WtCntsFG[i]);
	    }
	}
	fprintf(stderr, "\tBILD time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	exit(1);
	return;
#endif
	for(Int4 n=1; n < Hpt->NumSets(); n++){
	    sst_typ *sst=mcs->RtnCopyOfSST(n);
	    // sst_typ *sst=mcs->RtnCopyOfBestSST(n);
	    double d,Dd,FG,BG,DD,dD,**dd=mcs->GetResEvals(n);
	    bpps_typ *bpps=mcs->RtnBPPS(n);
	    // bpps->RtnWtFactor(); // assert(WtFact==100);
	    // UInt4	che[n]->GetResWtFG(i,r);
	    UInt4  **WtCntsFG=che[n]->GetResWtsFG();
	    UInt4  **WtCntsBG=che[n]->GetResWtsBG();
	    fprintf(stdout,"\n=========== %d.%s ==========\n",n,Hpt->SetName(n));
	    for(Int4 i=1; i <= mcs->RtnLengthMainCMSA(); i++){

FG=dms->bild(WtCntsFG[i]);
BG=dms->bild(WtCntsBG[i]);
double	TypScore=dms->Typical(WtCntsFG[i]);
	double totalFG=0.0,totalBG=0.0,ffg,fbg,Rt,rt;
	for(Int4 r=0; r <= nAlpha(AB); r++){
		WtCntsFB[r]=WtCntsFG[i][r] + WtCntsBG[i][r];
		totalFG += (double) WtCntsFG[i][r];
		totalBG += (double) WtCntsBG[i][r];
	} Dd=dms->bild(WtCntsFB); DD = (FG+BG) - Dd; Dd = Dd*(totalFG/(totalFG+totalBG));
		if(sst[i]==0){
		    if(DD < 5.0 || (FG-Dd) < DD) continue;
		    if(TypScore <= 0.0) continue;
		    IncdHist(TypScore, HG1);
		    fprintf(stdout,"  *** %d=%.2f[%.2f]{%.2f}: ",i,DD,FG-Dd,TypScore);
		    for(Int4 r=0; r <= nAlpha(AB); r++){
			fbg=((WtCntsBG[i][r]+1)/(totalBG+1));
			ffg=((WtCntsFG[i][r]+1)/(totalFG+1));
			Rt=ffg/fbg; 
			if(Rt > 1.0 && ffg > 0.04)
			  fprintf(stdout,"%c%.1f[%.2f] ",
				AlphaChar(r,AB),Rt,ffg);
			// fprintf(stdout,":%.1f:%.1f ",ffg,fbg);
		    }
		    for(Int4 r=0; r <= nAlpha(AB); r++){
			fbg=((WtCntsBG[i][r]+1)/(totalBG+1));
			ffg=((WtCntsFG[i][r]+1)/(totalFG+1));
			rt=fbg/ffg;
			if(rt > 1.0 && fbg > 0.04)
			  fprintf(stdout,"(%c%.1f[%.2f]) ",
				AlphaChar(r,AB),rt,fbg);
		    } fprintf(stdout,"\n"); continue;
		}

       		d=dms->sdotfs(WtCntsFG[i],WtCntsBG[i]);
		PutSST(stdout,sst[i],AB);
		for(Int4 r=1; r <= nAlpha(AB); r++){
		  if(MemSset(r,sst[i])){
		    IncdHist(TypScore, HG2);
		    fprintf(stdout,"%d=%.2f; P2P = %.2f[%.2f]{%.2f}: ",i,dd[i][r],DD,FG-Dd,TypScore);
		    for(Int4 r=0; r <= nAlpha(AB); r++){
			fbg=((WtCntsBG[i][r]+1)/(totalBG+1));
			ffg=((WtCntsFG[i][r]+1)/(totalFG+1));
			Rt=ffg/fbg; 
			if(Rt > 1.0 && ffg > 0.04)
			  fprintf(stdout,"%c%.1f[%.2f] ",AlphaChar(r,AB),Rt,ffg);
		    }
		    for(Int4 r=0; r <= nAlpha(AB); r++){
			fbg=((WtCntsBG[i][r]+1)/(totalBG+1));
			ffg=((WtCntsFG[i][r]+1)/(totalFG+1));
			rt=fbg/ffg;
			if(rt > 1.0 && fbg > 0.04)
			  fprintf(stdout,"(%c%.1f[%.2f]) ",AlphaChar(r,AB),rt,fbg);
		    } fprintf(stdout,"\n");
		    break;
		  }
		} free(dd[i]);
	    } free(dd); free(sst);
	} delete dms;
	PutHist(stdout,60,HG1); NilHist(HG1);
	PutHist(stdout,60,HG2); NilHist(HG2);
}

void	omc_typ::PrintOptions(FILE *fpmma, FILE *fphpt,FILE *ptrnfp)
{
	BooLean	ChangeNames=TRUE;
	// mcs->PutHyperPartition(stderr);
	switch (PrintMode) {
	   case 'B':
	     {
		FILE *sfp=open_file(infile,"_usr.sma","w");
        	mcs->PutSeqCntrb(0,0,sfp); fclose(sfp); 
	     } break;
	   case 'R':
	    {
		 mcs->SetFontSize(font_size); mcs->SetPageFormat(page_format);
		 if(NumHighlighted > 0) mcs->SetNumHighlighted(NumHighlighted);
		 if(NthSeqForDisplay > 0) mcs->SetNthSeqForDisplay(NthSeqForDisplay);
		 if(use_usr_sma) mcs->SetUserDisplaySet();
		 this->PrintRTF(FALSE,FALSE,KeyStart,KeyEnd,KeySetName);
	    } break;
	   case 'H':	// put map contributions.
		// mcs->PutHyperPartition(); 
		mcs->PutHyperPartition(stderr); 
	    break;
	   case 'h': mcs->PutSubHierarchyCntrb(stdout); break;
	   case 'C':	// put map contributions.
	     {
		ChangeNames=FALSE;
		if(NumHighlighted > 0) mcs->SetNumHighlighted(NumHighlighted);
		mcs->DoEvolve(); mcs->SampleColumns();
        	FILE *fp=open_file(infile,"_bst.out","w"); mcs->PutHyperPartition(fp); fclose(fp); 
		mcs->Put(FALSE);
		sprintf(str,"%s_C",infile); PutCheckPoint(str,ChangeNames);
		if(NthSeqForDisplay > 0) mcs->SetNthSeqForDisplay(NthSeqForDisplay);
		if(use_usr_sma) mcs->SetUserDisplaySet();
		this->PrintRTF(FALSE,FALSE,KeyStart,KeyEnd,KeySetName);
		// FILE *fp=
		// mcs->PutMapContributions(FILE *fp);
	     } break;
	   case 'c':	// Optimize and output contributions.
	     { 
		if(NumHighlighted > 0) mcs->SetNumHighlighted(NumHighlighted);
		mcs->DoEvolve(); mcs->SampleColumns(); 
		mcs->PutHyperPartition(stderr);
        	FILE *fp=0; 
		fp=open_file(infile,".tbl","w"); mcs->PutHyperPartition(fp); fclose(fp); 
		this->VerboseOff();
		this->Put(TRUE,fpmma,fphpt,ptrnfp);
#if 0
		if(NthSeqForDisplay > 0) mcs->SetNthSeqForDisplay(NthSeqForDisplay);
		if(use_usr_sma) mcs->SetUserDisplaySet();
		this->PrintRTF(FALSE,FALSE,KeyStart,KeyEnd,KeySetName);
#else
		ChangeNames=FALSE;
		sprintf(str,"%s_c",infile); PutCheckPoint(str,ChangeNames);
#endif
	     } break;
	   case 'O':	// Optimize and output contributions.
		ChangeNames=FALSE;
	   case 'o':	// Optimize and output contributions.
		ChangeNames=FALSE;	// changed this...
		// if(stable) mcs->SetStringency(MinimumLLR,MinimumSetSize,stable);
		if(PrintMode=='O'){ mcs->DoNotEvolve(); }
		else { mcs->DoEvolve(); mcs->SampleColumns(); }
		this->SetStringency(5); // don't delete anything!!!
#if 0
		mcs->ResetTemperature(300); mcs->Sample(2,2,2,2); mcs->RestoreBest(); 
		// Don't call CalcLLR(mcs);
		mcs->ResetTemperature(150); mcs->Sample(2,2,2,2); mcs->RestoreBest(); 
		mcs->ResetTemperature(50); mcs->Sample(2,2,3,2); mcs->RestoreBest(); 
#else
		// Don't call CalcLLR(mcs);
		// mcs->ResetTemperature(0.0); mcs->Sample(2,2,3,2); mcs->RestoreBest(); 
		mcs->ResetTemperature(0.0); mcs->Sample(2,2,2,2); mcs->RestoreBest(); 
#endif
		if(Hpt->SetName(1) == 0){
		   Int4 *P,i; 
		   if(Hpt->IsTree(P)){	// use the same names as before...
			assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
		   	for(i=1; i <= Hpt->NumBPPS(); i++){ Hpt->ReNameSet(i,Hpt->GrpName(i)); }
		   } free(P);
		} sprintf(str,"%s_o",infile); PutCheckPoint(str,ChangeNames);
	   case 'S':	// put sequence contributions.
	     {
        	FILE *fp=open_file(infile,"_sq.lpr","w");
		FILE *hfp=open_file(infile,"_sq.hst","w");
		FILE *sfp=open_file(infile,"_usr.sma","w");
#if 0	//****************** Compare BPPS with BILD scores. *******************
		this->CompareBPPS_BILD( );
#endif	//***************************************************************
        	mcs->PutSeqCntrb(fp,hfp,sfp);
        	fclose(fp); fclose(hfp); 
		if(sfp) fclose(sfp); 
        	fp=open_file(infile,"_bst.out","w"); mcs->PutHyperPartition(fp); fclose(fp); 
		if(PrintMode=='O') this->Put(TRUE);  // input fp not used!
		else mcs->Put(FALSE);
	     } break;
#if 0
	   case 's':	// print sequences along lineage to keynode.
	     {
		Int4	i,j,k,s,n,Ns;
		set_typ *set=this->RtnSets( );
		set_typ LineSet=Hpt->MkLineageSet(FocusNode);
		n=CardSet(LineSet);
		set_typ *SetX; NEW(SetX,n+3,set_typ);
		for(s=1,Ns=0; s < Hpt->NumSets(); s++){
		   if(MemberSet(s,LineSet)){
			Ns++; SetX[Ns]=Hpt->MkSubTreeSet(s);
		   }
		}
		for(s=1; s <= Ns; s++){
		   fprintf(stderr,"%d: ",s); PutSet(stderr,SetX[s]); fprintf(stderr,"\n");
		}
		exit(1);
	     } break;
#endif
	   case 't':	// remove failed nodes.
	     {
		double Temp=0.0;
		this->TrimNodes(stderr,Temp, 1);
		this->Put(TRUE);
	     } break;
	   case 'P':	// put sequence contributions.
	     {
		this->Put(TRUE);
	     } break;
	   case 'T':	// test option.
	     {
		mcs->ResetTemperature(300); mcs->Sample(2,2,2,2); mcs->RestoreBest(); 
		// mcs->PutHyperPartition(); // mcs->PutHyperPartition(stderr); 
		ResurrectRejects(stderr,300,1);
		// mcs->PutHyperPartition(); // mcs->PutHyperPartition(stderr); 
		// mcs->DoEvolve(); mcs->PutSetRelations(stderr); mcs->SampleColumns(); 
	     } break;
	   default:
	    break;
	}
}

