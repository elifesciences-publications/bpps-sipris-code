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

mcs_typ	*omc_typ::OptimizeInsertedNode(hsi_typ *hsi)
// Optimize final hsi based on above analysis...
// Returns a new hiearchy if sucsessful; otherwise returns null.
{
	FILE	*efp=0; efp=stderr;
        double	D,d;
	Int4	id,r,c,ID,pID,p,parent,*P=0,NumCols=mcs->RtnDefaultMaxCol();
	hpt_typ	*hpt=0;
	BooLean	changed;
	mcs_typ	*rtn_mcs=0,*xmcs=this->RtnCopy(); 

	// Extract hsi sets.
	set_typ	FG=hsi->RtnSetFG(), BG=hsi->RtnSetBG();
	set_typ	ChildFG=CopySet(FG),ChildBG=CopySet(BG),Changed=MakeSet(SetN(FG)); ClearSet(Changed);

	assert(Hpt->IsTree(P)); // Find the parent node...
	for(parent=0, c=1; c <= Hpt->NumBPPS(); c++){
		// if(MemberSet(Hpt->ItoSetID(c),ChildFG) && parent==0){ parent=P[c]; break; }
		if(MemberSet(Hpt->ItoSetID(c),ChildFG))
		   { if(parent==0) parent=P[c]; else assert(parent==P[c]);  }
	} free(P);	// pID is the Set identifier, not the index!
	assert(parent); pID=Hpt->ItoSetID(parent); assert(pID > 0 && pID < MaxNumNodes); 
	for(ID=1; ID <= MaxNumNodes; ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID
	assert(ID > 0 && ID <= MaxNumNodes);
	
     do {
     	if(CardSet(ChildFG) < 2 || CardSet(ChildBG) < 1){ if(rtn_mcs) DeleteMCS(rtn_mcs); rtn_mcs=0; break; }
	// PutSet(stderr,ChildFG); PutSet(stderr,ChildBG);
	// rtn_mcs=InsertNode(pID,ID,hsi,TRUE); // sample parent node seqs after insertion...
	rtn_mcs=InsertNode(pID,ID,hsi,FALSE); 
	hpt=rtn_mcs->GetHpt(); d=rtn_mcs->SampleColumns( );
        // if(efp) {rtn_mcs->PutMapContributions(efp,lpr); rtn_mcs->PutHyperPartition(efp); }

	//------------ Try moving child nodes between the parent (target) and insert nodes -----------
	char **HP=hpt->RtnHyperPartition(); 
	c=hpt->SetIDtoI(ID); p=hpt->SetIDtoI(pID);
	for(changed=FALSE,r=1; r < hpt->NumSets(); r++){
	   if(r==p) continue;
	   id=hpt->ItoSetID(r);
	   if(!MemberSet(id,ChildFG) && !MemberSet(id,ChildBG)) continue;
	   if(MemberSet(id,Changed)){
		if(efp && rtn_mcs->RtnContribLLR(r,c,d)){
			if(MemberSet(id,ChildFG)) fprintf(efp,"%d.Set%d: (+) %.2f\n",r,id,d);
			else fprintf(efp,"%d.Set%d: (-) %.2f\n",r,id,d);
		} continue;	// previously changed; don't modify again.
	   }
	   if(rtn_mcs->RtnContribLLR(r,c,d)){	// Calculate c-th BPPS subLLR (=d) for r-th subtree set.
		if(HP[r][c] == '+'){	// some of these are not direct children of parent!
		    if(efp) fprintf(efp,"%d.Set%d: (+) %.2f ",r,id,d);
		    if(d < 0.0){
		        if(efp) fprintf(efp," --> (-)");
			DeleteSet(id,ChildFG); AddSet(id,ChildBG); changed=TRUE; AddSet(id,Changed);
		    } if(efp) fprintf(efp,"\n");
		}
		// if(MemberSet(id,Changed)) continue;	// previously changed; don't modify again.
		if(HP[r][c] == '-'){
		    if(efp) fprintf(efp,"%d.Set%d: (-) %.2f",r,id,d);
		    if(d < 0.0){
		        if(efp) fprintf(efp," --> (+)");
			DeleteSet(id,ChildBG); AddSet(id,ChildFG); changed=TRUE; AddSet(id,Changed);
		    } if(efp) fprintf(efp,"\n");
		} assert(HP[r][c] != 'o');
	   }
	}
	if(changed){
		fprintf(stderr,"(hierarchy being modified)\n");
		assert(CardInterSet(ChildFG,ChildBG) == 0);
		ClearSet(TmpFG); ClearSet(TmpBG);
	        for(c=1; c <= Hpt->NumBPPS(); c++){ // Calculate sub-LLRs for each set.
		  	id=Hpt->ItoSetID(c);
			if(MemberSet(id,ChildFG)){ UnionSet(TmpFG,Set[c]); }
			if(MemberSet(id,ChildBG)){ UnionSet(TmpBG,Set[c]); }
		} 
       		sst_typ *xsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,D,NumCols,'M'); 
		fprintf(stderr,"       --- LPR=%.3f ---\n",D);
		hsi->Store(ChildFG,ChildBG,xsst,D,CardSet(TmpFG),CardSet(TmpBG));
     		DeleteMCS(rtn_mcs); rtn_mcs=0; 
	}
     } while(changed);
     NilSet(ChildFG); NilSet(ChildBG); NilSet(Changed); DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
     return rtn_mcs;
}

Int4	omc_typ::SampleInternalNode(Int4 pID, hsi_typ *subroot,double Temp)
#if 0	//*** STAGE 2: Sample each FG node (of insert node) for either merging or child assignment.
   Pass subroot to SampleInternalNode(); add all nodes at once.
   Start from higher up the tree and go down...
   for(s=1; s <= Hpt->NumBPPS(); s++){ if(hsi[s]) delete hsi[s]; } free(hsi);
   modify Hpt and sma, destroy mcs and create mcs_typ.
   generate internal node csq's for sma file.
/***************************************************************
     o                o(pID)
    /|\              /|\
   / | \  Insert-->   *(ID)	based on LPR calculations.
  /  |  \            /|\
 o   o   o          o o o
 ***************************************************************/
#endif	//************************************************************************
{
	Int4	x,p,id,i,j,n,R,r,c,I,pI,ID,*Parent;
	double	D,d; 
	char	*pttrn;
	set_typ	st;
	//--------------- Check input and get identifiers --------------
	assert(pID > 0 && pID < MaxNumNodes); assert(subroot != 0); 
	if(stable && MemberSet(pID,stable)) return 0; // This may not be necessary....already checked?
	for(ID=1; ID <= MaxNumNodes; ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID:
	assert(ID > 0 && ID <= MaxNumNodes);
	p = Hpt->SetIDtoI(pID); // Hpt->Put(stderr,FALSE); 
	fprintf(stderr,"================== parent %d (Set%d) ====================\n",p,pID);
	subroot->Put(stderr);

	//-----------------  Insert the new node ----------------------
	mcs_typ	*rtn_mcs=OptimizeInsertedNode(subroot); if(rtn_mcs==0) return 0;
	hpt_typ *hpt=rtn_mcs->GetHpt(); rtn_mcs->DoEvolve(); rtn_mcs->SaveBest=TRUE; 
	rtn_mcs->SampleColumns( ); 
	rtn_mcs->PutMapContributions(stderr); // rtn_mcs->PutHyperPartition(stderr); exit(1);

	//===================== Run a series of quick sampling checks =======================
	D=CalcLLR(mcs); d=rtn_mcs->CalcTotalLPR();
	if(d < 0.50*D){ DeleteMCS(rtn_mcs); return 0; }

	//=========  Require a minimum average subLLR contribution. ================
	Parent=0; assert(hpt->IsScrambledTree(Parent)); pI=hpt->SetIDtoI(pID); I=hpt->SetIDtoI(ID);
	double sum=0.0;
	for(n=0,i=2; i < hpt->NumSets(); i++){
		if(Parent[i]==I){ n++; assert(rtn_mcs->RtnContribLLR(i,I,d)); sum+=d; }
	} d=sum/(double)n; free(Parent); 
	fprintf(stderr," d = %.2f; cut = %.2f\n",d,MinimumLLR);
	if((d*1.1) < MinimumLLR){ DeleteMCS(rtn_mcs); return 0; }

	//=========  If the new is better than the old at this point then return. ================
	d=rtn_mcs->RtnBestLPR(); D = mcs->RtnBestLPR();
	if(d > D){  if(this->SampleMCS(rtn_mcs,ID)) { this->CalcLLR(mcs); return ID; } else return 0; }

	//========= Otherwise Sample over the parent and child sets. ================
	if(hpt->SetIDtoI(pID) == 1) st=0; 
	else {
	  st=MakeSet(hpt->NumSets()+1); ClearSet(st);
	  this->ReSetSubTree(st,hpt->SetIDtoI(pID),hpt);  // sample only over parent subtree.
	} 
	// Int4	modstart=rtn_mcs->TargetMod/4;
	Int4	modstart=rtn_mcs->TargetMod/8;
	modstart=MAXIMUM(Int4,modstart,rtn_mcs->ModStart);
	if(st) fprintf(stderr,"Card set 'st' = %d\n",CardSet(st));
	rtn_mcs->SampledSet=st; rtn_mcs->Sample(2,2,2,2,modstart); if(st) NilSet(st);
	rtn_mcs->SampledSet=0; rtn_mcs->RestoreBest(); d=rtn_mcs->CalcTotalLPR( );
	if(d > D) fprintf(stderr," Increase in LLR from %.2f to %.2f\n",D,d); 
	else fprintf(stderr," Decrease in LLR from %.2f to %.2f\n",D,d);
// 	rtn_mcs->PutHyperPartition(stderr);
	if(this->SampleMCS(rtn_mcs,ID)){ this->CalcLLR(mcs); return ID; } else return 0;
}

void    omc_typ::ReSetSubTree(set_typ subtree, Int4 c, hpt_typ *hpt)
{
        ClearSet(subtree);
        for(Int4 row=1; row < hpt->NumSets(); row++)
           { if(hpt->Cell(row,c) == '+') AddSet(row,subtree); }
}

void    omc_typ::ReSetIDSubTree(set_typ subtree, Int4 c)
{
        ClearSet(subtree);
        for(Int4 row=1; row < Hpt->NumSets(); row++){
            if(Hpt->Cell(row,c) == '+') AddSet(Hpt->ItoSetID(row),subtree);
        }
}

Int4	omc_typ::ForcedInternalNode(Int4 pID, hsi_typ *subroot)
#if 0	//*** STAGE 2: Sample each FG node (of insert node) for either merging or child assignment.
   Pass subroot to AddInternalNode(); add all nodes at once.
   Start from higher up the tree and go down...
   for(s=1; s <= Hpt->NumBPPS(); s++){ if(hsi[s]) delete hsi[s]; } free(hsi);
   modify Hpt and sma, destroy mcs and create mcs_typ.
   generate internal node csq's for sma file.
/***************************************************************
     o                o(pID)
    /|\              /|\
   / | \  Insert-->   *(ID)	based on LPR calculations.
  /  |  \            /|\
 o   o   o          o o o
 ***************************************************************/
#endif	//************************************************************************
{
	Int4	p,I,pI,ID;
	//--------------- Check input and get identifiers --------------
	assert(pID > 0 && pID < MaxNumNodes); assert(subroot != 0); 
	if(stable && MemberSet(pID,stable)) return 0; // This may not be necessary....already checked?
	for(ID=1; ID <= MaxNumNodes; ID++){ if(Hpt->SetIDtoI(ID) == 0) break; } // find an unused ID:
	assert(ID > 0 && ID <= MaxNumNodes);
	p = Hpt->SetIDtoI(pID); // Hpt->Put(stderr,FALSE); 
	fprintf(stderr,"================== parent %d (Set%d) ====================\n",p,pID);
	subroot->Put(stderr);
	//-----------------  Insert the new node ----------------------
	mcs_typ	*rtn_mcs=OptimizeInsertedNode(subroot); if(rtn_mcs==0) return 0;
	hpt_typ *hpt=rtn_mcs->GetHpt(); rtn_mcs->DoEvolve(); rtn_mcs->SaveBest=TRUE; 
	rtn_mcs->SampleColumns( ); 
	rtn_mcs->PutMapContributions(stderr); // rtn_mcs->PutHyperPartition(stderr); // exit(1);
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt();
	return 1;
}

Int4	omc_typ::InsertInternalNodes(FILE *fp, double Temp)
// Find related nodes that need a common (inserted) ancestral node.
{
    Int4  f,b,i,j,s,t,v,x,RtnHits=0,SetSize=mcs->GetSetSize();
    Int4  fID,bID,cID,pID,parent,*Parent=0,NumCols=mcs->RtnDefaultMaxCol();
    mcs->RestoreBest();	// WARNING: also need to UnRestoreBest()!!??
    // mcs->PutHyperPartition(); 
    double d,dd,D,DD;
    Int4  Start=1,End=Hpt->NumBPPS();

    if(NumForced > 0){ Start=ForcedP;  End=Start; }
    for(parent = Start; parent <= End; parent++) {
	if(Hpt->TypeOfSet(parent) != '?') continue;
	pID=Hpt->ItoSetID(parent); 
	if(stable && MemberSet(pID,stable)) continue;
	Parent=0; assert(Hpt->IsScrambledTree(Parent));

	ClearSet(ChildSetBG); ClearSet(ChildSetFG); ClearSet(ChildSetU); ClearSet(ChildSetBoth);
	ClearSet(SetFG); ClearSet(SetBG); 
	//============== Get subtree node sets. ====================
	set_typ	*subtree=0; NEW(subtree,Hpt->NumBPPS()+2,set_typ); 
	for(s=1; s <= Hpt->NumBPPS(); s++){ 	
		subtree[s]=MakeSet(MaxNumNodesPlus+2); // get set of sampled node subtree.
		// ReSetIDSubTree(subtree[s],s); // fprintf(stderr,"%d: ",s); PutSet(stderr,subtree[s]);
		ReSetSubTree(subtree[s],s); // fprintf(stderr,"%d: ",s); PutSet(stderr,subtree[s]);
	}
	//============== Get sequence subsets for each node. ====================
	NEW(Set,Hpt->NumSets() +3, set_typ);
	set_typ *set=mcs->CopyOfSeqSets();
	for(s=1; s <= Hpt->NumBPPS(); s++){
	    Set[s]=MakeSet(SetSize);  ClearSet(Set[s]); 
	    for(t=1; t <= Hpt->NumBPPS(); t++){ // Put subtree seqs into Set[s].
	        if(MemberSet(t,subtree[s])){ UnionSet(Set[s],set[t]); } 
	    } // fprintf(stderr,"CardSet(set[%d]) = %d (%d)\n",s,CardSet(set[s]),CardSet(Set[s])); 
	} for(s=1; s <= Hpt->NumBPPS(); s++){ NilSet(subtree[s]); } free(subtree);

	fprintf(stderr,"***************** parent %d **************************\n",parent);
#if 0	// Turned off: put parent in neither the FG nor BG, as it may be split between these.
	AddSet(pID,ChildSetBG); AddSet(pID,ChildSetBoth); AddSet(pID,ChildSetU);
	UnionSet(SetBG,set[parent]);	// Put parent node into BG sequence set.
#endif
	for(s=1; s <= Hpt->NumBPPS(); s++){	// Put children in Universal set.
	    if(Parent[s] == parent){ cID=Hpt->ItoSetID(s); AddSet(cID,ChildSetU); }
	} Int4 NumChild=CardSet(ChildSetU);

        hpt_typ *oHpt=Hpt->Copy();	// Save a copy of original Hpt.
	hsi_typ *hsi,*best_hsi=0; 

	if(CardSet(ChildSetU) > 3){	// parent + 2FG + 1BG = 4.
	  sst_typ	*xsst=0;
	  double	Lpr=0.0;
	  char		NdTyp;

//************** STEP 1: Search 3-set contrast alignments for candidate seed sets. *****************
	  if(NumForced > 0){
		for(f=1; f <= NumForced; f++){	// put subtree seqs for f into SetFG.
		    x=ForcedNode[f];
		    cID=Hpt->ItoSetID(x); AddSet(cID,ChildSetFG); UnionSet(SetFG,Set[x]);
	        }
		for(b=2; b <= Hpt->NumBPPS(); b++){	// b == background node.
		    bID=Hpt->ItoSetID(b); 
		    if(Parent[b] != parent || MemberSet(bID,ChildSetFG)) continue;
		    AddSet(bID,ChildSetBG); UnionSet(SetBG,Set[b]);  // put in the BG.
		}
        	xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M'); // free(xsst);
	        hsi= new hsi_typ(ChildSetFG,ChildSetBG,xsst,Lpr,CardSet(SetFG),CardSet(SetBG),lpr);
	        PrintPttrn(stderr,hsi->GetSST()); hsi->Put(stderr); fprintf(stderr,"\n");
		ForcedInternalNode(pID, hsi); delete hsi;
		// this->SampleMCS(mcs,ID);
		for(s=1; s<=oHpt->NumSets(); s++){ 
			if(set[s]) NilSet(set[s]); if(Set[s]) NilSet(Set[s]);
		} free(Set); Set=0; free(set); delete oHpt; free(Parent);
		return 1;
	  }
	  hmh_typ hmh(Hpt->NumBPPS());	// hsi mheap.
	  for(f=2; f <= Hpt->NumBPPS(); f++){	// don't use root node! f = foreground node.
	        hsi_typ *bst_hsi=0; 
		if(Parent[f] != parent) continue;	// must be a child of parent...
		cID=Hpt->ItoSetID(f); AddSet(cID,ChildSetFG);		
		UnionSet(SetFG,Set[f]);		// put subtree seqs for f into SetFG.
		for(b=2; b <= Hpt->NumBPPS(); b++){	// b == background node.
		    if(Parent[b] != parent || b == f) continue;
		    bID=Hpt->ItoSetID(b); AddSet(bID,ChildSetBG);
		    UnionSet(SetBG,Set[b]);		// put one Set in the BG.
        	    xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M'); free(xsst);
		    if(Lpr > 0.0){
	   	      for(Int4 f2=f+1; f2 <= Hpt->NumBPPS(); f2++){
	        	if(Parent[f2] != parent || f2 == b) continue;
	      		Int4 f2ID=Hpt->ItoSetID(f2); AddSet(f2ID,ChildSetFG);		
	        	UnionSet(SetFG,Set[f2]);		// put Set in the FG.
			// Insert component root node LLR == Lpr.
        	        xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M'); free(xsst);
		        if(Lpr > 0.0){
	    	          UnionSet3(ChildSetFG,ChildSetBG,ChildSetBoth);
			  hsi=CheckForBottom(3,3,0);  // return a three node hybrid hierarchy.
		          if(hsi){	
		            Lpr=hsi->Lpr();
			    if(preHG) IncdHist(Lpr,preHG);
		            if(bst_hsi == 0 || bst_hsi->Lpr() < hsi->Lpr()){
			      if(bst_hsi) delete bst_hsi; bst_hsi = hsi;
		            } else delete hsi;  hsi=0;
		          } 
			} IntersectNotSet(SetFG,Set[f2]); DeleteSet(f2ID,ChildSetFG);
		      }
	   	    } IntersectNotSet(SetBG,Set[b]);  DeleteSet(bID,ChildSetBG);
		} IntersectNotSet(SetFG,Set[f]); DeleteSet(cID,ChildSetFG);
		if(bst_hsi){   // put the best one for each node on the heap; will discard duplicates.
		    if(hmh.Insert(bst_hsi,bst_hsi->Lpr()) == 0){ delete bst_hsi; } bst_hsi=0;  
		}
	  }
//************** STEP 2: Do a full dfs search using seed sets. *****************
	  double key; Int4 Item;
	  hsi_typ *tmp_hsi=0;
	  hmh_typ hmhB(10);	// smaller min-max-heap of full partitions...
	  for(i=1; (hsi=hmh.DelMax(&key,&Item)) != 0; i++){
#if 0
	     fprintf(stderr," ++++++++++++ seed %d: %.2f nats +++++++++\n",Item,key);
	     PrintPttrn(stderr,hsi->GetSST()); hsi->Put(stderr); fprintf(stderr,"\n");
#endif
	     if(hmhB.IsSuperInHeap(hsi)){ delete hsi; continue; } // optimal config for hsi already in heap..
	     //----------- Reinstate sets as found for quick search... --------------
	     CopySet(ChildSetFG,hsi->RtnSetFG()); CopySet(ChildSetBG,hsi->RtnSetBG());
	     ClearSet(SetFG); ClearSet(SetBG); // CopySet(SetBG,set[parent]); // put parent in BG.
	     for(s=1; s <= Hpt->NumBPPS(); s++){
		x=Hpt->ItoSetID(s);
		if(MemberSet(x,ChildSetFG)) UnionSet(SetFG,Set[s]);
		if(MemberSet(x,ChildSetBG)) UnionSet(SetBG,Set[s]);
	     } UnionSet3(ChildSetFG,ChildSetBG,ChildSetBoth);
	     //--------------- Look all the way down ---------------
	     tmp_hsi=OnTheFlyLPR_dfs(2,INT_MAX,parent); delete hsi; hsi=tmp_hsi; tmp_hsi=0;
	     if(hsi){
		Lpr=hsi->Lpr();
		if(hmhB.Insert(hsi,Lpr) == 0){ 
fprintf(stderr," ========================= duplicate (%d) %.2f nats =====================\n\n",i,Lpr);
			delete hsi; // Will reject if hsi is a duplicate.
		} else {
fprintf(stderr," =========================== full (%d) %.2f nats =======================\n",i,Lpr);
PrintPttrn(stderr,hsi->GetSST()); hsi->Put(stderr); fprintf(stderr,"\n");
		} hsi=0;  
	     } 

	  }
//************** STEP 3: Sample using full search results. *****************
	  for(i=1; (hsi=hmhB.DelMax(&key,&Item)) != 0; i++){
	    fprintf(stderr,"%d: item = %d; %.2f nats; p = %d; pID = %d\n",
				i,Item,key,parent,Hpt->ItoSetID(parent));
            cID=this->SampleInternalNode(Hpt->ItoSetID(parent),hsi,Temp);
	    if(cID != 0){ 	//------> then Hpt and mcs are now new!!!
		fprintf(stderr," --------- cID = %d -------------\n",cID);
		// this->RestoreVeryBest(); RtnHits++; 
		mcs->RestoreBest(); RtnHits++; 
		// mcs->SaveBest=TRUE; mcs->StoreBest();
		this->CalcLLR(mcs); this->RestoreVeryBest(); 
		do { delete hsi; } while((hsi=hmhB.DelMax(&key,&Item)) != 0);  
		parent = 0; break;	//!!!!!!!!! start all over again... !!!!!!!!!!
	    } else { fprintf(stderr,"InsertNode() candidate failed.\n"); delete hsi; } hsi=0;
	  } fprintf(stderr," --------- Number of hits = %d -------------\n",RtnHits);
	} //--------  hmh & hmhB deleted when leaving this "if(CardSet(ChildSetU) > 3)" scope ------
	for(s=1; s<=oHpt->NumSets(); s++){ if(set[s]) NilSet(set[s]); if(Set[s]) NilSet(Set[s]);}
	free(Set); Set=0; free(set); delete oHpt; free(Parent);
    }
    // NilSet(ChildSetBG); NilSet(ChildSetFG); NilSet(ChildSetU); NilSet(ChildSetBoth);
    return RtnHits;
}

hsi_typ	*omc_typ::OnTheFlyLPR_dfs(Int4 depth, Int4 maxdepth, Int4 parent)
{
        double	Lpr,D,d;
    	Int4	NumCols=mcs->RtnDefaultMaxCol();
	Int4	i,s,t,c,cID,stage,N_neg=0;
	sst_typ	*xsst=0,*tsst=0;
	char	NdTyp;
	
	hsi_typ	*hsi=CheckForBottom(depth,maxdepth,parent);
	if(hsi) return hsi;
	dh_type dH=0;
	for(stage = 0; stage <= 1; stage++){
	  if(stage == 1 && CardSet(ChildSetFG) < 2) break;	// never use just one FG node...
	  for(s=2; s<=Hpt->NumBPPS(); s++){
	    cID=Hpt->ItoSetID(s);
	    if(!MemberSet(cID,ChildSetU)) continue;	// look at candidate child nodes only
	    if(MemberSet(cID,ChildSetFG) || MemberSet(cID,ChildSetBG)) continue;
	    //=============== Add another child node to FG or BG. =================
	    if(stage == 0) { AddSet(cID,ChildSetFG); UnionSet(SetFG,Set[s]); }	// into FG...
	    else { AddSet(cID,ChildSetBG); UnionSet(SetBG,Set[s]); }		// into BG..
	    UnionSet3(ChildSetFG,ChildSetBG,ChildSetBoth);
	    //=============== Compute the optimal pattern & the LLR. =================
            xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M');	// Get Pttrn & LPR.
	    if(Lpr > 0.0){
#if 1		//--------------- Calculate sub-LLRs for each set. ---------------
		for(N_neg=0,c=2; c <= Hpt->NumBPPS(); c++){   // ensure that each set matches pattern.
		    if(!MemberSet(Hpt->ItoSetID(c),ChildSetFG)) continue;
		    if(Hpt->TypeOfSet(c) == '!') NdTyp='L'; else NdTyp='M';
		    if(CalcSetvsPttrnLPR(0,Set[c],SetBG,xsst,NdTyp) < 0.0) N_neg++; 
		}
		if(N_neg == 0) 
#endif
		{
	          if(CardSet(ChildSetBoth) == CardSet(ChildSetU)){   // == hit bottom of search tree.
		    hsi_typ *rtn_hsi=OnTheFlyLPR_dfs(depth,maxdepth,parent);  // WILL return hsi_typ != 0.
		    if(hsi){
		       if(hsi->Lpr() < rtn_hsi->Lpr()){	delete hsi; hsi=rtn_hsi; } else delete rtn_hsi;
		    } else hsi=rtn_hsi;
	          } else {	//---- otherwise store node id on the heap.
			if(dH==0) dH=dheap(2*Hpt->NumBPPS()+1,3);
			if(stage == 1) insrtHeap(s+Hpt->NumBPPS(),-Lpr,dH); else insrtHeap(s,-Lpr,dH);
		  }
	        } 
	    } if(xsst){ free(xsst); xsst=0; }
	    if(stage == 0) { DeleteSet(cID,ChildSetFG); IntersectNotSet(SetFG,Set[s]); }
	    else { DeleteSet(cID,ChildSetBG); IntersectNotSet(SetBG,Set[s]); }
	  }
	}  //============ search deeper with the best configurations found at this level ==========
	if(dH==0) return hsi; 	// Bottom reached: dH should be null and hsi not null.
	else {			// otherwise hsi should be null.
	  hsi_typ *xhsi=0; assert(hsi==0);  // 
	  for(i=0; i < 1 && !emptyHeap(dH); i++){	// i < 1 --> look at best configuration only!!
	    Lpr= (double) -minkeyHeap(dH);
	    assert((s=delminHeap(dH)) != 0);
	    if(s > Hpt->NumBPPS()){ s=s-Hpt->NumBPPS(); stage=1; } else { stage = 0; }
	    cID=Hpt->ItoSetID(s);
	    if(stage == 0) { AddSet(cID,ChildSetFG); UnionSet(SetFG,Set[s]); }
	    else { AddSet(cID,ChildSetBG); UnionSet(SetBG,Set[s]); }
	    UnionSet3(ChildSetFG,ChildSetBG,ChildSetBoth);
	    //---------- search deeper starting with the best configuration. ------------
	    xhsi=OnTheFlyLPR_dfs(depth+1,maxdepth,parent);  // will return hsi_typ !
	    if(hsi == 0){ hsi=xhsi; }
	    else if(xhsi){
	    	if(xhsi->Lpr() > hsi->Lpr()){ delete hsi; hsi=xhsi; } else delete xhsi; xhsi=0;
	    } if(stage == 0) { DeleteSet(cID,ChildSetFG); IntersectNotSet(SetFG,Set[s]); }
	    else { DeleteSet(cID,ChildSetBG); IntersectNotSet(SetBG,Set[s]); }
	  } Nildheap(dH);
	} return hsi;
}

sst_typ	*omc_typ::ComputeHybridLLR(Int4 parent, double &Lpr)
// Compute the difference between the current state and the transition state LLRs
{
        double	d,dd,D,DD;
	char	NdTyp;
	Int4	t,NumCols=mcs->RtnDefaultMaxCol();
	sst_typ	*tsst,*xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M');
	if(Lpr <= 0.0) return xsst;
fprintf(stderr,"!! Insert component subLLR = %.3f\n",Lpr);
	for(dd=DD=0,t=2; t<=Hpt->NumBPPS(); t++){
	        if(!MemberSet(Hpt->ItoSetID(t),ChildSetFG)) continue;
		if(Hpt->TypeOfSet(t) == '!') NdTyp='L'; else NdTyp='M';

		IntersectNotSet(SetFG,Set[t],TmpBG); // TmpBG = SetFG intersect not Set[t].
		tsst=this->GetOptPttrnLPR(0,Set[t],TmpBG,d,NumCols,NdTyp); free(tsst);
		dd += d;

		IntersectNotSet(Set[parent],Set[t],TmpBG); // TmpBG = Set[parent] intersect not Set[t].
	        tsst=this->GetOptPttrnLPR(0,Set[t],TmpBG,D,NumCols,NdTyp); free(tsst);
		DD += D;
fprintf(stderr,"!! node %d: child insert = %.3f vs main %.3f; diff = %.3f\n",t,d,D,d-D);
	} Lpr = Lpr + dd - DD;
fprintf(stderr,"!! Net gain in LLR = %.3f\n",Lpr);
	return xsst;
}

hsi_typ	*omc_typ::CheckForBottom(Int4 depth, Int4 maxdepth, Int4 parent)
{
        double	Lpr;
	sst_typ	*xsst=0;
	hsi_typ	*hsi=0;
	
	if(depth >= maxdepth || (CardSet(ChildSetBoth) == CardSet(ChildSetU) 
		&& CardInterSet(ChildSetBoth,ChildSetU) == CardSet(ChildSetU))){
#if 1	//============================ transition state ================================
	   if(parent > 0) xsst=this->ComputeHybridLLR(parent,Lpr);
	   else xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,mcs->RtnDefaultMaxCol(),'M');
#else
	   Int4	NumCols=mcs->RtnDefaultMaxCol();
           xsst=this->GetOptPttrnLPR(0,SetFG,SetBG,Lpr,NumCols,'M');
#endif	//==============================================================================
	   if(dfsHG) IncdHist(Lpr,dfsHG);
	   hsi= new hsi_typ(ChildSetFG,ChildSetBG,xsst,Lpr,CardSet(SetFG),CardSet(SetBG),lpr);
	} return hsi;
}

