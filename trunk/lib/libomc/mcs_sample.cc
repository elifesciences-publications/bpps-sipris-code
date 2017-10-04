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

#include "mcs_typ.h"
#include "blosum62.h"

void	mcs_typ::UpdateDisplaySeqs( )
{
	Int4	n,s,len;
	for(n=1; n<= Hpt->NumBPPS(); n++){
	  assert(NumDisplayCMA == Hpt->NumBPPS());
	  e_type keyE=che[n]->KeyE( ); 
          len=LengthCMSA(1,DisplayCMA[n]);
	  e_type tE=TrueSeqCMSA(1,DisplayCMA[n]);
          assert(len == LenSeq(keyE)); assert(len == LenSeq(tE));
          for(s=1; s <=len; s++){
            unsigned char r = ResSeq(s,keyE);
            if(r != UndefAlpha(AB)){
		SetResidueCMSA(1,1,s,r,DisplayCMA[n]); // Fake seq...
		EqSeq(s,r,tE);		// True seq...
	    }
          } // PutCMSA(stderr,DisplayCMA[n]);
	}
}

void    mcs_typ::ReSetRelations(FILE *fp)
// Recompute FG and BG sets and relationships.
{
	Int4	x;
	for(x=1; x <= Hpt->NumBPPS(); x++) NilSet(SetBG[x]);  free(SetBG);
	for(x=1; x <= Hpt->NumBPPS(); x++) NilSet(SetFG[x]);  free(SetFG);
	for(x=1; x <= Hpt->NumBPPS(); x++) free(RelateFGs[x]);  free(RelateFGs);
	for(x=1; x <= Hpt->NumBPPS(); x++) free(RelateBGs[x]);  free(RelateBGs);
	RelateFGs=GetSetRelations("FG",Hpt->nGrpsFG(),Hpt->GrpsFG(),&SetFG);
       	RelateBGs=GetSetRelations("BG",Hpt->nGrpsBG(),Hpt->GrpsBG(),&SetBG);
	if(fp){
    	  fprintf(fp,"FG set:\n");
	  for(x=1; x <= Hpt->NumBPPS(); x++){ fprintf(fp,"%d. ",x); PutSet(fp,SetFG[x]); }
	  PutSetRelations(fp);
    	  fprintf(fp,"BG set:\n");
	  for(x=1; x <= Hpt->NumBPPS(); x++){ fprintf(fp,"%d. ",x); PutSet(fp,SetBG[x]); }
	  fflush(fp);
	}
}

BooLean	mcs_typ::SampleHpt(FILE *fp,Int4 SampledCol)
// WARNING: this routine assumes direct correspondence between column and rows!!
// So that the diagonal equates rows (sets) with columns (FD-categories)!!
/**********************************************************************
Sample Set5 into Set15?

HyperParTition:
!!!!!!!!!!!!!!!
+--oo--ooo----- 1.Set1?
++-oo--ooo----- 2.Set12!
+-+----ooo----- 3.Set11?
+-++---ooo----- 4.Set10!
+-+-+--ooo----- 5.Set9!
+--oo+-ooo----- 6.Set8!
+--oo-+-------- 7.Set15?        +--oo-+-------- 7.Set15?  (potential parent)
+--oo-++------- 8.Set16!
+--oo-+-+------ 9.Set14!
+--oo-+--+----- 10.Set7!
+--oo--ooo+---- 11.Set6!
+--oo--ooo-+--- 12.Set5!        +--oo--ooo-+--- 12.Set5!  (original)
                                   :        :       :
                                +--oo-+----+--- 12.Set5!  (new)
+--oo--ooo--+-- 13.Set4! 	(needs to be the same as parent except for diagonal).
+--oo--ooo---+- 14.Set3!
+--oo--ooo----+ 15.Set2!
-oooooooooooooo 16.Random=804.

/**********************************************************************/
// 3. Change that cell and recompute LPR.
// 4. If this improves the LPR then accept the change, else revert back... 
//************ Testing this for cdhBPPS routine: adding an internal node...
// WARNING: assumes a direct correspondence between columns and rows!!
{

	char	state,*State;
	double	olpr,nlpr,d;
	Int4	row_neg,row_pos,row_omit,col_neg,col_pos,col_omit;
	Int4	g,n,r,c,x,sq,*Parent=0,num_sq=NumSeqsCMSA(MainCMA);
	BooLean	IsChanged=FALSE;

	assert(Hpt->NumSets() == Hpt->NumBPPS() +1); // assumes direct correspondence between column and rows!!
	n = SampledCol;
	if(n <= 1) return FALSE;	// don't look at root node.
	if(Hpt->TypeOfSet(n) != '?') return FALSE; // Must be a misc. column to change within the hpt.
	if(Hpt->IsTree(Parent) == FALSE){ if(Parent) free(Parent); return FALSE; }
	char    **HP=Hpt->RtnHyperPartition();
        for(g=2; g < Hpt->NumSets(); g++){	// for each row... Don't look at Root or Rejects...
	    if(Hpt->Cell(g,n) == '+') continue;	// don't look at current child sets right now...
	    if(Parent[g] != Parent[n]) continue;  // look only at siblings of n.
	    state=Hpt->Cell(g,n);	// state == 'o' or '-'.
	    assert(state == '-');	// must be the case for siblings in a tree.
	    // State=AllocString(HP[g]);	// Save old HP states.
	    if(Hpt->TypeOfSet(g) == '?') ; // then need to bring along all of its child nodes as well.
#if 0
	    if(Hpt->Cell(g,n) == 'o') continue;	// don't look at omitted cells right now...
	    assert(Hpt->Cell(g,n) == '-'); // 
	    Hpt->RowsInColumn(n,row_neg,row_pos,row_omit);   // count # '-' cells in column n.
	    if(row_neg < 2) continue; 	// Need at least one background set? (how about parent node)
#endif
	    olpr = CalcTotalLPR(0,FALSE);
	    for(sq=1; sq <= num_sq; sq++){
		    if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); }
	    } nlpr=CalcTotalLPR(0,FALSE);
	    d=nlpr-olpr;
	    if(d < 10.0){	// then revert back to previous..
	      if(fp) fprintf(fp,"Failed to change HP[%d][%d] from '%c' to '+' (delta=%.3f).\n",
						g,n,state,d);
		for(sq=1; sq <= num_sq; sq++){
		     if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition(state,sq); } }
	    } else {	// accept change...
	      BooLean OkayToAdd=TRUE;
	      if(OkayToAdd){ 
	        if(fp) fprintf(fp,"Changed HP[%d][%d] from '%c' to '+' (delta=%.3f).\n",
						g,n,state,d);
	        Hpt->Change(state,'+',g,n); ReSetRelations(); assert(HP[g][n] == '+');
		IsChanged=TRUE;
	      } else { // revert back...
		for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('-',sq); } }
	      }
	    }
	} CalcTotalLPR( );
	if(Parent) free(Parent);
	return IsChanged;
}

BooLean	mcs_typ::SampleHpt(FILE *fp,Int4 Root, wdg_typ &Tree)
// WARNING: this routine assumes direct correspondence between column and rows!!
// So that the diagonal equates rows (sets) with columns (FD-categories)!!
{
    // 3. Change that cell and recompute LPR.
    // 4. If this improves the LPR then accept the change, else revert back... 

    Int4	g,n,r,c,x;
    char	state;
    double	olpr,nlpr;
    Int4	row_neg,row_pos,row_omit,col_neg,col_pos,col_omit;
    Int4	sq,num_sq = NumSeqsCMSA(MainCMA);

    wdg_typ NewTree=0;	// NewTree to merge with old tree...
    if(Tree){ NewTree=MkWdgraph(WdgraphN(Tree),WdgraphM(Tree)); }

    assert(Hpt->NumSets() == Hpt->NumBPPS() +1); // assumes direct correspondence between column and rows!!
    // Start at n=2 because don't want to modify main set...
    for(n=2; n<= Hpt->NumBPPS(); n++){	// 2. Find a column to work on.
	g=n; 	// assume direct correspondence between column and rows!!
	if(Hpt->TypeOfSet(g) != '?') continue; // find a misc. column 
    	// 1. Found miscellaneous column to change within the hpt.
	char    **HP=Hpt->RtnHyperPartition();
        for(g=1; g < Hpt->NumSets(); g++){		// for each row... Don't look at Rejects...
	    if(Hpt->Cell(g,n) == '+') continue;	// don't look at positive cells right now...
	    if(Hpt->Cell(g,n) == 'o') continue;	// don't look at omitted cells right now...
	    assert(Hpt->Cell(g,n) == '-');
	    Hpt->RowsInColumn(n,row_neg,row_pos,row_omit);   // count # '-' cells in column n.
	    if(row_neg < 2) continue; 
	    olpr = CalcTotalLPR(0,FALSE);
	    for(sq=1; sq <= num_sq; sq++){
		    if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); }
	    } nlpr=CalcTotalLPR(0,FALSE);
	    double d=nlpr-olpr;
	    if(d < 10.0){	// then revert back to previous..
	      if(fp) fprintf(fp,"Failed to change HP[%d][%d] from '-' to '+' (delta=%.3f).\n",g,n,d);
		for(sq=1; sq <= num_sq; sq++){
		     if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('-',sq); } }
	    } else {	// accept change...
	      BooLean OkayToAdd=TRUE;
	      if(Tree){		// Then add an edge from set in row g to the misc node n within the tree.
		// WARNING: This assumes that Hpt was first generated from a tree! Won't work otherwise.
		if(Hpt->AddEdgeToTree(n,g,Root,(Int4)d,Tree,NewTree) == 0) OkayToAdd=FALSE;
		// PrintNewickTreeWDG(stderr,Root, Tree); // need to pass in Root
	      } // Add edges to the input Tree.
	      if(OkayToAdd){ 
	        if(fp) fprintf(fp,"Changed HP[%d][%d] from '-' to '+' (delta=%.3f).\n",g,n,d);
	        Hpt->Change('-','+',g,n); ReSetRelations(); assert(HP[g][n] == '+');
	      } else { // revert back...
		for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('-',sq); } }
	      }
	      // Change the omitted columns for this set which is now included in misc set...
	      // 1. Find the rows to be omitted for this set in column g == row n set only.
	      for(olpr=nlpr,r=1; r < Hpt->NumSets(); r++){	// skip Rejected set..
		if(r == g) continue;
		if(HP[r][n] == '+'){	// if set 'r' also is in the FG for column n...
		    if(HP[g][r] == 'o'){   // & if set g is omitted from the BG for set r.
		      // then see whether this should be included in the BG.
		      for(sq=1; sq <= num_sq; sq++){
			if(MemberSet(sq,GrpSet[g])){ che[r]->SetPartition('-',sq); }
		      } nlpr=CalcTotalLPR(0,FALSE); d=nlpr-olpr;
		      if(d < 10.0){  // Then revert..
	      		if(fp) fprintf(fp,"Failed to change HP[%d][%d] from 'o' to '-' (delta=%.3f)\n",g,r,d);
			for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[g])){ che[r]->SetPartition('o',sq); } }
		      } else {
	      		if(fp) fprintf(fp,"Changed HP[%d][%d] from 'o' to '-' (delta=%.3f)\n",g,r,d);
	      		Hpt->Change('o','-',g,r); ReSetRelations(); assert(HP[g][r] == '-');
			olpr=nlpr;
		      }
		    }
		}
	      }
	      // Change '-' rows to 'o' in column g (== row g) if those rows == '-' for column n.
	      // 1. Find the rows to be changed to 'o' in column g.
	      for(r=1; r < Hpt->NumSets(); r++){	// skip Rejected set..
		if(r == g) continue;
		if(HP[r][n] == '-' && HP[r][g] == '-'){	
		      for(sq=1; sq <= num_sq; sq++){
			if(MemberSet(sq,GrpSet[r])){ che[g]->SetPartition('-',sq); }
		      } nlpr=CalcTotalLPR(0,FALSE); d=nlpr-olpr;
		      if(FALSE && d < 0.0){  // Then revert..
	      		if(fp) fprintf(fp,"Failed to change HP[%d][%d] from '-' to 'o' (delta=%.3f)\n",g,r,d);
			for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[g]->SetPartition('o',sq); } }
		      } else {
	      		if(fp) fprintf(fp,"Changed HP[%d][%d] from '-' to 'o' (delta=%.3f)\n",g,r,d);
	      		Hpt->Change('-','o',r,g); ReSetRelations(); assert(HP[r][g] == 'o');
			olpr=nlpr;
		      }
		}
	      }
	      // Change 'o' rows to '-' in column n (== row g) if those rows == 'o' for column g.
	      for(r=1; r < Hpt->NumSets(); r++){	// skip Rejected set..
		if(r == g) continue;
		if(HP[r][n] == '+' && HP[g][r] == 'o'){	
		      for(sq=1; sq <= num_sq; sq++){
			if(MemberSet(sq,GrpSet[g])){ che[r]->SetPartition('-',sq); }
		      } nlpr=CalcTotalLPR(0,FALSE); d=nlpr-olpr;
		      if(FALSE && d < 0.0){  // Then revert..
	      		if(fp) fprintf(fp,"Failed to change HP[%d][%d] from 'o' to '-' (delta=%.3f)\n",g,r,d);
			for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[g])){ che[r]->SetPartition('o',sq); } }
		      } else {
	      		if(fp) fprintf(fp,"Changed HP[%d][%d] from 'o' to '-' (delta=%.3f)\n",g,r,d);
	      		Hpt->Change('o','-',g,r); ReSetRelations(); assert(HP[g][r] == '-');
			olpr=nlpr;
		      }
		}
	      }
	    }
	}
    } 
    if(Tree){ Hpt->MergeOldTreeIntoNew(Root,Tree,NewTree); NilWdgraph(Tree); Tree=NewTree; }
    CalcTotalLPR( );
}

Int4	mcs_typ::RemoveSimilarSets( )
// Find sets that share a similar pattern and lpr with each other's patterns.
// remove the weaker of these and move sequences into the other set.
{
	Int4	  sq,g,h,i,j,k,n,lenI,lenJ,NumRemoved=0;
	bpps_typ  *ppsI,*ppsJ;
	sst_typ	  *sstI,*sstJ,*qst;
	double	  lpr,d,x,xII,xIJ,xJJ,xJI;

        Int4	num_sq = NumSeqsCMSA(MainCMA);
        char	**HP=Hpt->RtnHyperPartition();
	for(i=1; i < Hpt->NumSets(); i++){	// skip Rejected set..
	   if(Hpt->TypeOfSet(i) == '?') continue;
	   if(IsFailedSet[i]) continue;
	   ppsI=che[i]->BPPS();
	   lenI=ppsI->LenPattern( );
	   sstI=ppsI->RtnCopySST();
	   for(j=i+1; j < Hpt->NumSets(); j++){	// skip Rejected set..
		// if(j==i) continue; // shouldn't matter if loop from i+1...
	   	if(Hpt->TypeOfSet(j) == '?') continue;
	   	if(IsFailedSet[i]) continue;
		ppsJ=che[j]->BPPS();
	   	lenJ=ppsJ->LenPattern( );
		assert(lenI == lenJ);
	   	sstJ=ppsJ->RtnCopySST();
		lpr=CalcTotalLPR(0,FALSE); xII =Map[i]; xJJ=Map[j];

	// WARNING: this bypasses some restrictions for che_typ, such as min number column requirement...
		ppsI->ReplacePattern(sstJ); ppsJ->ReplacePattern(sstI);
		x=CalcTotalLPR(0,FALSE); xIJ=Map[i]; xJI=Map[j];
		ppsI->ReplacePattern(sstI); ppsJ->ReplacePattern(sstJ); // put pattern back...
		d=CalcTotalLPR(0,FALSE);
		if(x > 0.0 && xIJ > 0.0 && xJI > 0.0){
		  double fxJI=xJI/xJJ;
		  double fxIJ=xIJ/xII; 
		  char c=' ';
		  if(fxJI >= 0.50 && fxIJ >= 0.50) c='*';
		  fprintf(stderr,"%d vs %d swapped lpr = %.2f; lpr =%.2f (%.2f)\n",i,j,x,lpr,d);
		  fprintf(stderr,"   Map[%d] ptrn %d = %.3f (%.3f)%c\n",i,j,xIJ,xII,c);
		  fprintf(stderr,"   Map[%d] ptrn %d = %.3f (%.3f)%c\n",j,i,xJI,xJJ,c);
		  if(c=='*'){		// remove the weaker set.
		    if(xII < xJJ){ IsFailedSet[i]=TRUE; IsFailedBPPS[i]=TRUE;  g=i; h=j; }
		    else { IsFailedSet[j]=TRUE; IsFailedBPPS[j]=TRUE;  g=j; h=i; }
		    TransferAllSeqs(g,h); NumRemoved++; // move sequences from failed set g to h;
		    assert(g <= Hpt->NumBPPS()); n = g;
		    che[n]->SetMinNumColumns(0);	// this is required to remove all columns!!
    		    Int4 Length=LenSeq(che[n]->KeyE());
		    for(k=1; k <= Length; k++){ if(che[n]->PttrnPos(k)) che[n]->RemoveColumn(k); }
		  }
		} free(sstJ);
	   } free(sstI);
	} return NumRemoved;
}

BooLean	mcs_typ::ResurrectRejectSeq(Int4 sq)
// returns true if sq is assigned to root set or reject set else returns false
//  reassign sequence sq to root if it is in reject set.
{
  char  **HP=Hpt->RtnHyperPartition();
  Int4  n,RootG=1,RejectG=Hpt->NumSets();	
  if(sq == 0){
    for(sq = 1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
       if(!MemberSet(sq,GrpSet[RejectG])) continue; 
       DeleteSet(sq,GrpSet[RejectG]); // remove sq from reject set.
       AddSet(sq,GrpSet[RootG]); 
       for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition(HP[RootG][n],sq); }
    }
  } else {
    if(sq < 1 || sq > NumSeqsCMSA(TrueMainCMA)) print_error("ResurrectRejectSeq() input error");
    if(MemberSet(sq,Labeled)) return FALSE;
    if(MemberSet(sq,GrpSet[RootG])) return TRUE; // this sequence is already in root set.
    if(!MemberSet(sq,GrpSet[RejectG])) return FALSE; // this sequence is not in reject set.
    DeleteSet(sq,GrpSet[RejectG]); // remove sq from its set.
    AddSet(sq,GrpSet[RootG]); 
    for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition(HP[RootG][n],sq); }
  } this->StoreBest(); return TRUE; 
    
}

double	mcs_typ::SampleSeq(Int4 sq, double oLPR)
{
    const double NEG_INFINITY=-1.0e+100;
    if(MemberSet(sq,Labeled)) return oLPR;
    char	**HP=Hpt->RtnHyperPartition();
    double	nLPR,bestLPR;
    Int4	n,g,oG,bestG,lastG;
    // find the current set = oG.
    for(oG=0,g=1; g<= Hpt->NumSets(); g++){	
      if(MemberSet(sq,GrpSet[g])){ oG=g; break;} // this sequence is in set g.
    } lastG = oG;
    if(oG <= 0){
	fprintf(stderr,"sq = %d; random = %d..%d; oG = %d\n",
		sq,NumSeqsCMSA(TrueMainCMA)+1,NumSeqsCMSA(MainCMA),oG);
    	assert(oG > 0); // sequence must be in one of the sets!
    }
    DeleteSet(sq,GrpSet[lastG]); // remove sq from its set.

    // Obtain subLPRs for fast computing.
    double	*ColumnMap[256]; NEW(ColumnMap['o'],Hpt->NumBPPS()+3,double);
    NEW(ColumnMap['+'],Hpt->NumBPPS()+3,double); NEW(ColumnMap['-'],Hpt->NumBPPS()+3,double);
// update NullParameters before removing sequence...!
    double OmitMap=0.0;	// Map with sequence omitted from everything.
    for(n=1; n<= Hpt->NumBPPS(); n++){
	che[n]->SetPartition('o',sq);
	double *tmpdbl=che[n]->SubMap( ); ColumnMap['o'][n]=tmpdbl[0]; OmitMap+=tmpdbl[0]; 
	// can later speed up the below by only checking the, say, 5 best FG sets.
	che[n]->SetPartition('+',sq);
	tmpdbl=che[n]->SubMap( ); ColumnMap['+'][n]=tmpdbl[0]; 
	che[n]->SetPartition('-',sq);
	tmpdbl=che[n]->SubMap( ); ColumnMap['-'][n]=tmpdbl[0]; 
    }

    double	*SetLPR; NEW(SetLPR,Hpt->NumSets()+3,double);	// for storing possible LPRs.
    bestG=0; bestLPR=NEG_INFINITY; 
#if 0
    BooLean RestrictSets=FALSE;
    if(this->SampledSet && MemberSet(Hpt->NumSets(),SampledSet)) RestrictSets=TRUE;
#endif
    for(g=1; g<= Hpt->NumSets(); g++){	
        if(IsFailedSet[g]){ SetLPR[g]=NEG_INFINITY; continue; }
#if 0
	if(RestrictSets && !MemberSet(g,SampledSet)){ SetLPR[g]=NEG_INFINITY; continue; }
	this->TestSample(g); che[g]->PutSubLPRs(stderr); exit(1);
#endif
	double dp=0.0,dm=0.0,d_o=0.0;
        for(n=1; n<= Hpt->NumBPPS(); n++){
	   switch(HP[g][n]){	// Ratio of being in the FG vs BG at these positions.
	    case '+': dp += ColumnMap['+'][n]; break;
	    case '-': dm += ColumnMap['-'][n]; break;
	    case 'o': d_o += ColumnMap['o'][n]; break;
	    default: print_error("input error"); break;	
	   }
	}
	// double d0= dp + dm + d_o - OmitMap;
	double d0= dp + dm + d_o; SetLPR[g]=d0;
	if(d0 > bestLPR){ 
		// fprintf(stderr,"sq %d: LPR %.2f --> %.2f.\n",sq,bestLPR,nLPR);
		bestLPR=d0; bestG=g; 
	}
    }
    //************  Sample an elementary set for this sequence... ************
   Int4 sampledG=0;
   if(temperature > 0.0){	// then reset bestG based on sampling...
    assert(temperature <= MaxTemperature);
    double	*deltaLPR; NEW(deltaLPR,Hpt->NumSets()+3,double);
    double	sum=0.0,d,rand;
    double	D,K;
    if(temperature == 300.0) K=1.0; else K=300.0/temperature;
    for(g=1; g <= Hpt->NumSets(); g++){
        if(IsFailedSet[g]){ deltaLPR[g]=0.0; continue; }
	d=SetLPR[g]-bestLPR;	// Likelihood ratio vs best.
	D = d*K;
	if(D < -14){ deltaLPR[g]=0.0; } // d ~ 0.0: e^-14 < 1 millionth chance to sample. 
	else { 
	   d = exp(d);		// take the anti-log
	   if(temperature == 300.0){ deltaLPR[g]=d; sum += d; }
	   else {
		double y=(300.0/temperature);
		assert(!(d==0 && y <= 0)); assert(d >= 0);
		if(d==0) deltaLPR[g]=0;
		else {
		  d = pow(d,y);
		  // assert(isfinite(d));    // confirms that d is okay.
		  if(!isfinite(d)) assert(!"isfinite failed");
		  deltaLPR[g]=d; sum += d; 
		}
	   }
	}
    }
    rand = sum * ((double) Random()/(double) RANDOM_MAX);	// random number...
    for(sum=0.0, g=1; g <= Hpt->NumSets(); g++){
	sum += deltaLPR[g]; 
	if(sum >= rand){ sampledG=g; break; }
    } assert(sum >= rand); free(deltaLPR);
   } else sampledG=bestG; 

    //************* end of sampling section. ******************
    // SetPartition were probably not set right...
   //  DeleteSet(sq,GrpSet[lastG]); 
   AddSet(sq,GrpSet[sampledG]); 
   TotalLPR=SetLPR[sampledG];
   for(n=1; n<= Hpt->NumBPPS(); n++){
	che[n]->SetPartition(HP[sampledG][n],sq);
	Map[n]=ColumnMap[HP[sampledG][n]][n]; // <-- This should be the same as CalcTotalLPR( );
   } 
   free(SetLPR); 
   free(ColumnMap['o']); free(ColumnMap['+']); free(ColumnMap['-']); // move below later when use for map!
   
   BooLean StoreBestOK=FALSE;
   if(SaveBest) StoreBestOK=TRUE; else StoreBestOK=FALSE;
   for(n=1; n<= Hpt->NumBPPS(); n++){ if(che[n]->NumColumns() > MaxNumCol[n]) StoreBestOK = FALSE; }
   CheckValue(TotalLPR);
   if(StoreBestOK && TotalLPR > 0.0 && TotalLPR > BestLPR) StoreBest();
   return TotalLPR; // == return CalcTotalLPR(0);
}

#if 0	// Routine for Identifying the sets to be removed if LPR <= 0.

Int4	mcs_typ::FindSetForRemoval(Int4 failed_column)
// BPPS does not include a Misc Set.
{
	Int4 Nplus=0,Nmisc=0;
	Int2 *FailedPerRow;
	// 1. BPPS includes a Misc set.
	// 2. BPPS does not include a Misc Set.
	for(Int4 row=1; row <= Hpt->NumSets(); row++){
		char cell=Hpt->Cell(row,failed_column);
		if(cell == '+'){
		   Nplus++;
		   if(Hpt->TypeOfSet(row) == '?'){
			Nmisc++;
		   }
		}
	} fprintf(stdout,"Column %d: N+ = %d; N? = %d\n",failed_column,Nplus,Nmisc);
	return Nplus;
}
	
#endif

Int4	mcs_typ::RmUnfruitfulSets( )
// Find failed sets and remove the sequences and columns from those sets.
// For pmcBPPS program only unless mcBPPS used with -tree option.
// Then replace this with removing only the internal node for this category.
// Turn this back on to deal with redundant sets; this was used for the pmcBPPS procedure.
{
    Int4	sq,g,i,j,k,m,n,s,NumFailed=0;
    char	**HP=Hpt->RtnHyperPartition();

    assert(IsTreeHpt);
    assert((Hpt->NumSets() -1) == Hpt->NumBPPS()); // Must be true of a tree!
    //---------------- 1. Find failed sets: -----------------
    // Hpt->ChangeFailedCells(double *Map); // change '+' to '-' or 'o'.
    for(g=Hpt->NumSets() - 1; g > 1; g--){	// skip last set == Random.
	n=g; 	// look at column n == g. 
	if(Map[n] <= 0 && HP[g][n] == '+'){   // Found a Failed Set.
	      if(!IsFailedSet[g]){	     // if not previously found...
		 NumFailed++; IsFailedSet[g]=TRUE; 
	         // TransferAllSeqs(g,1); // afn 6/15/2012: transfer the sets to the main set.
	         assert(!IsFailedBPPS[n]); IsFailedBPPS[n]=TRUE;
                 for(j=1; j<= Hpt->NumSets(); j++){
		   if(j == g) continue;
		   if(HP[j][n] == '+'){	// column n & row j:
			Hpt->Change('+','-',j,n); // ReSetRelations(stderr); 
		   	assert(HP[j][n] == '-');
#if 0	// change BG for j 
			m = j;
                 	for(k=1; k <= Hpt->NumBPPS(); k++){
			   if(k == 
		        if(HP[k][m] == 'o' && HP[g]){
		   	}
#endif
		   }
	         } ReSetRelations(0);
	      } // RowsInColumn(Int4 col,Int4 &row_neg,Int4 &row_pos,Int4 &row_omit);
	}
    }
    //if(IsTreeHpt){ NumFailed += RemoveSimilarSets( ); }
    if(NumFailed == Hpt->NumSets()) print_error("SemiConvergent test failed...now exiting");

    //-------------------- 2. Move Failed Sets to Main set. -------------------------
    for(g=2; g <= Hpt->NumSets(); g++){	
      // WARNING: requires and assumes that the first set is the main set!!!!
      if(IsFailedSet[g] && CardSet(GrpSet[g]) > 0) TransferAllSeqs(g,1); 
    }
    //-------------------- 3. Remove all columns from Failed BPPS. -------------------------
    Int4 Length=LenSeq(che[1]->KeyE());
    for(n=1; n<= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n] && che[n]->NumColumns( ) > 0){
	  che[n]->SetMinNumColumns(0);	// this is required to remove all columns!!
	  for(k=1; k <= Length; k++){ if(che[n]->PttrnPos(k)) che[n]->RemoveColumn(k); }
	}
    } return NumFailed;
}

BooLean	mcs_typ::Sample(Int4 IterStart, Int4 IterEvolve, Int4 NumRounds,Int4 ColSampleStart,Int4 StartMod)
// returns TRUE if converged; else returns FALSE;
// 3. Perform Gibbs sampling on columns and sequences 
#if 0	//*****************************************************************
// Operations: remove or add columns; move sequences up or down.
set options:
   A & B != 0 (Intersecting but 	-->	
   	A < B	  (subset)	-->	only subpatterns allowed.
	A !< B	  (not subset)	-->	no common pattern positions allowed.
   A & B = 0 (disjoint)		-->	All common pattern positions allowed...
 					(As long as transitive non-existent)
	(A | C) & (B | D) != 0 	-->	C and D intersect. 
	Implies that FG BG sampling is coupled between sets.

 ========= HyperPartition: =========
      _Category_
Set:  1  2  3  4  5  6  7  8  9
  1:  +  +  +  -  o  o  -  -  -  RabA (4945)
  2:  +  +  -  +  +  -  -  -  -  Arl (1822)
  3:  +  +  -  +  -  +  -  -  -  GaITO (282)
  4:  +  +  o  o  o  o  -  -  -  OddRasLike (1811)
  5:  +  -  o  o  o  o  +  -  -  EFlike (8450)
  6:  +  -  o  o  o  o  -  +  -  ObgSF (10854)
  7:  +  -  o  o  o  o  -  -  +  SimibiA (6846)
  8:  +  o  o  o  o  o  o  o  o  OtherGPs (6773)
  9:  -  o  o  o  o  o  o  o  o  Random (24409)

#endif	//*****************************************************************
{
     double	nLPR,lastLPR,BestThisRound=-999999999.9;
     Int4	sq,g,i,j,k,m,n,s,t,x;
     BooLean	RmAllNeg=TRUE;

     assert(NumRounds >= 2);

     NumCalls++;	// how many times has this been called?
     // for(n=1; n<= Hpt->NumBPPS(); n++) che[n]->SpeakUp();
     assert(nAlpha(AB) == 20);
     this->CalcTotalLPR(0);  // Calculates all Map[n].
#if 0
// this->PutPttrnVsConsSeq(stderr,"debug 00");
CheckPttrnCsqMatch("mcs->Sample() debug 0");       // ZZZZZZZZZZZZ
this->PutHyperPartition(stderr);
// this->PutPttrnVsConsSeq(stderr,"debug XX");
#endif
#if 0
#elif 1
     DidRestoreBest=FALSE;	// if do sampling need to restore new Best.
#elif 0
     this->RevertToBest();
#else
     this->RestoreBest();
#endif
// this->PutPttrnVsConsSeq(stderr,"debug 11");
// CheckPttrnCsqMatch("mcs->Sample() debug 1");       // ZZZZZZZZZZZZ

     // 1. Find the optimal or nearly-optimal partition given the seed pattern.
  for(Int4 iter=IterStart; iter <= NumRounds; iter++){
     char	tmpstr[100];

     sprintf(tmpstr,"improvement in LPR (iter S%d)",iter);
     fprintf(stderr,"\n***************** iter S%d *******************\n",iter);
     // lastLPR=CalcTotalLPR(stderr);  // Calculates all Map[n].
     lastLPR=CalcTotalLPR(0);  // Calculates all Map[n].
     // if(iter %2 == 1) BestThisRound=lastLPR;
     BestThisRound=-999999999.9;
     Int4 num_sq = NumSeqsCMSA(TrueMainCMA);
     num_sq = NumSeqsCMSA(MainCMA);
// CheckPttrnCsqMatch("mcs->Sample() debug 2");       // ZZZZZZZZZZZZ

     nLPR=lastLPR;

     dh_type	dH = dheap(num_sq+3,4);
     if(SampledSet != 0){
       set_typ setX = 0;
       for(s=1; s < Hpt->NumSets(); s++){
	   if(MemberSet(s,SampledSet))
	     { if(setX == 0) setX = CopySet(GrpSet[s]); else UnionSet(setX,GrpSet[s]); }
       } for(sq=1; sq <= num_sq; sq++){
	  if(MemberSet(sq,setX) && !MemberSet(sq,Labeled)) insrtHeap(sq,(keytyp) Random(),dH); 
       } NilSet(setX);
       fprintf(stderr,"%d sequences from %d nodes being sampled.\n",ItemsInHeap(dH),CardSet(SampledSet));
     } else {
       for(sq=1; sq <= num_sq; sq++){
	  if(!MemberSet(sq,Labeled)) insrtHeap(sq,(keytyp) Random(),dH); 
       }
     }
     if(efp) fprintf(efp,"sampling sequences....\n");
     double lpr0=lastLPR; 
     Int4 ItmsInHeap=ItemsInHeap(dH);
     TargetMod=ItmsInHeap/10; TargetMod=MAXIMUM(Int4,TargetMod,5);

     // Int4 Mod,TargetMod,TempMod,toggle=0,Toggle=0;
     Int4 Mod,TempMod,toggle=0,Toggle=0;
     double TempInc,MinDelta,delta;
     TempMod=(Int4) ceil((double) ItemsInHeap(dH)/(double) 100); 
     TempMod=MAXIMUM(Int4,TempMod,5); TempInc=1.0; // lower temperature ~100 degrees.
     Mod=StartMod; MinDelta=0.1;
     if(iter >= IterEvolve) DoEvolve( ); // evolving csq here...update within sampling columns only...
// CheckPttrnCsqMatch("mcs->Sample() debug 3");       // ZZZZZZZZZZZZ
     if(iter >= ColSampleStart) SampleColumns(TRUE);
else if(Evolve) UpdateCSQ(stderr);
// CheckPttrnCsqMatch("mcs->Sample() debug 4");       // ZZZZZZZZZZZZ
     if(cfp) fprintf(cfp,"%d %.1f %.1f %d\n",Iteration,lastLPR,temperature,TotalColumns( ));
     for(i=1,j=1,t=1; (sq=delminHeap(dH)) != NULL; i++,j++,t++) {
	lpr0=nLPR; nLPR=SampleSeq(sq,nLPR);  Iteration++;
	if(nLPR > BestThisRound) BestThisRound=nLPR;
        if(nLPR > 0.0) SaveBest=TRUE;	// start saving the best configuration.
        if(iter >= ColSampleStart && j % Mod == 0){	// sample a new column configuration...
	   // need cfp here to avoid spike in LPR due to > MaxNumCols prior to Col. Sampling.
           if(cfp) fprintf(cfp,"%d %.1f %.0f %d\n",Iteration,nLPR,temperature,TotalColumns( ));
           if(ifp){
        	fprintf(ifp,"%d",Iteration);
		for(s=1; s< Hpt->NumSets(); s++){ fprintf(ifp,"\t%d",CardSet(GrpSet[s])); }
		fprintf(ifp,"\t%d\n",CardSet(GrpSet[s])-NumRandom);
	   }
     	   // fprintf(stderr,"sampling columns....\n");
	   lpr0=nLPR; 
	   if(0 && temperature >= 200.0) nLPR=SampleColumns(TRUE);	// Use neg. patterns to optimize LPR.
	   else nLPR=SampleColumns(FALSE);	// Add only patterns that contribute to the LPR.
	   if(nLPR > BestThisRound) BestThisRound=nLPR;
           if(nLPR > 0.0) SaveBest=TRUE;	// start saving the best configuration.
	   Iteration++;
// CheckPttrnCsqMatch("mcs->Sample() debug 5");       // ZZZZZZZZZZZZ
	}
else if(Evolve) UpdateCSQ();
	delta=100.0*(nLPR-lpr0)/lpr0;	// compute the % increase in the LPR. 
	if(t % TempMod == 0){
	   if(nLPR > 0.0 && (delta < MinDelta) && temperature > 0.0) {
		temperature -= TempInc;
     		if(temperature < MinTemperature) temperature = 0.0;
	   } t=0;
	}
	if(j % Mod == 0){
	   Toggle++; if(Toggle > 1) Toggle=0;
	   if(Toggle==0 && Mod < TargetMod) Mod++;
	   j=0; toggle++; if(toggle > PrintToggle) toggle=0;
	   if(toggle==0){
	     x=this->NumFailedNodes();
     	     fprintf(stderr,
		"%d(%d/%d): LPR = %.2f (Last = %.2f; gained %.3f%c)(%.1f K)(%d cols)%d:%d",
			iter,i,ItmsInHeap,nLPR,lpr0,delta,'%',temperature,TotalColumns(),Mod,TargetMod);
		if(x > 0) fprintf(stderr," (%d failed).\n",x); else fprintf(stderr,".\n");
	      // if(iter == 1) PutHyperPartition(stderr);
	    }
	}
     } Nildheap(dH);
     ModStart=Mod;
#if 0
this->PutHyperPartition(stderr);
CheckPttrnCsqMatch("mcs->Sample() debug 7");       // ZZZZZZZZZZZZ
#endif
     if(SampledColumn > 0) SampleHpt(stderr,SampledColumn);
#if 1	// Evolving csq here...
     if(iter >= IterEvolve) DoEvolve( );
#endif
     if(cfp) fflush(cfp); if(ifp) fflush(ifp); 
     // PutHyperPartition(outfp); fflush(outfp);
     // PutHyperPartition(stderr);

     double increase=(nLPR -lastLPR);
     fprintf(stderr,"%d: LPR = %.2f (Last = %.2f; gained %.3f)(%.1f K)\n",
			iter,nLPR,lastLPR,increase,temperature);
     fflush(stderr);
     if(iter == 1 && temperature > 300.0) temperature=300.0;
     // if(iter > 1 && (temperature < 100.0 || increase < 0.1)) 	//  increase less than 1/10th %.
     if(iter > 1 && (temperature < 100.0 && increase < 0.1)) 	//  increase less than 1/10th %.
     {
        if(iter >= (NumRounds-1) && BestThisRound <= BestLPR) return TRUE;
						// no improvement after 3 iters --> converged.
        increase=((nLPR -lastLPR)/lastLPR);
        double ppb = ((double)ppb_increase/1000000000.0);
        if(increase < ppb){ return TRUE; }
        // if(iter > 4 && temperature < 50.0 && increase < ppb){ return TRUE; }
     }
     fprintf(stderr,"%d: %d sequences labeled\n",iter,CardSet(Labeled)-NumRandom);
#if 1	// terminate after this iteration if failed nodes present...
     if(this->NoFailureMode && (this->NumFailedNodes() > 0 || nLPR <= 0.0)) return TRUE;
#endif
// CheckPttrnCsqMatch("mcs->Sample() debug 8");       // ZZZZZZZZZZZZ
    } // end of iter "for" loop 
    return FALSE;
}

