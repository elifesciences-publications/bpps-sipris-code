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

e_type  mcs_typ::RtnKeySeq(Int4 n)
// Return an array of consensus sequences.
{
        assert(n > 0 && n <= Hpt->NumBPPS());
        return che[n]->KeyE( );
}

e_type  *mcs_typ::RtnKeySeqs( )
// Return an array of consensus sequences.
{
        e_type  *KeySeq;
        NEW(KeySeq,Hpt->NumBPPS() +3, e_type);
        for(Int4 n=1; n<= Hpt->NumBPPS(); n++) KeySeq[n]=che[n]->KeyE( );
        return KeySeq;
}

char	*mcs_typ::RtnCopyOfPttrn(Int4 i)
{
	assert(i > 0 && i <= Hpt->NumBPPS());
	sst_typ *xsst=this->RtnCopyOfSST(i);
	char    *pttrn=SST2ArgStrAlpha(xsst,this->RtnLengthMainCMSA(),AB);
	free(xsst);
	return pttrn;
}

sst_typ	*mcs_typ::RtnCopyOfWorkingSST(Int4 n)
{
	assert(n > 0 && n <= Hpt->NumBPPS());
	Int4 i,len=this->RtnLengthMainCMSA();
	sst_typ *xsst; NEW(xsst,len+3,sst_typ);
	bpps_typ *bpps=che[n]->BPPS();
	sst_typ *tmpsst=bpps->RtnSST();
	Int4 k=bpps->LenPattern( );
	assert(k == len);
	for(i=1; i <= k; i++) xsst[i]=tmpsst[i];
	return xsst;
}

sst_typ	*mcs_typ::RtnCopyOfBestSST(Int4 n)
{
	assert(n > 0 && n <= Hpt->NumBPPS());
	Int4 i,len=this->RtnLengthMainCMSA();
	sst_typ *xsst; NEW(xsst,len+3,sst_typ);
        for(i=1; i <= len; i++) xsst[i]=best_sst[n][i];
	return xsst;
}

sst_typ	**mcs_typ::RtnCopyOfWorkingSSTs( )
{
	Int4 n,i,j,len=this->RtnLengthMainCMSA();
	sst_typ **xsst; NEWP(xsst,Hpt->NumBPPS()+3,sst_typ);
	for(n=1; n<= Hpt->NumBPPS(); n++) xsst[n]=this->RtnCopyOfWorkingSST(n);
	return xsst;
}

sst_typ	**mcs_typ::RtnCopyOfBestSSTs( )
{
	// assert(DidRestoreBest);
	Int4 n,i,j,len=this->RtnLengthMainCMSA();
	sst_typ **xsst; NEWP(xsst,Hpt->NumBPPS()+3,sst_typ);
	for(n=1; n<= Hpt->NumBPPS(); n++) xsst[n]=this->RtnCopyOfBestSST(n);
	return xsst;
}

sst_typ	*mcs_typ::RtnCopyOfSST(Int4 n)
// to be called after sampling is completed.
{
#if 0
	Int4 i,j,len=this->RtnLengthMainCMSA();
	sst_typ *xsst;
	
	assert(n > 0 && n <= Hpt->NumBPPS());
	NEW(xsst,len+3,sst_typ);
	if(DidRestoreBest==FALSE){
		bpps_typ *bpps=che[n]->BPPS();
		sst_typ *tmpsst=bpps->RtnSST();
		Int4 k=bpps->LenPattern( );
		assert(k == len);
		for(i=1; i <= k; i++) xsst[i]=tmpsst[i];
	} else {
                for(i=1; i <= len; i++) xsst[i]=best_sst[n][i];
	} return xsst;
#else
	if(DidRestoreBest) return RtnCopyOfBestSST(n);
	else return RtnCopyOfWorkingSST(n);
#endif
}

char	**mcs_typ::RtnCopyOfPttrns( )
{
	sst_typ **xsst=this->RtnCopyOfSSTs();
	char    **pttrn; NEWP(pttrn,Hpt->NumBPPS()+3,char);
	Int4	i;
	for(i=1; i <= Hpt->NumBPPS(); i++){
          pttrn[i]=SST2ArgStrAlpha(xsst[i],this->RtnLengthMainCMSA(),AB);
	} for(i=1; i <= Hpt->NumBPPS(); i++){ free(xsst[i]); } free(xsst);
	return pttrn;
}

sst_typ	**mcs_typ::RtnCopyOfSSTs( )
// to be called after sampling is completed.
{
	if(DidRestoreBest) return this->RtnCopyOfBestSSTs( );
	else return this->RtnCopyOfWorkingSSTs( );
#if 0
	Int4 n,i,j,len=this->RtnLengthMainCMSA();
	sst_typ **xsst;
	
	// assert(DidRestoreBest);
	NEWP(xsst,Hpt->NumBPPS()+3,sst_typ);
	for(n=1; n<= Hpt->NumBPPS(); n++){
		bpps_typ *bpps=che[n]->BPPS();
		sst_typ *tmpsst=bpps->RtnSST();
		Int4 k=bpps->LenPattern( );
		assert(k == len);
#if 0
		for(i=1; i <= k; i++){
		   if(tmpsst[i]){
			PutSST(stderr,tmpsst[i],AB); fprintf(stderr,"%d,",i);
		   }
		}
		char *s=SST2ArgStrAlpha2(tmpsst,k,AB);
		fprintf(stderr,"\n%d: %s\n",n,s);
#endif
		NEW(xsst[n],k+3,sst_typ);
		for(i=1; i <= k; i++) xsst[n][i]=tmpsst[i];
	    } else {
		NEW(xsst[n],len+3,sst_typ);
                for(i=1; i <= len; i++) xsst[n][i]=best_sst[n][i];
	    }
	} return xsst;
#endif
}

double  **mcs_typ::RtnCopyOfLPRs( )
// to be called after sampling is completed.
{
	  Int4 n,i;

	  assert(DidRestoreBest);
	  double total,**lpr; NEWP(lpr,Hpt->NumBPPS()+3,double);
	  for(total=0.0, n=1; n<= Hpt->NumBPPS(); n++){
		double *tmplpr=che[n]->SubMap( ); 
		Int4 k=che[n]->LengthPattern();
		NEW(lpr[n],k+3,double);
		for(i=0; i <= k; i++) lpr[n][i]=tmplpr[i];
		// lpr[n][0] = total lpr.
	  } return lpr;
}

double	*mcs_typ::ContributionsToLLR(BooLean global)
// compute the contribution of each sequence to the LLR.
{
     Int4 sq,n,g,num_sq=NumSeqsCMSA(TrueMainCMA);;
     char **HP=Hpt->RtnHyperPartition();
     double d,d_g,d_g0,*D; NEW(D,num_sq +3,double);
     for(sq=1; sq <= num_sq; sq++){
        for(g=1; g <= Hpt->NumSets(); g++){     // don't look at last (random) set!!
          if(GrpSet[g] && MemberSet(sq,GrpSet[g])){
             assert(Hpt->NumSets() == Hpt->NumBPPS()+1);
             if(!IsFailedSet[g] && g != Hpt->NumSets()){
		double  *Map0; NEW(Map0,Hpt->NumBPPS() +3, double);
                double lpr0 = CalcTotalLPR(0,FALSE);
		d_g0=Map[g];

                DeleteSet(sq,GrpSet[g]);
                for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition('o',sq); }

		d=lpr0-CalcTotalLPR(0,FALSE);
		d_g = d_g0 - Map[g];

                AddSet(sq,GrpSet[g]);
                for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition(HP[g][n],sq); }
             } else { d=-99999999.0; d_g=-99999999.0; }
	     break; // each sq is a member of only one set...
          } 
        } if(global) D[sq]=d; else D[sq]=d_g;
      } return D;
}

void	mcs_typ::UnLabelAllSeqs( )
// allow all sequences to be sampled between sets freely...
{
	if(Unlabeled) return;
	for(Int4 sq=1; sq <= SizeTrueMain; sq++){
		if(MemberSet(sq,Labeled)){ DeleteSet(sq,Labeled); } 
	} Unlabeled=TRUE;
}

void	mcs_typ::TransferAllSeqs(Int4 from, Int4 to)
{
	if(from == to) return;
	assert(from > 1 && from < Hpt->NumSets() && to > 0 && to <= Hpt->NumSets());
	Int4 n,sq,num_sq = NumSeqsCMSA(MainCMA);
	char **HP=Hpt->RtnHyperPartition();
	for(sq=1; sq <= num_sq; sq++){
	     if(MemberSet(sq,GrpSet[from])) {
		if(MemberSet(sq,Labeled)){ DeleteSet(sq,Labeled); } 
		// ^ remove these from Untouchable set.
		DeleteSet(sq,GrpSet[from]); AddSet(sq,GrpSet[to]);
		for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition(HP[to][n],sq); }
	     }
	}
}

set_typ	*mcs_typ::CopyOfPartitionSets_Private()
{
	// Copy and return each partition set within the hpt.
	Int4	c,r,g;
  	assert(DidRestoreBest);
	set_typ *RtnSet; NEW(RtnSet,Hpt->NumBPPS()+2,set_typ);
	set_typ Set=MakeSet(SetN(BestSet[1]));
	for(c=1; c<= Hpt->NumBPPS(); c++){
	   ClearSet(Set);
	   for(r=1; r<= Hpt->NumSets(); r++){
	     if(Hpt->RtnHyperPartition(r,c) == '+'){ UnionSet(Set,BestSet[r]); }
	   }
	   if(CardSet(Set) > 0){ RtnSet[c]=CopySet(Set);  }
	}  NilSet(Set);
	return RtnSet;
}

cma_typ	mcs_typ::RtnBstAsCMSA(Int4 n, char *name,sst_typ *xsst)
// return as a consensus seq. the best sequence corresponding to the FG for column n.
{
	unsigned char r;
	assert(n > 0 && n < Hpt->NumSets());
	// set_typ st=MakeSet(Hpt->NumSets()+1); this->ReSetSubTree(st,i,hpt);
	if(CardSet(BestSet[n]) < 1){
		fprintf(stderr,"mcs_typ::RtnBstAsCMSA( ): CardSet(BestSet[%d]) = %d.\n",n,CardSet(BestSet[n]));
		print_error("FATAL: failed to find a significant hierarchy.");
	}
	cma_typ cma=GetBestCsqAsCMSA(BestSet[n],TrueMainCMA);
	e_type  keyE=TrueSeqCMSA(1,cma); 
	assert(LenSeq(keyE)==LengthCMSA(1,TrueMainCMA));
	for(Int4 i=1; i <= LenSeq(keyE); i++){
	    r=ResSeq(i,keyE);
	    if(xsst && xsst[i]){ if(!MemSset(r,xsst[i])){ xsst[i] = SsetLet(r); } }
	} RenameCMSA(name,cma);
	return cma;
}

cma_typ	mcs_typ::RtnSeqAsCMSA(Int4 n, char *name,sst_typ *xsst)
{ return RtnSeqOrCsqAsCMSA('S',n,name,xsst); }

cma_typ	mcs_typ::RtnCsqAsCMSA(Int4 n, char *name,sst_typ *xsst)
{ return RtnSeqOrCsqAsCMSA('C',n,name,xsst); }

cma_typ	mcs_typ::RtnCsqSstAsCMSA(Int4 n, char *name,sst_typ *xsst)
{ return RtnSeqOrCsqAsCMSA('c',n,name,xsst); }

cma_typ	mcs_typ::RtnSeqOrCsqAsCMSA(char mode, Int4 n, char *name,sst_typ *xsst)
// return the consensus sequence corresponding to the FG for column n.
// if(xsst != 0) then redefine incompatible sets at pattern positions to match keyE.
{
	unsigned char r,x;
	e_type  keyE=0;
	if(mode == 'C' || mode == 'c'){
	   assert(n > 0 && n < Hpt->NumSets());
	   keyE=che[n]->KeyE( );
	} else if(mode == 'S'){
	   assert(n > 0 && n <= NumSeqsCMSA(TrueMainCMA));
	   keyE=GetSeqAsCsqCMSA(n,TrueMainCMA);
	   assert(xsst != 0);
	} else print_error("RtnSeqOrCsqAsCMSA( ) input error");
	FILE *fp=tmpfile(); 
	fprintf(fp,"[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",name);
        fprintf(fp,"(%d)",LenSeq(keyE));
        for(Int4 i=1; i <= LenSeq(keyE); i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",LenSeq(keyE),LenSeq(keyE));
        fprintf(fp,">%s consensus\n{()",name);
	for(Int4 i=1; i <= LenSeq(keyE); i++){
	    r=ResSeq(i,keyE);
	    if(xsst){
 		if(xsst[i]){
		    if(r==0 || !MemSset(r,xsst[i])){
	              for(r=0,x=1; x <= nAlpha(AB); x++) if(MemSset(x,xsst[i])){ r=x; break; }
		      assert(r != 0);
		    } // else if(r == 0) then fixed below; else r is a member of xsst[i].
		} else if(mode != 'c' && r != 0 && !MemSset(r,xsst[i])) xsst[i]=SsetLet(r); 
	    }
	    if(r==0){ r=AlphaCode('A',AB); }	// change 'X' residues to 'A'.
	    fprintf(fp,"%c",AlphaChar(r,AB));
	} fprintf(fp,"()}*\n\n_0].\n"); rewind(fp); 
	cma_typ	cma=ReadCMSA(fp,AB); fclose(fp);
	return cma;
}

set_typ	*mcs_typ::CopyOfBestTreeSets_Private( )
// Return significant sets only...
{
	set_typ *RtnSet;
	Int4	g;
#if 0
	if(!DidRestoreBest) RestoreBest();
  	// assert(DidRestoreBest);
#else
// fprintf(stderr,"** DidRestoreBest = %d\n",DidRestoreBest);
	if(DidRestoreBest==FALSE) return CopyOfSeqSets_Private();
#endif
	assert(IsTreeHpt);
	assert((Hpt->NumSets()-1) == Hpt->NumBPPS());
	CalcTotalLPR(0,FALSE);  // Calculates all Map[n].
	//  fprintf(outfp,"============ Best FD-table ============\n");
	// PutHyperPartition(outfp); fflush(outfp);
	// PutHyperPartition( ); 

	NEW(RtnSet,Hpt->NumSets()+2,set_typ);
	for(g=1; g<= Hpt->NumSets(); g++){
#if 0
	   if(g > 1 && g < Hpt->NumSets() && (Map[g] <= 5.0 || che[g]->NumColumns( ) < 4)){
		CopySet(RtnSet[1],BestSet[g]);  // copies BestSet[g] to MainSet.
	   } else {  // if(CardSet(BestSet[g]) > 0){ // allow empty sets for internal nodes...
		RtnSet[g]=MakeSet(SetN(BestSet[g])); ClearSet(RtnSet[g]); 
		CopySet(RtnSet[g],BestSet[g]);  // copies BestSet[g] to RtnSet[g].
	   } 
#else
	   RtnSet[g]=MakeSet(SetN(BestSet[g])); ClearSet(RtnSet[g]);  
	   CopySet(RtnSet[g],BestSet[g]);  // copies BestSet[g] to RtnSet[g].
#endif
	} return RtnSet;
}

set_typ	*mcs_typ::CopyOfSeqSets_Private()
{
	set_typ *RtnSet;
	Int4	g;
  	// assert(DidRestoreBest);
	NEW(RtnSet,Hpt->NumSets()+2,set_typ);
	for(g=1; g <= Hpt->NumSets(); g++){
	  if(DidRestoreBest){
	   RtnSet[g]=MakeSet(SetN(BestSet[g])); ClearSet(RtnSet[g]); 
	   if(CardSet(BestSet[g]) > 0){
		CopySet(RtnSet[g],BestSet[g]);  // copies BestSet[g] to RtnSet[g].
	   } 
	  } else {
	   RtnSet[g]=MakeSet(SetN(GrpSet[g])); ClearSet(RtnSet[g]); 
	   if(CardSet(GrpSet[g]) > 0) CopySet(RtnSet[g],GrpSet[g]);
	  }
	} return RtnSet;
}

