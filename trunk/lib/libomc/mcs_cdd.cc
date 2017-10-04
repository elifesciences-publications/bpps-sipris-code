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

BooLean mcs_typ::IsSeqEST_ENV(e_type E)
{
        char king,*result=strstr(E->info,"_EST");
        if(result==0){
		king=kingdomSeq(E); if(king == 0) return TRUE;
		if(strchr("MFEVBA",toupper(king)) == 0) return TRUE;
		return FALSE;
	} else {
                char *space=strchr(E->info,' ');
                if(space > result) return TRUE;
		else {
		  king=kingdomSeq(E); if(king == 0) return TRUE;
		  if(strchr(" MFEVBA",toupper(king)) == 0) return TRUE;
		  else return FALSE;
		}
        }
}

Int4	mcs_typ::FirstInSet(set_typ St,cma_typ cma)
// return the first good element in set St.
{
	double cut=0.0;
	for(Int4 i=1; i <= SetN(St); i++){
		if(MemberSet(i,St)){
		   e_type sE=TrueSeqCMSA(i,cma);
		   if(IsSeqEST_ENV(sE)) continue;
		   double prob=GetGappedProbCMSA(1,i,cma);
		   if(prob >= cut) return i; 
		}
	} return 0; 
}

Int4	*mcs_typ::RtnBestSeqs(set_typ *&GoodSet)
// for each Set get the best sequence.
{
    Int4	*rtn,sq,bst_sq,g,n,N=NumSeqsCMSA(TrueMainCMA);
    double	lpr,lpr0,diff,bst;
    char        **HP=Hpt->RtnHyperPartition();
    
    NEW(GoodSet,Hpt->NumSets()+3,set_typ);
    NEW(rtn,Hpt->NumSets()+3,Int4);
    for(g=1; g < Hpt->NumSets(); g++){	// skip Random set.
      GoodSet[g]=MakeSet(SetN(GrpSet[g])); ClearSet(GoodSet[g]);
      bst=-99999999; bst_sq=0;
      for(sq=1; sq <= N; sq++){
        if(MemberSet(sq,GrpSet[g])){
            lpr0 = CalcTotalLPR(0,FALSE);
            DeleteSet(sq,GrpSet[g]);
            for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition('o',sq); }
            lpr = CalcTotalLPR(0,FALSE);
            diff=lpr0-lpr;
	    e_type sE=TrueSeqCMSA(sq,TrueMainCMA);
	    // if(diff > 0.0 && !IsSeqEST_ENV(sE)) AddSet(sq,GoodSet);
	    if(diff > 0.0) AddSet(sq,GoodSet[g]);
#if 0   // adjust for sequence redundancy;  make sure this is a tree!
            UInt4   wt=che[g]->RtnAveSqIWt(sq),factor=che[g]->GetWtFactor();
            double Wt=(double)wt/(double)factor;
            assert(Wt <= 1.0 && Wt > 0.0);
            sqLLR[sq]=sqLLR[sq] * 1.0/Wt;
#endif
	    if(diff > bst && !IsSeqEST_ENV(sE)){ bst=diff; bst_sq=sq;  }
            AddSet(sq,GrpSet[g]);
            for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->SetPartition(HP[g][n],sq); }
        }
      } rtn[g]=bst_sq;
    } return rtn;
}

set_typ	mcs_typ::RtnPdbSet(Int4 g)
{
	assert(g > 0 && g < Hpt->NumSets());
    	Int4	sq,N=NumSeqsCMSA(TrueMainCMA);
	set_typ pdbSet=MakeSet(SetN(GrpSet[g])); ClearSet(pdbSet);
	for(sq=1; sq <= N; sq++){
           if(MemberSet(sq,GrpSet[g])){
            	if(PdbSeq(TrueSeqCMSA(sq,TrueMainCMA))){ AddSet(sq,pdbSet); }
	   }
	} return pdbSet;
}

Int4	mcs_typ::NoIndelsSet(set_typ NoIndels)
{
	cma_typ	cma=TrueMainCMA;
	assert(nBlksCMSA(cma) == 1);
	Int4	nins,inslen,ins,del,pos[3],sq,N=NumSeqsCMSA(cma);
	Int4	max,col,*Nino,*Ndel,*Nins,Len=LengthCMSA(1,cma);

	assert(SetN(NoIndels) >= N);
	ClearSet(NoIndels);

	NEW(Nino,Len+5,Int4); NEW(Nins,Len+5,Int4); NEW(Ndel,Len+5,Int4);
	for(col=0; col < Len; col++){
	    inslen=del=nins=0;
	    for(sq=1; sq <= N; sq++){
	        PosSiteCMSA(1, sq, pos, cma);
	        ins = InsertionCMSA(sq, pos[1]+col,cma);
	        if(ins > 0){ nins++; inslen+=ins; }
	        if(IsDeletedCMSA(sq,pos[1]+col, cma)){ del++; }
	    } Ndel[col+1]=del; 
	     if(col < (Len - 1)) Nino[col+1]=nins; Nins[col+1]=inslen; 
	}

	h_type HG=Histogram("number of residues inserted",0,Len,1);
	for(max=0,col=1; col < Len; col++){	// ignore extensions on the the end!!
	    if(Nins[col] > 0){
		IncdMHist(col,Nins[col],HG);
		if(max < Nins[col]) max = Nins[col];
	    }
	} PutHist(stdout,60,HG); NilHist(HG);

	// Int4 cut=10*N/Len;
	// Int4 cut=N/Len;
	Int4 cut=N/3;
	Int4 inc=1 + max/50;
	set_typ set=MakeSet(Len+5); ClearSet(set);
	HG=Histogram("residues inserted",0,max,inc);
	for(col=1; col < Len; col++){
	   if(Nins[col] > cut && col < (Len-9)) IncdHist(Nins[col],HG); 
	   else AddSet(col,set);
	}
	PutHist(stdout,60,HG); NilHist(HG);

	HG=Histogram("number of insertions",0,Len,1);
	for(col=1; col < Len; col++){
	    if(Nino[col] > 0) IncdMHist(col,Nino[col],HG);
	} PutHist(stdout,60,HG); NilHist(HG);

	HG=Histogram("number of deletions",0,Len,1);
	for(col=1; col <= Len; col++){
	    if(Ndel[col] > 0)  IncdMHist(col,Ndel[col],HG);
	} PutHist(stdout,60,HG); NilHist(HG);

	inc=1+N/100;
	HG=Histogram("numbers of indels per sq",0,Len,1);
	for(sq=1; sq <= N; sq++){
	    inslen=del=nins=0;
	    for(col=2; col <= Len-2; col++)
	    {
		if(!MemberSet(col,set)) continue;
	        PosSiteCMSA(1, sq, pos, cma);
	        ins = InsertionCMSA(sq, pos[1]+col-1,cma);
	        if(ins > 0){ nins++; inslen+=ins; }
	        if(IsDeletedCMSA(sq,pos[1]+col-1, cma)){ del++; }
	    } if(nins == 0 && del == 0) AddSet(sq,NoIndels);
	    IncdHist(nins+del,HG);
	 
	} PutHist(stdout,60,HG); NilHist(HG);
	free(Nino); free(Nins); free(Ndel); NilSet(set);
} 


void	mcs_typ::PutCDD(FILE *ofp)
// output sets for CDtree with first sequence for each descendant node of an internal node.
{
    Int4	n,g,N=NumSeqsCMSA(TrueMainCMA);
    char	tmp_str[200];
    Int4	max_size=100;
    set_typ	*GoodSet=0;
    e_type	sE=0;

    // ofp=open_file(infile,"_cdd.mma","w");
    Int4	*bstSq=RtnBestSeqs(GoodSet);
    set_typ	NoIndels=CopySet(GrpSet[1]); NoIndelsSet(NoIndels);
    for(g=1; g < Hpt->NumSets(); g++){	// skip Random set.
      if(CardSet(GrpSet[g]) > 0){	// output optimized minimal set if non-empty.
	set_typ TmpSet=CopySet(GrpSet[g]);
        // set_typ dst=Hpt->DescendantSet(g);  // all descendants of g.
	ClearSet(TmpSet);
	Int4 i,j=0,*P=0,first,ptr=1,*list=0;
	NEW(list,NumSeqsCMSA(MainCMA) +3, Int4);
	assert(Hpt->IsTree(P)); // first=FirstInSet(GrpSet[g],MainCMA);

	// best sequences.
	first=bstSq[g];
	if(first > 0){ j++; list[j] = first; AddSet(first,TmpSet);}
        for(Int4 gg=1; gg < Hpt->NumSets(); gg++){	// skip Random set.
           // if(MemberSet(gg,dst)) 	// if gg is a descendant of g.
	   if(P[gg] == g)	// if gg is a child of g.
	   {
		assert(GrpSet[gg] != 0);
		first=bstSq[gg]; // first=FirstInSet(GrpSet[gg],MainCMA);
		if(first > 0){ AddSet(first,TmpSet); j++; list[j] = first; }
	   }
	} free(P);

	// pdb sequences.
	set_typ pdbSet=RtnPdbSet(g);
	for(i=1; i <= N; i++){
	     if(MemberSet(i,TmpSet)) continue;	// if already on list.
	     if(MemberSet(i,pdbSet)){ AddSet(i,TmpSet); j++; list[j] = i; }
	} NilSet(pdbSet);

	// 'good' non-EST sequences.
	for(i=1; i <= N; i++){
	     if(j > max_size) break;
	     if(MemberSet(i,TmpSet)) continue;	// if already on list.
if(!MemberSet(i,NoIndels)) continue;
	     sE=TrueSeqCMSA(i,TrueMainCMA);
	     if(MemberSet(i,GoodSet[g]) && !IsSeqEST_ENV(sE))
			{ AddSet(i,TmpSet); j++; list[j] = i; }
	}

	// 'good' EST sequences.
	for(i=1; i <= N; i++){
	     if(j > max_size) break;
	     if(MemberSet(i,TmpSet)) continue;
if(!MemberSet(i,NoIndels)) continue;
	     sE=TrueSeqCMSA(i,TrueMainCMA);
	     if(MemberSet(i,GoodSet[g]) && !IsSeqEST_ENV(sE))
			{ AddSet(i,TmpSet); j++; list[j] = i; }
	}

	// worst non-EST sequences.
	for(i=1; i <= N; i++){
	     if(j > max_size) break;
	     if(MemberSet(i,TmpSet)) continue;
if(!MemberSet(i,NoIndels)) continue;
	     sE=TrueSeqCMSA(i,TrueMainCMA);
	     if(MemberSet(i,GrpSet[g]) && !IsSeqEST_ENV(sE))
		{ AddSet(i,TmpSet); j++; list[j] = i; }
	}

	// worst EST sequences.
	for(i=1; i <= N; i++){
	     if(j > max_size) break;
	     if(MemberSet(i,TmpSet)) continue;
if(!MemberSet(i,NoIndels)) continue;
	     if(MemberSet(i,GrpSet[g])) { j++; list[j] = i; }
	}

	ReNameCMSA(Hpt->ElmntSetName(g),MainCMA);
	PutSelectOneCMSA(ofp,0,list,MainCMA); free(list);
        // NilSet(dst);
	NilSet(TmpSet); 
      }
    } // fclose(ofp);
    NilSet(NoIndels);
    for(g=1; g < Hpt->NumSets(); g++){	
	if(GoodSet[g]) NilSet(GoodSet[g]);
    } free(GoodSet); free(bstSq);
}

