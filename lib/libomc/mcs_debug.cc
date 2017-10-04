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

void	mcs_typ::PutPttrnVsConsSeq(FILE *fp,char *msg)
{
	fprintf(fp,"============= %s =============\n",msg);
        for(Int4 i=1; i <= Hpt->NumBPPS(); i++){ this->PutPttrnVsConsSeq(fp,i); }
}

void	mcs_typ::PutPttrnVsConsSeq(FILE *fp,Int4 i)
{
	Int4    s,len=LengthCMSA(1,TrueMainCMA);
	unsigned char r;
        e_type xE=che[i]->KeyE( );
        sst_typ *xsst=RtnCopyOfWorkingSST(i);
	fprintf(fp,"-------- Working %d: --------\n",i);
        for(s=1; s <= len; s++){
		if(xsst[s] == 0) continue;
                r=ResSeq(s,xE); PutSST(fp,xsst[s],AB); 
		fprintf(fp,"%d\t%c",s,AlphaChar(r,AB));
		for(r = 1; r <= nAlpha(AB); r++){
		   if(MemSset(r,xsst[s])){
			double d=(double)che[i]->GetResWtFG(s,r)/(double) che[i]->GetWtFactor();
			fprintf(fp,"\t%c=%.1f",AlphaChar(r,AB),d);
		   }
		} fprintf(fp,"\n");
	} fprintf(fp,"\n"); free(xsst);

	fprintf(fp,"-------- Best %d: ---------\n",i);
        e_type bE=BestCsq[i];
	if(bE==0) return;
	sst_typ *zsst=best_sst[i];
	for(s=1; s <= len; s++){
		if(zsst[s] == 0) continue;
               	r=ResSeq(s,bE);
                PutSST(fp,zsst[s],AB); fprintf(fp,"%d\t%c\n",s,AlphaChar(r,AB));
	} fprintf(fp,"\n\n");
}

void	mcs_typ::CheckPttrnCsqMatch(char *msg)
// See whether patterns are consistent with Csqs
{
	Int4    iter,i,k,len=LengthCMSA(1,TrueMainCMA);
	unsigned char r;
        for(i=1; i <= Hpt->NumBPPS(); i++){
	  for(iter=1; iter <= 2; iter++){
	    e_type xE;
	    sst_typ *xsst=0;
	    if(iter==1){ xE=che[i]->KeyE( ); xsst=RtnCopyOfWorkingSST(i); } // working pattern & consensus.
	    else { xE=BestCsq[i]; xsst=best_sst[i]; if(xE == 0) return; } // best pattern & consensus.
            assert(LenSeq(xE)==len);
            for(k=1; k <= len; k++){
                if(xsst[k]==0) continue;
                r=ResSeq(k,xE);
                if(!MemSset(r,xsst[k])){
		   if(iter == 1) fprintf(stderr,"============= Working: %s =============\n",msg);
		   else fprintf(stderr,"============= Best: %s =============\n",msg);
                   PutSeq(stderr,xE,AB);
                   if(IsBestRestored()) fprintf(stderr,"best is restored\n");
                   fprintf(stderr,"%d.",i); PutSST(stderr,xsst[k],AB);
                   fprintf(stderr,"%d != '%c'\n",k,AlphaChar(r,AB));
		   this->PutHyperPartition(stderr);
		   FILE *fp=open_file("junk",".cma","w");
		   PutInSetCMSA(fp,GrpSet[i],TrueMainCMA); fclose(fp);
		   e_type csqE=GetSeqAsCsqCMSA(GrpSet[i],TrueMainCMA);
		   AlnSeqSW(stderr,11,1, csqE, xE, AB);
		   PutSeq(stderr,csqE,AB); PutSeq(stderr,xE,AB); 
		   fprintf(stderr,"%d: submap = %g\n",i,che[i]->SubMap( ));
#if 0
		   for(Int4 s=1; s <= len; s++){
			if(xsst[s] == 0) continue;
                	r=ResSeq(s,xE);
                   	PutSST(stderr,xsst[s],AB); fprintf(stderr,"%d\t%c\n",s,AlphaChar(r,AB));
		   } fprintf(stderr,"\n");
#else
		   this->PutPttrnVsConsSeq(stderr,i);
#endif
                   assert(MemSset(r,xsst[k]));
                }
            } if(iter==1) free(xsst);
	   }
         }
}

BooLean	mcs_typ::PutMapContributions2(FILE *fp,lpr_typ *xlpr)
// put the contributions to the Map for each set and intermediate node.
{
     assert(IsTreeHpt);	// make sure this is a tree...
     Int4	n,g,sq,num_sq = NumSeqsCMSA(TrueMainCMA);
     char	x,**HP=Hpt->RtnHyperPartition();
     BooLean	rtn=TRUE;

     double	lpr,lpr0,*subLLR;
     NEW(subLLR,Hpt->NumBPPS() +3,double); CalcTotalLPR(0,FALSE); 
     for(g=0; g <= Hpt->NumBPPS(); g++){ subLLR[g]=Map[g]; }
     fprintf(fp,"%d.%s(%d): (*) %.2f\n",1,Hpt->ElmntSetName(1),CardSet(GrpSet[1]),subLLR[1]);
     set_typ  *FullSet=this->RtnSubTreeSeqSet( );
#if 1	// Debug.
     set_typ	TmpFG=CopySet(GrpSet[1]),TmpBG=CopySet(GrpSet[1]),TmpSet=CopySet(GrpSet[1]);
     set_typ	FG=CopySet(GrpSet[1]),BG=CopySet(GrpSet[1]);
#endif
     for(g=2; g < Hpt->NumSets(); g++){	
#if 1	// Debugging lpr_typ:
	if(Hpt->TypeOfSet(g) != '?') continue;
	assert(FullSet[g]);
	che[g]->PutParametersBPPS(stderr);
	double D,d;
	// if(SetFG && SetFG[g]) PutSet(stderr,SetFG[g]);
	// if(SetBG && SetBG[g]) PutSet(stderr,SetBG[g]);
        ClearSet(TmpFG); ClearSet(TmpBG); ClearSet(TmpSet);
        for(n=1; n < Hpt->NumSets(); n++){ 
	    if(MemberSet(n,SetFG[g])){
		if(FullSet[n]){ UnionSet(TmpFG,FullSet[n]); }
		else { UnionSet(TmpFG,GrpSet[n]); }
	    }
	    if(MemberSet(n,SetBG[g])){
		if(FullSet[n]){ UnionSet(TmpBG,FullSet[n]); }
		else { UnionSet(TmpBG,GrpSet[n]); }
	    }
	} UnionSet3(TmpBG,GrpSet[g],BG); IntersectNotSet(TmpFG,GrpSet[g],FG);
#if 1
	sst_typ *csst=this->RtnCopyOfSST(g);
	// xlpr->PutParameters(stderr,TmpFG,TmpBG,csst,FALSE,'m');
	bpps_typ *bpps=che[g]->BPPS(); 
	e_type qE=bpps->RtnQuerySeq();
	xlpr->CompareBPPS(stderr,bpps,TmpFG,TmpBG,csst,FALSE,'m');
	xlpr->CompareSqWts(stderr,TmpFG,TmpBG,che[g]->GetResWtsFG(),che[g]->GetResWtsBG(),FALSE);
// continue;
	D=xlpr->CalcSetvsPttrnLPR(0,TmpFG,TmpBG,csst,FALSE,'m',qE); 
	d=xlpr->CalcSetvsPttrnLPR(0,FG,BG,csst,FALSE,'m',qE); free(csst); NilSeq(qE);
#else
	sst_typ *Xsst=xlpr->GetOptPttrnLPR(0,TmpFG,TmpBG,FALSE,D,DefaultMaxCol,'m');
	sst_typ *xsst=xlpr->GetOptPttrnLPR(0,FG,BG,FALSE,d,DefaultMaxCol,'m');
	free(xsst); free(Xsst);
#endif
	fprintf(stderr,"\n%d: lpr=%.1f-%.1f=%.1f (%d|%d)(%d|%d)\n",
		g,D,d,D-d,CardSet(TmpFG),CardSet(TmpBG),CardSet(FG),CardSet(BG)); 
	fprintf(stderr," che[%d]: lpr=%.1f (%d|%d).\n",
		g,subLLR[g],che[g]->RtnCardFG(),che[g]->RtnCardBG()); 
#endif
	assert(g <= Hpt->NumBPPS());
	if(IsFailedSet[g]) continue;
	// Int4 card=CardSet(GrpSet[g]);
	if(FullSet[g]) fprintf(fp,"%d.%s(%d|%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]),CardSet(FullSet[g]));
	else fprintf(fp,"%d.%s(%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]));
        for(n=1; n < g; n++){ 
#if 1
	   if(IsFailedBPPS[n]) continue;  // this should never be true;
	   if(HP[g][n] == '+'){
		assert(Hpt->TypeOfSet(n) == '?'); // this should correspond to an internal node.
		x='-';  // Put in background partition.
		// x='o';  // Omit from analysis.

		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," (%d) %.1f",n,lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); } 
		} 

		if(FullSet[g] == 0) continue;
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet[g])){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," [%.1f]",lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet[g])){ che[n]->SetPartition('+',sq); } 
		} 
	   }
#else
	   if(RtnContribLLR(g,n,lpr)){ fprintf(fp," (%d) %.1f",n,lpr); }
#endif
	} fprintf(fp," (*) %.1f\n",subLLR[n]);
     } 
     for(g=1; g < Hpt->NumSets(); g++) if(FullSet[g]) NilSet(FullSet[g]);
     fprintf(fp,"\n"); free(subLLR); free(FullSet);
     return rtn;
}

BooLean	mcs_typ::CheckInputSets()
{
	Int4	i,j;
        // 1. Make sure that input sequence sets are compatible with other input.
        assert(num_passed_in_sets > 0); assert(passed_in_sets); assert(passed_in_sets[1]);
        Int4 size_passed_in_set=SetN(passed_in_sets[1]);
        if((size_passed_in_set - 1) != (NumSeqsCMSA(TrueMainCMA) + NumRandom)){
		fprintf(stderr,"size_passed_in_set = %d; NumSeqsCMSA(TrueMainCMA) = %d; NumRandom = %d\n",
			size_passed_in_set,NumSeqsCMSA(TrueMainCMA),NumRandom);
        	assert((size_passed_in_set - 1) == (NumSeqsCMSA(TrueMainCMA) + NumRandom));
	}
        for(i=2; i <= num_passed_in_sets; i++){
                assert(SetN(passed_in_sets[i]) == size_passed_in_set);
        } set_typ SetU=MakeSet(size_passed_in_set); ClearSet(SetU);
        for(i=1; i <= num_passed_in_sets; i++){
           for(j=i+1; j <= num_passed_in_sets; j++){
                Int4 X=CardInterSet(passed_in_sets[i],passed_in_sets[j]);
                if(X != 0) print_error("mcs_typ input error: sequence sets not disjoint");
           } UnionSet(SetU,passed_in_sets[i]);
        } assert(!MemberSet(0,SetU));
#if 0
        if(NumDisplayCMA == num_passed_in_sets){
           if(CardSet(SetU) != NumSeqsCMSA(TrueMainCMA)){
                fprintf(stderr,"Card(SetU) = %d; NumSeqsCMSA(TrueMainCMA) = %d\n",
                        CardSet(SetU),NumSeqsCMSA(TrueMainCMA));
                assert(CardSet(SetU) == NumSeqsCMSA(TrueMainCMA));
           } assert(CardInterSet(SetU,InitSet[0]) == 0);     // no random seqs in other partitions & vice versa.
	} NilSet(SetU); return TRUE;
#else
        if(CardSet(SetU) != NumSeqsCMSA(MainCMA)){
                fprintf(stderr,"Card(SetU) = %d; NumSeqsCMSA(MainCMA) = %d\n",
                        CardSet(SetU),NumSeqsCMSA(MainCMA));
                assert(CardSet(SetU) == NumSeqsCMSA(MainCMA));
        } IntersectNotSet(SetU,RandomSet);
	if(CardSet(SetU) != NumSeqsCMSA(TrueMainCMA)){
                fprintf(stderr,"Card(SetU) = %d; NumSeqsCMSA(TrueMainCMA) = %d\n",
                        CardSet(SetU),NumSeqsCMSA(TrueMainCMA));
                assert(CardSet(SetU) == NumSeqsCMSA(TrueMainCMA));
	} NilSet(SetU); return TRUE;
#endif
}

BooLean mcs_typ::ChecksOut()
// check to see whether anything worthwhile was found.
{
        // new routine for defining previously identified sets.
        if(DidRestoreBest == FALSE) return FALSE;
        Int4    n,g;
        for(n=0,g=1; g<= Hpt->NumSets(); g++){
                if(Hpt->TypeOfSet(g) != '?' && g != Hpt->NumSets() && CardSet(BestSet[g]) > 0) n++;
        } if(n > 0) return TRUE; else return FALSE;
}

BooLean	mcs_typ::ConsistencyCheck()
{
	Int4	i,j,sq;
	set_typ SetI,SetJ,SetU,SetX,SetY,SetRm,FG,BG,Rm,U;
	// SetI=MakeSet(SetSize); SetJ=MakeSet(SetSize); 

	// 1. Make sure that GrpSets are disjoint.
	SetX=MakeSet(SetSize); SetY=0;
	SetU=MakeSet(SetSize); ClearSet(SetU);
	for(i=1; i<=Hpt->NumSets(); i++){
	   SetI=GrpSet[i];
	   UnionSet(SetU,SetI);	// SetU = SetU U SetI.
	   for(j=i+1; j<=Hpt->NumSets(); j++){
	   	SetJ=GrpSet[j];
		IntersectSet1(SetI,SetJ,SetX); // modifies SetX to SetI intersect SetJ 
		assert(CardSet(SetX) == 0);
	   }
	}
	// 2. Make sure that GrpSets include all of the sequences.
	if(CardSet(SetU) != SetSize-1){
		fprintf(stderr,"CardSet(SetU) = %d != %d\n",CardSet(SetU),SetSize-1);
		SetY=MakeSet(SetSize); FillSet(SetX); 
		IntersectNotSet(SetX,SetU,SetY); // SetY = SetX intersect not SetU
		PutSet(stderr,SetY);
		return FALSE;
	}
	
	// look at BestSets...
	U=MakeSet(SetSize); ClearSet(U);
	for(i=1; i<=Hpt->NumSets(); i++){
	   SetI=BestSet[i];
	   UnionSet(U,SetI);	// U = 'U' U SetI.
	   for(j=i+1; j<=Hpt->NumSets(); j++){
	   	SetJ=BestSet[j];
		IntersectSet1(SetI,SetJ,SetX); // modifies SetX to SetI intersect SetJ 
		assert(CardSet(SetX) == 0);
	   }
	}
	if(CardSet(U) != SetSize-1){
		fprintf(stderr,"CardSet(U) = %d != %d\n",CardSet(U),SetSize-1);
		if(SetY==0) SetY=MakeSet(SetSize); FillSet(SetX); 
		IntersectNotSet(SetX,U,SetY); // SetY = SetX intersect not U
		// PutSet(stderr,SetY);
	}

	UInt4    Nfg,Nbg,Nf,Nb,Nr,num_sq = NumSeqsCMSA(MainCMA);
	for(i=1; i<=Hpt->NumBPPS(); i++){
	   // if(IsFailedBPPS[i]) continue;  // These che[i] objects have been discarded and can be ignored.
	   FG=che[i]->RtnFG_Set(); BG=che[i]->RtnBG_Set(); Rm=che[i]->RtnRmSet();
           // set_typ che[i]->Gld=RtnGoldSet( );
	   for(Nf=Nb=Nr=0,sq=1; sq <= num_sq; sq++){
		assert(MemberSet(sq,FG) || MemberSet(sq,BG) || MemberSet(sq,Rm));
		assert(!(MemberSet(sq,FG) && MemberSet(sq,BG)));
		assert(!(MemberSet(sq,FG) && MemberSet(sq,Rm)));
		assert(!(MemberSet(sq,BG) && MemberSet(sq,Rm)));
		if(MemberSet(sq,FG)) Nf++; else if(MemberSet(sq,BG)) Nb++;
		else { assert(MemberSet(sq,Rm)); Nr++; }
	   } assert(Nf == che[i]->RtnCardFG()); assert(Nb == che[i]->RtnCardBG());
	   assert(Nf == CardSet(FG)); assert(Nb == CardSet(BG));
	   for(Nfg=Nbg=0,j=1; j <= Hpt->NumSets(); j++){
		// if(IsFailedSet[j]) continue;
		if(MemberSet(j,SetFG[i])) Nfg += CardSet(GrpSet[j]);
		if(MemberSet(j,SetBG[i])) Nbg += CardSet(GrpSet[j]);
	   } if(Nf != Nfg || Nb != Nbg){
		// PutSetRelations(stderr);
		fprintf(stderr,"FG set (%d): ",i); PutSet(stderr,SetFG[i]); 
		fprintf(stderr,"\nBG set (%d): ",i);PutSet(stderr,SetBG[i]);
		fprintf(stderr,"\n%d: Nf = %d; Nfg = %d; Nb = %d; Nbg = %d\n",i,Nf,Nfg,Nb,Nbg);
		assert(Nf==Nfg); assert(Nb==Nbg);
	   }
	}

	assert(CardSet(SetU) == SetSize-1);
	if(CardSet(U)  != SetSize-1){ fprintf(stderr,"CardSet(U) = %d; SetSize = %d\n",CardSet(U),SetSize); }
	assert(CardSet(U) == SetSize-1);
	if(SetY) NilSet(SetY); NilSet(SetX); NilSet(SetU); NilSet(U);
	return TRUE;
}

