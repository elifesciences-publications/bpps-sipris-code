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


void	omc_typ::FlattenHierarchy( )
{
        mcs->PutHyperPartition(stderr);
        mcs_typ *rtn_mcs=FlattenHiearchy(); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
        mcs->PutHyperPartition(stderr);
        rtn_mcs->NoFailureMode=TRUE;
        rtn_mcs->SaveBest=TRUE;
        // mcs->Sample(2,2,2,2);
        mcs->PutHyperPartition(stderr);
        // sprintf(str,"%s_flat",infile); PutCheckPoint(str,FALSE);
        FILE *fp=open_file(infile,"_flat.out","w"); mcs->PutHyperPartition(fp); fclose(fp); fp=0;
}

void	omc_typ::RandomStart(Int4 NumRandNodes)
// create and optimize an arbitrary hierarchy using real sequences.
// This is redundant with 'a' test option below; eliminate one of these eventually.
{
        Int4    id,i,j,n,N,s,iter=0,result=0,No,mode=0;
        double  best_lpr,last_lpr,d,D;
        mcs_typ *xmcs=0,*rmcs=0;

        SetStringency(stringency);
        assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
        mcs->NoFailureMode=TRUE;  // if node configuration subLPR <= 0.0, then reject it.
        mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
        if(Evolve) mcs->DoEvolve(); else mcs->DoNotEvolve();

	// 1. create a random hpt.
	Hpt->Put(stderr); // hpt->PutRandomize(stdout,Random);
        hpt_typ *rhpt=0;
	rhpt=Hpt->Randomize(NumRandNodes);
        rhpt->Put(stderr);
        // rhpt->PutSorted(stderr);
	// hpt_typ *shpt=rhpt->Sort( );
	// 2. rename sets and create sma files.
	cma_typ *sma; NEW(sma,rhpt->NumSets() + 3, cma_typ);
	set_typ *set; NEW(set,rhpt->NumSets() + 3, set_typ);
	sst_typ *osst=0;
	Int4	set_size=mcs->GetSetSize();
	for(i = 1; i < rhpt->NumSets(); i++){ set[i]=MakeSet(set_size); } set[i]=MakeSet(set_size);
	for(i = 1; i <= SizeTrueMainCMA; i++){
		s=random_integer(rhpt->NumSets()-1); s++;
		AddSet(i,set[s]);
	} s=rhpt->NumSets();
	for(; i <= SizeMainCMA; i++){ AddSet(i,set[s]); }
	sma[1]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA("Set1",sma[1]);
	for(i=2; i < rhpt->NumSets(); i++){
		cma_typ tcma=GetInSetCMSA(set[i],TrueMainCMA);
		sprintf(str,"Set%d",i); rhpt->ReNameSet(i,str); 
		sma[i]=MakeConsensusCMSA(tcma); RenameCMSA(str,sma[i]);
		osst=SST_FromSeq(TrueSeqCMSA(1,sma[i]));
		char *tmp[3]; tmp[0]=SST2ArgStrAlpha(osst,mcs->RtnLengthMainCMSA(),AB);
		// fprintf(stderr,"%d ('Set%d'): %s\n",i,i,tmp[0]); // debug...
		// rhpt->ReSetArgv(i,1,tmp); 
		free(tmp[0]); free(osst); TotalNilCMSA(tcma);
	} // rhpt->Put(stderr);
	// 4. Create new mcs.
	this->SetDefaultArguments( );
	mcs_typ *rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,rhpt->NumSets(),set,rhpt,sma,Argc,Argv);
	rhpt=rtn_mcs->GetHpt( );
	rtn_mcs->NoFailureMode=TRUE; rtn_mcs->SaveBest=TRUE; 
	rtn_mcs->PutHyperPartition(stderr); 
#if 0
	   for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,75);
		rtn_mcs->SetMaxNumColumns(i,50);
	   }
#endif
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->Sample(2,2,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->SampleColumns(TRUE);
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->LoadUpBestColumns(25);
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->Sample(2,2,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	rtn_mcs->PutHyperPartition(stderr); 
#if 0
	// rtn_mcs->UpdateCSQ();
	rtn_mcs->PutHyperPartition(stderr); 
	for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,50);
		rtn_mcs->SetMaxNumColumns(i,30);
	}
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->Sample(1,1,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->LoadUpBestColumns(25);
	for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,25);
		rtn_mcs->SetMaxNumColumns(i,20);
	}
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->Sample(1,1,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->SampleColumns(TRUE);
	rtn_mcs->PutHyperPartition(stderr); 
	rtn_mcs->Sample(2,2,2,2);
	rtn_mcs->PutHyperPartition(stderr); 
	// delete rhpt;
#endif
	sprintf(str,"%s_rand",infile); 
	FILE *fp=open_file(str,".out","w"); rtn_mcs->PutHyperPartition(fp); fclose(fp); 
	DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( ); PutCheckPoint(str,FALSE);
	// delete rtn_mcs;	// don't need to call DeleteMCS(); not made by CreateMCS();
}


BooLean omc_typ::CheckSSTvsCSQ(mcs_typ *xmcs)
// Make sure that the patterns are consistent with the consensus seqs.
{
	Int4    n,i,k,len=LengthCMSA(1,TrueMainCMA);
	unsigned char r;
	BooLean	okay=TRUE;
        hpt_typ *hpt=xmcs->GetHpt( );
        for(n=0,i=1; i < hpt->NumBPPS(); i++){
            e_type xE=xmcs->RtnKeySeq(i);
	    assert(LenSeq(xE)==len);
            sst_typ *xsst=xmcs->RtnCopyOfSST(i);
            for(k=1; k <= len; k++){
		if(xsst[k]==0) continue;
		r=ResSeq(k,xE);
                if(!MemSset(r,xsst[k])){
		   xmcs->PutHyperPartition(stderr);
		   xmcs->PutPttrnVsConsSeq(stderr,"omc->Sample() debug 0");
		   PutSeq(stderr,xE,AB);
		   if(xmcs->IsBestRestored()) fprintf(stderr,"best is restored\n");
                   fprintf(stderr,"%d.",i); PutSST(stderr,xsst[k],AB);
		   fprintf(stderr,"%d != '%c'\n",k,AlphaChar(r,AB));
		   assert(MemSset(r,xsst[k]));
		   okay=FALSE;
                }
            } free(xsst); 
         } return okay;
}

BooLean omc_typ::CheckBestCopy(mcs_typ *xmcs)
// Make sure that the copies are the same.
{
	double	R,d=best_mcs->CalcTotalLPR(),D;
	// const double tolerance=1e-6; // 1 part per million difference.
	const double tolerance=1e-4; // 1 part per 10,000.
        Int4    n,i,k,row,len=LengthCMSA(1,TrueMainCMA);
	BooLean	okay=CheckSSTvsCSQ(xmcs);
	if(!CheckSSTvsCSQ(best_mcs)) okay=FALSE;
	if(Best_LPR != d){
	   // check for rounding errors...
	   if(Best_LPR > d) R=(Best_LPR/d)-1.0; else R=(d/Best_LPR)-1.0; 
	   if(R > tolerance) okay=FALSE;
           hpt_typ *hpt=xmcs->GetHpt( );
           hpt_typ *bhpt=best_mcs->GetHpt( );
	   if(!bhpt->TheSame(hpt)){ fprintf(stderr,"WARNING: hpt objects differ!!\n"); okay=FALSE; }
           fprintf(stderr,"Best_LPR(%.2f) != best_mcs->LLR(%.2f); diff = %.2g%c\n",Best_LPR,d,R*100,'%');
           for(n=0,i=1; i <= hpt->NumBPPS(); i++){
	       D=best_mcs->RtnSubLPR(i); d=xmcs->RtnSubLPR(i);
	       if(D!=d){
		 if(D > d) R=(D/d)-1.0; else R=(d/D)-1.0;
		 if(R > tolerance){	
		    fprintf(stderr,"=========== SubLpr (%d) = %.3f != %.3f ============\n", i,D,d);
		    okay=FALSE; xmcs->PutPttrnLLR(stderr,i); best_mcs->PutPttrnLLR(stderr,i);
		 }
	       }
               for(row=1; row < hpt->NumSets(); row++){
	             if(xmcs->RtnContribLLR(row,i,d) && best_mcs->RtnContribLLR(row,i,D)){
			if(d != D){
				fprintf(stderr,"row %d: (%d) %.2f != %.2f\n",row,i,d,D); okay=FALSE;
			}
		     }
	       }
               e_type xE=xmcs->RtnKeySeq(i),bE=best_mcs->RtnKeySeq(i);
               if(!IdentSeqs(xE,bE)) AlnSeqSW(stderr,11,1,xE,bE,AB);
               sst_typ *xsst=xmcs->RtnCopyOfSST(i);
               sst_typ *bsst=best_mcs->RtnCopyOfSST(i);
               for(k=1; k <= LengthCMSA(1,TrueMainCMA); k++){
                        if(xsst[k] != bsst[k]){
                           fprintf(stderr,"%d.",i);
                           PutSST(stderr,xsst[k],AB); fprintf(stderr,"%d != ",k);
                           PutSST(stderr,bsst[k],AB); fprintf(stderr,"%d\n",k);
			   okay=FALSE; n++;
                        }
               } free(xsst); free(bsst);
           } // assert(Best_LPR == best_mcs->CalcTotalLPR( ));
         } return okay;
}

Int4    omc_typ::RecordChange(char *title)
{
        if(logfp){
          if(status > 0){
               fprintf(logfp,"%s: %d - %d = %d (%d seconds; %.1f nats; %d/%d failed; best=%.1f);",
                     title,status,t_stat,status - t_stat,time(NULL)-log_time,
                     this->CalcLLR(mcs),mcs->NumFailedNodes(),Hpt->NumBPPS(),Best_LPR);
               Hpt->PutAsTree(logfp);
          } else {
               fprintf(logfp,"%s: no change (%d seconds; %.1f nats; %d/%d failed; best=%.1f)\n",
                     title,time(NULL)-log_time,
                     this->CalcLLR(mcs),mcs->NumFailedNodes(),Hpt->NumBPPS(),Best_LPR);
          } fflush(logfp);
        } log_time=time(NULL);
        // if(status > 0) return 1; else return 0;
        return (status-t_stat);
}

void	omc_typ::CompareInput(FILE *fp,set_typ *set,hpt_typ *hpt,cma_typ *SMA)
{
	Int4	i,j,id;
#if 0
	FILE *fp=open_file(outfilename,".hpt","w"); hpt->Put(fp);  fclose(fp); fp=0;
	fp=open_file(outfilename,".sma","w");
	for(i=1; i <= hpt->NumBPPS(); i++) if(in_sma[i]) PutCMSA(fp,in_sma[i]); fclose(fp); fp=0;
#endif
	for(i=1; i <= hpt->NumSets(); i++){
	    if(i == hpt->NumSets()) { 
	    	fprintf(fp,"%d. %s(%d):",i,hpt->ElmntSetName(i),
			CardInterSetINotJ(set[i],RandomSet));
		fprintf(fp," (not present).\n"); continue; 
	    } fprintf(fp,"%d. %s(%d):",i,hpt->ElmntSetName(i),CardSet(set[i]));
	    id=hpt->ItoSetID(i);
	    if(SMA[i]){ fprintf(fp," %s.\n",NameCMSA(SMA[i]));
	    } else { fprintf(fp," (null).\n"); assert(SMA[i]); }
	}
	mcs->PutHyperPartition(stderr);
}

BooLean	omc_typ::CheckID4NewMCS(mcs_typ *xmcs, Int4 ID)
{
	hpt_typ *hpt=xmcs->GetHpt();
        Int4	i=hpt->SetIDtoI(ID);
        if(i < 1 || i >= hpt->NumSets()){
             fprintf(stderr,"ID = %d; i = %d; numsets=%d (new = %d)\n",ID,i,hpt->NumSets());
             hpt->Put(stderr,FALSE); assert(i > 0 && i < hpt->NumSets()); 
	     return FALSE;
	} return TRUE;
}

BooLean	omc_typ::TheSame(mcs_typ *xmcs)
{
	hpt_typ *hpt=xmcs->GetHpt();
	assert(Hpt->TheSame(hpt));
}

BooLean	omc_typ::OverlappingLineages(Int4 I, Int4 J)
{
	Int4	i,j,*P;
	BooLean	rtn=FALSE;
	assert(I > 0 && I < Hpt->NumSets());	// disallow Reject node.
	assert(J > 0 && J < Hpt->NumSets());	// disallow Reject node.
	Hpt->IsTree(P);
	set_typ	setI=MakeSet(Hpt->NumSets()+3); ClearSet(setI);
	set_typ	setJ=MakeSet(Hpt->NumSets()+3); ClearSet(setJ);
	for(i=I; P[i] != 0; i = P[i]){ AddSet(i,setI); }
	for(j=J; P[j] != 0; j = P[j]){ AddSet(j,setJ); }
	if(CardInterSet(setI,setJ) > 0) rtn=TRUE;
	NilSet(setI); NilSet(setJ); free(P);
	return rtn;
}

void	omc_typ::PutLineage(FILE *fp,Int4 node)
// move to hpt_typ evantually
{
	Int4	i,j,*P;
	char	*name;
	assert(node > 0 && node < Hpt->NumSets());	// disallow Reject node.
	Hpt->IsTree(P);
	for(i=node; P[i] != 0; i = P[i]){
		name=Hpt->ElmntSetName(i);	// checks that i is within range.
		fprintf(fp," %d (\"%s\") ",i,name);
		if(P[i] != 0) fprintf(fp," -> ");
	} fprintf(fp,"root\n"); free(P);
}

void	omc_typ::PutPatternOverlap(FILE *fp)
{
	Int4 i,j,m,n,s,x,len=mcs->RtnLengthMainCMSA(),*P;
	sst_typ *isst,*jsst,*tsst,**xsst=mcs->RtnCopyOfSSTs( );
	char	*name,*name2;
	Hpt->IsTree(P);

	h_type HG=Histogram("pattern matches between nodes",0,50,1.0);
	h_type sHG=Histogram("pattern subsets between nodes",0,50,1.0);
	for(i=1; i < Hpt->NumBPPS(); i++){
	    isst=xsst[i]; 
	    for(j=i+1; j <= Hpt->NumBPPS(); j++){
	 	jsst=xsst[j]; 
	        for(n=0,m=0,s=1; s <= len; s++){
		   BooLean overlap=FALSE;
		   if(isst[s] == 0 || jsst[s] == 0) continue;
		   // fprintf(fp,"{"); PutSST(fp,jsst[s],AB); fprintf(fp,"}\n");
		   if(isst[s] == jsst[s]){ n++; overlap=TRUE; }
		   else if(SubSset(isst[s],jsst[s])){ m++; overlap=TRUE; }
		   else if(SubSset(jsst[s],isst[s])){ m++; overlap=TRUE; }
		   // else if(IntersectSset(isst[s],jsst[s]) != 0){ is++; overlap=TRUE; }
		   if(overlap){
			fprintf(fp,"Node %d (\"%s\") ",i,Hpt->ElmntSetName(i));
			PutSST(fp,isst[s],AB);	
			fprintf(fp,"%d = node %d (\"%s\") ",s,j,Hpt->ElmntSetName(j));
			PutSST(fp,jsst[s],AB);	fprintf(fp,"%d\n  ",s); 
			char relFG=mcs->RtnRelateFG(i,j),relBG=mcs->RtnRelateBG(i,j);
			fprintf(fp,"  FGs: %c; BGs: %c\n  ",relFG,relBG); 
			PutLineage(fp,i); fprintf(fp,"  "); PutLineage(fp,j); 
			if(this->OverlappingLineages(i,j)) fprintf(fp," lineages overlap!\n"); else fprintf(fp,"\n");
		   }
	        } if(n > 0) IncdHist(n,HG); if(m > 0) IncdHist(m,sHG);
	        if(n >= 5 || m >= 25){
		  name=Hpt->ElmntSetName(i); name2=Hpt->ElmntSetName(j);
		  fprintf(fp,"Node %d (\"%s\") & node %d (\"%s\") = %d (%d sub)\n",
			i,name,j,name2,n,m);
	  	  char *pttrn=SST2ArgStrAlpha(isst,len,AB);
		  fprintf(fp,"%d: %s\n",i,pttrn); free(pttrn);
	  	  pttrn=SST2ArgStrAlpha(jsst,len,AB);
		  fprintf(fp,"%d: %s\n",j,pttrn);
		  NEW(tsst,len+5,sst_typ); 
	          for(s=1; s <= len; s++){
			if(isst[s] == 0) continue;
			if(isst[s] == jsst[s]) tsst[s]=isst[s];
		  }
	  	  pttrn=SST2ArgStrAlpha(tsst,len,AB);
		  fprintf(fp," common: %s\n\n",pttrn); free(tsst);
	        }
	    }
	}
	PutHist(fp,60,HG); NilHist(HG); 
	PutHist(fp,60,sHG); NilHist(sHG); 

	// mcs->PutSetRelations(stderr);
	//*************** check subset constraints **************
	set_typ subtree=MakeSet(Hpt->NumBPPS()+2); // get set for node's subtree.
	for(s=1; s <= len; s++){
	        for(i=1; i < Hpt->NumBPPS(); i++){
	 	  isst=xsst[i]; 
		  if(isst[s] == 0) continue;
		  ReSetSubTree(subtree, i, Hpt);
	          for(n=0,m=0,j=1; j <= Hpt->NumBPPS(); j++){
		     if(i == j) continue;
		     if(!MemberSet(j,subtree)) continue;
	 	     jsst=xsst[j]; 
		     if(jsst[s] == 0) continue;
		     name=Hpt->ElmntSetName(i); name2=Hpt->ElmntSetName(j);
		     if(isst[s] == jsst[s]){
			fprintf(fp,"ERROR: node %d (\"%s\") ",i,name);
			PutSST(fp,isst[s],AB); fprintf(fp,"%d == ",s);
			PutSST(fp,jsst[s],AB); fprintf(fp,"%d",s);
			fprintf(fp," node %d (\"%s\")!\n\n",j,name2);
		     } else if(SubSset(jsst[s],isst[s])){
			continue;
			fprintf(fp,"OKAY: node %d (\"%s\") ",i,name);
			PutSST(fp,isst[s],AB); fprintf(fp,"%d > ",s);
			PutSST(fp,jsst[s],AB); fprintf(fp,"%d",s);
			fprintf(fp," node %d (\"%s\")!\n\n",j,name2);
		     } else if(SubSset(isst[s],jsst[s])){	// is i < j?
			fprintf(fp,"ERROR: node %d (\"%s\") ",i,name);
			PutSST(fp,isst[s],AB); fprintf(fp,"%d < ",s);
			PutSST(fp,jsst[s],AB); fprintf(fp,"%d",s);
			fprintf(fp," node %d (\"%s\")!\n   (node)",j,name2);
			for(x=j; x != i && x > 0; x = P[x]){ 
				fprintf(fp," %d -> ",x);
			} fprintf(fp,"%d (parent)\n",i);
		     }
		  } 
	   }
	}
	for(i=1; i <= Hpt->NumBPPS(); i++) free(xsst[i]); free(xsst);
	free(P); NilSet(subtree);
}

#if 0
//================ for developing & testing 'B' option ====================

double  LogCumHypGeomProb(double N1,double N2, double n,double x)
#if 0   //****************************************************
N total balls with N1 red balls and N2 = N-N1 black balls.
Choose n balls at random.  The probability that the
group so chosen will contain x or more red balls is
given by: p=CumHypGeomProb(N1,N2,n,x).
#endif  //****************************************************
{
        double  p,K,end;

        if(x == 0) return 1.0;
        end = MINIMUM(double,N1,n);
        K = (lngamma(N1+1)+lngamma(N2+1)-lngamma(N2+N1+1)+lngamma(n+1)+lngamma(N2+N1-n+1));
        for(p=0.0; x <= end; x += 1.0){
           p += K-lngamma(x+1)-lngamma(N1-x+1)-lngamma(n-x+1)-lngamma(N2-n+x+1);
        } return p;
}

double	***BoltzScore;	// BoltzScore[n][j][r];
double	***BoltzStd;	// BoltzStd[n][j][r];
double	**Entropy;	// Entropy[n][j];

double	bpps_typ::PutBoltzmannLike(FILE *fp,double *StdFrq, Int4 j, UInt4 **CntBG,UInt4 **CntFG, Int4 n)
// for saving information on pattern residues and 
// read this in using csp->ReadBinomialBPPS(-cutoff,res_evals);
{
	sst_typ QstJ=qst[j];
	unsigned char QueryJ=query[j];
	Int4 i,r;
	double N,M,p,q,n1,n2,m1,m2,P,StdP;
//	fprintf(fp,"#             __Foreground___    __Background___     ____Information____    WtNumSeq\n");
//	fprintf(fp,"#   Pattern:  Match  Diverge     Match  Diverge \n");
	 
	double	map=0.0,map0=0.0,RefSubJ,StdRefSubJ,StdSubJ,d,D,DD,dd,RefMap;
	{ qst[j]=0; query[j]=0; }
	SubLPR(CntBG,CntFG); map0=subLPR[0];

#if 1
	UInt4 **RefCntFG,TotalCntFG,TotalCntBG,**StdCntFG,**StdCntBG;
	NEWP(RefCntFG,k+3,UInt4); NEWP(StdCntFG,k+3,UInt4); NEWP(StdCntBG,k+3,UInt4);
	for(TotalCntFG=0,TotalCntBG=0,r=0; r <= nAlpha(AB); r++)
		{ TotalCntFG += CntFG[j][r]; TotalCntBG += CntBG[j][r]; }
	for(i=1; i <= k; i++){
	    if(i != j){
		StdCntBG[i]=CntBG[i]; RefCntFG[i]=CntFG[i]; StdCntFG[i]=CntFG[i]; continue; 
	    }
	    NEW(RefCntFG[j],nAlpha(AB)+3,UInt4);
	    NEW(StdCntFG[j],nAlpha(AB)+3,UInt4); NEW(StdCntBG[j],nAlpha(AB)+3,UInt4);
	    for(r=0; r <= nAlpha(AB); r++){

		// m1=(double) CntBG[j][r]*wt_factor; M=(double) TotalCntBG*wt_factor;
		m1=(double) CntBG[j][r]; M=(double) TotalCntBG;
		d = m1/M;	// = the mean of beta distribution.
		D = floor( (d*(double)TotalCntFG) + 0.5); 
		if(D < 0.0) D = 0.0; RefCntFG[j][r]= (UInt4) D;

		d = StdFrq[r];		// = the overall frequency in the input alignment.
		D = floor( (d*(double)TotalCntFG) + 0.5); 
		if(D < 0.0) D = 0.0; StdCntFG[j][r]= (UInt4) D;

		D = floor( (d*(double)TotalCntBG) + 0.5); 
		if(D < 0.0) D = 0.0; StdCntBG[j][r]= (UInt4) D;

		// fprintf(fp,"%c=%.1f; ",AlphaChar(r,AB),wt_factor*(double)RefCntFG[j][r]);
	    } // fprintf(fp,"\n  Ref=%d; TotalFG=%d\n",C,TotalCntFG);
	}
#endif
	double	max=0,min=999999,entropy=0;
	double	maxFrq=0,minFrq=0;
	unsigned char maxR=0,minR=0;
	char	PttrnPos[3]; PttrnPos[0]=0;
	// for(unsigned char R=0; R <= nAlpha(AB); R++)
	for(unsigned char R=1; R <= nAlpha(AB); R++)
	{
	   if(R==0){ qst[j] = QstJ; query[j]=QueryJ; }
	   else { qst[j] = SsetLet(R); query[j]=R; } 

#if 1
	   SubLPR(StdCntBG,CntFG); RefMap=subLPR[0]; StdSubJ=subLPR[j];
	   SubLPR(StdCntBG,StdCntFG); RefMap=subLPR[0]; StdRefSubJ=subLPR[j];
	   SubLPR(CntBG,RefCntFG); RefMap=subLPR[0]; RefSubJ=subLPR[j];
#endif
	   SubLPR(CntBG,CntFG); map=subLPR[0];

	   m1 = MatchBG[j]; m2 = MisMatchBG[j]; M=m1+m2; p=(m1+1)/(M+2);
	   n1=MatchFG[j];  n2=MisMatchFG[j]; N=n1+n2; q=(n1+1)/(N+2);
	   // if(n1 <= 0.0) continue;
	   p = n1/N;
	   if(n1 > 0.0) entropy -= p*log(p); else n1=0.001;

#if 1
	   StdP=StdSubJ - StdRefSubJ; StdP=10*StdP/N;
	   P = subLPR[j] - RefSubJ; P=10*P/N;
#elif 1
	   P = subLPR[j]/(wt_factor*N);
#else
	   P = (map-map0)/(wt_factor*N);
#endif
#if 0
	   // P = LnCBP(n1,N,q) - LnCBP(n1,N,p); P = 100*P/N;
	   double Red=m1+n1,Black=m2+n2,mode,mean;
	   mean = (N+M)*Red/(Red+Black);
	   mode = n1*(N+M+1)*(Red+1)/(Red+Black+1);
	   dd=LogCumHypGeomProb(Red,Black,N,n1) - LogCumHypGeomProb(Red,Black,M,mode);
	   // dd=map-RefMap;
	   DD=-log((n1/N)/(m1/M));
#endif
#if 1
	   if(qst[j]){
	     char tmp[30],*ptr; ptr=tmp;
	     for(r=StartAlpha; r <= nAlpha(AB); r++){
		if(MemSset(r,qst[j])){ sprintf(ptr,"%c",AlphaChar(r,AB)); ptr++; }
	     }
	     if(0 || R == 0){
		PttrnPos[0]='!'; PttrnPos[1]=0;
	       sprintf(ptr,"%d",j);
	       fprintf(fp,"%12s: %5.0f %6.0f (%2.0f%c) %5.0f %6.0f (%2.0f%c)    %5.1f  %5.1f   %5.1f",
		  tmp,n1,n2, 100.0*n1/N,'%',m1,m2, 100.0*m1/M,'%', subLPR[j],P,StdP); 
	       if(R > 0) fprintf(fp,"      %5.1f\n", M+N);
	       else fprintf(fp,"      %5.1f\n---------\n", M+N);
	     } else {
	       fprintf(fp,"%s %.0f %.3f %.1f %.2f %.3f\n",tmp,n1,100*n1/N,-P,-StdP,100*m1/M); 
		BoltzScore[n][j][R]=P; BoltzStd[n][j][R]=StdP; 
	     }
	   }
	   // D=m1/M; d=n1/N;
	   if(!isfinite(P)) continue;
	   if(R > 0){
		//BoltzScore[n][j][R]=P; BoltzStd[n][j][R]=StdP; 
		if(-P > max){ max = -P;  maxR=R; maxFrq=100*m1/M; }
		if(-P < min){ min = -P;  minR=R; minFrq=100*n1/N; }
	   }
#else
	   if(!isfinite(P)) continue;
	   else if(R > 0){ 
		if(P > max){ max=P;  maxR=R; minR=R; maxFrq=n1/N; minFrq=m1/M; }
		// if(P > max){ max = P;  maxR=R; maxFrq=100*n1/N; }
		// if(P < min){ min = P;  minR=R; minFrq=100*n1/N; }
	   }
#endif
	}
	fprintf(fp,"Range(%d)%s=%.1f: %c%d(%.2f) --> %c%d(%.2f); entropy=%.3f\n\n",
		n,PttrnPos,max,AlphaChar(maxR,AB),j,maxFrq,AlphaChar(minR,AB),j,minFrq,entropy);
	Entropy[n][j]=entropy;
	qst[j] = QstJ; query[j]=QueryJ;
#if 0
	free(RefCntFG[j]); free(RefCntFG); 
	free(StdCntFG[j]); free(StdCntFG); free(StdCntBG[j]); free(StdCntBG); 
	return max-min;
#else
	// return max;
	return max*maxFrq;
#endif
}

double	che_typ::PutBoltzmannLike(FILE *fptr,double *frq,Int4 c, Int4 n)
{
        return pps->PutBoltzmannLike(fptr,frq,c,ResWt[2],ResWt[1],n);
}

double	mcs_typ::PutBoltzmannLike(FILE *fp,Int4 n,Int4 c)
{
	fprintf(stdout,"====== %d: %s (column %d) ========\n",n,Hpt->ElmntSetName(n),c);
	double  *freq = tFreqSeqSet(TrueDataCMSA(TrueMainCMA));
	// double  *freq = tFreqSeqSet(DataCMSA(TrueMainCMA));
        return che[n]->PutBoltzmannLike(fp,freq,c,n); 
}

double  GetSqProbCMSA(cma_typ cma, Int4 n, cma_typ CMA)
{
    assert(TotalLenCMSA(cma) == TotalLenCMSA(CMA));
    assert(nBlksCMSA(cma) == 1); assert(nBlksCMSA(CMA) == 1);
    ss_type             data=DataCMSA(CMA);
    st_type             sites = SitesCMSA(CMA);
    fm_type             *model=ModelsCMSA(cma);
    Int4                s,end;
    double              prob;
    unsigned char       *seq;

    if(n < 1 || n > NSeqsSeqSet(data)) print_error("GetSqProbCMSA( ) input error");
    end = SqLenSeqSet(n,data) + 1;
    end += 1 - LenFModel(model[1]);
    seq = SeqSeqSet(n,data);
    s = SitePos(1,n,1,sites);
    prob = (double) LikelihoodFModel(seq, s, model[1]);
    if(prob > 0.0) prob=log10(prob); else assert(prob > 0.0);
    return prob;
}
#endif

void	omc_typ::FindKeyPositions( )
//****************** Compare BPPS with BILD scores. *******************
// This version uses Altschul's approach: sum BILD - total BILD.
{
	hpt_typ *Hpt=mcs->GetHpt();
	char	dms_mode='F',wt_factor=100;	// 
	Int4	pernats=1000;
	dms_typ	*dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	che_typ **che=mcs->RtnChe( );	// need to create this function.
        mcs->PutHyperPartition(stdout); 
	UInt4  **RootWtCntsFG=che[1]->GetResWtsFG();
	double	dd,d,Dd,FG,BG,DD,dD,sumFG,total_d;
	double	max=-INT4_MAX,min=INT4_MAX; 
	double	Max=(double)-INT4_MAX,Min=(double)INT4_MAX; 
	Int4	i,n,x,*Parent; Hpt->IsTree(Parent);
	dh_type dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	dh_type ave_dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	UInt4	WtCntsNonRoot[50],WtCntsRootOnly[50];
	for(i=1; i <= mcs->RtnLengthMainCMSA(); i++){
	    double RootFG=dms->bild(RootWtCntsFG[i]);
	    fprintf(stdout,"\nRootFG[%d]=%.1f: ",i,RootFG); 
	    for(x=0; x <= nAlpha(AB); x++) WtCntsNonRoot[x]=0;
	    for(sumFG=FG=0.0,n=2; n < Hpt->NumSets(); n++){
	        if(Parent[n] != 1) continue; 
	        UInt4  **WtCntsFG=che[n]->GetResWtsFG();
		for(x=0; x <= nAlpha(AB); x++) WtCntsNonRoot[x] += WtCntsFG[i][x];
		dd=dms->bild(WtCntsFG[i]); sumFG += dd;
	        fprintf(stdout,"%d.%s=%.1f; ",n,Hpt->SetName(n),dd);
	    }
#if 0	    // sumFG computed at this point; need to add Root Only to sumFG.
	    for(x=0; x <= nAlpha(AB); x++) WtCntsRootOnly[x] = RootWtCntsFG[i][x] - WtCntsNonRoot[x];
	    sumFG += dms->bild(WtCntsRootOnly); 
	    DD = sumFG - RootFG;
#else	// ignore the root node.
	    FG = dms->bild(WtCntsNonRoot); 
	    DD = sumFG - FG;
#endif
	    insrtHeap(i,-(keytyp)DD,dH);
	    if(DD > Max) Max=DD;
	    if(DD < Min) Min=DD;
	} delete dms; fprintf(stdout,"\n\n");
	h_type	HG=Histogram("delta BILD scores",-50,100,5);
	while(!emptyHeap(dH)){
		dd=-minkeyHeap(dH); assert((i=delminHeap(dH)) != 0);
		Dd = 100*(dd - Min)/(Max - Min); IncdHist(Dd, HG);
		fprintf(stdout,"column %d: %.1f average dF (%.1f)\n",i,dd,Dd);
	} fprintf(stdout,"\n"); Nildheap(dH);
	fprintf(stdout,"Max=%.3f; Min=%.3f\n",Max,Min); 
	PutHist(stdout,60,HG); NilHist(HG);
}

void	omc_typ::FindKeyPositions2( )
//****************** Compare BPPS with BILD scores. *******************
{
	hpt_typ *Hpt=mcs->GetHpt();
	char	dms_mode='F',wt_factor=100;	// 
	Int4	pernats=1000;
	dms_typ	*dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	che_typ **che=mcs->RtnChe( );	// need to create this function.
	UInt4	*WtCntsFB; NEW(WtCntsFB,nAlpha(AB)+3,UInt4);
        mcs->PutHyperPartition(stdout); 
	UInt4  **RootWtCntsFG=che[1]->GetResWtsFG();
	double	dd,d,Dd,FG,BG,DD,dD,RootFG,aveFG;
	double	max=-INT4_MAX,min=INT4_MAX; 
	double	Max=-INT4_MAX,Min=INT4_MAX; 
	Int4	i,m,n,x,*Parent; Hpt->IsTree(Parent);
	dh_type dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	dh_type ave_dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
UInt4	WtCntsAll[50];
UInt4	WtCntsAllBut[50];
UInt4	WtCntsRootOnly[50];
	for(i=1; i <= mcs->RtnLengthMainCMSA(); i++){
	    RootFG=dms->bild(RootWtCntsFG[i]);
	    fprintf(stdout,"RootFG[%d]=%.1f: ",i,RootFG); 
// for(x=0; x <= nAlpha(AB); x++) WtCntsAll[x]=0;
for(x=0; x <= nAlpha(AB); x++) WtCntsAll[x]=0;
	    for(aveFG=FG=0.0,m=0,n=2; n < Hpt->NumSets(); n++){
	        if(Parent[n] != 1) continue;
	        UInt4  **WtCntsFG=che[n]->GetResWtsFG();
for(x=0; x <= nAlpha(AB); x++) WtCntsAll[x] += WtCntsFG[i][x];
		dd=dms->bild(WtCntsFG[i]);
	        fprintf(stdout,"%d.%s=%.1f; ",n,Hpt->SetName(n),dd);
		FG += dd; 

	        UInt4  **WtCntsBG=che[n]->GetResWtsBG();
		Dd=dms->bild(WtCntsBG[i]); 
#if 1
		aveFG += (dd + Dd) -RootFG; m++;
#else	//  sum bild scores only...
		aveFG += dd; m++;
#endif
	    } 
#if 0
	    for(x=0; x <= nAlpha(AB); x++) WtCntsRootOnly[x] = RootWtCntsFG[i][x] - WtCntsAll[x];
	    DD = FG + dms->bild(WtCntsRootOnly) - RootFG; 
#elif 1	// node FG + node BG - parent FG to compute delta-BILD score.
	    for(DD=0.0,n=2; n < Hpt->NumSets(); n++){
	        if(Parent[n] != 1) continue;
	        UInt4  **WtCntsFG=che[n]->GetResWtsFG();
	        // for(x=0; x <= nAlpha(AB); x++) WtCntsAllBut[x] = WtCntsAll[x]- WtCntsFG[i][x];
	        for(x=0; x <= nAlpha(AB); x++) WtCntsAllBut[x] = RootWtCntsFG[i][x]- WtCntsFG[i][x];
		dd=dms->bild(WtCntsFG[i]);
		BG=dms->bild(WtCntsAllBut);
		DD += dd + BG -RootFG;
	    }
#else
	    DD = FG - dms->bild(WtCntsAll); 
#endif
	    if(DD > max) max=DD;
	    if(DD < min) min=DD;
	    insrtHeap(i,-(keytyp)DD,dH);
	    fprintf(stdout," ==> %.1f\n",DD); 

#if 1
	    aveFG = aveFG/(double)m;
#else	// sum BILD scores only.
	    aveFG = aveFG - RootFG;
#endif
	    insrtHeap(i,-(keytyp)aveFG,ave_dH);
	    if(aveFG > Max) Max=aveFG;
	    if(aveFG < Min) Min=aveFG;
	} delete dms;
	h_type	HG1=Histogram("BILD scores",-50,100,5);
	while(!emptyHeap(dH)){
		dd=-minkeyHeap(dH); assert((i=delminHeap(dH)) != 0);
		Dd = 100*(dd - min)/(max - min); IncdHist(Dd, HG1);
		fprintf(stdout,"column %d: %.1f average dF (%.1f)\n",i,dd,Dd);
	} fprintf(stdout,"\n"); Nildheap(dH);
	fprintf(stdout,"max=%.3f; min=%.3f\n",max,min); 
	PutHist(stdout,60,HG1); NilHist(HG1);

	HG1=Histogram("aveBILD scores",-50,100,5);
	while(!emptyHeap(ave_dH)){
		dd=-minkeyHeap(ave_dH); assert((i=delminHeap(ave_dH)) != 0);
		Dd = 100*(dd - Min)/(Max - Min); IncdHist(Dd, HG1);
		fprintf(stdout,"column %d: %.1f average dF (%.1f)\n",i,dd,Dd);
	} fprintf(stdout,"\n"); Nildheap(ave_dH);
	fprintf(stdout,"Max=%.3f; Min=%.3f\n",Max,Min); 
	PutHist(stdout,60,HG1); NilHist(HG1);
}

Int4	omc_typ::TestSubRoutine(char Mode)
{
	Int4	id,i,j,n,N,s,iter=0,result=0,No,mode=0;
	double	best_lpr,last_lpr,d,D,dd,DD;
	mcs_typ *xmcs=0,*rmcs=0;

// mode=15;	
//  mode=18;	// ResurrectRejects( );	

	SetStringency(stringency);
	assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
	mcs->NoFailureMode=TRUE;  // if node configuration subLPR <= 0.0, then reject it.
	mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
	if(Evolve) mcs->DoEvolve(); else mcs->DoNotEvolve();
   switch (Mode){
     case 'K': { this->FindKeyPositions( ); return 0; } break;
     case 'k': {	// Find Cross-Conserved residue patterns.
		// this->FindCrossConserved(64); return 0; 
		this->FindCrossConserved(XC_line); return 0; 
	} break;
     case 'B': 
#if 0
	{ return this->TestHiHMM2( ); } break; // Boltzmann-like distributions delta-F scores.
#else
	{ return this->TestHiHMM( ); } break;
#endif
     case 'd': // debugging.
	{
	this->Sample(150,2,2,2,2,mcs,0.05);
        this->Sample(50,2,2,3,2,mcs,0.01);

	} break;
     case 'c': // test pattern sampling routines.
	{
	    TestSubRoutine('O');
	    mcs->SampleColumns( );
	    TestSubRoutine('O');
	} break;
     case 'O': // look for pattern overlap.
	{
	    mcs->DoEvolve();
	    mcs->SampleColumns( );
	    mcs->PutHyperPartition(stderr);
	    this->PutPatternOverlap(stderr);
	} break;
     case 'G': // Test GrowLeaves() routine.
	if(0) mcs->Sample(2,2,2,1); else mcs->SampleColumns( ); 
	mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	mcs->PutHyperPartition(stderr); 
#if 1
	// this->PutStringency( );
	this->AddLeaves(300,'a');
	// GrowLeaves(stderr,300,iter);
	// this->AddLeaves(300,'a');
#else
	GrowLeaves(stderr,300,iter);
#endif
	// InsertInternalNodes(stderr,300);
	mcs->SaveBest=TRUE; mcs->StoreBest();
	mcs->PutHyperPartition(stderr); 
	sprintf(str,"%s_G",infile); PutCheckPoint(str,FALSE);
	this->Put(TRUE);  // input fp not used!
      break;
     case 'A': // Test AddLeaf() routine.
	{
	  if(Hpt->NumSets() == 2) SproutLeaves(stderr,300,1);
	  mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	  // return result;
	  // xmcs = this->AddLeaf(3,12);
	  // this->SampleMCS(xmcs,12);
	  GrowLeaves(stderr,300,1);
	  mcs->PutHyperPartition(stderr); 
	  sprintf(str,"%s_A",infile); PutCheckPoint(str,FALSE);
	  FILE *fp=open_file(infile,".cntrb","w"); mcs->PutMapContributions(fp); fclose(fp);
	  mcs->PutLpr(infile);	// creates an <infile>.lpr output file.
	  // this->Put();
	  result=1;
      	} break;
     case 'D': // Remove (delete) internal nodes.
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	mcs->PutMapContributions(stderr);
	RmInternalNodes(stderr,300, 1);
	mcs->PutMapContributions(stderr);
      break;
     case 'S': // create a simulated hierarchical alignment.
	{	// have it mirror the hierarchy found...
	   // 1. Generate a random alignment using the TrueMainCMA overall model.
		// cma_typ GetInSetCMSA(set_typ set, cma_typ cma);
		mcs->DoEvolve();
		mcs->SampleColumns( ); 
	        mcs->Sample(2,2,2,2);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
		mcs->RestoreBest(); mcs->CalcTotalLPR();
	  	FILE *fp=open_file(infile,"_sim_source.out","w"); mcs->PutHyperPartition(fp); fclose(fp);
	  	sprintf(str,"%s_sim_source",infile); PutCheckPoint(str,FALSE);  // FALSE --> don't rename.
		PutSimulatedAln( ); 
	   // 2. 
	} break;
     case 'Z': // Create an optimized flat hierarchy as a poor starting point.
	GrowLeaves(stderr,300, 1, 'z');
	mcs->SampleColumns( ); 
	mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	mcs->PutHyperPartition(stderr); 
	GrowLeaves(stderr,300,2,'Z');
	// InsertInternalNodes(stderr,300);
	mcs->SaveBest=TRUE; mcs->StoreBest();
	mcs->PutHyperPartition(stderr); 
	sprintf(str,"%s_Z",infile); PutCheckPoint(str,FALSE);
	this->Put(TRUE);  // input fp not used!
      break;
     case 'z': // create and optimize an arbitrary hierarchy using real sequences.
     case 'a': // create and optimize an arbitrary hierarchy using real sequences.
	{
	   // 1. create a random hpt.
	   Int4 NumRandOperations=30;
	   Hpt->Put(stderr); // hpt->PutRandomize(stdout,Random);
           hpt_typ *rhpt=0;
	   fprintf(stderr,"Mode = '%c'\n",Mode);
	   if(Mode == 'a') rhpt=Hpt->Randomize(NumRandOperations);
	   else if(Mode == 'z') rhpt=Hpt->RandomFlat(NumRandOperations);
	   else print_error("This should not happen (omc_debug.cc)\n");
           rhpt->Put(stderr);
           // rhpt->PutSorted(stderr);
	   // hpt_typ *shpt=rhpt->Sort( );
	   // 2. rename sets and create sma files.
	   cma_typ *sma; NEW(sma,rhpt->NumSets() + 3, cma_typ);
	   set_typ *set; NEW(set,rhpt->NumSets() + 3, set_typ);
	   sst_typ *osst=0;
	   Int4	s,set_size=mcs->GetSetSize();
	   for(i = 1; i < rhpt->NumSets(); i++){ set[i]=MakeSet(set_size); } set[i]=MakeSet(set_size);
	   for(i = 1; i <= SizeTrueMainCMA; i++){
		s=random_integer(rhpt->NumSets()-1); s++;
		AddSet(i,set[s]);
	   } s=rhpt->NumSets();
	   for(; i <= SizeMainCMA; i++){ AddSet(i,set[s]); }
	   sma[1]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA("Set1",sma[1]);
	   for(i=2; i < rhpt->NumSets(); i++){
		cma_typ tcma=GetInSetCMSA(set[i],TrueMainCMA);
		sprintf(str,"Set%d",i); rhpt->ReNameSet(i,str); 
		sma[i]=MakeConsensusCMSA(tcma); RenameCMSA(str,sma[i]);
		osst=SST_FromSeq(TrueSeqCMSA(1,sma[i]));
		char *tmp[3]; tmp[0]=SST2ArgStrAlpha(osst,mcs->RtnLengthMainCMSA(),AB);
		// fprintf(stderr,"%d ('Set%d'): %s\n",i,i,tmp[0]); // debug...
		// rhpt->ReSetArgv(i,1,tmp); 
		free(tmp[0]); free(osst); TotalNilCMSA(tcma);
	   } // rhpt->Put(stderr);
	   // 4. Create new mcs.
	   this->SetDefaultArguments( );
	   mcs_typ *rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,rhpt->NumSets(),set,rhpt,sma,Argc,Argv);
	   rhpt=rtn_mcs->GetHpt( );
	   rtn_mcs->NoFailureMode=TRUE; rtn_mcs->SaveBest=TRUE; 
	   rtn_mcs->PutHyperPartition(stderr); 
#if 0
	   for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,75);
		rtn_mcs->SetMaxNumColumns(i,50);
	   }
#endif
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->Sample(2,2,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->SampleColumns(TRUE);
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->LoadUpBestColumns(25);
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->Sample(2,2,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	   rtn_mcs->PutHyperPartition(stderr); 
 // exit(1);
	
#if 0
	   // rtn_mcs->UpdateCSQ();
	   rtn_mcs->PutHyperPartition(stderr); 
	   for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,50);
		rtn_mcs->SetMaxNumColumns(i,30);
	   }
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->Sample(1,1,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->LoadUpBestColumns(25);
	   for(i=2; i < rhpt->NumSets(); i++){
		rtn_mcs->RmWorstColumn(i,25);
		rtn_mcs->SetMaxNumColumns(i,20);
	   }
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->Sample(1,1,2,5);	// IterStart,IterEvolve,NumRounds,ColSampleStart,ModStart
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->SampleColumns(TRUE);
	   rtn_mcs->PutHyperPartition(stderr); 
	   rtn_mcs->Sample(2,2,2,2);
	   rtn_mcs->PutHyperPartition(stderr); 
	   // delete rhpt;
#endif
	   DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
	   sprintf(str,"%s_a",infile); PutCheckPoint(str,FALSE);
	   // delete rtn_mcs;	// don't need to call DeleteMCS(); not made by CreateMCS();
	} break;
     case 'f': // flatten the hierarchy in order to test MoveDown() routine.
 	{	
	mcs->PutHyperPartition(stderr);
	mcs_typ *rtn_mcs=FlattenHiearchy(); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
	mcs->PutHyperPartition(stderr);
	rtn_mcs->NoFailureMode=TRUE; 
	rtn_mcs->SaveBest=TRUE; 
	// mcs->Sample(2,2,2,2);
	mcs->PutHyperPartition(stderr);
	sprintf(str,"%s_f",infile); PutCheckPoint(str,FALSE);
	FILE *fp=open_file(infile,"_init.out","w"); mcs->PutHyperPartition(fp); fclose(fp); fp=0;
	}
      break;
     case 'r': 	// Test Root node residue environments ('R' mode for rst_typ):
	mcs->PutHyperPartition(stderr); 
	mcs->Sample(2,2,2,2);
	mcs->PutHyperPartition(stderr);
	result=1;
	break;
     case 'x': 	// Test ResurrectRejects( );
	mcs->PutHyperPartition(stderr); 
	ResurrectRejects(stderr);
	result=1;
	break;
     case 'b': // Test "storing the best_mcs" procedure.
	for(Int4 z=1; z <= 100; z++){
	   Sample(300,2,2,2,2,mcs,1000000);
	   this->CalcLLR(mcs);
	   RestoreFinalBest();
	   RandomizeSeqAssign(stderr,10,z); 
	} result=1;
      break;
     case 'C': // Test RtnContribLLR() routine.
	for(i=1; i < Hpt->NumSets(); i++){
	   for(j=i+1; j < Hpt->NumSets(); j++){
		if(mcs->RtnContribLLR(j,i,d)){
			fprintf(stderr,"%d: (%d) %.2f nats\n",j,i,d);
		}
	   }
	}
	result=1;
      break;
     case 'L': // Test LowerBranches() routine.
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	// mcs->Sample(2,2,2,1);
	TrimNodes(stderr,300,1);
	mcs->PutHyperPartition(stderr);
	LowerBranches(stderr,300, iter);
	// MoveNodesDown(300);
	mcs->PutHyperPartition(stderr);
	// mcs->SaveBest=TRUE; 
	// this->Put( );
	sprintf(str,"%s_L",infile); PutCheckPoint(str,FALSE);
	result=1;
      break;
     case 'R': // Test RaiseBranches() routine.
	{
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	// mcs->Sample(2,2,2,1);
#if 0	// For cd00516 testing...dir_cd00516/OMC_BPPS/Test/90C2OaA 
{
set_typ subtree=MakeSet(Hpt->NumBPPS()+2);
Int4 gp,node=2,sibling=8,*Parent; 
assert(Hpt->IsScrambledTree(Parent)); gp=Parent[node];
ReSetSubTree(subtree,node); mcs->MoveDown(gp,sibling,node,subtree); NilSet(subtree); free(Parent);
mcs_typ *rtn_mcs=this->Operate('u',0,0); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );

#if 0
subtree=MakeSet(Hpt->NumBPPS()+2);
node=6;sibling=8; assert(Hpt->IsScrambledTree(Parent)); gp=Parent[node];
ReSetSubTree(subtree,node); mcs->MoveDown(gp,sibling,node,subtree); NilSet(subtree); free(Parent);
rtn_mcs=this->Operate('u',0,0); DeleteMCS(mcs); mcs=rtn_mcs; Hpt=mcs->GetHpt( );
#endif

mcs->SampleColumns();
sprintf(str,"%sA",infile); PutCheckPoint(str,FALSE);
FILE *fp=open_file(infile,"_init.out","w"); mcs->PutHyperPartition(fp); fclose(fp); fp=0;
}
#endif
	mcs->PutHyperPartition(stderr);
	// RaiseBranches(stderr,300, iter);
	MoveNodesUp(300);
	mcs->PutHyperPartition(stderr);
	FILE *fp=open_file(infile,".cntrb","w"); mcs->PutMapContributions(fp); fclose(fp);
	sprintf(str,"%s_R",infile); PutCheckPoint(str,FALSE);
	// this->Put();
	result=1;
        } break;
     case 'P': // Test pattern sampling routines.
	{
	  Int4	n;
	  mcs->PutSetRelations(stderr);
	  mcs->PutHyperPartition(stderr);
	  this->PutPatternOverlap(stderr);
	  // mcs->PutPttrns(stderr);
	  PrintTime(stderr);
	  mcs->SampleColumns();
	  PrintTime(stderr);
	  mcs->PutHyperPartition(stderr);
	  // mcs->PutPttrns(stderr);
	  this->PutPatternOverlap(stderr);
	  FILE *fp=open_file(infile,"_tstP_ptn.lpr","w"); mcs->PutPttrnLLRs(fp); fclose(fp);
	  fp=open_file(infile,"_tstP.lpr","w"); mcs->PutLpr(fp); fclose(fp);
	  // this->Put(TRUE);
	  // mcs->PutRTF(FALSE);
	}
      break;
     case 'T': // Test FuseNodes() routine.
	mcs->SampleColumns();
	TrimNodes(stderr,300,1);
      break;
     case 'F':
	if(FocusNode == 0 && FocusSeq == 0){ // Test FuseNodes() routine.
	  for(j=2; j < Hpt->NumSets(); j++){
	    if(mcs->SetCard(j) == 0){
		mcs_typ *tmp_mcs=DeleteNode(Hpt->ItoSetID(j)); assert(tmp_mcs);
		DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	    }
	  }
	  mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	  // mcs->Sample(2,2,2,1);
	  mcs->PutHyperPartition(stderr);
	  // MergeSiblings(stderr,300.0, 1);
	  this->FuseNodes( );
	  mcs->PutHyperPartition(stderr);
	  result=1;
	} else {			// Test -Focus option
	   mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
// cerr << "Focused Search 2\n";
	   if(GrowFocusedLeaf(stderr,300,1)){
	     mcs->SaveBest=TRUE; mcs->StoreBest();
	     mcs->PutHyperPartition(stderr); 
	     sprintf(str,"%s_F",infile); PutCheckPoint(str,FALSE);
	     this->Put(TRUE);  // input fp not used!
	   } 
	}
      break;
     case 'I': 	// Test Insert Nodes
	status=this->InsertInternalNodes(stderr,300);
        this->RestoreBest(); mcs->PutHyperPartition();
	mcs->SaveBest=TRUE; mcs->StoreBest();
	// mcs->PutMapContributions(stderr,lpr);
	this->CalcLLR(mcs); this->RestoreFinalBest(); // mcs->PutHyperPartition();
	mcs->PutHyperPartition(stderr); 
	sprintf(str,"%s_I",infile); PutCheckPoint(str,FALSE);
	result=1;
      break;
     case 'i': 	// Test Insert Nodes
	{
	// mcs->SampleColumns( ); 
	mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	mcs->PutHyperPartition(stderr); this->CalcLLR(mcs);
	// mcs->SampleColumns( ); mcs->Sample(2,2,2,1); 
	mcs->SaveBest=TRUE; mcs->StoreBest();
	// TrimNodes(stderr,300,iter); 
	// mcs->SampleColumns( ); mcs->Sample(2,2,2,1); 
	this->CalcLLR(mcs); this->PrintTime(stderr);
	sprintf(str,"preliminary LLRs for three node combinations");
	preHG=Histogram(str,0,10000,400);
	sprintf(str,"current vs hybrid LLR differences");
	dfsHG=Histogram(str,0,10000,200);
	AddInternalNodes(stderr,300,iter); this->CalcLLR(mcs);
	mcs->PutMapContributions(stderr,lpr);
	PutHist(stderr,60,preHG); NilHist(preHG); preHG=0;
	PutHist(stderr,60,dfsHG); NilHist(dfsHG); dfsHG=0;
	this->PrintTime(stderr);
	// mcs->PutHyperPartition(); 
	mcs->SaveBest=TRUE; mcs->StoreBest();
	// mcs->Sample(2,2,2,2); this->CalcLLR(mcs); mcs->StoreBest();
	// mcs->PutHyperPartition(); 
	this->CalcLLR(mcs); this->RestoreFinalBest(); // mcs->PutHyperPartition();
	mcs->PutHyperPartition(stderr); 
	sprintf(str,"%s_I",infile); PutCheckPoint(str,FALSE);
	result=1;
      } break;
#if 1
     case 'm':	// Merge all nodes below level 1.
      {
	mcs->PutHyperPartition(stderr);
	FILE *fp=open_file(infile,"_major.mma","w"); 
	mcs->PutMajorNodesMMA(fp); fclose(fp);
      } break;
#endif
     case 'M':	// Test PutMapContributions() routine.
#if 1
	mcs->DoEvolve();  
	mcs->SampleColumns( ); 
	mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
#endif
	mcs->PutHyperPartition(stderr);
	mcs->PutLpr(stdout);
	mcs->PutSetRelations(stdout);	// check FG & BG sets for IsConflict???
	// this->PutPatternOverlap(stderr);
	this->Put();
	this->PrintTime(stderr);
	// mcs->PutMapContributions(stdout,lpr);
	// this->MoveNodesUp(100.0);
	// this->MultiMoveNodesUp(100.0);
	result=1;
      break;
     case 't': 	// take a check point.
	// mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	mcs->PutHyperPartition(stderr); 
	sprintf(str,"%s_t",infile); PutCheckPoint(str);
	// fprintf(stderr,"file name = %s\n",str); 
	result=1;
       break;
     default:
    switch (mode){
     case 17: 	// Test check point.
	mcs->PutHyperPartition(stderr); this->CalcLLR(mcs);
	this->PrintTime(stderr);
	while(this->InsertInternalNodes(stderr,300) > 0){
		RaiseBranches(stderr,300, 1); No+=RecordChange("RaiseBranches");
		// RaiseMultiBranches(stderr,temp,iter); No+=RecordChange("RaiseMultiBranches");
		LowerBranches(stderr,300, 1); No+=RecordChange("LowerBranches");
	}
	this->PrintTime(stderr);
	mcs->SaveBest=TRUE; mcs->StoreBest();
	this->CalcLLR(mcs); this->RestoreFinalBest(); // mcs->PutHyperPartition();
	mcs->PutMapContributions(stderr);
	mcs->PutHyperPartition(stderr); 
	result=0;
	break;
     case 1:	// Test mad_typ object.
	// mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	// mcs->Sample(2,2,2,1);
	mcs->PutHyperPartition(stderr); this->CalcLLR(mcs);
	this->PrintTime(stderr);
#if 1
	result = this->InsertInternalNodes(stderr,300);
	this->PrintTime(stderr);
	mcs->SaveBest=TRUE; mcs->StoreBest();
	this->CalcLLR(mcs); this->RestoreFinalBest(); //  mcs->PutHyperPartition();
	mcs->PutMapContributions(stderr); mcs->PutHyperPartition(stderr); 
	result = 1;
#else
	while(this->InsertInternalNodes(stderr,300) > 0){
		this->CalcLLR(mcs); this->RestoreFinalBest();
		mcs->PutMapContributions(stderr); mcs->PutHyperPartition(stderr); 
		this->PrintTime(stderr);
	} result=0;
#endif
      break;
     case 3: // Test DeleteWorstNode() routine.
	mcs->DoEvolve();  mcs->Sample(2,3,2,1); TrimNodes(stderr,300,iter);
	mcs->PutHyperPartition(stderr);
	DeleteWorstNode(300); 
	mcs->PutHyperPartition(stderr);
	result=1;
      break;
     case 5: // Test RtnCopy() routine.
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	xmcs = this->RtnCopy( ); CalcLLR(xmcs);
	DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
	mcs->CalcTotalLPR(); this->AdjustSetIDs(); 
	result=1;
      break;
     case 6: // Test DeleteNode() routine.
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	id=Hpt->ItoSetID(2);
	xmcs = this->DeleteNode(id); DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
	mcs->CalcTotalLPR(); id=Hpt->ItoSetID(6);
	xmcs = this->DeleteNode(id); DeleteMCS(mcs); mcs=xmcs; Hpt=mcs->GetHpt();
	mcs->CalcTotalLPR();
	this->AdjustSetIDs(); 
	result=1;
      break;
     case 9: // test OptimizeDisplaySet() routine.
	mcs->Sample(2,2,2,1); mcs->StoreBest(); mcs->RestoreBest();
	mcs->PutHyperPartition(stderr); 
	rmcs=this->Optimize(); rmcs->PutHyperPartition(stderr);
	DeleteMCS(mcs); mcs=rmcs; Hpt=mcs->GetHpt(); mcs->SaveBest=TRUE; mcs->StoreBest(); 
	result=0;
      break;
     case 12: // Test RePartition() routine.
	mcs->Sample(2,2,2,1); mcs->StoreBest(); mcs->RestoreBest(); CalcLLR(mcs);
	RePartitionSeqs(stderr,iter); 
	RandomizeSeqAssign(stderr,10,iter); mcs->SampleColumns( );
	mcs->Sample(2,2,2,1);
	result=1;
      break;
     case 13: // Test RandomizeSeqAssign() routine.
	mcs->Sample(2,2,2,1); mcs->StoreBest(); mcs->RestoreBest(); CalcLLR(mcs);
	RandomizeSeqAssign(stderr,3,iter);
	mcs->Sample(2,2,2,1);
	result=1;
      break;
     case 14: //=================  end of debugging routines. ===================
#if 0
	mcs->SampleColumns( ); mcs->StoreBest(); mcs->RestoreBest(); mcs->CalcTotalLPR();
	// mcs->Sample(2,2,2,1);
	mcs->PutHyperPartition(stderr); 
	InsertIntermediateNodes(stderr,300,iter);
	// InsertInternalNodes(stderr,300,iter);
	mcs->SaveBest=TRUE; mcs->StoreBest();
	result=1;
	result=0;
#endif
       break;
     default:
	mcs->Sample(2,2,2,1); 
	mcs->StoreBest(); // mcs->RestoreBest(); 
	TrimNodes(stderr,100,1); CalcLLR(mcs); // Store best_mcs.
	this->RestoreFinalBest(); // mcs->PutHyperPartition(); 
	mcs->PutHyperPartition(stderr); 
	result=0;
     break;
    }
     break;
   } return result;
}


