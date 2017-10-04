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

void	mcs_typ::PutSubHierarchyCntrb(FILE *fp)
// Put the contributions of each subtree as a distince hierarchy.
{
     Int4   i,j,n,*P;
     double d;
     char   **HP=Hpt->RtnHyperPartition();
     assert(Hpt->IsTree(P)); free(P); 
     assert(IsTreeHpt);	// make sure this is a tree...
     assert(Hpt->NumSets() == Hpt->NumBPPS()+1);
     CalcTotalLPR(0,FALSE); 
     for(i=1; i <= Hpt->NumBPPS(); i++){
       for(d=0.0,j=1; j < Hpt->NumSets(); j++){
	   if(HP[j][i] == '+') d+=Map[j]; // row = j; col = i; 
       } fprintf(fp,"%d\t%s\t%.2f\n",i,Hpt->ElmntSetName(i),d);
     }
}

BooLean	mcs_typ::PutMapContributions(FILE *fp,Int4 g)
// put the contributions to the Map for each set and intermediate node.
{
     assert(IsTreeHpt);	// make sure this is a tree...
     Int4	n,sq,num_sq = NumSeqsCMSA(TrueMainCMA);
     char	x,**HP=Hpt->RtnHyperPartition();
     BooLean	rtn=TRUE;

     double	lpr,lpr0,*subLLR;
     if(IsFailedSet[g]) return FALSE;
     NEW(subLLR,Hpt->NumBPPS() +3,double); CalcTotalLPR(0,FALSE); 
     for(n=0; n <= Hpt->NumBPPS(); n++){ subLLR[n]=Map[n]; }
     set_typ  FullSet=this->RtnSubTreeSeqSet(g);
     // if(g==1) fprintf(fp,"%d.%s(%d|%d): (*) %.2f\n",1, Hpt->ElmntSetName(1),CardSet(GrpSet[1]),CardSet(FullSet),subLLR[1]);
     if(0) ; else {
	assert(g > 0 && g < Hpt->NumSets());
	// Int4 card=CardSet(GrpSet[g]);
	if(FullSet) fprintf(fp,"%d.%s(%d|%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]),CardSet(FullSet));
	else fprintf(fp,"%d.%s(%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]));
	Int4 end=g-1; if(Hpt->TypeOfSet(g) == '?') end = g;
        for(n=1; n <= end; n++){ 
	   if(IsFailedBPPS[n]) continue;  // this should never be true;
	   if(HP[g][n] == '+'){
		assert(Hpt->TypeOfSet(n) == '?'); // this should correspond to an internal node.
		x='o';  // Omit from analysis.
		x='-';  // Put in background partition.

#if 0	// subtract without Null model...
		// lpr=che[n]->CalcRawLLR( ); 
		lpr=che[n]->RawLLR( ); 
		for(sq=1; sq <= num_sq; sq++){	// remove these sequences from the foreground
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition(x,sq); } 
		} // lpr0=che[n]->CalcRawLLR( ); // lpr=subLLR[n];
		lpr0=che[n]->RawLLR( ); // lpr=subLLR[n];
		fprintf(fp," (%d) %.1f",n,lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); } 
		} 

		if(FullSet == 0) continue;
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet)){ che[n]->SetPartition(x,sq); } 
		} // lpr0=che[n]->CalcRawLLR( ); // lpr=subLLR[n];
		lpr0=che[n]->RawLLR( ); // lpr=subLLR[n];
		fprintf(fp," [%.1f]",lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet)){ che[n]->SetPartition('+',sq); } 
		}
#else
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," (%d) %.1f",n,lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); } 
		} 

		if(FullSet == 0) continue;
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet)){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," [%.1f]",lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet)){ che[n]->SetPartition('+',sq); } 
		} 
#endif
	   }
	} if(end < g) fprintf(fp," (%d) %.1f",g,subLLR[g]);
     } free(subLLR); if(FullSet) NilSet(FullSet);
     return rtn;
}

#if 1
double	mcs_typ::NodeDiversity(FILE *fp, Int4 i, set_typ Set)
// return the variablity of 
{
        Int4	j,N,num_names,namelen,NewNumber;
        double	RelStdDev,inc;

    	assert(i > 0 && i <= Hpt->NumSets());
        Int4 score;
	e_type cE=GetSeqAsCsqCMSA(Set,TrueMainCMA);
        for(score=0,j=1; j <= LenSeq(cE); j++){
              unsigned char r=ResSeq(j,cE);
              if(r) score += valAlphaR(r,r,AB);
        }
        inc=ceil((double)score/25.0);
        h_type  HG=Histogram("sequence scores",-200,score,inc);
        for(Int4 sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
	      if(!MemberSet(sq,Set)) continue;
              Int4 score=PseudoAlnScoreSqToCMSA(cE,sq,TrueMainCMA);
              IncdHist((double)score,HG);
        }
        RelStdDev=100.0*sqrt(VarianceHist(HG))/MeanHist(HG);
        if(fp) fprintf(fp,"%d.%s %.2f\n",i,Hpt->ElmntSetName(i),RelStdDev);
        NilSeq(cE); 
	if(fp) PutHist(fp,50,HG); NilHist(HG);
	return RelStdDev;
}
#endif

BooLean	mcs_typ::PutSeqCntrb(FILE *fp,FILE *hfp,FILE *sfp)
// fp = contribution file pointer; hfp = histogram; sfp = sma file pointer.
{
    Int4       *P;
    if(!Hpt->IsTree(P)){ free(P); return FALSE; } else free(P);
    Int4	i,j,b,n,g,sq,num_sq = NumSeqsCMSA(TrueMainCMA);
    char        **HP=Hpt->RtnHyperPartition(),str[300];
    double	**sqLLR,d,D,*RSD;
    Int4	*NumPhyla,*NumNull;
    char	x,**Percentile,c,*NumKing;
    char   *NameSaved=AllocString(NameCMSA(TrueMainCMA));
    if(hfp==0) hfp=stdout;
    FILE *vfp=0; // vfp=stderr;
     
    h_type	vHG=Histogram("Relative Standard Deviations (RSD) of seq-to-consensus scores",-200,100,2.0);
    NEWP(sqLLR,Hpt->NumBPPS()+3,double); 
    NEWP(Percentile,Hpt->NumBPPS()+3,char); 
    NEW(NumKing,Hpt->NumBPPS()+3,char); 
    NEW(NumPhyla,Hpt->NumBPPS()+3,Int4); 
    NEW(NumNull,Hpt->NumBPPS()+3,Int4); 
    NEW(RSD,Hpt->NumBPPS()+3,double);
    set_typ *ChkSet=0; NEW(ChkSet,Hpt->NumBPPS()+3,set_typ); 
    //******** 1. Get the FG sequence LLR contributions to each contrast alignment. **********
    dh_type dH=dH=dheap(num_sq+4,4);
    x='o';  // Omit from analysis.
    x='-';  // Put in background partition.
    for(n=1; n <= Hpt->NumBPPS(); n++){
      NEW(Percentile[n],num_sq+3,char); NEW(sqLLR[n],num_sq+3,double);
      ChkSet[n]=MakeSet(num_sq+3); ClearSet(ChkSet[n]);
      h_type HG=0; sprintf(str,"sequence contributions to subLLR %d",n);
      HG=Histogram(str,-30,100,1.0);
      for(g=1; g < Hpt->NumSets(); g++){	// don't look at last (random) set!!
	if(IsFailedSet[g]) continue;
	if(CardSet(GrpSet[g]) == 0){ continue; } // Eventually print consensus seqs for these?
	if(HP[g][n] != '+') continue;		// set g is a child of or is set n.
        for(sq=1; sq <= num_sq; sq++){
    	   if(MemberSet(sq,GrpSet[g])){
		assert(!MemberSet(sq,ChkSet[n])); AddSet(sq,ChkSet[n]);
		double lpr0 = che[n]->CalcLLR( );
		che[n]->SetPartition(x,sq);
		double lpr = che[n]->CalcLLR( );
		che[n]->SetPartition(HP[g][n],sq); 
		sqLLR[n][sq]=lpr0-lpr; 
#if 1	// adjust for sequence redundancy;  make sure this is a tree!
		// che[n]->IsPatternPos(Int4 j);
		UInt4	wt=che[n]->RtnAveSqIWt(sq),factor=che[n]->GetWtFactor();
		double Wt=(double)wt/(double)factor;
		assert(Wt <= 1.0 && Wt > 0.0);
		// sqLLR[n][sq]=sqLLR[n][sq] * (1.0/Wt);
		// sqLLR[n][sq]=sqLLR[n][sq] * sqrt(1.0/Wt); // muted because similarity might not involve
		sqLLR[n][sq]=sqLLR[n][sq] * pow((1.0/Wt),(1.0/1.5)); // muted because similarity might not involve
#endif
		insrtHeap(sq,(keytyp)sqLLR[n][sq],dH);
		IncdHist(sqLLR[n][sq], HG);
	   } 
    	}
      } 
      if(NumDataHist(HG) > 0){
        double max=MaxValHist(HG),min=MinValHist(HG);
        Int4   Max=(Int4)ceil(max),Min=(Int4)floor(min);
        double inc=ceil(((double)Max-(double)Min)/30.0); 
        NilHist(HG); 
        if(inc > 0.0){ HG=Histogram(str,Min,Max,inc); } else HG=0;
      } else { NilHist(HG);  HG=0; }
      Int4 NumInFG=CardSet(ChkSet[n]);
      for(sq=1; sq <= num_sq; sq++){
	  if(MemberSet(sq,ChkSet[n])){ if(HG) IncdHist(sqLLR[n][sq],HG); }
	  else Percentile[n][sq]=-1;
      } if(HG){ PutHist(hfp,60,HG); NilHist(HG);  HG=0; }
      for(i=1; !emptyHeap(dH); i++){
	sq=delminHeap(dH); 
	d=100*(double)i/(double)NumInFG;
	Percentile[n][sq]=(char)floor(d);
      }
      if(NumInFG > 0){
	d=NodeDiversity(vfp, n, ChkSet[n]); IncdHist(d,vHG); RSD[n]=d;
        if(vfp) fprintf(vfp,"RSD(%d) = %.2f\n",n,d);
        NumPhyla[n]=PutPhylaSeqSet(vfp,0, TrueDataCMSA(TrueMainCMA),NumKing[n],NumNull[n],ChkSet[n]);
      }
    } 
    if(hfp) PutHist(hfp,50,vHG); NilHist(vHG);
    //******** 2. Output each sequence's contribution to each contrast alignment. **********
    //******** 3. And output the best sequences from each phyla as an sma file. **********
    Int4   I,numPhyla=0,max_num_phyla=50;
    Int4   *phyla=GetPhylaSeqSet(0,&numPhyla,TrueDataCMSA(TrueMainCMA));
    BooLean *Printed=0; NEW(Printed,numPhyla+3,BooLean); 
    set_typ  sSet=MakeSet(num_sq+1);

    assert(emptyHeap(dH));
    assert(Hpt->NumSets() == Hpt->NumBPPS()+1);
    for(g=1; g < Hpt->NumSets(); g++){	// don't look at last (random) set!!
      if(IsFailedSet[g]) continue;
      if(GrpSet[g] == 0 || CardSet(GrpSet[g]) == 0){
#if 1	// Put internal node information regardless...test this...
        if(fp){
          fprintf(fp,"#"); PutMapContributions(fp,g);
          fprintf(fp," RSD=%.1f <%d(%d)> {%d}\n",RSD[g],NumPhyla[g],NumKing[g],NumNull[g]); 
	}
#endif
        if(sfp){
	  sst_typ *xsst=RtnCopyOfSST(g);
	  cma_typ cma=this->RtnCsqSstAsCMSA(g,Hpt->ElmntSetName(g),xsst); free(xsst);
	  if(NumSeqsCMSA(cma) > 0) PutCMSA(sfp,cma); TotalNilCMSA(cma);
	} continue;
      }
      if(fp){
        fprintf(fp,"#"); PutMapContributions(fp,g);
        fprintf(fp," RSD=%.1f <%d(%d)> {%d}\n",RSD[g],NumPhyla[g],NumKing[g],NumNull[g]); 
      }
      for(sq=1; sq <= num_sq; sq++){
	if(MemberSet(sq,GrpSet[g])){
 	   assert(MemberSet(sq,ChkSet[g]));
	   insrtHeap(sq,-(keytyp)sqLLR[g][sq],dH);
	   // fprintf(stderr,"sq=%d; g=%d; LLR=%g\n",sq,g,sqLLR[g][sq]);
	} 
      } ClearSet(sSet); 
      for(j=1; j <= numPhyla; j++) Printed[j]=FALSE;
      // while(!emptyHeap(dH))
      // while(sq=delminHeap(dH))	// seems problematic...sq gets set to 1!? boolean? due to gcc -O3 ?
      for(I=0,j=1; (sq=delminHeap(dH)) != 0; )
      {
        // sq=delminHeap(dH);
#if 0	// DEBUG...
	fprintf(stderr,"dH->n=%d; dH->h[1]=%d",dH->n,dH->h[1]); // PutHeap(stderr,dH);
        // sq=delminHeap(dH);
	// sq=rmHeap(dH->h[1],dH);
	fprintf(stderr," sq=%d; (%d %d)\n",sq,dH->n,dH->h[1]); // PutHeap(stderr,dH);
      } exit(1);
      while(0){
#endif
#if 1	// Get set of best sequences...
	if(I < max_num_phyla){
	   if(!Printed[phyla[sq]]) { AddSet(sq,sSet); Printed[phyla[sq]]=TRUE; I++; }
	}
#endif
	if(fp) fprintf(fp,"%5d: ",j); j++;
	e_type E=TrueSeqCMSA(sq,TrueMainCMA);
	Int4 s=FakeToRealCMA(sq,1,TrueMainCMA);
	Int4 e=FakeToRealCMA(sq,LengthCMSA(1,TrueMainCMA),TrueMainCMA);
	Int4 os=OffSetSeq(E);
    	for(i=0,n=1; n <= Hpt->NumBPPS(); n++){
	  if(!MemberSet(sq,ChkSet[n])) continue;
	  assert(Percentile[n][sq] >= 0);
	  if(fp){
	    if(g != n) fprintf(fp,"(%d) %.1f [%d] ",n,sqLLR[n][sq],Percentile[n][sq]);
	    else fprintf(fp,"(%d) %.1f ",n,sqLLR[n][sq]);
	  } 
	  // else fprintf(fp,"(%d) %.1f ",n,sqLLR[n][sq]);
	} c=KingdomSeq(E); if(c==0) c = 'U';
	if(fp){
	  if(PhylumSeq(E)) fprintf(fp," |%d:%d| <%s(%c)> ",s,e,PhylumSeq(E),c);
	  else fprintf(fp," |%d:%d| <%s(%c)> ",s,e,"unknown",'U');
	  PutSeqID(fp,E); fprintf(fp,"\n");
	}
	if(sqLLR[g][sq] > 0.0 && PdbSeq(TrueSeqCMSA(sq,TrueMainCMA))){ AddSet(sq,sSet); }
      } if(fp) fprintf(fp,"\n");
      ReNameCMSA(Hpt->ElmntSetName(g),TrueMainCMA);
      // fprintf(stderr,"CardSet(bst)=%d; CardSet(pdb)=%d\n",CardSet(sSet),CardSet(pdbSet));
      if(sfp){
	set_typ set[4]; set[1]=0; set[2]=sSet;
	sst_typ *xsst=RtnCopyOfSST(g);
	cma_typ cma[4]; cma[1]=this->RtnCsqSstAsCMSA(g,Hpt->ElmntSetName(g),xsst); free(xsst);
	cma[2]=TrueMainCMA; 
	// char *pttrn=RtnCopyOfPttrn(g);
	// fprintf(stderr,"%d.%s: %s\n",g,Hpt->ElmntSetName(g),pttrn);
	// if(g==9) PutCMSA(stderr,cma[1]);
	if(CardSet(sSet) < 1) PutCMSA(sfp,cma[1]); else PutMergedCMSA(sfp,2,set,cma,0); 
	TotalNilCMSA(cma[1]);
      }
    } Nildheap(dH); dH=0;

    //******** 4. Free up memory. **********
    free(phyla); NilSet(sSet); free(Printed); 
    for(n=1; n <= Hpt->NumBPPS(); n++){ free(sqLLR[n]); free(Percentile[n]); NilSet(ChkSet[n]); }
    free(sqLLR); free(Percentile);  free(ChkSet);
    ReNameCMSA(NameSaved,TrueMainCMA); free(NameSaved); 
    free(NumPhyla); free(NumKing); free(NumNull); free(RSD);
}


static sst_typ	*Str2SST(char *sst_str,cma_typ qcma)
// Use this to convert a pattern string into small sets
// for MainSet pattern (only upon output to file).
{
	sst_typ	*xsst=0,*qsst;
	a_type ab=AlphabetCMSA(qcma);
	Int4	i,s,r,nval,*pos;

	e_type Query=SeqSetE(1,DataCMSA(qcma));  // Fake sequence...

	NEW(qsst,LenSeq(Query)+2,sst_typ); NEW(xsst,strlen(sst_str)+2,sst_typ);
        NEW(pos,strlen(sst_str)+3,Int4);
        nval=ParseResidueSets(sst_str,pos,xsst,ab,"che_typ input error");
        for(i=1;i <= nval; i++){
		// s = RealToFakeCMSA(1,pos[i],qcma); // s == actual position in array.
		s = pos[i];
		// s == fake seq position for residue s in sq (with offset).
		if(s <= 0){
		   fprintf(stderr,"sst_str=%s\n",sst_str);
		   fprintf(stderr,"Position %d(%d):\n",pos[i],s);
		   print_error("che_typ input pattern option error");
		}
                if(s > LenSeq(Query)){ print_error("che_typ input option error"); }
                r=ResSeq(s,Query);
                // fprintf(stdout,"%c%d(%d)\n",AlphaChar(r,ab),pos[i],s);
		if(!MemSset(r,xsst[i])){
			fprintf(stderr,
			  "**************** Str2SST( ) FATAL ERROR! **************\n");
                	fprintf(stderr,"%c%d(%d)\n",AlphaChar(r,ab),pos[i],s);
			fprintf(stderr,"xsst[%d]: \"",pos[i]);
			PutSST(stderr,xsst[i],ab);
			fprintf(stderr,"\"\n");
			PutSeq(stderr,Query,ab);
			gsq_typ *gsq=gsqCMSA(1,qcma);
			gsq->Put(stderr,ab); gsq->Put(stderr,60,ab);
			fprintf(stderr,"%d(%d): %c\n",pos[i],s,AlphaChar(r,ab));
			fprintf(stderr,"Query offset = %d\n",OffSetSeq(Query));
			fprintf(stderr,"%s\n",sst_str);
			print_error("che_typ sst input error");
		} else qsst[s]=xsst[i];
        }
	free(xsst); free(pos); // exit(1);
	return qsst;
}

static void CheckCardinality(Int4 n, set_typ set, Int4 num, cma_typ cma)
{
	if(CardSet(set) != NumSeqsCMSA(cma)){
	  fprintf(stderr,"CardSet(fg_set[%d]=%d; NumSeqsCMSA(rtn_cma[%d]=%d\n",
                                n,CardSet(set),num,NumSeqsCMSA(cma));
	  assert(CardSet(set) == NumSeqsCMSA(cma));
	}
}

#if 0
void	mcs_typ::PutContinueFile(FILE *fp, BooLean Label)
// output a file for the next round of amcBPPS sampling...
{
  Int4	sq,g;

  RestoreBest();
  // 2. Label TrueMainCMA file...
  if(Label){
    for(g=1; g <= Hpt->NumSets(); g++){	
      if(IsFailedSet[g]){ assert(CardSet(GrpSet[g]) == 0); continue; }
      if(CardSet(GrpSet[g]) > 0){	// output optimized minimal set if non-empty.
	//*********************** Random Set ****************************
	if(Hpt->TypeOfSet(g) != '?' && strcmp(Hpt->ElmntSetName(g),"Random") != 0){
	    for(sq = 1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
		if(MemberSet(sq,GrpSet[g])){	// if member of the set...
		    e_type fE = FakeSeqCMSA(sq,TrueMainCMA);
		    e_type tE = TrueSeqCMSA(sq,TrueMainCMA);
		    LabelSeq(tE); LabelSeq(fE);
		}
	    }
	} 
      }
    }
  } PutCMSA(fp,TrueMainCMA);
}
#endif

void	mcs_typ::PutLpr(char *outfile)
{ FILE *fp=open_file(outfile,".lpr","w"); this->PutLpr(fp); fclose(fp); }

void	mcs_typ::PutLpr(FILE *fp)
{
     Int4	*P;
     BooLean	is_tree=Hpt->IsTree(P); free(P);
     for(Int4 n=Hpt->NumBPPS(); n > 0; n--){
	if(IsFailedBPPS[n]) continue;
	if(!Hpt->OutputThis(n)) continue;
	if(che[n]){
	  // fprintf(stderr,"Map[%d] = %g\n",n,Map[n]);
	  CalcTotalLPR(0,FALSE); // Calculates all Map[n]
	  // only don't put out if IsFailedBPPS is true.
	  if(Map[n] > 0){
		if(is_tree){
     		   fprintf(fp,":::::::::: BPPS category %d (\"%s%c\"): ::::::::::\n",
				n,Hpt->ElmntSetName(n),Hpt->TypeOfSet(n));
		} else {
     		   fprintf(fp,":::::::::: BPPS category %d: ::::::::::\n",n);
		} che[n]->PutPttrnFile(fp,n);
	  }
	}
     }
}

#if 0	// moved to mcs_cdd.cc
BooLean IsSeqEST(e_type E)
{
        char *result=strstr(E->info,"_EST");
        if(result==0) return FALSE;
        else {
                char *space=strchr(E->info,' ');
                if(space > result) return TRUE; else return FALSE;
        }
}

Int4	FirstInSet(set_typ St,cma_typ cma)
// return the first good element in set St.
{
	double cut=0.0;
	for(Int4 i=1; i <= SetN(St); i++){
		if(MemberSet(i,St)){
		   e_type sE=TrueSeqCMSA(i,cma);
		   if(IsSeqEST(sE)) continue;
		   double prob=GetGappedProbCMSA(1,i,cma);
		   if(prob >= cut) return i; 
		}
	} return 0; 
}
#endif

Int4	mcs_typ::PutMajorNodesMMA(FILE *mmafp)
{
     Int4	 n,g,d,N,*Parent;
     FILE	*ofp=0;
     cma_typ	mcma=0;

    if(mmafp) ofp=mmafp; else ofp=open_file(infile,"_major.mma","w");
    if(!Hpt->IsTree(Parent)) print_error("FATAL: hpt file does not correspond to a tree");
    for(N=0,g=1; g <= Hpt->NumSets(); g++) {
      if(IsFailedSet[g]){ continue; }
      if(g != Hpt->NumSets() && g != 1 && Parent[g] != 1) continue;
      set_typ setX,setT=Hpt->MkSubTreeSet(g);
      setX=CopySet(GrpSet[g]);
      for(n=2; n <= Hpt->NumSets(); n++){
	 if(g == 1 && CardSet(setX) == 0) break;
	 if(g != 1 && n != g && MemberSet(n,setT)){ UnionSet(setX,GrpSet[n]); }
      } ReNameCMSA(Hpt->ElmntSetName(g),MainCMA);
      if(strcmp(Hpt->ElmntSetName(g),"Random") == 0){
	    set_typ NonRandom=MakeSet(SetN(Labeled));
	    IntersectNotSet(GrpSet[g],Labeled,NonRandom);  // NonRandom = Set[g] && !Labeled.
	    if(CardSet(NonRandom) > 0){ PutInSetCMSA(ofp,NonRandom,MainCMA); N++; }
	    NilSet(NonRandom);
      } else { PutInSetCMSA(ofp,setX,MainCMA); N++; }
    } if(mmafp==0) fclose(ofp); else mmafp=ofp; free(Parent);
    return N;
}

double  mcs_typ::Put(BooLean put_pttrn, char **names,BooLean do_put_cdd, FILE *mmafp,FILE *ptrnfp)
{
     Int4	 n,g;
     FILE	*ofp=0;
     double	MinKeyFrq=0.5;
     char	tmp_str[200];
     cma_typ mcma=0;
     double	*fract_ignored,TotalMap=0.0;
     NEW(fract_ignored,25,double);

     Int4 TotalCHA=0;
     for(n=Hpt->NumBPPS(); n > 0; n--){ if(che[n] && Map[n] > 0){ TotalCHA++; } }
     Int4    NumberFound=0;
     FILE *pttrn_fp=0; if(put_pttrn) pttrn_fp=open_file(infile,".pttrns","w");
     FILE *ptninfo_fp=0; if(verbose) ptninfo_fp=open_file(infile,".lpr","w");
     FILE *info_fp=0; if(verbose) info_fp=open_file(infile,".info","w");
     for(n=Hpt->NumBPPS(); n > 0; n--){
	if(IsFailedBPPS[n]) continue;
	if(!Hpt->OutputThis(n)) continue;
	if(che[n]){
	  // fprintf(stderr,"Map[%d] = %g\n",n,Map[n]);
	  CalcTotalLPR(0,FALSE); // Calculates all Map[n]
     	  if(pttrn_fp){
	    fprintf(pttrn_fp,"%d: ",n); che[n]->PutBestPatterns(pttrn_fp,FALSE);
	  } else if(ptrnfp){
	    fprintf(ptrnfp,"%d: ",n); che[n]->PutBestPatterns(ptrnfp,FALSE);
	  }
	  // only don't put out if IsFailedBPPS is true.
	  if(Map[n] > 0){
#if 1
		bpps_typ *bpps=che[n]->BPPS();
	        mcma = che[n]->RtnMainSet();
	        double alpha,A0,B0;
	        alpha = bpps->Parameters(A0,B0);
	        SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,set_mode,mcma);
#endif
		TotalMap += Map[n];
		NumberFound++;
		if(PutIntermediateFiles){
			che[n]->ContribSeqMAP(n,SFBG[n],fract_ignored,MinNats,MinKeyFrq);
			// che[n]->PutAll(n,SFBG[n],MinNats,MinKeyFrq);
			che[n]->PutFgBgSeqs(n);
		}
		if(verbose){
		  che[n]->PutInfo(info_fp,n);
     		  // fprintf(pttrn_fp,"%d: ",n); che[n]->PutBestPatterns(pttrn_fp,FALSE);
     		  fprintf(ptninfo_fp,":::::::::: BPPS category %d: ::::::::::\n",n);
		  che[n]->PutPttrnFile(ptninfo_fp,n);
		}
	  }
	}
     } free(fract_ignored);	// set all to 0.0
     if(verbose){ fclose(info_fp); fclose(ptninfo_fp);} 
     if(pttrn_fp) fclose(pttrn_fp); 
    Int4 num_sq = NumSeqsCMSA(MainCMA);
    // sprintf(tmp_str,"_new.mma");
    if(mmafp) ofp=mmafp; else ofp=open_file(infile,"_new.mma","w");
    //  char	tmp_name[200];
    // sprintf(tmp_name,"%s",NameCMSA(MainCMA)); 
    for(g=1; g <= Hpt->NumSets(); g++){	
      // if(IsFailedSet[g]){ assert(CardSet(GrpSet[g]) == 0); continue; }
      if(IsFailedSet[g]){ continue; }
      if(PutIntermediateFiles){
       FILE *sfp=0;
       sprintf(tmp_str,"_set%d.seq",g);
       for(Int4 sq=1; sq <= num_sq; sq++){
        if(MemberSet(sq,GrpSet[g])){
          if(sfp==0) sfp=open_file(infile,tmp_str,"w");
	  PutSeq(sfp,TrueSeqCMSA(sq,MainCMA),AB);;
         }
       } fclose(sfp);
      }
      if(CardSet(GrpSet[g]) > 0){	// output optimized minimal set if non-empty.
#if 1
	 if(names && names[g]) ReNameCMSA(names[g],MainCMA);
	 else ReNameCMSA(Hpt->ElmntSetName(g),MainCMA);
#else
         ReNameCMSA(Hpt->ElmntSetName(g),MainCMA);
#endif
	 if(strcmp(Hpt->ElmntSetName(g),"Random") == 0){
	    set_typ NonRandom=MakeSet(SetN(Labeled));
	    IntersectNotSet(GrpSet[g],Labeled,NonRandom);  // NonRandom = Set[g] && !Labeled.
            // UnLabelPutInSetCMSA(ofp,NonRandom,MainCMA);
	    if(CardSet(NonRandom) > 0){
                PutInSetCMSA(ofp,NonRandom,MainCMA);	// keep labels as unputted
	    } NilSet(NonRandom);
#if 0	// added in order to lock sequences into a misc set on next iteration. 
	 } else if(Hpt->TypeOfSet(g) != '?'){
		LabelPutInSetCMSA(ofp,GrpSet[g],MainCMA);	// Label seqs in non-misc sets.
#endif
	 } else PutInSetCMSA(ofp,GrpSet[g],MainCMA);	// keep labels as inputted
      }
    } if(mmafp==0) fclose(ofp); else mmafp=ofp;
    // output sets for CDtree with first sequence for each descendant node of an internal node.
    if(do_put_cdd){	// see mcs_cdd.cc file...
      ofp=open_file(infile,"_cdd.mma","w"); this->PutCDD(ofp); fclose(ofp);
      ofp=open_file(infile,"_cdd.hpt","w"); Hpt->Put(ofp);  fclose(ofp);
      ofp=open_file(infile,"_cdd.ptn","w"); this->PutPttrnLLRs(ofp); fclose(ofp);
    }

#if 0	// skip this...
    sprintf(tmp_str,"_new.dma");
    ofp=open_file(infile,tmp_str,"w");
    for(i=1; i <= NumDisplayCMA; i++){
	// char *name=Hpt->ElmntSetName(i);
	// Hpt->ReNameFamilies(i,name);	// rename input cma for *.hpt file; fails if i disallowed.
	PutCMSA(ofp,DisplayCMA[i]); 	// then put display files
    }
    fclose(ofp);
#endif

#if 0
     if(PutIntermediateFiles){
       sprintf(tmp_str,"_new.hpt");
       ofp=open_file(infile,tmp_str,"w"); Hpt->Put(ofp); fclose(ofp);
     }
#endif
     return TotalMap;
}

void	mcs_typ::PutPttrns(FILE *fp)
{
     for(Int4 n=Hpt->NumBPPS(); n > 0; n--){
	if(IsFailedBPPS[n]) continue;
	if(!Hpt->OutputThis(n)) continue;
	if(che[n]){ fprintf(fp,"%d: ",n); che[n]->PutBestPatterns(fp,FALSE); }
     } 
}

BooLean	mcs_typ::PutEnhancedDisplaySet(FILE *sfp)
// sfp = sma file pointer.
{
    Int4       *P;
    if(!Hpt->IsTree(P)){ free(P); return FALSE; } else free(P);
    Int4	i,j,b,n,g,sq,num_sq=NumSeqsCMSA(TrueMainCMA);
    char        **HP=Hpt->RtnHyperPartition(),str[300];
    double	**sqLLR,d,D;
    char	x,c,*NameSaved=AllocString(NameCMSA(TrueMainCMA));
     
    NEWP(sqLLR,Hpt->NumBPPS()+3,double); 
    set_typ *ChkSet=0; NEW(ChkSet,Hpt->NumBPPS()+3,set_typ); 
    //******** 1. Get the FG sequence LLR contributions to each contrast alignment. **********
    x='o';  // Omit from analysis.
    x='-';  // Put in background partition.
    for(n=1; n <= Hpt->NumBPPS(); n++){
      NEW(sqLLR[n],num_sq+3,double);
      ChkSet[n]=MakeSet(num_sq+3); ClearSet(ChkSet[n]);
      for(g=1; g < Hpt->NumSets(); g++){	// don't look at last (random) set!!
	if(IsFailedSet[g]) continue;
	if(CardSet(GrpSet[g]) == 0){ continue; } // Eventually print consensus seqs for these?
	if(HP[g][n] != '+') continue;		// set g is a child of or is set n.
        for(sq=1; sq <= num_sq; sq++){
    	   if(MemberSet(sq,GrpSet[g])){
		assert(!MemberSet(sq,ChkSet[n])); AddSet(sq,ChkSet[n]);
		double lpr0 = che[n]->CalcLLR( );
		che[n]->SetPartition(x,sq);
		double lpr = che[n]->CalcLLR( );
		che[n]->SetPartition(HP[g][n],sq); 
		sqLLR[n][sq]=lpr0-lpr; 
#if 0
		// adjust for sequence redundancy;  make sure this is a tree!
		UInt4	wt=che[n]->RtnAveSqIWt(sq),factor=che[n]->GetWtFactor();
		double Wt=(double)wt/(double)factor;
		assert(Wt <= 1.0 && Wt > 0.0);
		sqLLR[n][sq]=sqLLR[n][sq] * pow((1.0/Wt),(1.0/1.5));
		// muted because similarity might not involve
#endif
	   } 
    	}
      } 
    } 
    //******** 2. Output the best sequences from each phyla as an sma file. **********
    Int4   I,numPhyla=0,max_num_phyla=30;
    Int4   *phyla=GetPhylaSeqSet(0,&numPhyla,TrueDataCMSA(TrueMainCMA));
    BooLean *Printed=0; NEW(Printed,numPhyla+3,BooLean); 
    set_typ  sSet=MakeSet(num_sq+1);
    assert(Hpt->NumSets() == Hpt->NumBPPS()+1);
    for(g=1; g < Hpt->NumSets(); g++){	// don't look at last (random) set!!
      if(IsFailedSet[g]) continue;
      if(GrpSet[g] == 0 || CardSet(GrpSet[g]) == 0){
        if(sfp){
	  sst_typ *xsst=RtnCopyOfSST(g);
	  cma_typ cma=this->RtnCsqSstAsCMSA(g,Hpt->ElmntSetName(g),xsst); free(xsst);
	  if(NumSeqsCMSA(cma) > 0) PutCMSA(sfp,cma); TotalNilCMSA(cma);
	} continue;
      }
      dh_type dH=dH=dheap(num_sq+4,4);
      for(sq=1; sq <= num_sq; sq++){
	if(MemberSet(sq,GrpSet[g])){
 	   assert(MemberSet(sq,ChkSet[g]));
	   if(sqLLR[g][sq] > 0.05){
#if 0	// DEBUG
		if(sqLLR[g][sq] < 1.0 || (g==13 && sq==462)) 
		  fprintf(stderr,"%s(%d): %.3f\n",Hpt->GrpName(g),g,sqLLR[g][sq]);
#endif
		insrtHeap(sq,-(keytyp)sqLLR[g][sq],dH);
	   }
	} 
      } ClearSet(sSet); 
      for(j=1; j <= numPhyla; j++) Printed[j]=FALSE;
      for(I=0,j=1; (sq=delminHeap(dH)) != 0; ) {
	if(I < max_num_phyla){	// Get set of best sequences...
	   if(!Printed[phyla[sq]]){ 
#if 0	// DEBUG...
		StrSeqID(str,50,TrueSeqCMSA(sq,TrueMainCMA));
		fprintf(stderr,"%s(%d).%d: %.3f = %s\n",Hpt->GrpName(g),g,sq,sqLLR[g][sq],str);
#endif
		AddSet(sq,sSet); Printed[phyla[sq]]=TRUE; I++; 
	   }
	} j++;
	e_type E=TrueSeqCMSA(sq,TrueMainCMA); c=KingdomSeq(E); if(c==0) c = 'U';
	if(sqLLR[g][sq] > 0.05 && PdbSeq(TrueSeqCMSA(sq,TrueMainCMA))){ AddSet(sq,sSet); }
      }

#if 1	// add additional sequences if too few of them...
      Int4 min_displayed=5,percent_ident=80;
      //******** 3. Add sequences to small sets. **********
      if((I=CardSet(sSet)) < min_displayed){
	while((sq=delminHeap(dH)) != 0) ;
	set_typ rpSet=0,inSet=CopySet(sSet); ClearSet(inSet);
        for(sq=1; sq <= num_sq; sq++){
	   if(MemberSet(sq,GrpSet[g]) && !MemberSet(sq,sSet)){
	      if(sqLLR[g][sq] > 0.05 ){
		 AddSet(sq,inSet); insrtHeap(sq,-(keytyp)sqLLR[g][sq],dH);
	      }
	   }
	}
	rpSet=RtnFastRepSetCMSA(0,percent_ident,inSet,TrueMainCMA);
        for(j=1; (sq=delminHeap(dH)) != 0; ) {
	   if(MemberSet(sq,rpSet)){ AddSet(sq,sSet); I++; }
	   if(I >= min_displayed) break;
	}
      }
#endif

      ReNameCMSA(Hpt->ElmntSetName(g),TrueMainCMA);
      if(sfp){
	set_typ set[4]; set[1]=0; set[2]=sSet;
	sst_typ *xsst=RtnCopyOfSST(g);
	cma_typ cma[4]; cma[1]=this->RtnCsqSstAsCMSA(g,Hpt->ElmntSetName(g),xsst); free(xsst);
	cma[2]=TrueMainCMA; 
	if(CardSet(sSet) < 1) PutCMSA(sfp,cma[1]); else PutMergedCMSA(sfp,2,set,cma,0); 
	TotalNilCMSA(cma[1]);
      } Nildheap(dH); dH=0;
    } 

    //******** 4. Free up memory. **********
    free(phyla); NilSet(sSet); free(Printed); 
    for(n=1; n <= Hpt->NumBPPS(); n++){ free(sqLLR[n]); NilSet(ChkSet[n]); }
    free(sqLLR); free(ChkSet);
    ReNameCMSA(NameSaved,TrueMainCMA); free(NameSaved); 
}

void	mcs_typ::PutRTF(BooLean updateCSQ, BooLean SaveChnFiles,Int4 KeyStart,Int4 KeyEnd,char *KeySetName, double **value)
//********************* output n-tier contrast hierarchical alignments. ********************
//********************* output n-tier contrast hierarchical alignments. ********************
//********************* output n-tier contrast hierarchical alignments. ********************
{
     Int4	*DisplayHits,NumTiers,g,q,i,n;
     double	MinKeyFrq=0.5;
     FILE	*ofp=0;
     char	tmp_str[200],tmp[200];
     FILE	*rtf_fp=0,*debug_fp=0;
     cma_typ	*DisplayedCMA=DisplayCMA;
#if 1
     if(!user_display_set){ 	// then create an alignment on the fly.
	// fprintf(stderr,"Creating a display alignment on the fly.\n");
	ofp = tmpfile(); this->PutEnhancedDisplaySet(ofp); rewind(ofp); 
	DisplayedCMA=MultiReadCMSA(ofp,&n,AB); fclose(ofp); assert(n == NumDisplayCMA);
     } else {
	fprintf(stderr,"User-based display set for contrast alignment.\n");
     }
#endif

     RestoreBest();
#if 1	// reset the display consensus to match the evolving consensus.
     if(updateCSQ){ this->UpdateDisplaySeqs( ); }
#endif
     if(cfp) fprintf(cfp,"%d %.1f %.0f %d\n",Iteration,TotalLPR,temperature,TotalColumns( ));
     NEW(DisplayHits,NumDisplayCMA+3,Int4);
     for(n=1; n <= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n]) continue;
	if(che[n]->NumColumns( ) < 1) continue;	// fixes core dump with "78: (454 0  seqs)(0)(0 cols)".
	for(g=1; g<=Hpt->nGrpsFG(n); g++){
           if(IsFailedSet[g]) continue;
	   q = Hpt->GrpsFG(n,g);
	   assert(q>0 && q<=NumDisplayCMA);
	   DisplayHits[q]++; 
	}
     }
     cma_typ mcma=0;
     for(q=1; q <= NumDisplayCMA; q++){
#if 1	// -print=R<int>:<int>=<name> option...
       if(KeySetName && strcmp(KeySetName,Hpt->ElmntSetName(q)) != 0) continue;
#endif
       if(IsFailedSet[q]) continue;
       if(Hpt->TypeOfSet(q) != '!') continue;	// Skip over sets not designated for output.
#if 1	// Fix for Array2SeqSet(sequence_type**, Int4, char*, alphabet_type*): Assertion `N > 0' failed.
       Int4 CardFG=1;
       for(n=Hpt->NumBPPS(); n >= 1; n--){
	  if(MemberSet(q,SetFG[n]) && CardSet(SetFG[n])==1){ // is this a leaf without seqs?
		CardFG = CardSet(che[n]->RtnFG_Set());  break;
	  }
       } if(CardFG < 1) continue;
	// SaveChnFiles=TRUE;
	// if(q==24) fprintf(stderr,"CardFG=%d\n",CardFG); // DEBUG
#endif
       BooLean	UsedSomewhere=FALSE;
       for(n=Hpt->NumBPPS(); n >= 1; n--){
	  if(IsFailedBPPS[n]) continue;
	  if(che[n]->NumColumns( ) < 1) continue;	// fixes core dump with "78: (454 0  seqs)(0)(0 cols)".
	  if(Hpt->OutputThis(n) && Hpt->RtnHyperPartition(q,n) == '+') UsedSomewhere=TRUE; 
       } if(!UsedSomewhere) continue;

	// 2. get sequence sets for fg and bg...
       // if(DisplayHits[q] > 1)	// this seemed to be working!?  Check into this later...
	// above assumes all subgroups modeled in the highest Misc category; not true if only 1 BPPS done.
       if(DisplayHits[q] > 0) {
	NumTiers=0;
	if(SaveChnFiles){ sprintf(tmp_str,"_set%d.chn",q); ofp=open_file(infile,tmp_str,"w"); }
	else ofp = tmpfile();
        set_typ	*fg_set,*bg_set; 
	NEW(fg_set, Hpt->NumBPPS()+3,set_typ); NEW(bg_set, Hpt->NumBPPS()+3,set_typ);
	hsw_typ HSW=0;
	if(passed_in_hsw){ HSW=passed_in_hsw; }
	else { HSW=chn[1]->RtnHSW(1); }	// Pass Henikoff weights on to other analyses
	hsw_typ *hsw;
	NEW(hsw, Hpt->NumBPPS()*2 + 3, hsw_typ);
	Int4	NumAnal=0;
#if 0
        sprintf(tmp_str,"_set%d.pttrn",q); FILE *pttrn_fp=open_file(infile,tmp_str,"w");
#else
	FILE *pttrn_fp=0;
#endif

	for(n=Hpt->NumBPPS(); n >= 1; n--){
	  if(IsFailedBPPS[n]) continue;
	  if(!Hpt->OutputThis(n)) continue;
          if(!MemberSet(q,SetFG[n])) continue;
	  if(che[n]->NumColumns( ) < 1) continue;	// fixes core dump with "78: (454 0  seqs)(0)(0 cols)".
          if(pttrn_fp){ fprintf(pttrn_fp,"%d: ",n); che[n]->PutBestPatterns(pttrn_fp,FALSE);}
	  mcma = che[n]->RtnMainSet();
          if(NumTiers==0){
	    NumTiers++;
	    PutCMSA(ofp,DisplayedCMA[q]);
	    // che[n]->PutMainSet(ofp,MinNats,MinKeyFrq);
	    sprintf(tmp,"%s",NameCMSA(mcma)); 
	    ReNameCMSA(Hpt->GrpName(n),mcma);
	    bpps_typ *bpps=che[n]->BPPS();
	    double alpha,A0,B0;
	    alpha = bpps->Parameters(A0,B0);
	    if(n==1) SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,'R',mcma);
	    else SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,set_mode,mcma);
#if 0		// Debug...
	    if(GlobalN){
		fprintf(stderr,"%d: GlobalN =%d; Contrast=%d\n",n,GlobalN,che[n]->Contrast);
		assert(GlobalN == che[n]->Contrast);
	    }
#endif
	    che[n]->PutFG(ofp,MinNats,MinKeyFrq);
	    fg_set[n]=che[n]->RtnFG_Set( );
	    NumAnal++; 	//********** rtf stuff **********
	    BooLean	*skip;
	    NEW(skip,NumSeqsCMSA(DisplayedCMA[q])+3,BooLean);
	    for(i=1; i<=NumSeqsCMSA(DisplayedCMA[q]); i++) skip[i]=TRUE;
	    for(i=1; i<=NumSeqsCMSA(DisplayedCMA[q]); i++){
		skip[i]=FALSE;
		PutSelectCMSA(ofp,skip,DisplayedCMA[q]);
		skip[i]=TRUE;
	    } free(skip); 
	    che[n]->PutBG(ofp,MinNats,MinKeyFrq); 
	    bg_set[n]=che[n]->RtnBG_Set( );
	    NumAnal++; 	//********** rtf stuff **********
	    ReNameCMSA(tmp,mcma);
	  } else {
	    NumTiers++;
	    sprintf(tmp,"%s",NameCMSA(mcma)); 
	    ReNameCMSA(Hpt->GrpName(n),mcma);
	    bpps_typ *bpps=che[n]->BPPS();
	    double alpha,A0,B0;
	    alpha = bpps->Parameters(A0,B0);
	    if(n==1) SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,'R',mcma);
	    else SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,set_mode,mcma);
#if 0		// Debug...
	    if(GlobalN){
		fprintf(stderr,"%d: GlobalN =%d; Contrast=%d\n",n,GlobalN,che[n]->Contrast);
		assert(GlobalN == che[n]->Contrast);
	    }
#endif
	    che[n]->PutFG(ofp,MinNats,MinKeyFrq);
	    fg_set[n]=che[n]->RtnFG_Set( );
	    NumAnal++; 	//********** rtf stuff **********
	    che[n]->PutBG(ofp,MinNats,MinKeyFrq);
	    bg_set[n]=che[n]->RtnBG_Set( );
	    NumAnal++; 	//********** rtf stuff **********
	    ReNameCMSA(tmp,mcma);
	  }
	}
        if(pttrn_fp){ fclose(pttrn_fp); pttrn_fp=0; }
	if(NumTiers > 1 && sst_str[0] != 0){	// put main set if pattern is specified.
	    cma_typ qcma = che[1]->RtnQuerySet();
	    sst_typ *msst = Str2SST(sst_str[0],qcma);
	    set_typ Set=MakeSet(NumSeqsCMSA(MainCMA)+1); 
	    FillSet(Set); DeleteSet(0,Set);	// WARNING: Don't count zero element.

	    sprintf(tmp,"%s",NameCMSA(MainCMA)); 
	    ReNameCMSA("MainSet",MainCMA);
	    PutInSetCMSA(ofp,Set,msst,MainCMA);
	    ReNameCMSA(tmp,MainCMA);
	    NilSet(Set);
	}
	if(SaveChnFiles){
		fclose(ofp); sprintf(tmp_str,"_set%d.chn",q);
		ofp=open_file(infile,tmp_str,"r"); 
	} else rewind(ofp); 
	{
	 char    *ArgV[100];
	 Int4     ArgC=0,NumCMA,num=0,NumHSW=0;
    	 // sprintf(tmp_str,"%s_set%d",infile,q);
    	 sprintf(tmp_str,"%s",Hpt->ElmntSetName(q));
	 ArgV[ArgC]=AllocString("chn_see"); ArgC++;
	 ArgV[ArgC]=AllocString(tmp_str); ArgC++;
	 // ArgV[ArgC]=AllocString("-method=U"); ArgC++;
	 if(0 && KeySetName) ArgV[ArgC]=AllocString("-S=L"); else ArgV[ArgC]=AllocString("-S=P"); ArgC++;
	 ArgV[ArgC]=AllocString("-concise");  ArgC++;
	 ArgV[ArgC]=AllocString("-sets=L"); ArgC++;
	 ArgV[ArgC]=AllocString("-sets_root=R"); ArgC++; // set root root to 'R' 
#if 1	// Add page formating for Kannan; 
    	 sprintf(tmp_str,"-S=%c",page_format);
	 ArgV[ArgC]=AllocString(tmp_str); ArgC++;
    	 sprintf(tmp_str,"-F=%d",font_size);
	 ArgV[ArgC]=AllocString(tmp_str); ArgC++;
#else
	 ArgV[ArgC]=AllocString("-F6"); ArgC++;
	 // ArgV[ArgC]=AllocString("-F10"); ArgC++;
#endif
	 if(KeySetName){
		sprintf(tmp_str,"-B=%d:%d",KeyStart,KeyEnd); ArgV[ArgC]=AllocString(tmp_str); ArgC++;
		ArgV[ArgC]=AllocString("-Nth=2"); ArgC++; // use first seq after concensus...
		// ArgV[ArgC]=AllocString("-C=BGyORMX"); ArgC++; // use first seq after concensus...
		ArgV[ArgC]=AllocString("-C=XBCGyORMX"); ArgC++; // use first seq after concensus...
	 }
	 if(NthSeqForDisplay > 0){ 
	    sprintf(tmp_str,"-Nth=%d",NthSeqForDisplay); ArgV[ArgC]=AllocString(tmp_str); ArgC++;
	 }
	 if(NumHighlighted > 0){
	    sprintf(tmp_str,"-N=%d",NumHighlighted); ArgV[ArgC]=AllocString(tmp_str); ArgC++;
	 }
	 //*********** passed in arguments...
	 for(Int4 a=0; a < Hpt->RtnArgC(q); a++){ ArgV[ArgC] = Hpt->SetArgV(q,a); ArgC++; }
	 char **Status=0; 
	 cma_typ *rtn_cma=MultiReadCMSA(ofp,&NumCMA,&Status,AB);
#if 1	// check for null alignments.
	 for(Int4 y=1; y <= NumCMA; y++) assert(NumSeqsCMSA(rtn_cma[y]) > 0);
#endif
double ***ResEvals; NEWPP(ResEvals,NumAnal +3, double);
	 // Some of these ^ input cma files have zero sequences; fix this...!!!
	 fclose(ofp); ofp=0;
	 {	// scope for chn_typ tmpchn.
	   for(NumHSW=num=0,n=Hpt->NumBPPS(); n >= 1; n--){
	        if(IsFailedBPPS[n]) continue;
	 	if(!Hpt->OutputThis(n)) continue;
         	if(!MemberSet(q,SetFG[n])) continue;
	  	if(che[n]->NumColumns( ) < 1) continue;	// fixes core dump with "78: (454 0  seqs)(0)(0 cols)".
	  	mcma = che[n]->RtnMainSet();
		if(num == 0){
		    num=2;
		    CheckCardinality(n,fg_set[n],num,rtn_cma[num]);	// FG cma file.
		    NumHSW++; hsw[NumHSW] = GetSubHSW(HSW,fg_set[n],rtn_cma[num],mcma); num++;
		    ResEvals[NumHSW] = this->GetResEvals(n);    // foreground only...
		    // Int4 NumSeqAln=NumSeqsCMSA(rtn_cma[1]); // rtn_cma[1] == display set.
		    for(Int4 sq=1; sq <= NumSeqsCMSA(rtn_cma[1]); sq++) num++; // skip dummy cmas
		    CheckCardinality(n,bg_set[n],num,rtn_cma[num]);	// BG cma file.
	 	    NumHSW++; hsw[NumHSW] = GetSubHSW(HSW,bg_set[n],rtn_cma[num],mcma); num++;
		    ResEvals[NumHSW] = ResEvals[NumHSW-1];
		} else {
		    CheckCardinality(n,fg_set[n],num,rtn_cma[num]);
		    NumHSW++; hsw[NumHSW] = GetSubHSW(HSW,fg_set[n],rtn_cma[num],mcma); num++;
		    ResEvals[NumHSW] = this->GetResEvals(n);
		    CheckCardinality(n,bg_set[n],num,rtn_cma[num]);
	 	    NumHSW++; hsw[NumHSW] = GetSubHSW(HSW,bg_set[n],rtn_cma[num],mcma); num++;
		    ResEvals[NumHSW] = ResEvals[NumHSW-1];
		}
	   } assert(NumHSW == NumAnal);
           chn_typ tmpchn(ArgC,ArgV,NumCMA,rtn_cma,hsw,Status,200,ResEvals);
           for(n=1; n <= NumHSW; n=n+2){
                for(Int4 s=1; s <= RtnLengthMainCMSA( ); s++) free(ResEvals[n][s]);
                free(ResEvals[n]);
           } free(ResEvals);
	   // chn_typ tmpchn(ArgC,ArgV,NumCMA,rtn_cma,0,Status,5); // computes hsw over again.
	   // percent_cut == 5 < 100 will remove sequences...
#if 0	// Add in values here...
	   double *val=0; NEW(val,RtnLengthMainCMSA( )+3,double);
	   for(Int4 z=1; z <= RtnLengthMainCMSA( ); z++) val[z]=-1.0;
           val[16]=25; val[24]=15; val[26]=10; val[35]=23;
	   tmpchn.AddXconserved(val);	// values will be freed by tmpchn.Free();
#else
	   if(value && value[q]){
	      for(Int4 z=1; z <= RtnLengthMainCMSA( ); z++){ 
		value[q][z]= 50*value[q][z]/value[q][0];	// 50 == rtf->hist_height
		if(value[q][z] < 1.0) value[q][z]=-1.0;
	      } tmpchn.AddXconserved(value[q]);	// values will be freed by tmpchn.Free();
	   }
#endif
	   assert(!tmpchn.OwnsCMAs());
	   if(PrintEachRTF){ sprintf(tmp_str,"%s_set%d",infile,q); tmpchn.PutHierarchicalAln(tmp_str); }
	   if(rtf_fp==0){
     	        sprintf(tmp_str,"_aln.rtf"); rtf_fp=open_file(infile,tmp_str,"w"); 
		tmpchn.PutHierarchicalAlnHead(rtf_fp);
	   } else {
		tmpchn.PutPageBreak(rtf_fp);	// print page break between alignments.
		tmpchn.PrintFileName(rtf_fp);	// print filename at top of page.
	   } tmpchn.PutHierarchAligns(rtf_fp);
	 }
	 for(n=1; n <= NumCMA; n++) if(rtn_cma[n]) TotalNilCMSA(rtn_cma[n]); free(rtn_cma);
	 for(n=1; n <= NumHSW; n++) NilHSW(hsw[n]); free(hsw);
	 for(n=0; n < ArgC; n++){ free(ArgV[n]); ArgV[n]=0; } ArgC=0;
	}	// chn_typ destructor called here...
        free(fg_set); free(bg_set); 
       }  // end of if(DisplayHits[q] > 1) // was: end of for(n=Hpt->NumBPPS(); n >= 1; n--) loop...
     }  // end of for(q=1; q <= NumDisplayCMA; q++) loop...
#if 1
     if(!user_display_set){ 	// then create an alignment on the fly.
        for(q=1; q <= NumDisplayCMA; q++) TotalNilCMSA(DisplayedCMA[q]); free(DisplayedCMA);
     }
#endif
     if(debug_fp) fclose(debug_fp);
     if(rtf_fp){  che[1]->PutTailRTF(rtf_fp); fclose(rtf_fp); }
     free(DisplayHits);
}

double  **mcs_typ::GetResEvals(Int4 n)
// Get res_evals and pvals for passing into tmpchn();
{
        assert(n > 0 && n<= Hpt->NumBPPS());
        bpps_typ *tmp_bpps=che[n]->BPPS();
        sst_typ *xsst=tmp_bpps->RtnSST( );
        Int4 s,r,len=tmp_bpps->LenPattern();
        double  **ResEvals,*tmp_pval=che[n]->SubMap( );
        NEWP(ResEvals,len+3,double);
        for(s=1; s <= len; s++){
            NEW(ResEvals[s],nAlpha(AB) +3, double);
            if(xsst[s]==0) continue;
            ResEvals[s][0]=tmp_pval[s];
            for(r=1; r <= nAlpha(AB); r++){
                if(MemSset(r,xsst[s])){ ResEvals[s][r] = tmp_pval[s]; }
            }
        } return ResEvals;
}

