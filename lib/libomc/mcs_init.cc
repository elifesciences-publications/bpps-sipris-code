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
#include "rst_typ.h"
#include "swt_typ.h"

char	**mcs_typ::GetSetRelations(char *Title,Int4 *nGrpsX, Int4 **GrpsX,set_typ **RtnSetX)
//**************** Identify relationships between sets. *******************
{
	Int4	m,n,i,j,f,g,b,r;
	char **RelateSets;
// subset[n1][n2]='<',superset='>',Intersect='+',Disjoint='0',Identical='='.
	set_typ *SetX;
	NEWP(RelateSets,Hpt->NumBPPS() + 5,char);
	NEW(SetX,Hpt->NumBPPS() + 5,set_typ);
        for(n=1; n<= Hpt->NumBPPS(); n++){
	   NEW(RelateSets[n],Hpt->NumBPPS() + 5,char);
	   SetX[n]=MakeSet(Hpt->NumSets()+1); ClearSet(SetX[n]);
	   for(i=1; i <= nGrpsX[n]; i++){ g = GrpsX[n][i]; AddSet(g,SetX[n]); }
	   RelateSets[n][n]='=';
	}
	Int4 ci,cn,cm;
	// set_typ Union=MakeSet(Hpt->NumSets()+1); 
        for(n=1; n < Hpt->NumBPPS(); n++){
	   cn=CardSet(SetX[n]); 
           for(m=n+1; m<= Hpt->NumBPPS(); m++){
		ci=CardInterSet(SetX[n],SetX[m]);
		if(ci == 0){ RelateSets[n][m]='0'; RelateSets[m][n]='0'; continue; } // disjoint.
		cm=CardSet(SetX[m]);
		// fprintf(stderr,"c%d=%d;c%d=%d;ci=%d\n",n,cn,m,cm,ci);
		if(cn == ci){	// n is subset; m is superset.
		    if(cn == cm){ RelateSets[n][m]='='; RelateSets[m][n]='='; continue; } // equal.
		    else { RelateSets[n][m]='<'; RelateSets[m][n]='>'; continue; } // subset/superset
		}
		if(cm == ci){	// m is proper subset; n is proper superset.
		    RelateSets[n][m]='>'; RelateSets[m][n]='<'; continue; // subset/superset
		}
		if(ci > 0){ RelateSets[n][m]='+'; RelateSets[m][n]='+'; continue; } // intersect
		else {
		   if(efp) fprintf(efp,"c%d=%d;c%d=%d;ci=%d\n",n,cn,m,cm,ci);
		   // fprintf(stderr,"cn=%d;cm=%d;ci=%d\n",cn,cm,ci);
		   print_error("mcs_typ::GetSetRelations( ): this should not happen");
		}
		// UnionSet3(SetX[n],SetX[m],Union); Int4 cu=CardSet(Union);
	   }
	}
   if(0){
        fprintf(stdout,"\n %s Set Relationships:\n     ",Title);
        for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(stdout,"%2d ",n); } fprintf(stdout,"\n");
        for(n=1; n <= Hpt->NumBPPS(); n++){
           fprintf(stdout,"%3d: ",n);
           for(m=1; m<= Hpt->NumBPPS(); m++){ fprintf(stdout," %c ",RelateSets[n][m]); }
	   fprintf(stdout,"\n");
	   // NilSet(SetX[n]);
	} fprintf(stdout,"\n");
   }
	// free(SetX);
	*RtnSetX = SetX;
	return RelateSets;
}

void    mcs_typ::PutSetRelations(FILE *fp)
{
	Int4	n,m;
        fprintf(fp,"\n FG Set Relationships:\n     ");
        for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(fp,"%2d ",n); } fprintf(fp,"\n");
        for(n=1; n <= Hpt->NumBPPS(); n++){
           fprintf(fp,"%3d: ",n);
           for(m=1; m<= Hpt->NumBPPS(); m++){ fprintf(fp," %c ",RelateFGs[n][m]); }
	   fprintf(fp,"\n");
	   // NilSet(SetX[n]);
	} fprintf(fp,"\n");
        fprintf(fp,"\n BG Set Relationships:\n     ");
        for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(fp,"%2d ",n); } fprintf(fp,"\n");
        for(n=1; n <= Hpt->NumBPPS(); n++){
           fprintf(fp,"%3d: ",n);
           for(m=1; m<= Hpt->NumBPPS(); m++){ fprintf(fp," %c ",RelateBGs[n][m]); }
	   fprintf(fp,"\n");
	   // NilSet(SetX[n]);
	} fprintf(fp,"\n");
}

void    mcs_typ::GetChnFiles( )
//************** get the chn files for the BPPS analyses. *****************
{
  Int4	i,j,n;
  FILE	*fp;

  //****** 3. create a DisplayCMA files for each analysis ***********************
  const Int4 max_num_TmpCMA=7000;
  cma_typ	TmpCMA[7004];
  Int4		NumTmpCMA;
  NEW(QryCMAs,Hpt->NumBPPS()+ 3,cma_typ);

  for(n=1; n<=Hpt->NumBPPS(); n++){
	if(Hpt->GrpsBG(n,1)==0){	// main set...all in foreground...
	   print_error("BG=(0) not yet implemented");
	} else if(passed_in_cma != 0 && IsTreeHpt){
	   assert(Hpt->NumBPPS() == NumDisplayCMA);
	   // QryCMAs[n]=CopyCMSA(DisplayCMA[n]);
	   fp = tmpfile(); PutCMSA(fp,DisplayCMA[n]);
	   rewind(fp); QryCMAs[n]=ReadCMSA(fp,AB); fclose(fp);
	} else {
	    if(Hpt->nGrpsFG(n) == 0){
		fprintf(stderr,"Tripartition %d category lacks a foreground set\n",n);
		print_error("Fatal hyperpartition syntax error.");
	    }
	    for(i=1,j=0; i <= Hpt->nGrpsFG(n); i++){
		if(Hpt->GrpsFG(n,i) > NumDisplayCMA){
			if(efp) fprintf(efp,"GrpsFG[%d][%d]=%d > %d.\n",
				n,i,Hpt->GrpsFG(n,i),NumDisplayCMA);
			print_error("Fatal: Number of display sets less than number of FG sets."); 
		}
		assert(Hpt->GrpsFG(n,i) > 0);
		j++; TmpCMA[j]=DisplayCMA[Hpt->GrpsFG(n,i)];
	    }
  	    fp = tmpfile(); PutMergedCMSA(fp,j,TmpCMA);
  	    rewind(fp); QryCMAs[n]=ReadCMSA(fp,AB); fclose(fp);
	    // use -nocsq option if want to use the csq provided in sma file (e.g., to rerun pmcBPPS).
	    if(NoCSQ==FALSE){	// then add a consensus sequence.
	       cma_typ ConsCMA = AddConsensusCMSA(QryCMAs[n]);
	       TotalNilCMSA(QryCMAs[n]); QryCMAs[n] = ConsCMA; ConsCMA=0;
	    }
	    RenameCMSA(Hpt->GrpName(n),QryCMAs[n]);
	    for(i=1; i <= Hpt->nGrpsFG(n); i++) TmpCMA[i]=0;
	}
  }

  //****** 4. create a dummy cma file for chn format (not used by anything right now). *****
  BooLean	*skip;
  if(passed_in_cma != 0) assert(NumSeqsCMSA(QryCMAs[1]) > 0);
  else assert(NumSeqsCMSA(QryCMAs[1]) > 1);
  NEW(skip,NumSeqsCMSA(QryCMAs[1])+3,BooLean);
  for(i=2; i<=NumSeqsCMSA(QryCMAs[1]); i++) skip[i]=TRUE;
  fp = tmpfile();
  PutSelectCMSA(fp,skip,QryCMAs[1]);
  rewind(fp); dummyCMA=ReadCMSA(fp,AB); fclose(fp); free(skip);
  // WriteCMSA("dummy.cma",dummyCMA);	// DEBUG.
  NEWP(chn,Hpt->NumBPPS() + 3, chn_typ);
  hsw_typ HSW=0;
  if(passed_in_hsw) HSW=passed_in_hsw;
  //****** 5. create a chn files for the analysis ***********************
  for(n=1; n<=Hpt->NumBPPS(); n++){
	if(NumSeqsCMSA(QryCMAs[n]) > max_num_TmpCMA){
		fprintf(stderr,"n=%d; Hpt->NumBPPS()=%d; NumSeqsCMSA(QryCMAs[n])=%d\n", n,Hpt->NumBPPS(),NumSeqsCMSA(QryCMAs[n]));
		print_error("FATAL: Input seed alignments contain too many sequences");
	}
	TmpCMA[1]=QryCMAs[n];
  	// fprintf(stderr,"n=%d; Hpt->NumBPPS()=%d; NumSeqsCMSA(QryCMAs[n])=%d\n", n,Hpt->NumBPPS(),NumSeqsCMSA(QryCMAs[n]));
	TmpCMA[2]=MainCMA;
	for(i=1,j=2; i <= NumSeqsCMSA(QryCMAs[n]); i++){
		// fprintf(stderr,"n=%d; i=%d; j=%d\n",n,i,j);
		j++; TmpCMA[j]=dummyCMA;
	} NumTmpCMA=j;
  	Int4	ArgC=0;
	char	*ArgV[10];
	ArgV[0]=AllocString("mcBPPS"); ArgV[1]=AllocString(infile); 
	ArgV[2]=AllocString("-Q"); ArgV[3]=AllocString("-concise"); ArgC=4;
	if(efp) fprintf(efp,"************ Analysis #%d... ************\n",n);
	if(HSW==0){
	    chn[n] = new chn_typ(ArgC,ArgV,NumTmpCMA,TmpCMA,200);
	    assert(!chn[n]->OwnsCMAs());
	    HSW=chn[n]->RtnHSW(1);	// Pass Henikoff weights on to other analyses
	} else {
  	    hsw_typ *hsw; NEW(hsw, NumTmpCMA+3, hsw_typ); // Make sure array is long enough..
	    hsw[1]=HSW;
	    chn[n] = new chn_typ(ArgC,ArgV,NumTmpCMA,TmpCMA,200,hsw);
	    assert(!chn[n]->OwnsCMAs()); free(hsw);
	}
#if 1   // Temporary fix for gaps ('-') in query sequence...need better fix later.
        {
	 cma_typ *IN_CMA=chn[n]->GetIN_CMSA();
         e_type fakeE1=FakeSeqCMSA(1,IN_CMA[1]);
         for(Int4 s=1; s<=LenSeq(fakeE1); s++){
           Int4 r=ResSeq(s,fakeE1);
           if(r==AlphaCode('X',AB)){ // set 'X' to 'A'in fake query seq.
             if(efp) fprintf(efp,"WARNING: setting 'X' at position %d in query to 'A'.\n",s);
             r=AlphaCode('A',AB); EqSeq(s,r,fakeE1);
           }
         }
        }
#endif
	for(i=0; i < ArgC; i++) free(ArgV[i]);
  }
}

Int4	mcs_typ::ReadArgFile()
{
	char      str[2003],*strp;
	Int4      line,i,j;
	FILE *tfp = 0;

	if(passed_in_sma){
	  // DisplayCMA=passed_in_sma; 
	  NumDisplayCMA=num_passed_in_sma;
	  NEW(DisplayCMA,NumDisplayCMA +3, cma_typ);
	  for(i=0; i <= NumDisplayCMA; i++) DisplayCMA[i]=passed_in_sma[i];
	} else {
	  tfp = open_file(infile,".sma","r");
	  DisplayCMA=MultiReadCMSA(tfp,&NumDisplayCMA,0,AB); fclose(tfp); tfp=0;
	}
	if(efp) fprintf(efp,"%d alignments in %s.sma\n",NumDisplayCMA,infile);

	// 1. Discard unused display CMA files && reorder consistent with Hyperpartition.
	Int4	NumFound=0;
	BooLean *found;
	cma_typ *DisplayCMA2;
	NEW(DisplayCMA2,Hpt->NumSets()+ 2, cma_typ);
	// Hpt->Put(stderr);	// DEBUG
	assert((Hpt->NumSets()-1) <= NumDisplayCMA);
	NEW(found,Hpt->NumSets() + 2, BooLean);
	for(i=1; i <= NumDisplayCMA; i++){
	   char *NameDisplay = NameCMSA(DisplayCMA[i]);
	   for(j=1; j <= Hpt->NumSets(); j++){	// includes Random?
		char *NameSet=Hpt->ElmntSetName(j);
		if(strcmp(NameDisplay,NameSet) == 0){
			if(found[j]){
				fprintf(stderr,"--> \"%s\" == \"%s\"\n",NameDisplay,NameSet);
				print_error("Fatal: duplicate name in *sma or *arg file");
			} found[j]=TRUE;  NumFound++;
			DisplayCMA2[j]=DisplayCMA[i]; DisplayCMA[i]=0;
		}
	   }
#if 1	   // change the order of the array passed in from calling environment.
	   if(DisplayCMA[i] != 0){ TotalNilCMSA(DisplayCMA[i]); DisplayCMA[i]=0; }
#else
	   if(passed_in_sma == 0 && DisplayCMA[i] != 0){ TotalNilCMSA(DisplayCMA[i]); DisplayCMA[i]=0; }
#endif
	}
	if(NumFound != (Hpt->NumSets()-1)){
	   fprintf(stderr,"NumFound = %d; Hpt->NumSets( ) = %d\n",
			NumFound,Hpt->NumSets());
	   for(j=1; j <= Hpt->NumSets(); j++){	// includes Random.
	     if(!found[j]) fprintf(stderr,"Not in sma file: \"%s\"\n",Hpt->ElmntSetName(j));
	   } print_error("*.sma and *.hpt files are inconsistent");
	} free(found);
#if 1	
	for(i=1; i <= Hpt->NumSets(); i++){	// original array will be freed elsewhere.
	    DisplayCMA[i]=DisplayCMA2[i]; // keep the original array; free DisplayCMA2.
	} NumDisplayCMA=NumFound;	// reset due to discarded files.
	free(DisplayCMA2);
#else
	if(passed_in_sma == 0) free(DisplayCMA); DisplayCMA=DisplayCMA2; NumDisplayCMA = NumFound;
#endif
	if(!passed_in_mcma){
           assert(NumDisplayCMA == num_passed_in_sets);
           MainCMA=AddRandomCMSA(TrueMainCMA,NumRandom);
        } else { MainCMA=passed_in_mcma; }
	assert(nBlksCMSA(MainCMA) == 1); assert(LengthMain==LengthCMSA(1,MainCMA));
	SizeMain= NumSeqsCMSA(MainCMA); SizeTrueMain= NumSeqsCMSA(TrueMainCMA);
	SetSize = SizeMain + 1;
	RandomSet=MakeSet(SetSize); ClearSet(RandomSet);
	for(i=SizeTrueMain+1; i <= SizeMain; i++) AddSet(i,RandomSet);
	if(PartitionRandomly){ RandomlyPartitionSets(); }
	else if(passed_in_sets){ PartitionByInputSetCMSA( ); CheckInputSets( ); }
	else PartitionBySeedAlnCMSA(NumDisplayCMA,DisplayCMA,Hpt->TypeOfSet());
	GetElmntSets( );
	Int4 N=ReadSeedPttrns( );
	// assert(Hpt->NumBPPS()==N);
	if(efp) fprintf(efp,"levels = %d; numBPPS = %d\n",Hpt->NumBPPS(),Hpt->NumBPPS());
	return N;
}

void	mcs_typ::InitDefaults( )
{
	user_display_set=FALSE;
	verbose=TRUE;
	GlobalA0=48.0; GlobalB0=2.0; GlobalRi=0.03;
#if 1	// new: allow more contamination in the root node; fraction of input seq set?
	MiscGlobalA0=42.0; MiscGlobalB0=8.0; MiscGlobalRi=0.3;
	RootGlobalA0=60.0; RootGlobalB0=40.0; RootGlobalRi=0.5;
	// RootGlobalA0=35.0; RootGlobalB0=15.0; RootGlobalRi=0.5;
#else
	MiscGlobalA0=40.0; MiscGlobalB0=10.0; MiscGlobalRi=0.4;
	RootGlobalA0=40.0; RootGlobalB0=10.0; RootGlobalRi=0.4;
#endif
	font_size=6; page_format='P';
	DefaultMaxCol=25; DefaultMinCol=5; // if DefaultMinCol too low == major problems
	MinimumSetSize=40; MinimumLLR=100; // Are always reset (or should be) from calling environment.
	MinimumNatsPerWtSeq=2.0;
	Global_rho=0.0001;
	// Global_rho=0.003;
	ppb_increase=100;
	PrintToggle=20;
	MaxNumberBPPS=1000;
	// SeedPttrnLen=15;
	SeedPttrnLen=25;	// Initial length of the pattern.
	// DefaultContrast=12;
	temperature=300; MinTemperature=50.0; MaxTemperature=3000;
	PutIntermediateFiles=FALSE;
	GlobalN=25;	// default contrast setting.
	MinNats=5;
	ModStart=5;
        IterStart_DF=1; IterEvolve_DF=3; NumRounds_DF=4;
	ColSampleStart_DF=2;
	del_as_random=1;
}

void	mcs_typ::InitFlags( )
{
	SemiConvergedState=0;
	IsTreeHpt=FALSE;
	Evolve=TRUE;
	Unlabeled=FALSE;
	PartitionRandomly=FALSE;
	Verbose=FALSE;
	StrictIndepend=FALSE;
	DidRestoreBest=FALSE;
	SaveSets=TRUE;
	NoSeeds=FALSE;
	PrintEachRTF=FALSE;
	NoCSQ=FALSE;
	NoFailureMode=FALSE;
	SaveBest=FALSE;		// 
}

void	mcs_typ::InitAsNull( )
{
	BestLPR = -99999999.8;
	TotalLPR = -99999999.9;
	SampledColumn=0;
	SampledSet=0;
	Skip=0;
	Iteration=0;
	NthSeqForDisplay=0;	// Use <int>th sequence for numbering rtf display constrast alignment.
	NumHighlighted=0;	// Number of pattern positions to highlight in constrast alignments.
	NumCalls=0;
	DisplayCMA=0;
  	QryCMAs=0;
	RandomSet=0;
	FileID=0;
	MaxNumCol=0;  SFBG=0;
	//======== Strings: ========
	program_name=0;	// what program is calling this object?
	RelateFGs=0; RelateBGs=0;
	//======== file pointers ========
	outfp=0; cfp=0; ifp=0; efp=0; // efp=stderr;
	//======== Sets: ========
	Labeled=0; SetFG=0; SetBG=0;
	//======== CMAs: ========
	MainCMA=0; dummyCMA=0; TrueMainCMA=0;
	//======== Other objects =======
	BestHpt=0;
	che=0; sst=0; chn=0;
	// initialize static arrays to zero to avoid core dumps
	for(Int4 g=0; g < MAX_NUM_ELMENTARY_SETS; g++){
		BestCsq[g]=0;
		SeedCMA[g] =0; 
		MinSeed2CsqScore[g]=0;
		sst_str[g]=0; best_sst[g]=0;
		IsFailedSet[g]=0; IsFailedBPPS[g]=0;
		InitSet[g]=0; BestSet[g]=0; 
		WorstToBest[g]=0;  Index1stBest[g]=0; GrpSet[g]=0;
		SubLPR[g]=0; Map[g]=0;
	}
}

void	mcs_typ::Init(Int4 argc, char *argv[])
// -N=20 -col=30:50 -A=48:2 -Am=40:10 -Ri=0.03 -Rim=0.4 -rho=0.003 -ppb=100 -set
{
	Int4	x,sq,m,n,i,j,g;

	assert(passed_in_cma);
	InitAsNull(); InitDefaults(); InitFlags();
	AB=AlphabetCMSA(passed_in_cma);
	ReadMainArg(argc, argv);

	if(passed_in_hpt) Hpt=passed_in_hpt;
	else Hpt = new hpt_typ(infile);	// create Arg file typ for output...
	{	// MDL adjustment for size of input set
		Int4	NumInternal,NumLeaves,NumNodes=Hpt->NumSets();
		NumInternal=Hpt->NumInternalNodes();
		NumLeaves=NumNodes-NumInternal-2;
		
#if 1	// don't favor internal nodes...
		double d=0.05/(double)(NumNodes-2);
		GlobalRi=d;	// was 0.03.
		MiscGlobalRi=d;	// was 0.3
		d=0.25;	// root.
		RootGlobalRi=d;	// was 0.5
		// random = 0.70 (i.e., the rest) by default.
#else	// favor internal nodes...
		double d=0.05/(double)NumLeaves;
		GlobalRi=d;	// was 0.03.
		d=0.25/(double)(NumInternal);
		MiscGlobalRi=d;	// was 0.3
		d=0.50/2.0;	// root + random.
		RootGlobalRi=d;	// was 0.5
#endif
#if 0	// debug...
		fprintf(stderr,"NumNodes=%d; NumInternal=%d; NumLeaves=%d\n",
			NumNodes,NumInternal,NumLeaves);
		fprintf(stderr,"GlobalRi=%.3g; MiscGlobalRi=%.3g; RootGlobalRi=%.3g\n",
			GlobalRi,MiscGlobalRi,RootGlobalRi);
#endif
	}
	if(efp) Hpt->Put(efp);
	TrueMainCMA = passed_in_cma;	// always passed in...
	NumRandom=ComputeNumRandom(NumSeqsCMSA(TrueMainCMA)); 
	assert(nBlksCMSA(TrueMainCMA) == 1); 
	LengthMain=LengthCMSA(1,TrueMainCMA);
	// TargetMod=(2*LengthMain) + 1; TargetMod=MAXIMUM(Int4,TargetMod,5);
	TargetMod=NumSeqsCMSA(TrueMainCMA)/10; TargetMod=MAXIMUM(Int4,TargetMod,5);
#if 1	// This accomodates lineage-specific input files...
	if(passed_in_sets){
		x=SetN(passed_in_sets[1]);
		n=NumSeqsCMSA(TrueMainCMA);
		m=x-n-1;
		if(m != NumRandom){
		    fprintf(stderr,"SetSize=%d != NumSeqs(%d) + NumRandom(%d) = %d\n",
			x,n,NumRandom,NumRandom+n+1);
		    fprintf(stderr,"Adjusting: SetSize=%d; NumSeqs = %d; NumRandom = %d\n",x,n,m);
		    if(m > NumRandom) NumRandom=m;
		    else print_error("mcs_typ::Init() Input set size is too small");
		}
	}
#endif
	Hpt->ChangeNumRandom(NumRandom);
	if(efp) Hpt->PutHyperPartition(efp);

	// 2a. Read in arguments and seed pattern strings for each run in <infile>.hpt
	ReadArgFile();	// creates MainCMA.

        for(g=1; g<= Hpt->NumSets(); g++){	
		IsFailedSet[g] = FALSE;		// include all sets initially...
		if(strcmp("Random",Hpt->ElmntSetName(g)) == 0){
			assert(g == Hpt->NumSets());
		} 
	}

	// 2c. Read main <file>.chn and create input cma files (chn_typ).
	GetChnFiles( );

	// 3. Create che objects, one for each search based on arguments + seed patterns.
	if(Hpt->NumBPPS() < 1) print_error("too few categories of constraints");
	if(Hpt->NumBPPS() > MaxNumberBPPS) print_error("too many categories of constraints");
	NEWP(che,MaxNumberBPPS+4,che_typ);
	NEWPP(sst,MaxNumberBPPS+4,sst_typ);
	NEW(SFBG,MaxNumberBPPS+4,char);
	NEW(MaxNumCol,MaxNumberBPPS+4,Int4);
	for(Int4 n=1; n <= MaxNumberBPPS; n++){ MaxNumCol[n]=DefaultMaxCol; }

#if 0	// make sure that first set is the main set...
	if(!Hpt->IsMiscSetForCol(1,1)) print_error("FATAL: FD-table format error");
#endif
	// 3a. Get seed patterns for each run using pattern strings. 
	// 3b. Get legal sets for each run (these will change during the search).
	// Need to initialize legal sets to be subsets of positions higher up the tree.
	for(Int4 n=1; n <= Hpt->NumBPPS(); n++){
		SetUpNthSrch(n, Hpt->nArg(n),Hpt->Argv(n)); IsFailedBPPS[n]=FALSE; 
	}

	// 3c. call che constructors for each run.
	// Use both telescoping & total set to determine subLPR ???
	//**************** Remove sequences that don't belong. *******************
	BooLean	*InFgBg;
	NEW(InFgBg, Hpt->NumSets()+3,BooLean);
	for(n=1; n<= Hpt->NumBPPS(); n++){	
	  for(g=1; g <= Hpt->NumSets(); g++) InFgBg[g]=FALSE;
	  for(i=1; i <= Hpt->nGrpsFG(n); i++){ g = Hpt->GrpsFG(n,i); InFgBg[g]=TRUE; }
	  for(i=1; i <= Hpt->nGrpsBG(n); i++){ g = Hpt->GrpsBG(n,i); InFgBg[g]=TRUE; }
	  for(g=1; g <= Hpt->NumSets(); g++){
	    if(!InFgBg[g]){	// i.e., if g is an omitted subgroup.
		for(sq=1; sq <= NumSeqsCMSA(MainCMA); sq++){
		  if(MemberSet(sq,InitSet[g])){
			BooLean rtn=che[n]->RemoveSeq(sq);
		  }
		}
	    }
	  }
	} free(InFgBg);
	//**************** Label sequences to be kept fixed. *******************
	Labeled=MakeSet(SetSize); ClearSet(Labeled); m=0;
        for(n=1; n<= Hpt->NumBPPS(); n++){
	   for(sq=1; sq <= SizeMain; sq++){
		if(che[n]->MemberGold(sq)){
		  if(!MemberSet(sq,Labeled)){ AddSet(sq,Labeled); m++; }
		}
	   }
	}
	if(efp) fprintf(efp,"%d labeled sequences found\n",m);
	// Add best matches within specific groups to Labeled set.
	Int4 Sq0;
  	// for(m=0,g=1; g<=Hpt->NumSets(); g++)
  	for(m=0,g=1; g < Hpt->NumSets(); g++) // don't label any sequences in the reject set...
	{
	    for(j=Index1stBest[g]; (Sq0=WorstToBest[g][j]) != 0; j++){
		if(MemberSet(Sq0,InitSet[g])){
			if(!MemberSet(Sq0,Labeled)){ AddSet(Sq0,Labeled); m++; }
		}
	    }
	}
	if(efp) fprintf(efp,"%d additional sequences labeled based on log-odds scores.\n",m);
	//**************** Create HyperPartition. *******************
	// for(g=1; g<= Hpt->NumSets(); g++){ HyperPartition[g] = Hpt->RtnHyperPartition(g); }
     	//**************** identify sequence subgroups. *******************
	for(g=1; g<= Hpt->NumSets(); g++){	
	   GrpSet[g]=MakeSet(SetSize); ClearSet(GrpSet[g]);
	   CopySet(GrpSet[g],InitSet[g]);
	}
	for(x=Hpt->NumSets(),sq=SizeTrueMain+1; sq <= SizeMain; sq++){
		AddSet(sq,Labeled);	// add random sequences to Labeled set.
	}
	RelateFGs=GetSetRelations("FG",Hpt->nGrpsFG(),Hpt->GrpsFG(),&SetFG);
	RelateBGs=GetSetRelations("BG",Hpt->nGrpsBG(),Hpt->GrpsBG(),&SetBG);

        //**************** Put sequences in their starting partitions. *******************
	// need to make sure that labeled sequences are in all the right partitions...
	// also label BG seqs in same hyperpartition:
     	for(g=1; g<= Hpt->NumSets(); g++){	
	  for(sq=1; sq <= SizeMain; sq++){
	   if(MemberSet(sq,GrpSet[g])){		// this sequence is in subgroup g.
             for(n=1; n<= Hpt->NumBPPS(); n++){
		switch (Hpt->RtnHyperPartition(g,n)){
		  case '+': 	
		    assert(!che[n]->RemovedSeq(sq));
	      	    if(!che[n]->MemberFG(sq)){ assert(che[n]->ChngPartition(sq)=='+'); }
		  break;
		  case '-': 
		    assert(!che[n]->RemovedSeq(sq));
	      	    if(!che[n]->MemberBG(sq)){ assert(che[n]->ChngPartition(sq)=='-'); }
		  break;
		  case 'o':
		   	assert(che[n]->RemovedSeq(sq));	// this was done above.
		  break;
		  default: print_error("mcs_typ Init( ) error");
		  break;
		}
	     }
           }
	  }
	}
	FILE *sfp=0;
	for(n=1; n<= Hpt->NumBPPS(); n++){	
	  if(sst_str[n] == 0){		// then find a seed pattern based on above partition.
		if(sfp==0 && PutIntermediateFiles){ sfp=open_file(infile,".seeds","w"); }
		sst_str[n]=FindSeedPattern(n); // resets the pattern to the seed pattern.
		if(PutIntermediateFiles) fprintf(sfp,"%d: %s\n",n,sst_str[n]);
	  }
     	} if(sfp) fclose(sfp);
	//********* For storing best ***********
	SaveBest=FALSE;		// 
	BestLPR=-999999999999.9;
     	for(g=1; g<= Hpt->NumSets(); g++){ BestSet[g]=MakeSet(SetSize); ClearSet(BestSet[g]); }
	Int4	Length=LenSeq(che[1]->KeyE());
	for(n=1; n<= Hpt->NumBPPS(); n++){ NEW(best_sst[n],Length+2,sst_typ); }
	//********* For storing best ***********
	// outfp = open_file(infile,".out","w");
	// PutHyperPartition(outfp); fflush(outfp);
#if 0	// debug...
     	// for(n=1; n<= Hpt->NumBPPS(); n++){ che[n]->BPPS()->PutParameters(stderr); }
	PutHyperPartition(stderr); exit(1);
#endif
}

/**************************** Global Variables ******************************/
Int4	mcs_typ::SetUpNthSrch(Int4 Level, Int4 argc,char *argv[])
{ 
	Int4	arg;
	// BooLean	UseGlobalSqWts=FALSE; // sets StartAlpha=1; which causes problems...
	BooLean	UseGlobalSqWts=TRUE;
	// BooLean	concise=FALSE;
	BooLean	concise=TRUE;
	char	compare='U';
	double	MinKeyFrq=0.5; // ,MaxGapFrq=0.5;
	Int4	min_nats=5;
	double	fract_ignored[20];
	double	A0=GlobalA0,B0=GlobalB0;
	char	Mode='R';	// can only use random order so that all sequences are in FG.
	// Int4	Contrast=12;
	Int4	Contrast=GlobalN,contrast=-1;
	char    sets_mode='L';
	// double	LnRho= 0.6931471805599452862;	// -log(0.5);
	// double	rho=0.1;
	// double	rho=0.003;
	double	rho=Global_rho;
	double	PriorRi=GlobalRi;
	Int4	min_num_col=DefaultMinCol;
	Int4	max_num_col=DefaultMaxCol;
	char    *pttrn_str=0;
	Int4	n=Level;
	SFBG[n]='B';

	if(del_as_random) Mode='r';	// this causes che_typ to treat deletions '-' as random background.
	fract_ignored[1]=0.0; fract_ignored[2]=0.0;

	char	type_node='L';
   {   //********************** for internal use. *************************
	if(Level == 1){
		type_node='R';
		A0=RootGlobalA0; B0=RootGlobalB0; PriorRi=RootGlobalRi; 
	} else for(Int4 row=1; row <= Hpt->NumSets(); row++){
		char cell=Hpt->Cell(row,Level);
		if(cell == '+'){
		    if(Hpt->TypeOfSet(row) == '?'){  // first positive row...
			type_node='M';
			A0=MiscGlobalA0; B0=MiscGlobalB0; PriorRi=MiscGlobalRi; 
		    } break;
		}	// for FGs with a miscellaneous set use more liberal settings.
	}
	for(arg = 0; arg < argc; arg++){
	   if(efp) fprintf(stderr,"argv[%d] = %s\n",arg,argv[arg]);
	   if(argv[arg][0] != '-') print_error(MCS_USAGE_START);
	   switch(argv[arg][1]) {
             case 'A': if(sscanf(argv[arg],"-A%lf:%lf",&A0,&B0)==2){
			if(A0 <= 0.0 || B0 <= 0.0) print_error(MCS_USAGE_START); 
			argv[arg][1] = ' '; 
		       } else print_error(MCS_USAGE_START); 
		break;
	     case 'B':
		if(sscanf(argv[arg],"-B=%c",&SFBG[n])==1){
		  if(!(SFBG[n]=='M' || SFBG[n]=='B' || SFBG[n]=='A')) print_error(MCS_USAGE_START);
		} else print_error(MCS_USAGE_START); 
		  argv[arg][1] = ' '; 
		break;
	     case 'c': 
              if(sscanf(argv[arg],"-col=%d:%d",&min_num_col,&max_num_col)==2){
			if(min_num_col < 2 || min_num_col > max_num_col){
				fprintf(stderr,"Min(%d)/Max(%d) number columns out of range\n",
					min_num_col,max_num_col);
				print_error(MCS_USAGE_START);
			} argv[arg][1] = ' '; 
	      } else if(strcmp("-concise",argv[arg]) == 0){
                        concise=TRUE;
	      } else if(sscanf(argv[arg],"-compare=%c",&compare)==1){
                        if(!isupper(compare)) print_error(MCS_USAGE_START);
		  	argv[arg][1] = ' '; 
	      } else {
		// NEW method based on fraction of poorest seqs to ignore.
		// fract_ignored=RealOption(argv[arg],'c',0.0,0.5000001,MCS_USAGE_START);
		// Int4    ParseReals(char *str, double *values, const char *msg);
		// for getting separate values for each partitions
		if(argv[arg][2] != '=') print_error(MCS_USAGE_START);
		Int4 n = ParseReals(argv[arg] + 3,fract_ignored,MCS_USAGE_START);
		if(n != 2) print_error(MCS_USAGE_START);
		fract_ignored[2]=fract_ignored[1];
		fract_ignored[1]=fract_ignored[0];
		if(fract_ignored[2] < 0.0 || fract_ignored[2] > 0.9) print_error(MCS_USAGE_START);
		if(fract_ignored[1] < 0.0 || fract_ignored[1] > 0.9) print_error(MCS_USAGE_START);
                argv[arg][1] = ' '; 
	      } break;
	     case 'g': 
	      if(strcmp("-global",argv[arg]) == 0) UseGlobalSqWts=TRUE;
	      else print_error(MCS_USAGE_START);
		break;
	     case 'm': 
		min_nats=IntOption(argv[arg],'m',0,1000,MCS_USAGE_START); 
                argv[arg][1] = ' '; 
		break;
	     case 'M': MinKeyFrq=RealOption(argv[arg],'M',0.0,1.0,MCS_USAGE_START);
                argv[arg][1] = ' '; break;
	     case 'N': 
             	if(sscanf(argv[arg],"-N=%d",&contrast)==1){
			if(contrast < 1) print_error(MCS_USAGE_START);
			argv[arg][1] = ' ';
                } else print_error(MCS_USAGE_START); 
		break;
             case 'P':
                if(argv[arg][2] == '=' && argv[arg][3] == 0) continue; // skip "-P= " strings.
                if(argv[arg][2] == '=' && isalpha(argv[arg][3])){
                        // argv[arg][1] = ' '; 
			pttrn_str=AllocString(argv[arg]+3);
                } else print_error(MCS_USAGE_START);
                break;
             case 'R':
		if(sscanf(argv[arg],"-Ri=%lf",&PriorRi)==1){
		   if(PriorRi <= 0.0 || PriorRi >= 1.0){
		      fprintf(stderr,"Ri (%.2f) out of range\n",PriorRi);
		      print_error(MCS_USAGE_START);
		   }
		   argv[arg][1] = ' ';
                } else print_error(MCS_USAGE_START); 
		break;
	     case 'r': 
	        if(sscanf(argv[arg],"-rho=%lf",&rho)==1){
		   if(rho <= 0.0 || rho > 0.5){
		        fprintf(stderr,"rho (%.2f) out of range\n",rho);
			print_error(MCS_USAGE_START);
		   }
		   argv[arg][1] = ' '; 
                } else print_error(MCS_USAGE_START); 
		break;
	     case 's': 
	        if(sscanf(argv[arg],"-sets=%c",&sets_mode)==1){
		   if(sets_mode != 'R' && sets_mode != 'T' && sets_mode != 'S' && sets_mode != 'M' 
			&& sets_mode != 'D' && sets_mode != 'O' && sets_mode != 'L'){
				print_error(MCS_USAGE_START);
		   }
		   // don't need to pass sets along to chn_typ, as it is not used within chn_pps
		} else print_error(MCS_USAGE_START);
		argv[arg][1] = ' '; break;
	     case 0: print_error(MCS_USAGE_START); break;
	     case ' ': break;	// ignore these...
	     default: 
		fprintf(stderr,"illegal input option (%c)\n",argv[arg][1]);
		print_error(MCS_USAGE_START); break; // do nothing.
	   }   // end of switch scope
	}	// end of argument loop.
GlobalN=DefaultMaxCol;
	if(contrast != -1) Contrast=contrast;
	else if(GlobalN > 0) Contrast=GlobalN; 
	else Contrast=max_num_col;
// fprintf(stderr,"Contrast=%d\n",Contrast);
    } // End of if(Hpt->RtnMode() == 'I' else 
	if(Level == 0) return 0;

	cma_typ *IN_CMA=chn[Level]->GetIN_CMSA();
        Int4 j,ri,ris;
	e_type Query=FakeSeqCMSA(1,IN_CMA[1]);
	set_mode=sets_mode;
	if(sst_str[n] && !concise) fprintf(stderr,"pattern %d: '%s'\n",n,sst_str[n]);
        if(sst_str[n]==0 && pttrn_str != 0) sst_str[n]=pttrn_str; // provide seed pattern.
	// should order the FG and BG sequences to match the input ordering.
	// Then should look for seed patterns.
// fprintf(stderr,"Creating che_typ for set %d\n",n);
#if 0	// pass in the sq_prior...
	double	sq_prior=this->GetPriorRi(n,0);
	che[n] = new che_typ(sst_str[n],chn[n],!concise,Mode,type_node,UseGlobalSqWts,sq_prior);
#elif 1	// set up above...
	// fprintf(stderr,"%d: PriorRi=%.3g\n",PriorRi);
	che[n] = new che_typ(sst_str[n],chn[n],!concise,Mode,type_node,UseGlobalSqWts,PriorRi);
     	// { che[n]->BPPS()->PutParameters(stderr); }
#else	// use old priorRi...
	che[n] = new che_typ(sst_str[n],chn[n],!concise,Mode,type_node,UseGlobalSqWts);
#endif
	e_type keyE=che[n]->KeyE( ); NEWP(sst[n],LenSeq(keyE)+3,sst_typ);
	bpps_typ *pps=che[n]->BPPS(); pps->UpdateSST(sst[n],keyE);
	che[n]->BeQuiet( );	// NOTE: !concise==verbose
	che[n]->SetMinNumColumns(min_num_col);
	MaxNumCol[n]=max_num_col;
#if 1
	che[n]->SetMaxNumColumns(LenSeq(Query)+2);
#else
	che[n]->SetMaxNumColumns(MaxNumCol[n]);
#endif
	if(Contrast > 0) che[n]->SetContrast(Contrast);
	return 1;
}

Int4	mcs_typ::ReadSeedPttrns( )
//************* read in seed patterns from <infile>.sp *************
{
	Int4	i,f,b;
	for(i=1; i <= Hpt->NumBPPS(); i++){
	//****************** Read partition information ***************
	   assert(Hpt->GrpName(i)!=0);
	   // check to make sure FG and BG don't overlap; later need to do within arg_typ only. 
	   for(f=1; f <= Hpt->nGrpsFG(i); f++){
	     for(b=1; b <= Hpt->nGrpsBG(i); b++){
		if(Hpt->GrpsFG(i,f)==Hpt->GrpsBG(i,b)){
		   fprintf(stderr,"Analysis #%d: group %d assigned to both FG & BG sets.\n",
				i,Hpt->GrpsFG(i,f));
		   print_error("Fatal: FG & BG set overlap disallowed.");
		}
	     }
	   }
	   if(NoSeeds) sst_str[i]=0; else sst_str[i]=Hpt->sst_str(i);
	} return Hpt->NumBPPS();
}

char	*mcs_typ::FindSeedPattern(Int4 Level)
{
	Int4 i,j;
	// char *pttrn=0;

	Int4 Length=LenSeq(che[Level]->KeyE());
	sst_typ *sst_Best; NEW(sst_Best, Length +4, sst_typ);
	dh_type dH = dheap(Length+5,4);
	e_type  cE1=che[Level]->KeyE();
	double	*FG=0,*BG=0;
	UInt4	wtFactor=che[Level]->GetWtFactor();
	NEW(FG, Length +4, double); NEW(BG, Length +4, double);
	for(j=1; j <= Length; j++){
	   unsigned char r; // ,r1=ResSeq(j,cE1);
	   sst_typ Bsst = 0;
	   double score,best_score=-99999.9,best_n,best_d;
	   double n,d;
           for(Int4 s=1; sst[Level][j][s]; s++){
		sst_typ xsst = sst[Level][j][s];
	        double	r1_f=0,r1_n=0,r2_f=0,r2_n=0;
		for(r=1; r <= nAlpha(AB); r++){
		   if(MemSset(r,xsst)){
			r1_f += che[Level]->GetResWtFG(j,r);
			r2_f += che[Level]->GetResWtBG(j,r);
		   } else {
			r1_n += che[Level]->GetResWtFG(j,r);
			r2_n += che[Level]->GetResWtBG(j,r);
		   }
		}
		// double p,q;
		n=((double)(r1_f+1*wtFactor)/(double)(r1_f+r1_n + 2*wtFactor));
		d=((double)(r2_f+1*wtFactor)/(double)(r2_f+r2_n + 2*wtFactor)); 
		score = n*log((n)/d) + (1-n)*log((1-n)/(1-d));
		if(score > best_score)
			{ best_score=score; Bsst=xsst; best_n=n; best_d=d; }
	   }
	   if(best_score > 0.0) {
	      sst_Best[j]=Bsst; insrtHeap(j,(keytyp)-best_score,dH);
	      FG[j]=best_n; BG[j]=best_d;
	   }
	}
	char tmp_str[3003];
	char str0[100];
	tmp_str[0]=0;
	// 4. Remove all columns and replace with new patterns
	Int4 MinCol=che[Level]->RtnMinNumColumns();
	che[Level]->SetMinNumColumns(0);
	// remove all patterns... 
	for(j=1; j <= Length; j++){ che[Level]->RemoveColumn(j); }

	// che[Level]->SpeakUp();
	che[Level]->BeQuiet( );
	tmp_str[3000]=0;
	for(i=0; !emptyHeap(dH); ){
		// double score=-(double)minkeyHeap(dH);
		j=delminHeap(dH); 
		if(che[Level]->AddColumn(j,sst_Best[j])) {
		  i++;
		  if(i > 1) strncat(tmp_str,",",200);
		  PutSST(str0,sst_Best[j],AB);  // copies pattern to str.
		  strncat(tmp_str,str0,22);
		  sprintf(str0,"%d",j);
		  strncat(tmp_str,str0,22);
		  assert(tmp_str[3000]==0);  // don't over-extend string

		}
#if 0		// turn this on to see pattern info...
		   fprintf(stderr,"%d: sst[%d]: ",i,j);
                   PutSST(stderr,sst_Best[j],AB);
                   if(FG) fprintf(stderr,"%d (%.3f)(%.3f/%.3f)\n",j,score,FG[j],BG[j]);
                   else fprintf(stderr,"%d (%.3f)\n",j,score);
#endif
		if(i >= SeedPttrnLen){ while(!emptyHeap(dH)) delminHeap(dH); }
	}

	che[Level]->SetMinNumColumns(MinCol);
	che[Level]->BeQuiet( );
	if(FG) free(FG); if(BG) free(BG); 
	free(sst_Best);
	Nildheap(dH);
	// fprintf(stderr,"%d: %s\n",Level,tmp_str);
	return AllocString(tmp_str);
}


