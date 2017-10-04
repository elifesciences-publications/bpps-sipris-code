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
#include "blosum62.h"

#define USAGE_START "Usage: omcBPPS <prefix> [options] \n\
     Input: <prefix>.mma and <prefix>.hsw as input\n\
           Note if <prefix>.hsw is not provided it will be created\n\
     Output: \n\
	   Checkpoint files: <prefix>_chk.sets, <prefix>_chk.sma, <prefix>_chk.hpt, <prefix>_chk.seed\n\
	   Editable display sequence file: <prefix>_usr.sma\n\
	   Hyperpartition: <prefix>_new.hpt\n\
	   Node sequence sets (in cma format): <prefix>_new.mma\n\
	   Hierarchy set and category summary: <prefix>_bst.out\n\
	   Sequence patterns for each category: <prefix>.pttrns\n\
    -A=<int1>:<int2> Attach node <int1> to its sibling node <int2>\n\
    -alp=<real>:<real>:<int> AddLeaf( ) parameters <real> = MinFrq, merge pcut: <int> = MinClique\n\
                     (default: -alp=0.4:0.0001:5)\n\
    -usr            use _usr.sma file with a checkpoint file and -run=R option\n\
    -D=<int>        Delete node <int> \n\
    -del            Don't treat deletions as random background\n\
    -evolve         Allow consensus sequences to evolve from the start\n\
    -flatten        Start sampler from a 'flattened' checkpoint file \n\
    -Focus=<int1>:<int2>  focus on sampling a new child node from the <int1>th node\n\
                    using the <int2>th sequence in the *.mma input file as a seed.\n\
                    Note that the <int2>th sequence must be assigned to the <int1>th node\n\
    -focus=<int>    focus sampling only on the subtree associated with the <int>th node\n\
                    (nodes will not be added or removed from other parts of the hierarchy)\n\
                    (this option requires *.hpt, *.sma and *.sets input files)\n\
    -I=<int1>:<int2>[,<int>]  Insert a new node between <int1> and <int2> list\n\
    -Lineage=<int>  Output assigned sequence sets for the lineage to <int>th node\n\
    -LLR=<int>      reset LLR parameters (default: 200 nats)\n\
    -log	    create a .log file of operations performed (with table)\n\
    -maxcol=<int>   Reset the maximum number of pattern positions per node (default: 25)\n\
    -maxdepth=<int> Set the maximum depth of the hierarchy; must be > 0 (default: 2)\n\
    -minsize=<int>  reset minimum set size (default: 50 squences)\n\
    -minnats=<real> reset minimum nats per (weighted) sequence (default: 10.0)\n\
    -names=<cmafile> Name subgroups based on this seed alignment\n\
    -Nth=<int>      Use <int>th sequence for numbering display alignment\n\
    -N=<int>        Maximum number of nodes in the tree (default: 1000)\n\
    -random=<int>   start with a random, optimized hpt of <int> > 2 nodes\n\
    -rtf            Output contrast alignments in rich txt format at end of run\n\
    -rtfN=<int>     Number of columns to highlight in rtf output file\n\
    -rtfP=<char>    Specify rtf page format (default: 'P')\n\
                    'l': 8.5\" x 11\" landscape\n\
                    'L': 11\" x 17\" landscape\n\
                    'p': 8.5\" x 11\" profile\n\
                    'P': 11\" x 17\" profile\n\
    -rtfF=<int>     Specify font size (range 4-24; default: 6)\n\
    -run=<char>     run optimization or output operations using checkpoint files\n\
                    'o': Optimize with column evolution and then run 'S' option\n\
                    'O': Optimize and then run 'S' option\n\
                    'c': Optimize columns \n\
                    'C': output checkpoint files only.\n\
                    'P': Call omc_typ->Put(TRUE)\n\
                    'R': output contrast alignment\n\
                    'S': output contributions to LLR\n\
    -run=R<int1>:<int2>=<name> \n\
		    output a contrast alignment for file <name> from position <int1> to <int2> only\n\
    -seed=<int>     provide seed for random number generator\n\
    -Seed           look for <filename>.seed input file\n\
    -strict         Impose strict pattern independence between categories\n\
    -test=<char>    Subroutine test mode\n\
                    <char> == 'k': find cross conserved residue patterns\n\
    -verbose        Generate lots of diagnostic files and stderr output\n\
    -XC=<int>       Look for cross conserved residues for the <int>th set in hpt.\n\
    -xc=<int>:<char>:<int> Set cross-conserved parameters: <heapsize>:<mode>:<weight>\n\
                      mode = X(cross-conserved),L(lineage), B(both) or N(none)\n\
                      weight = 1..200 (default: -xc=20:B:85)\n\
    -x              dummy\n\n"

#define PUBLIC_START "   Input: <prefix>.mma (multi-alignment in cma format)\n\
          To convert from fasta to cma format use the fa2cma program\n\
   Output: \n\
	<prefix>.cmd: command line file.\n\
	<prefix>.chk: checkpoint file (needed for step 2).\n\
	<prefix>_usr.sma: editable display sequence file.\n\
	<prefix>_bst.out: optimal hierarchy with summary of properties.\n\
	<prefix>_aln.rtf: MS-Word-readable contrast alignments in rich text format.\n\
   Options: \n\
    -usr            Use *_usr.sma input display file with -run=R option (requires checkpoint file)\n\
    -focus=<int>    Focus only on the subtree associated with the <int>th node.\n\
                    Nodes will not be added or removed from other parts of the hierarchy.\n\
                    This option requires a checkpoint *.chk file.\n\
    -maxcol=<int>   Set the maximum number of pattern positions per node (default: 25)\n\
    -maxdepth=<int> Set the maximum depth of the hierarchy; must be > 0 (default: 2)\n\
    -minsize=<int>  Set the minimum number of sequences required to create a node (default: 50 sequences)\n\
    -minnats=<real> Set minimum information content (in nats) per (downweighted) sequence (default: 10.0)\n\
    -rtfP=<char>    Specify rtf page format (default: 'P')\n\
                    'l': 8.5\" x 11\" landscape\n\
                    'L': 11\" x 17\" landscape\n\
                    'p': 8.5\" x 11\" profile\n\
                    'P': 11\" x 17\" profile\n\
    -rtfF=<int>     Specify rtf font size (range 4-24; default: 6)\n\
    -run=<char>     Run one of the following optimization or output operations (requires checkpoint file)\n\
                    'o': Optimize over sequences and patterns (without changing the hierarchy)\n\
                    'c': Optimize over patterns only\n\
                    'R': output contrast alignment\n\
    -seed=<int>     Provide a seed for the random number generator (ignored if *.chk file is present)\n\
  Note: If phylum and kingdom information is added to sequence deflines (in the input multiple alignment)\n\
      then contrast alignments will include one sequence from each phylum. Added information has\n\
      the following syntax: \">seq_id {<phylum(char)>}...\", where 'char' specifies the kingdom.\n\
      Kingdom codes: metazoan, M; fungi, F; protozoan, E; plants, V; bacteria, B; archaea, A.\n\
      Example: >RHOA_HUMAN {<chordata(M)>}Full=Transforming protein RhoA.\n\
  Reference:\n\
    Neuwald, A.F. 2014. A Bayesian sampler for optimization of protein domain hierarchies. \n\
      Journal of Computational Biology 21(3):269-286.\n\n"

void	omc_typ::PrintError(char *prgm_name)
{
	if(SecretCode > 0){
		fprintf(stderr,"Usage: %s %d <prefix> [options] \n",prgm_name,(int)SecretCode);
		print_error(PUBLIC_START);
	} else { print_error(USAGE_START); }
}

void	omc_typ::InitDefaults()
{
	MaxNumNodes=1000;
	// MaxDepthHierarchy=5;	// class -> subclass -> superfamily -> family -> subfamily
	MaxDepthHierarchy=3;	// class -> superfamily -> family (-> subfamily)
	font_size=6; page_format='P';
	Evolve=FALSE;
	XCuseX=TRUE; XCuseL=TRUE; XC_mode='B';
	XC_size=20; XC_wt=0.85;
	verbose=FALSE;
	flatten_hpt=FALSE;
	SaveSets=TRUE;
	DefaultMaxCol=0;
	random_start=0;
	MinimumNatsPerSeq=10.0;
	NumHighlighted=0;	// set by mcs automatically..
	del_as_random=1;
	AB=MkAlphabet(AMINO_ACIDS,GBLAST_BLOSUM62,SREL26_BLSM62);
	NEW(aafreq, nAlpha(AB)+3,double);
	for(Int4 r=0; r <= nAlpha(AB); r++) aafreq[r]=blosum62freq[r];	
	// make sure that order of AA residue freqs is consistent!
	// *MinFrq*      MaxGapFrq  *pcut*  Exact_pcut  *MinClique*  sets_mode
	// ^matching and 'X' freqs. 	    ^signif. Correlation     ^allowed patterns
	//                          ^merge_patterns     ^Min size for merging
	MP = new mp_type(0.4,0.5,0.0001,0.001,5,'L'); 	// Original
	// MP = new mp_type(0.2,0.5,0.0001,0.001,4,'L'); 
	// MP = new mp_type(0.2,0.5,0.002,0.001,4,'L'); 
	// MP = new mp_type(0.2,0.5,0.002,0.001,5,'L'); 
	// MP = new mp_type(0.2,0.5,0.002,0.001,4,'L'); 
	stringency=50;	// == minimum set size.
	MvNode=ToNode=0;
	ForcedP=0; NumForced=0; ForcedNode=0;
	DelNode=0;
}

void	omc_typ::InitAsNull()
{
	time1=time(NULL);
	NEWP(Argv,1000,char);
	outfilename=0;
	test_mode=0;
	Hpt=0;
	FocusSeq=0; FocusNode=0; stable=0;
	OutPutRTF=FALSE;
	use_usr_sma=FALSE;
	StrictIndepend=FALSE;
	Iteration=0; Set=0; SetFG=0; SetBG=0;
	NthSeqForDisplay=0;
	InitMode=' ';
	KeySetName=0;
	NameSeedAln=0;
	KeyStart=KeyEnd=0;
	XC_line=0;
	hsw=0;
	mcs=0;
	swt=0;
	Argc=0;
	PrintMode=0;
	// SMA=0; 
	DisplaySet=0; 
	best_mcs=0;
	goodHG=badHG=testHG=triesHG=0;
	logfp=0; LastNumNodes=0; status=0;
	dfsHG=0;
	preHG=0;
}

void	omc_typ::Init(Int4 argc, char *argv[])
{
	Int4	i,j,n,arg,x,y;
	double	d,D;
	UInt4	seed=18364592;
	Int4	seed_arg=0;
	BooLean	usingSeedOption=FALSE;
	Best_LPR=-999999999999.0;
	char	c;

	if(argc < 2){ PrintLicenseStatement();  PrintError(argv[0]); }
	TurnOffLicenseStatement();
	InitAsNull(); InitDefaults();
	goodHG=Histogram("Successful pattern LLRs",-100,500,25.0);
	badHG=Histogram("Failed pattern LLRs",-100,500,25.0);
	testHG=Histogram("Difference between heuristic and true LLRs",-1000,2000,100.0);
	triesHG=Histogram("Number of tries before success in AddLeaf()",0,20,1.0);
        for(arg = 2; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'A': 
		if(sscanf(argv[arg],"-A=%d:%d",&MvNode,&ToNode)!=2) PrintError(argv[0]);
		PrintMode='T'; 	// i.e., don't rename the nodes...
		break;
             case 'a': 
		if(sscanf(argv[arg],"-alp=%g:%g:%d",&D,&d,&x)!=3) PrintError(argv[0]);
		else if(d < 0.0 || d > 1.0 | D > 1.0 | D <= 0.0 | x < 3 | x > 10){
		   	PrintError(argv[0]);
		} else { delete MP; MP = new mp_type(d,0.5,D,0.001,x,'L'); }
		break;
             case 'D': 
		if(sscanf(argv[arg],"-D=%d",&DelNode)!=1) PrintError(argv[0]);
		PrintMode='T'; 	// i.e., don't rename the nodes...
		break;
             case 'd': 
		if(strcmp("-del",argv[arg]) == 0){ del_as_random=0; }
		else PrintError(argv[0]); break;
             case 'e': 
		if(strcmp("-evolve",argv[arg]) != 0) PrintError(argv[0]);
		Evolve=TRUE;
		break;
             case 'F': 
		   if(sscanf(argv[arg],"-Focus=%d:%d",&x,&y) != 2) PrintError(argv[0]);
		   else if(x < 1 || y < 1) PrintError(argv[0]);
		   else { FocusNode=x; FocusSeq=y; test_mode='F'; }
		break;
             case 'f': 
		   if(strcmp("-flatten",argv[arg]) == 0){
			flatten_hpt=TRUE;
		   } else if(sscanf(argv[arg],"-focus=%d",&x) != 1) PrintError(argv[0]);
		   else if(x < 1) PrintError(argv[0]);
		   else { FocusNode=x; FocusSeq=0; }
		break;
             case 'I': 
		if(sscanf(argv[arg],"-I=%d:%s",&ForcedP,str)!=2) PrintError(argv[0]);
		NEW(ForcedNode,strlen(str)+3,Int4);
		NumForced=ParseIntegers(str,ForcedNode, "-I option input error");
		PrintMode='T'; 	// i.e., don't rename the nodes...
		break;
             case 'l': 
		if(strcmp("-log",argv[arg]) != 0) PrintError(argv[0]);
		logfp=open_file(argv[1],".log","w");
		break;
             case 'L': 
#if 1
		   if(sscanf(argv[arg],"-Lineage=%d",&x) == 1){
		      print_error("omcBPPS -Linaeage option not yet implemented\n");
		      if(x < 1) PrintError(argv[0]);
		      else { FocusNode=x; PrintMode='s'; }
		   } 
#endif
		   if(sscanf(argv[arg],"-LLR=%d",&x) != 1) PrintError(argv[0]);
		   else if(x < 1) PrintError(argv[0]);
		   else { MinimumLLR = (double) x; }
		break;
             case 'm': 
		  if(sscanf(argv[arg],"-minnats=%lf",&d)==1){
			// fprintf(stderr,"%s: d=%.3f\n",argv[arg],d);
			if(d <= 0.0 || d > 100) print_error(MCS_USAGE_START);
			else MinimumNatsPerSeq=d;
		  } else if(sscanf(argv[arg],"-minsize=%d",&x)==1){
			if(x < 0) PrintError(argv[0]); else stringency=x; 
		  } else if(sscanf(argv[arg],"-maxdepth=%d",&x) == 1){
			if(x < 1) PrintError(argv[0]); else MaxDepthHierarchy=x+1; 
		  } else if(sscanf(argv[arg],"-maxcol=%d",&DefaultMaxCol) != 1){
                        PrintError(argv[0]);
                  } break;
             case 'N': 
		if(sscanf(argv[arg],"-Nth=%d",&x)==1){
		    if(x < 1) PrintError(argv[0]); NthSeqForDisplay=x;
		} else {
		    MaxNumNodes = IntOption(argv[arg],'N',1,100000,"cdhBPPS -N option input error");
		}
		break;
             case 'n': 
		if(sscanf(argv[arg],"-names=%s",str)==1){
		   NameSeedAln=AllocString(str);
		} break;
             case 'r': 
		if(sscanf(argv[arg],"-run=%c%d:%d=%s",&PrintMode,&KeyStart,&KeyEnd,str)==4){
		   KeySetName=AllocString(str);
		   if(KeyStart < 1 || KeyStart >= KeyEnd) PrintError(argv[0]); 
		} else if(sscanf(argv[arg],"-run=%c",&PrintMode)==1){ 
			if(!isalpha(PrintMode)) PrintError(argv[0]);
		} else if(sscanf(argv[arg],"-random=%d",&random_start)==1){
		    if(random_start <= 2){ PrintError(argv[0]); }
		} else if(sscanf(argv[arg],"-rtfN=%d",&x)==1){
		    if(x < 3) PrintError(argv[0]); NumHighlighted=x;
		} else if(sscanf(argv[arg],"-rtfF=%d",&x)==1){
		    if(x < 4 || x > 24) PrintError(argv[0]); font_size=x;
		} else if(sscanf(argv[arg],"-rtfP=%c",&c)==1){
		    if(strchr(" lLpP",c) == NULL) PrintError(argv[0]); page_format=c;
		} else if(strcmp("-rtf",argv[arg]) != 0) PrintError(argv[0]); else OutPutRTF=TRUE;
		break;
             case 'S': 
		{
		  if(strcmp("-Seed",argv[arg]) != 0) PrintError(argv[0]);
		  else {
		    FILE *sfp=open_file(argv[1],".seed","r");
		    if(fscanf(sfp,"-seed=%d",&seed) != 1) PrintError(argv[0]);
		    fclose(sfp); usingSeedOption=TRUE;
		  }
		} break;
             case 's': 
		if(strcmp("-strict",argv[arg]) == 0) StrictIndepend=TRUE;
#if 0
		else if(strcmp(argv[arg],"-sets")==0) SaveSets=TRUE; // else PrintError(argv[0]);
#endif
		else if(sscanf(argv[arg],"-seed=%d",&seed)!=1) PrintError(argv[0]);
		else seed_arg=arg;
		break;
             case 't': 
		if(sscanf(argv[arg],"-test=%c",&test_mode)==1){
			PrintMode='T';
		} else PrintError(argv[0]); break;
             case 'u': 
		if(strcmp("-usr",argv[arg]) == 0){ use_usr_sma=TRUE; }
		else PrintError(argv[0]); break;
             case 'v': 
		  if(strcmp("-verbose",argv[arg]) != 0) PrintError(argv[0]);
		  else this->VerboseOn(); 
		break;
             case 'X': 
		if(sscanf(argv[arg],"-XC=%d",&XC_line)==1){
		    if(XC_line < 1){ PrintError(argv[0]); }
		} else PrintError(argv[0]); 
		break;
             case 'x': 
		if(sscanf(argv[arg],"-xc=%d:%c:%d",&XC_size,&XC_mode,&x)==3){
		   if(XC_size < 1) PrintError(argv[0]); c=XC_mode;
		   if(c == 'X'){ XCuseX=TRUE; XCuseL=FALSE; }
		   else if(c == 'L'){ XCuseX=FALSE; XCuseL=TRUE; }
		   else if(c == 'B'){ XCuseX=XCuseL=TRUE; }
		   else if(c == 'N'){ XCuseX=XCuseL=FALSE; }
		   else PrintError(argv[0]);
		   if(x < 1 || x > 200) PrintError(argv[0]);
		   else XC_wt = (double)x/100.0;
		}
		break;
             default : PrintError(argv[0]);
           }
	  }else PrintError(argv[0]);
	}

	FILE *fp=0; chk_prefix=0;
	BooLean	ChkFileExists=FALSE;
	infile=AllocString(argv[1]);
	sprintf(str,"%s.chk",infile); fp=fopen(str,"r");
	if(fp != NULL){ usingSeedOption=TRUE; fclose(fp); ChkFileExists=TRUE; chk_prefix=AllocString(infile); }
	if(random_start > MaxNumNodes){
	  print_error("FATAL: random_start > maximum number of nodes");
	}
	if(use_usr_sma){
	    if(PrintMode != 'R' || usingSeedOption==FALSE){
		   print_error("ERROR 1: -usr option must be used with -run=R and *.chk files"); 
	    }
	} 
	//*********************** cma_typ *************************
   BooLean using_hierview=FALSE;
   if(mmafp){ 
	TrueMainCMA=ReadCMSA(mmafp,AB); fclose(mmafp); 
	assert(TrueMainCMA);
	if(TrueMainCMA==0) print_error("FATAL: Read input mmafp failed");
        swt = new swt_typ(TrueMainCMA,FALSE); hsw=swt->RtnHSW( ); own_hsw=FALSE;
   	using_hierview=TRUE;
   } else {
	sprintf(str,"%s.mma",infile); fp=fopen(str,"r");
	if(fp==NULL){
	   char *sp=strstr(infile,"_chk");
	   // if(sp == 0) sp=strstr(infile,"_ichk"); // needs more work...
	   if(sp && sp[4]==0){	// is this a check point file?
		FILE *sfp=open_file(infile,".seed","r");
		if(fscanf(sfp,"-seed=%d",&seed) != 1) print_error("seed file read error");
		fclose(sfp); usingSeedOption=TRUE;
	      chk_prefix=AllocString(infile); sp=strstr(chk_prefix,"_chk"); sp[0]=0;
	      fp=open_file(chk_prefix,".mma","r");
	      if(fp==NULL) print_error("FATAL: Opening of checkpoint mma file failed");
	   } else {
	        fprintf(stderr,"Could not open file \"%s\"\n",str);
		print_error("FATAL: Opening of input file failed");
	   }
	} TrueMainCMA=ReadCMSA(fp,AB); fclose(fp); fp=0; 
	if(TrueMainCMA==0) print_error("FATAL: Read input file failed");
	if(use_usr_sma && chk_prefix == 0){
	   if(!ChkFileExists) print_error("ERROR: -usr option must be used with -Seed -run=R and _chk files"); 
	}
	if(chk_prefix) sprintf(str,"%s.hsw",chk_prefix); else sprintf(str,"%s.hsw",infile);
	if((fp=fopen(str,"r")) == NULL){        // create file...
                swt = new swt_typ(TrueMainCMA,FALSE);
                hsw=swt->RtnHSW( );
		own_hsw=FALSE;
#if 1	// keep this for my own usage...
                fp = open_file(argv[1],".hsw","w");
                FWriteHSW(fp,hsw); fclose(fp);
		delete swt; swt=0;
        	fp=open_file(argv[1],".hsw","r");
                hsw=FReadHSW(fp,AB,TrueMainCMA); fclose(fp); own_hsw=TRUE;
#endif
        } else { hsw=FReadHSW(fp,AB,TrueMainCMA); fclose(fp); own_hsw=TRUE; }
   }

	sRandom(7061950);	// Always use the same seed for the random sequence set.
	NumRandom=1+(NumSeqsCMSA(TrueMainCMA)/3); 
#if 0	// This accomodates lineage-specific input files...
	sprintf(str,"%s.sets",argv[1]);
	if((fp=fopen(str,"r"))){
	   fclose(fp);
	   if((fp=fopen(str,"r"))){ 
		set_typ *set=ReadSets(fp,n); fclose(fp);
		assert(n > 0); i=SetN(set[1]);
		Int4 nrand=i-1-NumSeqsCMSA(TrueMainCMA);
		if(nrand > NumRandom) NumRandom=nrand;
		else if(nrand < NumRandom) print_error("omc_typ::Init(): input set size is too small");
		for(i=1; i <= n; i++) NilSet(set[i]); free(set);
	   }
	}
#else	// new routine...for hierview program.
	if(setfp){	// rewindable tmpfile() passed in.
		set_typ *set=ReadSets(setfp,n); 
		assert(n > 0); i=SetN(set[1]);
		Int4 nrand=i-1-NumSeqsCMSA(TrueMainCMA);
		if(nrand > NumRandom) NumRandom=nrand;
		else if(nrand < NumRandom) print_error("omc_typ::Init(): input set size is too small");
		for(i=1; i <= n; i++) NilSet(set[i]); free(set);
		rewind(setfp);
	} else { 
	  sprintf(str,"%s.sets",argv[1]); fp=fopen(str,"r");
	  if(fp){
		set_typ *set=ReadSets(fp,n); fclose(fp); 
		assert(n > 0); i=SetN(set[1]);
		Int4 nrand=i-1-NumSeqsCMSA(TrueMainCMA);
		if(nrand > NumRandom) NumRandom=nrand;
		else if(nrand < NumRandom) print_error("omc_typ::Init(): input set size is too small");
		for(i=1; i <= n; i++) NilSet(set[i]); free(set);
	  }
	}
#endif
	cma_typ rcma;	// Make the random sequence set...
	MainCMA=MkMainFileCMSA(TrueMainCMA,NumRandom,rcma);
	MainHSW=AddRandomHSW(hsw,TrueMainCMA,rcma,MainCMA);
	TotalNilCMSA(rcma); // destroy the temporary Random sequence alignment.

	if(using_hierview){ fp=0; }
	else {
	  fp=open_file(argv[1],".cmd","w");
	  for(arg=0; arg < argc; arg++) fprintf(fp,"%s ",argv[arg]);
	}
	if(seed == 18364592){  seed = (UInt4) time(NULL)/2; if(fp) fprintf(fp,"-seed=%d\n",seed); }
	else if(fp) fprintf(fp,"\n"); if(fp) fclose(fp); sRandom(seed);
	RandomSeed=seed;  // save for checkpoint

	SizeMainCMA = NumSeqsCMSA(MainCMA);
	SizeTrueMainCMA = NumSeqsCMSA(TrueMainCMA);
#if 1	// save swt from above...else swt==0;
	if(swt==0){ swt = new swt_typ(hsw);  own_hsw=TRUE; }
	lpr = new lpr_typ(TrueMainCMA,swt,del_as_random);
#else
	swt = new swt_typ(MainHSW); 
	lpr = new lpr_typ(MainCMA,swt);
#endif
	MaxNumNodesPlus=(MaxNumNodes*2);

	SetStringency(stringency);
	if(hptfp){ fp=hptfp; assert(setfp); FromSetsInit( ); }   // For hierview program call.
	else {		// from Checkpoint file...new format...
	  sprintf(str,"%s.chk",argv[1]);
	  if((fp=fopen(str,"r")) != NULL){
	      Int4    M,NumSMA=0,sd=0;
	      bpcp_typ bpcp; bpcp.Read(fp,AB); fclose(fp);
	      set_typ *chk_sets=bpcp.RtnSets(M); RandomSeed=bpcp.RtnSeed();
	      cma_typ *chk_sma=bpcp.RtnSMA(NumSMA); Hpt=bpcp.RtnHpt();
	      if(M != Hpt->NumSets()) print_error("FATAL: *.hpt and *.sets files are inconsistent"); 
#if 0
	      bpcp.Put(argv[1]);	// create old format files.
#endif
	      if(FocusNode > 0){ stable=MakeSet(MaxNumNodesPlus); ClearSet(stable); }
	      FromSetsInit(chk_sma,chk_sets); 
	  } else {	// use old checkpoint format...
	    sprintf(str,"%s.hpt",argv[1]);
	    if((fp=fopen(str,"r")) == NULL){
		if(XC_line > 0) print_error("ERROR: -XC_line option requires an input hpt file"); 
		if(use_usr_sma){
		   print_error("ERROR: -usr option must be used with -Seed -run=R and _chk files"); 
		}
		if(FocusNode > 1) PrintError(argv[0]); else AbInitioInit( ); 
	    } else { 	// opened hpt file okay.
		fclose(fp); 
		sprintf(str,"%s.sets",argv[1]);
		if((fp=fopen(str,"r")) == NULL){ 
		   if(use_usr_sma){
		     print_error("ERROR: -usr option must be used with -Seed -run=R and _chk files"); 
		   }
		   if(FocusNode > 0) PrintError(argv[0]); else FromFileInit( ); 
		} else {	// opened sets file okay...
		   if(FocusNode > 0){ stable=MakeSet(MaxNumNodesPlus); ClearSet(stable); }
		   fclose(fp); FromSetsInit( ); 
		}
	    }
	  }
	}
	mcs->SetStringency(MinimumLLR,MinimumSetSize,stable);
	Int4 SetSize=mcs->GetSetSize();
	RandomSet=MakeSet(SetSize); ClearSet(RandomSet);
	for(j=1,i=NumSeqsCMSA(TrueMainCMA)+1; j <= NumRandom; i++,j++) AddSet(i,RandomSet);
	SetFG=MakeSet(SetSize); SetBG=MakeSet(SetSize); TmpSet=MakeSet(SetSize);
	TmpFG=MakeSet(SetSize); TmpBG=MakeSet(SetSize); 

	ChildSetBG=MakeSet(MaxNumNodesPlus+1); ChildSetFG=MakeSet(MaxNumNodesPlus+1);
	ChildSetU=MakeSet(MaxNumNodesPlus+1);   // allow for added nodes.
	ChildSetBoth=MakeSet(MaxNumNodesPlus+1);   // allow for added nodes.
	free(chk_prefix); chk_prefix=0;
}

void	omc_typ::FromSetsInit(cma_typ *chk_sma, set_typ *chk_sets)
// hpt and sma files, but no sets file.
{
	InitMode='s';
	BooLean	rename=TRUE;
	if(PrintMode) rename=FALSE;
	if(FocusNode) rename=FALSE;	// WARNING: assumes input from previous run using Set%d notation!!!
	cma_typ *iSMA=RenameHptSMA(rename,chk_sma); MkFileIdHeap(); 
	Int4 id=0;
	SetUpFocusedSrch();
	if(PrintMode) ToggleOutFileName(TRUE); else ToggleOutFileName(); 
	this->SetDefaultArguments( );
        // Argv[Argc]=AllocString("-RandomInit"); Argc++;
	cma_typ *in_sma=OrderNodeSMAs(iSMA,rename); free(iSMA);
	Int4	n,i;
	set_typ *set=0;
	if(chk_sets == 0){
	   FILE *fp = 0;
	   if(setfp) fp=setfp; else fp=open_file(infile,".sets","r");
	   set=ReadSets(fp,n); fclose(fp); setfp=0;
	   if(n != Hpt->NumSets()) print_error("FATAL: *.hpt and *.sets files are inconsistent"); 
	} else { set=chk_sets; n= Hpt->NumSets(); }
	mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,n,set,Hpt,in_sma,Argc,Argv);
	mcs->UnLabelAllSeqs( );
	mcs->FileID=id; free(in_sma); mcs->SampledColumn=0; mcs->IsTreeHpt=TRUE;
	free(set);	// for(i=1; i <= n; i++) NilSet(set[i]); freed by mcs_typ...
	if(!Evolve) mcs->DoNotEvolve();	else mcs->DoEvolve(); // start out not evolving.
	// mcs->PutHyperPartition(stderr); // Hpt->PutHyperPartition(stderr);
}


cma_typ	*omc_typ::RenameHptSMA(BooLean rename, cma_typ *chk_sma)
// Order the input_sets to be consistent with Hpt.
// when using an input Hpt store the names and use Set%d instead.
// Also Store the sma correspondence for output later.
{
      Int4      a,s,i,id,N,NumInSMA;
      FILE *fp=0; 
      cma_typ *InSMA=0; 
      if(smafp){ assert(chk_sma == 0); assert(hptfp); fp=smafp; smafp=0; }
      else if(use_usr_sma){
	assert(chk_prefix || chk_sma); 
	sprintf(str,"%s_usr.sma",chk_prefix);
        if((fp=fopen(str,"r")) == NULL){	// retain _bst.sma files search...
	   sprintf(str,"%s_bst.sma",chk_prefix);
           if((fp=fopen(str,"r")) == NULL){
	      fprintf(stderr,"FATAL: %s_usr.sma file not found\n",chk_prefix); exit(1);
	      print_error("FATAL: usr.sma file not found");
	   }
	}
	if(chk_sma){ for(i=1; chk_sma[i]; i++) TotalNilCMSA(chk_sma[i]); free(chk_sma); }
      } else if(chk_sma){ 	// Use *.chk file as input
	 InSMA=chk_sma; for(N=0,i=1; chk_sma[i]; i++) N++; 
      } else fp = open_file(infile,".sma","r");
      if(fp){ InSMA=MultiReadCMSA(fp,&N,0,AB); fclose(fp); }
      cma_typ	*SMA; NEW(SMA,MaxNumNodesPlus+3,cma_typ);
      if(Hpt==0){
        if(hptfp) { Hpt = new hpt_typ(hptfp); fclose(hptfp); hptfp=0;}
        else { Hpt = new hpt_typ(infile); } 
      }
      // Hpt->Put(stderr);
      if(Hpt == 0 || Hpt->NumSets() > MaxNumNodes){
          fprintf(stderr,"Hpt->NumSets() = %d > %d = MaxNumNodesPlus\n",
                                        Hpt->NumSets(),MaxNumNodesPlus);
          assert(Hpt != 0 && Hpt->NumSets() <= MaxNumNodes); // 2 = add internal + leaf at same time.
      } NumInSMA=N;
      // InHpt=Hpt->Copy();  // make a copy...
      for(s=1; s < Hpt->NumSets(); s++){
	// InNameHpt[s]=AllocString(Hpt->ElmntSetName(s));
	BooLean found=FALSE;
        for(a=1; a <= N; a++){
	   if(InSMA[a] == 0) continue; 	// for rename == FALSE..
	   if(strcmp(Hpt->ElmntSetName(s),NameCMSA(InSMA[a])) == 0){
      		// then rename the SMA and Hpt file names to Set%d format...
		char *tmp_str=0;
		if(rename){
			tmp_str=AllocString(NameCMSA(InSMA[a])); sprintf(str,"Set%d",s); 
			ReNameCMSA(str,InSMA[a]); 
			SMA[s]=MakeConsensusCMSA(InSMA[a]); 
		} else { SMA[s]=InSMA[a]; InSMA[a]=0; }
		if(rename){ Hpt->ReNameSet(s,str); ReNameCMSA(tmp_str,InSMA[a]); free(tmp_str); }
		found=TRUE; break;
	   }
	} if(!found){
	    fprintf(stderr,"Seed sequence for file %s not found\n",Hpt->ElmntSetName(s));
	    print_error("hpt and sma files inconsistent");
	}
      } // Hpt->Put(stderr);
      for(i=1; i <= N; i++){ if(InSMA[i]) TotalNilCMSA(InSMA[i]); } free(InSMA);
      return SMA;
}

cma_typ	*omc_typ::OrderNodeSMAs(cma_typ *iSMA,BooLean renamed)
// Order the input_sets to be consistent with Hpt.
// NOTE:  Not really necessary as mcs_typ does this!!!!
{
      Int4      s,id;
      cma_typ	*out_sma=0;
      NEW(out_sma,MaxNumNodesPlus + 3,cma_typ);
      // assert(Hpt != 0 && Hpt->NumSets() <= MaxNumNodes + 1);
      if(Hpt == 0 || Hpt->NumSets() > MaxNumNodesPlus){
          fprintf(stderr,"Hpt->NumSets() = %d > %d = MaxNumNodesPlus\n",
                                        Hpt->NumSets(),MaxNumNodesPlus);
          assert(Hpt != 0 && Hpt->NumSets() <= MaxNumNodes); // 2 = add internal + leaf at same time.
      }
      for(s=1; s < Hpt->NumSets(); s++){
        char *name=Hpt->ElmntSetName(s);
	if(renamed){
         if(sscanf(name,"Set%d",&id) == 1){      // if name sma file  == name hpt set...
           assert(id > 0);
           if(id > (MaxNumNodes + 1)){  // +1 covers an additional, internal node added at the end...
                Hpt->PutHyperPartition(stderr); assert(id <= MaxNumNodes + 1);
           }
           assert(iSMA[id] != 0); out_sma[s]=iSMA[id]; 
         } else {
           fprintf(stderr,"error on: set = %d; id = %d; name = %s\n",s,id,name);
           print_error("omc_typ::OrderNodeSMAs() error 1");
         }
	} else {
		char *name2=NameCMSA(iSMA[s]);
		assert(strcmp(name,name2)==0);
		assert(iSMA[s] != 0); out_sma[s]=iSMA[s];
	}
      } return out_sma;
}

void    omc_typ::MkFileIdHeap( )
{ FileIdHeap=dheap(28,4); for(Int4 i=1; i <= 26; i++) insrtHeap(i,i,FileIdHeap); }

Int4    omc_typ::SetDefaultArguments(Int4 ppb)
{
	if(Argv[0]==0) Argv[0]=AllocString("omcBPPS");	// name of program.
        for(Int4 i=1; i < Argc; i++){ free(Argv[i]); Argv[i]=0; }
        Argv[1]=AllocString(outfilename); Argc=2;
        // Argv[Argc] = AllocString("-rho=0.003"); Argc++;
        Argv[Argc] = AllocString("-nocsq"); Argc++;    // Use the consensus sequences provided!
        sprintf(str,"-ppb=%d",ppb);        // default: "-ppb=1000000"
        Argv[Argc]=AllocString(str); Argc++;
        Argv[Argc]=AllocString("-iter=1:2:2"); Argc++;
        Argv[Argc] = AllocString("-col_start=2"); Argc++;
        Argv[Argc] = AllocString("-tree"); Argc++;
        if(!del_as_random){ Argv[Argc] = AllocString("-del"); Argc++; }
	if(StrictIndepend){ Argv[Argc] = AllocString("-strict"); Argc++; }
	if(MinimumNatsPerSeq > 0){ 
	    sprintf(str,"-minnats=%.2f",MinimumNatsPerSeq); Argv[Argc] = AllocString(str); Argc++; 
	}
	if(DefaultMaxCol >= 5){
		sprintf(str,"-maxcol=%d",DefaultMaxCol); Argv[Argc] = AllocString(str); Argc++; 
// fprintf(stderr,"%s\n",str);
	}
// Argv[Argc] = AllocString("-col=10:20"); Argc++;
// Argv[Argc] = AllocString("-A=10:20"); Argc++;
// Argv[Argc] = AllocString("-Am=10:40"); Argc++;
        return Argc;
}

#if 0
void	AbInitioReInit()
{
	assert(InitMode=='z');
}
#endif

void    omc_typ::AbInitioInit( )
{
	InitMode='z';
	if(flatten_hpt) print_error("FATAL: -flatten option requires an input hpt file");
#if 0	// 2 nodes...
	const char hpt_str[]="HyperParTition:\n!!\n+- 1.Set1?\n++ 2.Set2!\n-o 3.Random=20000.\n\n";
	FILE *fp=tmpfile(); fprintf(fp,"%s",hpt_str); rewind(fp); Hpt = new hpt_typ(fp); fclose(fp); 
	Hpt->Put(stderr);
	cma_typ *in_sma; NEW(in_sma,10,cma_typ);
	in_sma[1]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA("Set1",in_sma[1]);
	in_sma[2]=GetBestCsqCMSA(TrueMainCMA); RenameCMSA("Set2",in_sma[2]);
#else	// Root node only.
	const char hpt_str[]="HyperParTition:\n!\n+ 1.Set1?\n- 2.Random=20000.\n\n";
	FILE *fp=tmpfile(); fprintf(fp,"%s",hpt_str); rewind(fp); Hpt = new hpt_typ(fp); fclose(fp); 
	Hpt->Put(stderr);
	cma_typ *in_sma; NEW(in_sma,10,cma_typ);
	in_sma[1]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA("Set1",in_sma[1]);
#endif
	MkFileIdHeap(); 
	Int4 id=ToggleOutFileName(); this->SetDefaultArguments( );
	mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,Hpt,in_sma,Argc,Argv);
	mcs->UnLabelAllSeqs( );
	mcs->FileID=id; free(in_sma); mcs->SampledColumn=0; mcs->IsTreeHpt=TRUE;
	if(!Evolve) mcs->DoNotEvolve();	else mcs->DoEvolve(); // start as not evolving by default.
}

BooLean	omc_typ::SetUpFocusedSrch()
{
	Int4 i,id;
	if(FocusNode == 0) return FALSE;
	else {
	    assert(Hpt);
	    if(FocusNode < 1 || FocusNode >= Hpt->NumSets()) PrintError(program_name);
	    if(FocusNode == 1  && FocusSeq == 0) PrintError(program_name);
	    set_typ subtree=Hpt->MkSubTreeSet(FocusNode);
      	    for(i=1; i < Hpt->NumSets(); i++){
		if(MemberSet(i,subtree)){
		   if(FocusSeq==0) continue;	// in this case sample from entire subtree.
		   else if(FocusNode == i) continue;
		} id=Hpt->ItoSetID(i); AddSet(id,stable);
	    } 
	    // Hpt->Put(stderr); PutSet(stderr,subtree); PutSet(stderr,stable);
	    NilSet(subtree);
	} // exit(1);
	return TRUE;
}

void	omc_typ::FromFileInit( )
// hpt and sma files, but no sets file.
{
	InitMode='f';
	BooLean	rename=TRUE;
	if(PrintMode=='O' || PrintMode=='o' || PrintMode=='C') rename=FALSE;
	cma_typ *iSMA=RenameHptSMA(rename); MkFileIdHeap(); 
	Int4 id=0;
	// SetUpFocusedSrch();  // currently I require *.sets as input...
	if(PrintMode) ToggleOutFileName(TRUE); else ToggleOutFileName(); 
	this->SetDefaultArguments( );
        // Argv[Argc]=AllocString("-RandomInit"); Argc++;
	cma_typ *in_sma=OrderNodeSMAs(iSMA,rename); free(iSMA);
	mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,Hpt,in_sma,Argc,Argv);
	mcs->UnLabelAllSeqs( );
	mcs->FileID=id; free(in_sma); mcs->SampledColumn=0; mcs->IsTreeHpt=TRUE;
	if(!Evolve) mcs->DoNotEvolve();	else mcs->DoEvolve(); // start out not evolving.
}

void	omc_typ::PutStringency( )
{
	double	d,D,aveSS=0.0,aveLLR=0.0;
	for(Int4 x = 25; x <= 1000; x+=25){
		this->SetStringency(x);
		d=(double) MinimumSplitSize/(double) MinimumSetSize;
		D=(double) MinimumLLR/(double) MinimumSetSize;
		fprintf(stderr,"%d: %d %d %.1f %.3f %3f\n",
			x, MinimumSetSize,MinimumSplitSize, MinimumLLR,d,D);
		aveSS += d; aveLLR += D; 
	}
	aveSS = aveSS/16;
	aveLLR = aveLLR/16;
	fprintf(stderr,"\naveSS = %.3f; aveLLR = %.3f\n\n",aveSS,aveLLR);
	// aveSS = 2.286; aveLLR = 2.072 --> 2.3 and 2.0 ?
	exit(1);
}

void	omc_typ::SetStringency(Int4 x)
{
	MinimumSetSize=x; 
	MinimumSplitSize=(Int4) ceil(2.3*(double)MinimumSetSize);
	MinimumLLR=1.5*(double) MinimumSetSize;	  // 1.0 might be better.
	MinimumPreLLR = 0.25*MinimumLLR; 
	// MinimumPreLLR = 0.95*MinimumLLR; 
}

void	omc_typ::Free()
{
	Int4	i;
	for(i=0; i < Argc; i++) free(Argv[i]); free(Argv);
	if(mcs) DeleteMCS(mcs);  // WARNING: deletes old Hpt & closes output fp.
	if(best_mcs) DeleteMCS(best_mcs);  // WARNING: deletes old Hpt & closes output fp.
	if(FileIdHeap) Nildheap(FileIdHeap);
	if(swt) delete swt;
	if(TrueMainCMA) TotalNilCMSA(TrueMainCMA);
	delete lpr;
	if(MainCMA) TotalNilCMSA(MainCMA);
	if(hsw && own_hsw) NilHSW(hsw);
	if(MainHSW) NilHSW(MainHSW);
	if(DisplaySet){
	    for(i=1; i <= MaxNumNodesPlus; i++) if(DisplaySet[i]) NilSet(DisplaySet[i]);
	    free(DisplaySet);
	} NilSet(RandomSet);
	if(infile) free(infile);
	if(KeySetName) free(KeySetName);
	if(NameSeedAln) free(NameSeedAln);
	if(outfilename) free(outfilename);
	if(SetFG) NilSet(SetFG);
	if(SetBG) NilSet(SetBG);
	if(TmpFG) NilSet(TmpFG);
	if(TmpBG) NilSet(TmpBG);
	if(TmpSet) NilSet(TmpSet);
	NilSet(ChildSetBG); NilSet(ChildSetFG); NilSet(ChildSetU); NilSet(ChildSetBoth);
	if(goodHG) NilHist(goodHG);
	if(badHG) NilHist(badHG);
	if(testHG) NilHist(testHG);
	if(triesHG) NilHist(triesHG);
	if(stable) NilSet(stable);
	if(aafreq) free(aafreq);
	delete MP;
	NilAlpha(AB);
	if(logfp) fclose(logfp);
	if(ForcedNode != 0) free(ForcedNode);
	fprintf(stderr,"\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

