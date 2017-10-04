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

#include "chn_typ.h"

const char USAGE_START[]="USAGE: chn_see chn_fafile [options]\n\
   cmafafile = fafile with aligned sequences\n\
   options:\n\
     -a<int>               - allowed mismatch pattern positions (used with -P option)\n\
     -alpha=<int>:<int>=<real>   - alpha parameters for BPPS-like constraints.\n\
                              <int1> = match prior pseudocounts (range <int1>: 1-10000)\n\
                              <int2> = mismatch prior pseudocounts (range <int2>: 1-10000)\n\
                              <real> = posterior alpha from BPPS (range <real>: 0.0-1.0)\n\
     -B<int>:<int>         - set begin & end range (1st seq.) to be printed in alignment\n\
     -B=<int>:<int>        - set begin & end range (column position) to be printed in alignment\n\
     -b<int>               - minimum bars in histogram to show sidechain (default: 5)\n\
     -C<real>              - remove sequences in database alignment with percent\n\
                              identity >= <real> to the concensus orthologous sequence\n\
                              (range: 25..100)\n\
     -C=<str>              - Colors used for alignment categories (default: CyPBGRcmbgr)\n\
                              'C' = cyan          'c' = teal\n\
                              'y' = dk yellow\n\
                              'P' = pink\n\
                              'B' = blue          'b' = dk blue\n\
                              'G' = green         'g' = dk green\n\
                              'R' = red           'r' = dk red\n\
                              'M' = magenta       'm' = violet\n\
                              Note: need to input a string of length >= 6\n\
     -c<real>              - minimum fraction of conserved residues needed \n\
                              in order to show constraints in foreground set (default: 0.0)\n\
     -concise              - Don't create extra files\n\
     -D                    - Don't weight sequences\n\
     -Exp=<real>           - log10(Expected # of residue patterns) displayed (default: 0.0)\n\
     -F=<int>              - font size (4-24 points | default: 7 points)\n\
     -f<real>              - Set the maximum fraction of residues to highlight\n\
                              This option resets the linear-to-log parameter\n\
                              (default: 0.20 or 20%).\n\
     -g<int>,<int>         - gap opening and extension penalities for -M= and -m= options\n\
     -global               - Use global sequence weighting\n\
     -H                    - Don't show (i.e., hide) indel information\n\
     -hide=<int>           - hide insert regions >= <int> in length (<int> must be > 2)\n\
     -L=<int>              - maximum number of lines to display below alignment (range: 2-20)\n\
     -M                    - show marginal probability for questionable regions\n\
     -M<real>              - use marginal probability cutoff of <real> for display\n\
     -M=<cmafile>          - use background frequencies as specified by motif cmafile\n\
     -method=<char>        - set statistical set method for BPPS LPR(default: 'M')\n\
                        	'B' = Binomial model \n\
                        	'P' = Binomial BPPS model (cma input files must specify BPPS parameters)\n\
                        	'U' = used by cmc_typ only. \n\
     -m=<cmafile>          - use frequencies specified by motif cmafile, but\n\
                              don't use these as the background for the Main Set\n\
     -m                    - compute marginal probability weighted freq\n\
     -Nth=<int>            - Use <int>th sequence for numbering display alignment\n\
     -N=<int>              - Number of constrained positions to highlight\n\
                              This option resets the linear-to-log parameter\n\
                              (default: the number of pattern positions).\n\
     -O                    - Output a *.vsi file\n\
     -o                    - use pattern only mode\n\
     -P<chars><int>         - Define an alignment subset based on the occurance of\n\
                                 residue <char> at position <int> in the query sequence\n\
                                 (if lower case char are given the inverted set is used\n\
                                  e.g., -Pilv230 --> not I,L or V will be selected)\n\
     -P=<chars><int>        - same as above except <int> corresponds to column not query position\n\
     -p<real>              - cumulative binomial probability cutoff (for mimimum shading)\n\
                              (default: 0.01; range: 0.0 to 0.5)\n\
     -Q                    - Quiet (concise, non-verbose) rich text output\n\
     -R                    - Remove family alignment\n\
     -S                    - show standard alignment only\n\
     -S=<char>             - page setup option\n\
                              'P' = 11 x 17 portrait\n\
                              'p' = 8.5 x 11 portrait\n\
                              'L' = 11 x 17 landscape\n\
                              'l' = 8.5 x 11 landscape\n\
     -s<real>              - use fraction seq. aligned cutoff of <real> for display\n\
     -sets_root=<char>     - Residue sets option for root (default 'L')\n\
     -sets=<char>          - Residue sets options (default 'L')\n\
                        	'D' = believed to be the best sets\n\
                        	'O' = allow only one residue at pattern positions\n\
                        	'T' = tiny number of very similar residues\n\
                        	'S' = small number of similar residues\n\
                        	'M' = moderate number of residue sets\n\
                        	'L' = large number of residue sets\n\
                        	'R' = largest number of residue sets for root nodes\n\
                        	' ' = use alternative, blosum62-based mode\n\
     -see=<char>           - See residues in sets for each set option <char>\n\
     -U                    - use previous alignment as background instead of main set\n\
     -u                    - use ortholog freqs for pseudocounts\n\
     -V                    - Verbose rich text output\n\
     -2                    - show predicted secondary structure (use with -hide=500 or more)\n\
     -W                    - Write out the second and last *.cma alignment files\n\
     NOTE: a '^' may be placed in status array to ignore a selected position\n\
\n";

void	chn_typ::GetArg(Int4 argc,char *argv[])
{
	Int4	i,j,s,m,n,arg,r;
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
	     case 'a': 
		if(sscanf(argv[arg],"-alpha=%d:%d=%lf",&A0,&B0,&Alpha) == 3){
			if(Alpha < 0.0 || Alpha > 1.0) print_error(USAGE_START);
			if(A0 < 1 || B0 < 1) print_error(USAGE_START);
			if(A0 > 100000 || B0 > 100000) print_error(USAGE_START);
		} else MaxMisMatches=IntOption(argv[arg],'a',0,100,USAGE_START); 
		break;
             case 'B': 
		{
		  if(sscanf(argv[arg],"-B=%d:%d",&Begin,&End) == 2){
		     UseColPos=TRUE;
		     if(Begin < 1 || Begin > End) print_error(USAGE_START);
		  } else if(sscanf(argv[arg],"-B%d:%d",&Begin,&End) != 2)
                        print_error(USAGE_START); 
		  if(Begin < 1 || Begin > End) print_error(USAGE_START);
		  // Begin--; // start at zero not one. Done below...
		} break;
	     case 'b': minbars=IntOption(argv[arg],'b',1,24,USAGE_START); break;
	     case 'C': 
		if(argv[arg][2]=='='){
		  ColorString=AllocString(argv[arg]+2); ColorString[0]=' ';
		  if(strlen(ColorString) < 7) print_error(USAGE_START);
		} else {
		  percent_cut=RealOption(argv[arg],'C',0.0,100.0,USAGE_START); 
		}
		break;
	     case 'c':
		if(strcmp("-cha",argv[arg]) == 0){
			cha_as_input=TRUE;
		} else if(strcmp("-concise",argv[arg]) == 0){
			extra_files=FALSE;
		} else MinResFreqCutoff=RealOption(argv[arg],'c',0.0,1.0,USAGE_START); 
		break;
	     case 'D': NoWeights=TRUE; break;
	     case 'E': 
		if(sscanf(argv[arg],"-Exp=%lf",&ExpPatterns) != 1) print_error(USAGE_START); 
		// if(ExpPatterns > 0.0) print_error("-Exp value must be <= 0.0");
		break;
	     case 'F': fontsize=IntOption(argv[arg],'F',4,24,USAGE_START); 
			fontsize = fontsize*2;
		break;
	     case 'f': FractHighlight=RealOption(argv[arg],'f',0.00000,1.0,USAGE_START); 
			fprintf(stderr,"################ FractHighlight = %f ################\n",
					FractHighlight);
			InputFraction=TRUE;
		break;
             case 'g':
		if(strcmp("-global",argv[arg]) == 0){
			GlobalSqWts=TRUE;
			StartAlpha=0;
		} else if(sscanf(argv[arg],"-g%d,%d",&insert,&extend) != 2)
                                        print_error(USAGE_START); 
                     if(insert < 0 || extend < 0) print_error(USAGE_START);
                break;
	     case 'h': 	// hide inserts...
		if(sscanf(argv[arg],"-hide=%ld",&maxlen_gnull)==1) {
		 	if(maxlen_gnull < 3) print_error(USAGE_START); 
			// else fprintf(stderr,"maxlen_gnull = %d\n",maxlen_gnull);
		}
		// if(strcmp("-hide",argv[arg]) == 0) maxlen_gnull = 0;
		else print_error(USAGE_START);
                break;
	     case 'H': 
		if(argv[arg][2]==0){
		   HideIndelInfo=TRUE;
		} else print_error(USAGE_START);
		break;
	     case 'L': 
		// if(argv[arg][2] == '='){ print_error(USAGE_START); } else 
		MaxConcensusLines=IntOption(argv[arg],'L',2,20,USAGE_START);
	        break;
	     case 'M': 
		if(argv[arg][2] == 0){		// then show Marginal probabilities 
		  show_marg_prob=TRUE;
		} else if(argv[arg][2] == '='){ // then this is a BG_CMA file
		  if(argv[arg][3]==0) print_error(USAGE_START);
		  BG_CMA = ReadCMSA2(argv[arg]+3,AB);
		  ReNameCMSA("Motifs",BG_CMA);
		  use_bg_cma=TRUE;
		} else{
			mp_cutoff=RealOption(argv[arg],'M',0.0,1.0,USAGE_START); 
	        } break;
	     case 'm': 
		if(sscanf(argv[arg],"-method=%c",&ModeLPR)==1){
                        if(ModeLPR != 'B' && ModeLPR != 'P'){
				fprintf(stderr,"ModeLPR = %c\n",ModeLPR);
                                print_error(USAGE_START);
                        } 
		} else if(argv[arg][2] == '='){ // then this is a BG_CMA file
		  if(argv[arg][3]==0) print_error(USAGE_START);
		  BG_CMA = ReadCMSA2(argv[arg]+3,AB);
		  ReNameCMSA("motifs",BG_CMA);
		  use_bg_cma=FALSE;
		} else if(argv[arg][2] ==0){
			compute_mpwf=TRUE; 
		} else { print_error(USAGE_START); }
	        break;
	     case 'N': 
	      {
		Int4	n=0;
		if(sscanf(argv[arg],"-N=%d",&n)==1){
			if(n <= 0) print_error(USAGE_START);
			else InputFractionAsNumber=n;
		} else if(sscanf(argv[arg],"-Nth=%d",&Seq4numbering)==1){
			if(Seq4numbering < 1) print_error(USAGE_START);
	        } else print_error(USAGE_START);
		   //  } else Seq4numbering=IntOption(argv[arg],'Nth',1,100,USAGE_START); {
	      } break;
	     case 'O': 
		if(argv[arg][2] ==0) output_vsi=TRUE; 
		else print_error(USAGE_START); 
		break;
	     case 'o': 
		if(argv[arg][2] ==0) PatternOnlyMode=TRUE; 
		else print_error(USAGE_START); 
		break;
	     case 'P': 
		   {
#if 0
			if(Residue) print_error(USAGE_START);
			if(sscanf(argv[arg],"-P%[a-zA-Z]%d",residue_str,&Position) != 2){
                         	print_error(USAGE_START);
			} Residue = residue_str[0];
#elif 0		  // NEW multiple input pattern mode...
			NumResidues++;	// increment...
			if(NumResidues >= MAX_CHN_PATTERNS)
				print_error("Too many input patterns");
			if(sscanf(argv[arg],"-P%[a-zA-Z]%d",
				residue_str[NumResidues],
				&Position[NumResidues]) != 2){
                         		print_error(USAGE_START);
			} Residue[NumResidues] = residue_str[NumResidues][0];
#else		   // Even NEWER multiple input pattern mode... -PFVILM37,YFLM40,G54
		    char *Arg=0;
		    if(argv[arg][2]=='='){
		    	SelectColPos=TRUE;
		    	Arg=argv[arg]+3;
		    } else {
		    	Arg=argv[arg]+2;
		    }
		    do {
			NumResidues++;	// increment...
			if(NumResidues >= MAX_CHN_PATTERNS)
				print_error("Too many input patterns");
			if(sscanf(Arg,"%[a-zA-Z]%d",
				residue_str[NumResidues],
				&Position[NumResidues]) != 2){
                         		print_error(USAGE_START);
			} Residue[NumResidues] = residue_str[NumResidues][0];
			while(Arg[0] != ',' && Arg[0] != 0) Arg++;
			if(Arg[0] == ',') Arg++;
		    } while(Arg[0]);
#endif
		   } break;
	     case 'p': cbp_cut=RealOption(argv[arg],'p',0.0,0.5,USAGE_START); 
			break;
	     case 'Q': verbose=FALSE; Show2ndary=FALSE; break;
             case 'R': 
		if(argv[arg][2] == 0) RmFamilyAln=TRUE; 
		else print_error(USAGE_START);
		break;
	     case 'S': 
		if(argv[arg][2] == 0){ StdAlignmentOnly=TRUE; }
		else if(argv[arg][2] == '=' && argv[arg][4] == 0){
                        PageSetUp=argv[arg][3];
                        if(!isalpha(PageSetUp)) print_error(USAGE_START);
                } else print_error(USAGE_START);
		break;
             case 's':
                if(sscanf(argv[arg],"-see=%c",&sets_mode)==1){
			rst_typ *rst=new rst_typ(sets_mode);
			rst->Put(stdout); exit(0);
                } else if(sscanf(argv[arg],"-sets=%c",&sets_mode)==1){
                   if(sets_mode != 'M' && sets_mode != 'O' && sets_mode != 'L'
			&& sets_mode != 'G' && sets_mode != 'R'){
                                print_error(USAGE_START);
                   }
                } else if(sscanf(argv[arg],"-sets_root=%c",&root_sets_mode)==1){
		   if(strchr("TSMDOLR",root_sets_mode) == NULL) print_error(USAGE_START);
                } else fsa_cutoff=RealOption(argv[arg],'s',0.0,1.0,USAGE_START); 
		break;
	     case 'V': verbose=TRUE; Show2ndary=TRUE; break;
	     case 'u': use_ocma_pseudo=TRUE; break;
	     case 'U': use_sfbg=FALSE; break;
             case 'W': OutputAln=TRUE; break;
             // case 'w': OutputSeq=TRUE; break;
             case '2': Show2ndary=TRUE; break;
             case ' ': break; // ignore these...
	     default: ; print_error(USAGE_START);
	   }
	}
}

void	chn_typ::InitAsNull()
{
	Int4	i;
	for(i=0; i <= MAX_ALN_CHN_TYP; i++){ // these arrays go a bit beyond MAX_ALN_CHN_TYP.
		MargProb[i]=0; NullFreq[i]=0; FractSeqAln[0]=0; Obs[i]=0;
		rtf[i]=0; rtfQ[i]=0;
		SecondStrct[i]=0; 
		MCMA[i]=0; ColorCode[i]=0; AlignCode[i]=0; Freq[i]=0;
		NumSqMain[i]=0; BackGrnd[i]=0; Hist[i]=0; SuperAln[i]=0;
		WtFrq[i]=0; ObsIWt[i]=0; WtNumSq[i]=0; WtNumSeq[i]=0;
		merged_rasmol[i]=UCHAR_MAX;
	}
	for(i=0; i < MAX_CHN_PATTERNS; i++){
		InvertResSet[i]=FALSE;
		Residue[i]=0; Position[i]=Position0[i]=0;
	}
	A0s=0; B0s=0; Alphas=0; SetMode=0;
}

void	chn_typ::Init(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,double default_percent_cut,char *NameFReadSWT,hsw_typ *HSW)
{ 
	if(HSW != 0) assert(NameFReadSWT==0);
	// initialize parameters...
	// cbp_cut=0.00005;
	// GlobalSqWts=FALSE;
	Xconserved=0;
	RTFsCreated=FALSE;
	GlobalSqWts=TRUE;
	StartAlpha=1;	// use global only if specified on command line.
	sets_mode='L';		// use 'L' by default...
	root_sets_mode='L';		// use 'L' by default...
	cbp_cut=0.0;
	StdAlignmentOnly=FALSE;
	InputFraction=FALSE;
	Seq4numbering=0;
	status=0;	// need to free this? (Valgrind complaining...)
	// cbp_cut=0.001;
	mp_cutoff=0.0;fsa_cutoff=0.0;
	percent_cut=default_percent_cut;
	PageSetUp='P';
	char_per_line=0;fontsize=16; // font size is 8 points by default.
	verbose=FALSE;
	Alpha=-1.0; A0=0; B0=0;	// if set to zero then ModeLPR defaults to 'U';
	OwnCMAs=TRUE;

	// For IntegerWeights
	MainCMA=0;
	MainSqWt=0; SqIWtMain=0;

	extra_files=TRUE;
	ExpPatterns=0.0;
	NoWeights=FALSE;use_ocma_pseudo=FALSE; use_sfbg=TRUE;
	FractHighlight=0.2;	// set to 20% aligned.
	FractHighlightTOP=0.95;	// set to 95% aligned.
	OutputAln=FALSE;compute_mpwf=FALSE;
	minbars=5;
	Begin=0; End=INT4_MAX;
	UseColPos=FALSE;
	SelectColPos=FALSE;
	NumAnalysis=0; 
	Show2ndary=FALSE; 
	show_marg_prob=FALSE;
	RmFamilyAln=FALSE;
	PatternOnlyMode=FALSE;
	cha_as_input=FALSE;
	TAX_CMA=0;
	output_vsi=FALSE;
	BG_CMA=0;
	use_bg_cma=FALSE;
	NullFreqBG=0; ObsBG=0;
	Kingdom=0,Phylum=0;
	insert=18;extend=2;
	MinResFreqCutoff=0.0;	// show all residues by default..
	MaxMisMatches=0;
	ColorString=0;
	MaxConcensusLines=3;
	HideIndelInfo=FALSE;
	pcr = new pcr_typ;
	maxlen_gnull=INT4_MAX;
	if(cmc_res_vals) ModeLPR='U'; else ModeLPR='P';
	InputFractionAsNumber=0;
	Alphas=0; A0s=0; B0s=0;
	ownAB=0;

	Int4	i,j,s,m,n,arg,r;

	InitAsNull();
	// Information for subset construction.
	NumResidues=0; 
	// char	residue_str[30];

	merged_rasmol[0]=0; // i.e., set to FALSE;

	// Get input parameters and file names.
	if(argc < 2) print_error(USAGE_START);
	if(NumCMA > 0){ AB=AlphabetCMSA(CMA_ARRAY[1]); }
        else {
	  // AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // USING THIS BLOSUM IS IMPORTANT!!
	  AB=MkAlphabet(AMINO_ACIDS,GBLAST_BLOSUM62,SREL26_BLSM62); ownAB=AB;
	}
//******************* chn_typ::GetArg(Int4 argc,char *argv[]); *****************
//******************* chn_typ::GetArg(Int4 argc,char *argv[]); *****************
//******************* chn_typ::GetArg(Int4 argc,char *argv[]); *****************
	// sprintf(output_name,"%s",argv[1]); 
	output_name= AllocString(argv[1]);
	assert(isprint(output_name[0])); 
	// fprintf(stderr,"**** output_name(3)=%s ****\n",output_name);
	// for(arg = 0; arg < argc; arg++){ fprintf(stderr,"%s ",argv[arg]); } fprintf(stderr,"\n");
	this->GetArg(argc,argv); // Get the input arguments and set parameters accordingly.
	assert(isprint(output_name[0])); 
	if(cbp_cut == 0.0) cbp_cut=0.01;
	if(ColorString==0){
	  ColorString=AllocString(" CyPBGRcmbgr                                          ");
	}
	if(NumResidues > 0 && MaxMisMatches >= NumResidues){
		print_error("Fatal: MaxMisMatches >= # pattern positions");
	}

	assert(isprint(output_name[0])); 
//****************** Set up alignment format options ************************
//****************** Set up alignment format options ************************
//****************** Set up alignment format options ************************
	// Related pairs used for computing info by cha_typ.
	if(NumResidues > 0){
	  char tmp_name1[300];
	  char tmp_name2[300];
	  // afn: 1/17/08 name change for mismatches...
	 if(MaxMisMatches > 0 && NumResidues > 0){
	     sprintf(tmp_name2,"%s_%dm",output_name,MaxMisMatches); 
	     sprintf(tmp_name1,"%s",tmp_name2);
	 } else sprintf(tmp_name1,"%s",output_name);
	// skips this if NumResidues == 0.
	 for(i = 1; i <= NumResidues; i++){
	   sprintf(tmp_name2,"%s",tmp_name1);
	   if(i < 10){ // use only first patterns in file name.
		sprintf(tmp_name1,"%s_%s%d",tmp_name2,residue_str[i],Position[i]);
	   }
	   if(islower(Residue[i])){
		InvertResSet[i]=TRUE;
		Residue[i]=toupper(Residue[i]);
	   } else InvertResSet[i]=FALSE;
	   m = strlen(residue_str[i]);
	   if(InvertResSet[i]){
	     sst_typ NotThese=0;
	     for(j=0; j<m; j++){
		char aa=residue_str[i][j];
		sst_typ tmp_set=SsetLet(AlphaCode(aa,AB));
		NotThese= UnionSset(NotThese,tmp_set);
	     }
	     Residues[i] = 0;	// Empty set...
	     for(j=0; j <= nAlpha(AB); j++){
		if(!MemSset(j,NotThese)){
		    sst_typ tmp_set=SsetLet(j);
		    Residues[i] = UnionSset(Residues[i],tmp_set);
		}
	     }
	   } else {
	     Residues[i] = 0;	// Empty set...
	     for(j=0; j<m; j++){
		char aa=residue_str[i][j];
		sst_typ tmp_set=SsetLet(AlphaCode(aa,AB));
		Residues[i] = UnionSset(Residues[i],tmp_set);
	     }
	   }
#if 0	// DEBUG...
	   for(j=0; j<= nAlpha(AB); j++){
		if(MemSset(j,Residues[i])) printf("%c",AlphaChar(j,AB));
	   }
	   printf(" == %s%d\n",residue_str[i],Position[i]); // exit(1);
	   printf("\n");
#endif
	  free(output_name);
	  output_name = AllocString(tmp_name1);
	 }
	}
	assert(isprint(output_name[0])); 
	ReadChnAnalMultAln(argc,argv,NumCMA,CMA_ARRAY,USAGE_START);  // in chn_read.cc: read alignments
	assert(isprint(output_name[0])); 
	// GetSWT(NameFReadSWT);		// chn_swt.cc: calc. weighted numbers of sequences, residues
	GetSWT(HSW);		// chn_swt.cc: calc. weighted numbers of sequences, residues
	assert(isprint(output_name[0])); 
	if(!cmc_input){		// not sure why this was/is being created within constructor (afn:11-23-09).
		// Add extract Alphas option here...
		NEW(Alphas,NumberMCMA +3,double);
		NEW(A0s,NumberMCMA +3,Int4);
		NEW(B0s,NumberMCMA +3,Int4);
		NEW(SetMode,NumberMCMA +3,char);
		if(ModeLPR == 'P'){
		   for(Int4 anal=1; anal <= NumberMCMA; anal++){
			Alphas[anal]=GetBPPS_CMA(&A0s[anal],&B0s[anal],&SetMode[anal],MCMA[anal]);
		   }
		} CreateRTF( );	// in chn_rtf.cc: create rich text format objects.
	}	// rtf objects not needed for BPPS routine....only for chn_see.
	assert(isprint(output_name[0])); 
#if 1
//******************* Secondary structure prediction *************************
//******************* Secondary structure prediction *************************
//******************* Secondary structure prediction *************************
	double	predicted_accuracy,**ss_prob,max=0.0;
	Int4	aln,best_ss,t,left_flank=10,right_flank=10;

    for(aln=0; aln<= NumAnalysis; aln++) SecondStrct[aln]=0;
    if(Show2ndary) for(aln=1; aln<= NumAnalysis; aln++){
	if(CMA[aln]==0) continue;
	NEW(SecondStrct[aln],LengthCMSA(1,CMA[aln])+3,char);
	ss_prob=ComputeDSC(NULL,&predicted_accuracy,1,left_flank,right_flank,CMA[aln]);
	for(s=1; s <= LengthCMSA(1,CMA[aln]);s++){
	   for(max=0.0,t=1; t <= 3; t++){
		if(max < ss_prob[t][s]){ max=ss_prob[t][s]; best_ss=t; }
	   } SecondStrct[aln][s]=' ';
	   switch(best_ss){
	     case 1: if(max >= 0.75) SecondStrct[aln][s]='H'; 
		     else if(max >= 0.50) SecondStrct[aln][s]='h'; break;
	     case 2: if(max >= 0.75) SecondStrct[aln][s]='S';
		     else if(max >= 0.50) SecondStrct[aln][s]='s'; break;
	   }
	} for(t=1; t <= 3; t++) free(ss_prob[t]); free(ss_prob);
    }
#endif
	assert(isprint(output_name[0])); 
    	// fprintf(stderr,"**** output_name(4)='%s' ****\n",output_name);
}

void    chn_typ::Free( )
{
	Int4	i,s;

  if(!OutputAln){
	if(NoWeights){
	  for(i=0; i <= NumberMCMA; i++){
	   for(s=1; s<=LenSeq(qE); s++){ free(Obs[i][s]); } free(Obs[i]);
	  }
	} else {
	   for(s=1; s<=LenSeq(qE); s++){ free(Obs[0][s]); } free(Obs[0]); 
	}
        for(i=0; i <=NumAnalysis; i++) // afn: 11/10/07 ... make sure all is freed.
	{  if(status && status[i]) free(status[i]); }
	if(status) free(status);
        for(i=0; i <= MAX_ALN_CHN_TYP; i++){	// arrays are a bit longer than this.
	  if(SecondStrct[i]) free(SecondStrct[i]); 
	  if(NullFreq[i]){
	      if(NullFreq[i][0]){	// indicates that NullFreq[i] arrays should be freed.
		// This is probably never true right now...
		for(s=1; s<=LenSeq(qE); s++){ free(NullFreq[i][s]); }
	      } free(NullFreq[i]); 	
	  } 
	  if(rtf[i]) delete rtf[i]; 
	  if(rtfQ[i]) delete rtfQ[i]; 
	  if(WtNumSeq[i]) free(WtNumSeq[i]);
	} NilSMA(MA);
	if(swt[0]) delete swt[0];	// Integer Weights.
	for(i=1; i <=NumberMCMA; i++) if(swt[i]) delete swt[i]; 
	if(swt) free(swt);
	NilSeq(keyE); NilSeq(qE); NilSeq(full_csq);
	if(TAX_CMA) free(TAX_CMA);
	if(Kingdom) free(Kingdom);
	if(Phylum) free(Phylum);

	// for Integer weights (chn_swt.cc)
	if(SqIWtMain) free(SqIWtMain);
	if(MainCMA){
	  if(MainSqWt){
	    for(s=1; s<=LengthCMSA(1,MainCMA); s++) free(MainSqWt[s]); free(MainSqWt); 
	  }
	  if(OwnCMAs) TotalNilCMSA(MainCMA);
	}
	if(Alphas) free(Alphas);
	if(A0s) free(A0s);
	if(B0s) free(B0s);
	if(SetMode) free(SetMode);
  }
	if(Xconserved) free(Xconserved);
	if(output_name) free(output_name);
  	if(BG_CMA) TotalNilCMSA(BG_CMA);
	if(StdAlignmentOnly){
		if(OwnCMAs) NilCMSA(IN_CMA[1]);
		if(OwnCMAs) TotalNilCMSA(IN_CMA[2]);
	} else if(OwnCMAs) for(i=1; i<= Number; i++) TotalNilCMSA(IN_CMA[i]);
	free(IN_CMA); free(CMA);
	if(pcr) delete pcr;
	if(ColorString) free(ColorString);

	if(ownAB) NilAlpha(ownAB); 
}


