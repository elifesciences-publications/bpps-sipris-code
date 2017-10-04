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

#include "gsm_typ.h"

static char PROBE_USAGE[]= "\nUsage: probe fastafile database [options] \n\
  options (range:default):\n\
\t-A           - start by searching database with input file fastafile.msa \n\
\t-a<int>      - algebraic operations (0-9:1)\n\
\t               '0' = fast & dirty mode (no refinement)\n\
\t               '1' = fast mode (split & fuse once)\n\
\t               '2' = refine mode (do refinement steps)\n\
\t               '3' = rigorous mode (careful, time-consuming refinement)\n\
\t               '4' = experimental (time-consuming refinement)\n\
\t               '10'-'14' = experimental gapped version of probe.\n\
\t               '5'-'9' = not yet implemented.\n\
\t-B<char>     - residue frequencies to use for background model\n\
\t                m = use multiple alignment background counts.\n\
\t                s = use standard background counts.\n\
\t                d = use database background counts.\n\
\t                (default = use standard background counts).\n\
\t-b<int>      - initial mean number of blocks (1-10:5)\n\
\t-C<float>    - minimum 'field map' Cutoff for individual blocks (0.0)\n\
\t-D           - don't delete blocks before adding to cma population\n\
\t-d<float>    - maximum E-value for repeat domain detection\n\
\t-E<float>    - scan single-block E-value cutoff (0-50000:1)\n\
\t-e<float>    - scan E-value cutoff (0-50000:0.01)\n\
\t-f<int1:int2>- fix number of blocks to within the rangeint1-int2.\n\
\t-G<file>     - use guide sequence given in <file>\n\
\t-G=<int>,<int>- perform goscan gapped search with open & extend penalties (suggested: 18,2)\n\
\t-g<int>:<int>- specify gap opening and extension penalties (default: 10000,2000)\n\
\t-H           - create a histogram file (*.gaps) of gap lengths\n\
\t-h<int>      - heapsize for msa population (2-1000:10)\n\
\t-I<int>:<int>- save internal repeats with flanking regions of given lengths\n\
\t-i<int>      - purge score decrement (1-100:5)\n\
\t-L           - Don't mask low complexity regions\n\
\t-l<int>      - set rapid convergence limit (higher = longer to converge)\n\
\t-M<int>:<int>- minimum & maximum number of sequences in alignment\n\
\t-m<char>     - specify scan method (cCdDgGrR:c)\n\
\t               'g' or 'G' = Gribskov's (fragment or don't fragment)\n\
\t               'b' or 'B' = modified Gribskov's with motif background\n\
\t               'm' or 'M' = modified Gribskov's\n\
\t               'd' or 'D' = Dirichlet mixture priors\n\
\t               'r' or 'R' = product multinomial\n\
\t               'h' or 'H' = Henikoff's method\n\
\t-N<int>      - maximum number of repeats for scan step (1-1000:1)\n\
\t-n           - mask potential nonglobular regions using 'seg x 45 3.4 3.75'\n\
\t-n<int>      - minimum number of blocks with good matches in search phase (default: all blks)\n\
\t-O           - breed using input file fastafile.msa\n\
\t-O<int>      - Create an alignment from database search step to start Gibbs sampling\n\
\t                Purge the alignment at a cutoff of <int> percent identity\n\
\t-o<int>      - improve (optimize) by removing segments pulling down map\n\
\t-P<int>      - Patience in seconds for breed time (1-900000:600)\n\
\t-p<int>      - msa population size (2-10000:10)\n\
\t-R<int>      - maximum number of model refinement cycles (1-100:10)\n\
\t-r<int>      - minimum number of repeats to save a hit (1-1000:1)\n\
\t-S<int>      - purge maximum cutoff score (20-500:150)\n\
\t-s<int>      - random seed\n\
\t-T<real>     - Use simulated annealing starting with temperature T (0-300:nil)\n\
\t-t           - Use a template cmsa file to start alignment (0-300:nil)\n\
\t-trim=<real> - trim the cmsa file at <real> bits of information prior to a search\n\
\t-U           - use the repset option\n\
\t-U+          - use the repset option with priority on labeled sequences\n\
\t-u<char>     - use the scan mode specified by <char>\n\
\t                'g' = global with gap function\n\
\t                'G' = global without gap function (default)\n\
\t                'c' = local core with gap function\n\
\t                'C' = local core without gap function\n\
\t                'o' = local overlap with gaps\n\
\t                'O' = local overlap without gaps\n\
\t-v<int>      - maximum variation in map for convergence in optbreed (default = 1.0)\n\
\t-W           - output a *.see file to plot progress of sampler\n\
\t-W<int>      - Maximum starting motif width (default: 25)\n\
\t-w           - DON'T use sequence weighting in scan step\n\
\t-x           - dummy\n\n";

// \t-minlen=<int>- minimum length of sequences in alignment (used only with -G= option)\n\

gsm_typ::gsm_typ(int argc, char *argv[], a_type alphabet)
{ Init(0,argc, argv, alphabet); }

gsm_typ::gsm_typ(char *PRB_USAGE,int argc, char *argv[], a_type alphabet)
{ Init(PRB_USAGE,argc, argv, alphabet); }

gsm_typ::gsm_typ(char *PRB_USAGE,int argc, char *argv[], a_type alphabet, ss_type Data)
{ Init(PRB_USAGE,argc, argv, alphabet,Data); }

void    gsm_typ::InitDefaults()
{
	seed=18364592;
	minseq=5; number=-1; maxbreed=40; 
	min_rpt=1; breedVar=1.0; maxseq=10000; inc=5; maxseqleng = 0;

	// aveblk=20; avecol=6; minblk=12; maxblk=20; fix=TRUE;
	aveblk=5; avecol=15; minblk=3; maxblk=10; 
	mhpsz=3;
	TestMode=FALSE;
	MaxIter=4; // maximum of 5 rounds..

	SmplBlksCols=TRUE;
	maxrun=10; cutoff=150; 
	align_mode = 1; // default alignment method.
	maxrpts=1; maxcycle=1;
	left_flank=500;right_flank=500;
	// temperature = 300;  // 1/kT = 1.0 (T = 300 K)
	temperature = 150;  // 1/kT = 1/k*150 K = 2.0
	method='H';	// 'H' = Henikoff's not fragmented 
	patience=600;	// willing to wait 10 minutes by default
	// patience=1800;	// willing to wait 30 minutes by default
	pseudo=0.5;Ecut=1.0;ecut=0.01; minmap=0.0;
	bestmap=-9999999999.;
	mode='G'; 	// scan mode...
	limit=5;
	repeatEval=0.01;
	aafreq='s';

	pernats=1000.0; gapopen=10000; // = log(freq_open=0.05)*pernats = -30.
	gapextend=2000; 		//  = log(freq_extend=0.25)*pernats = -5
	indel_penalty=0.0;

	use_gseq=TRUE;
	MaxBadBlks=0;
	MinGoodBlks=1000;
	improve_cut=0.0;
	max_motif_width=25;
	COL=5; BLK=0;
}

void    gsm_typ::InitAsNull()
{
	time1=time(NULL); 
	//**************** new gapped search options *****************
	write=FALSE;
	StrtLens=0; StrtBlks=0;
	do_gapped_srch=FALSE;
	MinSeqLen=0;
	trim_info=0.0;
	//**************** end gapped search options *****************
	gss=0;
	AB=NULL;
	DBS_hits=NULL; 
	guide=NULL; Guide=NULL; WG=NULL; maH=NULL;
	cmsa=NULL; bestmsa=NULL; cmsa_in=NULL; data=NULL; in_data=NULL; input_data=FALSE;
	nsize=NULL;  counts=NULL;
	num=0;

	go=TRUE; weight=TRUE;

	report_gaps=FALSE; mask_nonglobular=FALSE;
	flag=FALSE;noseg=FALSE; combine=FALSE;
	repeats=FALSE; input_msa=FALSE; force=FALSE;fix=FALSE;
        useSA=FALSE; // flag to indicate use of sim. annealing option
	use_breed = FALSE; UseLabeled=FALSE; 
        UseRepset=FALSE; improve=FALSE; dont_delete=FALSE;

	cardB=0;oldcardB=0;
	cma_purge=0;	// don't use this by default.
	Run=0;
}

void	gsm_typ::Init(char *INPUT_USAGE,int argc, char *argv[], a_type alphabet, ss_type Data)
// Initialize all items for gsm_typ alignment and searches 
{
	Int4	i;
	char *PRB_USAGE=0;
	FILE	*fptr;

	Argv[0]=0;
	if(INPUT_USAGE == 0){ // for probe...
		PRB_USAGE=PROBE_USAGE;
    		if(argc < 3) print_error(PRB_USAGE);
		Argv[2] = AllocString(argv[2]);
    		if(!(argc == 3 || argv[3][0] == '-')) print_error(PRB_USAGE);
	} else {		// for gismo;
		PRB_USAGE=INPUT_USAGE;
    		if(argc < 2) print_error(PRB_USAGE);
    		if(!(argc == 2 || argv[2][0] == '-')) print_error(PRB_USAGE);
		Argv[2] = AllocString(argv[1]);
	}
	if(argc > 150) print_error(PRB_USAGE);
    	if(Data==0){
	   if((fptr=fopen(Argv[2],"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",Argv[2]);
		print_error("File does not exist!");
	   } fclose(fptr);
	}
	InitAsNull(); InitDefaults();	// initialize with default values...
	Argv[1] = AllocString(argv[1]);
        if(INPUT_USAGE){ limit=1; align_mode=1; }
	NEW(name,200,char); strcpy(name,argv[1]); AB=alphabet; 

// initialize with user input 
   BooLean	DoNotSeed=FALSE;
   Int4 arg_start=3;
   if(INPUT_USAGE) arg_start=2;
   for(i=arg_start; i < argc; i++){
     if(argv[i][0] != '-') print_error(PRB_USAGE);;
     switch(argv[i][1]) {
	case 'A': input_msa = TRUE; break;
	case 'a': align_mode=IntOption(argv[i],'a',0,19,PRB_USAGE); 
		  if(align_mode >= 10){
// fprintf(stderr,"FATAL: -a >10 option inactivated!\n"); print_error(PRB_USAGE);
			// input_msa = TRUE; 
			align_mode -= 10; use_gseq=TRUE; 
		  }
		  if(INPUT_USAGE && align_mode > 7) print_error(PRB_USAGE);
		  break;
	case 'b': aveblk=IntOption(argv[i],'b',1,100,PRB_USAGE); break;
	case 'B': if(!isalpha(aafreq=argv[i][2])) 
			print_error(PRB_USAGE); break;
	case 'C': minmap=RealOption(argv[i],'C',-5000,5000,PRB_USAGE); break;
	case 'D': 
		if(strcmp("-DoNotSeed",argv[i]) == 0) DoNotSeed=TRUE;
		else dont_delete=TRUE; break;
	case 'd': repeatEval=RealOption(argv[i],'d',0,500000,PRB_USAGE); break;
	case 'e': ecut=RealOption(argv[i],'e',0,5000,PRB_USAGE); break;
	case 'E': Ecut=RealOption(argv[i],'E',0,5000,PRB_USAGE); break;
	case 'F': force= TRUE; break;
	case 'f': fix = TRUE; 
		  if(sscanf(argv[i],"-f%d:%d",&minblk,&maxblk) != 2 ||
			minblk > maxblk || minblk < 1)
                        print_error(PRB_USAGE);
		  break;
	case 'H': report_gaps = TRUE; break;
	case 'h': mhpsz=IntOption(argv[i],'h',2,1000,PRB_USAGE); break;
	case 'G': 
        	   if(sscanf(argv[i],"-G=%d,%d",&Gap_o,&Gap_x) == 2){
			do_gapped_srch=TRUE;
		   } else if(argv[i][2] == 0) {
			print_error(PRB_USAGE);
		   } else {
		      guide=AllocString(argv[i]+2); 
		   } break;
        case 'g': if(sscanf(argv[i],"-g%d:%d",&gapopen,&gapextend) != 2)
                                        print_error(PRB_USAGE);
		  if(gapopen== 0 && gapextend==0)use_gseq=FALSE;
                  if(gapopen < 0 || gapextend < 0) print_error(PRB_USAGE); break;
	case 'i': inc=IntOption(argv[i],'i',1,100,PRB_USAGE); break;
	case 'I': repeats=TRUE; 
		  if(sscanf(argv[i],"-I%d:%d",&left_flank,&right_flank) != 2)
                        print_error(PRB_USAGE);
		  break;
	case 'L': noseg=TRUE; break;
	case 'l': limit=IntOption(argv[i],'l',1,500,PRB_USAGE); break;
	case 'M': if(sscanf(argv[i],"-M%d:%d",&minseq,&maxseq) != 2 
			|| minseq > maxseq || minseq < 2) print_error(PRB_USAGE);
		  break;
	case 'm': 
#if 0
		if(sscanf(argv[i],"-minlen=%d",&MinSeqLen) == 1){ 
			if(MinSeqLen < 2) print_error(PRB_USAGE);
		} else 
#endif
		if(!isalpha(method=argv[i][2])) print_error(PRB_USAGE); 
		break;
	case 'N': maxrpts=IntOption(argv[i],'N',1,1000,PRB_USAGE); break;
	case 'n': {
		 if(argv[i][2]==0) mask_nonglobular=TRUE; 
	         else MinGoodBlks=IntOption(argv[i],'n',0,1000,PRB_USAGE);
		} break;
	case 'O': use_breed = TRUE; 
		if(argv[i][2]) cma_purge=IntOption(argv[i],'O',10,100,PRB_USAGE);
		break;
	case 'o': improve = TRUE; 
		  improve_cut = RealOption(argv[i],'o',-10000.,10000.,PRB_USAGE); break;
	case 'p': maxcycle=IntOption(argv[i],'p',2,10000,PRB_USAGE); break;
	case 'P': patience=IntOption(argv[i],'P',1,9999999,PRB_USAGE); break;
	case 'R': maxrun=IntOption(argv[i],'R',0,100,PRB_USAGE); break;
	case 'r': min_rpt=IntOption(argv[i],'r',1,1000,PRB_USAGE); break;
	case 'S': cutoff=IntOption(argv[i],'S',20,5000,PRB_USAGE); break;
	case 's': seed=atoi(argv[i]+2); break;
	case 'T': useSA=TRUE; temperature=RealOption(argv[i],'T',0,5000.,PRB_USAGE); break;
	case 't': 
		if(sscanf(argv[i],"-trim=%lf",&trim_info) == 1){ 
			if(trim_info <= 0.0) print_error(PRB_USAGE);
		} else {
		   if(cmsa_in != NULL) break;
		   sprintf(str,"%s.cma",name); std::cerr << str; cmsa_in=ReadCMSA2(str,AB); 
		} break;
        case 'U': UseRepset=TRUE; 
		if(argv[i][2] == '+') UseLabeled=TRUE; 
		break;
        case 'u': if(!isalpha(mode=argv[i][2])) print_error(PRB_USAGE); 
                  if(argv[i][3] == '+') combine=TRUE;
                  else if(argv[i][3] == '\0') combine=FALSE;
                  else print_error(PRB_USAGE); break;
	case 'v': breedVar=RealOption(argv[i],'v',0.0,10000.0,PRB_USAGE); break;
	case 'w': weight=FALSE; break;
	case 'W': 
		if(argv[i][2] != 0){
			max_motif_width=IntOption(argv[i],'W',8,MAX_LENG_SITES,PRB_USAGE); 
		} else if(WG==NULL) WG=MakeWatchGibbs(100, 1000); 
		break;
	case ' ': break;  // ignore blanked out options.
	default: print_error(PRB_USAGE); break;
     }
   }
   if(minseq > maxseq) print_error(PRB_USAGE);
   if(cma_purge > 0 && maxrpts > 1) print_error("-O<int> option not yet implemented for repeats");
   maxcycle=maxcycle*mhpsz;
   maH=MkMSAHeap(mhpsz);

   // gapopen*=pernats; gapextend*=pernats; 

   // fptr = open_file(argv[1],".cmd","w");
   // for(i = 0; i < argc; i++) { if(argv[i][1] != ' ') fprintf(fptr,"%s ",argv[i]); }
   if(seed == 18364592) {  // not provided by user 
	seed = (UInt4) time(NULL); //fprintf(fptr,"-s%d\n",seed); 
   } // else fprintf(fptr,"\n"); fclose(fptr);

   if(Data != 0){ input_data=TRUE;  in_data=Data; }
   if(DoNotSeed == FALSE) sRandom(seed);
   if(guide != NULL) { Guide = MkGuide(guide, AB); PutGuide(stderr,Guide, AB); free(guide); }
}

void gsm_typ::read_input_data()
{
   // free up previous alignment at this point... (Stored as a file).
   if(cmsa != NULL){ NilCMSA(cmsa); cmsa= NULL; }
   if(bestmsa!=NULL) { NilCMSA(bestmsa);  bestmsa=NULL; }
   assert(data==NULL); // if(data!=NULL) { NilSeqSet(data);  data=NULL; }
   BLK=0;

   if(in_data){ data=in_data; in_data=0; }
   else if(!noseg) data = MkXSeqSet(name,AB); else data = MkSeqSet(name,AB);
   // gss.initialize(data,gapopen,gapextend,pernats,left_flank,right_flank);
   gss = new gss_typ(data,gapopen,gapextend,pernats,left_flank,right_flank);
   SetIndelPenaltySeqSet(indel_penalty,data);
   make_binomials(data);
}

void gsm_typ::make_binomials(ss_type sqset)
{
   maxlenModels = MinSeqSeqSet(sqset);
   // if(avecol == 0) avecol = (Int4) ((float) maxlenModels/(float)(aveblk*2))+1;
   // NEW(num,maxlenModels+2,Int4);
   if(num) free(num); NEW(num,aveblk*2+2,Int4);
   if(minseq > NSeqsSeqSet(sqset)){
	fprintf(stderr,"Only %d sequences in input set (at least %d required)\n",NSeqsSeqSet(sqset),minseq);
// assert(minseq <= NSeqsSeqSet(sqset));
	print_error("FATAL: input set has too few sequences");
   } bestmap=-9999999999.;
}

void gsm_typ::free_input_data()
{
     if(WG != NULL) { NilWatchGibbs(WG); WG=MakeWatchGibbs(100, 1000); }
     if(num) free(num); num=0; 
     if(gss) delete gss; gss=0;
}

gsm_typ::~gsm_typ()
{ 
    // afn: 11-13-2014
    if(gss) free_input_data();
    if(maH){	// afn: 10-24-2014.
      cma_typ xcma;
      while((xcma=DelMaxMSAHeap(&map, maH)) != 0) NilCMSA(xcma);  NilMSAHeap(maH); maH=NULL;
    }
    if(StrtLens) free(StrtLens);
    if(counts != NULL) free(counts); 
    if(nsize != NULL) free(nsize); 
    if(Guide != NULL) NilGuide(Guide);
    if(DBS_hits != NULL) NilSet(DBS_hits);
    free(name); 
    free(Argv[1]); free(Argv[2]); 
    if(WG != NULL) { NilWatchGibbs(WG); WG=NULL; }
    if(cmsa != NULL) NilCMSA(cmsa); 
    if(bestmsa != NULL) NilCMSA(bestmsa); 
    if(cmsa_in != NULL){
    	gss_typ *gssX=gssCMSA(cmsa_in); 	// not owned by cmsa_in!!!
    	gssX->~gss_typ();
	NilCMSA(cmsa_in); 
    }
    if(!input_data && data != NULL) NilSeqSet(data);
    // if(A != NULL) NilAlpha(A);
    double runtime=difftime(time(NULL),time1);
    fprintf(stderr,"\ttime: %.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
}

