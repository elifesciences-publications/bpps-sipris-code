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

#include "hat_typ.h"

#if 0
void    PutOneSeqAlign(FILE *fp, sap_typ sap, Int4 width, a_type AB)
// put only the first sequence's alignments
{
	UInt4   last_id=0;
	sap_typ	gsap=0;
        dsp_typ dsp=sap->segs; last_id=dsp->subject_id;
        fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE);
        fprintf(fp,"\t\tlength = %d\n\n",LenSeq(dsp->sE));
	for(gsap=sap; gsap!=NULL; gsap=gsap->next){
          dsp=gsap->segs;
          if(dsp->subject_id == last_id){
		PutGSeqAlign(fp,gsap,width,AB);
          } else break;
        } fprintf(fp,"\n");
}

sap_typ	MinEvalSAP(sap_typ head, double *MinEval)
// set MinEval to the minimum evalue for the first subject sequence in head.
// return the sap for the next sequence or null if end of list
{
	sap_typ sap=0;
	UInt4	id=0;
	double	eval,min=DBL_MAX;

	if(!head){ *MinEval=DBL_MAX; return 0; }
	dsp_typ dsp=head->segs; id = dsp->subject_id;
	for(sap=head; sap; sap=sap->next){
		dsp=sap->segs;
		if(dsp->subject_id == id){
			eval = EvalueGSAP(sap);
			if(eval < min) min = eval;
		} else { *MinEval = min; return sap; } // found another sap;
	}
	assert(sap == 0);
	*MinEval = min;	// reached end of list; sap == 0;
	return sap;
}
#endif

#define USAGE "usage: gapmap query_prefix dbsfile [options]\n\
      NOTE: Input multiple profile file = query.chk (created using mkmaps)\n\
            Input Template alignment file = query.cma \n\
      options:\n\
        -C=<int>  Cutoff score for eliminating overlaps (-P option)(default: 30)\n\
        -D        Don't SEG query sequence\n\
        -e=<real> E-value cutoff for showing scores and alignments\n\
        -prune=<char> Basis for input set pruning (default: 'T')\n\
                   'T': Template HSP trigger (slowest & most sensitive).\n\
                   'E': Template extensions (intermdiate).\n\
                   'H': Template hits (fastest & least sensitive).\n\
        -H=<int>  heapsize (default 500000)\n\
        -I<int>:<int> - flanking region lengths for sequence hits(default: 0:0)\n\
        -L        do low complexity masking of database sequences\n\
        -M        Output cma files for iterative recycling \n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -m=[<char>:<int>..<int>;]<char>:<int>..<int>.  generate output for a vsi file\n\
	          (e.g., -m=R:4..41,O:54..76,Y:81..104,G:112..137,B:141..156.)\n\
        -O        output sequence hits in separate files for each profile\n\
        -o=<real> output sequence hits in separate files only if E-value is above <real>\n\
        -P        Don't use PatchSAP() option to combine gapped HSPs into one\n\
        -r        Don't detect repetitive domains\n\
        -screen   pre-search database with file query0.cma to obtain searchable subset\n\
        -sense=<real>  sensitivity parameter for heuristic search \n\
		   range: 0.01...1.0; lower values: > sensitivity (default: 0.90)\n\
        -T=<int>  blast word hit threshold (default 11)\n\
        -t<real>  trim alignments at ends where marg. prob. is > <real>\n\
        -trigger=<real>  Number of bits to trigger gapping = <real> (default: 22.0)\n\
        -test=<int>  test mode (default: 0)\n\
                   test=0:  don't allow deletions next to insertions\n\
                   test=1:  allow deletions next to insertions\n\
        -w=<int>  alignment width (default 60)\n\
        -X=<int>  X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
	-z=<int>  set effective length of the database (use zero for the real size)\n\
        -z        see posMatrix\n\
\n"


#define MAX_NUM_CMA_FILES 1000
#if 0	//**********************************************************************
#endif  //**********************************************************************
int	gapmap_srch(int argc, char *argv[])
{
	Int4    T=11,gap_open=11,gap_extend=1;
	time_t	time1=time(NULL);
	Int4	left_flank=0,right_flank=0;
        double	x_parameter=25.0;
	Int4	width=60,x;
	UInt4	hpsz=1000000;
	Int4	printopt=0;
	double	Ethresh=0.001,Ecutoff=0.001;
	BooLean	see_smatrix=FALSE,use_chk_files_only=FALSE;
	BooLean	seg_seq=TRUE,success=FALSE,use_patch_sap=TRUE;
	BooLean	output_by_class=FALSE,out_all_class=FALSE;
	BooLean print_names_only=FALSE,find_repeats=TRUE;
	Int4	*IGNORE_THIS=0,num_bgd_files=0;
	char	mask_dbs=' ',*checkin=0,*checkout=0;
	FILE	*seqfile=NULL;
	char	*cmafile=0;
	double	misalncut=2.0;
	Int4	num_cma_files,overlap_cutoff=30,start_arg=3;
	UInt4	max_length_input_seqs=300000;
	Int4	NumRegions=0;
	Int4	*Start=0,*End=0,test=0;
	char	*colors;
	BooLean	SkipFile=FALSE, use_weak_cutoff= FALSE; 
	double	fractionMS=0.90,WeakEcutoff=0.0;
	char 	str[100],*Checkpoint,prune_mode='T';
	Int4	num_false_pos=0,first_false_pos;
	BooLean	OutputMAPS=FALSE;
	double	gap_trigger_bits=0.0;
	Int4	dbslen=0;
	Int4	Begin=1;	// Begin with template (II=0) or with first profile (II=1).
	// Leave this at 1 for now; can always add a copy of template to mpa file.
#if 0	// new false positive input files
	Int4	NumFalsePos=0,=TrueCMA_Files=0,FalseCMA_Files=0,TotalNumCMA_Files=0;
	cma_typ	*ALL_CMA,*TRUE_CMA,*FALSE_CMA;
	Int4	Num
#endif

	if(argc < 3) print_error(USAGE);
	//***************************** 1. Read input files *******************************
	start_arg=3;
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	cma_typ	TemplateCMA=0;

  	sprintf(str,"%s.tpl",argv[1]); 
	TemplateCMA=ReadCMSA2(str,A);
	num_false_pos=LevelCMSA(TemplateCMA);
	num_cma_files = NumSeqsCMSA(TemplateCMA) - 1;
	first_false_pos = num_cma_files - num_false_pos + 1;
	for(Int4 arg=0; arg < argc; arg++) { fprintf(stderr,"%s ",argv[arg]); }
	fprintf(stderr,"\n\n");
	// fprintf(stderr,"num_cma_files = %d\n",num_cma_files);
	//***************************** 2. Get input options *******************************
        for(Int4 arg = start_arg; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(USAGE);
           switch(argv[arg][1]) {
	     case 'C': overlap_cutoff=IntOption(argv[arg],'C',1,100,USAGE); break;
	     case 'D': seg_seq=FALSE; break;
             case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,USAGE); break;
	     case 'H': hpsz=IntOption(argv[arg],'H',1,100000000,USAGE); break;
             case 'h': Ethresh=RealOption(argv[arg],'h',0.0,10000,USAGE); break;
	     case 'I': // putsegs=TRUE;
		if(sscanf(argv[arg],"-I%d:%d",&left_flank,&right_flank) != 2)
                        print_error(USAGE); break;
             case 'L': 
		if(argv[arg][2] != '='){ mask_dbs='x'; }
		break;
	     case 'm': 	// R:7..41,O:54..76,Y:81..104,G:114..137,B:141..156.
		{
		  if(argv[arg][2] != '='){
			printopt=IntOption(argv[arg],'m',0,6,USAGE); 
		  } else {
			NumRegions=ParseInputRegionsVSI(argv[arg]+3,
				&Start,&End,&colors,USAGE);
			printopt=2;
		  }
		} break;
	     case 'M': OutputMAPS=TRUE; argv[arg][1] = ' '; break;
             case 'P': use_patch_sap=FALSE; argv[arg][1] = ' '; break;
	     case 'p': 
		if(sscanf(argv[arg],"-prune=%c",&prune_mode) != 1){
                        print_error(USAGE); 
		} else if(prune_mode != 'H' && prune_mode != 'E' && prune_mode != 'T'){
			print_error("prune_mode input error");
		} break;
             case 'r': find_repeats=FALSE; argv[arg][1] = ' '; break;
             case 's': 
             	if(sscanf(argv[arg],"-sense=%lf",&fractionMS) == 1){
			if(fractionMS < 0.01 || fractionMS > 1.0) print_error(USAGE);
			// else fprintf(stderr,"sensitivity = %f\n",fractionMS);
		} else if(strcmp(argv[arg],"-screen") == 0){
			SkipFile=TRUE;
		} else print_error(USAGE);
			break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,USAGE); break;
	     case 't': 
		if(sscanf(argv[arg],"-trigger=%lf",&gap_trigger_bits) == 1){
			if(gap_trigger_bits < 1.0 || gap_trigger_bits > 1000.0){
				print_error("gap_trigger out of permissible range");
			}
		} else if(sscanf(argv[arg],"-test=%d",&test) == 1){
		   if(test > 1 || test < 0) print_error(USAGE);
		} else misalncut=RealOption(argv[arg],'t',0.1,1.0,USAGE);
		break;
             case 'O': {
		  out_all_class=TRUE; output_by_class= TRUE; argv[arg][1] = ' '; 
		} break;
             case 'o': {
             		WeakEcutoff=RealOption(argv[arg],'o',0.0,10000,USAGE);
			// fprintf(stderr,"WeakEcutoff = %g\n",WeakEcutoff);
			output_by_class= TRUE; use_weak_cutoff= TRUE; argv[arg][1] = ' '; 
		} break;
	     case 'w': width=IntOption(argv[arg],'w',5,500000,USAGE); break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,USAGE); break;
             case 'Z': print_names_only=TRUE; break;
             case 'z': 
		if(sscanf(argv[arg],"-z=%ld",&dbslen) == 1){
			if(dbslen <= 0) print_error(USAGE);
			argv[arg][1] = ' ';
		} else see_smatrix= TRUE; break;
             default: print_error(USAGE);
           }
        }

	//***************************** 2. Initialize search *******************************
	e_type	qE,QueryE[MAX_NUM_CMA_FILES],trueE,fakeE;
	ss_type	data,Data;
	sap_typ	sap,SAP[MAX_NUM_CMA_FILES];
	// cma_typ	cma,CMA[MAX_NUM_CMA_FILES];
	cma_typ	cma,*CMA;
	Int4	II,maxrounds;

  if(!IGNORE_THIS) NEW(IGNORE_THIS,num_cma_files+3,Int4);
  for(Int4 n=first_false_pos; n <= num_cma_files; n++) IGNORE_THIS[n]=TRUE;
  for(II=0; II <= num_cma_files; II++){ QueryE[II]=0; SAP[II]=0; }
  // ReNumberSeqSet(Data);
  // Search database for hits using distinct queries for each cma file.
  QueryE[0]=TrueSeqCMSA(1,TemplateCMA); 
  for(II=1; II <= num_cma_files; II++){
	QueryE[II]=TrueSeqCMSA(II+1,TemplateCMA); 
	if(print_names_only){
		fprintf(stdout,"=====================> File %d: ",II);
		PutSeqID(stdout,QueryE[II]);
		fprintf(stdout,"\n");
	}
  }
  if(print_names_only){ fprintf(stdout,"\n"); exit(0); }
  Data=MakeSeqSet(argv[2],max_length_input_seqs,A);
  sprintf(str,"%s.map",argv[1]); Checkpoint = AllocString(str);

	//******************************** 3. Main search loop  *******************************
  FILE 	*chkfp=0;
  checkout=0; checkin=Checkpoint;
  chkfp=open_file(Checkpoint,"","r");
//***************************** TESTING Heuristics ********************************
  set_typ	SkipSet=0,UnionSetH=0,FalseSetH=0;
  set_typ	*SetM,*SetH,*SetE;
  NEW(SetM,num_cma_files+5,set_typ); NEW(SetH,num_cma_files+5,set_typ);
  NEW(SetE,num_cma_files+5,set_typ);
//***************************** TESTING Heuristics ********************************
  for(II=0; II <= num_cma_files; ){
	clock_t timeS=clock();
  // start main loop over all *cma files....
	cma=0;
	qE=CopySeq(QueryE[II]); 
	// PutSeq(stderr,qE,A);
	EqSeqI(0,qE); maxrounds=1;
	gpsi_type gpsi(qE,Data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,chkfp,dbslen);
	if(gap_trigger_bits > 0.0) gpsi.SetGapTriggerInBits(gap_trigger_bits);	// default: 22.0;
//***************************** TESTING ********************************
#if 1	// create skip file using query0.cma file with -screen option.
	if(II == 0){
	   // gpsi.CreateMarginalSet(); 
	} else if(II >= first_false_pos) {
		if(II == first_false_pos){
			fprintf(stderr,
			  "--(false +)-------------------------------------------------\n");
		}
		gpsi.ProvideSkipSet(UnionSetH);	 // only look at sequence in UnionSetH
	} else { gpsi.ProvideSkipSet(SkipSet); }
	gpsi.CreateMarginalSet(fractionMS); 
	gpsi.CreateExtendSet(); 
	gpsi.CreateHitSet(); 
	gpsi.KeepQuite();
#endif
//***************************** TESTING ********************************
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,A);
	Int4 **mtx=0;
	if(use_patch_sap) mtx=gpsi.CopyOFposMatrix();
	// fprintf(stderr,"II=%d; pass=%d\n",II,gpsi.ThisPass( ));
	   // if(gpsi.ThisPass( ) == 0) gpsi.SpeakUp(); 
	if(misalncut < 1.0) sap = gpsi.TrimmedSearch(T,misalncut,mask_dbs);
	else sap = gpsi.Search(T,mask_dbs);
	if(sap) success=TRUE; // else break;
	gpsi.ComputeMatrix(checkout);
	if(see_smatrix && gpsi.posMatrix){ 
		FILE *ofp = open_file(argv[1],".psm","w");
		gpsi.SeeMatrix(ofp); 
		fclose(ofp);
	}
//***************************** TESTING ********************************
	   char id_str[15];
	   if(II == 0){
		SetM[II]=gpsi.RtnMarginalSet();
	   	SetH[II] = gpsi.RtnHitSet(); SetE[II] = gpsi.RtnExtendSet();
		switch(prune_mode){
		  case 'H': SkipSet=SetH[II]; break;
		  case 'E': SkipSet=SetE[II]; break;
		  case 'T': SkipSet=SetM[II]; break;
		  default: print_error("prune_mode input error"); break;
		}
	        UnionSetH = MakeSet(SetN(SkipSet));
		CopySet(UnionSetH,SetH[II]);
	        FalseSetH = MakeSet(SetN(SkipSet));
		Int4 NumPruned=CardSet(SkipSet);
		if(NumPruned == 0){ fprintf(stderr,"no hits found\n"); exit(0); }
		fprintf(stderr,"   pruning input set from %d to %d sequences (%.1lf%%)\n",
			NSeqsSeqSet(Data),NumPruned,
				100.0*(double)NumPruned/(double)NSeqsSeqSet(Data));
	   	fprintf(stderr,"No.: %-10s     %-11s   %-11s   %-11s   (time)(TotalHits)\n",
			"name","Hits","Extend","Trigger");
	   	// fprintf(stderr,"%3d(%-2d): ",II,LevelCMSA(TemplateCMA[II])); 
	   	fprintf(stderr,"%3d: ",II); 
		// StrSeqID(id_str, 10, qE); 
	   	// fprintf(stderr,"%-10s     %-11d   %-11d   %-11d",id_str,
	   	fprintf(stderr,"%-10s     %-11d   %-11d   %-11d","template",
			CardSet(SetH[II]),CardSet(SetE[II]),CardSet(SetM[II]));
   	   	fprintf(stderr," (%0.2f s)\n",
				(double) (clock()-timeS)/(double)CLOCKS_PER_SEC);
	      if(Begin == 1){	// then print out and continue;
		II++; 
		continue; 	// don't use these hits in final alignment...
	      } 
	   } else {
		SetM[II] = gpsi.RtnMarginalSet();
	   	SetH[II] = gpsi.RtnHitSet(); SetE[II] = gpsi.RtnExtendSet();
	   	// PutSet(stderr,SkipSet); exit(1);
	   	// fprintf(stderr,"%3d(%-2d): ",II,LevelCMSA(TemplateCMA[II])); 
	   	fprintf(stderr,"%3d: ",II); 
	   	StrSeqID(id_str, 10, qE); 
	   	fprintf(stderr,"%-10s     %-11d   %-11d   %-11d",id_str,
			CardSet(SetH[II]),CardSet(SetE[II]),CardSet(SetM[II]));
#if 0
   	   	fprintf(stderr," (%0.2f s)\n",
				(double) (clock()-timeS)/(double)CLOCKS_PER_SEC);
#else
	   	if(II < first_false_pos){
	   		UnionSet(UnionSetH, SetH[II]);	// key track of true hits.
	   	} else {
	   		UnionSet(FalseSetH, SetH[II]);	// key track of true hits.
	   	}
   	   	fprintf(stderr," (%0.2f s)(%d)\n",
				(double) (clock()-timeS)/(double)CLOCKS_PER_SEC,
				CardSet(UnionSetH));
#endif 
	   }
//***************************** TESTING ********************************
	   sap_typ head;
	   if(use_patch_sap){ // fix SAP here...to patch together HSPs
		sap=gpsi.StealSAP();
		if(sap){
		   head=MakeSelfAlignedGSAP(qE);
		   head->next = sap;
		   SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data),qE,A); 
		
		   for(Int4 s=0; s < LenSeq(qE); s++) free(mtx[s]); free(mtx);
		} else SAP[II]=0;
	   } else {
		sap=gpsi.StealSAP();
		if(sap){
		  SAP[II]=MakeSelfAlignedGSAP(qE);
        	  SAP[II]->next=sap;
		} else SAP[II]=0;
	   }
	   II++;
   } // end loop over cma files
   if(chkfp) fclose(chkfp);

   FILE *fptr=0;
   NEW(CMA,num_cma_files+5,cma_typ);
   //************************* 4. process hits ************************
   // merge sap lists deleting weakest redundant hits. 
   //************************* 4a. check for hits ************************
   BooLean found_hits=FALSE,found_false=FALSE;
   for(II=Begin; II <= num_cma_files; II++){
	 if(IGNORE_THIS[II] && SAP[II]) found_false=TRUE; 
	 else if(SAP[II]) found_hits=TRUE; 
   }
   if(found_hits){
	//********************* 4b. Find the best matches and remove poor hits ********************
	sap_typ *saps;
	UInt4 *NUM_EACH_CMA=0;
	NEW(NUM_EACH_CMA,num_cma_files+3,UInt4);
	if(find_repeats){	// ******** finding less seqs for some reason...
	  FindBestOverlappingSAPs(SAP,num_cma_files,NSeqsSeqSet(Data),A);
	  RemoveRejectsSAP(SAP,num_cma_files,A);
	} else {
	  LabelBestSAPs(SAP,num_cma_files,NSeqsSeqSet(Data),use_patch_sap,A,NUM_EACH_CMA);
	}
   }
   //************************* 5. Discard empty hits ************************
   BooLean found_seqs=FALSE;
   Int4 s;
   for(II=Begin; II <= num_cma_files; II++){ 
	if(SAP[II] && !AtLeastOneLabeledSAPS(SAP[II]->next)){
		FreeGSeqAlignList(SAP[II]); SAP[II]=0;
	} else if(SAP[II]){
		found_seqs=TRUE;
	}
   }
   //************************* 6. Create cma files of best hits ************************
   FILE *cmafp=0;
   BooLean *Skip=0;
   if(OutputMAPS && found_seqs){
	NEW(Skip,num_cma_files+5,BooLean);
	sprintf(str,"%s_MAP",argv[1]);
   	cmafp=open_file(str,".cma","w");
   }
    // use this to output cma files for iterative cycling....
   for(II=Begin; II <= num_cma_files; II++){
       if(SAP[II]){
   	   // found_seqs=TRUE;
           fptr=tmpfile(); xSeqAlignToCMA(fptr,SAP[II],left_flank,right_flank,A);
           rewind(fptr); CMA[II]=ReadCMSA(fptr,A); fclose(fptr);
	   StrSeqID(str, 20, QueryE[II]); RenameCMSA(str,CMA[II]);
	   if(cmafp && II > 0){	// II > 1 eventually
		PutCMSA(cmafp,CMA[II]);
	   }
       } else { 
	    CMA[II]=0; 
	    if(Skip) Skip[II+1]=TRUE; 
	    if(II >= first_false_pos) num_false_pos--;
       }
	// index as II+1 because Template starts with consensus sequence as #1.
   }

#if 1
   if(cmafp){   
   	RenameCMSA("template",TemplateCMA);
   	SetLevelCMSA(num_false_pos,TemplateCMA);
   	PutSelectCMSA(cmafp,Skip,TemplateCMA);  free(Skip);
   	fclose(cmafp);
   }
#endif


   //************************* 7. output sequences ************************
   FILE *main_sq_fp=0;
   if(found_seqs){
	main_sq_fp=open_file(argv[2],"_aln.seq","w");
   }
   Int4 TotalTrue=0,TotalFalse=0;
   for(II=Begin; II <= num_cma_files; II++){ 
       	 if(SAP[II]){
	    if(II < first_false_pos) TotalTrue+= NumSeqsListGSAP(SAP[II]);
	    else TotalFalse += NumSeqsListGSAP(SAP[II]);
            if(printopt == 4) PutMultiGSeqAlign(stdout, SAP[II], width, A);
            else if(printopt == 2) PutGSeqAlignList(stdout, SAP[II]->next, width, A);
            else if(printopt == 1){
	        printf("===> File %d: %s.\n",II,NameCMSA(CMA[II]));
	        fprintf(stdout,"===> File %d: %s.\n",II,NameCMSA(CMA[II]));
		PutGSeqAlignList(stdout, SAP[II]->next, width, A);
	    }
   	    if(output_by_class){
   		FILE *sq_fp=0;
		sap=SAP[II]->next;
		// first sap is the consensus sequence, so need to use sap->next...
		if(out_all_class){
           	  sprintf(str,"%s_%d.seq",argv[2],II); sq_fp=open_file(str,"","w");
		  PutSeqsListGSAP(sq_fp, sap,A); fclose(sq_fp);
		}
#if 0
		if(use_weak_cutoff && EvalueGSAP(sap) > WeakEcutoff && !IGNORE_THIS[II]) {
           	  sprintf(str,"%s_%d.wsq",argv[2],II);
           	  sq_fp=open_file(str,"","w");
		  PutSeqsListGSAP(sq_fp,sap,A);
   		  fclose(sq_fp);
		  if(use_weak_cutoff){
		    sprintf(str,"%s_%d.weak",argv[2],II);
		    sq_fp=open_file(str,"","w");
		    // fprintf(stderr,"WeakEcut = %g; Evalue = %g\n",WeakEcutoff,EvalueGSAP(SAP[II]));
	            fprintf(sq_fp,"===> File %d: %s.\n",II,NameCMSA(CMA[II]));
		    PutGSeqAlignList(sq_fp, sap, width, A);
		    // PutGSeqAlign(stderr, sap, 60, A);
		    fclose(sq_fp);
		  }
		}
#else
		if(use_weak_cutoff && !IGNORE_THIS[II]) {
			sq_fp = 0;
			sap_typ tmp_sap;
			double MinEval=0.0;
			sap_typ head = sap;
			do { 
				tmp_sap= MinEvalSAP(head, &MinEval);
				if(MinEval > WeakEcutoff){
				   if(sq_fp == 0){ 
					sprintf(str,"%s_%d.weak",argv[2],II);
					sq_fp=open_file(str,"","w");
				   }
	            		   fprintf(sq_fp,"===> File %d: %s.\n",II,NameCMSA(CMA[II]));
				   PutOneSeqAlign(sq_fp, head, width, A);
				}
				head = tmp_sap;
			} while(head);
			if(sq_fp) fclose(sq_fp);
		}
#endif
	    }
	    if(main_sq_fp) PutSeqsListGSAP(main_sq_fp, SAP[II],A);
         }
   }
   if(main_sq_fp) fclose(main_sq_fp);

   //************************* 8. Use template to convert CMA files ************************
   //**************************** Key point in algorithm ***********************************
   cma_typ *OUT_CMA=0;
   if(test > 0) OUT_CMA=ConversionViaTemplateCMSA(TemplateCMA,CMA);
   else OUT_CMA=ConversionViaTemplateCMSA3(TemplateCMA,CMA);
   if(Begin == 0){ OUT_CMA[0]=CMA[0]; }
   //******************* CALL recursively for hierarchical alignments? **********************
   for(II=Begin; II <= num_cma_files; II++){	
	if(OUT_CMA[II]){
	   if(printopt==2){	// output regions in vsi format...(use -m= option)
		Int4 sq=2;
		PutVSIregionsCMSA(stdout,sq,colors,Start,End,NumRegions,OUT_CMA[II]);
	   }
	   RenameCMSA(NameCMSA(CMA[II]),OUT_CMA[II]); // Name output cma files as in template.
   	   if(output_by_class){
   		char str2[500];
		FILE *ofp;
		BooLean *skip=0; 
		// turn this on to omit concensus from the output alignment.
		// NEW(skip,NumSeqsCMSA(OUT_CMA[II])+3,BooLean); skip[1]=TRUE;
	        if(out_all_class) {
	           sprintf(str2,"%s_%d.cma",argv[2],II); 
     		   ofp = open_file(str2,"","w"); PutSelectCMSA(ofp,skip,OUT_CMA[II]);
		   fclose(ofp); 
		}
		sap=SAP[II]->next;	
#if 0	// this needs to be fixed...
		if(use_weak_cutoff && EvalueGSAP(sap) > WeakEcutoff && !IGNORE_THIS[II]){
	            sprintf(str2,"%s_%dw.cma",argv[2],II); 
     		    ofp = open_file(str2,"","w"); PutSelectCMSA(ofp,skip,OUT_CMA[II]);
		    fclose(ofp); 
		} 
#endif
		if(skip) free(skip);
	   }
	}
   } 

   //************************* 9. Summary of output files ************************
   printf("\nDatabase search results:\n");
   for(II=Begin; II <= num_cma_files; II++){
	if(found_false && II == first_false_pos){
		fprintf(stdout,
		  "--(false +)-------------------------------------\n");
	}
	if(OUT_CMA[II]){ 
	  RenameCMSA(NameCMSA(CMA[II]),OUT_CMA[II]);
	  if(IGNORE_THIS[II]){ printf("---------------------> File %d: %s (%d hits).\n",
                        II,NameCMSA(OUT_CMA[II]),NumSeqsCMSA(OUT_CMA[II])-1);
	  } else printf("=====================> File %d: %s (%d hits).\n",
                        II,NameCMSA(CMA[II]),NumSeqsCMSA(OUT_CMA[II])-1);
	} else {
	   StrSeqID(str, 20, QueryE[II]);
	   printf("--- File %d: %s (no hits).\n",
                        II,str);
	}
   }
   for(s=Begin-1,II=Begin; II <= num_cma_files; II++){	// compress array for output...
	if(!IGNORE_THIS[II] && OUT_CMA[II]!=0) {
		s++; CMA[s]=OUT_CMA[II]; 
		ReNameCMSA(NameCMSA(OUT_CMA[II]),CMA[s]);
	}
   }
   set_typ *set;
   NEW(set,s+3,set_typ);
   //************************* 10. output fasta sequence files ************************
   for(II=Begin; II <= s; II++){
	   set[II]=MakeSet(NumSeqsCMSA(CMA[II])+1);
	   FillSet(set[II]); DeleteSet(0,set[II]); DeleteSet(1,set[II]); 
   }
   //************************* 11. output Main alignment file ************************
   if(s > 0){
     fptr = open_file(argv[2],"_aln.cma","w");
     PutMergedCMSA(fptr,NameCMSA(TemplateCMA),s,set,CMA,NULL); 
     fclose(fptr);
   } else fprintf(stdout,"No hits found.\n");
#if 0
   Int4    *HitList=ListSet(UnionSetH);
   fptr = open_file(argv[2],"_aln.seq","w");
   // PutSet(fptr,UnionSetH);
   for(Int4 i=0; HitList[i] != -1; i++){
	PutSeq(fptr,SeqSetE(HitList[i],Data),A);
   } fclose(fptr);
#endif
   for(II=Begin; II <= s; II++){ if(set[II]) NilSet(set[II]); }
   if(Checkpoint) free(Checkpoint);
   if(TemplateCMA) TotalNilCMSA(TemplateCMA);	// merging ruins the seqset!!!!!
   for(II=0; II <= num_cma_files; II++){
	if(SetH[II]) NilSet(SetH[II]);
	if(SetE[II]) NilSet(SetE[II]);
	if(SetM[II]) NilSet(SetM[II]);
	if(OUT_CMA[II]) TotalNilCMSA(OUT_CMA[II]);	// merging ruins the seqset!!!!!
   } NilSeqSet(Data); NilAlpha(A);
   free(SetH); free(SetM); free(SetE);
   if(found_seqs){
      Int4 a = CardSet(UnionSetH);
      // Int4 f=CardInterSetINotJ(FalseSetH,UnionSetH);
      fprintf(stdout,"%d hits (true = %d; false = %d) in %d seqs.\n",
				TotalTrue+TotalFalse,TotalTrue,TotalFalse,a);
   } else fprintf(stdout,"no hits.\n");
   double runtime = difftime(time(NULL),time1);
   fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
   if(success) return 0;
   else return 1;
}

