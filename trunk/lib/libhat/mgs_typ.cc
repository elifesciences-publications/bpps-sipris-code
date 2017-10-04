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

#include "mgs_typ.h"

const char USAGE[]="usage: gapmap query_prefix dbsfile [options]\n\
      NOTE: Input multiple profile file = query.mpa (created using mkmaps)\n\
            Input Template alignment file = query.tpl \n\
      options:\n\
        -C=<int>  Cutoff score for eliminating overlaps (-P option)(default: 30)\n\
        -D        Don't SEG query sequence\n\
        -e=<real> E-value cutoff (default: 0.01/number of profiles)\n\
        -fract=<real> Fraction of residue identities for weak hits (default: 0.5).\n\
                    Used with the -o=<real> option.\n\
        -G        use global alignment (disallow deletions at ends)(not yet fully implemented)\n\
        -H=<int>  heapsize (default 500000)\n\
        -I=<int>:<int> - lengths of left & right flanking regions to be retained (default: 10:10)\n\
        -l        do low complexity masking of database sequences\n\
        -L        do low complexity and coiled coil masking of database sequences\n\
        -M        Output cma files for iterative recycling \n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -m=[<char>:<int>..<int>;]<char>:<int>..<int>.  generate output for a vsi file\n\
	          (e.g., -m=R:4..41,O:54..76,Y:81..104,G:112..137,B:141..156.)\n\
        -O        output sequence hits in separate files for each profile\n\
        -o=<real> output sequence hits in separate files only if E-value is above <real>\n\
        -P        Don't use PatchSAP() option to combine gapped HSPs into one\n\
        -prune=<char> Basis for input set pruning (default: 'T')\n\
                   'T': Template HSP trigger (slowest & most sensitive).\n\
                   'E': Template extensions (intermdiate).\n\
                   'H': Template hits (fastest & least sensitive).\n\
        -r        Don't detect repetitive domains\n\
        -sense=<real>  sensitivity parameter for heuristic search \n\
		   range: 0.01...1.0; lower values: > sensitivity (default: 0.05)\n\
	-size=<int>  set effective length of the database (use zero for the real size)\n\
        -T=<int>  blast word hit threshold (default 11)\n\
        -t<real>  trim alignments at ends where marg. prob. is > <real>\n\
        -tree     Use two-level search tree (names of non-terminal nodes start with T<int>_)\n\
        -trigger=<real>  Number of bits to trigger gapping = <real> (default: 22.0)\n\
        -test=<int>  test mode (default: 0)\n\
                   test=0:  don't allow deletions next to insertions\n\
                   test=1:  allow deletions next to insertions\n\
        -w=<int>  alignment width (default 60)\n\
        -weaklen=<int> Minimum number of matches required for weak hits (default: 0).\n\
                    Used with the -o=<real> option.\n\
        -X=<int>  X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
	-x        For mgapmap: compare _new.mma output files (instead of running mapgaps)\n\
\n";

void	mgs_typ::GetArg(int argc, char *argv[]) { GetArg(argc, argv,USAGE); }

void	mgs_typ::Init( )
{
	AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	// ************** files: ****************
	infile=0; database=0; Checkpoint=0;
	FDT=0;

	// ************** objects & arrays: ****************
	QueryE=0; SAP=0;
	TemplateCMA=0; ALN_CMA=0; SRCH_CMA=0;
	IGNORE=0;
   	found_false=FALSE;	// mgs class boolean...

	// ************** psi-blast parameters: ****************
	T=11;gap_open=11;gap_extend=1; x_parameter=25.0;
	gap_trigger_bits=0.0;
	left_flank=10;right_flank=10;
	printopt=0; width=60;
	hpsz=1000000;
	Ethresh=0.001;Ecutoff=0;
	seg_seq=TRUE;
	mask_dbs=' ';
	misalncut=2.0;
	dbslen=0;
	OwnData=FALSE;

	// ************** database parameters: ****************
	max_length_input_seqs=300000;

	// ************** MAPGAPS parameters: ****************
	use_patch_sap=TRUE;
	output_by_class=FALSE;out_all_class=FALSE;
	print_names_only=FALSE;find_repeats=TRUE;
	num_cma_files;overlap_cutoff=30;
	use_weak_cutoff= FALSE; 
	fractionMS=0.05;WeakEcutoff=0.0;
	prune_mode='T';
	num_false_pos=0; first_false_pos=0;
	OutputMAPS=FALSE;
   	found_seqs=FALSE;
	SeqsTrue=0,SeqsFalse=0;
	FalseData=0; FalseCheckPoint=0;
	fract_IDs=0.50;
	MinWeakAlnLen=0;
	Use2LevelTree=FALSE;
	DoGlobalAlign = FALSE;

	// ************** Test parameters ****************
	test=0;
	Begin=1;	// Begin with template (II=0) or with first profile (II=1).

	// ************** VSI parameters: ****************
	Start=0;End=0;
	NumRegions=0;

	//******************** tree ************************
	dft=0;
}

void	mgs_typ::GetArg(int argc, char *argv[],const char *usage)
{
	if(argc < 3) print_error(usage);
	//***************************** 1. Read input files *******************************
	a_type	A=AB;

	for(Int4 arg=0; arg < argc; arg++) { fprintf(stderr,"%s ",argv[arg]); }
	fprintf(stderr,"\n\n");
	// fprintf(stderr,"num_cma_files = %d\n",num_cma_files);
	//***************************** 2. Get input options *******************************
        for(Int4 arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(usage);
           switch(argv[arg][1]) {
	     case 'C': overlap_cutoff=IntOption(argv[arg],'C',1,100,usage); break;
	     case 'D': seg_seq=FALSE; break;
             case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,usage); break;
	     case 'f': 
		if(sscanf(argv[arg],"-fract=%lf",&fract_IDs) != 1)
                        print_error(usage); break;
	     case 'G': DoGlobalAlign = TRUE; break;
	     case 'H': hpsz=IntOption(argv[arg],'H',1,100000000,usage); break;
             case 'h': Ethresh=RealOption(argv[arg],'h',0.0,10000,usage); break;
	     case 'I': // putsegs=TRUE;
		if(sscanf(argv[arg],"-I=%d:%d",&left_flank,&right_flank) != 2)
                        print_error(usage); break;
             case 'l': if(argv[arg][2] != '='){ mask_dbs='x'; } break;
             case 'L': 
		if(argv[arg][2] != '='){ mask_dbs='b'; }  // mask out both low complex & coiled coils.
		break;
	     case 'm': 	// R:7..41,O:54..76,Y:81..104,G:114..137,B:141..156.
		{
		  if(argv[arg][2] != '='){
			printopt=IntOption(argv[arg],'m',0,6,usage); 
		  } else {
			NumRegions=ParseInputRegionsVSI(argv[arg]+3,
				&Start,&End,&colors,usage);
			printopt=2;
		  }
		} break;
	     case 'M': OutputMAPS=TRUE; argv[arg][1] = ' '; break;
             case 'O': {
		  out_all_class=TRUE; output_by_class= TRUE; argv[arg][1] = ' '; 
		} break;
             case 'o': {
             		WeakEcutoff=RealOption(argv[arg],'o',0.0,10000,usage);
			// fprintf(stderr,"WeakEcutoff = %g\n",WeakEcutoff);
			output_by_class= TRUE; use_weak_cutoff= TRUE; argv[arg][1] = ' '; 
		} break;
             case 'P': use_patch_sap=FALSE; argv[arg][1] = ' '; break;
	     case 'p': 
		if(sscanf(argv[arg],"-prune=%c",&prune_mode) != 1){
                        print_error(usage); 
		} else if(prune_mode != 'H' && prune_mode != 'E' && prune_mode != 'T'){
			print_error("prune_mode input error");
		} break;
             case 'r': find_repeats=FALSE; argv[arg][1] = ' '; break;
             case 's': 
             	if(sscanf(argv[arg],"-sense=%lf",&fractionMS) == 1){
			if(fractionMS < 0.01 || fractionMS > 1.0) print_error(usage);
			// else fprintf(stderr,"sensitivity = %f\n",fractionMS);
		} else if(sscanf(argv[arg],"-size=%ld",&dbslen) == 1){
			if(dbslen <= 0) print_error(usage);
			argv[arg][1] = ' ';
		} else print_error(usage);
			break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,usage); // argv[arg][1] = ' ';
		break;
	     case 't': 
		if(strcmp("-tree",argv[arg]) == 0){ 
			Use2LevelTree=TRUE; 
		} else if(sscanf(argv[arg],"-trigger=%lf",&gap_trigger_bits) == 1){
			if(gap_trigger_bits < 1.0 || gap_trigger_bits > 1000.0){
				print_error("gap_trigger out of permissible range");
			}
		} else if(sscanf(argv[arg],"-test=%d",&test) == 1){
		   if(test > 1 || test < 0) print_error(usage);
		} else misalncut=RealOption(argv[arg],'t',0.1,1.0,usage);
		break;
	     case 'w': 
		if(argv[arg][2] == 0){
			width=IntOption(argv[arg],'w',5,500000,usage);
		 } else if(sscanf(argv[arg],"-weaklen=%ld",&MinWeakAlnLen) != 1){
                        print_error(usage); 
		 } break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,usage); break;
             case 'Z': print_names_only=TRUE; break;
             case 'x':  break; // ignore this
             case ' ':  break; // ignore this
             default: print_error(usage);
           }
        }
	infile=argv[1];
	database=argv[2];
}

void	mgs_typ::SetUpSearch(ss_type in_data)
{
  Int4 II,j;
  ALN_CMA=0;
  //***************************** 2. Initialize search *******************************
  UnionSetH=0;FalseSetH=0; TrueSetH=0;
  sprintf(str,"%s.tpl",infile); 
  TemplateCMA=ReadCMSA2(str,AB);
  if(nBlksCMSA(TemplateCMA) != 1) print_error("Template file input error");
  if(!IsOkayTemplateCMA(TemplateCMA)) print_error("Template file formating error");
  // if(MinWeakAlnLen==0) MinWeakAlnLen=LengthCMSA(1,TemplateCMA);
  // ^ not what I want; set this to 0 by default.
  num_cma_files = NumSeqsCMSA(TemplateCMA) - 1;
  if(Ecutoff == 0){		// Adjustment is less conservative as profiles are not independent
#if 0	// AFN: 6/6/2017.
	if(num_cma_files < 10){			// i.e., 2..10
		Ecutoff = 0.001;
	} else if(num_cma_files < 100){		// i.e., 10...100
		Ecutoff = 0.0004;
	} else if(num_cma_files < 1000){	// i.e., 100...1000
		Ecutoff = 0.0001;
	} else {
		Ecutoff = 0.00001;
	} 
#else
	Ecutoff=0.01/(double) num_cma_files;
#endif
  }
  first_false_pos = num_cma_files + 1;
  FILE	  *fp;
  BooLean IncludeFalsePos=FALSE;
  sprintf(str,"%s.xup",infile);
  if((fp=fopen(str,"r")) != 0){
	fclose(fp); sprintf(str,"%s.xpq",infile);
	if((fp=fopen(str,"r")) != 0){ fclose(fp); IncludeFalsePos=TRUE; }
  }
  if(IncludeFalsePos){
  	sprintf(str,"%s.xup",infile); FalseCheckPoint = AllocString(str);
  	sprintf(str,"%s.xpq",infile); FalseData=MakeSeqSet(str,max_length_input_seqs,AB);
  	num_false_pos=NSeqsSeqSet(FalseData);
	num_cma_files += num_false_pos;
  } else {
	fprintf(stderr,"False positive files (%s.xpq, %s.xup) not provided.\n\n",
		infile,infile);
	num_false_pos=0;
  }
#if 1	// dft tree file?
  sprintf(str,"%s.dft",infile);
  if((fp=fopen(str,"r")) != 0){
	fclose(fp); 
	dft = new dft_typ(infile);
	if(!dft->IsValid()){ 
		delete dft;
		dft=0;
	}
  } else {
	fprintf(stderr,"Pruning tree file (%s.dft) not provide...\n",infile);
	fprintf(stderr,"   performing a single-level search instead.\n");
	dft=0;
  }
#endif
  NEW(QueryE,num_cma_files+3,e_type);
  NEW(SAP,num_cma_files+3,sap_typ);
  // for(II=0; II <= num_cma_files; II++){ QueryE[II]=0; SAP[II]=0; }
  QueryE[0]=TrueSeqCMSA(1,TemplateCMA); 
  for(j=1,II=1; II <= num_cma_files; II++){
     if(II >= first_false_pos){ QueryE[II] = SeqSetE(j,FalseData); j++; }
     else {    
	QueryE[II]=TrueSeqCMSA(II+1,TemplateCMA); 
	if(1 && print_names_only){
		fprintf(stdout,"=====================> File %d: ",II);
		PutSeqID(stdout,QueryE[II]);
		fprintf(stdout,"\n");
	}
     }
  }
  if(print_names_only){ fprintf(stdout,"\n"); exit(0); }
  if(in_data) { OwnData=FALSE; Data=in_data; }
  else { OwnData=TRUE; Data=MakeSeqSet(database,max_length_input_seqs,AB); }
  // ReNumberSeqSet(Data);	// this should be done automatically.
  // sprintf(str,"%s.map",infile); Checkpoint = AllocString(str);
  sprintf(str,"%s.mpa",infile); Checkpoint = AllocString(str);
  SeqsTrue=0,SeqsFalse=0;

  if(!IGNORE) NEW(IGNORE,num_cma_files+3,BooLean);
  for(Int4 n=first_false_pos; n <= num_cma_files; n++) IGNORE[n]=TRUE;
  // ReNumberSeqSet(Data);
  // Search database for hits using distinct queries for each cma file.
}

Int4	*mgs_typ::RtnTree( )
// creates a two level tree for testing...
{
  Int4    II;
  Int4	  *Level;
  char id_str[15];
  NEW(Level,num_cma_files+3,Int4);	// Level[0]=0;
  for(II=1; II <= num_cma_files; II++){
	cma_typ cma=0;
  	e_type  qE=CopySeq(QueryE[II]); 
	if(II < first_false_pos){
#if 0
		Int4	x;
		StrSeqID(id_str, 10, qE); 
		if(sscanf(id_str,"T%d_",&x) == 1){
   	   		// fprintf(stderr,"...shift point...\n");
			Level[II]=1;
			// Threshold=11;
		} else {
			Level[II]=2;
			// Threshold=18;
		}
#else
		Level[II]=1;
#endif
	} else if(II >= first_false_pos){
	        Level[II]=1;
	}
   } return Level;
}

BooLean	mgs_typ::Search(FILE *efp)
// Takes a checkpoint file and array of corresponding Query sequences as input
// returns an array of SAP[II], one for each profile.
// also returns UnionSetH, FalseSetH.
{
#if 1	// testing trees...
  if(dft){
	Int4  *Level=dft->RtnLevels( );
	Int4  N=dft->NumNodes( );
	if(N != num_cma_files + 1){
		fprintf(stderr,"num_cma_files=%d; num_false_pos=%d; Tree Nodes=%d\n",
			num_cma_files,num_false_pos,N);
		print_error("dft file input error");
	}
  	return TreeBasedSearch(Level);
  } else {
  	Int4	*Level=RtnTree( );
  	return TreeBasedSearch(Level);
  }
#endif
  BooLean success=FALSE;
  sap_typ sap;
  Int4    II,maxrounds;
  set_typ *SetM=0,*SetH=0,*SetE=0;
  set_typ SkipSet;
  e_type  qE;
  FILE *chkfp=open_file(Checkpoint,"","r");
//***************************** TESTING Heuristics ********************************
  NEW(SetM,num_cma_files+5,set_typ); NEW(SetH,num_cma_files+5,set_typ);
  NEW(SetE,num_cma_files+5,set_typ);
//***************************** TESTING Heuristics ********************************
  char id_str[15];
  set_typ ShiftPointSet=0;
  clock_t timeTot=clock();
  for(II=0; II <= num_cma_files; ){
	if(II == first_false_pos){
		fclose(chkfp); 
		chkfp=open_file(FalseCheckPoint,"","r");
	}
	clock_t timeS=clock();
  // start main loop over all *cma files....
	cma_typ cma=0;
	qE=CopySeq(QueryE[II]); 
	// PutSeq(stderr,qE,AB);
	EqSeqI(0,qE); maxrounds=1;
	gpsi_type gpsi(qE,Data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,chkfp,dbslen);
	if(gap_trigger_bits > 0.0) gpsi.SetGapTriggerInBits(gap_trigger_bits);	// default: 22.0;
	// else if(II==0) gpsi.SetGapTriggerInBits(15.0);  // more sensitive for first search;
//***************************** TESTING ********************************
	Int4 Threshold=T;
	if(II == 0){
	   // gpsi.CreateMarginalSet(); 
	} else if(ShiftPointSet && II < first_false_pos){  // when Use2LevelTree == TRUE.
		Int4	x;
		StrSeqID(id_str, 10, qE); 
		if(sscanf(id_str,"T%d_",&x) == 1){
   	   		// fprintf(stderr,"...shift point...\n");
			Threshold=T;
			gpsi.ProvideSkipSet(SkipSet);	// pop back up to root node...
			
		} else {
			// Threshold=18;
			Threshold=T;
			gpsi.ProvideSkipSet(ShiftPointSet);
		}
	} else if(II >= first_false_pos) {
	   if(II == first_false_pos){
		fprintf(stderr,
		  "--(excluded profiles)-------------------------------------------------\n");
	   }
	   gpsi.ProvideSkipSet(UnionSetH);	 // only look at sequence in UnionSetH
	} else {
	   gpsi.ProvideSkipSet(SkipSet); 
	}
	gpsi.CreateMarginalSet(fractionMS); 
	gpsi.CreateExtendSet(); 
	gpsi.CreateHitSet(); 
	gpsi.KeepQuite();
//***************************** TESTING ********************************
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,AB);
	Int4 **mtx=0;
	if(use_patch_sap) mtx=gpsi.CopyOFposMatrix();
	// fprintf(stderr,"II=%d; pass=%d\n",II,gpsi.ThisPass( ));
	   // if(gpsi.ThisPass( ) == 0) gpsi.SpeakUp(); 
	if(misalncut < 1.0) sap = gpsi.TrimmedSearch(Threshold,misalncut,mask_dbs);
	else sap = gpsi.Search(Threshold,mask_dbs);
	if(sap) success=TRUE; // else break;
	gpsi.ComputeMatrix(0);
//***************************** TESTING ********************************
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
		// CopySet(UnionSetH,SetH[II]);
	        FalseSetH = MakeSet(SetN(SkipSet));
		Int4 NumPruned=CardSet(SkipSet);
		if(NumPruned == 0){ fprintf(stderr,"no hits found\n"); exit(0); }
		fprintf(stderr,"   pruning input set from %d to %d sequences (%.1lf%%)\n",
			NSeqsSeqSet(Data),NumPruned,
				100.0*(double)NumPruned/(double)NSeqsSeqSet(Data));
	   	fprintf(stderr,"No.: %-10s     %-11s   %-11s   %-11s   (time)(Hits) level\n",
			"name","Hits","Extend","Trigger");
	   	fprintf(stderr,"%3d: ",II); 
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
	   	fprintf(stderr,"%3d: ",II); 
	   	StrSeqID(id_str, 10, qE); 
	   	fprintf(stderr,"%-10s     %-11d   %-11d   %-11d",id_str,
			CardSet(SetH[II]),CardSet(SetE[II]),CardSet(SetM[II]));
	   	if(II < first_false_pos){
	   		UnionSet(UnionSetH, SetH[II]);	// key track of true hits.
	   	} else {
	   		UnionSet(FalseSetH, SetH[II]);	// key track of false hits.
	   	}
   	   	fprintf(stderr," (%0.2f s)(%d)",
				(double) (clock()-timeS)/(double)CLOCKS_PER_SEC,
				CardSet(UnionSetH));
	   }
	   if(Use2LevelTree){
	       Int4	x;
	       StrSeqID(id_str, 10, qE); 
	       if(sscanf(id_str,"T%d_",&x) == 1){
   	   		fprintf(stderr,"pruning node");
			// ShiftPointSet=SetH[II];
			ShiftPointSet=SetE[II];
	       } 
	   }
	   fprintf(stderr,"\n");
//***************************** TESTING ********************************
	   sap_typ head;
	   if(use_patch_sap){ // fix SAP here...to patch together HSPs
		sap=gpsi.StealSAP();
		if(sap){
		   head=MakeSelfAlignedGSAP(qE);
		   head->next = sap;
		   SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data),qE,AB); 
		
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
   fprintf(stderr," total search time = %0.2f s.\n",
				(double) (clock()-timeTot)/(double)CLOCKS_PER_SEC);
   for(II=0; II <= num_cma_files; II++){
	if(SetH && SetH[II]) NilSet(SetH[II]);
	if(SetE && SetE[II]) NilSet(SetE[II]);
	if(SetM && SetM[II]) NilSet(SetM[II]);
   }
   if(SetH){ free(SetH); SetH=0; }
   if(SetM){ free(SetM); SetM=0; }
   if(SetE){ free(SetE); SetE=0; }
   if(chkfp) fclose(chkfp);
   return success;
}

BooLean	mgs_typ::TreeBasedSearch(Int4 *Level)
// Takes a checkpoint file and array of corresponding Query sequences as input
// returns an array of SAP[II], one for each profile.
// also returns UnionSetH, FalseSetH.
{
  BooLean success=FALSE;
  sap_typ sap;
  Int4    II,maxrounds;
  set_typ *SetM=0,*SetH=0,*SetE=0;
  set_typ SkipSet;
  e_type  qE;
  FILE *chkfp=open_file(Checkpoint,"","r");

  // New tree-based items:
  // tree is a depth-first-traversal integer-array representation of the profile tree.
  stk_typ stk(num_cma_files + 3);
  stk.Push(0);	// push the template node onto the stack
  Int4	  Parent=0,MaxLevel=6;
  for(II=1; II <= num_cma_files; II++){
	MaxLevel=MAXIMUM(Int4,MaxLevel,Level[II]);
  }

//***************************** TESTING Heuristics ********************************
  NEW(SetM,num_cma_files+5,set_typ); NEW(SetH,num_cma_files+5,set_typ);
  NEW(SetE,num_cma_files+5,set_typ);
//***************************** TESTING Heuristics ********************************
  char id_str[15];
  set_typ ShiftPointSet=0;
  clock_t timeTot=clock();
  for(II=0; II <= num_cma_files; ){
	if(II == first_false_pos){
		fclose(chkfp); 
		chkfp=open_file(FalseCheckPoint,"","r");
	}
	clock_t timeS=clock();
  // start main loop over all *cma files....
	cma_typ cma=0;
	qE=CopySeq(QueryE[II]); 
	// PutSeq(stderr,qE,AB);
	EqSeqI(0,qE); maxrounds=1;
	gpsi_type gpsi(qE,Data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,chkfp,dbslen);
	if(gap_trigger_bits > 0.0) gpsi.SetGapTriggerInBits(gap_trigger_bits);	// default: 22.0;
//***************************** TESTING ********************************
	Int4 Threshold=T;
	if(II == 0){
	   // gpsi.CreateMarginalSet(); 
	} else if(II < first_false_pos){  // when Use2LevelTree == TRUE.
		switch (Level[II]){
		  case 0: print_error("mgs_typ::TreeBasedSearch() - this should not happen"); break;
		  case 1: Threshold=11; 
			ShiftPointSet=SetM[Parent];	// then shift set...
			break;
		  case 2: Threshold=13; 
			ShiftPointSet=SetM[Parent];	// Parent is at level 1...
			break;
		  case 3: Threshold=15; 
			ShiftPointSet=SetE[Parent];	// Parent is at level 2...
			break;
		  case 4: Threshold=17;
			ShiftPointSet=SetE[Parent];	// then shift set...
			break;
		  default: Threshold=18;
			ShiftPointSet=SetE[Parent];	// then shift set...
			break;
		}
		gpsi.ProvideSkipSet(ShiftPointSet);
	} else if(II >= first_false_pos) {
	   if(II == first_false_pos){
		fprintf(stderr,
		  "--(excluded profiles)-------------------------------------------------\n");
	   }
	   gpsi.ProvideSkipSet(UnionSetH);	 // only look at sequence in UnionSetH
	} else {
	   gpsi.ProvideSkipSet(SkipSet); 
	}
	gpsi.CreateMarginalSet(fractionMS); 
	gpsi.CreateExtendSet(); 
	gpsi.CreateHitSet(); 
	gpsi.KeepQuite();
//***************************** TESTING ********************************
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,AB);
	Int4 **mtx=0;
	if(use_patch_sap) mtx=gpsi.CopyOFposMatrix();
	// fprintf(stderr,"II=%d; pass=%d\n",II,gpsi.ThisPass( ));
	   // if(gpsi.ThisPass( ) == 0) gpsi.SpeakUp(); 
	if(misalncut < 1.0) sap = gpsi.TrimmedSearch(Threshold,misalncut,mask_dbs);
	else sap = gpsi.Search(Threshold,mask_dbs);
	if(sap) success=TRUE; // else break;
	gpsi.ComputeMatrix(0);
//***************************** TESTING ********************************
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
		// CopySet(UnionSetH,SetH[II]);
	        FalseSetH = MakeSet(SetN(SkipSet));
		Int4 NumPruned=CardSet(SkipSet);
		if(NumPruned == 0){ fprintf(stderr,"no hits found\n"); exit(0); }
		fprintf(stderr,"   pruning input set from %d to %d sequences (%.1lf%%)\n",
			NSeqsSeqSet(Data),NumPruned,
				100.0*(double)NumPruned/(double)NSeqsSeqSet(Data));
	   	fprintf(stderr,
		       "ID:  %-10s     %-11s   %-11s   %-11s   (time)(Hits) Level (parent)\n",
			"name","Hits","Extend","Trigger");
	   	fprintf(stderr,"%3d: ",II); 
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
	   	fprintf(stderr,"%3d: ",II); 
	   	StrSeqID(id_str, 10, qE); 
	   	fprintf(stderr,"%-10s     %-11d   %-11d   %-11d",id_str,
			CardSet(SetH[II]),CardSet(SetE[II]),CardSet(SetM[II]));
	   	if(II < first_false_pos){
	   		UnionSet(UnionSetH, SetH[II]);	// key track of true hits.
	   	} else {
	   		UnionSet(FalseSetH, SetH[II]);	// key track of false hits.
	   	}
   	   	fprintf(stderr," (%0.2f s)(%d)",
				(double) (clock()-timeS)/(double)CLOCKS_PER_SEC,
				CardSet(UnionSetH));
	   }
	   { Int4 i;
	    for(i=1; i <= Level[II]; i++) fprintf(stderr," ");
	    fprintf(stderr," %d",Level[II]);
	    for(; i <= MaxLevel; i++) fprintf(stderr," ");
           } fprintf(stderr," (%d)\n",Parent);
	   if(II < first_false_pos && Level){ //************* Use tree here ***************
	       Int4	x;
	       if(II > 0 && II < num_cma_files){
	         assert(1 <= Level[II+1] && Level[II+1] <= Level[II]+ 1); //  for II >= 1.
	         // (adds an ith node as a child along the lineage leading to the i -1 node) 
	         x = Level[II] - Level[II+1];
	         if(x == -1){			// II is a parent of II + 1; 
		    Parent = II;
		    stk.Push(II);
	         } else if(x == 0){		// II and II+1 are siblings.
			// no need to do anything.
	         } else if(x > 0){		// Pop back x generations to common ancestoral node.
		    while(x >= 0){ Parent=stk.Pop(); x--; }
		    stk.Push(Parent);
	         } else {
			print_error("This should not happen");
	         }
	       }
	   } // stk.Put(stderr);	// check stack...
//***************************** TESTING ********************************
	   sap_typ head;
	   if(use_patch_sap){ // fix SAP here...to patch together HSPs
		sap=gpsi.StealSAP();
		if(sap){
		   head=MakeSelfAlignedGSAP(qE);
		   head->next = sap;
		   SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data),qE,AB); 
		
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
   fprintf(stderr," total search time = %0.2f s.\n",
				(double) (clock()-timeTot)/(double)CLOCKS_PER_SEC);
   for(II=0; II <= num_cma_files; II++){
	if(SetH && SetH[II]) NilSet(SetH[II]);
	if(SetE && SetE[II]) NilSet(SetE[II]);
	if(SetM && SetM[II]) NilSet(SetM[II]);
   }
   if(SetH){ free(SetH); SetH=0; }
   if(SetM){ free(SetM); SetM=0; }
   if(SetE){ free(SetE); SetE=0; }
   if(chkfp) fclose(chkfp);
   return success;
}

#if 0	// use this to make the alignment global...Moved to my_ncbi.cc

void	MkGlobalSeqAlign(sap_typ sap,a_type AB)
// make the Subject globally aligned with query...
// WARNING: This assumes that sap is for a mapgaps alignment; will cause major problems for MCMC alignments.
{
	Int4		r,i,j,s,e,q,qs,ss,len,length,index,total,m,start,end;
	Int4		qs1,qs2,ss1,ss2,len1,len2;
	Int4		lenQ,lenS;
	dsp_typ		dsp = sap->segs; 
	e_type		qE=dsp->qE,sE=dsp->sE;
	
	// if(dsp->numseg < 2) return;

	// 1. first check last segment (easiest to fix).
	i=dsp->numseg - 1;	// last segment.
	qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1]; len = dsp->lens[i];
	// qs = 0...lenQ-1; ss = 0...lenS-1
	lenQ=LenSeq(qE)-qs; // length of full query seq.
	lenS=LenSeq(sE)-ss; // length of full query seq.
	if(lenQ > len && lenS > len)
	{
		// fprintf(stderr,"FOUND DELETION ON END\n");
		if(lenS >= lenQ) dsp->lens[i] = lenQ;
		else dsp->lens[i] = lenS;
		// PutGSeqAlignList(stderr, sap, 60, AB);
	}

	// 2. next check first segment.
	i=0;	// first segment.
	qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1]; len = dsp->lens[i];
	if(qs > 0 && ss > 0){
		// fprintf(stderr,"FOUND DELETION AT BEGIN\n");
		if(ss >= qs){	// e.g., if ss = 12, qs = 5 then lengthen by qs = 5.
			dsp->starts[2*i] = 0;
			dsp->starts[2*i+1] = ss - qs;	// e.g., ss = 12 - 5 = 7.
			dsp->lens[i] = len + qs;	// e.g., len = len + 5.
		} else {	// e.g., if ss = 5, qs = 12 then lengthen by ss = 5
			dsp->starts[2*i] = qs - ss;	// e.g., qs = 12 - 5 = 7.
			dsp->starts[2*i+1] = 0;  
			dsp->lens[i] = len + ss;
		}
		// PutGSeqAlignList(stderr, sap, 60, AB);
	}
#if 0
	i=0; qs1 = dsp->starts[2*i]; ss1 = dsp->starts[2*i+1]; len1 = dsp->lens[i];
	i=1; qs2 = dsp->starts[2*i]; ss2 = dsp->starts[2*i+1]; len2 = dsp->lens[i];
	if(ss1 == -1 && qs1 != -1 && ss2 != -1 && qs2 != -1){
fprintf(stderr,"FOUND DELETION ON END\n");
		// if ss is deleted at end and matches at next segment.
		// then fix this as much as possible.
	    start = ss2+1;  // start of subject seq. in alignment
	    if(len1 <= start){  
		// if the first segment can be filled with the subject seq.
		// then collapse first two segments into one.
		i=0; qs1 = dsp->starts[2*i]; ss1 = dsp->starts[2*i+1]; len1 = dsp->lens[2*i];
		j=1; dsp->starts[2*j]=dsp->starts[2*i];  // query start = beginning of alignment.
		     dsp->starts[2*j+1]=start-len1;    // back up subject start by len1.
		     dsp->lens[j] = len1+len2;
		
		// delete first segment.
        	for(i=0,j=1; j < dsp->numseg; i++,j++){
		  dsp->starts[2*i] = dsp->starts[2*j] ;
		  dsp->starts[2*i+1] = dsp->starts[2*j+1] ;
		  dsp->lens[i] = dsp->lens[j];
		} dsp->numseg--;
	    } else if(start > 1){
	    }
	} 	// else do nothing...
#endif
}

void	MakeGlobalSeqAlign(sap_typ head,a_type AB)
{
	Int4 i=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
		// fprintf(stderr,"Globalizing align %d\n",++i);
		MkGlobalSeqAlign(gsap,AB);
	}
}

#endif

void	mgs_typ::ProcessHits( )
{
   FILE *fptr=0;
   Int4    II;
   NEW(SRCH_CMA,num_cma_files+5,cma_typ);
   //************************* 4. process hits ************************
   // merge sap lists deleting weakest redundant hits. 
   //************************* 4a. check for hits ************************
   BooLean found_hits=FALSE; 
   found_false=FALSE;	// mgs class boolean...
   for(II=Begin; II <= num_cma_files; II++){
	 if(IGNORE[II] && SAP[II]) found_false=TRUE; 
	 else if(SAP[II]) found_hits=TRUE; 
   }
   if(found_hits){
	//********************* 4b. Find the best matches and remove poor hits ********************
	UInt4 *NUM_EACH_CMA=0;
	NEW(NUM_EACH_CMA,num_cma_files+3,UInt4);
	if(find_repeats){	// ******** finding less seqs for some reason...
	  // also Labels Best SAPs.
	  FindBestOverlappingSAPs(SAP,num_cma_files,NSeqsSeqSet(Data),AB);
#if 0	// Debug...
   	  for(II=Begin; II <= num_cma_files; II++){ 
		if(SAP[II] && II == 13) PutGSeqAlignList(stderr, SAP[II], 60, AB);
		if(SAP[II] && II == 14) PutGSeqAlignList(stderr, SAP[II], 60, AB);
	  }
#endif
	  RemoveRejectsSAP(SAP,num_cma_files,AB); // not called from FindBestOverlappingSAPs();
	} else {
	  LabelBestSAPs(SAP,num_cma_files,NSeqsSeqSet(Data),use_patch_sap,AB,NUM_EACH_CMA);
#if 0	// Debug...
   	  for(II=Begin; II <= num_cma_files; II++){ 
		if(SAP[II] && II == 13) PutGSeqAlignList(stderr, SAP[II], 60, AB);
		if(SAP[II] && II == 14) PutGSeqAlignList(stderr, SAP[II], 60, AB);
	  }
#endif
	  // RemoveRejectsSAP(SAP,num_cma_files,AB); // called from LabelBestSAPs();
	}
   }
   //************************* 5. Discard empty hits ************************
   found_seqs=FALSE;
   // Int4 s;
   for(II=Begin; II <= num_cma_files; II++){ 
#if 0	// Debug...
		if(SAP[II] && II == 13) PutGSeqAlignList(stderr, SAP[II], 60, AB);
		if(SAP[II] && II == 14) PutGSeqAlignList(stderr, SAP[II], 60, AB);
#endif
	if(SAP[II] && !AtLeastOneLabeledSAPS(SAP[II]->next)){
		FreeGSeqAlignList(SAP[II]); SAP[II]=0;
	} else if(SAP[II]){
		found_seqs=TRUE;
	}
   }
   //************************* 6. Create cma files of best hits ************************
   for(II=Begin; II <= num_cma_files; II++){
       if(SAP[II]){
   	   // found_seqs=TRUE;
#if 1
	   if(DoGlobalAlign){
		MakeGlobalSeqAlign(SAP[II],AB);
	   }
#endif
           fptr=tmpfile(); xSeqAlignToCMA(fptr,SAP[II],left_flank,right_flank,AB);
           rewind(fptr); SRCH_CMA[II]=ReadCMSA(fptr,AB); fclose(fptr);
	   StrSeqID(str, 20, QueryE[II]); RenameCMSA(str,SRCH_CMA[II]);
       } else { 
	    SRCH_CMA[II]=0; 
	    if(II >= first_false_pos) num_false_pos--;
       }
	// index as II+1 because Template starts with consensus sequence as #1.
   }
}

void	mgs_typ::ConvertViaTemplate( )
// create converted cma files using the template. 
{
   ALN_CMA=0;
   if(test > 0) ALN_CMA=ConversionViaTemplateCMSA(TemplateCMA,SRCH_CMA);
   else ALN_CMA=ConversionViaTemplateCMSA3(TemplateCMA,SRCH_CMA);
   if(Begin == 0){ ALN_CMA[0]=SRCH_CMA[0]; }
}

void	mgs_typ::PutMAPInputFile( )
// use this to output cma files for iterative cycling....
{
   Int4    II;
   FILE *cmafp=0;
   BooLean *Skip=0;
   if(OutputMAPS && found_seqs){
     NEW(Skip,num_cma_files+5,BooLean);
     sprintf(str,"%s_MAP",infile);
     cmafp=open_file(str,".cma","w");
     for(II=Begin; II <= num_cma_files; II++){
       if(SRCH_CMA[II]){ // II > 1 eventually
	   if(II > 0){ PutCMSA(cmafp,SRCH_CMA[II]); }
       } else {
	   if(Skip) Skip[II+1]=TRUE; 
       }
       // index as II+1 because Template starts with consensus sequence as #1.
     }
     RenameCMSA("template",TemplateCMA);
     PutSelectCMSA(cmafp,Skip,TemplateCMA);  free(Skip);
     fclose(cmafp);
   }
}

void	mgs_typ::PutSeqs(BooLean putseq)
{
   //************************* 7. output sequences ************************
   Int4    II,sq;
   FILE *tsq_fp=0,*fsq_fp=0;
   FILE *wk_fp=0,*wsq_fp=0;
   e_type sE;
   sap_typ sap;
   if(!found_seqs) return;
   HSPsTrue=HSPsFalse=0;
   TrueSetH=0; FalseSetH=0;
   for(II=Begin; II <= num_cma_files; II++){ 
       	 if(SAP[II]){
       	    if(!SAP[II]->next) continue;
	    sap_typ head=SAP[II]->next;
	    // first sap is the consensus sequence, so need to use sap->next...
	    // MakeSelfAlignedGSAP(qE) used above to add this..
	    if(II < first_false_pos){
	       	if(tsq_fp==0){
   		  TrueSetH = MakeSet(SetN(UnionSetH));
		  if(putseq) tsq_fp=open_file(database,"_aln.seq","w");
		  else tsq_fp=stderr; 
	       	}
		for(sap=head; sap!=NULL; sap=sap->next){
		  sE = SubjectSeqGSAP(sap); sq = SeqI(sE);
		  if(!MemberSet(sq,TrueSetH)){
			if(putseq) PutSeq(tsq_fp,sE,AB); 
			AddSet(sq,TrueSetH); 
		  }
		}
	    }
	    if(II >= first_false_pos){
	    	if(fsq_fp==0){
   		  FalseSetH = MakeSet(SetN(UnionSetH));
		  fsq_fp=open_file(database,"_aln.fsq","w");
		}
		for(sap=head; sap!=NULL; sap=sap->next){
		   sE = SubjectSeqGSAP(sap); sq = SeqI(sE);
		   if(!MemberSet(sq,FalseSetH)){ PutSeq(fsq_fp,sE,AB); AddSet(sq,FalseSetH); }
		}
	    }
	    if(II < first_false_pos) HSPsTrue+= NumHSPsListGSAP(head);
	    else HSPsFalse += NumHSPsListGSAP(head);
            if(printopt == 4) PutMultiGSeqAlign(stderr, SAP[II], width, AB);
            else if(printopt == 2) PutGSeqAlignList(stdout, head, width, AB);
            else if(printopt == 1){
	        printf("===> File %d: %s.\n",II,NameCMSA(SRCH_CMA[II]));
	        fprintf(stderr,"===> File %d: %s.\n",II,NameCMSA(SRCH_CMA[II]));
		PutGSeqAlignList(stderr, head, width, AB);
	    }
   	    if(output_by_class){
		if(out_all_class){
		  FILE *sq_fp=0;
           	  sprintf(str,"%s_%d.seq",database,II); sq_fp=open_file(str,"","w");
		  PutSeqsListGSAP(sq_fp, head,AB); fclose(sq_fp);
		}
		if(use_weak_cutoff && !IGNORE[II]) {
			double MinEval=0.0;
			sap_typ tmp_sap;
			for(sap = head; sap; sap= tmp_sap){ 
				tmp_sap= MinEvalSAP(sap, &MinEval);
				if(MinEval > WeakEcutoff){
				   Int4 alnlen, Ident,Inserts,minlen;
				   minlen = IdentitiesGSAP(&alnlen, &Ident,&Inserts,sap);
				   double fraction=(double)Ident/(double)alnlen;
				   if(fraction <= fract_IDs && alnlen >= MinWeakAlnLen){
				     if(wk_fp == 0){ 
           	  			// sprintf(str,"%s_%d.wsq",database,II);
           	  			sprintf(str,"%s.wsq",database);
           	  			wsq_fp=open_file(str,"","w");
					// sprintf(str,"%s_%d.weak",database,II);
					sprintf(str,"%s.weak",database);
					wk_fp=open_file(str,"","w");
				     }
	            		     fprintf(wk_fp,"===> File %d: %s.\n",
							II,NameCMSA(SRCH_CMA[II]));
				     PutOneSeqAlign(wk_fp, sap, width, AB);
				     e_type sE = SubjectSeqGSAP(sap);
				     PutSeq(wsq_fp,sE,AB);
				   }
				}
			}
		}
	    }
         }
   }
   if(wk_fp) fclose(wk_fp);
   if(wsq_fp) fclose(wsq_fp);
   if(tsq_fp  && putseq) fclose(tsq_fp);
   if(fsq_fp) fclose(fsq_fp);
}

void	mgs_typ::PutAln( )
{
   Int4    II;
   sap_typ sap;
   //************************* 8. Use template to convert CMA files ************************
   //******************* CALL recursively for hierarchical alignments? **********************
   // only looks at the CMA files that are in template...
   for(II=Begin; II <= num_cma_files; II++){	
	if(SRCH_CMA[II]){
	   if(printopt==2){	// output regions in vsi format...(use -m= option)
		Int4 sq=2;
		PutVSIregionsCMSA(stdout,sq,colors,Start,End,NumRegions,ALN_CMA[II]);
	   }
	   if(II < first_false_pos){
		char *tmp_str=AllocString(NameCMSA(SRCH_CMA[II]));
	   	RenameCMSA(tmp_str,ALN_CMA[II]); free(tmp_str);
	   } // Name output cma files as in template.
   	   if(output_by_class){
   		char str2[500];
		FILE *ofp;
		BooLean *skip=0; 
		// turn this on to omit concensus from the output alignment.
		// NEW(skip,NumSeqsCMSA(ALN_CMA[II])+3,BooLean); skip[1]=TRUE;
	        if(out_all_class) {
	           sprintf(str2,"%s_%d.cma",database,II); 
		   if(II < first_false_pos){
     		      ofp = open_file(str2,"","w"); PutSelectCMSA(ofp,skip,ALN_CMA[II]);
		   } else {
     		      ofp = open_file(str2,"","w"); PutSelectCMSA(ofp,skip,SRCH_CMA[II]);
		   }
		   fclose(ofp); 
		}
		sap=SAP[II]->next;	
#if 0	// this needs to be fixed...
		if(use_weak_cutoff && EvalueGSAP(sap) > WeakEcutoff && !IGNORE[II]){
	            sprintf(str2,"%s_%dw.cma",database,II); 
     		    ofp = open_file(str2,"","w"); PutSelectCMSA(ofp,skip,ALN_CMA[II]);
		    fclose(ofp); 
		} 
#endif
		if(skip) free(skip);
	   }
	}
   } 
}


void	mgs_typ::PutSummary( )
//************************* 9. Summary of output files ************************
{
   Int4    II;
   printf("\nDatabase search results:\n");
   for(II=Begin; II <= num_cma_files; II++){
	if(found_false && II == first_false_pos){
		fprintf(stdout,
		  "--(excluded profiles)-------------------------------------\n");
	}
	if(SRCH_CMA[II]){ 
	  if(IGNORE[II]){ printf("---------------------> File %d: %s (%d hits).\n",
                        II,NameCMSA(SRCH_CMA[II]),NumSeqsCMSA(SRCH_CMA[II])-1);
	  } else {
	  	char *tmp=AllocString(NameCMSA(SRCH_CMA[II]));
	  	RenameCMSA(tmp,ALN_CMA[II]); free(tmp);
		printf("=====================> File %d: %s (%d hits).\n",
                        II,NameCMSA(ALN_CMA[II]),NumSeqsCMSA(ALN_CMA[II])-1);
	  }
	} else {
	   StrSeqID(str, 20, QueryE[II]);
	   printf("--- File %d: %s (no hits).\n",
                        II,str);
	}
   }
#if 0	// Output file names only for mcBPPS... 4/01/10 afn.
   for(II=Begin; II <= num_cma_files; II++){
	if(II < first_false_pos){
	  if(TRUE || SRCH_CMA[II]){ 	// found some hits...
	    StrSeqID(str, 50, QueryE[II]);
	    if(str[0] == 'T' && isdigit(str[1])) printf("\n");
	     printf("%s,",str);
	  }
	} 
   } printf("\n");
#endif
   if(found_seqs){
      Int4 a = CardSet(UnionSetH);
      if(TrueSetH) SeqsTrue = CardSet(TrueSetH);
      else SeqsTrue = 0;
      if(FalseSetH) SeqsFalse = CardSet(FalseSetH);
      else SeqsFalse = 0;
      // Int4 f=CardInterSetINotJ(FalseSetH,UnionSetH);
#if 0
      fprintf(stderr,"%d domains in %d seqs (%d relevant in %d sq; %d irrelevant in %d sq).\n",
		HSPsTrue+HSPsFalse,a,HSPsTrue,SeqsTrue,HSPsFalse,SeqsFalse);
#endif
      fprintf(stdout,"%d domains in %d seqs (%d relevant in %d sq; %d irrelevant in %d sq).\n",
		HSPsTrue+HSPsFalse,a,HSPsTrue,SeqsTrue,HSPsFalse,SeqsFalse);
   } else fprintf(stderr,"no hits.\n");
}

void	mgs_typ::PutMainAln(FILE *ofp)
{
   Int4	s,II;
   FILE	*fptr=0;
   for(s=Begin-1,II=Begin; II <= num_cma_files; II++){	// compress array for output...
	if(!IGNORE[II] && ALN_CMA[II]!=0) {
		s++; SRCH_CMA[s]=ALN_CMA[II]; 
#if 0	// DEBUG...
		char	str[100]; sprintf(str,"_aln%d.cma",s);
		fptr = open_file(database,str,"w");
		PutCMSA(fptr,SRCH_CMA[s]); fclose(fptr);
#endif
		char *tmp=AllocString(NameCMSA(ALN_CMA[II]));
		ReNameCMSA(tmp,SRCH_CMA[s]); free(tmp);
	}
   }
   set_typ *set;
   NEW(set,s+3,set_typ);
   //************************* 10. output fasta sequence files ************************
   for(II=Begin; II <= s; II++){
	   set[II]=MakeSet(NumSeqsCMSA(SRCH_CMA[II])+1);
	   FillSet(set[II]); DeleteSet(0,set[II]); DeleteSet(1,set[II]); 
   }
   //************************* 11. output Main alignment file ************************
   if(s > 0){
     if(ofp) PutMergedCMSA(ofp,database,s,set,SRCH_CMA,NULL);
     else {
       fptr = open_file(database,"_aln.cma","w");
       // PutMergedCMSA(fptr,NameCMSA(TemplateCMA),s,set,SRCH_CMA,NULL); 
       PutMergedCMSA(fptr,database,s,set,SRCH_CMA,NULL); 
       fclose(fptr);
     }
   } else fprintf(stderr,"No hits found.\n");
   for(II=Begin; II <= s; II++){ if(set[II]) NilSet(set[II]); }
}

void	mgs_typ::CleanUpSearch()
{
   if(Checkpoint) { free(Checkpoint); Checkpoint=0; }
   // merging ruins the seqset!!!!!
   for(Int4 II=0; II <= first_false_pos; II++){
	if(ALN_CMA && ALN_CMA[II]) TotalNilCMSA(ALN_CMA[II]);	// merging ruins the seqset!!!!!
   } 
   for(Int4 II=0; II <= num_cma_files; II++){
	// if(SRCH_CMA && SRCH_CMA[II]) TotalNilCMSA(SRCH_CMA[II]);
	// merging ruins the seqset!!!!!
   }
   if(ALN_CMA){ free(ALN_CMA); ALN_CMA=0; }
   if(SRCH_CMA){ free(SRCH_CMA); SRCH_CMA=0; }
}

void	mgs_typ::Free( )
{
   CleanUpSearch();
   if(TemplateCMA){ TotalNilCMSA(TemplateCMA); TemplateCMA=0; }
   if(SRCH_CMA){
     for(Int4 II=0; II <= num_cma_files; II++){
   	// if(SRCH_CMA[II]){ NilCMSA(SRCH_CMA[II]); }
	// This appears to be shared by other files...
     } free(SRCH_CMA); SRCH_CMA=0; 
   }
   if(QueryE) free(QueryE);
   if(SAP) free(SAP);
   if(OwnData) NilSeqSet(Data);
   if(FDT) delete FDT;
   NilAlpha(AB);
}

//************************** Compare Routines ***************************8
void	mgs_typ::StoreSummary(Int4 R)
//************************* 9. Summary of output ************************
{
   assert(FDT != 0);
   Int4    II;
   for(II=Begin; II <= num_cma_files; II++){
	if(found_false && II == first_false_pos){  // (excluded profiles)

	}
	if(SRCH_CMA[II]){ 
	  if(!IGNORE[II]){ 
	  	RenameCMSA(NameCMSA(SRCH_CMA[II]),ALN_CMA[II]);
		FDT->AssignHits(R,II+1,NumSeqsCMSA(ALN_CMA[II])-1);
	  }
	}
   }
   if(found_seqs){
      Int4 a = CardSet(UnionSetH);
      if(TrueSetH) SeqsTrue = CardSet(TrueSetH); else SeqsTrue = 0;
      if(FalseSetH) SeqsFalse = CardSet(FalseSetH); else SeqsFalse = 0;
      // FDT->AssignHits(R,0,HSPsTrue);
      fprintf(stderr,"%3d: %d domains in %d seqs (%d relevant in %d sq; %d irrelevant in %d sq).\n",
		R,HSPsTrue+HSPsFalse,a,HSPsTrue,SeqsTrue,HSPsFalse,SeqsFalse);
   } 
}

BooLean	mgs_typ::Compare(Int4 Ni,char **iR,char **iName,char *TypeOfSetI,
					Int4 No, char **oR,char **oName,char *TypeOfSetO)
{
	Int4 Number,f,r;
	Begin=0;
  	FDT= new fdt_typ(Ni,iR,iName,TypeOfSetI,No, oR,oName,TypeOfSetO);
	FILE *fp=open_file(database,"_new.mma","r");
	cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	if(Number > No){
		fprintf(stderr,"Number of cmafiles (%d) inconsistent with Hpt rows (%d)\n",Number,No);
		exit(1);
	} else if(Number < No){	// Make sure all Number are in Hpt & in the right order.
	     for(r=f=1; f <= Number; r++,f++){
		if(r > No) {
			fprintf(stderr,"cmafile names (%d) inconsistent with Hpt (%d)\n",Number,No);
			exit(1);
		}
		if(strcmp(NameCMSA(IN_CMA[f]),oName[r]) !=0){
			fprintf(stderr,"%3d: NameRow=%s missing\n",r,oName[r]);
			f--; 
		} else fprintf(stderr,"%3d: NameOut=%s; NameRow=%s\n",r,NameCMSA(IN_CMA[f]),oName[r]);
	     }
 	}
	for(r=f=1; f <= Number; r++,f++){
	   ss_type data=TrueDataCMSA(IN_CMA[f]);
	   // fprintf(stderr,"NameOut=%s; NameRow=%s\n",NameSeqSet(data),oName[f]);
	   if(strcmp(NameCMSA(IN_CMA[f]),oName[r]) !=0){ f--; continue; }
	   assert(NSeqsSeqSet(data) > 0);
	   // fprintf(stderr,"DEBUG: seqs = %d\n",NSeqsSeqSet(data));
	   FDT->AssignHits(r,0,NSeqsSeqSet(data));
       	   SetUpSearch(data);
	   if(Search( )) {
       		ProcessHits( );
       		ConvertViaTemplate( );
       		// if(OutputMAPS) PutMAPInputFile( );
       		// PutSeqs( );
       		// PutAln( ); 
		// PutSummary( ); 
		StoreSummary(r); 
		// PutMainAln( );
       		CleanUpSearch();
	   }
	} 
	fprintf(stdout,"\n%s vs %s:\n",database,infile);
	FDT->Put(stdout); fprintf(stdout,"\n\n");
	for(f=1; f <= Number; f++) TotalNilCMSA(IN_CMA[f]); free(IN_CMA);
	return TRUE;
}

//************************** CompareCMSAs Routines ***************************8
BooLean	mgs_typ::CheckInput(Int4 N_cma, Int4 N_rows, cma_typ *CMSA,char **Name_row)
{
	Int4	f,r;
	if(N_cma > N_rows){
		fprintf(stderr,"Number of cmafiles (%d) inconsistent with Hpt rows (%d)\n",N_cma,N_rows);
		exit(1);
	} else if(N_cma <= N_rows){	// Make sure all CMAs and Hpt rows are in the right order.
	     for(r=f=1; f <= N_cma; r++,f++){
		if(r > N_rows) {
			fprintf(stderr,"cmafile names (%d) inconsistent with Hpt (%d)\n",N_cma,N_rows);
			exit(1);
		}
		if(strcmp(NameCMSA(CMSA[f]),Name_row[r]) !=0){
			fprintf(stderr,"%3d: NameRow=%s missing\n",r,Name_row[r]);
			f--; 
		} else fprintf(stderr,"%3d: NameCMSA=%s; NameRow=%s\n",r,NameCMSA(CMSA[f]),Name_row[r]);
	     }
 	} return TRUE;
	
}

#if 0
Int4	mgs_typ::Intersection(ss_type P1, ss_type P2)
// subseq ...
{
	Int4	n,m,Num=0;
        for(n = 1; n <= NSeqsSeqSet(P1); n++){
            e_type E1 = SeqSetE(n,P1); 
	    BooLean found=FALSE;
            for(m=1; m <= NSeqsSeqSet(P2); m++){
                e_type E2 = SeqSetE(m,P2);
                if(IsSubSeq(E1,E2)){ found=TRUE; break; } 
            } if(found) Num++; 
        } return Num;
}
#else

Int4	mgs_typ::Intersection(ss_type P1, ss_type P2)
{
	Int4	n,m,Num=0;
        for(n = 1; n <= NSeqsSeqSet(P1); n++){
            e_type E1 = SeqSetE(n,P1);
	    BooLean found=FALSE;
            for(m=1; m <= NSeqsSeqSet(P2); m++){
                e_type E2 = SeqSetE(m,P2);
		if(IdentSeqs(E1,E2)){ found=TRUE; break; } 
            } if(found) Num++; 
        } return Num;
}
#endif

BooLean	mgs_typ::CompareCMSAs(Int4 Ni,char **iR,char **iName, char *TypeOfSetI,
					Int4 No, char **oR,char **oName,char *TypeOfSetO)
{
	Int4 Number,f,g,r,s,Number2;
  	FDT= new fdt_typ(Ni,iR,iName,TypeOfSetI,No, oR,oName,TypeOfSetO);
	FILE *fp=open_file(database,"_new.mma","r");
	cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	CheckInput(Number,No,IN_CMA,oName);

	fp=open_file(infile,"_new.mma","r");
	// fp=open_file(infile,".sma","r");
	cma_typ *CMPR_CMA=MultiReadCMSA(fp,&Number2,AB); fclose(fp);
	CheckInput(Number2,Ni,CMPR_CMA,iName);

	for(r=f=1; f <= Number; r++,f++){
	   ss_type data=TrueDataCMSA(IN_CMA[f]);
	   // fprintf(stderr,"NameOut=%s; NameRow=%s\n",NameSeqSet(data),oName[f]);
	   if(strcmp(NameCMSA(IN_CMA[f]),oName[r]) !=0){ f--; continue; }
	   assert(NSeqsSeqSet(data) > 0);
	   // fprintf(stderr,"DEBUG: seqs = %d\n",NSeqsSeqSet(data));
	   FDT->AssignHits(r,0,NSeqsSeqSet(data));
	   for(s=g=1; g <= Number2; s++,g++){
	      ss_type data2=TrueDataCMSA(CMPR_CMA[g]);
	      if(strcmp(NameCMSA(CMPR_CMA[g]),iName[s]) !=0){ g--; continue; }
	      Int4 hits=Intersection(data,data2);
	      FDT->AssignHits(r,s,hits);
	   }
	} 
	fprintf(stdout,"\n%s vs %s:\n",database,infile);
	FDT->Put(stdout); fprintf(stdout,"\n\n");
	for(f=1; f <= Number; f++) TotalNilCMSA(IN_CMA[f]); free(IN_CMA);
	return TRUE;
}



