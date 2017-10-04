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

#define USAGE "usage: curated_srch query dbsfile [options]\n\
      NOTE: Put concatenated cma files in query.cma \n\
            with the last being the Template for the rest.\n\
      options:\n\
        -C<int>   Cutoff score for eliminating overlaps (-P option)(default: 30)\n\
        -D        Don't SEG query sequence\n\
        -e<real>  E-value cutoff for showing scores and alignments\n\
        -H<int>   heapsize (default 500000)\n\
        -I<int>:<int> - flanking region lengths for sequence hits(default: 0:0)\n\
        -L        do low complexity masking of database sequences\n\
        -L=<int>  Omit the last <int> cma files from output cma file\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -m=[<char>:<int>..<int>;]<char>:<int>..<int>.  generate output for a vsi file\n\
	          (e.g., -m=R:4..41,O:54..76,Y:81..104,G:112..137,B:141..156.)\n\
        -O=<str>  cma files for families omitted from output cma file (sample str: \"3,5-7,9,11-17\")\n\
        -o        output sequence hits in separate files for each input cma file\n\
        -P        Don't use PatchSAP() option to combine gapped HSPs into one\n\
        -r        Don't detect repetitive domains\n\
        -S        Save *.chk files for a future run\n\
        -s        output subsequence instead of full sequences\n\
        -T<int>   blast word hit threshold (default 11)\n\
        -t<real>  trim alignments at ends where marg. prob. is > <real>\n\
        -test=<int>  test mode (default: 0)\n\
                   test=0:  don't allow deletions next to insertions\n\
                   test=1:  allow deletions next to insertions\n\
        -U        Use *.chk files instead of cma files\n\
        -w<int>   alignment width (default 60)\n\
        -X<int>   X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
        -z        see posMatrix\n\
\n"


#define MAX_NUM_CMA_FILES 1000
#if 0	//**********************************************************************
	Curated cma search:

	1. Read in template cmafile and obtain query sequences (from TrueSeq).

	2. 

#endif  //**********************************************************************
int	curated_srch(int argc, char *argv[])
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
	BooLean	output_by_class=FALSE;
	BooLean output_subseqs=FALSE,print_names_only=FALSE,find_repeats=TRUE;
	Int4	*IGNORE_THESE=0,num_bgd_files=0;
	char	mask_dbs=' ',*checkin=0,*checkout=0;
	FILE	*seqfile=NULL;
	char	*cmafile=0;
	double	misalncut=2.0;
	Int4	num_cma_files,overlap_cutoff=30,start_arg=3;
	UInt4	max_length_input_seqs=300000;
	Int4	NumRegions=0;
	Int4	*Start=0,*End=0,test=0;
	char	*colors;
	BooLean	RemoveChkFiles=TRUE;

	if(argc < 3) print_error(USAGE);
	//***************************** 1. Read input files *******************************
	start_arg=3;
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	FILE    *fp=open_file(argv[1],".cma","r");
	cma_typ *IN_CMA=MultiReadCMSA(fp,&num_cma_files,A);
	fclose(fp);

	if(num_cma_files < 2 || num_cma_files >= MAX_NUM_CMA_FILES) print_error(USAGE);
	IN_CMA[0]=IN_CMA[num_cma_files]; // move Template to front.
	if(NumSeqsCMSA(IN_CMA[0]) != num_cma_files) {
		fprintf(stderr,"\n*********** Template CMA input error! ***********\n\n");
		fprintf(stderr,"NumSeqsCMSA(IN_CMA[0]) = %d; num_cma_files = %d\n\n",
			NumSeqsCMSA(IN_CMA[0]),num_cma_files);
		print_error(USAGE);
	}
	IN_CMA[num_cma_files]=0; num_cma_files--; 

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
		if(argv[arg][2] != '='){
		  mask_dbs='x'; 
		} else {
		  Int4 omit_at_end=0;
		  if(sscanf(argv[arg],"-L=%d",&omit_at_end) != 1){
                        print_error(USAGE);
		  } else if(omit_at_end > 0 && omit_at_end < num_cma_files){
		     NEW(IGNORE_THESE,num_cma_files+3,Int4);
		     while(omit_at_end > 0){
			omit_at_end--;
			IGNORE_THESE[num_cma_files - omit_at_end] = TRUE;
		     }
		  } else { print_error(USAGE); } argv[arg][1] = ' ';
		}
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
	     case 'O':  {
		num_bgd_files++;
		if(argv[arg][2] != '=') print_error(USAGE);
		char *str=argv[arg]+3;
		Int4 *values,j;
		NEW(values,num_cma_files+20,Int4);
		// WARNING: may get a segmentation error if values > num_cma_files!!!
		Int4 n=ParseIntegers(str, values, "cma_gblastpgp: -O option input error");
		if(n > num_cma_files) print_error(USAGE);
		NEW(IGNORE_THESE,num_cma_files+3,Int4);
                for(Int4 i=1; i<= n; i++){
                    if((j=values[i]) > 0 && j <= num_cma_files){
			if(IGNORE_THESE[j] > 0) print_error("-O option error: overlapping input sets\n");
			IGNORE_THESE[j]=num_bgd_files; 
                    } else print_error(USAGE);
                } free(values);
	       } break;
             case 'P': use_patch_sap=FALSE; argv[arg][1] = ' '; break;
             case 'r': find_repeats=FALSE; argv[arg][1] = ' '; break;
             case 's': output_subseqs=TRUE; argv[arg][1] = ' '; break;
             case 'S': RemoveChkFiles=FALSE; argv[arg][1] = ' '; break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,USAGE); break;
	     case 't': 
		if(sscanf(argv[arg],"-test=%d",&test) == 1){
		   if(test > 1 || test < 0) print_error(USAGE);
		} else misalncut=RealOption(argv[arg],'t',0.1,1.0,USAGE);
		break;
             case 'o': output_by_class= TRUE; argv[arg][1] = ' '; break;
             case 'U': use_chk_files_only= TRUE; RemoveChkFiles=FALSE; argv[arg][1] = ' '; break;
	     case 'w': width=IntOption(argv[arg],'w',5,500000,USAGE); break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,USAGE); break;
             case 'Z': print_names_only=TRUE; break;
             case 'z': see_smatrix= TRUE; break;
             default: print_error(USAGE);
           }
        }

	//***************************** 2. Initialize search *******************************
	e_type	qE,QueryE[MAX_NUM_CMA_FILES],trueE,fakeE;
	ss_type	data,Data[MAX_NUM_CMA_FILES];
	sap_typ	sap,SAP[MAX_NUM_CMA_FILES];
	cma_typ	cma,CMA[MAX_NUM_CMA_FILES];
	char 	str[100],*Checkpoint[MAX_NUM_CMA_FILES];
	Int4	II,maxrounds;

  if(IGNORE_THESE || print_names_only){
   for(Int4 s=1; s <= num_cma_files; s++){ 
     if(IGNORE_THESE && IGNORE_THESE[s]){
	  fprintf(stderr,"---------------------- %d: %s (%d)(ignore).\n",
		s,NameCMSA(IN_CMA[s]),IGNORE_THESE[s]);
     } else fprintf(stderr,"================================ %d: %s.\n", s,NameCMSA(IN_CMA[s])); 
   } 
  }
  if(!IGNORE_THESE) NEW(IGNORE_THESE,num_cma_files+3,Int4);
  if(print_names_only) exit(0);
  for(II=0; II <= num_cma_files; II++){
	Checkpoint[II]=0; QueryE[II]=0; Data[II]=0; SAP[II]=0; CMA[II]=0;
  }
  Data[0]=MakeSeqSet(argv[start_arg-1],max_length_input_seqs,A);
  // Search database for hits using distinct queries for each cma file.
  for(II=1; II <= num_cma_files; II++){
	CMA[II]=IN_CMA[II];
	cma=CMA[II];
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		// need to renumber sequence ids in alignment for blast.
		trueE = TrueSeqCMSA(sq,cma);
		fakeE = FakeSeqCMSA(sq,cma);
		EqSeqI(sq,trueE); EqSeqI(sq,fakeE);
	}
	Data[II]=TrueDataCMSA(cma);
  	// ReNumberSeqSet(Data[II]);
	trueE=TrueSeqCMSA(1,cma); fakeE=FakeSeqCMSA(1,cma);
	if(!IdentSeqs(trueE,fakeE) || LenSeq(trueE) != LengthCMSA(1,cma) 
			|| nBlksCMSA(cma) != 1){
		   fprintf(stderr,"file %s_%d.cma: \n\n",argv[1],II);
		   fprintf(stderr,"trueE = %d aa; cma = %d aa \n\n",
					LenSeq(trueE),LengthCMSA(1,cma));
		   PutSeq(stderr,trueE,A);
		   PutSeq(stderr,fakeE,A);
		   print_error("incompatible query sequence");
	}
	QueryE[II]=CopySeq(trueE);
	ChangeInfoSeq(NameCMSA(cma),QueryE[II]); 
	sprintf(str,"%s_%d.chk",argv[1],II);
	Checkpoint[II] = AllocString(str);
  }

	//*********************************** 3. Main search loop  ***********************************
  BooLean	fake_cma_srch;
  if(use_chk_files_only) fake_cma_srch=FALSE; else fake_cma_srch=TRUE;
  for(II=1; II <= num_cma_files; ){
  // start main loop over all *cma files....
	// TESTING FakeSearchCMSA( );
	if(fake_cma_srch){
		checkin=0; 
		checkout=Checkpoint[II]; 
		cma=CMA[II]; data=Data[II]; 
		qE=CopySeq(QueryE[II]); 
		EqSeqI(0,qE);
		maxrounds=2;
	} else {
		cma=0; data=Data[0]; qE=CopySeq(QueryE[II]); 
		EqSeqI(0,qE);
		checkout=0; 
		checkin=Checkpoint[II];
		maxrounds=1;
	}
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,A);
	Int4 **mtx=0;
	if(checkin && use_patch_sap) mtx=gpsi.CopyOFposMatrix();
        do { 
	   sap=0;
// fprintf(stderr,"II=%d; fake_cma_srch=%d; pass=%d\n",II,fake_cma_srch,gpsi.ThisPass( ));
	   if(fake_cma_srch){
		gpsi.KeepQuite();
		// gpsi.SpeakUp( );
	   	sap = gpsi.FakeSearchCMSA(T,cma); 
		fake_cma_srch=FALSE;
	   } else {
		if(gpsi.ThisPass( ) == 0) gpsi.SpeakUp();
		if(misalncut < 1.0) sap = gpsi.TrimmedSearch(T,misalncut,mask_dbs);
	   	else sap = gpsi.Search(T,mask_dbs);
		// if(gpsi.ThisPass( ) == 1 && sap)PutGSeqAlignList(stdout, sap, 60, A);
	   }
	   if(sap) success=TRUE; else break;
	   gpsi.ComputeMatrix(checkout);
	   if(see_smatrix && gpsi.posMatrix){ 
		FILE *ofp = open_file(argv[1],".psm","w");
		gpsi.SeeMatrix(ofp); 
		fclose(ofp);
	   }
	} while(gpsi.NotConverged( ));
	if(checkin){
	   sap_typ head;
	   if(use_patch_sap){ // fix SAP here...to patch together HSPs
		sap=gpsi.StealSAP();
		if(sap){
		   head=MakeSelfAlignedGSAP(qE);
		   head->next = sap;
		   SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data[0]),qE,A); 
		
		   // fprintf(stderr,"DEBUG 1\n");
		   for(Int4 s=0; s < LenSeq(qE); s++) free(mtx[s]); free(mtx);
		   // fprintf(stderr,"DEBUG 2\n");
		} else SAP[II]=0;
		II++;
	   } else {
		sap=gpsi.StealSAP();
		if(sap){
		  SAP[II]=MakeSelfAlignedGSAP(qE);
        	  SAP[II]->next=sap;
		} else SAP[II]=0;
		II++;
	   }
  	   if(!use_chk_files_only) fake_cma_srch=TRUE;
	} 
   } // end loop over cma files

   //************************* 4. process hits ************************
   FILE *fptr=0;
   // merge sap lists deleting weakest redundant hits. 
   //************************* 4a. check for hits ************************
// fprintf(stderr,"DEBUG 1\n");
   for(II=1; II <= num_cma_files; II++){
	sap=SAP[II];
	if(sap){
	  fptr=stdout;
#if 0
          if(printopt == 4) PutMultiGSeqAlign(stderr, sap, width, A);
          else if(printopt == 1){
	        printf("===> File %d: %s.\n",II,NameCMSA(IN_CMA[II]));
		PutGSeqAlignList(stderr, sap->next, width, A);
	  }
#endif
        } // else printf("================== File %d(%s): no hits!\n",II,NameCMSA(CMA[II]));
   }
   if(fptr){
// fprintf(stderr,"DEBUG 2\n");
	//********************* 4b. Find the best matches and remove poor hits ********************
	sap_typ *saps;
	UInt4 *NUM_EACH_CMA=0;
	NEW(NUM_EACH_CMA,num_cma_files+3,UInt4);
	if(find_repeats){	// ******** finding less seqs for some reason...
	  FindBestOverlappingSAPs(SAP,num_cma_files,NSeqsSeqSet(Data[0]),A);
#if 0
   	  for(II=1; II <= num_cma_files; II++){ 
       	    if(SAP[II]){
	        printf("DEBUG: ===> File %d: %s.\n",II,NameCMSA(IN_CMA[II]));
		PutGSeqAlignList(stderr, SAP[II]->next, width, A);
            }
	  }
	  printf("*********** entering RemoveRejectsSAP( ) *************\n");
#endif
	  RemoveRejectsSAP(SAP,num_cma_files,A);
	} else {
#if 0		// passing argv did not do anything, so I deleted it...
	  if(output_by_class){
           LabelBestSAPs(argv[1],SAP,num_cma_files,NSeqsSeqSet(Data[0]),
			use_patch_sap,A,NUM_EACH_CMA);
	  } else {
	   LabelBestSAPs(0,SAP,num_cma_files,NSeqsSeqSet(Data[0]),
			use_patch_sap,A,NUM_EACH_CMA);
	  }
#else
	   LabelBestSAPs(SAP,num_cma_files,NSeqsSeqSet(Data[0]),use_patch_sap,
			A,NUM_EACH_CMA);
#endif
	}
   }
   //************************* 5. Discard empty hits ************************
   Int4 s;
   // X. Create CMA list for merging into Main cma file...
// fprintf(stderr,"DEBUG 3\n");
   for(II=1; II <= num_cma_files; II++){ 
	if(SAP[II] && !AtLeastOneLabeledSAPS(SAP[II]->next)){
		// free SAP[II] here...
		FreeGSeqAlignList(SAP[II]); SAP[II]=0;
	} CMA[II]=0; 
   }
   //************************* 6. Create cma files of best hits ************************
// fprintf(stderr,"DEBUG 4\n");
   BooLean found_seqs=FALSE;
   for(II=1; II <= num_cma_files; II++){
       if(SAP[II]){
   	   found_seqs=TRUE;
           sprintf(str,"%s_%d.cma",argv[1],II);
           fptr=open_file(str,"","w");
           xSeqAlignToCMA(fptr,SAP[II],left_flank,right_flank,A);
           fclose(fptr);
#if 0
	   if(IGNORE_THESE[II]==0) { CMA[II]=ReadCMSA2(str,A); }
       	   else { CMA[II]=0; }
#else
	   CMA[II]=ReadCMSA2(str,A);
	   RenameCMSA(NameCMSA(IN_CMA[II]),CMA[II]);
#endif
       } else { CMA[II]=0; }
   }
   // print out best hits only:
   //************************* 7. output sequence alignments? ************************
// fprintf(stderr,"DEBUG 5\n");
   FILE *main_sq_fp=0;
   if(found_seqs){
	main_sq_fp=open_file(argv[1],"_Main.seq","w");
   }
   for(II=1; II <= num_cma_files; II++){ 
       	 if(SAP[II]){
            if(printopt == 4) PutMultiGSeqAlign(stderr, SAP[II], width, A);
            else if(printopt == 2) PutGSeqAlignList(stdout, SAP[II]->next, width, A);
            else if(printopt == 1){
	        printf("===> File %d: %s.\n",II,NameCMSA(IN_CMA[II]));
		PutGSeqAlignList(stderr, SAP[II]->next, width, A);
	    }
#if 1
   	    if(output_by_class){
   		FILE *sq_fp=0;
           	sprintf(str,"%s_%d.seq",argv[1],II);
           	sq_fp=open_file(str,"","w");
		PutSeqsListGSAP(sq_fp, SAP[II],A);
   		fclose(sq_fp);

	    }
	    if(main_sq_fp) PutSeqsListGSAP(main_sq_fp, SAP[II],A);
#endif
         }
   }
   if(main_sq_fp) fclose(main_sq_fp);

   //************************* 8. Use template to convert CMA files ************************
   //**************************** Key point in algorithm ***********************************
// fprintf(stderr,"DEBUG 6 (test=%d)\n",test);
   cma_typ *OUT_CMA=0;
   if(test > 0) OUT_CMA=ConversionViaTemplateCMSA(IN_CMA[0],CMA);
   else OUT_CMA=ConversionViaTemplateCMSA3(IN_CMA[0],CMA);
   //******************* CALL recursively for hierarchical alignments? **********************
// fprintf(stderr,"DEBUG 7\n");
   for(II=1; II <= num_cma_files; II++){	
	if(OUT_CMA[II]){
	    char str2[500];
	    sprintf(str2,"%s_%d.cma",argv[1],II); 
#if 1	// output regions in vsi format...
	// -m=R:6..20;O:32..45;Y55..68.
	if(printopt==2){
	   Int4 sq=2;
	   PutVSIregionsCMSA(stdout,sq,colors,Start,End,NumRegions,OUT_CMA[II]);
	}
#endif
#if 1	// Rename output cma files according to input files names.
	    RenameCMSA(NameCMSA(IN_CMA[II]),OUT_CMA[II]);
#endif
	    WriteCMSA(str2,OUT_CMA[II]);
	}
   } 
   //************************* 9. Summary of output files ************************
// fprintf(stderr,"DEBUG 8\n");
   printf("\nDatabase search results:\n");
   for(II=1; II <= num_cma_files; II++){
	if(OUT_CMA[II]){ 
	  if(IGNORE_THESE[II]) printf("---------------------> File %d: %s (%d hits).\n",
                        II,NameCMSA(IN_CMA[II]),NumSeqsCMSA(OUT_CMA[II])-1);
	  else printf("=====================> File %d: %s (%d hits).\n",
                        II,NameCMSA(IN_CMA[II]),NumSeqsCMSA(OUT_CMA[II])-1);
	} else {
	   printf("--- File %d: %s (no hits).\n",
                        II,NameCMSA(IN_CMA[II]));
	}
   }
// fprintf(stderr,"DEBUG 9\n");
   for(s=0,II=1; II <= num_cma_files; II++){	// compress array for output...
	if(!IGNORE_THESE[II] && OUT_CMA[II]!=0) {
		s++; CMA[s]=OUT_CMA[II]; 
		ReNameCMSA(NameCMSA(IN_CMA[II]),CMA[s]);
	}
   }
   set_typ *set;
   NEW(set,s+3,set_typ);
   //************************* 10. output fasta sequence files ************************
// fprintf(stderr,"DEBUG 10\n");
   for(II=1; II <= s; II++){
	   set[II]=MakeSet(NumSeqsCMSA(CMA[II])+1);
	   FillSet(set[II]); DeleteSet(0,set[II]); DeleteSet(1,set[II]); 
#if 0
	   if(output_by_class){
		ss_type data0=TrueDataCMSA(CMA[II]);
		for(Int4 sq=2; sq <= NSeqsSeqSet(data0); sq++){
			PutSeqSetE(sq_fp,sq, data0);
		}
	   }
#endif
   }
   //************************* 11. output Main alignment file ************************
// fprintf(stderr,"DEBUG 11\n");
   if(s > 0){
     fptr = open_file(argv[1],"_Main.cma","w");
     PutMergedCMSA(fptr,s,set,CMA,NULL); 
     fclose(fptr);
   } else fprintf(stderr,"No hits found.\n");
// fprintf(stderr,"DEBUG 12\n");
   for(II=1; II <= s; II++) NilSet(set[II]);
   for(II=0; II <= num_cma_files; II++){
	if(IN_CMA[II]) TotalNilCMSA(IN_CMA[II]);	// merging ruins the seqset!!!!!
	if(OUT_CMA[II]) TotalNilCMSA(OUT_CMA[II]);	// merging ruins the seqset!!!!!
	if(Checkpoint[II]){
	  if(RemoveChkFiles){
	    char cmd_str[200];
	    sprintf(cmd_str,"rm -f %s",Checkpoint[II]);
	    system(cmd_str);
	  } free(Checkpoint[II]);
	}
	if(QueryE[II]) NilSeq(QueryE[II]); 
   } NilSeqSet(Data[0]); NilAlpha(A);
   double runtime = difftime(time(NULL),time1);
   fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
   if(success) return 0;
   else return 1;
}


