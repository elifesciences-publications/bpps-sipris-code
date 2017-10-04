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

#define USAGE "usage: matblast template_alignment database [options]\n\
      NOTE: The template alignment file must be in fasta format \n\
      options:\n\
        -e<real>  E-value cutoff for showing scores and alignments\n\
        -False=<str> input the 'false positive' cma file so that number can be computed\n\
        -false=<int> input the number of 'false positive' cma files included\n\
                  (these must be placed between the true cma files and the template)\n\
        -H=<int>  heapsize (default 500000)\n\
	-l        do low complexity masking of database sequences\n\
	-L        do low complexity and coiled coil masking of database sequences\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -sense=<real> heuristic sensitivity (default: 0.90)\n\
        -X<int>   X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
        -z        see posMatrix\n\
\n"

// omit below for now...
//        -O        output sequence hits in separate files for each profile\n\


#if 0	//**********************************************************************
#endif  //**********************************************************************
int	matblast_srch(int argc, char *argv[])
{ return matblast_srch(argc, argv,USAGE); }

int	matblast_srch(int argc, char *argv[],const char *usage)
{
	Int4    T=11;
	time_t	time1=time(NULL);
        double	x_parameter=25.0;
	Int4	x;
	UInt4	hpsz=500000;
	Int4	printopt=0,overlap_cutoff=30;
	double	Ethresh=0.001,Ecutoff=0.001;
	BooLean	see_smatrix=FALSE,use_patch_sap=TRUE;;
	BooLean	success=FALSE,*IGNORE_THIS,find_repeats=TRUE;
	BooLean print_names_only=FALSE,SeparateFiles=FALSE;
	char	*checkin=0,*checkout=0;
	Int4	num_cma_files,start_arg=3,num_false_pos=0,FirstFalse;
	char	FalsePosFile[200];
	char    mask_dbs=' ';
	cma_typ	cma,true_cma;
	FILE	*fp;
	double  fractionMS=0.90;
	sap_typ	sap,*SAP=0;
	Int4    left_flank=0,right_flank=0;

	if(argc < 3) print_error(usage);
	//***************************** 1. Read input files *******************************
	start_arg=3;
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	fp=open_file(argv[1],".tpl","r");
	cma=ReadCMSA(fp,A); fclose(fp);

	num_cma_files = NumSeqsCMSA(cma);

	for(Int4 arg=0; arg < argc; arg++) { fprintf(stderr,"%s ",argv[arg]); }
	fprintf(stderr,"\n\n");
	// fprintf(stderr,"num_cma_files = %d\n",num_cma_files);
	//***************************** 2. Get input options *******************************
        for(Int4 arg = start_arg; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(usage);
           switch(argv[arg][1]) {
	     case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,USAGE); break;
	     case 'F': {	
		cma_typ fcma=0;
		  if(sscanf(argv[arg],"-False=%s",FalsePosFile) == 1){
			fcma=ReadCMSA2(FalsePosFile,A);
			num_false_pos = NumSeqsCMSA(fcma);
			if(num_false_pos < 1 || num_false_pos >= num_cma_files){
				print_error(usage);
			} else NilCMSA(fcma);
		  } else print_error(usage);
		} break;
	     case 'f': 	
		if(sscanf(argv[arg],"-false=%d",&num_false_pos) == 1){
			if(num_false_pos < 1 || num_false_pos >= num_cma_files){
				print_error(usage);
			}
		} else print_error(usage);
		break;
	     case 'H': 	
			hpsz=IntOption(argv[arg],'H',1,1000000,usage); 
		break;
	     case 'l': mask_dbs='x'; break;
	     case 'L': mask_dbs='b'; break;
	     case 'm': 	
		{
		  if(argv[arg][2] != '='){
			printopt=IntOption(argv[arg],'m',0,6,usage); 
		  } else print_error(usage);
		} break;
	     case 'O': 
		  if(argv[arg][2] == 0) SeparateFiles=TRUE; 
		  else print_error(usage);
		break;
	     case 's': 	
                if(sscanf(argv[arg],"-sense=%f",&fractionMS) == 1){
                        if(fractionMS < 0.1 || fractionMS > 1.0) print_error(usage);
                } else print_error(usage);
		break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,usage); break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,usage); break;
             case 'Z': print_names_only=TRUE; break;
             case 'z': see_smatrix= TRUE; break;
             default: print_error(usage);
           }
        }

	//***************************** 2. Initialize search *******************************
	e_type	qE,*QueryE,trueE,fakeE;
	ss_type	data; 
	char 	str[100],*Checkpoint;
	Int4	II,maxrounds;

  NEW(SAP,num_cma_files+3,sap_typ);
  NEW(QueryE,num_cma_files+3,e_type);
  NEW(IGNORE_THIS,num_cma_files+3,BooLean);
  for(II=0; II <= num_cma_files; II++){ QueryE[II]=0; SAP[II]=0; }
  // check to make sure that input file is okay.
  // SetLevelCMSA(num_false_pos,cma);

  FirstFalse=num_cma_files+1;
  if(num_false_pos > 0){	// then create a template profile lacking the false positives.
	BooLean *skip=0;
	Int4 N = NumSeqsCMSA(cma); 
	NEW(skip,N+3,BooLean);
	FirstFalse=N-num_false_pos+1;
	for(Int4 i=N - num_false_pos + 1; i <= N; i++){
		skip[i]=TRUE; IGNORE_THIS[i] = TRUE;
	}
	FILE *tmpfp=tmpfile(); PutSelectCMSA(tmpfp,skip,cma); free(skip);
	rewind(tmpfp); true_cma=ReadCMSA(tmpfp,A); fclose(tmpfp);
  } else true_cma = cma;
  
  trueE=TrueSeqCMSA(1,cma); fakeE=FakeSeqCMSA(1,cma);
  QueryE[0]=CopySeq(trueE);
  if(!IdentSeqs(trueE,fakeE) || LenSeq(trueE) != LengthCMSA(1,cma)
			|| nBlksCMSA(cma) != 1){
	fprintf(stderr,"file %s_%d.cma: \n\n",argv[1],II);
	fprintf(stderr,"trueE = %d aa; cma = %d aa \n\n",
					LenSeq(trueE),LengthCMSA(1,cma));
	PutSeq(stderr,trueE,A); PutSeq(stderr,fakeE,A);
	print_error("incompatible query sequence");
  }

  for(II=1; II <= num_cma_files; II++){
	trueE=TrueSeqCMSA(II,cma); 
	QueryE[II]=CopySeq(trueE);
	// ChangeInfoSeq(NameCMSA(cma),QueryE[II]); 
  }
  sprintf(str,"%s.mpa",argv[2]); Checkpoint = AllocString(str);

  //******************************** 3. Main search loop  *******************************
// fprintf(stderr,"DEBUG 1\n");
  FILE *chkfp=0;
  set_typ SkipSet=0,SetM,SetH,SetE;
  // chkfp=open_file(Checkpoint,"","w");
  checkout=Checkpoint; checkin=0;
  data=MakeSeqSet(argv[2],200000,A); ReNumberSeqSet(data);
  ss_type Data=0;
  {
     // Data=TrueDataCMSA(true_cma);
	// not working with true_cma...FIND OUT WHY!
	// missing one sequence...in output file below...
     Data=TrueDataCMSA(cma);
     ReNumberSeqSet(Data);
     qE=CopySeq(QueryE[0]); maxrounds=2;
     gpsi_type gpsi(qE,Data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
// fprintf(stderr,"DEBUG 1a\n");
	
     // gpsi.SpeakUp( );
     gpsi.KeepQuite();
     // sap = gpsi.FakeSearchCMSA(T,true_cma); // do this with template alignment...
     sap = gpsi.FakeSearchCMSA(T,cma); // do this with template alignment...
// fprintf(stderr,"DEBUG 1b\n");
     if(sap) success=TRUE; // else break;
     // else print_error("This should not happen");
     // gpsi.ComputeMatrix(chkfp, checkout);
     gpsi.ComputeMatrix(checkout);
  }
// fprintf(stderr,"DEBUG 2\n");
  // checkout=0; 
  // if(checkout) chkfp=open_file(Checkpoint,"","a+b");
  // 2. Save the sequences found in this search to search against the remaining profiles
  BooLean *skip;
  NEW(skip,num_cma_files+5,BooLean);
  fprintf(stderr,"\nNo.: %-10s     %-11s   %-11s   %-11s   (time)(TotalHits)\n",
                        "name","Hits","Extend","Trigger");
  Int4 DatabaseSize=NSeqsSeqSet(data);
  for(II=1; II <= num_cma_files; II++){ // start main loop over all *cma files....
	clock_t timeS=clock();
	if(II == 1){ checkin=Checkpoint; checkout=0;
	} else if(II == 2){
		chkfp=open_file(Checkpoint,"","a+b");
		checkin=0; checkout=Checkpoint;
	}
	qE=CopySeq(QueryE[II]); 
	EqSeqI(0,qE);
	maxrounds=2;
	// gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
        if(II > 1){ gpsi.ProvideSkipSet(SkipSet); }
	gpsi.CreateHitSet(); gpsi.CreateExtendSet(); gpsi.CreateMarginalSet(fractionMS);
	// gpsi.SpeakUp( );
	gpsi.KeepQuite();
           // if(misalncut < 1.0) sap = gpsi.TrimmedSearch(T,misalncut,mask_dbs); else 
	sap = gpsi.Search(T,mask_dbs);
	if(sap){ 
		success=TRUE; // else break;
	   	if(!gpsi.ComputeMatrix(chkfp, checkout)){
			skip[II]=TRUE;	// probably query only found itself...
		}
	} else { skip[II]=TRUE;	}
	SetM=gpsi.RtnMarginalSet(); SetH=gpsi.RtnHitSet(); SetE=gpsi.RtnExtendSet();
	fprintf(stderr,"%3d: ",II);
	char id_str[15]; StrSeqID(id_str, 10, qE);
	fprintf(stderr,"%-10s     %-11d   %-11d   %-11d",id_str,
                        CardSet(SetH),CardSet(SetE),CardSet(SetM));
	fprintf(stderr," (%0.2f s)\n",(double) (clock()-timeS)/(double)CLOCKS_PER_SEC);
	if(II == 1){
	   SkipSet=SetM; SetM=0;
#if 0
	   Int4 NumPruned=CardSet(SkipSet);
	   if(NumPruned == 0){ fprintf(stderr,"no hits found\n"); exit(0); }
		fprintf(stderr,"   pruning input set from %d to %d sequences (%.1lf%%)\n",
                        NSeqsSeqSet(Data),NumPruned,
                                100.0*(double)NumPruned/(double)NSeqsSeqSet(Data));
#endif
	} else if(SetM) NilSet(SetM); 
// fprintf(stderr,"DEBUG 3\n");
	if(SetH) NilSet(SetH); if(SetE) NilSet(SetE);
//********************************* GET SAPs ***********************************
        sap_typ head;
#if 0
        if(use_patch_sap){ // fix SAP here...to patch together HSPs
                sap=gpsi.StealSAP();
                if(sap){
                   head=MakeSelfAlignedGSAP(qE);
                   head->next = sap;
                   SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data),qE,A);
                   for(Int4 s=0; s < LenSeq(qE); s++) free(mtx[s]); free(mtx);
                } else SAP[II]=0;
        } else 
#endif
	{
                sap=gpsi.StealSAP();
                if(sap){
                  SAP[II]=MakeSelfAlignedGSAP(qE);
                  SAP[II]->next=sap;
                } else SAP[II]=0;
        }
   } // end loop over sequences in Template alignment 
   if(chkfp) fclose(chkfp);


   //******************************** Output CMA files ******************************
   //******************************** Output CMA files ******************************
   //******************************** Output CMA files ******************************
   //************************* 4. process hits ************************
   // merge sap lists deleting weakest redundant hits.
   //************************* 4a. check for hits ************************
   BooLean found_hits=FALSE,found_false=FALSE;
   for(II=1; II <= num_cma_files; II++){
         if(IGNORE_THIS[II] && SAP[II]) found_false=TRUE;
         else if(SAP[II]) found_hits=TRUE;
   }
   if(found_hits){
        //********************* 4b. Find the best matches and remove poor hits ********************
        sap_typ *saps;
        UInt4 *NUM_EACH_CMA=0;
        NEW(NUM_EACH_CMA,num_cma_files+3,UInt4);
        if(find_repeats){       // ******** finding less seqs for some reason...
          FindBestOverlappingSAPs(SAP,num_cma_files,NSeqsSeqSet(data),A);
          RemoveRejectsSAP(SAP,num_cma_files,A);
        } else {
          LabelBestSAPs(SAP,num_cma_files,NSeqsSeqSet(data),use_patch_sap,A,NUM_EACH_CMA);
        }
   }
   //************************* 5. Discard empty hits ************************
   Int4 s;
   for(II=1; II <= num_cma_files; II++){
        if(SAP[II] && !AtLeastOneLabeledSAPS(SAP[II]->next)){
                FreeGSeqAlignList(SAP[II]); SAP[II]=0;
        }
   }
   //************************* 6. Create cma files of best hits ************************
   BooLean found_seqs=FALSE;
   cma_typ *CMA;
   NEW(CMA,num_cma_files+5,cma_typ);

   fp=open_file(argv[2],".tpl","w");
   PutSelectCMSA(fp,skip,cma); fclose(fp);

   BooLean *Skip;
   NEW(Skip,num_cma_files+5,BooLean);
   sprintf(str,"%s_map",argv[2]);
   FILE *cmafp=open_file(str,".cma","w");
   FILE	*tsq_fp=0;
   set_typ InputSeqSet=MakeSet(DatabaseSize+1); ClearSet(InputSeqSet);
   for(II=1; II <= num_cma_files; II++){
       if(SAP[II] && !skip[II]){
           found_seqs=TRUE;
           fp=tmpfile(); 
#if 0	   // fix deletions at ends of sequence.
	   MakeGlobalSeqAlign(SAP[II]);  
#endif
	   xSeqAlignToCMA(fp,SAP[II],left_flank,right_flank,A);
           rewind(fp); CMA[II]=ReadCMSA(fp,A); fclose(fp);
           StrSeqID(str, 20, QueryE[II]); RenameCMSA(str,CMA[II]);
	   if(II > 1){
		PutCMSA(cmafp,CMA[II]);
#if 1	// output sequence hits...
	   	sap_typ sap,head=SAP[II]->next;
		if(tsq_fp==0){
                  tsq_fp=open_file(argv[2],"_map.seq","w");
                }
		for(sap=head; sap!=NULL; sap=sap->next){
                  e_type sE = SubjectSeqGSAP(sap);
		  Int4 sq = SeqI(sE);
                  if(!MemberSet(sq,InputSeqSet)){
                        if(tsq_fp) PutSeq(tsq_fp,sE,A);
                        AddSet(sq,InputSeqSet);
                  }
                }
#endif
#if 0	// put as separate cma files in one file right now...
		if(SeparateFiles){
           	  FILE *tmpfp=open_file(str,".cma","w");
   		  fclose(tmpfp);
	        }
#endif
	   }
       } else {
	  Skip[II]=TRUE; CMA[II]=0; 
	  if(II >= FirstFalse) num_false_pos--;
       }
   }
   NilSet(InputSeqSet);

   RenameCMSA("template",cma);
   SetLevelCMSA(num_false_pos,cma);
   PutSelectCMSA(cmafp,Skip,cma);  free(Skip);
   fclose(cmafp);

   //********************************* Run mkmaps *****************************
   //********************************* Run mkmaps *****************************
   //********************************* Run mkmaps *****************************

   free(skip);
   //************************* Free up memory ************************
   for(II=0; II <= num_cma_files; II++){
	if(QueryE[II]) NilSeq(QueryE[II]);
   }
   TotalNilCMSA(cma); NilAlpha(A);
   double runtime = difftime(time(NULL),time1);
   fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
   if(success) return 0;
   else return 1;
}


