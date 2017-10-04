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

sap_typ	*MergeSAP(char *outfile,sap_typ *HEAD, Int4 num_sap, Int4 num_seqs, e_type qE,
	BooLean get_all_hits,a_type AB,BooLean output_subseqs,
	UInt4 *NUM_EACH_CMA,Int4 num_bgd_files, Int4 *IGNORE_THESE)
// merge the SAPs into one list
// Move this to my_ncbi because of private information.
// ScoreGSAP(sap); EvalueGSAP(sap); NextGSAP(sap); 
// SubjectSeqGSAP(sap); SeqI(E); SubjectSeqGSAP(sap);
// sap_typ ConcatenateGSAP(sap_typ head, sap_typ sap);
{
	sap_typ	*head,*BEST_SAP,tail,sap,gsap;
	Int4	len,score,sq,s,*BEST_CMA,*BEST_SCORE,min_len;
	UInt4	J,start,end,n,**NUM_RPTS;

	e_type	sE;
	NEW(head,num_bgd_files+4,sap_typ);
	min_len = LenSeq(qE)/2;
	NEW(BEST_SAP,num_seqs+3,sap_typ);
	NEW(BEST_CMA,num_seqs+3,Int4);
	NEW(BEST_SCORE,num_seqs+3,Int4);
	NEWP(NUM_RPTS,num_sap+3,UInt4);
	for(s=1; s <= num_sap; s++){
	    NEW(NUM_RPTS[s],num_seqs+3,UInt4);
	    n=0;
	    for(gsap=HEAD[s]; gsap!=NULL; gsap=gsap->next){
		sE=SubjectSeqGSAP(gsap); sq=SeqI(sE); score=ScoreGSAP(gsap);
		if(sq > num_seqs) print_error("sequence id error in MergeSAP()");
		len = GetStartEndGSAP(gsap, &start,&end);
		if(len < min_len) continue;
		NUM_RPTS[s][sq]++;
		if(score > BEST_SCORE[sq]){ 
		   BEST_SCORE[sq]=score; BEST_CMA[sq]=s; BEST_SAP[sq]=gsap; 
		} 
	    } 
	}
	tail=head[0]=MakeSelfAlignedGSAP(qE);
	tail->next=0; 
	for(Int4 i=1; i <= num_bgd_files; i++){
	  tail=head[i]=MakeSelfAlignedGSAP(qE);
	  tail->next=0; 
	}
	UInt4 total_hits=0;
	sap_typ	*FINAL_SAP; 

   if(get_all_hits){ 		// get all saps from best cma...
	for(sq=1; sq <= num_seqs; sq++){
	   for(s=1; s <= num_sap; s++){ if(BEST_CMA[sq]==s) total_hits+=NUM_RPTS[s][sq]; }
	   if(NUM_EACH_CMA && BEST_CMA[sq] > 0){
	     assert(BEST_CMA[sq] <= num_sap);
	     NUM_EACH_CMA[BEST_CMA[sq]]++;
	   }
	}
	NEW(FINAL_SAP,total_hits +3,sap_typ);

	for(J=0,s=1; s <= num_sap; s++){
	    for(gsap=HEAD[s]; gsap!=NULL; gsap=gsap->next){
		sE=SubjectSeqGSAP(gsap); sq=SeqI(sE);
		len = GetStartEndGSAP(gsap, &start,&end);
		if(len < min_len) continue;
	    	if(s==BEST_CMA[sq]){ J++; FINAL_SAP[J]=gsap; }
	    }
	} assert(J == total_hits);

     if(outfile){
      if(output_subseqs){    //*********** output subseqs only ************
	Int4 left=20,right=20;
	char str[100];
	for(s=1; s <= num_sap; s++){
	  sprintf(str,"_%d.seq",s);
          FILE *sfp=0;
	  for(J=1; J <= total_hits; J++){
	    sE=SubjectSeqGSAP(FINAL_SAP[J]);
	    assert(sE); sq=SeqI(sE);
	    if(s==BEST_CMA[sq]){
	      if(sfp==0) sfp=open_file(outfile,str,"w"); 
	      PutGSubSeq(sfp,FINAL_SAP[J],left,right,AB);
	    }
	  }
	  if(sfp) fclose(sfp);
	}
      } else { 		     //************* OLD FULL SEQS ******************
	e_type **SEQS;
	Int4	*SeqPt;
	NEWP(SEQS,num_sap+3,e_type); NEW(SeqPt,num_sap+3,Int4);
	for(s=1; s <= num_sap; s++) NEW(SEQS[s],total_hits+3,e_type);
	for(J=1; J <= total_hits; J++){
		sE=SubjectSeqGSAP(FINAL_SAP[J]);
		assert(sE); sq=SeqI(sE);
		s=BEST_CMA[sq];
		assert(s > 0 && s <= num_sap);
		SEQS[s][SeqPt[s]] = sE; SeqPt[s]++;
	}
	for(s=1; s <= num_sap; s++){
	   char str[100];
	   if(SeqPt[s] > 0){
	     sprintf(str,"_%d.seq",s);
             FILE *fptr=open_file(outfile,str,"w"); 
	     for(sq=0; sq < SeqPt[s]; sq++) PutSeq(fptr,SEQS[s][sq],AB);
	     fclose(fptr);
	   }
	   free(SEQS[s]);
	} free(SEQS); free(SeqPt);
      }
     }
	// DEstroy all SAPS where s != BEST_CMA[sq] right HERE...

	sap_typ *bg_sap=0;
	if(num_bgd_files > 0){
		NEW(bg_sap,num_bgd_files+4,sap_typ);
		for(Int4 i=1; i <= num_bgd_files; i++) bg_sap[i]=head[i];
	}
	gsap=head[0];
	for(J=1; J <= total_hits; J++){
#if 1	// eliminate all hits on ignore list...
	    sE=SubjectSeqGSAP(FINAL_SAP[J]);
	    assert(sE); sq=SeqI(sE);
	    s=BEST_CMA[sq];
	    if(num_bgd_files > 0 && IGNORE_THESE[s]){	// Then put these into the background set.
		Int4 i=IGNORE_THESE[s];
		bg_sap[i]->next=FINAL_SAP[J]; 
		bg_sap[i]=bg_sap[i]->next;
	    } else {
		gsap->next=FINAL_SAP[J]; 
		gsap=gsap->next;
	    }
#else
		gsap->next=FINAL_SAP[J]; 
		gsap=gsap->next;
#endif
	} gsap->next=0;
	for(Int4 i=1; i <= num_bgd_files; i++) bg_sap[i]->next=0;
	free(bg_sap);
   } else {			// single best sap...
	dh_type	dH=dheap(num_seqs+3,4);;
	for(sq=1; sq <= num_seqs; sq++) insrtHeap(sq,(keytyp)-BEST_SCORE[sq],dH);
	while((sq=delminHeap(dH)) != 0){
	   if(BEST_SAP[sq]){ 
		tail->next=BEST_SAP[sq]; tail=tail->next; tail->next=0; 
	   }
	} Nildheap(dH); 
   }
	for(s=1; s <= num_sap; s++) free(NUM_RPTS[s]); free(NUM_RPTS);
	free(BEST_SAP); free(BEST_SCORE);
	return head;
}

#define USAGE "usage: cma_gblastpgp query dbsfile [options]\n\
      NOTE: put concatenated cma files in query.cma\n\
      options:\n\
        -C<int>   Cutoff score for eliminating overlaps (-P option)(default: 30)\n\
        -D        Don't SEG query sequence\n\
        -e<real>  E-value cutoff for showing scores and alignments\n\
        -H<int>   heapsize (default 100000)\n\
        -I<int>:<int> - flanking region lengths for sequence hits(default: 0:0)\n\
        -L        do low complexity masking of database sequences\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -O=<str>  cma files for families omitted from output cma file (sample str: \"3,5-7,9,11-17\")\n\
        -o        output sequence hits in separate files for each input cma file\n\
        -P        Don't use PatchSAP() option to combine gapped HSPs into one\n\
        -s        output subsequence instead of full sequences\n\
        -T<int>   blast word hit threshold (default 11)\n\
        -t<real>  trim alignments at ends where marg. prob. is > <real>\n\
        -U        Use *.chk files instead of cma files\n\
        -w<int>   alignment width (default 60)\n\
        -X<int>   X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
        -z        see posMatrix\n\
\n"

#define MAX_NUM_CMA_FILES 100
int	cma_gblastpgp(int argc, char *argv[])
// Int4	cma_gblastpgp(Int4 argc, char *argv[])
{
	Int4    T=11,gap_open=11,gap_extend=1;
	time_t	time1=time(NULL);
	Int4	left_flank=0,right_flank=0;
        double	x_parameter=25.0;
	Int4	width=60,x;
	UInt4	hpsz=300000,printopt=0;
	double	Ethresh=0.001,Ecutoff=0.001;
	BooLean	see_smatrix=FALSE,use_chk_files_only=FALSE;
	BooLean	seg_seq=TRUE,success=FALSE,use_patch_sap=TRUE;
	BooLean	output_by_class=FALSE;
	BooLean output_subseqs=FALSE;
	BooLean print_names_only=FALSE;
	Int4	*IGNORE_THESE=0,num_bgd_files=0;
	char	mask_dbs=' ',*checkin=0,*checkout=0;
	FILE	*seqfile=NULL;
	char	*cmafile=0;
	double	misalncut=2.0;
	Int4	num_cma_files,overlap_cutoff=30,start_arg=3;
	UInt4	max_length_input_seqs=250000;

#if 0
	if(argc < 4) print_error(USAGE);
	start_arg=4;
	// a_type		A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	if(sscanf(argv[2],"%d",&num_cma_files) != 1) print_error(USAGE);
	if(num_cma_files < 1 || num_cma_files >= MAX_NUM_CMA_FILES) print_error(USAGE);
#else
	if(argc < 3) print_error(USAGE);
	start_arg=3;
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	FILE    *fp=open_file(argv[1],".cma","r");
	cma_typ *IN_CMA=MultiReadCMSA(fp,&num_cma_files,A);
	fclose(fp);
	if(num_cma_files < 1 || num_cma_files >= MAX_NUM_CMA_FILES) print_error(USAGE);
#endif
	for(Int4 arg=0; arg < argc; arg++) { fprintf(stderr,"%s ",argv[arg]); }
	fprintf(stderr,"\n\n");
	fprintf(stderr,"num_cma_files = %d\n",num_cma_files);
	NEW(IGNORE_THESE,num_cma_files+3,Int4);
        for(Int4 arg = start_arg; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(USAGE);
           switch(argv[arg][1]) {
	     case 'C': overlap_cutoff=IntOption(argv[arg],'C',1,100,USAGE); break;
             case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,USAGE); break;
	     case 'H': hpsz=IntOption(argv[arg],'H',1,1000000,USAGE); break;
             case 'h': Ethresh=RealOption(argv[arg],'h',0.0,10000,USAGE); break;
	     case 'I': // putsegs=TRUE;
		if(sscanf(argv[arg],"-I%d:%d",&left_flank,&right_flank) != 2)
                        print_error(USAGE); break;
             case 'L': mask_dbs='x'; break;
	     case 'm': printopt=IntOption(argv[arg],'m',0,6,USAGE); break;
	     case 'O':  {
		num_bgd_files++;
		if(argv[arg][2] != '=') print_error(USAGE);
		char *str=argv[arg]+3;
		Int4 *values,j;
		NEW(values,num_cma_files+20,Int4);
		// WARNING: may get a segmentation error if values > num_cma_files!!!
		Int4 n=ParseIntegers(str, values, "cma_gblastpgp: -O option input error");
		if(n > num_cma_files) print_error(USAGE);
                for(Int4 i=1; i<= n; i++){
                    if((j=values[i]) > 0 && j <= num_cma_files){
			if(IGNORE_THESE[j] > 0) print_error("-O option error: overlapping input sets\n");
			IGNORE_THESE[j]=num_bgd_files; fprintf(stderr,"ignore file %d\n",j);
                    } else print_error(USAGE);
                } free(values);
	       } break;
             case 'P': use_patch_sap=FALSE; argv[arg][1] = ' '; break;
             case 's': output_subseqs=TRUE; argv[arg][1] = ' '; break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,USAGE); break;
	     case 't': misalncut=RealOption(argv[arg],'t',0.1,1.0,USAGE); break;
             case 'o': output_by_class= TRUE; argv[arg][1] = ' '; break;
             case 'U': use_chk_files_only= TRUE; argv[arg][1] = ' '; break;
	     case 'w': width=IntOption(argv[arg],'w',5,50000,USAGE); break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,USAGE); break;
             case 'Z': print_names_only=TRUE; break;
             case 'z': see_smatrix= TRUE; break;
             default: print_error(USAGE);
           }
        }

	e_type	qE,QueryE[MAX_NUM_CMA_FILES],queryE=0,trueE,fakeE;
	ss_type	data,Data[MAX_NUM_CMA_FILES];
	sap_typ	sap,SAP[MAX_NUM_CMA_FILES];
	cma_typ	cma,CMA[MAX_NUM_CMA_FILES];
	char 	str[100],*Checkpoint[MAX_NUM_CMA_FILES];
	Int4	II,maxrounds;

  for(II=0; II <= num_cma_files; II++){
	Checkpoint[II]=0; QueryE[II]=0; Data[II]=0; SAP[II]=0; CMA[II]=0;
  }
  Data[0]=MakeSeqSet(argv[start_arg-1],max_length_input_seqs,A);
  queryE=ReadSeqFA(argv[1],0,A);
  QueryE[0]=CopySeq(queryE);
  for(II=1; II <= num_cma_files; II++){
#if 0
	sprintf(str,"%s_%d.cma",argv[1],II);
	CMA[II]=ReadCMSA2(str,A); 
#else
	CMA[II]=IN_CMA[II];
#endif
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
	QueryE[II]=CopySeq(queryE);
	if(!IdentSeqs(trueE,fakeE) || !IdentSeqs(trueE,queryE) ||
		LenSeq(queryE) != LengthCMSA(1,cma) || nBlksCMSA(cma) != 1){
		   fprintf(stderr,"file %s_%d.cma: ",argv[1],II);
		   print_error("incompatible with query sequence");
	}
	sprintf(str,"%s_%d.chk",argv[1],II);
	Checkpoint[II] = AllocString(str);
  }
  for(Int4 s=1; s <= num_cma_files; s++){ 
     if(IGNORE_THESE[s]) fprintf(stderr,"-------------------------------- %d: %s (%d).\n",
		s,NameCMSA(CMA[s]),IGNORE_THESE[s]);
     else fprintf(stderr,"================================ %d: %s.\n", s,NameCMSA(CMA[s])); 
  } if(print_names_only) exit(0);
  BooLean	fake_cma_srch;
  if(use_chk_files_only) fake_cma_srch=FALSE; else fake_cma_srch=TRUE;
  for(II=1; II <= num_cma_files; ){
  // start main loop over all *cma files....
	// TESTING FakeSearchCMSA( );
	if(fake_cma_srch){
	    	sprintf(str,"%s_%d.cma",argv[1],II);
		checkin=0; checkout=Checkpoint[II]; 
		cma=CMA[II]; data=Data[II]; qE=QueryE[II]; 
		qE=CopySeq(queryE);
		maxrounds=2;
	} else {
		cma=0; data=Data[0]; qE=QueryE[0]; 
		qE=CopySeq(queryE);
		checkout=0; checkin=Checkpoint[II];
		maxrounds=1;
	}
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,A);
	Int4 **mtx=0;
	if(checkin && use_patch_sap) mtx=gpsi.CopyOFposMatrix();
        do { 
	   sap=0;
	   if(fake_cma_srch){
	   	sap = gpsi.FakeSearchCMSA(T,cma); 
		fake_cma_srch=FALSE;
	   } else {
		if(misalncut < 1.0) sap = gpsi.TrimmedSearch(T,misalncut,mask_dbs);
	   	else sap = gpsi.Search(T,mask_dbs);
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
	   if(use_patch_sap){ // fix SAP here...to patch together HSPs

		sap_typ	head=MakeSelfAlignedGSAP(qE);

            	head->next = gpsi.StealSAP();
		SAP[II] = PatchSAPs(overlap_cutoff,mtx,head,NSeqsSeqSet(Data[0]),queryE,A); 
		II++;
// fprintf(stderr,"DEBUG 1\n");
		for(Int4 s=0; s < LenSeq(qE); s++) free(mtx[s]); free(mtx);
// fprintf(stderr,"DEBUG 2\n");
	   } else {
#if 1		// test this procedure...
		SAP[II]=MakeSelfAlignedGSAP(qE);
        	SAP[II]->next=gpsi.StealSAP(); II++;
#else
        	SAP[II]=gpsi.StealSAP(); II++;
#endif
	   }
  	   if(!use_chk_files_only) fake_cma_srch=TRUE;
	} 
   } // end loop over cma files
   FILE *fptr=0;
   // merge sap lists deleting weakest redundant hits. 
   for(II=1; II <= num_cma_files; II++){
	sap=SAP[II];
	if(sap){
	  fptr=stdout;
#if 0
	  printf("================== File %d: %s =======================\n",
			II,NameCMSA(CMA[II]));
			// II,NameCMSA(CMA[II]),NumSeqsListGSAP(sap));
			// II,NameCMSA(CMA[II]),NumHSPsListGSAP(sap));
#endif
	  if(printopt == 4) PutMultiGSeqAlign(stderr, sap, width, A);
	  else if(printopt == 1) PutGSeqAlignList(stderr, sap->next, width, A);
	} // else printf("================== File %d(%s): no hits!\n",II,NameCMSA(CMA[II]));
   }
   if(fptr){
	sap_typ *saps;
	UInt4 *NUM_EACH_CMA=0;
	NEW(NUM_EACH_CMA,num_cma_files+3,UInt4);
       if(output_by_class){
          saps=MergeSAP(argv[1],SAP,num_cma_files,NSeqsSeqSet(Data[0]),
			queryE,use_patch_sap,A,output_subseqs,NUM_EACH_CMA,num_bgd_files,IGNORE_THESE);
       } else {
          saps=MergeSAP(0,SAP,num_cma_files,NSeqsSeqSet(Data[0]),queryE,
			use_patch_sap,A,output_subseqs,NUM_EACH_CMA,num_bgd_files,IGNORE_THESE);
       }
       for(Int4 s=1; s <= num_cma_files; s++){ 
	   if(NUM_EACH_CMA[s] > 0){
		if(IGNORE_THESE[s]){
	  	   printf("-------------------------------- %d: %s (%d)[%d].\n",
			 	s,NameCMSA(CMA[s]),NUM_EACH_CMA[s],IGNORE_THESE[s]);
		} else {
	  	   printf("================================ %d: %s (%d).\n",
			 	s,NameCMSA(CMA[s]),NUM_EACH_CMA[s]);
		}
	   } else printf(" %d: %s no hits.\n",
			 	s,NameCMSA(CMA[s]));
       }
       // NOTE: MergeSAP() will eventually destroy all saps not returned...

       fptr=open_file(argv[1],"_Main.cma","w"); 
       xSeqAlignToCMA(fptr,saps[0],left_flank,right_flank,A);
       // tSeqAlignToCMA(fptr,saps[0],left_flank,right_flank,A);
       // gpsi.PutSeqsAsCMA(stdout,left_flank,right_flank); 
       fclose(fptr); 

       for(Int4 bg_file=1; bg_file <= num_bgd_files; bg_file++){
        if(saps[bg_file] && saps[bg_file]->next){
           sprintf(str,"%s_BG%d",argv[1],bg_file);
           fptr=open_file(str,".cma","w");
           xSeqAlignToCMA(fptr,saps[bg_file],left_flank,right_flank,A);
           fclose(fptr);
        }
       }

       // FreeGSeqAlignList(sap);

#if 0
	{
         Int4 Nset;
         sprintf(str,"%s.cma",argv[1]);
         cma_typ cma=ReadCMSA2(str,A);

         sprintf(str,"%s.purge",argv[1]);
         FILE *fp = open_file(str,".cma","w");
         PutRepSetCMSA(fp,purge,&Nset,cma); fclose(fp);

         sprintf(str,"%s.purge.cma",argv[1]);
         cma_typ cma2=ReadCMSA2(str,A);
         sprintf(str,"%s.purge",argv[1]);
         WriteMtfCMSA(str, cma2, NULL);
         TotalNilCMSA(cma2);
#endif
   } else fprintf(stderr,"No hits!\n\n");
   for(II=0; II <= num_cma_files; II++){
	// if(CMA[II]) TotalNilCMSA(CMA[II]); // merging ruins the seqset!!!!!!!!!
	if(Checkpoint[II]) free(Checkpoint[II]);
	if(QueryE[II]) NilSeq(QueryE[II]); 
   } NilSeqSet(Data[0]); NilAlpha(A);
   double  runtime=difftime(time(NULL),time1);
   fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
   if(success) return 0;
   else return 1;
}


