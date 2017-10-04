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
	
#if 0
void	chn_typ::ReadChnAnalMultAln(Int4 argc, char *argv[],const char *USAGE)
{ ReadChnAnalMultAln(argc, argv[],0,0,USAGE); }
#endif

void	chn_typ::ReadChnAnalMultAln(Int4 argc, char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,const char *USAGE)
//****************** Read in CHAIN analysis multiple alignments **********************
//****************** Read in CHAIN analysis multiple alignments **********************
//****************** Read in CHAIN analysis multiple alignments **********************
{
	Int4    i,j,s,m,n,r;
        FILE    *fp=0;
        if(!OutputAln && extra_files && !StdAlignmentOnly){
          // fp = open_file(argv[1],".cmd","w");
          fp = open_file(output_name,".ccmd","w");
          for(i = 0; i < argc; i++) { if(argv[i][1] != ' ') fprintf(fp,"%s ",argv[i]); }
          fprintf(fp,"\n"); fclose(fp); fp=0;
        }
        // fprintf(stderr,"cbp_cut=%g\n",cbp_cut);
        cbp_cut = log10(cbp_cut);
        // if(FractHighlight > 1.0){ FractHighlight=0.50; }
        double  total_units;
        switch(PageSetUp){
                case 'L': total_units=22650; break;
                case 'P': total_units=14000; break;
                case 'l': total_units=14000; break;
                case 'p': total_units=10400; break;
                default: print_error(USAGE); break;
        }
        char_per_line=(Int4) floor((total_units/(double)(fontsize*6)));
        char_per_line = char_per_line - 21;  // subtract (16+5)from each side of alignment

        // 1. gblastpgp keyfa against orthogs to get *.oma file.

	Number=0; fp=0;
	if(cmc_rtf_mode){
	  if(verbose) fprintf(stderr,"input CMC files passed into chn_typ constructor.\n");
	} else if(cmc_input){
	  // fprintf(stderr,"input CMA files passed into chn_typ constructor.\n");
	  // fp=open_file(argv[1],".mma","r");
	} else if(cha_as_input){
	  fprintf(stderr,"Reading input file (%s.cha).\n",argv[1]);
	  fp=open_file(argv[1],".cha","r");
	} else if(StdAlignmentOnly){
	  fprintf(stderr,"Reading input file (%s.cma).\n",argv[1]);
	  fp=open_file(argv[1],".cma","r");
	} else {
	  fprintf(stderr,"Reading input file (%s.chn).\n",argv[1]);
	  fp=open_file(argv[1],".chn","r");
	}
	// Retrieve information regarging selected columns.
	char	**Status=0;
	if(fp != 0){   // if fp opened above...
		assert(CMA_ARRAY == 0); 
		IN_CMA=MultiReadCMSA(fp,&Number,&Status,AB); 
		OwnCMAs=TRUE; fclose(fp); fp=0;
	} else {	// borrow cma from above...
		assert(CMA_ARRAY); NEW(IN_CMA,NumCMA + 5,cma_typ); Number=NumCMA;
		if(cmc_rtf_mode) Status=cmc_rtf_status; else Status=0;
		OwnCMAs=FALSE;
		for(Int4 j=1; j <= NumCMA; j++) IN_CMA[j]=CMA_ARRAY[j];
	} if(fp) fclose(fp);
#if 1	// New: get rid of Quarkjunk script...
	if(StdAlignmentOnly){
		if(Number != 1) print_error("chn_see input error for -S option");
		cma_typ tmp_cma=IN_CMA[1]; free(IN_CMA); 
		NEW(IN_CMA,4,cma_typ);
		IN_CMA[1]=tmp_cma; IN_CMA[2]=CopyCMSA(tmp_cma);
		char *tmp_status = 0;
#if 0
		tmp_status=Status[1]; free(Status);
#else	// the following needs more work.
		if(Status){
		     tmp_status=Status[1]; free(Status);
		} else {
		     assert(nBlksCMSA(tmp_cma) == 1);
		     NEW(tmp_status, LengthCMSA(1,tmp_cma) +5, char);
		     for(Int4 j=0; j < LengthCMSA(1,tmp_cma); j++) tmp_status[j]='*';
		}
#endif
		NEWP(Status,4,char);
		Status[1] = tmp_status;
		Status[2] = AllocString(tmp_status);
		Number=2;
	}
#endif
#if 1
	if(InputFractionAsNumber > 0){
		FractHighlight=(double) InputFractionAsNumber/(double) LengthCMSA(1,IN_CMA[1]);
                fprintf(stderr,"################ FractHighlight = %f ################\n",
                                        FractHighlight);
                InputFraction=TRUE;
	}
#endif
	if(cmc_input){
		NEWP(Status, Number + 3, char);
	} else for(i=1; i <= Number; i++){
		if(strstr(Status[i]+1,"!")){
		  if(extra_files) fprintf(stdout,"%d: %s\n",i,Status[i]+1); 
		} else if(i == 2 && NumResidues > 1){  // multiple selected residues...
		  Int4 k;
	          for(j=1; j <= NumResidues; j++){
	            if(SelectColPos){
			k = Position[j];
		    } else {
	            	k = RealToFakeCMSA(1, Position[j], IN_CMA[1]);
		    }
	            if(k < 1 || k > LengthCMSA(1,IN_CMA[1]))
				print_error("-P option input error 1");
	            fprintf(stderr,"Real = %d; Fake = %d\n",Position[j],k);
		    Status[i][k]='!';
		  }
		} else { free(Status[i]); Status[i]=0; }
	} if(extra_files && verbose) fprintf(stdout,"(%d cma files)\n",Number);
	// NOTE: NEED TO FREE UP THESE STRINGS LATER!!!
	if(Number < 2) { 
	   // if(Number==1)PutCMSA(stdout,IN_CMA[1]);
	   fprintf(stderr,"%d chn cmafiles read\n",Number); 
	   print_error("chn file input error");
	}
	cma_typ	cma=IN_CMA[1]; 
	if(NumSeqsCMSA(cma) < Seq4numbering){
		fprintf(stderr,"You chose seq #%d for numbering, ",Seq4numbering);
		fprintf(stderr,"but the input file contains only %d sequences.\n",
			NumSeqsCMSA(cma));
		print_error("FATAL ERROR: exiting...\n\n");
		// print_error(USAGE_START); break;
	}
	full_csq=CopySeq(TrueSeqCMSA(1,IN_CMA[2]));
	// 2. Get Consensus sequence from the ortholog alignment.
#if 0	// OLD METHOD...
	// keyE=CopySeq(SeqSetE(1,TrueDataCMSA(cma)));
	keyE=CopySeq(SeqSetE(1,DataCMSA(cma))); // This eliminates overhangs...
#elif 0	// NEWER METHOD THAT ACCOMODATES TRIMMING
	Int4 start = TruePosCMSA(1,1,cma);
	Int4 end = start+LengthCMSA(1,cma)-1;
	if(LenSeq(SeqSetE(1,TrueDataCMSA(cma))) != LengthCMSA(1,cma)){
		keyE=CopySeq(FakeSeqCMSA(1,cma)); // This eliminates overhangs...???
		// keyE=MkSubSeq(start,end,SeqSetE(1,DataCMSA(cma)));
		// get sub of FakeSeq...
		if(Seq4numbering == 0) Seq4numbering = 1;
		// print_error("Input error: first sequence must not contain gaps");
	} else {
		keyE=MkSubSeq(start,end,SeqSetE(1,TrueDataCMSA(cma)));
	}
#else	// afn:10-30-2008 NEWEST METHOD THAT ACCOMODATES BOTH TRIMMING AND GAPS
	Int4 start = TruePosCMSA(1,1,cma);	// ignores offset...
	Int4 end = start+LengthCMSA(1,cma)-1;	// will jump over gaps...
#if 0
	if(start < 1){ 			// accomodate gaps on front end
		e_type tmpE=SeqSetE(1,TrueDataCMSA(cma));
		SetOffSetSeq(1,tmpE);
		tmpE=FakeSeqCMSA(1,cma);
		SetOffSetSeq(1,tmpE);
	}
#elif 0
	if(start < 1){ 			// accomodate gaps on front end
	   if(LenSeq(SeqSetE(1,TrueDataCMSA(cma))) != LengthCMSA(1,cma)){
		keyE=CopySeq(FakeSeqCMSA(1,cma)); // This eliminates overhangs...???
		// keyE=MkSubSeq(start,end,SeqSetE(1,DataCMSA(cma)));
		// get sub of FakeSeq...
		// if(Seq4numbering == 0) Seq4numbering = 1;
		if(start < 1) SetOffSetSeq(1,keyE);
	   } else {
		// print_error("Input error: first sequence must not contain gaps");
		keyE = CopySeq(SeqSetE(1,DataCMSA(cma)));
	   }
	} else 
#endif
        {	// Make sure that there are no deletions at ends of alignment
	     e_type fakeE = SeqSetE(1,DataCMSA(cma));
	     if(start < 1 || IsDeletedCMSA(1, start, cma)){
		print_error("FATAL: query sequence appears to start with a deletion"); 
	     }
	     if(IsDeletedCMSA(1, end, cma)){
		print_error("FATAL: query sequence appears to end with a deletion"); 
	     }
      
	     keyE=MkSubSeq(start,end,SeqSetE(1,DataCMSA(cma)));  // Fake sequence...
	// ^ This jumps over extensions at end of block.
	}
/**********************************************************************************
fprintf(stderr,"DEBUG1: start=%d; end = %d\n",start,end);

 Quarkjunk (chn_see with -S option) is exiting here for the following input (which has 
the first seq. 209==209 but where start == 0 due to gap): 

[0_(1)=RfcLike(1){go=19,gx=2,pn=5.0,lf=0,rf=0}:
(209)**************************************************************************!*!*!*************!!*******************!******!****************!********!****!*!********!!!**!***********!****!!****!!***********!*****

$1=209(209):
>gi|67470376|ref|XP_651156.1|  {|6(100)|<Entamoebidae(E)>}Activator 1 40 kDa subunit [Entamoeba histolytica HM-1:IMSS]
{()-IPWVEKYRPKLLDEIIGNVDIIKTLKSFRDSKQFPHLLLCGQPGIGKTTSIHCLAHELLKDRYKDAVLELNASDERGIETIRSTIKSFCEKKLVLPDNlPKIVILDEADSMTTAAFQALRRTMEIHSKTTRFVLACNTPEKIIEPIQSRCARLTFRPLGEEELMNRIKEIAHCEGVDIEDDAVKALEIVSEGDMRKAINALQTCAIIQG()}*

_0].

error message:
  MkSubSeq(): start=0; end = 208
  Seq: SubSeq out of range

Need to fix this...

 **********************************************************************************/
	assert(LenSeq(keyE) == LengthCMSA(1,cma));
	// fprintf(stderr,"start = %d; end = %d\n",start,end); 
	// PutSeq(stderr,keyE,AB);
	// correcting for offset...
	if(Begin==0){	// then use default settings
	  Begin=1; End = LenSeq(keyE);
	} else if(Seq4numbering > 0){
	  if(!UseColPos){
		Begin = RealToFakeCMSA(Seq4numbering,Begin,cma);
	  	Int4 input_end=End;
	  	End = RealToFakeCMSA(Seq4numbering,End,cma);
	  	Int4 true_end=TruePosCMSA(Seq4numbering,LengthCMSA(1,cma), cma);
	  	// find true end of sequence.
	  	fprintf(stderr,"End = %d; true_end=%d; input_end=%d\n", End,true_end,input_end);
	  	if(End == 0) print_error("Fatal: Input beyond end!");
	  }
	} else {
	  if(!UseColPos){
	      fprintf(stderr,"End = %d; Begin=%d\n", End,Begin);
	      Begin = RealToFakeCMSA(1,Begin,cma); End = RealToFakeCMSA(1,End,cma); 
	      if(End == 0){
	  	fprintf(stderr,"End = %d; Begin=%d; offset = %d; fakeLen=%d\n", End,Begin,
			OffSetSeq(TrueSeqCMSA(1,cma)),LenSeq(FakeSeqCMSA(1,cma)));
		PutSeq(stderr,TrueSeqCMSA(1,cma),AB);
		PutSeq(stderr,FakeSeqCMSA(1,cma),AB);
		gsq_typ	*gsq=gsqCMSA(1,cma);
		gsq->Put(stderr,AB);
		print_error("Fatal: Input either within an insertion or beyond end of sequence!");
	      }
	  }
	}
#if 1	// AFN: 8/11/05: try using this option always...
	if(Seq4numbering == 0) Seq4numbering = 1;
#endif
	if(Begin < 1){ 
	  fprintf(stderr,"Begin = %d; End = %d\n",Begin,End);
	  print_error("-B option out of permissible range...aborting...");
	}
	Begin--;		// start at zero not one.
#endif
	// But may throw off residue numbering.
	// WARNING: move keyE to first position in the PSI-BLAST alignment!!
	qE=MkConsensusCMSA(cma);
	// fp=open_file(argv[1],".csq","w"); PutSeq(fp,qE,AB); fclose(fp);
	// move this to chn_aln procedure.
	if(extra_files && !StdAlignmentOnly){
	  // fp=open_file(argv[1],".seq","w"); 
	  fp=open_file(argv[1],".dsq","w"); 
	  PutSeqSetEs(fp,TrueDataCMSA(cma)); fclose(fp);
	}

	NEW(CMA,2*Number+4,cma_typ);
	MCMA[0]=CMA[0]=CMA[1]=cma;
	assert(nBlksCMSA(CMA[1]) == 1); // just one block for now...
    	assert(cbp_cut < 0.0);
	NumSeqAln=NumSeqsCMSA(cma);
	if(Number > 2){
	  NEW(TAX_CMA,NumSeqAln+3,cma_typ);
	  for(i=3,j=1; j <= NumSeqAln; j++,i++){ TAX_CMA[j]=IN_CMA[i]; }
	} else TAX_CMA=0;
	status = 0;
	// fprintf(stderr,"Checking integrity of input alignment.\n");
	if(Number > 2) { 
	   NumberMCMA = Number - NumSeqAln - 1;
	   NumAnalysis=1; m=1;
	   if(extra_files && verbose) fprintf(stderr,"%d chn cmafiles read\n",Number);
	   if(Number < (NumSeqAln+2)){ // should be one alignment for each sequence 
		fprintf(stderr,"Number = %d < NumSeqAln+2=%d\n",Number,(NumSeqAln+2));
		print_error("inconsistent number of subfamily alignment files");
	   } else if(Number >= (NumSeqAln+2)){  // NumSeqAln + 2 == 1 background model.
		NEWP(status,2*Number+4,char);
		if(Number >= (NumSeqAln+MAX_ALN_CHN_TYP)){
		   fprintf(stderr,"Number = %d; NumSeqAln=%d\n",Number,NumSeqAln);
		   fprintf(stderr,"Maximum alignments = %d < %d - %d = %d\n",
			MAX_ALN_CHN_TYP-1,Number,NumSeqAln,Number-NumSeqAln);
		   print_error("Too many input alignments");
		} else {	// extra hierarchical alignments
		  NumAnalysis++;
		  CMA[NumAnalysis]=MCMA[1]=IN_CMA[2];	
		  if(Status[2]) status[NumAnalysis]=Status[2];
		  for(j=2,i=(NumSeqAln+3) ; i <= Number; i++){
	   		m++; NumAnalysis++;
			CMA[NumAnalysis]=MCMA[m]=IN_CMA[i]; 
		  	if(Status[i]) status[NumAnalysis]=Status[i];
			if(LengthCMSA(1,MCMA[m]) != LenSeq(qE)) {
				fprintf(stderr,"Query and main sets of unequal lengths!\n");
				fprintf(stderr,"Main set = %d; query = %d\n",
					LengthCMSA(1,MCMA[m]),LenSeq(qE));
				print_error("chn file input error");
			}
		  }  
		}
	   } else { CMA[2]=MCMA[1]=IN_CMA[2]; }
	   assert(m == NumberMCMA);
	} else if(Number == 2){
	   // NumberMCMA=1; NumAnalysis++; MCMA[1]=IN_CMA[2];
	   NumberMCMA=1; NumAnalysis=2; MCMA[1]=IN_CMA[2];
	} else print_error("query alignment file input error");
	if(OutputAln){ PutAln(output_name,'W'); exit(0); } 
	if(RmFamilyAln){
		FILE *ofp = open_file(argv[1],"_new.chn","w");
		PutCMSA(ofp,cma);
		assert(MCMA[2]);
		PutCMSA(ofp,MCMA[2]);
		for(j=1; j<=NumSeqAln; j++){
			PutCMSA(ofp,TAX_CMA[j]);
		}
		for(j=3; j<=NumberMCMA; j++) PutCMSA(ofp,MCMA[j]);
		fclose(ofp);
		exit(0);
	}
#if 1	// Subset option (e.g., -PW207)   
	if(NumResidues > 0){
	  assert(!cmc_input);
	  if(NumAnalysis != 2) 
		print_error("-P option cannot be used with multiple alignments");
	  assert(NumberMCMA == 1);

	  BooLean	*skipsqS,*skipsqN;
	  cma_typ	cmaS,cmaN;
	  unsigned char	r1;
	   
	  for(i=1; i <= NumResidues; i++){
#if 1
	      if(SelectColPos){
	        Position0[i] = Position[i];
	      } else Position0[i] = RealToFakeCMSA(1, Position[i], CMA[1]);
	      if(extra_files) fprintf(stderr,"Real = %d; Fake = %d\n",Position[i],Position0[i]);
	      if(Position0[i] < 1 || Position0[i] > LenSeq(qE))
				print_error("-P option input error 1");
#else
	      e_type firstE = TrueSeqCMSA(1,CMA[1]);
	      Position0[i] = Position[i] - OffSetSeq(firstE);
	      if(Position0[i] < 1 || Position0[i] > LenSeq(qE))
				print_error("-P option input error 1");
#endif
#if 1	// NOT SURE THAT THIS IS NECESSARY...Add FORCE option later...
	      // r1 = ResidueCMSA(1,1,Position0[i],MCMA[1]); 
	      r1 = ResidueCMSA(1,1,Position0[i],CMA[1]); 
	      if(!MemSset(r1,Residues[i])){
		fprintf(stderr,"%c%d in seq 1 of query set fails to match input pattern %s%d\n",
			AlphaChar(r1,AB),Position0[i],residue_str[i],Position[i]);
		print_error("-P option input error 2");
	      }
#endif
	  }
	  NEW(skipsqS,NumSeqsCMSA(MCMA[1])+2,BooLean);
	  NEW(skipsqN,NumSeqsCMSA(MCMA[1])+2,BooLean);
	  skipsqS[1]=skipsqN[1]=FALSE;
	  //************* Find sequences matching input pattern ************
#if 0
	  for(Int4 sq=2; sq <=NumSeqsCMSA(MCMA[1]); sq++)
#else
	  for(Int4 sq=1; sq <=NumSeqsCMSA(MCMA[1]); sq++)
#endif
	  {
	    skipsqS[sq]=FALSE; skipsqN[sq]=TRUE;   // put in foreground by default..
	    Int4 mismatches=0;
	    for(i = 1; i <= NumResidues; i++){
	      r1 = ResidueCMSA(1,sq,Position0[i],MCMA[1]);
	      // If any residue fails to match the pattern put seq in the background set.
	      if(!MemSset(r1,Residues[i])){ mismatches++; }

	      // NOTE: Fix problem with removing query family sequences from background set
	      // Remove these from superfamily only; not from Non-matching set!
	      // Should eventually simply retain the original main set!!!!
	    }
	    if(mismatches > MaxMisMatches){ skipsqS[sq]=TRUE; skipsqN[sq]=FALSE; }
	  }

// fprintf(stderr,"DEBUG 1\n");
	  //************* Output cma files of sequence matching pattern ************
          fp=tmpfile(); PutSelectCMSA(fp,skipsqS,MCMA[1]); 
          rewind(fp); cmaS=ReadCMSA(fp,AB); fclose(fp);

          fp=tmpfile(); PutSelectCMSA(fp,skipsqN,MCMA[1]); 
          rewind(fp); cmaN=ReadCMSA(fp,AB); fclose(fp);
#if 1
	  // fp=open_file(output_name,"_fg.cma","w");
	  // PutSelectCMSA(fp,skipsqS,MCMA[1]); fclose(fp);

	  fp=open_file(output_name,"_bg.cma","w");
	  PutSelectCMSA(fp,skipsqN,MCMA[1]); fclose(fp);
#endif
	  free(skipsqS); free(skipsqN);

// fprintf(stderr,"DEBUG 2\n");

          TotalNilCMSA(MCMA[1]);	// used for pattern matching...
	   // if(InvertResSet){ cma_typ cmaTmp=cmaS; cmaS = cmaN; cmaN = cmaTmp; }
	  MCMA[1]=CMA[2]=cmaS;
#if 1	// output selected sequences
fp=open_file(output_name,".chn","w");
PutCMSA(fp,IN_CMA[1]); PutCMSA(fp,cmaS);
for(j=3; j <= Number; j++) PutCMSA(fp,IN_CMA[j]);
fclose(fp);
// PutSeqSetEs(fp,TrueDataCMSA(cmaS));
#endif
	  MCMA[2]=CMA[3]=cmaN;
	  IN_CMA[2] = cmaS; Number++; IN_CMA[Number] = cmaN;
	  // IN_CMA has a few extra fields to accomodate additions...
	  NumberMCMA++; NumAnalysis++; 
	  assert(isprint(output_name[0])); assert(!OutputAln);
	  PutAln(); // SAVE SEQS. in output file.
	}
#endif
	for(m=0; m <= NumberMCMA; m++){
	   char	*AlignName=NameCMSA(MCMA[m]); AlignCode[m]= AlignName[0];
	   if(LengthCMSA(1,MCMA[m]) != LenSeq(qE)){
#if 0
		fprintf(stderr,"LengthCMSA(1,MCMA[%d]) == %d != %d\n",
					m,LengthCMSA(1,MCMA[m]),LenSeq(qE));
#else
		fprintf(stderr,"Query and main sets of unequal lengths!\n");
		fprintf(stderr,"Main set = %d; query = %d\n",
					LengthCMSA(1,MCMA[m]),LenSeq(qE));
#endif
		print_error("chn file input error");
	   }
	}
	free(Status);
}

