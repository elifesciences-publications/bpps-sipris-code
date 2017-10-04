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

e_type	**gsm_typ::search(Int4 inso, Int4 insx, Int4 del, BooLean segmask, char *msafile)
{
	e_type	**ListE;
	sma_typ MA;
	Int4	j,n,total,*mtfcnts=NULL;
	double	*Freq;

	fprintf(stderr,"creating profile...\n");
	// if(bestmap <= 1. && !input_msa) print_error("failed to find a significant motif");
	// input_msa=FALSE;
	if(counts) free(counts); if(nsize) free(nsize);
        if(DBS_hits) NilSet(DBS_hits);
	number = GetFastaInfo(Argv[2], MAX_IN_SEQS, &counts, &nsize, AB);
	for(maxseqleng=0,j=1; j<= number; j++) 
		maxseqleng = MAXIMUM(Int4,nsize[j],maxseqleng);
	DBS_hits = MakeSet(number+1);
	// oldcardB = cardB = 0;
	for(total=0, j=0; j<=nAlpha(AB); j++) total += counts[j];
        NEW(Freq,nAlpha(AB)+2,double);
        if(aafreq == 'm') {
          NEW(mtfcnts, nAlpha(AB)+2, Int4);
          MA=ReadSMA(Argv[2]); CountsSMA(mtfcnts, MA); NilSMA(MA);
          for(n=0, j=0; j<=nAlpha(AB); j++) n += mtfcnts[j];
          for(j=0; j<=nAlpha(AB); j++) Freq[j]= (double)mtfcnts[j]/(double)n;
          free(mtfcnts);
        } else if(aafreq == 'd') {
          for(j=0; j<=nAlpha(AB); j++) Freq[j] = (double)counts[j]/(double)total;
        } else { 	// aafreq = 's'
          for(j=0; j<= nAlpha(AB); j++) Freq[j] = blosum62freq[j];
        }

	FILE	*fptr = open_file(Argv[2],"","r");
	gsn_typ	F;
	F=MakeGOScan(msafile,AB,maxrpts,method,mode,ecut,Ecut,maxseqleng,
		weight,minmap,pseudo,Freq);
#if 1
	// double tmp_d = (double) TotLenPrtnModel(PrtnModelGOScan(F)) * 0.90;
	MinSeqLen = TotLenPrtnModel(PrtnModelGOScan(F));
#endif
        if(!segmask) NoMaskGOScan(F);
	if(mask_nonglobular) MaskNonGlobularGOScan(F);
	SetRptEvalGOScan(repeatEval,F);
	fprintf(stderr,"scanning database...\n");
	FILE *ofp=open_file(name,".gscn","w");
#if 1
	ListE = GapSeqGOScanScan(ofp,fptr,inso,insx,left_flank,right_flank,
					min_rpt,number,total,nsize,F);
#else	// code below is a work in progress; want to speed up purge by starting with a cma. afn: 7_17_12
#if 0
Int4    GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
        Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total,
        unsigned short *nsize, char ***Operation, e_type **RtnE, e_type **FullE,
        unsigned short **FullR, Int4 **Start, char Mode, gsn_typ F,
        UInt4 minlen, UInt4 maxlen)
        Int4 NumHits=GapSeqGOScanScan(fptr,a,b,gapo,gapx,left,right,min_rpt,number,
                total,nsize,&operation,&ListE,&FullE,&FullR,&start,Test,F,minlen,maxlen);
#endif
	// Taken from GOScanToCMSA() within goscan.cc
	// need a,b,operation,&ListE,&FullE,&FullR,&start,Test,minlen,maxlen.
	Int4		a=18,b=2;
        char            **operation;
        e_type          *ListE,*FullE;
        unsigned short  *FullR;
        Int4            *start;
	cma_typ         cma=0;
	Int4		minlen=0,maxlen=INT4_MAX;
	ptm_typ         PM=PrtnModelGOScan(F);
	Int4            nblks=NumModelsPrtnModel(PM);
	BooLean         *skip;

        Int4 NumHits=GapSeqGOScanScan(fptr,a,b,inso,insx,left,right,min_rpt,number,total,
		nsize,&operation,&ListE,&FullE,&FullR,&start,Test,F,minlen,maxlen);
        if(NumHits > 0){
          cma=MakeCMSA(ListE,NumHits,operation,start,nblks,PrtnModelLengths(PM),
                        gapo,gapx,pernats,left,right,argv[1],PrtnModelA(PM),FullE,FullR);
                        // 10000,2000,1000,left,right,argv[1],PrtnModelA(PM),FullE,FullR);
          ss_type FullSeq = FullSeqCMSA(cma);
          if(!AddFull) RmFullCountsCMSA(cma);
          if(verbose) { PutAlnCMSA(stdout,cma); }
	} // WARNING: NOT DEALLOCATING ALL MEMORY WHEN NumHits == 0! Fix here later...
        free(operation); free(start);
	// free(counts); free(nsize); 
#endif
	fclose(ofp); NilGOScan(F); fclose(fptr); free(Freq);
	return ListE;
}

void    gsm_typ::gapped_srch()
// NEW: AFN 7/8/09.
// output 
{
	Int4 del = 0;
	BooLean segmask=TRUE;
	if(trim_info > 0.0){
          sprintf(str,"%s.trim.msa",name);
	} else {
          sprintf(str,"%s.msa",name);
	}
	e_type	**ListE=search(Gap_o, Gap_x, del, segmask, str);
	FILE *sfp=open_file(name,".grpts","w");
        for(Int4 s=1; ListE[s]; s++){
	   Int4 r,NumRpts=0;
           for(r=1; ListE[s][r]; r++){ }
	   r--;
	   if(r >= min_rpt){
             for(r=1; ListE[s][r]; r++){ 
		e_type tE=ListE[s][r];
		if(LenSeq(tE) >= MinSeqLen) PutSeq(sfp,tE,AB); 
		NilSeq(tE); 
	     }
           } free(ListE[s]);
        } free(ListE);
	fclose(sfp);
}

void gsm_typ::search( )
{
	Int4	j,n,total;
	FILE	*sfp,*fptr,*fp2;
	gsn_typ	F;
	snh_typ	sH;
	h_type	HG;
	sma_typ MA;
	Int4	*mtfcnts=NULL;

	fprintf(stderr,"creating profile...\n");
	if(bestmap <= 1. && !input_msa) print_error("failed to find a significant motif");
	input_msa=FALSE;
	if(number == -1){ /*** run GetFastaInfo( ) on first iteration ***/
	   number = GetFastaInfo(Argv[2], MAX_IN_SEQS, &counts, &nsize, AB);
	   for(maxseqleng=0,j=1; j<= number; j++) 
		maxseqleng = MAXIMUM(Int4,nsize[j],maxseqleng);
	   DBS_hits = MakeSet(number+1);
	   oldcardB = cardB = 0;
	}
	fptr = open_file(Argv[2],"","r");
	// Read cma file and create trimmed cma and msa files...
	fprintf(stderr,"***************** trim_info = %.3f ****************\n",trim_info);
	if(trim_info > 0.0){
	   sprintf(str,"%s.cma",name);
           cma_typ tmp_cma=ReadCMSA2(str,AB);
           Int4 *RmLeft,*RmRight,*TrimLimit;
           NEW(RmLeft,nBlksCMSA(tmp_cma)+3,Int4);
           NEW(RmRight,nBlksCMSA(tmp_cma)+3,Int4);
           NEW(TrimLimit,nBlksCMSA(tmp_cma)+3,Int4);
           cma_typ cma2 = TrimCMSA(trim_info,TrimLimit,RmLeft,RmRight,tmp_cma);
           sprintf(str,"%s.trim",name); PutAlnCMSA(str,cma2,NULL);
           NilCMSA(cma2);
           free(RmLeft); free(RmRight); free(TrimLimit);
           TotalNilCMSA(tmp_cma);
           sprintf(str,"%s.trim.msa",name);
        } else {
           sprintf(str,"%s.msa",name);
	} 

	for(total=0, j=0; j<=nAlpha(AB); j++) total += counts[j];
        NEW(freq,nAlpha(AB)+2,double);
        if(aafreq == 'm') {
          NEW(mtfcnts, nAlpha(AB)+2, Int4);
          MA=ReadSMA(Argv[2]); CountsSMA(mtfcnts, MA); NilSMA(MA);
          for(n=0, j=0; j<=nAlpha(AB); j++) n += mtfcnts[j];
          for(j=0; j<=nAlpha(AB); j++) freq[j]= (double)mtfcnts[j]/(double)n;
          free(mtfcnts);
        } else if(aafreq == 'd') {
          for(j=0; j<=nAlpha(AB); j++) freq[j] = (double)counts[j]/(double)total;
        } else { 	// aafreq = 's'
          for(j=0; j<= nAlpha(AB); j++) freq[j] = blosum62freq[j];
        }

	F=MakeGOScan(str,AB,maxrpts,method,mode,ecut,Ecut,maxseqleng,
		weight,minmap,pseudo,freq);
	if(mask_nonglobular) MaskNonGlobularGOScan(F);
	SetRptEvalGOScan(repeatEval,F);
	fprintf(stderr,"scanning database...\n");
	ptm_typ PM=PrtnModelGOScan(F);
	MaxBadBlks=NumBlksPrtnModel(PM)-MinGoodBlks;
	if(MaxBadBlks < 0) MaxBadBlks=0;
	sH=GOScanScan(fptr,number, total, nsize,1,MaxBadBlks,F);
	fprintf(stderr,"output results...\n");
        if(nScanHeap(sH) > 0){
	     fp2 = open_file(name,".scn","w");
	     fprintf(fp2,"goscan %s %s.msa -c -e%f -E%f\n",
		Argv[2],name,ecut,Ecut);

             HG=HistScanHeap(sH);
             PutHist(fp2,60,HG);
             PutInfoScanHeap(fp2,sH,AB);
             for(j=1; j<=nblksScanHeap(sH); j++){
                  HG=histScanHeap(j,sH); PutHist(fp2,60,HG);
             } fclose(fp2);
	     if(report_gaps){
	     	fp2 = open_file(name,".gaps","w");
	     	PutGapsScanHeap(fp2, sH); fclose(fp2);
	     }
	     sfp = open_file(name,".seq","w");
             PutSeqScanHeap(sfp, sH, AB);
	     AddSetScanHeap(sH, DBS_hits, TRUE);
	     cardB = CardSet(DBS_hits);
	     fclose(sfp); 
	     if(NumFailSeqScanHeap(sH) > 0){
	     	sfp = open_file(name,".fsq","w");
             	PutFailSeqScanHeap(sfp, sH, AB);
	     	fclose(sfp); 
	     }
    	     if(repeats){ 
	        sfp=open_file(name,".rpts","w");
		// PutRptsScanHeap(sfp,left_flank,right_flank, sH, AB); // OLD
		if(min_rpt < 2) PutRptsScanHeap(sfp,left_flank,right_flank,sH,AB);
		else PutMinRptsScanHeap(sfp,left_flank,right_flank,min_rpt,sH,AB);
		fclose(sfp);
	     }
        }
	fclose(fptr); 

#if 1	//*********************** Gapped repeat search *****************************
	if(do_gapped_srch) gapped_srch();
#endif

	fprintf(stdout,"total # of sequences detected: %d\n", cardB);
	if(cardB > oldcardB){ oldcardB = cardB; }
	else if(!force) maxrun = Run;	// i.e., quit 

	if(Run < maxrun){  
	  if(cma_purge > 0){ // Then add step to output scanned sequences as *cma file.
                assert(maxrpts == 1);
	        Run++;
	        cma_typ	tmpcma1,tmpcma2;
                cma_typ tmpcma = ScanHeap2CMSA(sH,0,name,AB,FALSE);
                sH = NULL;          /** WARNING: ScanHeap2CMSA( ) destroys sH ***/
                // if(tmpcma!= NULL){ WriteMtfCMSA(name,tmpcma,NULL); }

	// Purge cma file 
                Int4 Nset;
	        sprintf(name,"%s%d",Argv[1],Run);
                sprintf(str,"%s.purge",name);
                FILE *fp = open_file(str,".cma","w");
                PutRepSetCMSA(fp,cma_purge,&Nset,tmpcma); fclose(fp);

                sprintf(str,"%s.purge.cma",name);
                tmpcma1=ReadCMSA2(str,AB);

	        if(tmpcma != NULL) TotalNilCMSA(tmpcma);
	        if(sH) NilScanHeap(sH); NilGOScan(F); free(freq);
		tmpcma=cmsa; cmsa=0; 
	        tmpcma2 = optbreed(tmpcma1);
		if(tmpcma2) TotalNilCMSA(tmpcma2);
		cmsa=tmpcma; // This gets freed in ~gsm_typ().
         } else {
	  NilScanHeap(sH); 
	  NilGOScan(F); 
	  free(freq);
	  fprintf(stderr,"purging sequence set...\n");
	  if(do_gapped_srch) sprintf(str,"%s.grpts",name);
	  else if(repeats) sprintf(str,"%s.rpts",name);
	  else sprintf(str,"%s.seq",name);
	  sprintf(name,"%s%d",Argv[1],Run);
	  if(UseRepset){  // then first generate representative set...
	        // fp2 = open_file(name,".rep","w");
	        fp2 = open_file(name,"","w");
		// n= RepSetCluster(fp2, str, 0, !noseg, 'b', cutoff, AB);
		// n=RepSetCluster(fp2,str,0,!noseg,'b',cutoff,UseLabeled,AB);
		n=RepSetCluster(fp2,str,0,!noseg,'S',cutoff,UseLabeled,AB);
		// this should call gapped blast with score option...
		fclose(fp2);
	        // sprintf(str,"%s.rep",name);
	  } n =PurgeFiles(str,name,cutoff,inc,minseq,maxseq,AB);
	  if(n < minseq) go = FALSE;
	 }
	} else { NilScanHeap(sH); NilGOScan(F); free(freq); }
}

