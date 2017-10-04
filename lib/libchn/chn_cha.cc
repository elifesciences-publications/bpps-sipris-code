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

void    chn_typ::PutCHA_RTNS(FILE *fptr,unsigned short numRtns)
{
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3");
        for(unsigned short i=1; i <=numRtns; i++) fprintf(fptr,"\n\\par");
        fprintf(fptr," }\n");
}

void	chn_typ::PutCHA_TAIL(FILE *fptr)
{ fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1  {\\f3\n\\par}}\n"); }

void	chn_typ::PutCHA_HEADER(FILE *fptr)
{
    fprintf(fptr,"{\\rtf1\\ansi \\deff4\\deflang1033");
    fprintf(fptr,"{\\fonttbl{\\f0\\froman\\fcharset0\\fprq2 Tms Rmn;}\n");
    fprintf(fptr,"{\\f1\\fswiss\\fcharset0\\fprq2 Arial;}\n");
    fprintf(fptr,"{\\f2\\fswiss\\fcharset0\\fprq2 Arial Narrow;}\n");
    fprintf(fptr,"{\\f3\\fmodern\\fcharset0\\fprq1 Courier New;}\n");
    fprintf(fptr,"{\\f4\\froman\\fcharset0\\fprq2 Times New Roman;}\n");
    fprintf(fptr,"{\\f5\\fswiss\\fcharset0\\fprq2 System;}\n");
    fprintf(fptr,"{\\f6\\fmodern\\fcharset0\\fprq1 Courier New;}}\n");
    fprintf(fptr,"{\\colortbl;\\red0\\green0\\blue0;\\red0\\green0\\blue255;");
    fprintf(fptr,"\\red0\\green255\\blue255;\\red0\\green255\\blue0;");
    fprintf(fptr,"\\red230\\green5\\blue190;\\red255\\green0\\blue0;");
    fprintf(fptr,"\\red255\\green255\\blue0;\\red255\\green255\\blue255;");
    fprintf(fptr,"\\red0\\green0\\blue128;\\red0\\green128\\blue128;");
    fprintf(fptr,"\\red0\\green128\\blue0;\\red128\\green0\\blue128;");
    fprintf(fptr,"\\red128\\green0\\blue0;\\red128\\green128\\blue0;");
    /** 15 dk grey; 16 lt grey (original) **/
    fprintf(fptr,"\\red128\\green128\\blue128;\\red192\\green192\\blue192;");
    /** 17 pink; 18 orange; **/
#if 0	// two new colors: dark orange; dark cyan
    fprintf(fptr,"\\red255\\green102\\blue255;\\red225\\green180\\blue0;}");
#else
    fprintf(fptr,"\\red255\\green102\\blue255;\\red225\\green180\\blue0;");
    fprintf(fptr,"\\red225\\green110\\blue10;\\red57\\green154\\blue181;}");
#endif
    fprintf(fptr,"{\\stylesheet{\\widctlpar \\f4\\fs%d \\snext0 Normal;}",fontsize);
    fprintf(fptr,"{\\*\\cs10 \\additive Default Paragraph Font;}}\n");
#if 1	// Test header
{
    fprintf(fptr,"{\\header \\pard\\plain \\qr \\widctlpar \\b\\f2\\fs%d \n",fontsize+2);
    time_t time0=time(0);
    struct tm *tp;
    if(tp=localtime(&time0)){
       char min_str[10];
       if(tp->tm_min < 10) sprintf(min_str,":0%d",tp->tm_min);
       else sprintf(min_str,":%d",tp->tm_min);
       fprintf(fptr,"{%s (%d%s - %d/%d/%d) page }\n", output_name,
                tp->tm_hour,min_str,tp->tm_mon+1,tp->tm_mday,tp->tm_year+1900);
    } else fprintf(fptr,"{ %s page }\n",output_name); // NEW: filename at top.
    // fprintf(fptr,"{\\field{\\*\\fldinst { PAGE }}}{ \\par }}\n",fontsize+2);
    fprintf(fptr,"{\\field{\\*\\fldinst { PAGE }}}{ \\par }}\n");
}
#endif

    switch(PageSetUp){
	  case 'l': // 8.5" x 11" landscape.
	   fprintf(fptr,"\\paperw15840\\paperh12240\\margl720\\margr720");
    	   fprintf(fptr,"\\margt1080\\margb1080 ");
    	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   // fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\headery360\\endnhere \n"); // header high up
    	   // page_width=12000.0;
	   break;
	  case 'L': // 11" x 17" landscape
	   fprintf(fptr,"\\paperw24480\\paperh15840\\margl720\\margr720");
	   fprintf(fptr,"\\margt1080\\margb1080 ");
	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   // fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\headery360\\endnhere \n");
	   // print_error("PutCHA_HEADER( ) 11\" x 17\" landscape not yet implemented");
	   break;
	  case 'p': // 8.5" x 11" portrait.
	   fprintf(fptr,"\\paperw12240\\paperh15840\\margl720\\margr720");
    	   fprintf(fptr,"\\margt1080\\margb1080 ");
    	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
           // fprintf(fptr,"\\linex0\\endnhere\\sectdefaultcl \n");
           fprintf(fptr,"\\linex0\\headery360\\endnhere\\sectdefaultcl \n");
	   // page_width=7000.0;
	   break;
	  case 'P': // 11" x 17" portrait.
	   fprintf(fptr,"\\paperw15840\\paperh24480\\margl720\\margr720");
           fprintf(fptr,"\\margt1080\\margb1080 ");
           fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
           // fprintf(fptr,"\\linex0\\endnhere\\sectdefaultcl \n");
           fprintf(fptr,"\\linex0\\headery360\\endnhere\\sectdefaultcl \n");
            // page_width=12000.0;
	   break;
	  default: print_error("PutCHA_HEADER( ) PageSetUp input error");
    }
    fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs%d \n",fontsize);
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
#if 1
    if(StdAlignmentOnly) fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 %s }{\\f3 \n\\par }\n",
		fontsize+4,NameCMSA(IN_CMA[1])); // top.
    else fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 %s }{\\f3 \n\\par }\n",
		fontsize+4,output_name); // NEW: filename at top.
#else
    time_t time0=time(0);
    struct tm *tp;
    if(tp=localtime(&time0)){
       fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 %s (%d:%d - %d/%d/%d) }{\\f3 \n\\par }\n",
		fontsize+2,output_name,
                tp->tm_hour,tp->tm_min,tp->tm_mon+1,tp->tm_mday,tp->tm_year+1900);
    } else 
    fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 %s }{\\f3 \n\\par }\n",
		fontsize+2,output_name); // NEW: filename at top.
#endif
}

char    chn_typ::FractionToCharPBC(double fract)
{
    const char rtn[]="0123456789!!!";
    Int4 f = (Int4) floor(10.0*fract);
    if(f >= 0 && f <= 10) return rtn[f]; else return '?';
}

void	chn_typ::PutCHA(FILE *fptr, const char *ColorFont)
// Create an alignment figure for publication:
// Use this routine when have only one gapped 'block'.
{
	Int4	color,colorB,i,j,n,N;
	char	r,c,str[50],state,laststate,*ss_struct=0;
	char	*color_font=0;

    assert(ntypSMA(MA)==1);
// PutClustalwSMA(stdout, 1, MA); exit(1);
    char	*gnull=gnullSMA(1,MA);
    Int4	next_start,start;
    Int4	gstart=0,gend;
    Int4	*sq_index[MAX_ALN_CHN_TYP_PLUS],Analysis,analysis;
    Int4	*pos_key_seq=0;
    char	**conserved=0;
    // Int4	fontsize=rtf[1]->FontSize();
    // assert(fontsize > 4 && fontsize < 50);
    // Int4	gapsize=fontsize-4;
    Int4	gapsize=fontsize;

    Int4	anal_index,AnalysisList[MAX_ALN_CHN_TYP_PLUS],NumAnalIndex=NumAnalysis;
    assert(NumAnalysis < MAX_ALN_CHN_TYP);
    if(BG_CMA){ 
      AnalysisList[0]=0; anal_index=NumAnalysis+1; 
      for(Analysis=1; Analysis <= NumAnalysis+1; Analysis++){
	   AnalysisList[Analysis] = anal_index; anal_index--;
      } NumAnalIndex=NumAnalysis+2; AnalysisList[NumAnalIndex]=1;
      assert(Analysis==NumAnalysis+2);
    } else {
      AnalysisList[0]=0; anal_index=NumAnalysis; 
      for(Analysis=1; Analysis <= NumAnalysis; Analysis++){
	   AnalysisList[Analysis] = anal_index; anal_index--;
      } NumAnalIndex=NumAnalysis+1; AnalysisList[NumAnalIndex]=1;
    }
    //**************** Classify residue positions for 3-tier hierarchies ****************
    //**************** Classify residue positions for 3-tier hierarchies ****************
    //**************** Classify residue positions for 3-tier hierarchies ****************
    if(NumAnalysis!=3 && status){	// 'status' is Generated by BPPS procedure...
     for(Int4 na=1; na <= NumAnalysis; na++){
      if(status[na]){
	Int4 NumCols = rtf[na]->FixSetStatus(status[na],rtfQ[na]);
	double FractHighlite;
	if(PatternOnlyMode) FractHighlite = 1.01*(double)NumCols/(double) rtf[na]->Length();
	else FractHighlite = 1.1*(double)NumCols/(double) rtf[na]->Length();
	if(FractHighlite > 1.0) FractHighlite = 1.0;
	if(FractHighlite <= 0.0) FractHighlite = 0.10;
	if(InputFraction) FractHighlite=FractHighlight;
	rtf[na]->ReInit(FractHighlite,-cbp_cut);
	if(extra_files) fprintf(stdout,"fixed status = %s\n",status[na]+1);
      }
     }
    } else 
    if(NumAnalysis==3 && status && status[2]){	// 'status' is Generated by BPPS procedure...
      Int4 NumCols = rtf[2]->FixSetStatus(status[2],rtfQ[2]);
      double FractHighlite;
      if(PatternOnlyMode) FractHighlite = 1.01*(double)NumCols/(double) rtf[2]->Length();
      else FractHighlite = 1.1*(double)NumCols/(double) rtf[2]->Length();
      if(FractHighlite > 1.0) FractHighlite = 1.0;
      if(FractHighlite <= 0.0) FractHighlite = 0.10;
      if(InputFraction) FractHighlite=FractHighlight;
      rtf[2]->ReInit(FractHighlite,-cbp_cut);
      if(extra_files) fprintf(stdout,"fixed status = %s\n",status[2]+1);
      for(anal_index=0; anal_index <= NumAnalIndex; anal_index++){
	   Analysis=AnalysisList[anal_index];
	   if(status[2] && status[Analysis] == 0){
	    if(Analysis == 3){ // main set.
	        // use ave col_prob for '!' in status[2][j] to find
	        //  cutoff for main set histogram.
	        // Get mainset status using columns in 
	        // status[Analysis] = rtf[Analysis]->DeriveStatus(status[2],rtfQ[Analysis],'M');
	        status[Analysis] = rtf[Analysis]->DeriveStatus(status[2],rtf[Analysis],'M');
		if(status[Analysis] !=0){
		   // fprintf(stderr,"DEBUG: status = %s\n",status[Analysis]+1);
		   NumCols = rtf[Analysis]->FixSetStatus(status[Analysis],rtfQ[Analysis]);
#if 0
		   fprintf(stderr,"Fixed status = %s\n",status[Analysis]+1);
		   fprintf(stderr,"NumCols = %d\n",NumCols);
#endif
      		   if(PatternOnlyMode) FractHighlite = 
				1.01*(double)NumCols/(double) rtf[Analysis]->Length();
      		   else FractHighlite = 1.1*(double)NumCols/(double) rtf[Analysis]->Length();
      		   if(FractHighlite > 1.0) FractHighlite = 1.0;
      		   if(FractHighlite <= 0.0) FractHighlite = 0.10;
      		   if(InputFraction) FractHighlite=FractHighlight;
      		   rtf[Analysis]->ReInit(FractHighlite,-cbp_cut);
		}
	    } else if(Analysis == 1){ // Family set.
	        status[Analysis] = rtf[Analysis]->DeriveStatus(status[2],rtfQ[Analysis],'F');
		if(status[Analysis] !=0){
		   NumCols = rtf[Analysis]->FixSetStatus(status[Analysis],rtfQ[Analysis]);
      		   FractHighlite = 1.5*(double)NumCols/(double) rtf[Analysis]->Length();
      		   if(FractHighlite > 1.0) FractHighlite = 1.0;
      		   if(FractHighlite <= 0.0) FractHighlite = 0.10;
      		   if(InputFraction) FractHighlite=FractHighlight;
      		   rtf[Analysis]->ReInit(FractHighlite,-cbp_cut);
		}
	    }
	   }
      }
      if(status[4]){ free(status[4]); status[4]=0; }
      if(BG_CMA && use_bg_cma){	// then determine status based on Main Set (Analysis=3).
	        // status[4] = rtf[4]->DeriveStatus(status[3],rtfQ[4],'M');
	        status[4] = rtf[4]->DeriveStatus(status[3],rtf[4],'M');
#if 0
	        for(i=1; i <= rtf[0]->Length(); i++){
	           if(status[4][i] == '?') { status[4][i]='!'; }
	       	}
	fprintf(stderr,"Raw status = %s\n",status[4]+1);
#endif
		// fprintf(stderr,"status = %s\n",status[4]+1);
		NumCols = rtf[4]->FixSetStatus(status[4],rtfQ[4]);
#if 0
	fprintf(stderr,"Fixed status = %s\n",status[4]+1);
	fprintf(stderr,"NumCols = %d\n",NumCols);
#endif
      		if(PatternOnlyMode) FractHighlite = 
				1.01*(double)NumCols/(double) rtf[4]->Length();
      		else FractHighlite = 1.1*(double)NumCols/(double) rtf[4]->Length();
      		if(FractHighlite > 1.0) FractHighlite = 1.0;
      		if(FractHighlite <= 0.0) FractHighlite = 0.10;
      		if(InputFraction) FractHighlite=FractHighlight;
      		rtf[4]->ReInit(FractHighlite,-cbp_cut);
      }
      if(status[0] == 0){	// Intermediate Categories...
	char	*tmp_status;
        Int4    i,n;
        NEW(tmp_status,rtf[0]->Length()+3,char);
        for(i=1; i <= rtf[0]->Length(); i++){
	   tmp_status[i]='*';
           for(n=1; n<=4; n++){
                if(status[n] && status[n][i] != '*') { tmp_status[i]='!'; break; }
	   }
        }
        if(extra_files) fprintf(stdout,"tmp status = %s\n",tmp_status+1);
	status[0] = rtf[0]->DeriveStatus(tmp_status,rtfQ[0],'I');
        if(extra_files) fprintf(stdout,"intermediate status = %s\n",status[0]+1);
	free(tmp_status);
      }
    }

    //***************************************************************************
    NEWP(conserved,lengthSMA(1,MA)+3,char);
    // Loop for fancy formating of page...
#if 1	// NEW: 10_22_04: AFN.
    if(pos_key_seq ==0 && Seq4numbering > 0){
	if(nseqSMA(MA) < Seq4numbering) print_error("Sequence for numbering doesn't exist");
        NEW(pos_key_seq,glengthSMA(1,MA)+3,Int4); // residue number at each position.
	Int4 current_pos=startSMA(1,Seq4numbering,MA); // block 1, sequence 'Seq4numbering'.
        char *gseq=gseqSMA(1,Seq4numbering,MA); 
// fprintf(stderr,"seq %d(%d): %s\n",Seq4numbering,current_pos,gseq); exit(0);
	if(!isalpha(gseq[0])) current_pos++; // gap at start --> true start is 1 ahead;
	for(i=0; i < glengthSMA(1,MA); i++){
	   if(isalpha(gseq[i])){ pos_key_seq[i]=current_pos; current_pos++; }
	   else pos_key_seq[i]=current_pos;
        }
    }
#endif
  if(maxlen_gnull < INT4_MAX){

    gnull_insrt_len = gnullInsrtLenSMA(1,MA);

#if 0
    for(j=0; j < glengthSMA(1,MA); j++){
	fprintf(stderr,"gnull[%d] = %c %d\n",j,gnull[j],gnull_insrt_len[j]);
    }
#endif
  } else gnull_insrt_len=0;
    Begin = MINIMUM(Int4,Begin,glengthSMA(1,MA));
    next_start=Begin; 
    for(gstart=0,j=0; j < Begin; gstart++){ 
	if(gnull[gstart] != '_') j++; 
	else if(gstart > glengthSMA(1,MA)) return; // range outside of seq.
    }
#if 1	// FIX PROBLEM WITH GAPS AT START...WITH -B OPTION
    while(gnull[gstart] == '_'){ gstart++; }
#endif
    for(anal_index=0; anal_index <= NumAnalIndex; anal_index++){
        NEW(sq_index[anal_index],nseqSMA(MA)+2,Int4);
        for(n=1; n<= nseqSMA(MA); n++) {
	   sq_index[anal_index][n]=startSMA(1,n,MA);
           char *gseq= gseqSMA(1,n,MA);
#if 1	// Repair problem with gaps at beginning...
	   // If gap at start then true start is one position ahead of startSMA();
	   if(!isalpha(gseq[0])) sq_index[anal_index][n]++; // 
#endif
#if 1	// FIX PROBLEM WITH GAPS AT START...WITH -B OPTION
	   // for(i=0,j=0; j < (Begin + leng_start_gap); i++)
	   for(i=0; i<gstart; i++) if(isalpha(gseq[i])){ sq_index[anal_index][n]++; }
#else
	   for(i=0,j=0; j < Begin; i++){
		if(gnull[i] != '_') j++;
		if(isalpha(gseq[i])){ sq_index[anal_index][n]++; }
	   }
#endif
	}
    }
#if 1
// fprintf(stderr,"NumAnalysis=%d; color=%s\n",NumAnalysis,ColorFont);
// use 2,3,4,..,(NumAnalysis+1)/2;
if(ColorFont == 0){
   n=(NumAnalysis+1)/2;
   switch(n){
#if 0
     case 2: color_font=AllocString(" LM   "); break;
     case 3: color_font=AllocString(" LRM   "); break;
     case 4: color_font=AllocString(" LoRM   "); break;
     case 5: color_font=AllocString(" LyoRM   "); break;
     case 6: color_font=AllocString(" LgyoRM   "); break;
     case 7: color_font=AllocString(" LngyoRM   "); break;
     case 8: color_font=AllocString(" LBngyoRM   "); break;
     case 9: color_font=AllocString(" LmBngyoRM   "); break;
     case 10: color_font=AllocString(" LDmBngyoRM   "); break;
     case 11: color_font=AllocString(" WLDmBngyoRM   "); break;
     default: color_font=AllocString(" CyPBGRcmbgr                               "); break;
#else
     case 2: color_font=AllocString(" Ly   "); break;
     case 3: color_font=AllocString(" LRy   "); break;
     case 4: color_font=AllocString(" LoRy   "); break;
     case 5: color_font=AllocString(" LMoRy   "); break;
     case 6: color_font=AllocString(" LgMoRy   "); break;
     case 7: color_font=AllocString(" LCgMoRy   "); break;
     case 8: color_font=AllocString(" LBCgMoRy   "); break;
     case 9: color_font=AllocString(" LPBCgMoRy   "); break;
     case 10: color_font=AllocString(" LDPBCgMoRy   "); break;
     case 11: color_font=AllocString(" WLDPBCgMoRy   "); break;
     default: color_font=AllocString(" CyPBGRcmbgr                               "); break;
#endif
   }
} else color_font=AllocString(ColorFont);
#endif
    //**************** Main loop for each 'page' of alignment *******************
    //**************** Main loop for each 'page' of alignment *******************
    //**************** Main loop for each 'page' of alignment *******************
    Int4 loop=0;
    Int4 lines_on_page=0;
    Int4 loops_on_page=0;
    char last_gnull=0;
    for( ; gstart < glengthSMA(1,MA) && next_start < End; gstart=gend,loop++){
     gend=MINIMUM(Int4,glengthSMA(1,MA),gstart+char_per_line);
#if 1	// Fix Kannan's bug: afn: 8/21/10.
     start=next_start;
     { 
        Int4 ii,jj;
        for(ii=gstart,jj=start; ii < gend && jj < End; ii++){ if(gnull[ii] != '_'){ jj++; } }
	next_start=jj;
	// NOTE: if there is a deletion at the start of file
	// then jj index may be incorrect as that is not counted as a residue...
	// can check for '-' characters using MA array.
     }
#endif
#if 1	// add page breaks.... only for tabloid portrait and fontsize '8' (default) for now.
     // Need to add a routine with a table for various page setups and fontsizes.
     // Need exact accounting of lengths left on page given size minus margins...
     if(PageSetUp == 'P'){
       if((fontsize <= 16 && fontsize >= 12) && loops_on_page > 0 && lines_on_page > 0){
	Int4 total_lines;
	if(fontsize == 16) total_lines = 160 - loops_on_page*10;
	else if(fontsize == 14) total_lines = 185 - loops_on_page*10;
	else if(fontsize == 12) total_lines = 215 - loops_on_page*10;
	else total_lines =100;	// this setting should not be inputable.
	Int4 lines_left = total_lines-lines_on_page;
	Int4 lines_in_loop=lines_on_page/loops_on_page;
	if(verbose){
	  fprintf(stderr,"*************** lines_on_page=%3d *****************\n",lines_on_page);
	  fprintf(stderr,"*************** lines_in_loop=%3d *****************\n",lines_in_loop);
	  fprintf(stderr,"***************** lines_left=%3d ******************\n",lines_left);
	}
	if(lines_left <= lines_in_loop){
	  fprintf(fptr,"\n\\par \\page \n");
     	  lines_on_page=0; loops_on_page=0;
	}
       }
     } else if(PageSetUp == 'L'){
       if(fontsize == 12 && loops_on_page > 0 && lines_on_page > 0){
	Int4 total_lines;
	if(fontsize == 12) total_lines = 180;
	else total_lines =100;	// this setting should not be inputable.
	Int4 lines_left = total_lines-lines_on_page;
	Int4 lines_in_loop=lines_on_page/loops_on_page;
	if(verbose){
	  fprintf(stderr,"*************** lines_on_page=%3d *****************\n",lines_on_page);
	  fprintf(stderr,"*************** lines_in_loop=%3d *****************\n",lines_in_loop);
	  fprintf(stderr,"***************** lines_left=%3d ******************\n",lines_left);
	}
	if(lines_left <= lines_in_loop){
	  fprintf(fptr,"\n\\par \\page \n");
     	  lines_on_page=0; loops_on_page=0;
	}
       }
     } loops_on_page++;
#endif
     // loop for various types of alignments....
     // 0=Query SubFamily A alignment with standard background frequencies.
     // 1=Family A alignment with Family B background frequencies.
     // 2=Family B alignment with Family C background frequences.
     // 3=Family C alignment with standard background frequences.
     //**************** Loop over each level of the hierarchy *******************
     //**************** Loop over each level of the hierarchy *******************
     //**************** Loop over each level of the hierarchy *******************
     char *AlignName=0;
     for(anal_index=0; anal_index <= NumAnalIndex; anal_index++){
	// if(StdAlignmentOnly) if(anal_index >=1) break;
	Analysis=AnalysisList[anal_index];
	analysis=BackGrnd[Analysis];
#if 0	// This doesn't appear to be working; possible read/write errors!!!!
	if(NumAnalysis == 3 && Analysis == 4){
		if(start > EndBG_CMA) continue; // skip motif alignment if no motifs.
		if(next_start < StartBG_CMA) continue; // skip motif alignment if no motifs.
		// NOTE: next_start will have been calculated for the top alignment?
		// Need to also deal with regions between motifs...
	}
#endif
	double *fract_seqs=FractSeqAln[Analysis];;
	if(fract_seqs && Analysis > 0) {
	   if(MCMA[Analysis-1]){
		AlignName=NameCMSA(MCMA[Analysis-1]); // These are off by one...
	   	if(strcmp(AlignName,"BackGround")==0) continue; // antiquated convention...
	   	// if(AlignName[0] == 'B') continue;	// Old code...
	   } else AlignName=0; 
	   color=get_font_color_code(color_font[ColorCode[Analysis]]);
	   if(use_sfbg && Analysis == 1){	// 'U' option.
		colorB = get_font_color_code(color_font[ColorCode[AnalysisList[1]]]);
	   } else if(Analysis == NumAnalysis) colorB = get_font_color_code('L'); 
	   else if(ColorCode[Analysis+1] == 0){ colorB=color; }
	   else colorB = get_font_color_code(color_font[ColorCode[Analysis+1]]);
	} else{
	   AlignName=0;
	   color=get_font_color_code(color_font[ColorCode[Analysis]]);
	   if(use_sfbg && Analysis == 1){	// 'U' option.
		colorB = get_font_color_code(color_font[ColorCode[AnalysisList[1]]]);
	   } else if(Analysis == NumAnalysis) colorB = get_font_color_code('L'); 
	   // if(Analysis == NumAnalysis) colorB = get_font_color_code('L'); 
	   else if(Analysis == 0) colorB = get_font_color_code('L'); 
	   else colorB = get_font_color_code(color_font[ColorCode[Analysis+1]]);
	}

	for(i=gstart,j=start; i < gend && j < End; i++){
	  if(gnull[i] != '_'){
		j++;
		if(conserved[j]) free(conserved[j]);
		double *res_evals=rtf[Analysis]->ResEvalues(j);
#if 0	// Revert back to initial method, as AlignNames are inconsistent.
#if 0
	if(Analysis==0) conserved[j]=rtf[analysis]->Conserved(res_evals);
	else switch(AlignCode[Analysis-1])
	{	// Note: MargProb[anal] actually corresponds to anal-1;
	   case 'C': conserved[j]=rtfQ[analysis]->Conserved(res_evals); break;
	   case 'S': case 'F': case 'Q': 
	   default: conserved[j]=rtf[analysis]->Conserved(res_evals); break;
	}
#else	// AFN 12/22/05 - New...
	if(Analysis==0) conserved[j]=rtf[analysis]->Conserved(j);
	else switch(AlignCode[Analysis-1])
	{	// Note: MargProb[anal] actually corresponds to anal-1;
	   case 'C': conserved[j]=rtfQ[analysis]->Conserved(j); break;
	   case 'S': case 'F': case 'Q': 
	   default: conserved[j]=rtf[analysis]->Conserved(j); break;
	}
#endif
#else
  		// conserved[j]=rtfQ[analysis]->Conserved(res_evals); 
  		// conserved[j]=rtf[analysis]->Conserved(res_evals); 
  		conserved[j]=rtf[analysis]->Conserved(j); // AFN: 7/27/07 (fixes bug).
#endif
	  }
	}
        double	**wtfreq=Freq[Analysis];	// background freq (= WtFrq[] )
	ss_struct=SecondStrct[Analysis];
	Int4	tab1=15*fontsize*6;
	Int4	tab2=tab1+(fontsize*6);
	Int4	tab3;
	// NEW: find end tab site.
	for(i=gstart,j=start; i < gend && j < End; i++){ if(gnull[i] != '_') j++; }
	tab3=tab2+((i-gstart+1)*fontsize*6)+(fontsize*6);

	if(Analysis==0){
// fprintf(stderr,"**************** DEBUG **************** Analysis=%d \n",Analysis);
	  // next_start=j; // WARNING: assumes always call Analysis==0 first!!!
          fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
	  lines_on_page++;
          fprintf(fptr,"\\widctlpar\\tqr\\tx%d\\tx%d\\tx%d\\adjustright \\fs%d\\cgrid ",
                tab1,tab2,tab3,fontsize);
	}
	//************************ output tick marks **************************
	// This prints out before the next file, so need to run one 'dummy' run after done printing
	BooLean	PrintTickMarks=FALSE;
	if(anal_index==1) PrintTickMarks=TRUE;
#if 1	// print tick marks below all alignments: AFN 6-9-2017.
	else if(!BG_CMA && anal_index > 1 && anal_index < (NumAnalIndex-1)) PrintTickMarks=TRUE;
#endif
	else if(!BG_CMA && anal_index > NumAnalysis && NumAnalysis >= 2) PrintTickMarks=TRUE;
	// BG_CMA is a block based motif file....
	else if(BG_CMA && anal_index > (NumAnalysis+1) && NumAnalysis >= 2) PrintTickMarks=TRUE;
	else if(anal_index == 2){	// if anal_index==1 was a background alignment...
		Int4 tmpAnal=AnalysisList[1]; // look at previous
		assert((MCMA[tmpAnal -1]) != 0 && tmpAnal > 0);
		AlignName=NameCMSA(MCMA[tmpAnal -1]); // These are off by one...
		assert(AlignName);
	   	if(strcmp(AlignName,"BackGround")==0) PrintTickMarks=TRUE;; // antiquated convention...
	}
	if(PrintTickMarks){
	 if(Seq4numbering > 0){	  // checked for possible segmentation fault above...
	  assert(Seq4numbering <= nseqSMA(MA));
	  fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 position",fontsize);
	  fprintf(fptr,"\\tab  \\tab }"); 
	  Int4 site,site0=0;
	  Int4 spaces=0;
          char *gseq2=gseqSMA(1,Seq4numbering,MA); // pointer set to [0,1,...] array
	  site=sq_index[anal_index][1];
	  last_gnull=0;
	  for(spaces=0,i=gstart,j=start; i < gend && j < End; i++){
		assert(i < glengthSMA(1,MA));
		if(pos_key_seq) site=pos_key_seq[i];
		if(gnull[i] != '_'){
			j++;
			if(anal_index > NumAnalysis) sq_index[anal_index][1]++;
		}
#if 1
		if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]  && !isalpha(gseq2[i])){
		   if(gnull[i] == '_' && last_gnull != '_'){
			// spaces += 2 + (Int4) ceil(log10(gnull_insrt_len[i]));
			spaces += 2 + (Int4) ceil(log10((double)gnull_insrt_len[i]+0.00001));
		   }
		} else if(!isalpha(gseq2[i])) spaces++;
#else
		if(!isalpha(gseq2[i])) spaces++;
#endif		
		else {
		   if(site%10 == 0){	// print out position.
		   	if(spaces > 4){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 ",gapsize);
			  for(Int4 sp=5; sp <= spaces; sp++) fprintf(fptr," ");
			  fprintf(fptr,"}\n"); spaces=4;
			}
		   	if(spaces == 4){
			  assert(site < 100000);
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %5d}\n",gapsize,site);
			  spaces=0;
		   	} else if(site < 10000 && spaces == 3){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %4d}\n",gapsize,site);
			  spaces=0;
		   	} else if(site < 1000 && spaces == 2){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %3d}\n",gapsize,site);
			  spaces=0;
		   	} else if(site < 100 && spaces == 1){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %2d}\n",gapsize,site);
			  spaces=0;
		   	} else if(site < 10 && spaces == 0){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %d}\n",gapsize,site);
			  spaces=0;
			} else spaces++;
		   } else if(site%5 == 0){
		   	if(spaces > 0){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 ",gapsize);
			  for(Int4 sp=1; sp <= spaces; sp++) fprintf(fptr," ");
			  fprintf(fptr,"}\n");
			  spaces=0;
			} fprintf(fptr,"{\\f3\\fs%d\\cf1 .}\n",gapsize);
		   } else spaces++;
	   	   site++;
		} last_gnull=gnull[i];
	  }
	  if(spaces > 0){
		fprintf(fptr,"{\\f3\\fs%d\\cf1 ",gapsize);
		for(Int4 sp=1; sp <= spaces; sp++) fprintf(fptr," "); fprintf(fptr,"}\n");
	  } fprintf(fptr,"\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
#if 1	// Get rid of this!!!!!!!!!!!
	 } else {	// NOTE: eventually merge this with the above!!!!
	  if(anal_index==1) fprintf(fptr,"{\\b\\f2\\fs%d\\cf1 position",fontsize);
	  // else fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d position",fontsize,color);
	  else fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d SubFamily",fontsize,color);
	  fprintf(fptr,"\\tab  \\tab }"); 
	  Int4 site,site0=0;
	  last_gnull=0;
	  for(site=sq_index[anal_index][1],i=gstart,j=start; i < gend && j < End; i++){
		if(gnull[i] == '_'){
#if 1	// hide null option.
		  if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]){
		   if(i > gstart && last_gnull != '_'){
			 fprintf(fptr,"{\\f3\\fs%d\\cf1    }\n",gapsize);
		   }
		  } else {
#endif
		   if(i > gstart) fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		   // don't increment at start because true start site is a position ahead.
		  }
		} else {
		   j++; site0++;
		   if(site%10 >= 6){	// don't print anything yet.
		   } else if(site%10 == 0){
		   	if(site0 < 5){
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 ",gapsize);
			  for(Int4 sp=1; sp <= site0; sp++) fprintf(fptr," ");
			  fprintf(fptr,"}\n");
			} else {
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 %5d}\n",gapsize,site);
			}
		   } else if(site%5 == 0) fprintf(fptr,"{\\f3\\fs%d\\cf1 .}\n",gapsize);
		   else fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		   if(anal_index > NumAnalysis) sq_index[anal_index][1]++;
		   site++;
		}
		last_gnull=gnull[i];
	  } fprintf(fptr,"\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
#endif
	 }
	}
	//**************************** StdAlignment *********************************
	if(StdAlignmentOnly){
	  if(anal_index >=1){ 
	    for(i=gstart,j=start; i < gend && j < End; i++){
		if(gnull[i] != '_'){
		   j++; 
		   sq_index[anal_index][1]++;
		}
	    } break;
	  }
	}
	if(!BG_CMA && anal_index > NumAnalysis){ continue; }
	if(BG_CMA && anal_index > NumAnalysis+1){ continue; }
	if(!SuperAln[Analysis] && anal_index > 0) continue;	// AFN 5_7_10; omite last alignment from output.
	//****************** output histogram at top of figure **********************
	fprintf(fptr,"{\\f3 \n\\par }"); lines_on_page++;
   if(!verbose){  // then put histogram of superalignment above subalignment.
	// if(rtf[Analysis]->Printable() && Hist[Analysis]) 
// NOTE: I may want to delete this option and revert back to original...
#if 0
	if(rtfQ[Analysis]->Printable()){
	    rtfQ[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,
					tab3,maxlen_gnull,gnull_insrt_len);
	}
#else
	// this prints out the histogram above each alignment.
	if(rtf[Analysis]->Printable()){
#if 1	// for cross conserved values...
	    if(Analysis==0 && Xconserved){
		double *value=Xconserved; 
		double *pval=0; NEW(pval,rtf[Analysis]->Length()+3,double);
		double *old_value=rtf[Analysis]->SwapValue(value,pval);
	        rtf[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,
					tab3,maxlen_gnull,gnull_insrt_len);
		value=rtf[Analysis]->SwapValue(old_value,pval); 
		if(status && status[Analysis]==0){
		   NEW(status[Analysis],rtf[Analysis]->Length()+3,char);
		   for(Int4 z=1; z <= rtf[Analysis]->Length(); z++){
			if(value[z] > 0) status[Analysis][z]='!';
			else status[Analysis][z]='*';
		   }
		} free(pval);
	    } else rtf[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,
					tab3,maxlen_gnull,gnull_insrt_len);
#else
	    rtf[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,
					tab3,maxlen_gnull,gnull_insrt_len);
#endif
	}
#endif
   } else {	  // otherwise put sub-alignment histogram.
	if(rtfQ[Analysis] && rtfQ[Analysis]->Printable() && Hist[Analysis]){
	    fprintf(fptr,"\n{\\f3 \n\\par }");
	    lines_on_page++;
#if 0
	    rtfQ[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,tab3);
#else
	    rtfQ[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,tab3,
				maxlen_gnull,gnull_insrt_len);
#endif
	}
   }
	//*************** Print out dots for BPPS selected columns ******************
      if(NumAnalysis!=3){ // getting this will trigger the next conditional statement...
	//*************** Print out dots for BPPS selected columns ******************
	if(status && status[Analysis]){	// BPPS generated hierarchical alignment.
	  // fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Subfamily",fontsize,color);
	  char family_str[200];
#if 1
	  // AlignName=NameCMSA(MCMA[Analysis-1]); // These are off by one...
	  if(Analysis > 0){
	    // fprintf(stderr,"%d.%s\n",Analysis,status[Analysis]+1);
	    strcpy(family_str,NameCMSA(MCMA[Analysis-1]));
	  } else 
#endif
	  switch(Analysis){
	   // case 0: strcpy(family_str,"Intermediate"); break;
	   // case 0: strcpy(family_str,"Cross-Conserved"); break;
	   case 0: strcpy(family_str,"CrossConserved"); break;
	   case 1: strcpy(family_str,"Family"); break;
	   case 2: strcpy(family_str,"Superfamily"); break;
	   case 3: strcpy(family_str,"Main set"); break;
	   case 4: strcpy(family_str,"motifs"); break;
	   default: strcpy(family_str,"Unknown"); break;
	  }
	  fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s",fontsize,1,family_str); // black
	  fprintf(fptr,"\\tab  \\tab }"); 
	  Int4 site,site0=0;
	  last_gnull=0;
	  for(site=sq_index[anal_index][1],i=gstart,j=start; i < gend && j < End; i++){
		if(gnull[i] == '_'){
		  if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]){
		   if(last_gnull != '_'){
#if 0
			  fprintf(fptr,"{\\f3\\fs%d\\cf1    }\n",gapsize);   // always 3 spaces...
#else		// fix problem with -hide option and dots '*' = '3f below.
			  Int4 spaces = 2 + (Int4) ceil(log10((double) gnull_insrt_len[i]+0.00001));
			  fprintf(fptr,"{\\f3\\fs%d\\cf1 ",gapsize);
			  for(Int4 sp=1; sp <= spaces; sp++) fprintf(fptr," ");
			  fprintf(fptr,"}\n");
#endif
		   }
		  } else {
			fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		  }
		} else {
		   j++; site0++;
		   if(status[Analysis][j] == '!'){	// dots printed out here!
			fprintf(fptr,"{\\f3\\fs%d\\cf1 \\u9679\\'3f}\n",gapsize);
		   } else if(status[Analysis][j] == '?'){
			fprintf(fptr,"{\\f3\\fs%d\\cf%d \\u9679\\'3f}\n",
				gapsize,get_font_color_code('L'));
		   } else {	// status == '*' or "." or "^"...
			fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		   } if(anal_index > NumAnalysis) sq_index[anal_index][1]++;
		   site++;
		}
		last_gnull=gnull[i];
	  } fprintf(fptr,"\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
	}
// XXX EDITED UP TO HERE...
      } else if(NumAnalysis==3){ // this accesses the next conditional statement...
	// Print out dots for selected columns...
	if(status && status[Analysis]){	// BPPS generated hierarchical alignment.
	  // fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Subfamily",fontsize,color);
	  // fprintf(stderr,"DEBUG: fs = %d\nAnalysis=%d\nSTATUS=%s\n",fontsize,Analysis,status[Analysis]+1);
	  switch(Analysis){
	   case 0: // conventional alignment.
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Intermediate",fontsize,1); // black
	     break;
	   case 1: // family
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Query set",fontsize,1); // black
	     break;
	   case 2: // superfamily
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Foreground",fontsize,1); // black
	     break;
	   case 3: // main set
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Background",fontsize,1); // black
	     break;
	   case 4: // motifs
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d motifs",fontsize,1); // black
	     break;
	   default: // unknown alignment.
	     fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Unknown",fontsize,1); // black
	     break;
	  } 
	  fprintf(fptr,"\\tab  \\tab }"); 
	  Int4 site,site0=0;
#if 1	// crs file HERE...
	  unsigned char *qseq = seqSMA(1,1,MA); // query seq
	  char crs_color=' ';
	  // ColorString=" CyPBGRcmbgr";
	  switch(Analysis){
		case 0: crs_color='G'; break;
		case 1: crs_color='B'; break;
		case 2: crs_color='Y'; break;
		case 3: crs_color='M'; break;
	   	default: // unknown alignment.
	     	break;
	  }
#endif
	  last_gnull=0;
	  for(site=sq_index[anal_index][1],i=gstart,j=start; i < gend && j < End; i++){
		if(gnull[i] == '_'){
		  if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]){
		   if(last_gnull != '_'){
			 fprintf(fptr,"{\\f3\\fs%d\\cf1    }\n",gapsize);
		   }
		  } else {
			 fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		  }
		} else {
		   j++; site0++;
	// OUTPUT *crs file HERE...
	// (Later make a heap structure to output only the strongest sites...)
		// new format for VSI...
		   if(output_vsi && status[Analysis][j] == '!'){
		     Int4 bar_height = rtf[Analysis]->BarHeight(j);
		     if(bar_height >= 4){
		       if(pos_key_seq){
			 pcr->append(AlphaChar(qseq[j],AB),pos_key_seq[i],crs_color,bar_height);
		       } else {
			 pcr->append(AlphaChar(qseq[j],AB),site,crs_color,bar_height);
		       }
		     }
		   }
		   if(status[Analysis][j] == '!'){
			// fprintf(fptr,"{\\f3\\fs%d\\cf1 *}\n",gapsize);
			fprintf(fptr,"{\\f3\\fs%d\\cf1 \\u9679\\'3f}\n",gapsize);
		   } else if(status[Analysis][j] == '?'){
			fprintf(fptr,"{\\f3\\fs%d\\cf%d \\u9679\\'3f}\n",
				gapsize,get_font_color_code('L'));
		   } else {	// status == '*' or "." or "^"...
			fprintf(fptr,"{\\f3\\fs%d\\cf1  }\n",gapsize);
		   }
		   if(anal_index > NumAnalysis) sq_index[anal_index][1]++;
		   site++;
		}
		last_gnull=gnull[i];
	  } fprintf(fptr,"\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
	} // End of "if(status[Analysis])" section.
      }  // End of "else if(NumAnalysis==3)..." section.
	//**************** output secondary structure... *******************
	if(ss_struct && gnull_insrt_len){
	  fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d pred. structure",fontsize,color);
	  fprintf(fptr,"\\tab  \\tab }"); 
	  strcpy(str,"\\b");
	  for(i=gstart,j=start; i < gend && j < End; i++){
		if(gnull[i] == '_')
        		fprintf(fptr,"{%s\\f3\\fs%d\\cf1  }\n",str,gapsize);
		else {
			j++;
        		fprintf(fptr,"{%s\\f3\\fs%d\\cf1 %c}\n",str,gapsize,ss_struct[j]);
		}
	  } fprintf(fptr,"\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
	}
	//**************** output pvalue ... *******************
	if(verbose && gnull_insrt_len){
	 rtfQ[Analysis]->PutBarHeights(fptr,start,End,gstart,gend,gnull,
						color,colorB,rtfQ[analysis]);
	 rtfQ[Analysis]->PutLineBinomial(fptr,start,End,gstart,gend,gnull,
						color,colorB,rtfQ[analysis]);
	 // NEW: output Relative Entropy
	 if(!fract_seqs) rtf[Analysis]->PutLineRelEntropy(fptr,start,End,gstart,gend,gnull,
						color,colorB,rtf[analysis]);
	}
	//**************** SEQUENCE ALIGNMENT... *******************
	//**************** SEQUENCE ALIGNMENT... *******************
	//**************** SEQUENCE ALIGNMENT... *******************
	for(N=0, n=1; n<= nseqSMA(MA); n++){
		char id_str[30];
		Int4 res_printed=0; 
		Int4 colorK=color;
		char tmp_str[200],kingdom,king_color;

		if(n==1 && Seq4numbering == 2) continue;
		if(Analysis==0){
		   if(Kingdom){
                        switch(Kingdom[n]){
			   case 'b': case 'B': king_color='m'; break;
			   case 'a': case 'A': king_color='B'; break;
			   case 'f': case 'F': king_color='y'; break;
			   case 'd': case 'D': king_color='y'; break;
			   case 'e': case 'E': case 'p': case 'P': king_color='C'; break;
			   case 'm': case 'M': king_color='R'; break;
			   case 'v': case 'V': case 'g': case 'G': king_color='G'; break;
			   default: king_color='D'; break;
                        } colorK=get_font_color_code(king_color);
		   } 
		   if(StdAlignmentOnly){
			strncpy(id_str,seqidSMA(n,MA),28);
			// sprintf(id_str,"%s",seqidSMA(n,MA));
		   } else if(Phylum && Phylum[n]){
			strncpy(id_str,Phylum[n],28);
			//sprintf(id_str,"%s",Phylum[n]); 
		   } else sprintf(id_str,"unknown"); 
		} else {
			strncpy(id_str,seqidSMA(n,MA),28);
			// sprintf(id_str,"%s",seqidSMA(n,MA)); 
		}
		id_str[14]=0; // end at 14
		if(n < nseqSMA(MA)){
		  fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s }",fontsize,colorK,id_str);
		  fprintf(fptr,"{\\f2\\fs%d\\cf1 \\tab %d\\tab }",
				fontsize,sq_index[anal_index][n]);
		  strcpy(str,"\\b");
	        } else {	// underline last sequence...
		  fprintf(fptr,"{\\b\\f2\\fs%d\\ul\\cf%d %s }",fontsize,colorK,id_str);
		  fprintf(fptr,"{\\f2\\fs%d\\ul\\cf1 \\tab %d\\tab }",
						fontsize,sq_index[anal_index][n]);
	 	  strcpy(str,"\\b\\ul");
		}
        	char *gseq= gseqSMA(1,n,MA);
		unsigned char *seq = seqSMA(1,n,MA);
	        last_gnull=0;
		for(laststate=0,i=gstart,j=start; i < gend && j < End; i++){
		   if(gnull[i] == '_'){
		     state='i';
		     c = gseq[i];
		     if(isalpha(c)){ res_printed++; sq_index[anal_index][n]++; }
		     else { assert(isprint(c)); }

		     if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]){
		       if(last_gnull != '_'){
			 // Int4 numdots = 2 + (Int4) ceil(log10(gnull_insrt_len[i]));
			 Int4 numdots = 2 + (Int4) ceil(log10((double) gnull_insrt_len[i]+0.00001));
			 if(state == laststate){  // always false because gnull[i] != last_gnull
			    assert(!"this should not happen");
#if 0	// state == laststate) is always false because gnull[i] != last_gnull
			    if(isalpha(c)){
				// fprintf(fptr,"(x)");
				fprintf(fptr,"(%d)",gnull_insrt_len[i]);
			    } else {
				fprintf(fptr,"...");
				// fprintf(fptr,"---");
			    }
#endif
			 } else {
			   if(laststate != 0) fprintf(fptr,"}\n");
			   sprintf(tmp_str,"{%%s\\f3\\fs%%d\\cf15 (%%%dd)",numdots-2);
			   if(isalpha(c)){
			     Int4 num_Res=0;
			     for(Int4 col=i+gnull_insrt_len[i]-1; col >= i; col--){
				if(isalpha(gseq[col])) num_Res++;
			     }
			     fprintf(fptr,tmp_str,str,gapsize,num_Res);  
			   } else {
			     // fprintf(fptr,"{%s\\f3\\fs%d\\cf16 ...",str,gapsize);  
			     fprintf(fptr,"{%s\\f3\\fs%d\\cf15 ",str,gapsize);  
			     for(Int4 d = 1; d <= numdots; d++) fprintf(fptr,".");
			   }
			   // fprintf(fptr,"{%s\\f3\\fs%d\\cf8 (x)",str,gapsize);  
			 }
		       } // else skip over the rest...
		     } else {
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i'){
			    	  fprintf(fptr,"{%s\\f3\\fs%d\\cf15 %c",str,gapsize,c);  
			    } else {
				  fprintf(fptr,"{%s\\f3\\fs%d\\cf8 %c",str,gapsize,c);  
			    }
			  }
		     } laststate = state;

		   } else { 	// i.e., (gnull[i] != '_')
			j++;
			r = seq[j]; 
#if 0
			assert(!( r!=0 && gseq[i] == '-'));
#else
			if(r!=0 && gseq[i] == '-'){
			  Int4 jx=MAXIMUM(Int4,1,j-5);
			  fprintf(stderr,"context: '");
			  for( ; jx <= j; jx++){
				fprintf(stderr,"%c",AlphaChar(seq[jx],AB));
			  }
			  fprintf(stderr,"'\n %s(%s)(%d): %c%d(%d)\n",
				seqidSMA(n,MA),id_str,n,AlphaChar(r,AB),i,j);
			  fprintf(stderr,"%s\n",gseq);
			  fprintf(stderr,"j=%d; start=%d; i=%d; gstart=%d\n",j,start,i,gstart);
			  fprintf(stderr,"gnull=%s\n",gnull);
			  fprintf(stderr,"glength=%d\n",glengthSMA(1,MA));
			  // fprintf(fptr,"%s", MA->gseq[n][1]);
			  FILE *dbfp=open_file("debug_out",".msa","w");
			  PutSMA(dbfp,MA); fclose(dbfp);
			  dbfp=open_file("debug_out",".rtf","w");
			  SMA2RTF(dbfp,-2.0,1.8,2.5,MA); fclose(dbfp);
			  // MA created by:
			  // PutAlnCMSA(fp,CMA[1]); rewind(fp); MA=ReadSMA(fp); fclose(fp);
			  // need to make sure that the array is not over flow; appears to be bug.
#if 0
    for(j=0; j < glengthSMA(1,MA); j++){
	fprintf(stderr,"gnull[%d] = %c %d\n",j,gnull[j],gnull_insrt_len[j]);
    }
#endif
			  assert(!( r!=0 && gseq[i] == '-'));
			}
#endif
			if(r==0 && gseq[i] == '-'){ c = '-'; }
			else { c=AlphaChar(r,AB); res_printed++; sq_index[anal_index][n]++; }
			// New multiple Residues...
			if(NumResidues==1 && Residue[1] && j == Position0[1]){
			   if(MemSset(AlphaCode(c,AB),Residues[NumResidues])) state = '!';
			   else state = ' '; 
#if 1			// Added feature to allow first status == '^' position to be ignored
			} else if(rtf[analysis]->IgnoreResPos(j)){
				state = '!';
#endif
			} else {
			   state = rtf[analysis]->ConservedState(conserved[j], r,c); 
			}
			rtf[analysis]->PutConservedState(fptr, state, laststate, c, str);
			laststate = state;
		   }
		   last_gnull=gnull[i];
		} 
		// if(n==1) next_start=j; // turning this off was causing a core dump!
		// but I haven't checked this out rigorously yet...
		if(n < nseqSMA(MA)){
		  fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
		  if(res_printed == 0) fprintf(fptr,"\\tab  }{\\f3 \n\\par }");
		  else {
		    if(Seq4numbering > 0 && Seq4numbering == n){ 
		      fprintf(fptr,"\\tab %d*}{\\f3 \n\\par }",sq_index[anal_index][n]-1);
		    } else fprintf(fptr,"\\tab %d}{\\f3 \n\\par }",sq_index[anal_index][n]-1);
		  }
		} else {
		  fprintf(fptr,"}\n{\\f2\\fs%d\\ul ",fontsize);
		  if(res_printed == 0) fprintf(fptr,"\\tab  }{\\f3 \n\\par }");
		  else {
		    if(Seq4numbering > 0 && Seq4numbering == n){ 
			fprintf(fptr,"\\tab %d*}{\\f3 \n\\par }",sq_index[anal_index][n]-1);
		    } else fprintf(fptr,"\\tab %d}{\\f3 \n\\par }",sq_index[anal_index][n]-1);
		  }
		} lines_on_page++;
	}  /********************* END SEQUENCE ALIGNMENT *******************************/
	   /********************* END SEQUENCE ALIGNMENT *******************************/
	   /********************* END SEQUENCE ALIGNMENT *******************************/
	if(SuperAln[Analysis]){
	 if(show_marg_prob){
	   if(AlignName && AlignName[0]=='C'){
	     if(MargProb[Analysis]){
	  	put_marg_prob(fptr,gstart,start,gend,color,gapsize,gnull,Analysis);
	     }
	   }
	 }
   	 if(verbose){ 	// then put another histogram above super-alignment consensus.
	  if(AlignName && (AlignName[0]=='C' || AlignName[0]=='S' || AlignName[0]=='F')
		&& rtf[Analysis] && rtf[Analysis]->Printable() && Hist[Analysis])
	  { // Print another histogram above consensus...
	    fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s }{\\f3 \n\\par }\n",
						fontsize,color,AlignName);
	    lines_on_page++;
	    // fprintf(fptr,"{\\f3 \n\\par }\n"); // one carriage 
	    // fprintf(fptr,"{\\f2\\fs%d %s}\n",fontsize,AlignName);
	    // fprintf(fptr,"{\\f3 \n\\par }\n");
	    // rtf[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,tab3);
	    rtf[Analysis]->PutLine(fptr,start,End,gstart,gend,gnull,tab1,tab2,tab3,
							maxlen_gnull,gnull_insrt_len);
	  } 
	  rtf[Analysis]->PutBarHeights(fptr,start,End,gstart,gend,gnull,
						color,colorB,rtf[analysis]);
   	 }
	 //*********** output super-alignment info at bottom of figure **********
	 //*********** output super-alignment info at bottom of figure **********
	 //*********** output super-alignment info at bottom of figure **********

         //**************** output res_evals=ResEvals[Analysis] ****************
         //**************** output res_evals=ResEvals[Analysis] ****************
         //**************** output res_evals=ResEvals[Analysis] ****************
	 // next_start=rtf[Analysis]->PutResEvals(fptr,start,End,gstart,gend,gnull,
	 rtf[Analysis]->PutResEvals(fptr,start,End,gstart,gend,gnull,
			wtfreq,WtNumSeq[Analysis],color,rtf[Analysis],verbose,
			WtNumSq[Analysis],NumSqMain[Analysis],PatternOnlyMode,
			MaxConcensusLines,maxlen_gnull,gnull_insrt_len); 
	 lines_on_page+=MaxConcensusLines;
	}
	//**************** output fraction of sequences aligned... *******************
	if(fract_seqs){
#if 1	// NEW: output Relative Entropy
	 if(verbose){
	  rtf[Analysis]->PutLineRelEntropy(fptr,start,End,gstart,gend,gnull,
						color,colorB,rtf[analysis]);
	 }
#endif
#if 1	// NEW: PRINT OUT INSERTIONS AND DELETION INFORMATION HERE:
   if(!HideIndelInfo && gnull_insrt_len == 0){
	if(Analysis > 0 && MCMA[Analysis-1]){
	 for(char XX=1; XX >=0; XX--){
	  // if(!verbose && XX==0) break;  // always output deletions.
	  if(XX){
	    fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d insertions\\tab \\tab }",
		fontsize,color);
	    strcpy(str,"\\b");
	  } else {
	    fprintf(fptr,"{\\b\\f2\\fs%d\\ul\\cf%d deletions\\tab \\tab }",
		fontsize,color);
	    strcpy(str,"\\b\\ul");
	  }
	  char	indel;
	  cma_typ tmp_cma=MCMA[Analysis-1];
	  for(laststate=0,i=gstart,j=start; i < gend && j < End; i++){
		   if(gnull[i] == '_'){
			  state=' ';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		   } else {
			j++;
			Int4 hits=0;
			for(Int4 n_0=1; n_0 <= NumSeqsCMSA(tmp_cma); n_0++){
				if(XX){ if(InsertionCMSA(1,n_0,j,tmp_cma) > 0) hits++; }
				else { if(IsDeletedCMSA(1,n_0,j,tmp_cma)) hits++; }
			}
			state='X';
			double fract_hits=(double)hits/(double)NumSeqsCMSA(tmp_cma);
#if 0
			if(fract_hits < 0.10){
				indel=FractionToCharPBC(10.0*fract_hits);
				if(indel < '5') state = 'L';
				else state = 'D';
				if(indel=='0') indel=' ';
				else if(indel=='?') indel='!';
			} else {
				indel=FractionToCharPBC(fract_hits);
				if(indel < '6') state = 'R';
				else state = 'Y'; 
				if(indel=='?') indel='!';
			}
#else
			if(fract_hits < 0.10){
				indel=FractionToCharPBC(10.0*fract_hits);
				state = 'L';
				if(indel=='0') indel=' ';
				else if(indel=='?') indel='!';
			} else {
				indel=FractionToCharPBC(fract_hits);
				state = 'B'; 
				if(indel=='?') indel='!';
			}
#endif
			if(state == laststate) fprintf(fptr,"%c",indel);
			else { 
			   if(laststate != 0) fprintf(fptr,"}\n");
	        	   if(state == 'R') fprintf(fptr,"{%s\\f3\\cf6 %c",str,indel);
	        	   else if(state == 'Y') fprintf(fptr,
					"{%s\\f3\\cf6\\highlight7 %c",str,indel);
	        	   else if(state == 'B') fprintf(fptr,"{%s\\f3\\cf1 %c",str,indel);
	        	   else if(state == 'D') fprintf(fptr,"{%s\\f3\\cf15 %c",str,indel);
	        	   else fprintf(fptr,"{%s\\f3\\cf15 %c",str,indel);
			} laststate = state;
		   }
	  } fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
	  if(XX) fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  else fprintf(fptr,"}\n");	// for outputing res diversity.
	  lines_on_page++;
	 // else fprintf(fptr,"}\n{\\f2\\fs%d \\tab }",fontsize);	// for outputing res diversity.
	 }
	} else if(BG_CMA && Analysis == NumAnalysis+1){
	  //*********************** Special case for BG motif model... ********************
	  //*********************** Special case for BG motif model... ********************
	  //*********************** Special case for BG motif model... ********************
	  fprintf(fptr,"{\\b\\f2\\fs%d\\ul\\cf%d deletions\\tab \\tab }",
		fontsize,color);
	  strcpy(str,"\\b\\ul");
	  char	indel;
	  cma_typ tmp_cma=MCMA[Analysis-1];
	  for(laststate=0,i=gstart,j=start; i < gend && j < End; i++){
		   if(gnull[i] == '_'){
			  state=' ';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		   } else {
			j++;
			Int4 hits=0;
			state='X';
			if(NullFreqBG[j]==NullFreqBG[0]){  // then deleted...
				indel=FractionToCharPBC(0.99);
				state = 'B'; 
				if(indel=='?') indel='!';
			} else {
				indel=FractionToCharPBC(0.0);
				state = 'L';
				if(indel=='0') indel=' ';
				else if(indel=='?') indel='!';
			}
			if(state == laststate) fprintf(fptr,"%c",indel);
			else { 
			   if(laststate != 0) fprintf(fptr,"}\n");
	        	   if(state == 'B') fprintf(fptr,"{%s\\f3\\cf1 %c",str,indel);
			   // else fprintf(fptr,"{%s\\f3\\cf16 %c",str,indel);
			   else fprintf(fptr,"{%s\\f3\\cf15 %c",str,indel);
			} laststate = state;
		   }
	  } fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
	  fprintf(fptr,"}\n");	// for outputing res diversity.
	}
    }	// end HideIndelInfo if statement.
#endif 	// END: PRINT OUT INSERTIONS AND DELETION INFORMATION HERE:
#if 0
	  fprintf(fptr,"{\\b\\f2\\fs%d\\ul\\cf%d %% seqs. aligned \\tab %d\\tab }",
		fontsize,color,NumSqMain[Analysis]);
	  strcpy(str,"\\b\\ul");
	  char	fract_sqaln;
	  for(laststate=0,i=gstart,j=start; i < gend && j < End; i++){
		   if(gnull[i] == '_'){
			  state=' ';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		   } else {
			j++;
			state='X';
			fract_sqaln=FractionToCharPBC(fract_seqs[j]);
			if(fract_sqaln =='!') fract_sqaln='9';
			if(state == laststate) fprintf(fptr,"%c",fract_sqaln);
			else { 
			   if(laststate != 0) fprintf(fptr,"}\n");
	        	   fprintf(fptr,"{%s\\f3\\cf1 %c",str,fract_sqaln);
			} laststate = state;
		   }
	  } fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
	  if(FALSE && WtNumSq[Analysis] > 0) // print with wt_res_freq instead...
		fprintf(fptr,"}\n{\\f2\\fs%d\\ul\\cf%d \\tab %d}{\\f3 \n\\par }",
			fontsize,color,(Int4) WtNumSq[Analysis]);
	  else if(Analysis > 0 && MCMA[Analysis-1])  // skipped 1st analysis...
		fprintf(fptr,"}\n{\\f3\\fs%d \\tab }{\\f3 %.1f\n\\par }",fontsize,
		   ResidueDiversityCMSA(0,0,MCMA[Analysis-1]));
	  else fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
	  // NOTE: PERHAPS SHOULD PRINT OUT CMA DIVERSITY AT END OF THIS LINE!
#else 
	if(!HideIndelInfo){
	  fprintf(fptr,"\n{\\f2\\fs%d ",fontsize);
	  if(Analysis > 0 && MCMA[Analysis-1])  // skipped 1st analysis...
		fprintf(fptr,"}\n{\\f3\\fs%d \\tab }{\\f2 %.1f\n\\par }",fontsize,
		   ResidueDiversityCMSA(0,0,MCMA[Analysis-1]));
	  else fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	  lines_on_page++;
	}
	  // NOTE: PERHAPS SHOULD PRINT OUT CMA DIVERSITY AT END OF THIS LINE!
#endif
#if 1	// put BackGround below alignment...
	 // if(anal_index > 1)
	 if(anal_index > 2)	// no background shown for root node...
	 {
	   Int4 tmpAnal=AnalysisList[anal_index-1]; // look at previous
	   assert((MCMA[tmpAnal -1]) != 0 && tmpAnal > 0);
	   AlignName=NameCMSA(MCMA[tmpAnal-1]); // These are off by one...
	   assert(AlignName);
	   if(strcmp(AlignName,"BackGround")==0){
             double	**tmp_wtfreq=Freq[tmpAnal];	// background freq (= WtFrq[] )
	     color = 15; // dk grey
	     // color = 16; // lt grey
	     // ignore returned next_start here as this was set above....
	     // fprintf(fptr,"{\\f3 \n\\par }");
	     rtf[Analysis]->PutResEvals(fptr,start,End,gstart,gend,gnull,
			tmp_wtfreq,WtNumSeq[tmpAnal],color,rtf[tmpAnal],verbose,
			WtNumSq[tmpAnal],NumSqMain[tmpAnal],PatternOnlyMode,
			MaxConcensusLines,maxlen_gnull,gnull_insrt_len,TRUE,FALSE); 
	     lines_on_page+=MaxConcensusLines;
	   }
	 }
#endif
	}
	//**************** output marginal probability... *******************
	if(MargProb[Analysis] && SuperAln[Analysis] 
			&& !(AlignName && AlignName[0]=='C')){
	   put_marg_prob(fptr,gstart,start,gend,color,gapsize,gnull,Analysis);
	} // end of output marginal probability.
     } // end of Analysis loop
     PutCHA_RTNS(fptr,4);
   } // end of 'fancy-formating-of-page' loop...
   free(color_font);
   // PutCHA_RTNS(fptr,2);
   if(pos_key_seq) free(pos_key_seq);
   for(j=1; j<=lengthSMA(1,MA); j++) if(conserved[j]) free(conserved[j]);
   free(conserved);
   for(anal_index=0; anal_index <= NumAnalIndex; anal_index++) free(sq_index[anal_index]);
   // fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3\n\\par\n\\par }\n");
}

