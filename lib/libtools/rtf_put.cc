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

#include "rtf_typ.h"

char    rtf_typ::IntegerToChar(Int4 x)
{
    const char rtn[]=" 123456789!!!";
    if(x >=0 && x <= 10) return rtn[x];
    else return '?';
}

char    rtf_typ::FractionToChar(double fract)
{
    const char rtn[]="0123456789!!!";
    // Int4 f = (Int4) ceil((10.0*fract - 0.5));
    Int4 f = (Int4) floor((10.0*fract));
    if(f >= 0 && f <= 10) return rtn[f]; else return '?';
}

Int4	rtf_typ::PutResEvals(FILE *fptr,Int4 start,Int4 end,Int4 gstart,
	   Int4 gend,char *gnull,double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,
	   BooLean verbose,double WtNumSq,Int4 RawNumSeqs, BooLean PatternOnlyMode,
	   Int4 MaxConcensusLines)
{ return PutResEvals(fptr,start,end,gstart,gend,gnull,wtfreq,wtnsq,color,rtf,
	verbose,WtNumSq,RawNumSeqs, PatternOnlyMode,MaxConcensusLines,INT4_MAX,0); }

#define PUT_HIST_EVAL_WTFRQ 0

Int4	rtf_typ::PutResEvals(FILE *fptr,Int4 start,Int4 end,Int4 gstart,
	   Int4 gend,char *gnull,double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,
	   BooLean verbose,double WtNumSq,Int4 RawNumSeqs, BooLean PatternOnlyMode,
	   Int4 MaxConcensusLines,Int4 maxlen_gnull,Int4 *gnull_insrt_len)
{ return PutResEvals(fptr,start,end,gstart,gend,gnull,wtfreq,wtnsq,color,rtf,verbose,
		WtNumSq,RawNumSeqs,PatternOnlyMode,MaxConcensusLines,maxlen_gnull,gnull_insrt_len, FALSE,TRUE); }

Int4	rtf_typ::PutResEvals(FILE *fptr,Int4 start,Int4 end,Int4 gstart,
	   Int4 gend,char *gnull,double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,
	   BooLean verbose,double WtNumSq,Int4 RawNumSeqs, BooLean PatternOnlyMode,
	   Int4 MaxConcensusLines,Int4 maxlen_gnull,Int4 *gnull_insrt_len, BooLean UnderLine, BooLean IsFG)
// find consensus residues and corresponding pvalues.
{
	char	*consensus[27];
	double	*cons_pval[27];
	Int4	index,high_index,max_index=25;
	double	*tmp_pval;
	Int4	factor=2,x,i,j;
	Int4	gapsize=fontsize;
	char	c,r,state,laststate,str[50];
	//============= First sort by Pvalues ================
	double sum_pvals=0.0;
	Int4	num_pvals=0,next_start=0;
	dh_type dH=dheap(max_index+2,4);
	char	**conserved=0;

	NEWP(conserved,length+3,char);
        for(i=gstart,j=start; i < gend && j < end; i++){
          if(gnull[i] != '_'){
                // j++; conserved[j]=rtf->Conserved(res_evals[j]);
                j++; conserved[j]=rtf->Conserved(j);
          }
        }
	for(x=1; x <= max_index; x++){
	    	MEW(consensus[x],length+3,char);
	    	MEW(cons_pval[x],length+3,double);
	} high_index=0;
#if PUT_HIST_EVAL_WTFRQ
	h_type	HG=Histogram("res_evals",-100,100,1.0);
	h_type	HG2=Histogram("wtfreq",0,5,0.05);
#endif
	for(i=gstart,j=start; i < gend && j < end; i++){
	    if(gnull[i] != '_'){
		j++;
	        tmp_pval=res_evals[j];
	        for(index=0,r=1; r<=nAlpha(AB); r++){
		   if(tmp_pval[r] <= -2.0){ sum_pvals+=tmp_pval[r]; num_pvals++; }
		   if(!wtfreq || (wtfreq[j] && wtfreq[j][r] > 0.1)) 
		      insrtHeap(r,((keytyp)(tmp_pval[r])),dH);
#if PUT_HIST_EVAL_WTFRQ
		   IncdHist(tmp_pval[r],HG);
		   if(wtfreq && (wtfreq[j] && wtfreq[j][r])) IncdHist(wtfreq[j][r],HG2); 
#endif
		   cons_pval[r][j] = 0.0;
		   consensus[r][j] = ' ';
	        }
	        while((r=delminHeap(dH)) != 0){
#if 1
			if(tmp_pval[r] > ExpPttrns) break;
			// ExpPttrns is set by chn_typ from calling environment 
			// corresponds to ExpPatterns of chn_typ
			// if(tmp_pval[r] > 20.0) break;
			// setting to 20.0 = expect to show ~20 chance residue positions...
			if((wtfreq[j] &&  wtfreq[j][r] < 0.1)) break;
#elif 0 
			if(tmp_pval[r] > -1.3) break;
			// if(tmp_pval[r] > -2.0) break;
			// if((wtfreq[j] &&  wtfreq[j][r] < 0.1)) break;
#else 
			if(tmp_pval[r] > -2.0 && 
			    (!wtfreq || (wtfreq[j] &&  wtfreq[j][r] <= 0.1))) break;
#endif
			index++;
			cons_pval[index][j] = tmp_pval[r];
			consensus[index][j] = AlphaChar(r,AB);
			if(index >= max_index) break;
		}
		while(delminHeap(dH));	// empty heap
		if(index > high_index) high_index = index;
	    }
	} Nildheap(dH);
#if PUT_HIST_EVAL_WTFRQ
	PutHist(stderr,60,HG); NilHist(HG);
	PutHist(stderr,60,HG2); NilHist(HG2);
#endif
	double	target_pval=5.0;
	double	ave_pval=target_pval;
	if(num_pvals > 0){ ave_pval = -sum_pvals/(double)num_pvals; }
	factor = (Int4) floor((ave_pval/target_pval) + 0.5); 
	factor = MAXIMUM(Int4,1,factor);
#if 1	// change factor
	if(MaxPval <= 0.0){ factor = 1; }
	else {
	  factor = (Int4) ceil((MaxPval/9.0));
	  factor = MAXIMUM(Int4,1,factor);
	}
#endif
	//============= output conserved residues ================
	char last_gnull=0;
	if(high_index > MaxConcensusLines){
		high_index=MaxConcensusLines;
	} 
	if(high_index > 10){
		fprintf(stderr,"********* high_index=%d; ExpPttrns=%g *********\n",
								high_index,ExpPttrns);
	}
	for(index=1; index <= high_index; index++){
		char id_str[30]; 
		if(index==1) {
			if(PatternOnlyMode) sprintf(id_str,"%s (%d):","pattern",RawNumSeqs);
			else if(IsFG) sprintf(id_str,"%s (%d):","foreground",RawNumSeqs);
			else sprintf(id_str,"%s (%d):","background",RawNumSeqs);
			// sprintf(id_str,"%s(%.2f):","conserved",rtf->GetLinearToLog());
		} else sprintf(id_str," ");
		// fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s \\tab %d }",fontsize,color,id_str,RawNumSeqs);
		fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s }",fontsize,color,id_str);
#if 0
		if(index==1) fprintf(fptr,"{\\f2\\fs%d\\cf1 \\tab %d\\tab }", fontsize,start+1);
		else 
#endif
		fprintf(fptr,"{\\f2\\fs%d\\cf1 \\tab \\tab }", fontsize);
		// fprintf(fptr,"{\\f2\\fs%d\\cf1 \\tab }", fontsize);
		strcpy(str,"\\b");
		last_gnull=0;
		for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
		   if(gnull[i] == '_'){
		     state='i';
		     c = ' ';
		     if(gnull_insrt_len && gnull_insrt_len[i] >= maxlen_gnull){
			if(last_gnull != '_'){		// first insert position
			  if(state != laststate){	
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i'){
			    	  // fprintf(fptr,"{%s\\f3\\fs%d\\cf16 ",str,gapsize);  
			    	  fprintf(fptr,"{%s\\f3\\fs%d\\cf15 ",str,gapsize);  
			    } else {
				  fprintf(fptr,"{%s\\f3\\fs%d\\cf8 ",str,gapsize);  
			    }
			  }
			  // Int4 numspace = 2 + (Int4) ceil(log10(gnull_insrt_len[i]));
			  Int4 numspace = 2 + (Int4) ceil(log10((double) gnull_insrt_len[i]+0.00001));
			  for(Int4 d = 1; d <= numspace; d++) fprintf(fptr," ");
			  laststate = state;
			}
		     } else {
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i'){
			    	  // fprintf(fptr,"{%s\\f3\\fs%d\\cf16 %c",str,gapsize,c);  
			    	  fprintf(fptr,"{%s\\f3\\fs%d\\cf15 %c",str,gapsize,c);  
			    } else {
				  fprintf(fptr,"{%s\\f3\\fs%d\\cf8 %c",str,gapsize,c);  
			    }
			  } laststate = state;
		     }
		   } else {
			j++;
			c = consensus[index][j]; r = AlphaCode(c,AB);
			// if(Observed && !Observed[j]) state='W'; else 
			if(c == ' '){
			  state='i';
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    // fprintf(fptr,"{%s\\f3\\fs%d\\cf16 %c",str,gapsize,c);  
			    fprintf(fptr,"{%s\\f3\\fs%d\\cf15 %c",str,gapsize,c);  
			  }
			} else {
			  if(UnderLine){
		            state='z';
			    if(state == laststate){ fprintf(fptr,"%c",c); }
                            else {
                              if(laststate != 0) fprintf(fptr,"}\n");
	     		      fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
			    }
			  } else {
			    state = rtf->ConservedState(conserved[j], r,c);
			    char tmp_str[3]; tmp_str[0]=state; tmp_str[1]=0;
			    if(PatternOnlyMode && strstr("xyu*",tmp_str)){ 
				  state='W'; // white on white...  
			    } rtf->PutConservedState(fptr,state,laststate,c,str);
			  }
			} laststate = state;
		   } last_gnull=gnull[i];
		} next_start=j;	// WARNING: TWO DISTINCT PATHS TO SET next_start...
		fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
		// if(index==1) fprintf(fptr,"\\tab %d",j);
		if(index==1){
#if 1
			double log_two_fold=log10(pow(2.0,1.0/(1.0-rtf->GetLinearToLog())));
			// double log_two_fold=pow(2.0,1.0/(1.0-rtf->GetLinearToLog()));
			fprintf(fptr,"\\tab %.1f",log_two_fold);
#else
			fprintf(fptr,"\\tab %.0f",pow(100,1.0-rtf->GetLinearToLog()));
#endif
		} else fprintf(fptr,"\\tab ");
		fprintf(fptr,"}{\\f3 \n\\par }");
	}
	//============= output binomial tail probabilities ================
   if(verbose && gnull_insrt_len == 0){
	for(index=1; index <= high_index; index++){
	    
	    if(index==1) fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d log-#trials (%dx)",
				fontsize,color,factor);
	    else fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d ",fontsize,color);
	    fprintf(fptr,"\\tab  \\tab }"); 
	    strcpy(str,"\\b");
	    for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
		   if(gnull[i] == '_'){
			  state='x';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		   } else {
			j++;
			char	pchar=rtf->PvalueToChar(1.0/(double)factor,
						cons_pval[index][j]);
			if(pchar == '0') pchar = ' ';
			// state = rtf->ShowInfo(fptr,pchar,laststate,str);
			state = ShowInfo(fptr,pchar,laststate,str);
			laststate = state;
		   }
	    } fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	}
   }
	//============= output weighted residue frequencies ================
	if(!PatternOnlyMode && wtfreq){
	   for(index=1; index <= high_index; index++){
	    if(index==1){
		fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d wt_res_freqs (%d):\\tab \\tab }",
				fontsize,color,(Int4)ceil(WtNumSq-0.5));
	    } else if(UnderLine && index==high_index){
		 fprintf(fptr,"{\\b\\ul\\f2\\fs%d\\cf%d \\tab  \\tab }",fontsize,color);
	    } else fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d \\tab  \\tab }",fontsize,color);
	    if(index == high_index && UnderLine){ strcpy(str,"\\b\\ul"); } else strcpy(str,"\\b");
	    last_gnull=0;
	    for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
		   if(gnull[i] == '_'){
		     state='x';
		     if(gnull_insrt_len && maxlen_gnull <= gnull_insrt_len[i]){
			if(last_gnull != '_'){
#if 0
			  if(state == laststate) fprintf(fptr,"   ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1    ",str,gapsize);
                          }
#else
			  if(state != laststate){
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1 ",str,gapsize);
                          }
			  // Int4 numspace = 2 + (Int4) ceil(log10(gnull_insrt_len[i]));
			  Int4 numspace = 2 + (Int4) ceil(log10((double)gnull_insrt_len[i] + 0.00001));
			  for(Int4 d = 1; d <= numspace; d++) fprintf(fptr," ");
#endif
			  laststate = state;
		        }
		     } else {
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		     }
		   } else {
			j++;
			char pchar=' ';;
			c = consensus[index][j];
			if(c != ' '){
	                  r = AlphaCode(c,AB);
#if 0		// I turned this off: fractions don't add up...
			  // NEW...print out freq without PSEUDOCOUNTS: 
			  double wtfrq;
		   	  if(wtfreq[j]){
			    if(wtnsq){ 
				// WARNING: assumes 1 pseudocount each for 20 residues
				wtfrq = (wtfreq[j][r] * (wtnsq[j] + 20.0)) - 1.0;
				wtfrq = wtfrq /wtnsq[j];
				if(wtfrq < 0.0) wtfrq=0.0;
			    } else {
				// wtfrq = wtfreq[j][r];
				wtfrq = 0.0;
			    }
		 	    if(wtfrq > 0.1){
			       pchar=FractionToChar(wtfrq);
			       if(pchar=='?') pchar='!';
			    } else pchar=' ';
			  } else pchar=' ';
#else 	// OLD...
		   	  if(wtfreq[j] && wtfreq[j][r] > 0.1) pchar=FractionToChar(wtfreq[j][r]);
			  else pchar=' ';
#endif
			  if(pchar=='0') pchar = ' ';
			  else if(pchar=='!') pchar = '9';
			}
			if(UnderLine){		// UnderLine == Background patterns & freqs.
			  // for background frequencies don't highlight.
		          state='z';
			  if(state == laststate){ fprintf(fptr,"%c",pchar); }
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
	     		    fprintf(fptr,"{%s\\f3\\cf15 %c",str,pchar);
			  }
			} else {
				// state = rtf->ShowInfo(fptr,pchar,laststate,str);
				state = ShowInfo(fptr,pchar,laststate,str);
			} laststate = state;
		   } last_gnull=gnull[i];
	    } fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
	   }
	}
	for(x=1; x <= max_index; x++){ free(cons_pval[x]); free(consensus[x]); }
	for(j=1; j<=length; j++) if(conserved[j]) free(conserved[j]);
	free(conserved);
	return next_start;
}

void    rtf_typ::PutLineBinomial(FILE *fptr,Int4 start,Int4 end, Int4 gstart,
		Int4 gend, char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2)
{
	char	str[100],state,laststate;
	Int4	i,j,gapsize=fontsize;
	Int4	factor;     // factor to bring PValue into range 0..10
	
	assert(start >= 0 && start < length);
	if(MaxPval <= 0.0){ factor = 1; }
	else {
	  factor = (Int4) ceil((MaxPval/9.0));
	  factor = MAXIMUM(Int4,1,factor);
	}
	fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d log-#trials (%dx)\\tab}",
		fontsize,color,factor);
	// fprintf(fptr,"{\\b\\f3\\cf%d\\highlight%d *}\n",colorB,colorB); 

	if(color==colorB) fprintf(fptr,"{\\b\\f3\\cf8\\highlight%d x}\n",colorB); 
	else fprintf(fptr,"{\\b\\f3\\cf%d\\highlight%d x}\n",colorB,colorB); 

	fprintf(fptr,"{\\b\\f2\\fs%d \\tab}", fontsize);
        strcpy(str,"\\b");
        for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
            if(gnull[i] == '_'){
                 state='x';
                 if(state == laststate) fprintf(fptr," ");
                 else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                 } laststate = state;
            } else {
                j++;
		assert(j <= length);
                // char pval_char=PvalueToChar(0.5,-PValue[j]);
                char pval_char;
	        if(j == Position) pval_char='*';
		else pval_char=PvalueToChar(1.0/(double)factor,-PValue[j]);
#if 0
		if(rtf2){	// use first rtf.
		  state = rtf2->ShowInfo(fptr,pval_char,laststate,str);
		} else {
                  state = ShowInfo(fptr,pval_char,laststate,str);
		}
#else 
                state = ShowInfo(fptr,pval_char,laststate,str);
#endif
                laststate = state;
            }
        } 
#if 0	// Not sure want this! I DON'T SEEM TO BE USING rtf2!!!!!!!!!!!!!!!
	if(rtf2) fprintf(fptr,"}\n{\\f2\\fs%d \\tab %.2f}{\\f3 \n\\par }",
			fontsize,rtf2->LinearToLog());
	else
#endif
#if 1
	double log_two_fold=log10(pow(2.0,1/(1.0-LinearToLog)));
	// double log_two_fold=pow(2.0,1.0/(1.0-LinearToLog));
	fprintf(fptr,"}\n{\\f2\\fs%d \\tab %.1f}{\\f3 \n\\par }",fontsize,log_two_fold);
#else
	fprintf(fptr,"}\n{\\f2\\fs%d \\tab %.0f}{\\f3 \n\\par }",fontsize,
				pow(100,1.0-LinearToLog));
#endif
}

void    rtf_typ::PutLine(FILE *fptr,Int4 start,Int4 end, Int4 gstart,Int4 gend,char *gnull,
		Int4 tab1,Int4 tab2,Int4 tab3)
{ PutLine(fptr,start,end, gstart,gend,gnull,tab1,tab2,tab3,INT4_MAX,0); }

void    rtf_typ::PutLine(FILE *fptr,Int4 start,Int4 end, Int4 gstart,Int4 gend,char *gnull,
		Int4 tab1,Int4 tab2,Int4 tab3,Int4 maxlen_gnull,Int4 *gnull_insrt_len)
// Output an histogram of selective constraint Values for conserved positions.
{
        char    laststate,state,last_gnull=0;
        Int4    i,j;

        fprintf(fptr,
            "\\pard \\sl14\\slmult1\\widctlpar\\tqr\\tx%d\\tx%d\\tx%d\\adjustright \\fs%d ",
                tab1,tab2,tab3,fontsize);  // NOTE: sl14 == crunches down lines...
#if 1  // add a line at the top
        fprintf(fptr,"{\\f3\\cf1\\up%d\\cgrid \\tab \\tab }\n",hist_up);
        fprintf(fptr,"{\\f3\\up%d\\cgrid \\par }\n",hist_up);
#endif
        // for(Int4 line=hist_height; line >= 1; line--)
        for(Int4 line=hist_height; line >= 0; line--)
	{
           fprintf(fptr,"{\\f3\\cf1\\up%d\\cgrid \\tab \\tab }\n",hist_up);
           for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
                if(gnull[i] != '_'){
                  j++;
                  // Int4 I = GetHistHeight(IgnorePos[j],Value[j]); // AFN 12/22/05
                  Int4 I = GetHistHeight(Value[j]);

		  if(line == 0){		// line on bottom for significance.
		    if(I > 0) state = 'R';
		    else {
			// if(PValue[j] >= 3.0) state = 'D';
		    	// else if(PValue[j] >= 2.0) state = 'G';
		    	// if(PValue[j] >= 4.0) state = 'G';
		    	if(PValue[j] >= MinPval) state = 'G';
		    	else state = 'W';
		    }
		  } else if(I >= line){ 
			state = 'R';	// use only red for now...
                  } else state = 'W';

                  if(state == laststate) {
                        switch(state){
			  case 'D': case 'G':
			  case 'M': case 'R': case 'B': fprintf(fptr,"_"); break;
                          case 'W': fprintf(fptr," "); break;
                          default: print_error("This should not happen...");
                                break;
                        }
                  } else {	// state != laststate
                        if(laststate) fprintf(fptr,"}\n");
                        switch(state){
                          case 'G': // lt gray.
                            // fprintf(fptr,"{\\b\\f3\\ulth\\cf16\\up%d\\cgrid _",hist_up); break;
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf15\\up%d\\cgrid _",hist_up); break;
			    // get ride of bold on first line...
                            // fprintf(fptr,"{\\f3\\ulth\\cf15\\up%d\\cgrid _",hist_up); break;
                          case 'D': // dark gray.
                            // fprintf(fptr,"{\\f3\\ulth\\cf15\\up%d\\cgrid _",hist_up); break;
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf15\\up%d\\cgrid _",hist_up); break;
                          case 'M': // magenta 
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf5\\up%d\\cgrid _",hist_up); break;
                          case 'R': // red
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf6\\up%d\\cgrid _",hist_up); break;
                          case 'B': // brown
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf13\\up%d\\cgrid _",hist_up); break;
                          case 'W': // white
                            fprintf(fptr,"{\\b\\f3\\ulth\\cf8\\up%d\\cgrid  ",hist_up); break;
                          default: print_error("This should not happen..."); break;
                        }
                  }
                } else {		// gnull[i] == '_'
                  state='W';
		  if(gnull_insrt_len && gnull_insrt_len[i] >= maxlen_gnull){
							// use hide inserts option
		    if(gnull[i] != last_gnull){	// i.e., for the first '_' seen...
                        if(state != laststate){
                          if(laststate) fprintf(fptr,"}\n");
                          fprintf(fptr,"{\\b\\f3\\ulth\\cf8\\up%d\\cgrid ",hist_up);
			}
			// Int4 numspace = 2 + (Int4) ceil(log10(gnull_insrt_len[i]));
			Int4 numspace = 2 + (Int4) ceil(log10((double) gnull_insrt_len[i] + 0.00001));
			for(Int4 d = 1; d <= numspace; d++) fprintf(fptr," ");
		    }
		  } else {
                    if(state == laststate) fprintf(fptr," ");
                    else {
                        if(laststate) fprintf(fptr,"}\n");
                        fprintf(fptr,"{\\b\\f3\\ulth\\cf8\\up%d\\cgrid  ",hist_up);
                    }
		  }
                } laststate=state; last_gnull=gnull[i];
           } fprintf(fptr,"}\n{\\f3\\up%d\\cgrid \\tab }",hist_up);
           fprintf(fptr,"{\\f3\\up%d\\cgrid \\par }\n",hist_up);
        } fprintf(fptr,"{\\f3\\up%d\\cgrid \\par }\n",hist_up);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
        fprintf(fptr,"\\widctlpar\\tqr\\tx%d\\tx%d\\tx%d\\adjustright \\fs%d\\cgrid \n",
                tab1,tab2,tab3,fontsize);
}

void    rtf_typ::PutBarHeights(FILE *fptr,Int4 start,Int4 end, Int4 gstart,
		Int4 gend, char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2)
{
	char	str[100],state,laststate;
	Int4	i,j,gapsize=fontsize;
	
	assert(start >= 0 && start < length);
	fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d bar height \\tab}", fontsize,color);

	// fprintf(fptr,"{\\b\\f3\\cf%d\\highlight%d *}\n",colorB,colorB); 
	if(color==colorB) fprintf(fptr,"{\\b\\f3\\cf8\\highlight%d x}\n",colorB); 
	else fprintf(fptr,"{\\b\\f3\\cf%d\\highlight%d x}\n",colorB,colorB); 

	fprintf(fptr,"{\\b\\f2\\fs%d \\tab}", fontsize);
        strcpy(str,"\\b");
        for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
            if(gnull[i] == '_'){
                 state='x';
                 if(state == laststate) fprintf(fptr," ");
                 else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                 } laststate = state;
            } else {
                j++;
		assert(j <= length);
		Int4 I = GetHistHeight(Value[j]);
		// Int4 I = GetHistHeight(IgnorePos[j],Value[j]); // AFN 12/22/05
                char pval_char;
	        if(j == Position) pval_char='*';
		else if(I<=0) pval_char=' ';
		else if(I > hist_height) pval_char='!';
		else if(I <= 9) pval_char=IntegerToChar(I);
		else if(I <= 35) pval_char='a'+(char)(I-10); // assumes height ranges from 10..35
		else if(I <= 61) pval_char='A'+(char)(I-36); // assumes height ranges from 36..61
	 	else print_error("Histogram Height is set too high");
		// WARNING: temperary fix for '?' problem...
	        if(pval_char == '?') pval_char = '!';
		// find out what is wrong (probably pseudocounts)
#if 0
		if(rtf2){	// use first rtf.
		  state = rtf2->ShowInfo(fptr,pval_char,laststate,str);
		} else {
                  state = ShowInfo(fptr,pval_char,laststate,str);
		}
#else 
                  state = ShowInfo(fptr,pval_char,laststate,str);
#endif
		laststate = state;
            }
        } 
#if 1
	fprintf(fptr,"}\n{\\f2\\fs%d \\tab %.2f}{\\f3 \n\\par }",fontsize,LinearToLog);
#else
	fprintf(fptr,"}\n{\\f2\\fs%d \\tab %.0f}{\\f3 \n\\par }",
				fontsize,pow(100,1.0-LinearToLog));
#endif
}

//*************************** from RMA_TYP ****************************

char	rtf_typ::ConservedState(char *conserved, char r,char c)
{
	char state;

	if(c=='-') state='u';
	else if(c=='?') state='u';
	else if(conserved==NULL || !conserved[r]) state='u';
	else {
	  switch (conserved[r]){
	    case 'w': 	// weak...
		switch(c){
		  case 'G': case 'P': state = 'x'; 
		    break;
		  case 'C': 
		  case 'A': case 'I': case 'L': case 'V': case 'M': 
		  case 'F': case 'W': case 'Y': 
			state = 'y'; break;
		  case 'H': case 'K': case 'R': 
		  case 'N': case 'Q': case 'S': case 'T': 
		  case 'D': case 'E': state = '*';
		    break;
		  default: print_error("get_conserved_state(w): This should not happen!");
		    break;
		}
	     break;
	    case 'm': case 's':
		switch(c){
		  case 'G': state = 'g'; break;
		  case 'P': state = 'p'; break;
		  case 'C': state = 'c'; break;
		  case 'A': case 'I': case 'L': case 'V': case 'M': 
			state = 'n'; break;
		  case 'F': case 'W': case 'Y': state = 'a'; break;
		  case 'H': state = 'h'; break;
		  case 'K': case 'R': state = 'b'; break;
		  case 'N': case 'Q': case 'S': case 'T': 
			state = 'o'; break;
		  case 'D': case 'E': state = 'd'; break;
		  default: print_error("get_conserved_state(ms): This should not happen!");
		   break;
		}
		if(conserved[r]=='s') state = toupper(state);
	     break;
	    default: print_error("get_conserved_state(?): This should not happen!");
	     break;
	  }
	} return state;
}

void	rtf_typ::PutConservedState(FILE *fptr,char state,char laststate,char c,char *str)
#if 0	//************************* states *********************************
   	group:		weak:		moderate:	Strong:
Turns:
	G:		'x'		'g'		'G'
	P:		'x'		'p'		'P'
Nonpolar:
	C:		'y'		'c'		'C'
	AILVM:		'y'		'n'		'N'
	FWY:		'y'		'a'		'A'
polar:
	H:		'*'		'h'		'H'
	KR:		'*'		'b'		'B'
	NQST:		'*'		'o'		'O'
	DE:		'*'		'd'		'D'
Defines set:
	selected residue: '!'
	other residue:    ' '
unconserved:              'u';
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
#endif	//******************************************************************
{
	if(state == laststate) fprintf(fptr,"%c",c);
	else {
	   if(laststate != 0) fprintf(fptr,"}\n");
	   switch(state){
	     case '!': fprintf(fptr,"{%s\\f3\\cf1\\highlight6 %c",str, c); break;  // black on red
	     case ' ': fprintf(fptr,"{%s\\f3\\cf6\\highlight1 %c",str, c); break;  // red on black
	     case 'G': fprintf(fptr,"{%s\\f3\\cf8\\highlight4 %c",str, c); break;  // white on green
	     case 'P': fprintf(fptr,"{%s\\f3\\cf8\\highlight1 %c",str, c); break;  // white on black

	     case 'C': fprintf(fptr,"{%s\\f3\\cf2\\highlight7 %c",str, c); break;  // blue on yellow
	     case 'N': fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str, c); break;  // red on yellow
	     case 'A': fprintf(fptr,"{%s\\f3\\cf5\\highlight7 %c",str, c); break;  // magenta on yellow

	     case 'H': fprintf(fptr,"{%s\\f3\\cf8\\highlight2 %c",str, c); break;  // white on blue
	     case 'B': fprintf(fptr,"{%s\\f3\\cf8\\highlight3 %c",str, c); break;  // white on cyan
	     case 'O': fprintf(fptr,"{%s\\f3\\cf8\\highlight5 %c",str, c); break;  /** white on magenta **/
	     case 'D': fprintf(fptr,"{%s\\f3\\cf8\\highlight6 %c",str, c); break;  /** white on red **/

	     case 'g': fprintf(fptr,"{%s\\f3\\cf4 %c",str,c); break;  /** green on white **/
	     case 'p': fprintf(fptr,"{%s\\f3\\cf1 %c",str,c); break;  /** black on white **/

	     case 'c': case 'n': 
	     case 'a': fprintf(fptr,"{%s\\f3\\cf14\\highlight7 %c",str, c);
		 break;  /** dk yellow on yellow **/

	     case 'h': fprintf(fptr,"{%s\\f3\\cf2 %c",str,c); break;  /** blue on white **/
	     case 'b': fprintf(fptr,"{%s\\f3\\cf3 %c",str,c); break;  /** cyan on white **/
	     case 'o': fprintf(fptr,"{%s\\f3\\cf5 %c",str,c); break;  /** magenta on white **/
	     case 'd': fprintf(fptr,"{%s\\f3\\cf6 %c",str,c); break;  /** red on white **/

	     case 'y': case '*':
	     case 'x': fprintf(fptr,"{%s\\f3\\cf15 %c",str,c); break;  /** dk grey on white **/
	     case 'W': fprintf(fptr,"{%s\\f3\\cf8 %c",str,c); break;  /** white on white **/
	     case 'u':
	     default : 
		  // fprintf(fptr,"{%s\\f3\\cf16 %c",str,c);
		// break;  /** light grey on white **/
		  fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
		break;  /** dk grey on white **/
	   } 
	}
}

char	rtf_typ::PvalueToChar(double factor, double pval)
// convert a pvalue to a symbolic character; scaling factor is used.
{
	char	pval_show;
	Int4	tmp_d = (Int4) floor((-pval*factor)+0.5);
	if(tmp_d < 0) tmp_d=0;
	switch(tmp_d){
		case 0: pval_show = '0'; break;
		case 1: pval_show = '1'; break;
		case 2: pval_show = '2'; break;
		case 3: pval_show = '3'; break;
		case 4: pval_show = '4'; break;
		case 5: pval_show = '5'; break;
		case 6: pval_show = '6'; break;
		case 7: pval_show = '7'; break;
		case 8: pval_show = '8'; break;
		case 9: pval_show = '9'; break;
		default: pval_show = '!'; break;
	} return pval_show;
}

char	rtf_typ::ShowInfo(FILE *fptr, char info_show, char laststate, char *str)
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
{
	char	state;
	switch (info_show){
	  case '*': state = '*'; break;
	  case '0': info_show=' ';
	  case '1': case '2': state='g'; break;
	  case '3': case '4': state='G'; break;
	  case '5': case '6': state='B'; break;
	  case '7': case '8': state='r'; break;
	  case '9': case '!': state='R'; break;
	  case ' ': state='u'; break;
	  case '-': state='u'; break;
	  default: 
		char sstr[3]; sstr[0]=info_show; sstr[1]=0;
		if(strstr("ABC",sstr)) state = 'g';
		else if(strstr("DEF",sstr)) state='G';
		else if(strstr("GHIJK",sstr)) state='B';
		else if(strstr("LMNOPQR",sstr)) state='r';
		else if(strstr("STUVWXYZ",sstr)) state='R';
		else state='u'; break;
	}
	if(state == laststate) fprintf(fptr,"%c",info_show);
	else {
	   if(laststate != 0) fprintf(fptr,"}\n");
	   switch (state) {
	     case '*': // black on red
	        fprintf(fptr,"{%s\\f3\\cf1\\highlight6 *",str); break;
	     case 'w': // white
	        fprintf(fptr,"{%s\\f3\\cf8 %c",str,info_show); break;
	     case 'g': // light grey --> dark gray
	        fprintf(fptr,"{%s\\f3\\cf15 %c",str,info_show);
	        // fprintf(fptr,"{%s\\f3\\cf16 %c",str,info_show);
		break;
	     case 'G': // grey
	        fprintf(fptr,"{%s\\f3\\cf15 %c",str,info_show); break;
	     case 'b': // blue
	        fprintf(fptr,"{%s\\f3\\cf2 %c",str,info_show); break;
	     case 'B': // black
	        fprintf(fptr,"{%s\\f3\\cf1 %c",str,info_show); break;
	     case 'r': // red
	        fprintf(fptr,"{%s\\f3\\cf6 %c",str,info_show); break;
	     case 'R':  
		fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str,info_show); break;  // red on yellow
	     case ' ': 
	     default: fprintf(fptr,"{%s\\f3\\cf1 %c",str,info_show); break;
	   }
	} return state;
}

Int4	get_font_color_code(char ColorCode)
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
{
	Int4	color;

    switch (ColorCode){		// get color for sequence name...
	case 'R': color = 6; break;	// red
	case 'B': color = 2; break;	// blue
	case 'G': color = 4; break;	// green
	case 'M': color = 5; break;	// magenta
	case 'C': color = 3; break;	// cyan
	case 'y': color = 14; break;	// dk yellow
	case 'r': color = 13; break;	// dk red
	case 'b': color = 9; break;	// dk blue
	case 'g': color = 11; break;	// dk green
	case 'm': color = 12; break;	// violet
	case 'c': color = 10; break;	// teal
	// case 'L': color = 16; break;	// lt grey
	case 'L': color = 15; break;	// dk grey; get rid of lt grey
	case 'D': color = 15; break;	// dk grey
	case 'W': color = 8; break;	// white
	case 'Y': color = 7; break;	// yellow
	case 'P': color = 17; break;	// pink
	case 'O': color = 18; break;	// orange
	case 'o': color = 19; break;	// dk orange
	case 'n': color = 20; break;	// dk cyan 
	default: color = 1; break;       // black
    } return color;
}

void    rtf_typ::PutLineRelEntropy(FILE *fptr,Int4 start,Int4 end, Int4 gstart,
                Int4 gend, char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2)
{
        char    str[100],state,laststate;
        Int4    i,j,gapsize=fontsize;

        assert(start >= 0 && start < length);
        fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d Info (3rd bits)\\tab}", fontsize,color);

	if(color==colorB) fprintf(fptr,"{\\b\\f3\\cf8\\highlight%d x}\n",colorB); 
	else fprintf(fptr,"{\\b\\f3\\cf%d\\highlight%d x}\n",colorB,colorB); 

	fprintf(fptr,"{\\b\\f2\\fs%d \\tab}", fontsize);
#if 0
        // NEW: miniture Venn diagram showing background family color.
        fprintf(fptr,"{\\b\\f2\\cf%d\\highlight%d\n",color,colorB);
        fprintf(fptr,"{\\field{\\*\\fldinst SYMBOL 183 \\\\f \"Symbol\" \\\\s %d}\n",
                fontsize/2);
        fprintf(fptr,"{\\fldrslt\\f2\\fs%d}}}\n",fontsize);
        // END NEW.
        fprintf(fptr,"{\\b\\f2\\fs%d \\tab}", fontsize);
#endif
        strcpy(str,"\\b");
        for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
            if(gnull[i] == '_'){
                 state='x';
                 if(state == laststate) fprintf(fptr," ");
                 else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                 } laststate = state;
            } else {
                j++;
                assert(j <= length);
                char pval_char=PvalueToChar(1.0,-(3*RelEntropy[j]));
#if 0
                if(rtf2){       // use first rtf.
                  state = rtf2->ShowInfo(fptr,pval_char,laststate,str);
                } else {
                  state = ShowInfo(fptr,pval_char,laststate,str);
                }
#else
                state = ShowInfo(fptr,pval_char,laststate,str);
#endif
		laststate = state;
            }
        } fprintf(fptr,"}\n{\\f2\\fs%d \\tab }{\\f3 \n\\par }",fontsize);
}

Int4	rtf_typ::PutResEvalsXC(FILE *fptr,Int4 start,Int4 end,Int4 gstart,
	   Int4 gend,char *gnull,double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,
	   Int4 MaxConcensusLines,Int4 maxlen_gnull,Int4 *gnull_insrt_len, BooLean UnderLine)
// Put consensus residues for cross-conserved patterns: need to pass in sst_typ 
{
	char	*consensus[27];
	double	*cons_pval[27];
	Int4	index,high_index,max_index=25;
	double	*tmp_pval;
	Int4	factor=2,x,i,j;
	Int4	gapsize=fontsize;
	char	c,r,state,laststate,str[50];
	//============= First sort by Pvalues ================
	double sum_pvals=0.0;
	Int4	num_pvals=0,next_start=0;
	dh_type dH=dheap(max_index+2,4);
	char	**conserved=0;

print_error("PutResEvalsXC() not yet implemented");
	NEWP(conserved,length+3,char);
        for(i=gstart,j=start; i < gend && j < end; i++){
          if(gnull[i] != '_'){ j++; conserved[j]=rtf->Conserved(j); }
        }
	for(x=1; x <= max_index; x++){
	    	MEW(consensus[x],length+3,char); MEW(cons_pval[x],length+3,double);
	} high_index=0;
	for(i=gstart,j=start; i < gend && j < end; i++){
	    if(gnull[i] != '_'){
		j++;
	        tmp_pval=res_evals[j];
	        for(index=0,r=1; r<=nAlpha(AB); r++){
		   if(tmp_pval[r] <= -2.0){ sum_pvals+=tmp_pval[r]; num_pvals++; }
		   if(!wtfreq || (wtfreq[j] && wtfreq[j][r] > 0.1)) 
		      insrtHeap(r,((keytyp)(tmp_pval[r])),dH);
		   cons_pval[r][j] = 0.0; consensus[r][j] = ' ';
	        }
	        while((r=delminHeap(dH)) != 0){
			if(tmp_pval[r] > ExpPttrns) break;
			if((wtfreq[j] &&  wtfreq[j][r] < 0.1)) break;
			index++;
			cons_pval[index][j] = tmp_pval[r];
			consensus[index][j] = AlphaChar(r,AB);
			if(index >= max_index) break;
		}
		while(delminHeap(dH));	// empty heap
		if(index > high_index) high_index = index;
	    }
	} Nildheap(dH);
	double	target_pval=5.0;
	double	ave_pval=target_pval;
	if(num_pvals > 0){ ave_pval = -sum_pvals/(double)num_pvals; }
	factor = (Int4) floor((ave_pval/target_pval) + 0.5); 
	factor = MAXIMUM(Int4,1,factor);
	if(MaxPval <= 0.0){ factor = 1; }
	else { factor = (Int4) ceil((MaxPval/9.0)); factor = MAXIMUM(Int4,1,factor); }
	//============= output conserved residues ================
	char last_gnull=0;
	if(high_index > MaxConcensusLines){ high_index=MaxConcensusLines; } 
	if(high_index > 10){
		fprintf(stderr,"********* high_index=%d; ExpPttrns=%g *********\n",
								high_index,ExpPttrns);
	}
	for(index=1; index <= high_index; index++){
		char id_str[30]; 
		if(index==1) { sprintf(id_str,"%s :","pattern"); }
		else sprintf(id_str," ");
		fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d %s }",fontsize,color,id_str);
		fprintf(fptr,"{\\f2\\fs%d\\cf1 \\tab \\tab }", fontsize);
		strcpy(str,"\\b");
		last_gnull=0;
		for(laststate=0,i=gstart,j=start; i < gend && j < end; i++){
		   if(gnull[i] == '_'){
		     state='i';
		     c = ' ';
		     if(gnull_insrt_len && gnull_insrt_len[i] >= maxlen_gnull){
			if(last_gnull != '_'){		// first insert position
			  if(state != laststate){	
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i'){
			    	  fprintf(fptr,"{%s\\f3\\fs%d\\cf15 ",str,gapsize);  
			    } else {
				  fprintf(fptr,"{%s\\f3\\fs%d\\cf8 ",str,gapsize);  
			    }
			  }
			  Int4 numspace = 2 + (Int4) ceil(log10((double) gnull_insrt_len[i]+0.00001));
			  for(Int4 d = 1; d <= numspace; d++) fprintf(fptr," ");
			  laststate = state;
			}
		     } else {
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i'){
			    	  fprintf(fptr,"{%s\\f3\\fs%d\\cf15 %c",str,gapsize,c);  
			    } else {
				  fprintf(fptr,"{%s\\f3\\fs%d\\cf8 %c",str,gapsize,c);  
			    }
			  } laststate = state;
		     }
		   } else {
			j++;
			c = consensus[index][j]; r = AlphaCode(c,AB);
			if(c == ' '){
			  state='i';
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    fprintf(fptr,"{%s\\f3\\fs%d\\cf15 %c",str,gapsize,c);  
			  }
			} else {
			  if(UnderLine){
		            state='z';
			    if(state == laststate){ fprintf(fptr,"%c",c); }
                            else {
                              if(laststate != 0) fprintf(fptr,"}\n");
	     		      fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
			    }
			  } else {
			    state = rtf->ConservedState(conserved[j], r,c);
			    char tmp_str[3]; tmp_str[0]=state; tmp_str[1]=0;
			    if(strstr("xyu*",tmp_str)){ 
				  state='W'; // white on white...  
			    } rtf->PutConservedState(fptr,state,laststate,c,str);
			  }
			} laststate = state;
		   } last_gnull=gnull[i];
		} next_start=j;	// WARNING: TWO DISTINCT PATHS TO SET next_start...
		fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
		// if(index==1) fprintf(fptr,"\\tab %d",j);
		if(index==1){
			double log_two_fold=log10(pow(2.0,1.0/(1.0-rtf->GetLinearToLog())));
			fprintf(fptr,"\\tab %.1f",log_two_fold);
		} else fprintf(fptr,"\\tab ");
		fprintf(fptr,"}{\\f3 \n\\par }");
	}
	//============= output weighted residue frequencies ================
	for(x=1; x <= max_index; x++){ free(cons_pval[x]); free(consensus[x]); }
	for(j=1; j<=length; j++) if(conserved[j]) free(conserved[j]);
	free(conserved);
	return next_start;
}


