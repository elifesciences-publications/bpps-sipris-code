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

#include "rma_typ.h"

rma_typ::rma_typ(Int4 maxrpts,cma_typ cma)
// -G200..1300,20..100/0..50:500,35..250
{ init(maxrpts,"-P200..1300,20..100/0..50:500,35..250", cma); }

rma_typ::rma_typ(Int4 maxrpts,char *pssm_arg, cma_typ cma)
// input as "-P100..1100,20..120:500,40..400"
// or as -P150..1400,15..70/5..70:500,40..400  (Best settings??)
{ init(maxrpts,pssm_arg,cma); }

void	rma_typ::init(Int4 maxrpts,char *pssm_arg,cma_typ cma)
{
}

void	rma_typ::Put(FILE *fp)
{
}

void	rma_typ::Free( )
{
}

#if 0
char	FractionToCharRMA(double fract)
{
    char c;
    Int4 f = (Int4) ceil(10.0*fract);
    switch(f){
      case 0: c='0'; break;
      case 1: c='1'; break;
      case 2: c='2'; break;
      case 3: c='3'; break;
      case 4: c='4'; break;
      case 5: c='5'; break;
      case 6: c='6'; break;
      case 7: c='7'; break;
      case 8: c='8'; break;
      case 9: c='9'; break;
      case 10: c='!'; break;
      default: c='?'; break;
    } return c;
}
#endif

char	*ConservedRMA(double *observed, double cutoff, double *freq, a_type A,
	double *pvalue,double delta_p, double *res_evals)
/* return the conserved residues in column "observed". */
// scale is the difference 
{
	Int4	r,r2,nres,n;
	double	total,p,q,f,low_cut,high_cut,min_p;
	double	P,adjust;
	char	*conserved;
	BooLean	flag,debug=0,use_temp=FALSE;

	assert(cutoff <= 0.0);
	low_cut=cutoff;
	cutoff = low_cut - delta_p;  // reset from minimum to middle value
	high_cut=cutoff - delta_p;
	*pvalue=min_p=0.0;
	if(observed != NULL){
	    nres=0;
	    NEW(conserved,nAlpha(A)+2,char); // conserved[r]=0; ---> 'u';
	    for(total=0.0, r=1; r<=nAlpha(A); r++){ total += (double)observed[r]; }
	    if(debug){
		fprintf(stderr,"cutoff = %.2f; low_cut = %.2f\n",cutoff,low_cut);
		for(r = 1; r <= nAlpha(A); r++){
		   if(observed[r] > 0){
			q = freq[r];
			fprintf(stderr,"%c: %d/%d (p=%.1f)\n", 
				AlphaChar(r,A), (Int4) observed[r], (Int4)total,
				Log10CBP(observed[r], total, q));
   		   }
		} fprintf(stderr,"\n");
	    }
	    /** 1. find clearly conserved residues **/
	    for(f=0.0, n=0, r=1; r <= nAlpha(A); r++){
		if(observed[r] > 2){
		   q = freq[r];
		   p = Log10CBP(observed[r], total, q);
		   if(res_evals) res_evals[r]=p;
		   P = p;
		   if(p < min_p) min_p = p;
	           if(P <= low_cut) {
			if(debug) fprintf(stderr," %c: (%.1f)[%d/%d]", 
				AlphaChar(r,A), 
				p,(Int4) observed[r], (Int4)total);
			n += observed[r]; f += q; 
			if(P <= cutoff){
		   		if(P <= high_cut) conserved[r] = 's'; 
				else  conserved[r] = 'm';
			} else conserved[r] = 'w';
			nres++;
		   }
		} else if(res_evals) res_evals[r]=0.0;
	    }
            /** 2. find weakly conserved related residue pairs **/
            for(n=0, r=1; r < nAlpha(A); r++){
                if(!conserved[r] && observed[r] > 2){
                  for(r2 = 1; r2 <= nAlpha(A); r2++){
                        if(!conserved[r2] && r2 != r && observed[r2] > 2 
                                        && valAlphaR(r,r2,A) > 0){
                           q = freq[r] + freq[r2];
		  	   p = Log10CBP(observed[r] + observed[r2], total, q);
		   	   P = p;
			   if(p < min_p) min_p=p;
	                   if(P <= low_cut) {
			     if(debug) fprintf(stderr," %c %c (%.1f)", 
					AlphaChar(r,A), AlphaChar(r2,A), p);
                             nres+=2;
			     n += observed[r] + observed[r2];
			     f += freq[r] + freq[r2];
			     if(P <= cutoff){
		   		if(P <= high_cut){
                                  conserved[r]='s'; conserved[r2]='s';
				} else { conserved[r]='m'; conserved[r2]='m'; }
			     } else { conserved[r]='w'; conserved[r2]='w'; }
                           }
                        }
                  }
                }
            } 
	    /** 3. find weakly conserved residues related to clear ones **/
            do {
	      flag = FALSE;
	      for(r=1; r <= nAlpha(A); r++){
		if(conserved[r]){
	          for(r2 = 1; r2 <= nAlpha(A); r2++){
		   if(!conserved[r2] && observed[r2] >= 2){
			if(valAlphaR(r,r2,A) >= 0){
		   	   q = freq[r2];
		           p = Log10CBP(observed[r2], (total-n), q/(1.0 - f));
		   	   P = p;
			   if(p < min_p) min_p=p;
	                   if(P <= low_cut) {
				if(debug) fprintf(stderr," %c (%.1f)", 
						AlphaChar(r2,A), p);
				flag = TRUE; nres++;
				conserved[r2] = 'w';
			   }
			}
		   }
		  }
		}
	      }
	    } while(flag);
	    if(nres==0){ free(conserved); return NULL; }
#if 0
	    // 3. reset to total conserved residues if more highly conserved
	    for(q=0.0,r2=0, r = 1; r <= nAlpha(A); r++){
	      if(conserved[r]){ r2 += observed[r]; q += freq[r]; }
	    }
	    p = Log10CBP(r2, total, q);
	    if(use_temp) P = pow(p,temperature); else P = p;
	    if(P <= low_cut) {
	      for(r=1; r <= nAlpha(A); r++){
		if(conserved[r]){
	   	  if(P <= cutoff){
		     if(P <= high_cut) conserved[r] = 's';
		     else if(conserved[r] != 's') conserved[r] = 'm';
		  } // else if(conserved[r] != 'm' && conserved[r] != 's') conserved[r] = 'm';
		  // else conserved[r] = 'w';
		}
	      }
	    }
#endif
	    *pvalue = min_p;
	    return conserved;
	} else return NULL;
}

float	**GetInfoCMSA(double purge_cut, double ***freq, sma_typ MA, cma_typ cma,
	char mode, double ***Observed)
// Create an alignment figure for publication:
{
	Int4	t,n,s,length;
	double	*observed,info,sub_info,total,f,q;
	float	**RtnInfo;
	a_type	A=AlphabetCMSA(cma);
	char	r;
	BooLean	*use;

    use=purge_sma_cma(purge_cut,MA,cma);
    NEW(observed,nAlpha(A)+2,double);
    NEWP(RtnInfo,ntypSMA(MA)+2,float);
    for(t=1; t <= ntypSMA(MA); t++){
	length = lengthSMA(t,MA);
        NEW(RtnInfo[t],length+2,float);
	for(s=1; s <= length; s++){
	    if(Observed){
	      if(Observed[t][s]){
	        for(r=0; r<=nAlpha(A); r++) observed[r]=Observed[t][s][r];
	      } else for(r=0; r<=nAlpha(A); r++) observed[r]=0.0;
	    }else {
	      for(r=0; r<=nAlpha(A); r++) observed[r]=0.0;
	      for(n=1; n<=nseqSMA(MA); n++){
		if(use[n]){ r=residueSMA(t,n,s,MA); observed[r]+=1.0; }
	      }
	    } // use column freq from alternative cma here for subfamily models...

	    if(Observed){
              for(total=0.0,r=1; r <= nAlpha(A); r++){ 
		total += (double) observed[r]+1.0; // fully uninformed prior...
	      }
	    } else {
              for(total=0.0,r=1; r <= nAlpha(A); r++){ 
		total += (double) observed[r]+nAlpha(A)*freqSMA(r,MA); 
	      }
	    }
            for(sub_info=info=0.0,r=1; r <= nAlpha(A); r++){
		if(observed[r] > 0){
		   if(Observed){
			f = (observed[r]+1.0)/total;
                        q = freqSMA(r,MA);
                        if(q > 0.0) info += f*log(f/q);
			if(freq[t][s] != 0){ q = freq[t][s][r]; }
                        if(q > 0.0) sub_info += f*log(f/q);
		   } else {
                        f = (observed[r]+nAlpha(A)*freqSMA(r,MA))/total;
                        // f = (double) (observed[r]+1)/total;
                        q = freqSMA(r,MA);
                        if(q > 0.0) info += f*log(f/q);
			if(freq[t][s] != 0){ q = freq[t][s][r]; }
			if(q > 0.0) sub_info += f*log(f/q);
		   }
                }
	    } info /= log(2.0); sub_info /= log(2.0);
	    RtnInfo[t][s]=sub_info;
	}
    } free(observed); free(use);
    return RtnInfo;
}

float	***SubtractInfoCMSA(double purge_cut, Int4 rpts, char *pssm_arg, cma_typ main_cma, 
	cma_typ sub_cma, sma_typ subMA, char mode,Int4 JackCut)
// mode = 'n' for net (subtract main from sub); 'm' for main; 's' for subfamily
{
	// assert(nBlksCMSA(main_cma) == 1); // just one block for right now...
	assert(nBlksCMSA(sub_cma) == 1); // just one block for now...

	Int4	s;
	double	***Observed=0;
	double	***freq;	// freq[blk][col][res];
	float	***rtninfo;

     	NEWPP(rtninfo,5,float);
	freq=SuperModelFreqsCMSA(rpts,pssm_arg,sub_cma,main_cma,&Observed,JackCut);

     	// 1. Find alignment relating main family to subfamily via consensus.
	// net information...
	rtninfo[1]=GetInfoCMSA(purge_cut,freq,subMA,sub_cma,mode,0);
	for(s=1; s <=LengthCMSA(1,sub_cma); s++){ if(freq[1][s]) free(freq[1][s]); }
	free(freq[1]); free(freq);

	NEWPP(freq,4,double);
	NEWP(freq[1],LengthCMSA(1,sub_cma)+3,double);
	rtninfo[2]=GetInfoCMSA(purge_cut,freq,subMA,sub_cma,mode,0);
	rtninfo[3]=GetInfoCMSA(purge_cut,freq,subMA,sub_cma,mode,Observed);
	for(s=1; s <=LengthCMSA(1,sub_cma); s++){
		if(Observed[1][s]) free(Observed[1][s]); 
	} free(Observed[1]); free(Observed);
	free(freq[1]); free(freq);
	return rtninfo;
}

double  ***SuperModelFreqsCMSA(Int4 rpts, char *pssm_arg,cma_typ cma,
	cma_typ main_cma,double ****Observed,Int4 JackCut)
#if 0	//*********************************************************************
	Get position-specific residue frequencies based on a superfamily
	alignment (main_cma) for use in a subfamily alignment (cma).
	If Observed != NULL then also return the observed number of 
	each type of residue at each position (Observed[t][s][r]).
	NOTE: This code only works with one subfamily block (t==1) right now...
#endif	//*********************************************************************
{
        double  ***freq;
	a_type	A=AlphabetCMSA(cma);

	assert(TotalLenCMSA(cma) >= TotalLenCMSA(main_cma));
	// assert(nBlksCMSA(main_cma) == 1); // just one block for right now...
	assert(nBlksCMSA(cma) == 1); // just one block for right now...

	// 1. Get a consensus sequence for the subfamily alignment
	e_type  E=MkConsensusCMSA(cma);	// implemented....

	// 2. Align consensus against the main model...
	// 2a. Get rma_typ for alignment
        Int4 oper_len,score,start;
        char *operation;
	HMM_typ *hmm=new HMM_typ(rpts,pssm_arg,main_cma,200,0);
        operation=hmm->Align(stderr,E,1,&score,&start,&oper_len);

// std::cerr << operation; std::cerr << std::endl;

	// ========== 5. Create a gapped sequence. ===========
	Int4    *newpos,s,b;
	NEW(newpos,rpts*nBlksCMSA(main_cma)+2,Int4);  // for new sites.
        gsq_typ *gsq; gsq = new gsq_typ[1];
	
        gsq->initialize(5,5,operation,oper_len,start,E,newpos);
	// gsq->Put(stderr, A);
	e_type fakeE = gsq->FakeSeq( );

	Int4	lenM,begin;
	double	***mfreq;
	BooLean *skip=0;
	double	***observed,***Obs;

	unsigned char *seq=SeqPtr(fakeE);
        NEWPP(freq,4,double);
        NEWP(freq[1],LengthCMSA(1,cma)+3,double); // Leave null at insertions.
	NEWPP(mfreq,nBlksCMSA(main_cma)+3,double);
	NEWPP(observed,nBlksCMSA(main_cma)+3,double);
	if(Observed){ NEWPP(Obs,4,double); NEWP(Obs[1],LengthCMSA(1,cma)+3,double); }
	for(Int4 r=1; r <= rpts; r++){
	  begin = (r-1)*nBlksCMSA(main_cma);
	  fprintf(stderr,"newpos[%d] = %d\n",r,newpos[begin+1]);
	  if(Observed){
	    if(JackCut < 100){
	      skip=FindCloseCMA(JackCut,fakeE,newpos+begin,main_cma);
	      for(b=1; b <= nBlksCMSA(main_cma); b++){	
	      	mfreq[b]=ColResFreqsCMSA(b,skip,&observed[b],main_cma);
	      } free(skip);
	    } else {
	      for(b=1; b <= nBlksCMSA(main_cma); b++){	
		mfreq[b]=ColResFreqsCMSA(b,&observed[b],main_cma);
	      }
	    }
	  } else {
	    if(JackCut < 100){
	      skip=FindCloseCMA(JackCut,fakeE,newpos+begin, main_cma);
	      for(b=1; b <= nBlksCMSA(main_cma); b++){	
	         mfreq[b]=ColResFreqsCMSA(b,skip,main_cma);
	      } free(skip);
	    } else {
	      for(b=1; b <= nBlksCMSA(main_cma); b++){	
		mfreq[b]=ColResFreqsCMSA(b,main_cma);
	      }
	    }
	  }
	  for(b=1; b <= nBlksCMSA(main_cma); b++){
	    begin = (r-1)*nBlksCMSA(main_cma) + b;
	    lenM=LengthCMSA(b,main_cma);
            for(s=1; s <=lenM; s++){
	     Int4 fakesite = newpos[begin]+s-1;
	     if(!gsq->IsDeleted(fakesite)){	// Skip over deletions.
	        Int4 truesite=gsq->FakeToReal(fakesite);
		if(Observed){ Obs[1][truesite]=observed[b][s]; }
		freq[1][truesite]=mfreq[b][s];   
		double	totfreq=0.0;
		for(Int4 res=0; res <=nAlpha(A); res++) totfreq+=mfreq[b][s][res];
#if 0
	  	fprintf(stderr," %c: true = %d; fake = %d; s = %d; totfreq=%.3f\n",
			AlphaChar(seq[fakesite],A),truesite,fakesite,s,totfreq);
#endif
	     } else { free(mfreq[b][s]); if(Observed) free(observed[b][s]); }
	   }
	  }
	}
	for(b=1; b <= nBlksCMSA(main_cma); b++){
		free(mfreq[b]); if(Observed) free(observed[b]);
	}
	free(mfreq); mfreq=0; if(Observed) free(observed); observed=0;
	if(Observed){ *Observed = Obs; }
	free(newpos); free(operation); delete []gsq; NilSeq(E);
	delete hmm;
	return freq;
}

//*************************** RTF OUTPUT ROUTINES ******************************

void	NewCMA2RTF_HEADER(FILE *fptr,char PageSetUp)
{ NewCMA2RTF_HEADER(fptr,PageSetUp,12); }

void	NewCMA2RTF_HEADER(FILE *fptr,char PageSetUp,Int4 fontsize)
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
    fprintf(fptr,"\\red255\\green0\\blue255;\\red255\\green0\\blue0;");
    fprintf(fptr,"\\red255\\green255\\blue0;\\red255\\green255\\blue255;");
    fprintf(fptr,"\\red0\\green0\\blue128;\\red0\\green128\\blue128;");
    fprintf(fptr,"\\red0\\green128\\blue0;\\red128\\green0\\blue128;");
    fprintf(fptr,"\\red128\\green0\\blue0;\\red128\\green128\\blue0;");
    fprintf(fptr,"\\red128\\green128\\blue128;\\red192\\green192\\blue192;}");
    /** 15 dk grey; 16 lt grey (original) **/
    fprintf(fptr,"{\\stylesheet{\\widctlpar \\f4\\fs%d \\snext0 Normal;}",fontsize);
    fprintf(fptr,"{\\*\\cs10 \\additive Default Paragraph Font;}}\n");

    switch(PageSetUp){
	  case 'l': // 8.5" x 11" landscape.
	   fprintf(fptr,"\\paperw15840\\paperh12240\\margl720\\margr720");
    	   fprintf(fptr,"\\margt1080\\margb1080 ");
    	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
    	   // page_width=12000.0;
	   break;
	  case 'L': // 11" x 17" landscape
	   fprintf(fptr,"\\paperw24480\\paperh15840\\margl720\\margr720");
	   fprintf(fptr,"\\margt1080\\margb1080 ");
	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
	   // print_error("NewCMA2RTF_HEADER( ) 11\" x 17\" landscape not yet implemented");
	   break;
	  case 'p': // 8.5" x 11" portrait.
	   // page_width=7000.0;
	   break;
	  case 'P': // 11" x 17" portrait.
	    fprintf(fptr,"\\paperw15840\\paperh24480\\margl720\\margr720");
            fprintf(fptr,"\\margt1080\\margb1080 ");
            fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
            fprintf(fptr,"\\linex0\\endnhere\\sectdefaultcl \n");
            // page_width=12000.0;
	   break;
	  default: print_error("NewCMA2RTF_HEADER( ) PageSetUp input error");
    }
    fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs%d \n",fontsize);
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
}

char get_conserved_state(char *conserved, char r,char c)
{
	char state;

	if(conserved==NULL || !conserved[r]) state='u';
	else {
	  switch (conserved[r]){
	    case 'w': 	// weak...
		switch(c){
		  case 'G': case 'P': state = 'x'; 
		    break;
		  case 'C': 
		  case 'A': case 'I': case 'L': case 'V': case 'M': 
		  case 'F': case 'W': case 'Y': 
			state = 'y'; 
		    break;
		  case 'H': case 'K': case 'R': 
		  case 'N': case 'Q': case 'S': case 'T': 
		  case 'D': case 'E': state = '*';
		    break;
		  default: print_error("This should not happen!");
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
		  default: print_error("This should not happen!");
		   break;
		}
		if(conserved[r]=='s') state = toupper(state);
	     break;
	    default: print_error("This should not happen!");
	     break;
	  }
	} return state;
}

void	print_conserved_state(FILE *fptr, char state, char laststate, 
			char c, char *str)
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
unconserved: 'u';
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
#endif	//******************************************************************
{
	if(state == laststate) fprintf(fptr,"%c",c);
	else {
	   if(laststate != 0) fprintf(fptr,"}\n");
	   switch(state){
	     case 'G': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight4 %c",str, c);
		 break;  // white on green
	     case 'P': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight1 %c",str, c);
		 break;  // white on black

	     case 'C': 
		fprintf(fptr,"{%s\\f3\\cf2\\highlight7 %c",str, c);
		 break;  // blue on yellow
	     case 'N': 
		fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str, c);
		 break;  // red on yellow
	     case 'A': 
		fprintf(fptr,"{%s\\f3\\cf5\\highlight7 %c",str, c);
		 break;  // magenta on yellow

	     case 'H': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight2 %c",str, c);
		 break;  // white on blue
	     case 'B': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight3 %c",str, c);
		 break;  // white on cyan
	     case 'O': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight5 %c",str, c);
		 break;  /** white on magenta **/
	     case 'D': 
		fprintf(fptr,"{%s\\f3\\cf8\\highlight6 %c",str, c);
		 break;  /** white on red **/

	     case 'g':
		  fprintf(fptr,"{%s\\f3\\cf4 %c",str,c);
		break;  /** green on white **/
	     case 'p':
		  fprintf(fptr,"{%s\\f3\\cf1 %c",str,c);
		break;  /** black on white **/

	     case 'c': case 'n': case 'a':
		fprintf(fptr,"{%s\\f3\\cf14\\highlight7 %c",str, c);
		 break;  /** dk yellow on yellow **/

	     case 'h':
		  fprintf(fptr,"{%s\\f3\\cf2 %c",str,c);
		break;  /** blue on white **/
	     case 'b':
		  fprintf(fptr,"{%s\\f3\\cf3 %c",str,c);
		break;  /** cyan on white **/
	     case 'o':
		  fprintf(fptr,"{%s\\f3\\cf5 %c",str,c);
		break;  /** magenta on white **/
	     case 'd':
		  fprintf(fptr,"{%s\\f3\\cf6 %c",str,c);
		break;  /** red on white **/

	     case 'y': case '*':
	     case 'x':
		  fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
		break;  /** dk grey on white **/
#if 0 
	     case 'y':
		  fprintf(fptr,"{%s\\f3\\cf14 %c",str,c);
		break;  /** dk yellow on white **/
	     case '*':
		  fprintf(fptr,"{%s\\f3\\cf12 %c",str,c);
		break;  /** violet on white **/
#endif
	     case 'W':
		  fprintf(fptr,"{%s\\f3\\cf8 %c",str,c);
		 break;  /** white on white **/
	     case 'u':
	     default : 
		  fprintf(fptr,"{%s\\f3\\cf16 %c",str,c);
		 break;  /** light grey on white **/
	   } 
	}
}

char	PvalueToCharRMA(double factor, double pval)
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

double *ComputRelEntropyRMA(Int4 length, char **conserved, char *info_show, Int4 t,
	double *observed, double ***Observed, double ***freq, sma_typ MA, a_type A,
	double cbp_cut,BooLean *use,char *pval_show, double *binomial_tail,
	double delta_p,double **ResEvals) 
/******* compute relative entropy of postions *******/
{
	Int4	s,r,tmp_d,n; 
	double	total,sub_info,info,d,f,q;
	double	*info_real,pval,*res_evals;
        h_type HG2=Histogram("Kyte-Doolittle hydrophobicity",0,100,0.2);
	NEW(info_real,length+3,double);
	for(s=1; s <= length; s++){
	    if(Observed){
	      if(Observed[t][s]){
	        for(r=0; r<=nAlpha(A); r++) observed[r]=Observed[t][s][r];
	      } else for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	    } else {
	      for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	      for(n=1; n<=nseqSMA(MA); n++){
		if(use[n]){ r=residueSMA(t,n,s,MA); observed[r]++; }
	      }
	    }
	    // use column freq from alternative cma here for subfamily models...
	    if(ResEvals) res_evals = ResEvals[s]; else res_evals = 0;
	    if(freq[t][s] != 0)
	      conserved[s]=ConservedRMA(observed,cbp_cut,freq[t][s],A,&pval,
			delta_p,res_evals);
	    else conserved[s]=ConservedRMA(observed,cbp_cut,FreqSMA(MA),A,&pval,
			delta_p,res_evals);
	    double factor=0.5;
	    pval_show[s] = PvalueToCharRMA(factor,pval);
	    binomial_tail[s] = -pval;

	    double kd_score = 0.0;
	    tmp_d=0;
            // for(total=0.0,r=1; r <= nAlpha(A); r++){ total += observed[r]; }
	    if(Observed){
              for(total=0.0,r=1; r <= nAlpha(A); r++){ 
		total += observed[r]+1.0; // uninformed prior...
	      }
	    } else {
              for(total=0.0,r=1; r <= nAlpha(A); r++){ 
		// total += observed[r]+freqSMA(r,MA); 
		total += observed[r]+nAlpha(A)*freqSMA(r,MA); 
	      }
	    }
	    if(total < 4){ info=sub_info=0.0; }
	    else {
// double total_q=0.0,total_f=0.0;
              for(sub_info=info=0.0,r=1; r <= nAlpha(A); r++){
		// if(observed[r] > 0.0){
			kd_score += observed[r]*blsm62kyte_doolittle[r];
		   if(Observed){
			f = (observed[r]+1.0)/total;
                        q = freqSMA(r,MA);
                        if(q > 0.0) info += f*log(f/q);
			if(freq[t][s] != 0){ q = freq[t][s][r]; }
                        if(q > 0.0) sub_info += f*log(f/q);
#if 0
fprintf(stdout,"%5d  %c  %.4f %.4f\n",s,AlphaChar(r,A),f,q);
	total_q+=q; total_f+=f;
	assert(f >= 0.0 && f <= 1.0); assert(q >= 0.0 && q <= 1.0);
#endif
		   } else {
                        f = (observed[r]+nAlpha(A)*freqSMA(r,MA))/total;
                        // f = (observed[r]+1)/total;
                        q = freqSMA(r,MA);
                        if(q > 0.0) info += f*log(f/q);
			if(freq[t][s] != 0){ q = freq[t][s][r]; }
			if(q > 0.0) sub_info += f*log(f/q);
#if 0
fprintf(stdout,"%5d  %c  %.4f %.4f\n",s,AlphaChar(r,A),f,q);
	total_q+=q; total_f+=f;
	assert(f >= 0.0 && f <= 1.0); assert(q >= 0.0 && q <= 1.0);
#endif
		   }
                // }
              } // info /= log(2.0); sub_info /= log(2.0);
#if 0
	assert(total_f > 0.999 && total_f  < 1.001);
	assert(total_q > 0.999 && total_q  < 1.001);
#endif
	    }
	    info_real[s]=sub_info;
	    assert(sub_info >= 0.0);
	    tmp_d = (Int4) floor((5*(sub_info)));
	    // if(freq[t][s] == 0){ info_show[s]='-'; } 
	    // else if(Observed && Observed[t][s]==0){ info_show[s]='-'; }
	    if(Observed && Observed[t][s]==0){ info_show[s]='-'; }
	    else if(tmp_d >= 10) info_show[s] = '!';
	    else {
	      switch(tmp_d){
		case 1: info_show[s] = '1'; break;
		case 2: info_show[s] = '2'; break;
		case 3: info_show[s] = '3'; break;
		case 4: info_show[s] = '4'; break;
		case 5: info_show[s] = '5'; break;
		case 6: info_show[s] = '6'; break;
		case 7: info_show[s] = '7'; break;
		case 8: info_show[s] = '8'; break;
		case 9: case 10: info_show[s] = '9'; break;
		default: info_show[s] = '0'; break;
	      }
	    }
#if 0
fprintf(stderr,"info[%d][%d] = %.2f; ",t,s,info);
fprintf(stderr,"sub_info = %.2f; ",sub_info);
fprintf(stderr,"diff = %.1f --> %c (%.1f)\n-----\n",10.*(info-sub_info),info_show[s],total);
#endif
	    kd_score /= total;
            IncdHist(kd_score,HG2); 
	}
	if(!Observed) ; // PutHist(stderr,60,HG2);  
	NilHist(HG2);
	return info_real;
}

char	put_info_rma_typ(FILE *fptr, char info_show, char laststate, char *str)
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
{
	char	state;
	switch (info_show){
	  case '0':  info_show=' ';
	  case '1': case '2': state='g'; break;
	  case '3': case '4': state='G'; break;
	  case '5': case '6': state='B'; break;
	  case '7': case '8': state='r'; break;
	  case '9': case '!': state='R'; break;
	  case ' ': state='u'; break;
	  case '-': state='u'; break;
	  default: state='u'; break;
	}
	if(state == laststate) fprintf(fptr,"%c",info_show);
	else {
	   if(laststate != 0) fprintf(fptr,"}\n");
	   switch (state) {
	     case 'w': // white
	        fprintf(fptr,"{%s\\f3\\cf8 %c",str,info_show);
		break;
	     case 'g': // light grey
	        fprintf(fptr,"{%s\\f3\\cf16 %c",str,info_show);
		break;
	     case 'G': // grey
	        fprintf(fptr,"{%s\\f3\\cf15 %c",str,info_show);
		break;
	     case 'b': // blue
	        fprintf(fptr,"{%s\\f3\\cf2 %c",str,info_show);
		break;
	     case 'B': // black
	        fprintf(fptr,"{%s\\f3\\cf1 %c",str,info_show);
		break;
	     case 'r': // red
	        fprintf(fptr,"{%s\\f3\\cf6 %c",str,info_show);
		break;
	     case 'R':  
		fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str,info_show);
		 break;  // red on yellow
	     case ' ': 
	     default:
	        fprintf(fptr,"{%s\\f3\\cf1 %c",str,info_show);
		break;
	   }
	} return state;
}

void	NewCMA2RTF_TAIL(FILE *fptr)
{ fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1  {\\f3\n\\par}}\n"); }

void    NewCMA2RTF_RETURNS(FILE *fptr,unsigned short numRtns)
{
	fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3");
	for(unsigned short i=1; i <=numRtns; i++) fprintf(fptr,"\n\\par");
	fprintf(fptr," }\n");
}

void	NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, double ***freq, sma_typ MA,
	cma_typ cma,double ***Observed,char ColorCode)
{ NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,freq,MA,
		cma,Observed,0,'X',ColorCode); }

void	put_rasmol_rma(FILE *fp_ras, sma_typ MA, Int4 t,char su, a_type A, char *info_show,
	double *info_real)
// Rasmol output file option...
{
	char	*gnull = gnullSMA(t,MA);
        const char    *residue[ ] = { "nil", /** X CGASTNDEQK RHWYFVILMP **/
                "cys","gly","ala","ser","thr","asn","asp","glu","gln","lys",
                "arg","his","trp","tyr","phe","val","ile","leu","met","pro"};
	char	*resname;
	Int4    stickwidth=100,i,j;
	unsigned char	*sq;
	char	sidechain_str[200],*sidechain_color,*gseq;
	Int4	s0,s,r;

	if(!fp_ras) return;
	assert(t == 1);
	s = startSMA(t,1,MA);
	sq = seqSMA(t,1,MA);
	gseq = gseqSMA(t,1,MA);
	for(s0=1,i=j=0; j < lengthSMA(t,MA); i++){
	 if(gnull[i] != '_'){
	  j++;
       	  if(su == ' ') fprintf(fp_ras,"select %d\n",s);
       	  else fprintf(fp_ras,"select %d%c\n",s,su);
       	  fprintf(fp_ras,"strands off\n");
       	  fprintf(fp_ras,"color [220,220,255]\n");
       	  fprintf(fp_ras,"cartoons\n");
	  sidechain_color=sidechain_str;
	  unsigned char Red=255,Blue,Green;
	  switch (info_show[j]){
 	    case '0': sidechain_color=0; break;
	    case '1': sidechain_color=0; break;
	    // case '1': Blue=225; Green=225; break;
	    case '2': sidechain_color=0; break;
	    // case '2': Blue=200; Green=200; break;
	    case '3': sidechain_color=0; break;
	    // case '3': Blue=175; Green=175; break;
	    case '4': Blue=150; Green=150; break;
	    case '5': Blue=125; Green=125; break;
	    case '6': Blue=100; Green=100; break;
	    case '7': Blue=75; Green=75; break;
	    case '8': Blue=50; Green=50; break;
	    case '9': Blue=25; Green=25; break;
	    case '!': Blue=0; Green=0; break;
	    case ' ': case '-': sidechain_color=0; break;
	    default: sidechain_color=0; break;
	  }
	  r = sq[s0];
	  assert(r <= 20);
	  if(sidechain_color){
	   sprintf(sidechain_str,"[%d,%d,%d]",Red,Blue,Green);
	   resname = (char *) residue[r];
	          char res_char = AlphaChar(r,A);
        	  if(su!=' '){
                    if(res_char == 'G'){
                       fprintf(fp_ras, "select %s%d%c.ca\nspacefill 400\n",
                               resname,s,su);
                    } else {
                       fprintf(fp_ras,
                          "select %s%d%c.ca,(%d%c and sidechain)\nwireframe %d\n",
                               resname,s,su,s,su,stickwidth);
                       fprintf(fp_ras,"select %d%c and sidechain\n", s,su);
                    }
                  } else {
                    if(res_char == 'G'){
                       fprintf(fp_ras, "select %s%d.ca\nspacefill 400\n",
                               resname,s);
                    } else {
                       fprintf(fp_ras,"select %s%d.ca,(%d and sidechain)\nwireframe %d\n",
                               resname,s,s,stickwidth);
                       fprintf(fp_ras,"select %d and sidechain\n", s);
                    }
                  } fprintf(fp_ras,"color %s\n",sidechain_color);
	  } s0++;
	  if(!(r==0 && gseq[i] == '-')) s++;
	 }
	}
}

void	NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, double ***freq, sma_typ MA, cma_typ cma,
	double ***Observed,FILE *fp_ras,char su,char ColorCode)
{ NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,freq,MA,
		cma,Observed,fp_ras,su,0,0,ColorCode); }

Int4	get_font_color_code2(char ColorCode)
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow; 
// 15 dk grey; 16 lt grey...
{
	Int4	color;

    switch (ColorCode){		// get color for sequence name...
	case 'B': color = 2; break;	// blue
	case 'b': color = 9; break;	// dk blue
	case 'C': color = 3; break;	// cyan
	case 'c': color = 10; break;	// teal
	case 'G': color = 4; break;	// green
	case 'g': color = 11; break;	// dk green
	case 'M': color = 5; break;	// magenta
	case 'm': color = 12; break;	// violet
	case 'Y': color = 7; break;	// yellow
	case 'y': color = 14; break;	// dk yellow
	case 'R': color = 6; break;	// red
	case 'r': color = 13; break;	// dk red
	case 'L': color = 16; break;	// lt grey
	case 'D': color = 15; break;	// dk grey
	case 'W': color = 8; break;	// white
	default: color = 1; break;       // black
    } return color;
}

void	NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, double ***freq, sma_typ MA, cma_typ cma,
	double ***Observed,FILE *fp_ras,char su,double **MargProb,
	char *title,char ColorCode)
// Create an alignment figure for publication:
{
	Int4	t,i,j,n,s,N,length,gap,ntyp;
	double	*observed;
	char	*pval_show,*info_show,r,c,str[50],state,laststate,**conserved;
	a_type	A=AlphabetCMSA(cma);
	BooLean	*use;
    	Int4	color=get_font_color_code2(ColorCode);
	double	*binomial_tail;

    assert(cbp_cut < 0.0 && infoLO >= 0.1 && infoLO < infoHI);
    ntyp = ntypSMA(MA);
    if(Observed != 0) use=purge_sma_cma(purge_cut,MA,cma);
    else use=purge_sma_cma(200.0,MA,cma);  // don't remove any from orthologs
    NEW(observed,nAlpha(A)+2,double);
    if(title != 0){
	fprintf(fptr,"{\\f2\\fs16 %s\\tab }{\\f3 \n\\par }",title);
	// fprintf(fptr,"%s",title);
	// fprintf(fptr,"\\tab }{\\f3 \n\\par }");
    }
    for(t=1; t <= ntyp; t++){
	length = lengthSMA(t,MA);
	NEWP(conserved, length+2, char);
	NEW(info_show, length+2, char); NEW(pval_show, length+2, char); 
	NEW(binomial_tail,length+2,double);

        double *info_real=ComputRelEntropyRMA(length,conserved,info_show,t,observed,
		Observed,freq,MA,A,cbp_cut,use,pval_show,binomial_tail,2.0,0);
	free(binomial_tail);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
	BooLean	*show_position=0;
	if(key_seq){
	  NEW(show_position,glengthSMA(t,MA) + 2, BooLean);
	  char *gnull = gnullSMA(t,MA);
	  for(i=j=0; j < lengthSMA(t,MA); i++){
	    if(gnull[i] != '_'){ j++; show_position[i]=TRUE; }
	    else for(n=1; n<= nseqSMA(MA); n++){
	      char *gseq = gseqSMA(t,n,MA);
    	      if(IsInSeqSet(key_seq,seqidSMA(n,MA))){
		c = gseq[i];
		if(c != '.') show_position[i] = TRUE;
	      }
	    }
	  }
	}
	char *gnull,*gseq;
	for(N=0, n=1; n<= nseqSMA(MA); n++){
    		if(key_seq && !IsInSeqSet(key_seq,seqidSMA(n,MA))) continue;
#if 0	// new option to show lowercase kingdom sequences only.
    		e_type trueE=TrueSeqCMSA(n,cma);
		fprintf(stderr,"seq n kingdom=%c\n",kingdomSeq(trueE));
		PutSeqInfo(stderr,trueE);
		if(isupper(kingdomSeq(trueE))) continue;
#endif
	   	if(t==1){
		   fprintf(fptr,"{\\f2\\fs16\\cf%d ",color);
		   fprintf(fptr,"%s",seqidSMA(n,MA));
		   fprintf(fptr,"\\tab %d\\tab }",startSMA(1,n,MA));
	   	}
		if(t < ntypSMA(MA)) {
			gap = startSMA(t+1,n,MA) - endSMA(t,n,MA) - 1;
		}
		if(IsProbSMA(MA)){
		  // if(probSMA(n,t,MA) < 0.0) strcpy(str,"\\b\\i\\strike"); 
		  // else if(probSMA(n,t,MA) < 1.3) strcpy(str,"\\b\\i");
		  if(probSMA(n,t,MA) < 1.3) strcpy(str,"\\b\\i");
		  else strcpy(str,"\\b");
		} else strcpy(str,"\\b");
		gnull = gnullSMA(t,MA);
	        gseq = gseqSMA(t,n,MA);
		unsigned char *seq = seqSMA(t,n,MA);
		for(laststate=0,i=j=0; j < lengthSMA(t,MA); i++){
		   if(gnull[i] == '_'){
			if(!show_position || show_position[i]){
			  if(Observed && !Observed[t][j]) state='w';
			  else state='i';
			  c = gseq[i];
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    if(state=='i')
			      fprintf(fptr,"{%s\\f3\\fs14\\cf16 %c",str,c);  // gap color
			    else fprintf(fptr,"{%s\\f3\\fs14\\cf8 %c",str,c);  // gap color
			  } laststate = state;
			}
		   } else {
			j++;
			r = seq[j]; 
			if(r==0 && gseq[i] == '-') c = '-';
			else c = AlphaChar(r,A); 
			if(Observed && !Observed[t][j]) state='W';
			else state = get_conserved_state(conserved[j], r,c);
			print_conserved_state(fptr, state, laststate, c, str);
			laststate = state;
		   }
		}
		fprintf(fptr,"}\n{\\f2\\fs16 ");
		if(t==ntyp) fprintf(fptr,"\\tab %d", endSMA(t,n,MA));
		else if(gap <=0) fprintf(fptr,"\\tab   ");
		else fprintf(fptr,"\\tab (%d)", gap);
#if 1
		if(IsProbSMA(MA)) {
			fprintf(fptr,"\\tab [%.1f]",probSMA(n,t,MA));
		}
#endif
		fprintf(fptr,"\\tab }{\\f3 \n\\par }");
	}  /** end sequences **/
	//**************** output information... *******************
	if(t==1){
		fprintf(fptr,"{\\f2\\fs16 ");
		fprintf(fptr,"information");
		fprintf(fptr,"\\tab  \\tab }");
	}
	strcpy(str,"\\b");
	gnull = gnullSMA(t,MA);
	for(laststate=0,i=j=0; j < lengthSMA(t,MA); i++){
		   if(gnull[i] == '_'){
			if(!show_position || show_position[i]){
			  state='x';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs14\\cf1  ",str);
                          } laststate = state;
			}
		   } else {
			j++;
			state = put_info_rma_typ(fptr, info_show[j],laststate,str);
			laststate = state;
		   }
	} fprintf(fptr,"}\n{\\f2\\fs16 ");
	fprintf(fptr,"\\tab }{\\f3 \n\\par }");
	//**************** output pvalue ... *******************
	if(t==1){
		fprintf(fptr,"{\\f2\\fs16 ");
		fprintf(fptr,"binomial tail");
		fprintf(fptr,"\\tab  \\tab }");
	}
	strcpy(str,"\\b");
	gnull = gnullSMA(t,MA);
	for(laststate=0,i=j=0; j < lengthSMA(t,MA); i++){
		   if(gnull[i] == '_'){
			if(!show_position || show_position[i]){
			  state='x';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs14\\cf1  ",str);
                          } laststate = state;
			}
		   } else {
			j++;
			state = put_info_rma_typ(fptr,pval_show[j],laststate,str);
			laststate = state;
		   }
	} fprintf(fptr,"}\n{\\f2\\fs16 ");
	fprintf(fptr,"\\tab }{\\f3 \n\\par }");
#if 1
	//**************** output marginal probability... *******************
	if(MargProb){
	  if(t==1){
		fprintf(fptr,"{\\f2\\fs16 ");
		// fprintf(fptr,"marginal prob");
		fprintf(fptr,"%% seqs. aligned");
		fprintf(fptr,"\\tab  \\tab }");
	  }
	  strcpy(str,"\\b");
	  gnull = gnullSMA(t,MA);
	  char	marg_prob;
	  Int4	prob;
	  for(laststate=0,i=j=0; j < lengthSMA(t,MA); i++){
		   if(gnull[i] == '_'){
			if(!show_position || show_position[i]){
			  state=' ';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs14\\cf1  ",str);
                          } laststate = state;
			}
		   } else {
			j++;
			state='X';
			prob = (Int4) floor(10.0*MargProb[j][j]);
			switch(prob){
			  case 0: marg_prob='.'; break;
			  case 1: marg_prob='1'; break;
			  case 2: marg_prob='2'; break;
			  case 3: marg_prob='3'; break;
			  case 4: marg_prob='4'; break;
			  case 5: marg_prob='5'; break;
			  case 6: marg_prob='6'; break;
			  case 7: marg_prob='7'; break;
			  case 8: marg_prob='8'; break;
			  case 9: marg_prob='9'; break;
			  case 10: marg_prob='*'; break;
			  default: marg_prob='?'; break;
			}
			if(state == laststate) fprintf(fptr,"%c",marg_prob);
			else { 
			   if(laststate != 0) fprintf(fptr,"}\n");
	        	   fprintf(fptr,"{%s\\f3\\cf1 %c",str,marg_prob);
			} laststate = state;
		   }
	  } fprintf(fptr,"}\n{\\f2\\fs16 ");
	  fprintf(fptr,"\\tab }{\\f3 \n\\par }");
	}
	//*******************************************************************
#endif
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3\n\\par\n\\par }\n");

	// if(fp_ras) put_rasmol_rma(fp_ras, MA, t,su,A,info_show,info_real);
	if(fp_ras) put_rasmol_rma(fp_ras, MA, t,su,A,pval_show,info_real);
	free(info_real);
	if(show_position) free(show_position);
	for(s=1; s <= length; s++) if(conserved[s] != NULL) free(conserved[s]);
	free(conserved); free(info_show); free(pval_show);
    } free(observed); free(use);
}

void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, double ***freq, sma_typ MA,
	cma_typ cma,FILE *fp_ras, char su,char Color)
{ NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,freq,MA,cma,0,fp_ras,su,Color); }

void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, cma_typ cma,double ***Observed,char ColorCode)
{
	FILE	*fp=tmpfile();
	double	***freq;
	Int4	b;
	
	PutAlnCMSA(fp,cma); rewind(fp);
        sma_typ MA=ReadSMA(fp); fclose(fp);

	NEWPP(freq,nBlksCMSA(cma)+3,double);
	for(b=1; b <= nBlksCMSA(cma); b++){
	   NEWP(freq[b],LengthCMSA(b,cma)+3,double);
	   //for(Int4 s=1; s <=LengthCMSA(b,cma); s++){ freq[b][s]=FreqSMA(MA);	}
	}
	NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,freq,MA,cma,Observed,ColorCode);
	for(b=1; b <= nBlksCMSA(cma); b++){ free(freq[b]); } free(freq);
	NilSMA(MA);
}

void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, cma_typ cma, char ColorCode)
{ NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,cma,0,ColorCode); }

double *HenikoffWeights2CMSA(cma_typ cma,double **wCnts)
{
        double  *w;
        Int4    i,n,t,st,p,l,res,pos[5],col,nCol,posit,sum;
        Int4    **rawCnts,**blPos,*blLen,*nResTyp;
        Int4    nSeq = NumSeqsCMSA(cma);
        Int4    nBlks = nBlksCMSA(cma);
        a_type  A = AlphabetCMSA(cma);
        ss_type SeqSet = DataCMSA(cma);
        NEW(w,nSeq+1,double); NEWP(blPos,nBlks+1,Int4);
        for(i=0;i<=nBlks;i++) NEW(blPos[i],nSeq+1,Int4);
        NEW(blLen,nBlks+1,Int4);
        for(sum=0,i=1;i<=nBlks;i++) {
                blLen[i] = LengthCMSA(i,cma);
                sum += blLen[i];
        }
        nCol = sum;
        NEWP(rawCnts,nCol+1,Int4);
        for(t=0;t<=nCol;t++){ NEW(rawCnts[t],nAlpha(A)+1,Int4);}
        for(n=1;n<=nSeq;n++){
           for(t=1;t<=nBlks;t++){
               st = PosSiteCMSA(t,n,pos,cma);
               if(st != 0) blPos[t][n] = pos[1]; else blPos[t][n] = 0;
           }
        }
        for(n=1;n<=nSeq;n++){
                col = 0;
                for(t=1;t<=nBlks;t++){
                        p = blPos[t][n];
                        if(p != 0){
                                for(l=1;l<=blLen[t];l++){
                                        col += 1;
                                        posit = p+l-1;
                                        res=SeqP(n,posit,SeqSet);
                                        rawCnts[col][res] += 1;
                                }
                        }
                }
        }
        NEW(nResTyp,nCol+1,Int4);
        for(t=1;t<=nCol;t++){
                for(i=1;i<=nAlpha(A);i++){
                        if (rawCnts[t][i] != 0) nResTyp[t] += 1;
                }
        }
        for(col=0,t=1;t<=nBlks;t++){
                for(l=1;l<=blLen[t];l++){
                        col += 1;
                        for(n=1;n<=nSeq;n++){
                                p = blPos[t][n];
                                if(p !=0){
                                        posit = p+l-1;
                                        res = SeqP(n,posit,SeqSet);
                                        w[n] += 1./(rawCnts[col][res]*nResTyp[col]);
                                }
                        }
                }
        }
        double  maxw;
        for(maxw=0.,n=1;n<=nSeq;n++){ if(w[n] > maxw) { maxw = w[n]; } }
        for(n=1;n<=nSeq;n++) { w[n] /= maxw; }
        for(n=1;n<=nSeq;n++){
                col = 0; 
                for(t=1;t<=nBlks;t++){
                        p = blPos[t][n];
                        if(p != 0){
                                for(l=1;l<=blLen[t];l++){
                                        col += 1;
                                        posit = p+l-1;
                                        res = SeqP(n,posit,SeqSet);
                                        wCnts[col][res] += w[n];
                                }
                        }
                }
        }
        for(t=0;t<=nCol;t++){ free(rawCnts[t]); } 
        free(rawCnts); free(nResTyp);
        for(t=0;t<=nBlks;t++){ free(blPos[t]); } 
        free(blPos); free(blLen);
        return w;
}

double 	WeightedCountsCMSA(cma_typ cma, Int4 *NoSeq, double *wNoSeq, double **wCnts)
//returns weighted number of sequences
//NoSeq = raw number of seqs at each position
//wNoSeq = weighted number of seqs at each position
//wCnts = weighted counts
{
        Int4            i,j,totNseq=NumSeqsCMSA(cma),*bound1,*bound2,lenSeq;
        gss_typ         *gss = gssCMSA(cma);
        double          nseq=0.;
        e_type          E;
        unsigned char   *seq;

        double *w = HenikoffWeights2CMSA(cma,wCnts);
        for(i=1;i<=totNseq;i++) nseq += w[i];
        NEW(bound1,totNseq+1,Int4); NEW(bound2,totNseq+1,Int4);
        for(i=1;i<=totNseq;i++){
                E = gss->FakeSeq(i);
                seq = SeqPtr(E);
                lenSeq = LenSeq(E);
                j=1; while(seq[j] == 0) j++;
                bound1[i] = j;
                j=lenSeq; while(seq[j] == 0) j--;
                bound2[i] = j;
        }
        for(j=1;j<=TotalLenCMSA(cma);j++){
                for(i=1;i<=totNseq;i++){
                        if(j>=bound1[i] && j<=bound2[i]){
                                NoSeq[j] +=1;
                                wNoSeq[j] += w[i];
                        }
                }
        } free(bound1); free(bound2); free(w);
        return nseq;
}

//******************* END NEW CMA2RTF ***************************

void    SubtractInfoCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, 
	double infoLO, double infoHI, ss_type key_seq, char PageSetUp, 
	Int4 rpts, char *pssm_arg,cma_typ main_cma, cma_typ sub_cma,
	char mode,Int4 JackCut)
{ SubtractInfoCMA2RTF(fptr, "1XXXX",'X',cbp_cut, purge_cut, infoLO, infoHI, key_seq, 
		PageSetUp, rpts, pssm_arg,main_cma, sub_cma, mode,JackCut,0); }

void    SubtractInfoCMA2RTF(FILE *fptr,char *filename, char subunit, double cbp_cut, 
	double purge_cut, double infoLO, double infoHI, ss_type key_seq, char PageSetUp, 
	Int4 rpts, char *pssm_arg,cma_typ main_cma, cma_typ sub_cma,
	char mode,Int4 JackCut,double **MargProb)
{
	// assert(nBlksCMSA(main_cma) == 1); // just one block for right now...
	assert(nBlksCMSA(sub_cma) == 1); // just one block for now...

	FILE	*fp_ras=0;
        FILE    *fp=tmpfile();
	double	***Observed=0;
	double	***freq;	// freq[blk][col][res];
	a_type	A=AlphabetCMSA(main_cma);

        PutAlnCMSA(fp,sub_cma); rewind(fp);
        sma_typ MA=ReadSMA(fp); fclose(fp);

#if 1
	if(MargProb==0){
	   Int4		i,j,*NoSeq,Len=LengthCMSA(1,sub_cma);
	   double	*wNoSeq, **wCnts;
	   NEW(NoSeq,Len+5,Int4); NEW(wNoSeq,Len+5,double);
	   NEWP(wCnts,Len+5,double); NEWP(MargProb,Len+3,double);
           for(i = 1; i <= Len; i++){
		NEW(wCnts[i],nAlpha(A)+5,double); 
		NEW(MargProb[i],Len+3,double);
	   }
	   double effective_no_seq = WeightedCountsCMSA(main_cma, NoSeq, wNoSeq, wCnts);
	   fprintf(stderr,"effective_no_seq = %.3f\n",effective_no_seq);
           for(i = 1,j=1; i <= Len; i++){
                double max_p=0.0;
                Int4    npos=0;
                MargProb[i][i] = wNoSeq[j]/effective_no_seq;
       fprintf(stderr,"effective_num_seq=%.2f;",wNoSeq[i]);
       fprintf(stderr,"actual number=%d\n",NoSeq[i]);
                j++;
           }
	}
#endif
	if(mode=='t'){ // new method...
	  // 1. Find alignment relating main family to subfamily via consensus.
	  freq=SuperModelFreqsCMSA(rpts,pssm_arg,sub_cma,main_cma,&Observed,JackCut);
	  // freq=SuperModelFreqsCMSA(rpts,pssm_arg,sub_cma,main_cma,&Observed,JackCut,MargProb);
#if 1	  // print rasmol header.
	  char su=subunit;
	  fp_ras= open_file(filename,".ras","w");
	  if(fp_ras){
	        // fprintf(fp_ras,"load pdb %s.pdb\n",filename);
	        fprintf(fp_ras,"load pdb %s.pdb\n","1XXXX");
        	fprintf(fp_ras,"set background black\n");
        	fprintf(fp_ras,"color [220,220,255]\nwireframe off\n");
        	if(su == 'X') {
                  fprintf(fp_ras,"select protein\nset strands 1\nstrands\n");
                  fprintf(fp_ras,"strands\n");
        	} else {
                  fprintf(fp_ras,"select *%c\nset strands 1\nstrands\n",su);
                  fprintf(fp_ras,"select not *%c and protein\n",su);
                  fprintf(fp_ras,"color [220,220,255]\n");
                  fprintf(fp_ras,"set strands 1\nstrands\n");
        	}
	  }
#endif
PutSMA(stdout, MA);
	  NewCMA2RTF_HEADER(fptr,PageSetUp);
	  NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,sub_cma,' ');
	  NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,freq,MA,sub_cma,
			0,fp_ras,su,MargProb,"Difference",'R');
	  NewCMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,sub_cma,Observed,'B');
	  NewCMA2RTF_TAIL(fptr);
#if 1	  // print rasmol tail.
	  if(fp_ras){
        	fprintf(fp_ras,"select dna or rna\nwireframe 80\n");
        	fprintf(fp_ras,"select all\n"); fclose(fp_ras);
	  }
#endif
	  for(Int4 s=1; s <=LengthCMSA(1,sub_cma); s++){
		if(Observed[1][s]) free(Observed[1][s]); 
	  } free(Observed[1]); free(Observed);
	} else {	// old method
	  freq=SuperModelFreqsCMSA(rpts,pssm_arg,sub_cma,main_cma,0,JackCut);
	  CMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,PageSetUp,freq,MA,sub_cma);
	}
	for(Int4 s=1; s <=LengthCMSA(1,sub_cma); s++){ 
		if(freq[1][s]) free(freq[1][s]);
	} free(freq[1]); free(freq);
        NilSMA(MA);
}


