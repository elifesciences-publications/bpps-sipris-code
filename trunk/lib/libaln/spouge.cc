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

/* spouge.c - john spouge statistics */
#include "spouge.h"

#if 0	// does not compile at TIGR
#include "sls_pssm.h"

using namespace Sls;

static pssm_typ *spouge_pssm=0;
#endif

void InitPsiGap (js_type S)
{

	S->initialize = FALSE;
#if 0
    Int4  *scoreGap[2];

	/*** TEST SPOUGE SUGGESTION ***/
	scoreGap[0] = S->scoreGap[0]; S->scoreGap[0] = NULL;
	scoreGap[1] = S->scoreGap[S->num_mtfs];
	S->scoreGap[S->num_mtfs]=NULL;
	/*** TEST SPOUGE SUGGESTION ***/
#endif
    switch(S->mode) {
#if 0	// shut off until I can fix this....
     case 'a' :
	fprintf(stderr,"compute alternative Gap function mode\n");
#if 0
	for(Int4 x=0; x <= S->num_mtfs; x++) printf("S->scoreGap[%d]=%d\n",
			x,(Int4)S->scoreGap[x]);
#endif
	spouge_pssm = new pssm_typ(0,S->maxLength,S->num_aa,S->freq,S->num_mtfs,
			S->lengthMotif0,S->scoreAmino,S->maxLength,S->scoreGap,S->maxLength,
			S->scoreGap,9); 
			// <upperGapLength_>,<penalty_function_>,
			// <alternative_upperGapLength>,<alternative_penalty_function>,9);
	fprintf(stderr,"Done computing alternative Gap function mode.\n");
#if 0
	S->neuwald = (void *) Vec_NewGlobalNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, NULL, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);

	spouge_pssm = new pssm_typ(		//new constructor
              const Int4 dimAminos_,//#(amino     acids)
              const double *freqAminos_, // amino freq's; [0...dimAminos_-1]
              const Int4 dimMotif_,//     #(motifs)
              const Int4 *lengthMotif_,// lengthMotif_ [0...dimMotif_-1]
	      Int4 ***scoreAmino_,// position-wise amino acid scores 
			// [motif   #(0...dimMotif_-1)]
			// [position(0..lengthMotif_-1)]
			// [aminoacid(0...dimAminos_-1)]; 
			// dimensions of scoreAmino_ must be adjusted with dimMotif_, 
			// lengthMotif_ and dimAminos_; else     constructor     produces error
              const Int4 upperGapLength_,// upper bound on all gaps; 
				// must be >=0; and     =-1     if undefined
              Int4 **penalty_function_,//penalty function [0...dimMotif_][0...upperGapLength_]
              const Int4 alternative_upperGapLength_,// upper bound on all alternative gaps;
				// must     be >=0; and     =-1     if undefined
              Int4 **alternative_penalty_function_,//alternative penalty function 
				// [0...dimMotif_][0...alternative_upperGapLength_]
              const Int4 max_penalty_number_of_digits_=9);  
				//maximum number of digits for internal calculations

#endif
     break;
#endif	// above is broken.
     case 'd' :
     case 'D' :
	fprintf(stderr,"compute GapLocalDisjointThreshold\n");
	S->neuwald = (void *) Vec_NewLocalDisjointNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, NULL, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);
	if(S->neuwald == NULL) 
		print_error("LocalDisjoint : failure");
     break;
#if 0
     case 'o' :
	print_error("GapLocalOverlapThreshold not yet implemented");
     break;
     case 'O' :
	fprintf(stderr,"compute LocalOverlapThreshold\n");
	S->neuwald = (void *) Vec_NewLocalOverlapNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, NULL, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);
#if 0
	if(!Vec_TestSentinel(S->num_mtfs,S->lengthMotif,S->sentinel)){
		print_error("LocalCoreOverlap : sentinel failure");
	}
	S->neuwald = (void *) Vec_NewLocalCoreDisjointNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, S->sentinel, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);
#endif
	if(S->neuwald == NULL) print_error("LocalOverlap : failure");
     break;
#endif
     case 'g':
     case 'G':
     case 'O':
     case 'o':
     case 'e':
	fprintf(stderr,"compute GapGlobalThreshold\n");
	S->neuwald = (void *) Vec_NewGlobalNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, NULL, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);
	if(S->neuwald == NULL) print_error("Global : failure");
     break;
     case 'c':
     case 'C':
	fprintf(stderr,"compute GapLocalCoreThreshold\n");
	if(!Vec_TestSentinel(S->num_mtfs,S->lengthMotif,S->sentinel)){
		print_error("LocalCore : sentinel failure");
	}
	S->neuwald = (void *) Vec_NewLocalCoreDisjointNeuwald(
		S->maxLength, S->num_aa, S->freq, S->num_mtfs, 
		S->lengthMotif, S->sentinel, S->scoreAmino, S->dimscoreGap,
		S->scoreGap);
	if(S->neuwald == NULL) print_error("LocalCore : failure");
     break;
     default:	
		print_error("InitPsiGap( ) Input error: failure");
     break;
    }
    fprintf(stderr,"done computing\n");
#if 0
	/*** TEST SPOUGE SUGGESTION ***/
	S->scoreGap[0] = scoreGap[0];
	S->scoreGap[S->num_mtfs] = scoreGap[1];
	/*** TEST SPOUGE SUGGESTION ***/
#endif
}

Int4	SpougeThreshold(js_type S,double target, Int4 length)
{ 
	if(S->initialize) InitPsiGap(S);
#if 0
    if(S->mode == 'o' || S->mode == 'O')
	return Vec_NeuwaldThreshold(S->neuwald,length+S->sum_extends,target); 
    else return Vec_NeuwaldThreshold(S->neuwald,length,target); 
#endif 
    if(S->mode == 'a'){ 
	print_error("spouge mode 'a' not working");
#if 0	// **************** shut off: broken ******************
	fprintf(stderr,"compute pssm threshold\n");
	bool calculation_is_successful;
	// double score= spouge_pssm->threshold_function(2,length,target,calculation_is_successful,0);
	// target = -log(1.0 - target);

	double score= spouge_pssm->threshold_function(0,length,target,calculation_is_successful);
	if(!calculation_is_successful) fprintf(stderr,"Computational error!\n");
	else fprintf(stderr,"Computed threshold score successful.\n");
	fprintf(stderr,"Score = %d; target = %g; length=%d\n",(Int4) score, target,length);
// score = -99;
	return (Int4) score;
#if 0
	// need to add this to new routines.
double pssm_typ::threshold_function(   
	  //returns score t for which P(T>=t)~threshold_; current object must be correct
  const Int4 maximum_number_of_alternative_penalties_,
	  //maximum number of alternative penalties; must be between 0 and d_dimMotif+1
  const Int4 sequence_length_,
					//sequence length; must be >=sum of motifs lengths
  const double threshold_,	//threshold for the tail probability
        bool &calculation_is_successful_, //return true iff claculation is successful
        double eps_=0);
		//accuracy of the result; 
		// if this parameter does not defined then celected automatically (recommended)

#endif
#endif	//********************************************

    } else {
	Int4 score=Vec_NeuwaldThreshold(S->neuwald,length,target); 
	fprintf(stderr,"Score = %d; target = %g\n",score, target);
	return score;
    }
}

/*** NeuwaldTail returns -ln(P-value); convert to log10(P-value); ***/
double	SpougePvalue(js_type S,UInt8 dbs_leng, Int4 length, Int4 score)
{ 
	double	p;

	if(S->initialize) InitPsiGap(S);
#if 0
    if(S->mode == 'o' || S->mode == 'O') {
	p = Vec_NeuwaldTail(S->neuwald,length+S->sum_extends,score)/2.3025851; 
    } else p = Vec_NeuwaldTail(S->neuwald, length, score)/2.3025851; 
#endif
    bool calculation_is_successful;
    double threshold=1000.0;
    if(S->mode == 'a'){
	print_error("spouge mode 'a' not working");
#if 0	//*************** shut off *****************
	fprintf(stderr,"compute pssm pvalue\n");
	p = spouge_pssm->calculate_tail_probability(0,length,score,
			calculation_is_successful);
			// calculation_is_successful,threshold);
	// p = 1 - exp(-p);
	if(calculation_is_successful) fprintf(stderr,"successful pvalue calculation.\n");
	else fprintf(stderr,"Unsuccessful pvalue calculation.\n");
#endif
    } else p = Vec_NeuwaldTail(S->neuwald, length, score)/2.3025851; 
	p *= (double)dbs_leng/(double)length;
	return p; 
}

void	PutHistSpouge(FILE *fp, Int4 t, double inc, js_type S)
{
	h_type  H;
	Int4	g;
	char	str[100];

	if(t >= 0 && t <= S->num_mtfs && (S->maxgap[t] > 0)) {
	 if(S->scoreGap[t] != NULL){
	   sprintf(str,"motifs %d-%d gap penalties", t,t+1);
	   H=Histogram(str,0,S->maxgap[t],inc);
	   for(g=0; g<= S->maxLength; g++){
		if(S->scoreGap[t][g] > SHRT_MIN) IncdMHist(g, -S->scoreGap[t][g], H);
		// if(S->scoreGap[t][g] > SHRT_MIN) IncdMHist(g, S->scoreGap[t][g], H);
	   }
	   PutHist(fp,60,H);
           NilHist(H);
	 }
	}
}

Int4	MaxGapScoreSpouge(js_type S)
/*** returns the maximum gap score that can be obtained ***/
{
	Int4	t,g,maxscore,max;

	if(S->scoreGap == NULL) return 0;
	for(maxscore=0,t=1; t< S->num_mtfs; t++){
	   max = INT4_MIN;
	   for(g=0; g<= S->maxLength; g++) {
		max = MAXIMUM(Int4,max, S->scoreGap[t][g]);
	   }
	   maxscore += max;
	}
	return maxscore;
}

void	MedianSpouge(double *median, js_type S)
/*** WARNING: assumes that scoreGap[t][g] = the number of gaps observed! ***/
{
	Int4	t,g,total,half;

	median[0]=median[S->num_mtfs] = 0;
	for(t=1; t< S->num_mtfs; t++){
	   for(total=g=0; g<= S->maxLength; g++) total+=S->scoreGap[t][g];
	   half = total/2;
	   for(total=g=0; g <= S->maxLength; g++){
		total+=S->scoreGap[t][g];
		if(total > half) { median[t] = (double) g; break; }
	   }
	}
}

void	SmoothGapScoresSpouge2(Int4 nl, Int4 nr, Int4 m, js_type S)
/*********************************************************************
  WARNING: ans must have dimensions [1...2*n] and n MUST be an integer
  power of two!
/*********************************************************************/
{
	double	*data,*ans,*respns,total,pseudo,v,factor=5.0,expect;
	Int4	*lscore,adj,mode,max_g,s,num_g,ninety;
	double	sum,ave;
	UInt4	np,n;
	Int4	cutoff,t,g,num;
	BooLean	flag;

	n = S->maxLength +1; np = 1;
	/*** set n = to the smallest integer power of two >= maxLength +1. ***/
	do { np *= 2; } while(np < n); n = np;
	np = nl+nr+1;
	NEW(data,n+3,double); NEW(ans,2*n+3,double); NEW(respns,n+3,double);
	NEW(lscore,n+3,Int4);
	S->dimscoreGap[0] = 0; S->dimscoreGap[S->num_mtfs] = 0;
	for(g=0; g<= S->maxLength; g++) 
		S->scoreGap[0][g] = S->scoreGap[S->num_mtfs][g]=0;
	// fprintf(stderr,"\n******************** input gaps ********************\n");
#if 1	//******************************************************
	for(t=1; t< S->num_mtfs; t++){
	    S->maxgap[t]=100;
	    PutHistSpouge(stdout, t, 1.0, S);
	    S->maxgap[t]=0;
	    h_type  H= Histogram("Gonnet gap function",0,S->maxLength,1);
	    Int4 extend=1;
	    double tot;
	    double *pen = new double [S->maxLength+2];
	    for(g=0; g<= S->maxLength; g++) pen[g]=0.0;
	    for(tot=0,0,g=0; g<= S->maxLength; g++) {
		s = S->scoreGap[t][g]; 
		if(s > 0){
	// s=1;
	tot+=s;
		  // for(Int4 gp=0; gp <=S->maxLength; gp++)
		  for(Int4 gp=g+1; gp <=S->maxLength; gp++)
		  {
		   if(gp != g){
			pen[gp]+=s*log(pow(abs(gp-g),1.7));  // Gonnet
			pen[gp]+=s*abs(gp-g);   // affine
		   }
		  }
		}
	    }  
	    double min_pen;
	    for(min_pen=DBL_MAX,g=0; g<= S->maxLength; g++){ 
		pen[g]/=tot; 
		if(min_pen > pen[g]) min_pen=pen[g];
	    }
#if 1
	    for(g=0; g<= S->maxLength; g++){ 
		IncdMHist(g,(Int4)((pen[g]-min_pen)*extend), H);
	    }
	    // PutHist(stdout,60,H); NilHist(H);
#endif
	    delete [ ] pen;
	// } exit(1);
#endif	//******************************************************
	// for(t=1; t< S->num_mtfs; t++){
	   data++; ninety = -1;
	   for(total=0.0,max_g=g=num_g=0; g<= S->maxLength; g++){
		s = S->scoreGap[t][g];
		if(s > 0){ max_g = g; }
		num_g += s;
		data[g] = (double) s;
		total += data[g];
	   }
	   expect = 2.0*total/(double) n;
	   // fprintf(stderr,"np =%d; maxLength = %d\n",np,n);
	   data--;
	   savgol(respns, np, nl, nr, 0, 2);
	   convlv(data, n, respns, np, 1, ans);
	   ans++;
	   for(flag=FALSE,s=0,cutoff=-1,g=0; g<= S->maxLength; g++){
		ans[g] = ans[g]/expect;
		s += S->scoreGap[t][g];
// fprintf(stderr,"ans[%d] = %g; s = %d; num_g = %d\n",g,ans[g],s,num_g); 
		if(ans[g] < 0.0) { 
		    if(flag) { cutoff = g; break; }
		} else if(ans[g] > 0.0){
		    if(((double)s/(double)num_g) > 0.9) {	/** > 90% **/
			 flag = TRUE;
		    }
		}
	   }
	   // cutoff = MINIMUM(Int4, cutoff*2, S->maxLength); 
	   cutoff = MINIMUM(Int4, cutoff+10, S->maxLength); 
	   cutoff = MAXIMUM(Int4, cutoff, max_g); 
	   // fprintf(stderr,"Motif %d: gap of 0-%d; max_g = %d\n",t,cutoff,max_g);
	   S->maxgap[t] = max_g;
	 if(cutoff > 0) {
/****** NEW: redo with more stringent settings ****/
	   ans--;
	   savgol(respns, np, nl, nr, 0, m);
	   convlv(data, n, respns, np, 1, ans);
	   ans++;
/****** NEW: redo with more stringent settings ****/
	   for(sum=0.0, num=g=0; g<= cutoff; g++){
		ans[g] = ans[g]/expect;
		if(ans[g] < 1.0) ans[g] = 0.0;
		else { ans[g] = factor*log(ans[g]); num++; }
		sum += ans[g];
		lscore[g] = (Int4) floor(ans[g] + 0.4999); 
	   } ave = sum/(double) num;
	   adj = (Int4) floor(ave + 0.4999);
	   // fprintf(stderr,"sum = %g; num = %d; ave = %g; adj = %d\n",sum,num,ave,adj);
	   for(g=0; g<= S->maxLength; g++){
		if(g > cutoff) S->scoreGap[t][g] = SHRT_MIN;
		else S->scoreGap[t][g] = lscore[g] - adj;
	   }
#if 1	   // see if this fixes -uo option problem...
	   S->scoreGap[t][S->maxLength+1] = SHRT_MIN;
#endif
	   S->dimscoreGap[t] = cutoff+1;
	   PutHistSpouge(stderr, t, 1.0, S); 
#if 1
	    // for(g=0; g<= S->maxLength; g++){ 
		// if(g < 100) IncdMHist(g,(Int4)-S->scoreGap[t][g], H);
	    // } 
	    PutHist(stdout,60,H); NilHist(H);
#endif
	 } else { /*** ignore the gaps in these regions ***/
		for(g=0; g<= S->maxLength; g++) S->scoreGap[t][g] = 0;
		S->dimscoreGap[t] = 0; 
	 }
	 ans--;
	}
	free(data); free(ans); free(respns); free(lscore);
	/****** set end gaps == 0 ***/
	for(g=0; g<= S->maxLength; g++){
		S->scoreGap[0][g] = S->scoreGap[S->num_mtfs][g] = 0;
	}
#if 1	// NEW: set maximum scores to zero...
	Int4	max;
	for(t=1; t< S->num_mtfs; t++){
	    for(max=g=0; g<= S->maxLength; g++){
		if(max < S->scoreGap[t][g]) max = S->scoreGap[t][g];
	    }
	    for(g=0; g<= S->maxLength; g++) S->scoreGap[t][g] -= max;
	}
#endif
}


void	SmoothGapScoresSpouge(Int4 nl, Int4 nr, Int4 m, js_type S)
/*********************************************************************
  WARNING: ans must have dimensions [1...2*n] and n MUST be an integer
  power of two!
/*********************************************************************/
{
	double	*data,*ans,*respns,total,pseudo,v,factor=5.0,expect;
	Int4	*lscore,adj,mode,max_g,s,num_g,ninety;
	double	sum,ave;
	UInt4	np,n;
	Int4	cutoff,t,g,num;
	BooLean	flag;

	// SmoothGapScoresSpouge2(nl,nr,m,S); return;
	n = S->maxLength +1; np = 1;
	/*** set n = to the smallest integer power of two >= maxLength +1. ***/
	do { np *= 2; } while(np < n); n = np;
	np = nl+nr+1;
	NEW(data,n+3,double); NEW(ans,2*n+3,double); NEW(respns,n+3,double);
	NEW(lscore,n+3,Int4);
	S->dimscoreGap[0] = 0; S->dimscoreGap[S->num_mtfs] = 0;
	for(g=0; g<= S->maxLength; g++) 
		S->scoreGap[0][g] = S->scoreGap[S->num_mtfs][g]=0;
	// fprintf(stderr,"\n******************** input gaps ********************\n");
	for(t=1; t< S->num_mtfs; t++){
	   data++; ninety = -1;
	   for(total=0.0,max_g=g=num_g=0; g<= S->maxLength; g++){
		s = S->scoreGap[t][g];
		if(s > 0){ max_g = g; }
		num_g += s;
		data[g] = (double) s;
		total += data[g];
	   }
	   expect = 2.0*total/(double) n;
	   // fprintf(stderr,"np =%d; maxLength = %d\n",np,n);
	   data--;
	   savgol(respns, np, nl, nr, 0, 2);
	   convlv(data, n, respns, np, 1, ans);
	   ans++;
	   for(flag=FALSE,s=0,cutoff=-1,g=0; g<= S->maxLength; g++){
		ans[g] = ans[g]/expect;
		s += S->scoreGap[t][g];
// fprintf(stderr,"ans[%d] = %g; s = %d; num_g = %d\n",g,ans[g],s,num_g); 
		if(ans[g] < 0.0) { 
		    if(flag) { cutoff = g; break; }
		} else if(ans[g] > 0.0){
		    if(((double)s/(double)num_g) > 0.9) {	/** > 90% **/
			 flag = TRUE;
		    }
		}
	   }
	   // cutoff = MINIMUM(Int4, cutoff*2, S->maxLength); 
	   cutoff = MINIMUM(Int4, cutoff+10, S->maxLength); 
	   cutoff = MAXIMUM(Int4, cutoff, max_g); 
	   // fprintf(stderr,"Motif %d: gap of 0-%d; max_g = %d\n",t,cutoff,max_g);
	   S->maxgap[t] = max_g;
	 if(cutoff > 0) {
/****** NEW: redo with more stringent settings ****/
	   ans--;
	   savgol(respns, np, nl, nr, 0, m);
	   convlv(data, n, respns, np, 1, ans);
	   ans++;
/****** NEW: redo with more stringent settings ****/
	   for(sum=0.0, num=g=0; g<= cutoff; g++){
		ans[g] = ans[g]/expect;
		if(ans[g] < 1.0) ans[g] = 0.0;
		else { ans[g] = factor*log(ans[g]); num++; }
		sum += ans[g];
		lscore[g] = (Int4) floor(ans[g] + 0.4999); 
	   } ave = sum/(double) num;
	   adj = (Int4) floor(ave + 0.4999);
	   // fprintf(stderr,"sum = %g; num = %d; ave = %g; adj = %d\n",sum,num,ave,adj);
	   for(g=0; g<= S->maxLength; g++){
		if(g > cutoff) S->scoreGap[t][g] = SHRT_MIN;
		else S->scoreGap[t][g] = lscore[g] - adj;
	   }
#if 1	   // see if this fixes -uo option problem...
	   S->scoreGap[t][S->maxLength+1] = SHRT_MIN;
#endif
	   S->dimscoreGap[t] = cutoff+1;
	   // PutHistSpouge(stderr, t, 1.0, S); 
	 } else { /*** ignore the gaps in these regions ***/
		for(g=0; g<= S->maxLength; g++) S->scoreGap[t][g] = 0;
		S->dimscoreGap[t] = 0; 
	 }
	 ans--;
	}
	free(data); free(ans); free(respns); free(lscore);
	/****** set end gaps == 0 ***/
	for(g=0; g<= S->maxLength; g++){
		S->scoreGap[0][g] = S->scoreGap[S->num_mtfs][g] = 0;
	}
#if 1	// NEW: set maximum scores to zero...
	Int4	max;
	for(t=1; t< S->num_mtfs; t++){
	    for(max=g=0; g<= S->maxLength; g++){
		if(max < S->scoreGap[t][g]) max = S->scoreGap[t][g];
	    }
	    for(g=0; g<= S->maxLength; g++) S->scoreGap[t][g] -= max;
	    // PutHistSpouge(stderr, t, 1.0, S); 
	}
#endif
#if 1	// NEW: set tail to log-linear...
	for(t=1; t< S->num_mtfs; t++){
	   // PutHistSpouge(stderr, t, 1.0, S); 
	   for(g=S->maxLength; g >= 0; g--) if(S->scoreGap[t][g] > SHRT_MIN) break;
	   for(max=S->scoreGap[t][g]; g >= 0; g--){
		if(S->scoreGap[t][g] > max) break;
	   } g++;
// std::cerr << g; std::cerr << " ********************************\n"; std::cerr << max; 
//	   std::cerr << " ********************************\n";
	   Int4 extend=1,gp;
	   for(gp=1; g <= S->maxLength; g++,gp++){
		S->scoreGap[t][g]-=gp*factor;   // affine
	   } // PutHistSpouge(stderr, t, 1.0, S); 
//	   std::cerr << " ********************************\n";
	}
#endif
}

js_type	MkSpouge(Int4 num_aa, double *freq, Int4 num_mtfs, smx_typ *sM,
        Int4 maxLength, Int4 **scoreGap, char mode)
{
	js_type S;
	Int4	i,m,j,r,len;

	NEW(S,1,spouge_type);
	S->num_aa = num_aa;
	S->maxLength = maxLength;
	S->mode = mode;
	S->initialize = TRUE;
	S->neuwald = NULL;
	S->sum_extends=0;
	NEW(S->freq,num_aa+2,double);
	for(i=0; i< num_aa; i++){ S->freq[i] = freq[i]; }
	S->num_mtfs = num_mtfs;
	NEW(S->lengthMotif, num_mtfs+1, size_t);
	NEW(S->lengthMotif0, num_mtfs+1, Int4);
	NEW(S->maxgap, num_mtfs+2, Int4);
	for(i=0; i< num_mtfs; i++){
		S->lengthMotif[i] = LenSMatrix(sM[i+1]);
		S->lengthMotif0[i] = LenSMatrix(sM[i+1]);
	}

/********** sentinel [0...dimMotif-1][2] ******************
0 <= sentinel [i][0] <= sentinel [i][1] + 1 <= lengthMotif [i]
 **********************************************************/
	// if(mode == 'c' || mode == 'C' || mode == 'o' || mode == 'O'){
	if(mode == 'c' || mode == 'C'){
	    NEWP(S->sentinel, num_mtfs+2, size_t);
	    for(i=0; i< num_mtfs; i++){	
		NEW(S->sentinel[i],3,size_t);
		len = S->lengthMotif[i];
		if(len > 8){
		    S->sentinel[i][0] = 2;
		    S->sentinel[i][1] = len - 3;
		    S->sum_extends += 4;
		} else if(len > 5){
		    S->sentinel[i][0] = 1;
		    S->sentinel[i][1] = len - 2;
		    S->sum_extends += 2;
		} else {
		    S->sentinel[i][0] = 0;
		    S->sentinel[i][1] = len-1;
		}
#if 0	/** TEST: should be equivalent to global method **/
		    S->sentinel[i][0] = 0;
		    S->sentinel[i][1] = len-1;
		    S->sum_extends = 0;
#endif
	
	    }
	} else S->sentinel=NULL;

/********** sentinel [0...dimMotif-1][2] ******************/
/********** GAPS ******************/
  if(scoreGap != NULL && islower(mode)){
	NEW(S->dimscoreGap, num_mtfs+1, size_t);
	NEWP(S->scoreGap, num_mtfs+1, Int4);
	for(i=0; i<= num_mtfs; i++){
	   NEW(S->scoreGap[i], maxLength+2, Int4);
	   for(j=0; j<= maxLength; j++){
	   	S->scoreGap[i][j] = scoreGap[i][j]; 
// if(scoreGap[i][j]!=0)
//	fprintf(stderr,"scoreGap[%d][%d] = %d\n",i,j,scoreGap[i][j]);
	   }
	}
	SmoothGapScoresSpouge(12,12,4,S);
  } else S->scoreGap = NULL;
/********** GAPS ******************/

	NEW(S->threshold, maxLength +1, Int4);
	NEWPP(S->scoreAmino, num_mtfs +2, Int4);
	for(m=0; m< num_mtfs; m++){
	   NEWP(S->scoreAmino[m], S->lengthMotif[m]+2, Int4);
	   for(j=0; j< S->lengthMotif[m]; j++){
		NEW(S->scoreAmino[m][j], num_aa + 2, Int4);
		for(r=0; r< S->num_aa; r++){
			S->scoreAmino[m][j][r] = ValSMatrix(j+1,r,sM[m+1]);
		}
	   }
	}
#if 0	/**************** SPOUGE III ************************/
	fprintf(stderr,"max = %d; %d motifs; %d amino acids.\n",
		maxLength,num_mtfs,num_aa);
	// InitPsiGap (S); // do this only when necessary...
#endif /**************** SPOUGE III ************************/

	return S;
}

void	NilSpouge(js_type S)
{
	Int4	m,j;

	for(m=0; m< S->num_mtfs; m++){
	   for(j=0; j< S->lengthMotif[m]; j++){ free(S->scoreAmino[m][j]); }
	   free(S->scoreAmino[m]);
	}
	if(S->scoreGap!=NULL){
	  free(S->dimscoreGap);
	  for(j=0; j<= S->num_mtfs; j++){
	 	if(S->scoreGap[j] != NULL) free(S->scoreGap[j]);
	  }
	  free(S->scoreGap);
	}
	if(S->sentinel != NULL){
	    for(m=0; m< S->num_mtfs; m++) free(S->sentinel[m]);
	    free(S->sentinel);
	}
	if(S->neuwald != NULL) Vec_FreeNeuwald(S->neuwald);
	free(S->scoreAmino);
	free(S->maxgap);
	free(S->freq);
	free(S->threshold);
	free(S->lengthMotif);
	free(S->lengthMotif0);
	free(S);
#if 0
	if(spouge_pssm) delete spouge_pssm;
#endif
}

