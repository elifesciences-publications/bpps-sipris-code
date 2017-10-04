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

rtf_typ::rtf_typ(double FractHighlight,Int4 FontSize,Int4 Length, double cutoff,
	double **Obs,double **freq, 
	Int4 MinBars, e_type qE, a_type A, Int4 position, double min_res_freq,char SetsMode)
{ Construct(FractHighlight,FontSize,Length,cutoff,Obs,freq,0,0,MinBars,qE,A,position,
		min_res_freq,SetsMode,'B',-1.0,0,0,0); }

rtf_typ::rtf_typ(double FractHighlight,Int4 FontSize,Int4 Length, double cutoff,
	double **Obs,double **freq,
	Int4 MinBars,e_type qE,a_type A, Int4 position, double min_res_freq,char SetsMode,
	double **ObsBG,double **freqFG,char ModeLPR,double Alpha,Int4 A0, Int4 B0)
{ Construct(FractHighlight,FontSize,Length,cutoff,Obs,freq,ObsBG,freqFG,MinBars,qE,A,position,
		min_res_freq,SetsMode,ModeLPR,Alpha,A0,B0,0); }

// for cmc_typ == ModeLPR is 'U'.
rtf_typ::rtf_typ(double FractHighlight,Int4 FontSize,Int4 Length, double cutoff,
	double **Obs,double **freq,
	Int4 MinBars,e_type qE,a_type A, Int4 position, double min_res_freq,char SetsMode,
	double **ObsBG,double **freqFG,double **cmc_res_evals)
{ Construct(FractHighlight,FontSize,Length,cutoff,Obs,freq,ObsBG,freqFG,MinBars,qE,A,position,
		min_res_freq,SetsMode,'U',-1.0,0,0,cmc_res_evals); }


#define OUTPUT_BALL_IN_URN_PVALUES  0

void	rtf_typ::Construct(double FractHighlight,Int4 FontSize,Int4 Length,double cutoff,
	double **Obs,double **freq,
	double **ObsBG,double **freqFG,
	Int4 MinBars,e_type qE,a_type A,Int4 position, double min_res_freq,char SetsMode,
	char ModeLPR,double alpha, Int4 A0, Int4 B0,double **cmc_res_evals)
{ 
	length=Length;
	AB=A;
	ExpPttrns=1.0;
        MinResFreqCutoff=min_res_freq;
        // MinBlosum62Cutoff=-0.200; 
        MinBlosum62Cutoff=0.0; 
        AveBlosumScoreCutoff=0.10; 
	hist_height=50;
	IgnorePosition=FALSE;
	// hist_height=35;
	// hist_height=25;
	Position=position;
	Query=qE;
	minbars=MinBars;
	double  pval;
	static Int4 NumCalls=0;
	Int4	s,i;
	// below: New add residue set data type (use this for cbp_dfs routine)
	if(SetsMode == ' ') rst=0; else rst=new rst_typ(SetsMode);
	assert(minbars < hist_height);
        NEWP(res_evals,length+3, double); NEW(binomial, length+2, double);
	NEW(RelEntropy,length+2, double);
	dh_type dH=dheap(length+2,4);
        for(i=0; i <= length; i++){ NEW(res_evals[i],nAlpha(AB)+3,double); }
	for(s=1; s <= length; s++){
	   if(cmc_res_evals){
	      pval = -cmc_res_evals[s][0]; 
	      for(Int4 r=1; r <= nAlpha(AB); r++){ res_evals[s][r]=-cmc_res_evals[s][r]; }
	   } else if(ObsBG && freqFG){ // NEW: Cumulative binomial models...
	      // ModeLPR == 'P' ???
	      csp_typ *csp=new csp_typ(Obs[s],ObsBG[s],freqFG[s],freq[s],AB,min_res_freq,rst);
	      if(alpha > 0.0) csp->SetModeBPPS(A0,B0,alpha);
	      // if(csp->ComputeBinomialBPPS(-cutoff,&pval,res_evals[s],psst[s],subLPR[s]));
		// ^^ planned routine once I read BPPS subLPRs from file above.
	      if(csp->ComputeBinomialBPPS(-cutoff,&pval,res_evals[s])){
#if OUTPUT_BALL_IN_URN_PVALUES
	      for(Int4 r=1; r <= nAlpha(AB); r++){
		if(pval < -10 && res_evals[s][r] <= pval){
	           fprintf(stdout,"%c%d: Obs=%g;,ObsBG =%g; freqFG=%g; freq=%g; p=%g\n",
			AlphaChar(r,AB),s,Obs[s][r],ObsBG[s][r],freqFG[s][r],freq[s][r],
			res_evals[s][r]); // rst->Put(stdout);
		}
	      }
#endif
              if(TRUE || NumCalls == 0) insrtHeap(s,((keytyp) pval),dH);
	      } delete csp;
	   } else {	//  ModeLPR == 'B' ???
	      csp_typ *csp=new csp_typ(Obs[s],0,0,freq[s],AB,min_res_freq,rst);
	      csp->ComputeBinomial(-cutoff,&pval,res_evals[s]);
#if OUTPUT_BALL_IN_URN_PVALUES
	      double total=0.0,key=0.0;
	      for(Int4 r=1; r <= nAlpha(AB); r++){
		if(pval < -10 && res_evals[s][r] <= pval){
	           fprintf(stdout,"%c%d: Obs=%g;,freqBG=%g; p=%g (%g)\n",
			AlphaChar(r,AB),s,Obs[s][r],freq[s][r],res_evals[s][r],pval);
		   // rst->Put(stdout);
		   key+=Obs[s][r];
		} else total+=Obs[s][r];
	      } if(key > 0.0) fprintf(stdout," %g vs %g (%.2f\%)\n",key,total,100.0*key/(key+total));
#endif
	      delete csp;
	   }
	   RelEntropy[s] = ComputeRelEntropy(Obs[s],freq[s]);
	   // eventually replace this with use of Position within routines.
	   if(Query){
	      if(s != position){
		Int4 r = ResSeq(s,Query);
		binomial[s]=-res_evals[s][r];
	      } else binomial[s]=0.0;
	   } else {
	      if(s != position) binomial[s]=-pval;
	      else binomial[s]=0.0;
	   }
	}  //****************  end of for(s=1; ...) loop ********************88
#if 0	// print results as for BPPS procedure (*.pttrn file)...
	// if(ObsBG && freqFG && NumCalls==0)
	FILE *ofp=0; ofp=stdout;
	if(ObsBG && freqFG)
        {
	   NumCalls++;	// called for rtf, rtfQ, etc.
	   Int4	hits=0;
	   char	str[30];
	   pval = -minkeyHeap(dH);
	   while((s=delminHeap(dH)) != 0){
	     if(pval > 5.0){
		if(ofp && hits == 0){
			fprintf(ofp,"Method = '%c'\n", ModeLPR);
			fprintf(ofp,"print results as for BPPS procedure (*.pttrn file)...\n");
	        	fprintf(ofp,"%10s %6s %6s (%2s\%) %6s %6s (%s\%) %6s %6s\n",
				"Pattern","n1","n2"," ","m1","m2"," ","pval","WtNumSq");
		}
		hits++;
		double n1,n2,m1,m2;
		n1=n2=m1=m2=0.0;
	        Int4 nc=0;
		for(Int4 r=1; r <= nAlpha(AB); r++){
		   if(res_evals[s][r] <= -pval) {
			// fprintf(stderr,"%c%d eval=%lf\n",AlphaChar(r,AB),s,res_evals[s][r]);
			str[nc]=AlphaChar(r,AB); nc++; 
			n1+=Obs[s][r]; m1+=ObsBG[s][r]; 
		   } else { n2 += Obs[s][r]; m2+=ObsBG[s][r]; }
	        } 
		str[nc]=0; 
	      if(hits <= 25){ // DEBUG...
if(ofp && nc > 0){	// problem here but don't bother fixing this now...
		if(nc == 0){
		 for(Int4 r=1; r <= nAlpha(AB); r++){
		  fprintf(stderr,"res_evals[%d][%c] = %lf; pval=%f\n", s,AlphaChar(r,AB),res_evals[s][r],pval);
		 }
		}
	        fprintf(ofp,"%5s%-5d %6.0f %6.0f (%2.0f\%) %6.0f %6.0f (%2.0f\%) %6.1f %6.1f\n",
			str,s,n1,n2,100*n1/(n1+n2),m1,m2,100*m1/(m1+m2),pval,n1+n2+m1+m2);
		if(hits % 10 == 0) fprintf(ofp,"\n");
}
	      }
	     }
	     pval = -minkeyHeap(dH);
	   }
	}
#endif
	Nildheap(dH); fontsize=FontSize; 
	init(FractHighlight,binomial,cutoff); 
}

Int4    rtf_typ::get_unit_height(Int4 fontsize)
{
        switch(fontsize){       // NOTE: 13 (26) and 9 (18) are problem fonts...
           case 8: return 2;   // May be able to fix by changing paragraph crunch.
	   case 9: return 2; // NEW...Halfpoints...
           case 10: return 3;
	   case 11: return 3; // NEW...
           case 12: return 4;
           case 13: case 15: case 23: return 5; // NEW...
           case 14: case 16: case 22: return 5;
           case 17: case 19: case 21: return 6; // NEW...
           case 18: case 20: case 24: return 6;
           case 25: case 27: return 7; // NEW...
           case 26: case 28: return 8;
           case 29: case 31: case 33: case 35: case 43: case 45: return 9; // NEW...
           case 30: case 32: case 34: case 36: case 44: case 46: return 10;
           case 37: case 39: case 41: case 47: return 11; // NEW...
           case 38: case 40: case 42: case 48: return 12;
        } print_error("get_unit_height( ) input error");
        return 0;
}

void    rtf_typ::Free( )
{
	Int4 s;
	for(s=0; s<=length; s++) free(res_evals[s]); free(res_evals);
	free(binomial); free(PValue); free(Value); free(RelEntropy); // free(IgnorePos);
	for(s=1; s<=length; s++) free(ResValue[s]); free(ResValue);
	delete rst;
}

double	rtf_typ::GetLineToLog(Int4 Length, double *PvalReal, double MinPvalue,
	double Fraction)
#if 0
 Convert from the Fraction of residues to be highlighted to LineToLog.
 Given an array of Length p-values, find the value of LineToLog that will 
 highlight the inputed Fraction cutoff
 WARNING: Pvalues are expressed as -log10(p-value)!
 WARNING: Purify sometimes complaining about an 'Uninitialized memory read' in this routine.
 AFN: This apparently was the L variable returned below; This has been fixed.
#endif
{
	double	MaxPvalue,OneBelow;
	Int4	i,s,TargetNum;

	if(Fraction >= 1.0) return 1.0;		// use straight log scale 
	else if(Fraction <= 0.0) return 0.0;	// use straight linear scale 
	TargetNum = (Int4) floor((Length*Fraction) + 0.5);

        dh_type dH=dheap(Length+2,4);

	for(s=1; s<=Length; s++) insrtHeap(s,((keytyp)-PvalReal[s]),dH);
	MaxPvalue = -minkeyHeap(dH);
	TargetPval=OneBelow=MinPvalue;
	// fprintf(stderr,"Fraction = %f; TargetNum = %d\n",Fraction,TargetNum);
	for(i=0; (s=delminHeap(dH)) != 0; ){
	   if(PvalReal[s] < MinPvalue || i >= TargetNum){
		double tmp_d;
		while(!emptyHeap(dH) && (tmp_d=-minkeyHeap(dH)) == TargetPval){
			s=delminHeap(dH); i++;
		}
		OneBelow=PvalReal[s]; break;
	   } i++; TargetPval = PvalReal[s]; 
	   // fprintf(stderr,"%d: Pvalue[%d]=%g\n",i,s,PvalReal[s]);
	} Nildheap(dH);
	Int4 TargetShow=i;
	if(i <= 0) return 0.0;	// It doesn't matter.
	else if(TargetPval <= MinPvalue) return 1.0; // highlight everything possible
/********************* Finding the Target LineToLog L: ************
	Note that Value = pow(Trials,K)/K where K = 1.0 - L.
	 (K varies from 1 to 0 as L varies from 0 to 1.)
	Let's standardize the number of trials by dividing them all by 
	  MinTrials, so that MinTrials = 1, T = TargetTrials/MinTrials,
	  and  M = MaxTrials/MinTrials.
	This implies that 1/H = (pow(T,K)-1)/(pow(M,K)-1)
	 where H = hist_height, and 1 = pow(MinTrials==1,K).
	But 1/H = (pow(T,K)-1)r)/(pow(M,K)-1).
	   --> pow(M,K) -1 = H * (pow(T,K)-1)
	   --> pow(M,K) -1 = H * pow(T,K) - H
	Let Vm = pow(M,K)/K = LogToLinear value for MaxTrials 
	And Vt = pow(T,K)/K = LogToLinear value for TargetTrials.
	This routine chooses K to be the smallest value (in 1/1000
	   increments) such that Vt/Vm > 1/H or Vm/Vt < H.
	NOTE: -log(Pvalue) = log(Trials)
	NOTE: R = Ratio MaxTrials to TargetTrials
		= pow(pow(10,MaxPvalue)/pow(10,MinPvalue),K)/
	      		pow(pow(10,TargetPval)/pow(10,MinPvalue),K);
		= pow(10,K*(MaxPvalue-TargetPval));
/******************************************************************/
	double	L,H,K,R;	// K == 1 implies linear.
	H = (double) hist_height;
	// fprintf(stderr,"TargetNum = %d; TargetShow = %d\n",TargetNum,TargetShow);
	MaxRatio=1;
	for(L=1.0,K=0.001; K < 1.0; K += 0.001){
		R = pow(10,K*(MaxPvalue-TargetPval));
		if(R > H) break;
// fprintf(stderr,"K=%g; R = %g\n",K,R);
		MaxRatio=R;
		L = 1.0 - K;
	} return L;
}

double	rtf_typ::GetValue(double pvalue)
{
	if(pvalue < TargetPval) return -1.0;
	double V;
	switch(mode){
           case '0': V=(double) hist_height*
			pow(10,(pvalue-TargetPval))/pow(10,(MaxPval-TargetPval)); 
		break;
           case 'i': V=(double) hist_height*
			pow(10,(1.0-LinearToLog)*(pvalue-TargetPval))/MaxRatio; 
		break;
           case '1': // use logarithmic scale as it is.
		V = (double) hist_height*((pvalue-TargetPval)/(MaxPval-TargetPval));
		break;
	   default: print_error("This should not happen");
	} 
	if(V > (double) hist_height){
// fprintf(stderr,"\n{V = %.1f > hist_height (%d) '%c' Target=%.1f; MaxP=%.1f; pval=%.1f}\n\t",
//			V,hist_height,mode,TargetPval,MaxPval,pvalue);
		V = (double) hist_height;
	}
	return V;
}

void    rtf_typ::init(double FractHighlight,double *PvalReal,double MinPvalue)
// LinearToLog approach; seeks to give visual effect
// of true significance (i.e., actual number of trials needed
// to get observed elevation by chance alone, but not quite).
{
    Int4	r,i;

    assert(FractHighlight > 0.0 && FractHighlight <= 1.0);
    // if(FractHighlight ==0.0) FractHighlight = 0.1;
    LinearToLog=GetLineToLog(length,PvalReal,MinPvalue,FractHighlight);
    if(LinearToLog == 0.0) mode='0';		// linear
    else if(LinearToLog < 1.0) mode='i';	// intermediate
    else mode='1';				// log
    hist_up=get_unit_height(fontsize);
    NEW(PValue,length+2,double); 
    NEW(Value,length+2,double);
    MaxValue=MaxPval=0.0;
    NEWP(ResValue,length+2,double);
    for(i=1; i<=length; i++){
	PValue[i]=PvalReal[i];
	if(PValue[i] > MaxPval) MaxPval=PValue[i];	// MaxPval needed by GetValue()
    }
    for(i=1; i<=length; i++){
	NEW(ResValue[i],nAlpha(AB)+3,double);
	Value[i]=GetValue(PValue[i]);
#if 0	// for debugging related to histogram heights: afn 11/13/07.
        if(Value[i] > 1.0){
	  fprintf(stderr,"%-5d:  Val=%.1f; Pva =%.1f; Pttn=\"",i,Value[i],PValue[i]);
	  for(r=1; r <= nAlpha(AB); r++){
		ResValue[i][r]=GetValue(-res_evals[i][r]);
		// if(ResValue[i][r] >= PValue[i])
		if(ResValue[i][r] >= 1.0) {
		   fprintf(stderr,"%c(%.1f)",AlphaChar(r,AB),ResValue[i][r]);
		}
	  }
	  fprintf(stderr,"\"\n");
        } else for(r=1; r <= nAlpha(AB); r++) ResValue[i][r]=GetValue(-res_evals[i][r]);
#else
        for(r=1; r <= nAlpha(AB); r++) ResValue[i][r]=GetValue(-res_evals[i][r]);
#endif
	if(Value[i] > MaxValue) MaxValue=Value[i];
   } MinPval=MinPvalue;
   // fprintf(stderr,"\n(target=%.1f; MaxRatio=%.1f; mode=%c)\n\n",
   // 		TargetPval,MaxRatio,mode);	// DEBUG: rtf histogram heights

#if 0
     // h_type HG=Histogram("E-Values",0,10,0.05);
     h_type HG=Histogram("-log10(E-Values)",0,100,2.0);
     for(i=1; i <= length; i++){
                // if(PValue[i] > 0.0) IncdHist((double)length*20.0*pow(10,-PValue[i]),HG);
                if(PValue[i] > 0.0) IncdHist(PValue[i]-log10(length*20.0),HG);
     } 
     // PutHist(stdout,60,HG); 
     NilHist(HG);
     // fprintf(stdout,"MinPval = %.3f\n",MinPval);
#endif
    bars_s_cut=5;
}

char	*rtf_typ::Conserved(Int4 j)
{
	assert(j > 0 && j <= length);
	// return Conserved(IgnorePos[j],res_evals[j]); 	// AFN 12/22/05
	return Conserved(res_evals[j]); 
}

char	*rtf_typ::Conserved(double *ResEvals)
{ return Conserved(FALSE,ResEvals); } 

char	*rtf_typ::Conserved(BooLean Ignore,double *ResEvals)
{
	char    *conserved;
	NEW(conserved,nAlpha(AB)+3,char);
	for(Int4 r=1; r<=nAlpha(AB); r++){
		double V=-ResEvals[r];
		// if(V < TargetPval){ continue; } // unconserved.
		if(V < MinPval){ continue; } // unconserved = 0.
		// if(Ignore) { conserved[r]='w'; continue; } 	// AFN 12/22/05
		V = GetValue(V);
		Int4 bars=GetHistHeight(V);
		if(bars < 1) conserved[r]='w'; 
		else if(bars < bars_s_cut) conserved[r]='m'; 
		else conserved[r]='s'; 
	} return conserved;
}

Int4	rtf_typ::FixSetStatus(char *status,rtf_typ *rtfQ)
// Add '?' and '!' at positions missed by BPPS procedure.
// 
{
     double	ave,sum,max=-99999.,min=DBL_MAX;
     double	aveQ,sumQ,maxQ=-99999.,minQ=DBL_MAX;
     Int4	i,n,nQ,N;
     // 1. find the average Value 
     for(N=0,sum=0.0,n=0,sumQ=0.0,nQ=0,i=1; i <= length; i++){
	if(status[i] == '!'){
#if 1	
		// if PValue is not significant, then ignore these positions...
		// These can correspond to positions selected to be in the
		// superfamily by the BPPS procedure, but that don't align
		// with anything in the main set.
		if(PValue[i] <= 1.0) { status[i] = '*'; continue; }
#endif
		N++;
		sum += PValue[i]; n++; 
		sumQ += rtfQ->PValue[i]; nQ++; 
		if(max < PValue[i]) max=PValue[i];
		if(min > PValue[i]) min=PValue[i];
		if(maxQ < rtfQ->PValue[i]) maxQ=rtfQ->PValue[i];
		if(minQ > rtfQ->PValue[i]) minQ=rtfQ->PValue[i];
	}
     }
     // fprintf(stderr,"min = %g; minQ = %g\n",min,minQ);
     if(n == 0 || nQ == 0) return 0;
     ave = sum/(double)n;
     aveQ = sumQ/(double)nQ;
     for(i=1; i <= length; i++){
	if(status[i] != '!'){
#if 0
	  if(PValue[i] >= min && rtfQ->PValue[i] >= minQ){
		status[i]='!'; N++;
	  } else if(PValue[i] > min){ status[i]='?'; N++; }
#else
	  // if(PValue[i] > min){ status[i]='?'; N++; }
#if 0
	  if(PValue[i] >= ave){ status[i]='?'; N++; }
	  // if(PValue[i] >= aveQ){ status[i]='?'; N++; }
#else	// add feature to ignore a position... // AFN: 2/3/2009
	  if(Position == 0 && status[i] == '^'){
		// status[i]='*'; 
		IgnorePosition=TRUE;
		Position = i; binomial[i]=0.0;
	  } else if(PValue[i] >= ave){ status[i]='?'; N++; }
	  // if(PValue[i] >= aveQ){ status[i]='?'; N++; }
#endif
#endif
	}
     } return N;
}

char	*rtf_typ::DeriveStatus(char *statusSF,rtf_typ *rtfQ,char category)
// Determine which columns to select for the main-set category
// based on those selected from the superfamily category.
// '*' or anything else == ignored; '!' == key column.
{
     double	ave,sum,max=-99999.,min=DBL_MAX;
     double	aveQ,sumQ,maxQ=-99999.,minQ=DBL_MAX;
     double	stdev,stdevQ,d;
     Int4	i,n,nQ;
     char	*status;
     // 1. find the average Value 
     stdev=stdevQ=0.0;
     // Determine mean and standard deviation for Normal approximation.
     for(sum=0.0,n=0,sumQ=0.0,nQ=0,i=1; i <= length; i++){
	if(statusSF[i] == '!'){
		if(PValue[i] <= 1.0) { continue; } // see comments above...
		d = PValue[i]; stdev += d*d; sum += d; n++; 
		if(max < d) max=d;
		if(min > d) min=d;
		d = rtfQ->PValue[i]; stdevQ += d*d; sumQ += d; nQ++; 
		if(maxQ < d) maxQ=d;
		if(minQ > d) minQ=d;
	}
     }
     if(n == 0 || nQ == 0){
	NEW(status,length+3,char);
	for(i=1; i <= length; i++) status[i]='*';
	return status;
     }
     ave = sum/(double)n;
     stdev=sqrt((stdev/(double)n)-ave*ave);
     aveQ = sumQ/(double)nQ;
     stdevQ=sqrt((stdevQ/(double)nQ)-aveQ*aveQ);
     double OneSDs,OneSDsQ,M1SDs,M1SDsQ,M2SDs,M2SDsQ;
     OneSDs=ave+stdev; OneSDsQ=aveQ+stdevQ;
     M1SDs=ave-stdev; M1SDsQ=aveQ-stdevQ;	// 33 % below...
     M2SDs=ave-2*stdev; M2SDsQ=aveQ-2*stdevQ;	// 5 % below...
     NEW(status,length+3,char);
     if(category == 'M'){		//************** Main Set...
      for(i=1; i <= length; i++){
	if(statusSF[i] == '*' && PValue[i] > ave && rtfQ->PValue[i] > aveQ){
	// if(statusSF[i] == '*' && PValue[i] > OneSDs && rtfQ->PValue[i] > OneSDsQ){
		status[i]='!';
	} else status[i]='*';
      }
     } else if(category == 'I'){	//************** Intermediate Set...
      for(i=1; i <= length; i++){
	// if(statusSF[i] == '*' && PValue[i] > ave && rtfQ->PValue[i] > aveQ){
	// if(statusSF[i] == '*' && PValue[i] >= min && rtfQ->PValue[i] >= minQ){
	if(statusSF[i] == '*' && PValue[i] >= M1SDs && rtfQ->PValue[i] >= M1SDsQ){
	// if(statusSF[i] == '*' && PValue[i] >= M2SDs && rtfQ->PValue[i] >= M2SDsQ){
		status[i]='!';
		// fprintf(stdout,"%d: %g %g\n",i,PValue[i],rtfQ->PValue[i]);
	} else status[i]='*';
      }
     } else {				//************** Family Set...
       min=DBL_MAX; minQ=DBL_MAX;
       for(i=1; i <= length; i++){
	// if(statusSF[i] != '!' && PValue[i] > ave && rtfQ->PValue[i] > aveQ){
	if(statusSF[i] != '!' && PValue[i] > OneSDs && rtfQ->PValue[i] > OneSDsQ){
		status[i]='!';
		if(min > PValue[i]) min=PValue[i];
		if(minQ > rtfQ->PValue[i]) minQ=rtfQ->PValue[i];
	} else status[i]='*';
       }
     }
     return status;
}

double	rtf_typ::ComputeRelEntropy(double *observed, double *freq)
/* return the conserved residues in column "observed". */
// scale is the difference 
{
	Int4	r;
	double	total,p,q,RE=0.0,pseudo=0.01;

	if(observed != NULL){
	    for(total=0.0, r=1; r<=nAlpha(AB); r++){
		total += (double)observed[r] + pseudo;
	    }
	    for(r=1; r <= nAlpha(AB); r++){
		if(observed[r] > 0.01){
		   p = (observed[r]+pseudo)/total;
		   q = freq[r];
		   RE+= p*log(p/q);
		} 
	    } return RE;
	} else return 0.0;
}

