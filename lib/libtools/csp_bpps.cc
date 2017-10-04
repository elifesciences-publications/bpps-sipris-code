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

#include "csp_typ.h"
#include "blosum62.h"

#define TRUNCATED_HIGHLIGHT_LOG10P  (-3.0)

double csp_typ::CalcBiNomBPPS(double n, double m)
// This model corresponds to Jun's BPPS model...
// n = matching FG; m = matching BG
// N-n = mis-match FG; M-m = mismatch BG.
{
	assert(N >= n && M >= m);
	double		b = 1.0;	// pseudocounts
	// double		m1=m+b,n1=n+b, m2=M-m+b,n2=N-n+b;
	double		m1=m,n1=n, m2=M-m,n2=N-n;

/************************* updateSumXI **************************/

	// double	theta_x=m1/(m1+m2);
	double alpha0=(double)A0/(double)(A0 + B0);
	double	theta_x=(m1+b)/(m1+m2+b+b);
	double	SumXI= n1*(alpha0/(alpha0+((1-alpha0)*theta_x))); // from Equation [3].
	assert(SumXI < n1);
/****************************************************************/
/************************ updateAlpha ***************************/

/****************************************************************/
/************************ updateTheta ***************************/
	// get theta_j
	double	mat = n1 + m1 + b;
        double  sum = mat + n2 + m2 + b;
	theta_x = mat/sum; 	// for null map below
	assert(mat > SumXI);
	mat -= SumXI;	// subtract out FG non-contamination which doesn't count.
	sum -= SumXI;	// subtract out FG non-contamination which doesn't count.
	double	theta_j = mat/sum; 
	assert(theta_j > 0.0 && theta_j < 1.0);
#if 0
	fprintf(stderr,"theta=%f; SumXI=%.2f; n1=%.2f; n2 = %.2f; m1=%.2f; m2 = %.2f; mat=%f; sum=%f; ",
		theta_j,SumXI,n1,n2,m1,m2,mat,sum);
#endif
/****************************************************************/
	double	submap=0.0;
	submap += n1*log((1-alpha)*theta_j + alpha); // delta_j == 1.
	submap += n2*log((1-alpha)*(1.0-theta_j));	  // delta_j == 0.
	submap += m1*log(theta_j);	
	submap += m2*log((1.0-theta_j));
	// submap += logRho;
	// put everything in background for NullMAP...
	double nullsubmap = ((n1+m1)*log(theta_x) + (n2+m2)*log((1.0-theta_x)));
	// submap += logOneMinusRho;
	double lpr=submap-nullsubmap;
#if 0
	if(n == N || lpr > -802 && lpr < -801){
	  fprintf(stderr,
		"SumXI=%g; n1=%.1f; n2=%.1f; m1=%.1f; m2=%.1f;map=%g-%g(th_x=%g; th_j=%g)\n",
		SumXI,n1,n2,m1,m2,submap,nullsubmap,theta_x,theta_j);
	}
	fprintf(stderr," lpr=%f\n",lpr);
#endif
	if(isfinite(lpr)) return -lpr;
	else return 0.0;
}

void	csp_typ::bnml_BPPS_dfs(Int4 rs, Int4 *ss, sst_typ sset, double n, 
		double m, double *min_Eval,double *min_p)
// bnml_BPPS_dfs = BPPS log-probability ratio - depth first search
// Recursively find 'related' residue sets that are significant via the BPPS model.
{
	Int4	rx,r,sx;
	double	p,q,P;

	// PutDSets(stderr,sets);
	for(rx=rs+1; rx <= nAlpha(AB); rx++){
           if(ObsFG[rx] >= 2){	// need to see at least two of this residue type.
	     BooLean flag=TRUE;
	     // Similar residue sets are determined here...
	     sst_typ sst=SsetLet(rx);
	     sst_typ nsst=UnionSset(sst,sset);
	     //if(!rst->IsLegalSet(nsst)){ flag=FALSE; break; }
	     if(!rst->IsLegalSet(nsst)){ flag=FALSE; }
	     if(flag){
	          sst_typ	sst=SsetLet(rx);
	          sst=UnionSset(sst,sset);	// add residue rx to sset for next stage in dfs.
		  q=(n+ObsFG[rx]+1.0)/(N + 2.0);
		  p=(m+ObsBG[rx]+1.0)/(M + 2.0);
		  assert(alpha > 0.0);
		  // P=CalcBiNomBPPS(n+ObsFG[rx], m+ObsBG[rx]);
		  if(q > p) P=CalcBiNomBPPS(n+ObsFG[rx],m+ObsBG[rx]); else P = 0.0;
#if 0		// Attempt to change -c option...
		  double fract_not_set = (n+ObsFG[rx])/N;
 		  if(fract_not_set < MinResFreqCutoff){
			if(P < TRUNCATED_HIGHLIGHT_LOG10P) P = TRUNCATED_HIGHLIGHT_LOG10P;
		  }
#endif
	          for(r=1; r <= rx; r++){	// make sure that prob is the best found yet...
		    if(MemSset(r,sst) && P >= ResEvals[r]){ flag=FALSE; break; }
	          }
		  if(flag){	// if this set has the best prob found yet then save this group 
		   sx = findDSets(rx,sets);
		   if(sx != *ss){ *ss = linkDSets(*ss,sx,sets); }
		   if(P < *min_p) *min_p=P;
	           for(r=1; r <= rx; r++){ if(MemSset(r,sst) && P < min_Eval[r]) min_Eval[r]=P; }
		  }
		  if(0) fprintf(stderr,"Entering DFS(%d) with %c (n=%g;m=%g)(p=%g)\n",(UInt4)sst,
		  		AlphaChar(rx,AB),n+ObsFG[r],m+ObsBG[r],*min_p);
		  bnml_BPPS_dfs(rx,ss,sst,n+ObsFG[rx],m+ObsBG[rx],min_Eval,min_p);
		  if(0) fprintf(stderr,"Exiting DFS(%d)\n",(UInt4)sst);
	     }
	   }
	} 
}

BooLean	csp_typ::ComputeBinomialBPPS(double cutoff, double *pvalue,double *resEvals)
/* return the conserved residues in this column. */
// scale is the difference 
{
	Int4	r,r2,s1,s2;
	double	n,m,P,min_p,p,q;
	double	min_Eval[30];	// minimum E-values for residues.
	BooLean	OkayToUse=FALSE;

	assert(cutoff <= 0.0);
	assert(ObsBG); 
        if(MinResFreqCutoff > 0.0){ 
		print_error("MinResFreqCutoff not implemented for ComputeBinomialBPPS( )");
	}
	*pvalue=min_p=0.0;
	ResEvals=resEvals;
	BooLean found_set[30];
	sets = DSets(nAlpha(AB));
	//	PutDSets(stderr,sets);
	for(N=M=0.0, r=1; r<=nAlpha(AB); r++){
		found_set[r] = FALSE;
		N+= ObsFG[r]; M+= ObsBG[r]; 
		ResEvals[r]=0.0;
	}
	/** 1. find single conserved residues **/
	// fprintf(stderr,"DEBUG 0\n");
	for(r=1; r <= nAlpha(AB); r++){
		if(ObsFG[r] >= 2) {
		  // fprintf(stderr,"DEBUG 0a\n");
		   q=(ObsFG[r]+1.0)/(N + 2.0);
		   p=(ObsBG[r]+1.0)/(M + 2.0);
		   double P;
		   assert(alpha > 0.0);
		   // P=CalcBiNomBPPS(ObsFG[r],ObsBG[r]);
		   if(q > p) P=CalcBiNomBPPS(ObsFG[r],ObsBG[r]); else P = 0.0;
		   // ResEvals == the evalue for a set with only one residue 
		   if(ResEvals) ResEvals[r]=min_Eval[r]=P;
		   if(P < min_p) min_p = P;
		} else if(ResEvals){ ResEvals[r]=min_Eval[r]=0.0; }
	}
	if(N == 0.0){ *pvalue=0.0; NilDSets(sets); return FALSE; }
	if(M == 0.0){ *pvalue=0.0; NilDSets(sets); return FALSE; }
        /************************ 2. find conserved related residue pairs *******************/
        for(r=1; r <= nAlpha(AB); r++){
		if(ObsFG[r] > 2) {
		  s1 = findDSets(r,sets); // PutDSets(stderr,sets);
		  if(0) fprintf(stderr,"Entering DFS0 with %c (%g; %g)\n",
					AlphaChar(r,AB),ObsFG[r],ObsBG[r]);
		  bnml_BPPS_dfs(r,&s1,SsetLet(r),ObsFG[r],ObsBG[r],min_Eval,&min_p);
		  if(0) fprintf(stderr,"Exiting DFS0\n");
		  OkayToUse=TRUE;
                }
	}
	BooLean debug_it=FALSE;
	// if(min_p < -1600.0 && min_p > -1610.0) debug_it=TRUE;
        dh_type dH=dheap(32,4);
	// if(debug_it) PutDSets(stderr,sets);
	// PutAlpha(stderr,AB);
        for(r=1; r <= nAlpha(AB); r++){
	      if(debug_it) 
		   fprintf(stderr,
			"Entering findDSets 2 with r= %c; log-eval = %g = %g; log-p=%g\n",
				AlphaChar(r,AB),ResEvals[r],min_Eval[r],min_p);
	      // if(!found_set[r] && ResEvals[r] < 0.0)
	      // if(!found_set[r] && ObsFG[r] >= 2.0)
	      if(!found_set[r] && min_p < 0.0)
	      {
		s1 = findDSets(r,sets);
		for(r2=1; r2 <= nAlpha(AB); r2++){
		  s2 = findDSets(r2,sets);	
		  // find all residues in the same disjoint set with r.
		  if(s1 == s2){
		        if(debug_it) fprintf(stderr,"r2=%c (%g)\n",AlphaChar(r2,AB),ResEvals[r]);
			insrtHeap(r2,((keytyp)ResEvals[r2]),dH); 
		  }
		}
		Int4	old_r2=0;
                for(double tiny=0.0; ((r2=delminHeap(dH)) != 0); tiny+=0.00001){
		    if(old_r2 == r2){
			PutDSets(stderr,sets);
			PutHeap(stderr, dH);
		    	fprintf(stderr,"r2=%c (%g)\n",AlphaChar(r2,AB),ResEvals[r]);
			assert(old_r2 != r2);
		    }
		    ResEvals[r2]=min_Eval[r2]-tiny; 
		    found_set[r2] = TRUE;
		    old_r2=r2;
		}
	      }
	      if(debug_it) fprintf(stderr,"Exiting findDSets 2\n");
	} Nildheap(dH); NilDSets(sets);
	*pvalue=min_p;
	if(min_p > cutoff){ return FALSE; }
	else return OkayToUse;
}

