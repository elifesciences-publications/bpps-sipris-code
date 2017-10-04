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

csp_typ::csp_typ(double *obsFG,double *obsBG, double *frqFG, double *frqBG, a_type A,
						double min_res_freq,rst_typ *RST)
{ Init(obsFG,obsBG,frqFG,frqBG,A,min_res_freq,RST); }

void	csp_typ::Init(double *obsFG,double *obsBG, double *frqFG, double *frqBG,
	a_type A,double min_res_freq,rst_typ *RST)
{ 
	AB=A;
        MinResFreqCutoff=min_res_freq;
        // MinBlosum62Cutoff=-0.200; 
        MinBlosum62Cutoff=0.0; 
        AveBlosumScoreCutoff=0.10; 
	alpha=0.95;
	A0=9.0; B0=1.0;

	assert(RST); rst=RST;
	assert(obsFG); assert(frqBG);
	ObsFG=obsFG; ObsBG=obsBG;
	FreqFG=frqFG; // currently not used.
	FreqBG=frqBG;
	sets=0;
	ResEvals=0;
}

void	csp_typ::Free()
{
}

#define TRUNCATED_HIGHLIGHT_LOG10P  (-3.0)

//*******************************************************************************

void	csp_typ::cbp_dfs0(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
	double *obs, double q, double *freq, double total,double *min_Eval,
	double *min_p)
// cbp_dfs = cumulative binomial probability - depth first search
// Recursively find 'related' residue sets that are statistically significant via the Urn model.
{
	Int4	rx,r,sx;
	for(rx=rs+1; rx <= nAlpha(AB); rx++){
           if(obs[rx] > 2){	// need to see at least three of this residue type.
	     BooLean flag=TRUE;
	     // Similar residue sets are determined here...
	     if(rst){ //********** afn: 7/27/07 new residue set specification mode **********
		   sst_typ sst=SsetLet(rx);
		   sst_typ nsst=UnionSset(sst,sset);
		   if(!rst->IsLegalSet(nsst)){ flag=FALSE; break; }
	     } else {	//****** old method...
	       for(r=1; r < rx; r++){
		if(MemSset(r,sset) && blosum62[r][rx] <= MinBlosum62Cutoff){
			flag=FALSE; break; 
		}
	       }
	     }
	     if(flag){
	          sst_typ	sst=SsetLet(rx);
	          sst=UnionSset(sst,sset);	// add residue rx to sset for next stage in dfs.
		  double	p=Log10CBP(Obs+obs[rx], total, q+freq[rx]);
#if 1		// Attempt to change -c option...
		  double fract_in_set = (Obs+obs[rx])/total;
 		  if(fract_in_set < MinResFreqCutoff){
			if(p < TRUNCATED_HIGHLIGHT_LOG10P) p = TRUNCATED_HIGHLIGHT_LOG10P;		
		  }
#endif
	          for(r=1; r <= rx; r++){	// make sure that prob is the best found yet...
		    if(MemSset(r,sst) && p >= ResEvals[r]){ flag=FALSE; break; }
	          }
		  if(flag){	// if this set has the best prob found yet then save this group 
		   sx = findDSets(rx,sets);
		   if(sx != *ss){ *ss = linkDSets(*ss,sx,sets); }
		   if(p < *min_p) *min_p=p;
	           for(r=1; r <= rx; r++){
		     if(MemSset(r,sst) && p < min_Eval[r]) min_Eval[r]=p;
		   }
		  }
		  cbp_dfs0(rx,ss,sst,Obs+obs[rx],obs,q+freq[rx],freq,total,min_Eval,min_p);
	     }
	   }
	} 
}

// OLD cbp_dfs() return; pre-debug...
BooLean	csp_typ::cbp_dfs2(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
	double *obs, double q, double *freq, double total,double *min_Eval,
	double *min_p)
// cbp_dfs = cumulative binomial probability - depth first search
// Recursively find 'related' residue sets that are statistically significant via the Urn model.
// #define MIN_BLOSUM62_CUTOFF	(0.0)
// #define MIN_BLOSUM62_CUTOFF	(-0.050) // adds in NT, AT, VT as well.
// #define MIN_BLOSUM62_CUTOFF	(-0.200) // adds in QED and more...
{
	Int4	rx,r,sx;
	BooLean	rtnflag=FALSE;
	assert(MinResFreqCutoff > 0.0);
	for(rx=1; rx <= nAlpha(AB); rx++){
	   if(MemSset(rx,sset)) continue;
           if(obs[rx] > 2){
	     BooLean flag=TRUE;
	     double ave_score=0.0;
	     Int4 size=0;
	     if(rst){ //********** afn: 7/27/07 new residue set specification mode **********
		  sst_typ sst=SsetLet(rx);
		  sst_typ nsst=UnionSset(sst,sset);
		  if(!rst->IsLegalSet(nsst)){ flag=FALSE; break; }
		  // ignore size and ave_score in this mode...
	     } else {	//****** old method...
	       for(r=1; r <= nAlpha(AB); r++){
	        if(r == rx) continue;
		// If any edge from rx into the current group of related residues
		// (any r element of sset) is below the cutoff then don't use residue rx.
		if(MemSset(r,sset)){
		   if(blosum62[r][rx] <= MinBlosum62Cutoff){ flag=FALSE; break; }
		   else { ave_score += blosum62[r][rx]; size++; }
		}
	       }
	     }
	     if(flag){
	        if(size > 0) ave_score = ave_score/(double)size;
		else ave_score = -9.99;
		if(TRUE || ave_score > AveBlosumScoreCutoff){ // go further only if rx is good match to current set...
			// WARNING: THIS MAY NOT BE TRUE LATER ON...
	          sst_typ	sst=SsetLet(rx);
	          sst=UnionSset(sst,sset);	// add residue rx to sset for next stage in dfs.
		  double	p=Log10CBP(Obs+obs[rx], total, q+freq[rx]);
	          for(r=1; r <= nAlpha(AB); r++){		
		    // make sure that p is the best found yet...including where fract_in_set < MinResFreqCutoff
		    if(r == rx) continue;	// min_Eval[rx] is not yet set...
		    if(MemSset(r,sst) && p >= ResEvals[r]){ flag=FALSE; break; }
		    if(MemSset(r,sst) && p >= min_Eval[r]){ flag=FALSE; break; }
	          }
		  double fract_in_set = (Obs+obs[rx])/total;
		  if(flag){	// set ResEval[r] = p for checking at lower depth...
				// reset later by ComputeBinomial()
#if 1		    // Attempting to eliminate highlight but show res_freq in *csp file.
 		    if(fract_in_set < MinResFreqCutoff){
		      	// fprintf(stderr,"Log10CBP() = %f\n",p);
		      	if(p < TRUNCATED_HIGHLIGHT_LOG10P) p = TRUNCATED_HIGHLIGHT_LOG10P;
		    }
#endif
#if 0
		    for(r=1; r <= nAlpha(AB); r++){ if(MemSset(r,sst)) ResEvals[r]=p; }
#else
		    for(r=1; r <= nAlpha(AB); r++){
			if(MemSset(r,sst)) ResEvals[r]=p; 
		    }
#endif
		  } // ^This determines whether enlarging residue set 
			// actually helps improve the (rejected) E-value
		  if(flag){	// if current set has the best prob found yet 
				// then save this group 
		      // Found a set that meets the criteria for acceptance.
		      sx = findDSets(rx,sets);
		      if(sx != *ss){ *ss=linkDSets(*ss,sx,sets); }
		      rtnflag=TRUE;
		      if(p < *min_p) *min_p=p;
	              for(r=1; r <= nAlpha(AB); r++){
		        if(MemSset(r,sst) && p < min_Eval[r]) min_Eval[r]=p;
		      }
		  }
		  if(0) fprintf(stderr,"%c (obs=%d)",AlphaChar(rx,AB),(Obs+obs[rx]));
		  if(cbp_dfs2(rx,ss,sst,Obs+obs[rx],obs,q+freq[rx],freq,total,min_Eval,min_p)) {
		     // fprintf(stderr,"************** FOUND AFTER THE FACT ****************\n");
		     // this IS used a lot...
		     // found to be significant afterwards, so merge sets.
		     sx = findDSets(rx,sets);
		     if(sx != *ss){ *ss=linkDSets(*ss,sx,sets); }
		     rtnflag=TRUE;
		  }
		  if(0) fprintf(stderr,"Exiting DFS \n");
		}
	     }
	   }
	} return rtnflag;
}

// NEW cbp_dfs() routine; debugged...


#define RTF_TYP_USE_STD_BACKGROUND	1

BooLean	csp_typ::ComputeBinomial(double cutoff, double *pvalue,double *resEvals)
{ 
	ResEvals=resEvals;
	if(MinResFreqCutoff > 0.0){
		return ComputeBinomial2(ObsFG, cutoff,FreqBG, pvalue); 
	} else return ComputeBinomial2(ObsFG, cutoff,FreqBG, pvalue); 
}

BooLean	csp_typ::ComputeBinomial0(double *observed, double cutoff, double *freq, double *pvalue)
/* return the conserved residues in column "observed". */
// scale is the difference 
{
	Int4	r,r2,s1,s2;
	double	total,p,q,min_p,obs[30];
	double	min_Eval[30];	// minimum E-values for residues.
	BooLean	OkayToUse=FALSE;

	assert(cutoff <= 0.0);
	// cutoff = cutoff - 2.0;  
	*pvalue=min_p=0.0;
	if(observed != NULL){
	    BooLean not_set[30];
	    sets = DSets(nAlpha(AB));
#if RTF_TYP_USE_STD_BACKGROUND
	    for(total=0.0, r=1; r<=nAlpha(AB); r++){
		not_set[r] = TRUE;
		total += observed[r]; 
		obs[r]=observed[r];
	    }
#else	// WARNING: THIS CODE IS CURRENTLY MESSED UP AND WON'T WORK!!!
	    for(total=0.0, r=1; r<=nAlpha(AB); r++){
		not_set[r] = TRUE;
		total += observed[r]; 
	    }
	    double frq[30],sum_frq=0.0;
	    for(r=1; r<=nAlpha(AB); r++){
		frq[r] = (total*freq[r])+((double)(nAlpha(AB))*(obs[r]/total));
		// frq[r] = (total*freq[r])+(obs[r]/total); // one pseudocount.
		sum_frq+=frq[r];
	    }
	    for(r=1; r<=nAlpha(AB); r++){
		frq[r] = frq[r]/sum_frq;
		obs[r]=observed[r]+frq[r];  // pseudocounts...
	    }
#endif
	    /** 1. find conserved residues **/
	    for(r=1; r <= nAlpha(AB); r++){
		if(obs[r] > 0.01) {
#if RTF_TYP_USE_STD_BACKGROUND
		   q = freq[r];
#else
		   q = frq[r];
#endif
		   p = Log10CBP(obs[r], total, q);
		   if(MinResFreqCutoff > 0.0){	//************************************
		     // ResEvals == the evalue for a set with only one residue 
		     if((obs[r]/total) >= MinResFreqCutoff){
		   	if(ResEvals) ResEvals[r]=min_Eval[r]=p;
		   	if(p < min_p) min_p = p;
		     } else {
#if 1
			if(p < TRUNCATED_HIGHLIGHT_LOG10P) p = TRUNCATED_HIGHLIGHT_LOG10P; 
		   	if(ResEvals) ResEvals[r]=p; min_Eval[r]=p;
#endif
		     }
		   } else {			//************************************
		     // ResEvals == the evalue for a set with only one residue 
		     if(ResEvals) ResEvals[r]=min_Eval[r]=p;
		     if(p < min_p) min_p = p;
		   }				//**********************************
		} else if(ResEvals){ ResEvals[r]=min_Eval[r]=0.0; }
	    }
	    if(total == 0.0){ *pvalue=0.0; NilDSets(sets); return FALSE; }
            /** 2. find conserved related residue pairs **/
            // for(r=1; r < nAlpha(AB); r++){	// old: r < nAlpha(AB)??!!!
            for(r=1; r <= nAlpha(AB); r++){
                if(obs[r] > 2){
		  s1 = findDSets(r,sets);
#if RTF_TYP_USE_STD_BACKGROUND
		  if(MinResFreqCutoff > 0.0){ //*************** NEW ****************
			if(0) fprintf(stderr,"Entering DFS with %c (%.3f)",
					AlphaChar(r,AB),MinResFreqCutoff);
			BooLean UsePos=cbp_dfs2(r,&s1,SsetLet(r),obs[r],obs,freq[r],
				freq,total,min_Eval,&min_p);
#if 1
			if(UsePos) OkayToUse=TRUE;
#endif
			if(0) fprintf(stderr,"\n");
		  } else {
			cbp_dfs0(r,&s1,SsetLet(r),obs[r],obs,freq[r],freq,total,
				min_Eval,&min_p);
			OkayToUse=TRUE;
		  }
#else
		  cbp_dfs2(r,&s1,SsetLet(r),obs[r],obs,frq[r],frq,
			total,min_Eval,&min_p);
#endif
                }
	    }
            dh_type dH=dheap(32,4);
            for(r=1; r <= nAlpha(AB); r++){
	     //********************** NEW OPTION *********
	     if(MinResFreqCutoff > 0.0){
	      if(not_set[r]){
	       if(min_Eval[r] < 0.0){
		s1 = findDSets(r,sets);
		 for(r2=1; r2 <= nAlpha(AB); r2++){
		   s2 = findDSets(r2,sets);	
		   // find all residues in the same disjoint set with r.
		   if(s1 == s2){ insrtHeap(r2,((keytyp)min_Eval[r2]),dH); }
		 }
                 for(double tiny=0.0; ((r2=delminHeap(dH)) != 0); tiny+=0.00001){
		    ResEvals[r2]=min_Eval[r2]+tiny; 
		    not_set[r2] = FALSE;
		 }
	       } else { ResEvals[r]=0.0; }
	      }
	     } else { //********************************************
	      if(not_set[r] && ResEvals[r] < 0.0){
		s1 = findDSets(r,sets);
		for(r2=1; r2 <= nAlpha(AB); r2++){
		  s2 = findDSets(r2,sets);	
		  // find all residues in the same disjoint set with r.
		  if(s1 == s2){ insrtHeap(r2,((keytyp)ResEvals[r2]),dH); }
		}
                for(double tiny=0.0; ((r2=delminHeap(dH)) != 0); tiny+=0.00001){
		    ResEvals[r2]=min_Eval[r2]+tiny; 
		    not_set[r2] = FALSE;
		}
	      }
	     } 	// ***** End of NEW OPTION if/else statement *********
	    } Nildheap(dH); NilDSets(sets);
	    *pvalue=min_p;
	    if(min_p > cutoff){ return FALSE; }
	    else return OkayToUse;
	    // else return TRUE;
	} else return FALSE;
}

BooLean	csp_typ::cbp_dfs(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
	double *obs, double q, double *freq, double total,double *min_Eval,
	double *min_p)
// cbp_dfs = cumulative binomial probability - depth first search
// Recursively find 'related' residue sets that are statistically significant via the Urn model.
{
	Int4	rx,r,sx;
	BooLean	rtnflag=FALSE;
	// assert(MinResFreqCutoff > 0.0);
	for(rx=1; rx <= nAlpha(AB); rx++){
	   if(MemSset(rx,sset)) continue;
           if(obs[rx] > 2){
	     BooLean flag=TRUE;
	     Int4 size=0;
	     assert(rst);
	     sst_typ sst=SsetLet(rx);
	     sst_typ nsst=UnionSset(sst,sset);
	     if(!rst->IsLegalSet(nsst)){ flag=FALSE; }
	     if(flag){
	          sst_typ	sst=SsetLet(rx);
	          sst=UnionSset(sst,sset);	// add residue rx to sset for next stage in dfs.
		  double	p=Log10CBP(Obs+obs[rx], total, q+freq[rx]);
	          for(r=1; r <= nAlpha(AB); r++){		
		    // make sure that p is the best found yet
		    // ...including where fract_in_set < MinResFreqCutoff
		    if(r == rx) continue;	// min_Eval[rx] is not yet set...
		    if(MemSset(r,sst) && (p >= ResEvals[r] || p >= min_Eval[r])){
				flag=FALSE; break; 
		    }
	          }
		  if(flag){	// set ResEval[r] = p for checking at lower depth...
				// reset later by ComputeBinomial()
		    double fract_in_set=(Obs+obs[rx])/total;
 		    if(fract_in_set < MinResFreqCutoff){
		        // Attempting to eliminate highlight but show res_freq in *csp file.
		      	// fprintf(stderr,"Log10CBP() = %f\n",p);
		      	if(p < TRUNCATED_HIGHLIGHT_LOG10P) p = TRUNCATED_HIGHLIGHT_LOG10P;
		    }
		    for(r=1; r <= nAlpha(AB); r++){ if(MemSset(r,sst)) ResEvals[r]=p; }
		  } // ^This determines whether enlarging residue set 
			// actually helps improve the (rejected) E-value
		  if(flag){	// if current set has the best prob found yet then save this group 
		      // Found a set that meets the criteria for acceptance.
		      sx = findDSets(rx,sets);
		      if(sx != *ss){ *ss=linkDSets(*ss,sx,sets); }
		      rtnflag=TRUE;
		      if(p < *min_p) *min_p=p;
	              for(r=1; r <= nAlpha(AB); r++){
		        if(MemSset(r,sst) && p < min_Eval[r]) min_Eval[r]=p;
		      }
		  }
		  if(cbp_dfs(rx,ss,sst,Obs+obs[rx],obs,q+freq[rx],freq,total,min_Eval,min_p)) {
		     // fprintf(stderr,"************** FOUND AFTER THE FACT ****************\n");
		     // this IS used a lot...
		     // found to be significant afterwards, so merge sets.
		     sx = findDSets(rx,sets);
		     if(sx != *ss){ *ss=linkDSets(*ss,sx,sets); }
		     rtnflag=TRUE;
		  }
	     }
	   }
	} return rtnflag;
}

BooLean	csp_typ::ComputeBinomial2(double *observed, double cutoff, double *freq, double *pvalue)
/* return the conserved residues in column "observed". */
// scale is the difference 
{
	Int4	r,r2,s1,s2;
	double	total,p,min_p,obs[30];
	double	min_Eval[30];	// minimum E-values for residues.
	BooLean	OkayToUse=FALSE;

	assert(ResEvals);
	// assert(MinResFreqCutoff > 0.0);
	assert(cutoff <= 0.0);
	// cutoff = cutoff - 2.0;  
	*pvalue=min_p=0.0;
	if(observed != NULL){
	    BooLean not_set[30];
	    sets = DSets(nAlpha(AB));
	    for(total=0.0, r=1; r<=nAlpha(AB); r++){
		not_set[r] = TRUE;
		total += observed[r]; 
		obs[r]=observed[r];
	    }
	    if(total == 0.0){ *pvalue=0.0; NilDSets(sets); return FALSE; }
	    /** 1. find conserved residues **/
	    for(r=1; r <= nAlpha(AB); r++){
		if(obs[r] > 0.01) {
		   p = Log10CBP(obs[r],total,freq[r]);
		   // ResEvals == the evalue for a set with only one residue 
		   if((obs[r]/total) >= MinResFreqCutoff){
		   	ResEvals[r]=min_Eval[r]=p;
		   	if(p < min_p) min_p = p;
		   } else {
			if(p < TRUNCATED_HIGHLIGHT_LOG10P) p = TRUNCATED_HIGHLIGHT_LOG10P; 
		   	ResEvals[r]=p; min_Eval[r]=p;
		   }
		} else { ResEvals[r]=min_Eval[r]=0.0; }
	    }
            /** 2. find conserved related residue pairs **/
            // for(r=1; r < nAlpha(AB); r++){	// old: r < nAlpha(AB)??!!!
            for(r=1; r <= nAlpha(AB); r++){
                if(obs[r] > 2){
		  s1 = findDSets(r,sets);
		  BooLean UsePos=cbp_dfs(r,&s1,SsetLet(r),obs[r],obs,freq[r],
				freq,total,min_Eval,&min_p);
		  if(UsePos) OkayToUse=TRUE;
                }
	    }
            dh_type dH=dheap(32,4);
            for(r=1; r <= nAlpha(AB); r++){
	     if(not_set[r]){
	       if(min_Eval[r] < 0.0){
		s1 = findDSets(r,sets);
		 for(r2=1; r2 <= nAlpha(AB); r2++){
		   s2 = findDSets(r2,sets);	
		   // find all residues in the same disjoint set with r.
		   if(s1 == s2){ insrtHeap(r2,((keytyp)min_Eval[r2]),dH); }
		 }
                 for(double tiny=0.0; ((r2=delminHeap(dH)) != 0); tiny+=1e-20){
		    ResEvals[r2]=min_Eval[r2]+tiny; 
		    not_set[r2] = FALSE;
		 }
	       } else { ResEvals[r]=0.0; }
	     }
	    } Nildheap(dH); NilDSets(sets);
	    *pvalue=min_p;
	    if(min_p > cutoff){ return FALSE; }
	    else return OkayToUse;
	} else return FALSE;
}


