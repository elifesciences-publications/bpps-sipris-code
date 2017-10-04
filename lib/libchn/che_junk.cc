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

#include "che_typ.h"

BooLean che_typ::SSetOkay(Int4 s, sst_typ sst,FILE *fp)
// If the FG set is less than the BG set for any residues return FALSE.
{
        double  d,bg,fg;
        UInt4   f,b,tb,tf;
        if(fp) fprintf(fp,"%d: ",s);
        for(Int4 r=1; r <= nAlpha(AB); r++){
            if(MemSset(r,sst)){
                tf = TotalResWt[1][s]; tb = TotalResWt[2][s];
                assert(tf > 0 && tb > 0);
                f=ResWt[1][s][r]; b=ResWt[2][s][r];
                if(fp==0){
                   if(f == 0) return FALSE;
                   if(b == 0) continue;
                   fg=(double)f/(double)tf; assert(fg >= 0.0 && fg <= 1.0);
                   bg=(double)b/(double)tb; assert(bg >= 0.0 && bg <= 1.0);
                   if((fg/bg) <= 0.20) return FALSE;
                } else {
                   if(f == 0) d=0.0;
                   else if(b == 0) d=1000.0;
                   else {
                     fg=(double)f/(double)tf; assert(fg >= 0.0 && fg <= 1.0);
                     bg=(double)b/(double)tb; assert(bg >= 0.0 && bg <= 1.0);
                     d=fg/bg;
                   } fprintf(fp,"%c%.2f ",AlphaChar(r,AB),d);
                }
            }
        } if(fp) fprintf(fp,"\n");
        return TRUE;
}

BooLean che_typ::ResetColumnJ(Int4 j, Int4 r, double prior_rho, Int4 s,sst_typ *sstJ)
// Update the csq used to restrain the sst.
{
        assert(j >= 1 && j <= LenSeq(Query));
        EqSeq(j,r,Query);
        if(XedSeq(Query)) EqXSeq(j,r,Query); 
	double *RhoJ=GetJthRhoCategoricalPriors(0, prior_rho, sstJ, AB);
        bpps_typ *pps=this->BPPS();  pps->UpdateRhoJ(j,RhoJ); free(RhoJ);
	if(s > 0){ 
		assert(sstJ[s] != 0);
		assert(MemSset(r,sstJ[s])); pps->AddColumn(j,sstJ[s]);
		if(sstJ[s] == 0) pps->RmColumn(j); 
	} else pps->RmColumn(j); 
}

BooLean	che_typ::SampleAnyPattern(Int4 j,sst_typ **sst, BooLean	**skip, unsigned char *NumResSets, 
			double prior_rho, double temperature)
// Sample a pattern irrespective of Query constraint...
// Call as: if(che[n]->SampleAnyPattern(j,rst->LegalResSets(),rst->EqualsPrevious,
//			rst->NumResSets,temperature)...
// resSSet[r][rset] == for pattern rset for residue r 
{
	double	map,*submap,**deltaLPR,bst_map=-99999999999999.9; 
	Int4	s,num_pttrns=0,old_s=0,bst_s;
	sst_typ	bst_sst=0,*sstJ,old_sst;
	unsigned char r,r_bst=0,r_old=ResSeq(j,Query);
	BooLean	**Skip;
	for(s=1; sst[r_old][s]; s++){ if(sst[r_old][s] == old_sst) old_s=s; break; }

	if(pps->NumColumns() >= MaxNumColumns) return FALSE;
	old_sst = pps->RtnSST(j); 
	NEWP(deltaLPR, nAlpha(AB)+ 3,double);
	NEWP(Skip, nAlpha(AB)+ 3,BooLean);
 	for(r=1; r <= nAlpha(AB); r++){
	   sstJ=sst[r];
	   for(s=1; sstJ[s]; s++) num_pttrns++;
	   NEW(deltaLPR[r], num_pttrns + 3,double);
	   NEW(Skip[r], num_pttrns + 3,BooLean);
	   for(s=1; sstJ[s]; s++){
		if(skip[r][s]){ Skip[r][s]=TRUE; continue; } // don't resampling same pattern.
		if(!SSetOkay(j,sstJ[s])){Skip[r][s]=TRUE;  continue; }  // This pattern can't work.
		this->ResetColumnJ(j,r,prior_rho,s,sstJ);
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
		map = submap[0];
		deltaLPR[r][s]=submap[j];
		if(submap[j] > bst_map) // Insist on positive sub_map at column position...
		   { bst_map=submap[j]; bst_sst=sstJ[s]; r_bst=r; bst_s=s; }
	   }
	}
	if(bst_sst==0){ 
		// this->ResetColumnJ(j,r_old,prior_rho,old_s,sst[r_old]);
		this->ResetColumnJ(j,r_old,prior_rho,0,sst[r_old]);
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
 		for(r=1; r <= nAlpha(AB); r++) free(deltaLPR[r]); 
    	  	free(deltaLPR); return FALSE; 
	} else {
	  //**************** Sample a pattern ****************
	  if(temperature > 0.0){		// else use bst_sst;
	    double d,y,sum=0.0;
	    assert(temperature <= 300.0);
 	    for(r=1; r <= nAlpha(AB); r++){
	       sstJ=sst[r];
	       for(s=1; sstJ[s]; s++){
		   if(Skip[r][s]) continue;
		   d=deltaLPR[r][s]-bst_map; 
		   assert(d <= 0.0);
		   if(d < -100){ d = 0.0; deltaLPR[r][s]=0.0; }
		   else {
			d = exp(d);
			if(temperature == 300.0){ deltaLPR[r][s]=d; sum += d; }
			else {
			  y=(300.0/temperature);
			  assert(!(d==0 && y <= 0)); assert(d >= 0);
			  if(d == 0) deltaLPR[r][s]=0;
			  else {
			    d = pow(d,y);
			    // assert(isfinite(d)); this is failing!!!
			    if(!isfinite(d)) assert(!"isfinite failed");
			    deltaLPR[r][s]=d; sum += d;
			  }
			}
		   }
	       }
	    }
	    double rand = sum * ((double) Random()/(double) RANDOM_MAX);        // random number...
 	    for(sum=0.0,r_bst=0,r=1; r <= nAlpha(AB); r++){
	      sstJ=sst[r];
	      for(s=1; sstJ[s]; s++){
		   if(Skip[r][s]) continue;
		   sum += deltaLPR[r][s];
		   if(sum >= rand){ bst_sst=sstJ[s]; bst_map=deltaLPR[r][s]; r_bst=r; bst_s=s; break; }
	      } if(r_bst > 0) break;
	    } assert(sum >= rand);
    	  }
 	  for(r=1; r <= nAlpha(AB); r++){ free(deltaLPR[r]); free(Skip[r]);  }
	  free(deltaLPR); free(Skip);
	  //**************** End of pattern sampling ****************
	  if(Verbose){
	   fprintf(stderr,"Redefined pattern set: '");
	   PutSST(stderr,old_sst,AB); fprintf(stderr,"%d' -> '",j);
	   PutSST(stderr,bst_sst,AB); fprintf(stderr,"%d'\n",j); 
	  }
	  this->ResetColumnJ(j,r_bst,prior_rho,bst_s,sst[r_bst]);
	  submap = pps->SubLPR(ResWt[2],ResWt[1]);
	  return TRUE; 
	}
}

