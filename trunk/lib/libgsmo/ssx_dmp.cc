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

#include "ssx_typ.h"
#include <math.h>
#include "blosum62.h"
// #include "dirichlet.h"

double  ssx_typ::DirichletLLR(char mode,FILE *efp)
/**********************************************************************
 DMS->Alpha == con
 DMS->alpha == con * P
 DMS->weight == mix.
 DMS->nDC == 30 here...
 **********************************************************************/
{
    Int4        n,l,i,j,t,T,s,ncol,b,x;
    Int4	blk,nblks=nTypeSites(SitesCMSA(CMA));
    double      MAP,weight=0,v,dd,d,total,*obs,*nonsite;
    UInt4	*wtcnt;

    /******** compute weight for number of blocks??? ********/
    for(weight=0.0, s=1; s<=N; s++){
        if(nBlksCMSA(CMA) == 1 && nSites(1,s,SitesCMSA(CMA)) == 0) continue;
#if 0   // new for computing column contributions.
        if(Blk && Col==0){ // for computing block contributions.
          T = LenSeqCMSA(s,CMA) + nblks - 1;
          for(t=1; t <= nblks; t++){
                if(t == Blk) continue;
                T -= LengthCMSA(t,CMA);
          } d=(double)Wt[s]/(double)wt_factor;
          if(T > (nblks-1)) weight += d*lnbico(T, nblks-1);
        } else
#endif
        {
          T = LenSeqCMSA(s,CMA) + nblks;
          for(t=1; t <= nblks; t++) T -= LengthCMSA(t,CMA);
          dd=(double)Wt[s]/(double)wt_factor;
          if(T > nblks) weight += dd*lnbico(T, nblks);
        }
    }
    if(efp) fprintf(efp,"weight = %.2f\n",weight);
    /******** end compute weight for number of blocks ********/
    /******** end NEW compute weight for number of blocks ********/
    NEW(obs,nAlpha(AB)+3, double); NEW(nonsite,nAlpha(AB)+3, double);
    for(b=0;b<=nAlpha(AB);b++){
	 nonsite[b]=((double)TtlWtCnts[b] + WtPseudo)/(double)wt_factor; 
    }
    /***** NULL MAP ******/
    for(MAP=0.0,total=0.0, b=1;b<= nAlpha(AB); b++) {
	MAP -= lngamma(nonsite[b]); total += nonsite[b];
    } MAP += lngamma(total);
    if(efp) fprintf(efp,"Null LogLikelihood = %.2f\n",MAP); dd=MAP;

    /***** SITES MAP *******/
    for(blk=1; blk <= nblks; blk++) {
	for(i=1; i <= LengthCMSA(blk,CMA); i++){
	    wtcnt=WtCnts[blk][i];
	    for(total=0.0,b=1;b<= nAlpha(AB); b++){
		obs[b]=(double)wtcnt[b]/(double)wt_factor; total += obs[b];
	        nonsite[b] -= obs[b];
	    } // Dirichlet mixture priors...
#if 1
	    d = DMS->ComputeLogLike(obs,efp);
	    MAP += d;
	    if(efp) fprintf(efp,"%d: MAP + %.3f = %.2f\n",i,d,MAP);
#else
	    double  tPrb=0.0,mix,con,a;
	    for(x=1;x<=30; x++){
		mix=Dirichlet30Mix[x]; con=Dirichlet30Conc[x];
		v = log(mix) + lgamma(con) - lgamma(total+con);
	        for(b=1;b <= nAlpha(AB); b++) {
		    a=con*Dirichlet30P[x][b]; 
		    v += lgamma(obs[b]+a) - lgamma(a);
		} // tPrb += mix*exp(v); 
		if(x==1) tPrb=v;
		else if(tPrb > v) tPrb +=log(1.0 + exp(v-tPrb));
	        else tPrb = v + log(1.0+exp(tPrb-v));
            	if(efp) fprintf(efp,"tPrb[%d][%d] = %.2f; v=%.2f; mix=%.2f\n",i,x,tPrb,v,mix);
	    } if(efp) fprintf(efp,"tPrb[%d] = %.2f\n",i,tPrb);
	    MAP += tPrb;
#endif
	}
    } 
    if(efp) fprintf(efp,"Sites logLikelihood = %.2f\n",MAP-dd);  dd=MAP;
    /***** NONSITES MAP *******/
    for(total=0.0,b=1;b<= nAlpha(AB); b++) {
	MAP += lngamma(nonsite[b]); total += nonsite[b];
    } MAP -= lngamma(total);
    if(efp) fprintf(efp,"Non-sites logLikelihood = %.2f\n",MAP-dd);
    free(nonsite); free(obs);

    /**************** return MAP *************/
    dd=NDL->IndelPenalty(0,mode,Wt);
    if(efp) fprintf(efp,"Indel Penalty = %.2f\n",dd); 
    return (MAP - weight + dd);
    // return MAP + dd;
}

double  ssx_typ::AdjstDirichletLLR(double targetWt,FILE *efp)
/**********************************************************************
 Adjust for differences in sequence weights; gp not included!
 **********************************************************************/
{
    assert(nBlksCMSA(CMA) == 1);
    Int4        n,l,i,j,t,T,s,r,ncol,b,blk=1;
    double      MAP,weight=0,v,dd,d,total,*obs,*nonsite;
    if(ReCalcSMX) CalcSMX();

    //============ Compute adjusted WtCnts =================
    double      sqwt=this->RtnWtNumSeqs(),Ratio=targetWt/sqwt;
    UInt4	x,*wtcnt,**aWtCnts,*aTtlWtCnts;
    NEWP(aWtCnts,LengthCMSA(1,CMA)+5,UInt4); NEW(aTtlWtCnts,nAlpha(AB)+5,UInt4);
    for(i=1; i <= LengthCMSA(1,CMA); i++){
      NEW(aWtCnts[i],nAlpha(AB)+5,UInt4);
      for(r=0; r<=nAlpha(AB); r++){
	x=WtCnts[1][i][r];
        if(x!=0){
                d=(double)x*Ratio; x=(UInt4) ceil(d-0.49);
        } aWtCnts[i][r]=x; aTtlWtCnts[r] += x;
     }
    }
#if 0	// add weights...
    for(weight=0.0, s=1; s<=N; s++){
        if(nBlksCMSA(CMA) == 1 && nSites(1,s,SitesCMSA(CMA)) == 0) continue;
        {
          T = LenSeqCMSA(s,CMA) + nblks;
          for(t=1; t <= nblks; t++) T -= LengthCMSA(t,CMA);
          dd=(double)Wt[s]/(double)wt_factor;
          if(T > nblks) weight += dd*lnbico(T, nblks);
        }
    }
    if(efp) fprintf(efp,"weight = %.2f\n",weight);
#endif
    //=============== Compute nonsite counts ===================
    NEW(obs,nAlpha(AB)+3, double); NEW(nonsite,nAlpha(AB)+3, double);
    for(b=0;b<=nAlpha(AB);b++){
	 nonsite[b]=((double)TtlWtCnts[b] + WtPseudo)/(double)wt_factor; 
    }
    /***** NULL MAP ******/
    for(MAP=0.0,total=0.0, b=1;b<= nAlpha(AB); b++) {
	MAP -= lngamma(nonsite[b]); total += nonsite[b];
    } MAP += lngamma(total);
    if(efp) fprintf(efp,"Null LogLikelihood = %.2f\n",MAP); dd=MAP;

    /***** SITES MAP *******/
    for(i=1; i <= LengthCMSA(blk,CMA); i++){
	    wtcnt=aWtCnts[i];
	    for(total=0.0,b=1;b<= nAlpha(AB); b++){
		obs[b]=(double)wtcnt[b]/(double)wt_factor; total += obs[b];
	        nonsite[b] -= obs[b];
	    } // Dirichlet mixture priors...
	    d = DMS->ComputeLogLike(obs,efp); MAP += d;
	    if(efp) fprintf(efp,"%d: MAP + %.3f = %.2f\n",i,d,MAP);
    } free(obs);
    if(efp) fprintf(efp,"Sites logLikelihood = %.2f\n",MAP-dd);  dd=MAP;
    /***** NONSITES MAP *******/
    for(total=0.0,b=1;b<= nAlpha(AB); b++) {
	MAP += lngamma(nonsite[b]); total += nonsite[b];
    } MAP -= lngamma(total);
    if(efp) fprintf(efp,"Non-sites logLikelihood = %.2f\n",MAP-dd);
    for(i=1; i <= LengthCMSA(blk,CMA); i++) free(aWtCnts[i]);
    free(aWtCnts); free(aTtlWtCnts); free(nonsite);

    /**************** return MAP *************/
    return (MAP);
}



