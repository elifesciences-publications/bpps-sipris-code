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

// C Code from Stephen Altchsul...
/******************************************************

   Enclosed below is some C code to produce scores the way PSI-BLAST did at one point.  You need to provide a vector of counts for the 20 amino acids, a vector of amino acid background frequencies, and a standard integral-valued substitution matrix.  You also need to provide a scale parameter and a pseudocount parameter.

   The program is initialized by first calling mcalc() with the vector of background frequencies, and a standard square integral substitution matrix such as BLOSUM-62.  The program calculates M[][], a square transition matrix which is used to calculate PSI-BLAST like scores.

   Next, any time you want to convert a vector of counts to a vector of scores, you call psiscore().  You need to provide the M[][] calculated before, as well as background frequencies, the amino acid counts, a pseudocount parameter, and a scale parameter.  Pseudocounts equal to 11 seem to work reasonably well.  A scale parameter 2 will yield scores in half-bits, 3 in third-bits, etc.
   Let me know whether these instructions are clear enough, and whether you have any difficulty getting the program to run.

-Stephen
 ******************************************************/

double	ssx_typ::ColScoreP2P(ssx_typ *that,Int4 i,Int4 j)
// return the P2P scoring object for doing an alignment.
{
	assert(nBlksCMSA(CMA) == 1);
	assert(i > 0 && i <= LengthCMSA(1,CMA));
	cma_typ cma=that->CMA; 
	assert(nBlksCMSA(cma) == 1);
	assert(j > 0 && j <= LengthCMSA(1,cma));
#if 0
	return DMS->sdotfs(that->WtCnts[1][j],this->WtCnts[1][i]);
#else
	double FG,BG,DD;
	FG=DMS->bild(that->WtCnts[1][j]);
	BG=DMS->bild(this->WtCnts[1][i]);
	UInt4   *WtCntsFB; NEW(WtCntsFB,nAlpha(AB)+3,UInt4);
	for(Int4 r=0; r <= nAlpha(AB); r++){
                WtCntsFB[r]=WtCnts[1][j][r] + WtCnts[1][i][r];
        } DD=DMS->bild(WtCntsFB); free(WtCntsFB);
	// return ((FG+BG) - DD);
	return (DD-(FG+BG));
#endif
}

/************************
this->mcalc(blosum62freq,AlphaR(AB),Mtrx);
 ***********************/
// double	**ssx_typ::mcalc(double *back, int **pam, double **M)
// double	**ssx_typ::mcalc(char **pam)
double	**ssx_typ::mcalc(Int4 **pam)
// Calculates the implied probability of co-occurence of pairs of residues.
 //	double	back[];		/*	Background frequencies		*/
 // 	int	pam[][NAA];	/*	Standard substituion matrix	*/
 //	double	M[][NAA];	/*	Transition matrix		*/
{
	Int4	NAA=nAlpha(AB);
 	int	i,j,k,low,high,range;
 	double	up,beta,ftemp,New,lambda,sum;
 	double	*p,**M; 
	NEWP(M,NAA+5,double);
	for (i=1;i <= NAA;++i) NEW(M[i],NAA+5,double);

 	for (low=high=0,i=1;i <= NAA;++i){
	    for (j=1;j <= NAA;++j) {
 		if ((k=pam[i][j])>high) high=k;
 		if (k<low) low=k;
	    }
 	}
 	range=high-low; NEW(p,range+5,double);
 	// for (i=0;i <= range;++i) p[i]=0;
 	for (i=1;i <= NAA;++i){
		for (j=1;j <= NAA;++j) p[pam[i][j]-low]+=BackGrnd[i]*BackGrnd[j]; 
	} up=0.5;
 	do {
 		up*=2;
 		beta=exp(up);
 		ftemp=exp(up*(low-1));
 		for (sum=i=0;i<=range;++i) sum+= p[i] * (ftemp*=beta);
 	} while(sum < 1.0);
 	for (lambda=j=0;j<25;++j) {
 		New=(lambda+up)/2.0;
 		beta=exp(New);
 		ftemp=exp(New*(low-1));
 		for (sum=i=0;i<=range;++i) sum+= p[i] * (ftemp*=beta);
 		if (sum>1.0) up=New; else lambda=New;
 	}
 	for (i=1;i <= NAA;++i){
		for (j=1;j <= NAA;++j){ M[i][j]=BackGrnd[j]*exp(lambda*pam[i][j]); } 
	}
#if 1
	fprintf(stderr,"done with mcalc\n");
	double D=0.0;
 	for (i=1;i <= NAA;++i) {
	   double total=0.0;
 	   for (j=1;j <= NAA;++j){
		// fprintf(stderr,"M[%c][%c]=%.4f\n",AlphaChar(i,AB),AlphaChar(j,AB),M[i][j]);
		total+=M[i][j];
	   } 
#if 0	// Don't normalize
	   fprintf(stderr," total = %.4f\n",total); D+=total;
#else	// normalize.
 	   for (j=1;j <= NAA;++j) M[i][j]=M[i][j]/total;
 	   for (total=0.0,j=1;j <= NAA;++j){
		// fprintf(stderr,"M[%c][%c]=%.4f\n",AlphaChar(i,AB),AlphaChar(j,AB),M[i][j]);
		total+=M[i][j];
	   } // fprintf(stderr," total = %.4f\n",total); D+=total;
	
#endif
	} // fprintf(stderr,"  grand total=%.4f\n",D);
	// exit(1);
#endif
	free(p);
	return M;
}

// void ssx_typ::psiscore(double *count,double *score)
void ssx_typ::psiscore(UInt4 *count,double *score)
//  	double	pseudo;		/*	Pseudocount parameter	*/
// 	double	scale;		/*	Scale parameter		*/
// 	double	count[];	/*	Count vector		*/
//  	double	score[];	/*	Calculated score vector	*/
//  	double	back[];		/*	Background frequencies	*/
// 	double	M[][NAA];	/*	Transition matrix	*/
{
	Int4	NAA=nAlpha(AB);
 	int	i,j;
 	double	sum;
 	double	scale;
 	double	f[NAA];
 	double	q[NAA];
 	double	alpha;

 	for (sum=0,i=1;i <= NAA;++i) sum+=(double)count[i];
 	for (i=1;i <= NAA;++i) f[i]=(double)count[i]/sum;

// assert(pernats==1000);

 	alpha=this->Pseudo/(this->Pseudo+sum-1);
 	// scale=this->Scale/log(2.0);  // Scale is already in nats.

 	for (i=1;i <= NAA;++i) {
 		for (q[i]=0,j=1;j <= NAA;++j) q[i]+=Matrix[j][i]*f[j];
 		q[i]=alpha*q[i]+(1-alpha)*f[i];
 		// score[i]=scale*log(q[i]/back[i]);
 		score[i]=pernats*log(q[i]/BackGrnd[i]);
 	} score[0]=0;
	// fprintf(stderr,"done with psiscore\n");
}

void    ssx_typ::calc_lngamma()
// static variables are guaranteed to be initialized to zero.
// May want to use this here??? (Instead of cmsa array).
// Need to recalculate this every time SeqWts() is called.
{
	Int4    time1=time(NULL),sq;
        UInt4   i;
        double  d;
	VrtlN=0;     // Virtual number of sequences * 100 after down weighting.
        // for(WtN=0,sq=1; sq <= N; sq++){ WtN += Wt[sq]; }
	VrtlN=WtN + TotalWtPseudo + N;	// Add on N to account for rounding up!!!
	VrtlN=2*VrtlN;	// double size to adjust for conversion errors.
	if(LGM != 0) delete LGM;
	// LGM = new lgm_typ(VrtlN,wt_factor);
	LGM = new lgm_typ(4*VrtlN,wt_factor);
#if 0	// Debug...
	fprintf(stderr,"WtN=%d; VrtlN=%d; %.2f Mbytes.\n",WtN,VrtlN,(double)(8*VrtlN)/1000000.0);
	fprintf(stderr,"\tcalc_lgamma time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
#endif
}

double	ssx_typ::FractionDeleted(Int4 blk, Int4 site)
{
	Int4	n,i,j,sq,nblks=nBlksCMSA(CMA);
	double	d,w,total;
	assert(blk > 0 && blk <= nBlksCMSA(CMA));
	assert(site > 0 && site <= LengthCMSA(blk,CMA));
	for(total=0,d=0,sq=1; sq <= NumSeqsCMSA(CMA); sq++){
		w=(double)Wt[sq];
                if(IsDeletedCMSA(blk,sq,site,CMA)) d += w;
		total += w;
	} return (d/total);
}

double	**ssx_typ::FractionDeleted( )
{
	Int4	n,blk,i,j,sq,nblks=nBlksCMSA(CMA);
	double	**D,d,w,total;
	NEWP(D,nblks +3, double);
	for(blk=1; blk <= nblks; blk++) {
	   NEW(D[blk],LengthCMSA(blk,CMA) +3, double);
	   for(i=1; i <= LengthCMSA(blk,CMA); i++){
	      for(total=0,d=0,sq=1; sq <= NumSeqsCMSA(CMA); sq++){
		w=(double)Wt[sq]/(double)wt_factor;
                if(IsDeletedCMSA(blk,sq,i,CMA)) d += w;
		total += w;
	      } D[blk][i]=d/total;
	   }
	} return D;
}

#define debug_vrtl_map 0

double  ssx_typ::VirtualMap(char mode,Int4 Blk,Int4 Col)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: This version doesn't require pseudocounts as these are already 
 included within VrtlCnts.
 **********************************************************************/
{
    Int4        n,len,i,j,t,T,s,ncol,off_col,on_col,nblks=nBlksCMSA(CMA);
    UInt4 	cnts,b,*nonsite,*site,total;
    double      D,d,MAP,weight;
    st_type	S=SitesCMSA(CMA);

    //=========== 1. Compute weight for number of blocks ==============
    for(weight=0.0, s=1; s<=N; s++){
	if(nBlksCMSA(CMA) == 1 && nSites(1,s,S) == 0) continue;
#if 1	// new for computing column contributions.
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
	  d=(double)Wt[s]/(double)wt_factor;
	  if(T > nblks) weight += d*lnbico(T, nblks);
	}
    } // weight=0;

    //============ 2. Compute log-probability for conserved blocks. ==============
    NEW(nonsite,nAlpha(AB) +3,UInt4);
    for(b=1;b<= nAlpha(AB); b++){
#if 0
	if(nonsite[b] >= TtlWtCnts[b]) nonsite[b] = WtPseudo;
	else nonsite[b] = TtlWtCnts[b] - nonsite[b];
        // fprintf(stderr,"non-site[%c] = %u/%u.\n",AlphaChar(b,AB),nonsite[b],TtlWtCnts[b]);
#else
	nonsite[b] = TtlWtCnts[b] + WtPseudo;
#endif
    }
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len=LengthCMSA(t,CMA);
#if 1	// new for computing column contributions.
	if(t == Blk){
	  if(Col==0) continue;
	  ncol = len; on_col += ncol;
	  for(D=0.0,i=1; i <= len; i++) {
	    if(i == Col){ on_col--; continue; }
	    site = VrtlCnts[t][i];
	    for(d=0.0,b=1;b<= nAlpha(AB); b++) {
		d += LGM->LGamma(site[b]); nonsite[b] -= site[b];
	    } D += d; D -= LGM->LGamma(TtlVrtlCnts[t][i]); 
#if debug_vrtl_map
	    UInt4 x=0,y;
	    for(b=1;b<= nAlpha(AB); b++){
		x += site[b];
		fprintf(stderr,"   LnGamma(site[%c] = %.3f) = %.3f\n",
			AlphaChar(b,AB),(double)(site[b])/(double)wt_factor,
			LGM->LGamma(site[b]));
	    }
	    double dd=LGM->LGamma(TtlVrtlCnts[t][i]);
	    fprintf(stderr,"column %d.%d: %.2f - %.3f = %.2f (%.2f =?= %.2f cnts)\n",
			t,i,d,dd,d-dd,(double)x/(double)wt_factor,
			(double)(y+TtlVrtlCnts[t][i])/(double)wt_factor);
#endif
	  } MAP+=D;
	} else 
#endif
	{
#if 0	// no fragmentation right now...use all positions!
          ncol = nColsFModel(model[t]); on_col += ncol;
	  off_col += LenFModel(model[t]) - ncol;
#else
	  ncol = len; on_col += ncol;
#endif
	  for(D=0.0,i=1; i <= len; i++) {
	    site = VrtlCnts[t][i];
	    for(b=1;b<= nAlpha(AB); b++) {
		D += LGM->LGamma(site[b]); nonsite[b] -= site[b];
	    } 
	    D -= LGM->LGamma(TtlVrtlCnts[t][i]);
#if debug_vrtl_map
	    UInt4 x=0,y=0;
	    for(d=0.0,b=1;b<= nAlpha(AB); b++){
		x += site[b];
		d += LGM->LGamma(site[b]); 
		fprintf(stderr,"   LnGamma(site[%c] = %.2f) = %.3f\n",
			AlphaChar(b,AB),(double)(site[b])/(double)wt_factor,
			LGM->LGamma(site[b]));
	    }
	    double dd=LGM->LGamma(TtlVrtlCnts[t][i]);
	    fprintf(stderr,"   LnGamma(site[all] = %.3f) = %.3f\n",
			(double)TtlVrtlCnts[t][i]/(double)wt_factor, dd);
	    fprintf(stderr,"column %d.%d: %.2f - %.3f = %.2f (%.2f =?= %.2f cnts) SumLLR=%.2f\n",
			t,i,d,dd,d-dd,(double)x/(double)wt_factor,
			(double)(y+TtlVrtlCnts[t][i])/(double)wt_factor,MAP+D);
#endif
	  } MAP+=D;
	}
        // fprintf(stderr,"  Running Map=%.2f (change=%.2f)\n",MAP,D);
        // if(debug_vrtl_map && t==4) fprintf(stderr,"  block 4 LLR = %.2f\n",D); 
    }
    if(debug_vrtl_map) fprintf(stderr,"Map=%.2f; weight = %.2f\n",MAP,weight);

    //============ 3. Compute log-probability for non-aligned regions. ============
#if 0
    if(CMA->FullSeq){	// check for problems due to flanking regions...
	for(b=1;b<= nAlpha(AB); b++) assert(nonsite[b] >= 0);
    }
#endif
    for(D=0.0,total=0, b=1;b<= nAlpha(AB); b++) {
	if(nonsite[b] > 0){ D += LGM->CalcLnGamma(nonsite[b]); total += nonsite[b]; }
        // fprintf(stderr,"Running non-site (%c:%u/%u) = %.2f.\n",AlphaChar(b,AB),nonsite[b],TtlWtCnts[b],D);
    } D -= LGM->CalcLnGamma(total);
    if(debug_vrtl_map) fprintf(stderr,"Non-conserved contribution (%.2f) + MAP (%.2f) = %.2f.\n",D,MAP,MAP+D);
    MAP += D; 

    //============ 4. Subtract Null Log-Proability =============
    for(D=0.0,total = 0, b=1;b<= nAlpha(AB); b++) {
	if(TtlWtCnts[b] > 0){ D -= LGM->CalcLnGamma(TtlWtCnts[b]); total += TtlWtCnts[b]; }
    } D += LGM->CalcLnGamma(total);
    if(debug_vrtl_map) fprintf(stderr,"Null Map=%.2f; Map=%.2f.\n",-D,MAP);
    MAP += D; 

    //============= 5. Add weight for column transfers =================
#if 1	// don't use because doesn't change? But might later if implement column sampling.
    if(nblks > 1){
	if(Blk > 0 && Col==0) weight +=  lnbico(on_col-nblks-2,nblks-2);
	else weight +=  lnbico(on_col-nblks-1,nblks-1);
	// if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }
#endif
    // D=JLH->IndelPenalty(0,mode,Wt);
    D=NDL->IndelPenalty(0,mode,Wt);
#if debug_vrtl_map
    fprintf(stderr,"LLR=%.2f (aln) - %.2f (blks) - %.2f (indel penalty) = %.2f; mode='%c'.\n",
		MAP,weight,-D,MAP-weight+D, mode);
#endif
    free(nonsite);
    return (MAP - weight + D);
}

double  ssx_typ::VirtualMap2(char mode)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: This version doesn't require pseudocounts as these are already 
 included within VrtlCnts.
 **********************************************************************/
{
    Int4        n,len,i,j,t,T,s,ncol,off_col,on_col,nblks=nBlksCMSA(CMA);
    UInt4 	cnts,b,*nonsite,*site,total,Ps=WtPseudo/2,tPs=TotalWtPseudo/2;
    double      D,d,MAP,weight;
    st_type	S=SitesCMSA(CMA);

    //=========== 1. Compute weight for number of blocks ==============
    for(weight=0.0, s=1; s<=N; s++){
	if(nBlksCMSA(CMA) == 1 && nSites(1,s,S) == 0) continue;
	T = LenSeqCMSA(s,CMA) + nblks;
	for(t=1; t <= nblks; t++) T -= LengthCMSA(t,CMA);
	d=(double)Wt[s]/(double)wt_factor;
	if(T > nblks) weight += d*lnbico(T, nblks);
    }

    //============ 2. Compute log-probability for conserved blocks. ==============
    NEW(nonsite,nAlpha(AB) +3,UInt4);
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len=LengthCMSA(t,CMA);
	ncol = len; on_col += ncol;
	for(D=0.0,i=1; i <= len; i++) {
	    site = VrtlCnts[t][i];
	    for(b=1;b<= nAlpha(AB); b++) {
		D += LGM->LGamma(site[b]+Ps); nonsite[b] += site[b];
	    } D -= LGM->LGamma(TtlVrtlCnts[t][i]+tPs); 
	} MAP+=D;
        // fprintf(stderr,"  Running Map=%.2f (change=%.2f)\n",MAP,D);
        // if(debug_vrtl_map && t==4) fprintf(stderr,"  block 4 LLR = %.2f\n",D); 
    }
    if(debug_vrtl_map) fprintf(stderr,"Map=%.2f; weight = %.2f\n",MAP,weight);

    //============ 3. Compute log-probability for unaligned regions. ============
    for(b=1;b<= nAlpha(AB); b++){
        // fprintf(stderr,"non-site[%c] = %u/%u.\n",AlphaChar(b,AB),nonsite[b],TtlWtCnts[b]);
	if(nonsite[b] >= TtlWtCnts[b]) nonsite[b] = 0;	// Ps added below...
	else nonsite[b] = TtlWtCnts[b] - nonsite[b];
    }
    for(D=0.0,total=0, b=1;b<= nAlpha(AB); b++) {
	if(nonsite[b] > 0){ D += LGM->CalcLnGamma(nonsite[b]+Ps); total += nonsite[b]; }
        // fprintf(stderr,"Running non-site (%c:%u/%u) = %.2f.\n",AlphaChar(b,AB),nonsite[b],TtlWtCnts[b],D);
    } D -= LGM->CalcLnGamma(total+tPs);
    if(debug_vrtl_map) fprintf(stderr,"Non-conserved contribution (%.2f) + MAP (%.2f) = %.2f.\n",D,MAP,MAP+D);
    MAP += D; 

    //============ 4. Subtract Null Log-Proability =============
    for(D=0.0,total = 0, b=1;b<= nAlpha(AB); b++) {
	if(TtlWtCnts[b] > 0){ D -= LGM->CalcLnGamma(TtlWtCnts[b]+Ps); total += TtlWtCnts[b]; }
    } D += LGM->CalcLnGamma(total+tPs);
    if(debug_vrtl_map) fprintf(stderr,"Null Map=%.2f; Map=%.2f.\n",-D,MAP);
    MAP += D; 

    //============= 5. Add weight for column transfers =================
#if 1	// don't use because doesn't change? But might later if implement column sampling.
    if(nblks > 1){
	weight +=  lnbico(on_col-nblks-1,nblks-1);
	// if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }
#endif
    // D=JLH->IndelPenalty(0,mode,Wt);
    D=NDL->IndelPenalty(0,mode,Wt);
#if debug_vrtl_map
    fprintf(stderr,"LLR=%.2f (aln) - %.2f (blks) - %.2f (indel penalty) = %.2f; mode='%c'.\n",
		MAP,weight,-D,MAP-weight+D, mode);
#endif
    free(nonsite);
    return (MAP - weight + D);
}
 
              
double  ssx_typ::RelMap(char mode)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 **********************************************************************/
{
    cma_typ	cma=this->CMA; // return UnGappedRelMapCMSA(cma);
    Int4        n,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(cma);
    fm_type	*model = ModelsCMSA(cma);    
    Int4	nblks=nBlksCMSA(cma);
    a_type	A;
    double	indel_penalty=0.0;
    ss_type	truedata = TrueDataCMSA(cma);

    if(cma->FullSeq){	// then use full_counts and freq;
	counts = CntsSeqSet(cma->FullSeq); freq=tFreqSeqSet(cma->FullSeq);
    } else { // OLD: prior to full_counts.
       counts = CntsSeqSet(truedata); freq=tFreqSeqSet(truedata); 
    }
    for(len=0, t=1; t <= nblks; t++) {
	n = MaxLenFModel(model[t]); if(n > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    /******** compute weight for number of blocks ********/
    st_type	S=SitesCMSA(cma);
    if(cma->FullSeq){ 	// NEW: full_counts.
      for(weight=0.0, s=1; s<=NSeqsSeqSet(cma->FullSeq); s++){
	if(cma->FullRpts[s] > 0){
	  n = cma->FullRpts[s]*nBlksCMSA(cma);
	  T = SqLenSeqSet(s,cma->FullSeq) + n;
	  T -= n*TotalLenCMSA(cma);
	  if(T < n) T = n;  // allows for gaps...
	  weight += lnbico(T,n);  // number of positions...
	}
      }
    } else {
      for(weight=0.0, s=1; s<=N; s++){
	if(nBlksCMSA(cma) ==1 && nSites(1,s,S) == 0) continue;
	T = SqLenSeqSet(s,data) + nblks;
	for(t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	if(T < nblks) T = nblks;  // allows for gaps...
	weight += lnbico(T, nblks);
      }
    }
    /******** end compute weight for number of blocks ********/
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(AB);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	/*************** NEW: column width weight ******************/
	MAP -= lnbico(len - 2, ncol - 2);
	/***********************************************************/
	/***** MAP *******/
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(AB); b++) {
		MAP += cma->lngamma[site[b]][b]; 
		nonsite[b] -= site[b];
	    }
	    MAP -= lgamma((double)(totsite-site[0]) + npseudo);
	    // MAP -= lgamma((double)totsite + npseudo);
	  } 
	}
    }
    v = lgamma(npseudo);

    /***** NONSITES MAP *******/
    if(cma->FullSeq){	// check for problems due to flanking regions...
	for(b=1;b<= nAlpha(AB); b++) assert(nonsite[b] >= 0);
    }
    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
	   if(Ps[b] > 0){
		MAP += lgamma((double) nonsite[b] + Ps[b]);
		total += (double) nonsite[b];
        	v -= lgamma(Ps[b]);   /** term pulling down MAP **/
	   }
    }
    MAP += v *(double) on_col; 
    MAP -= lgamma((double) total + npseudo);

    /***** subtract NULL MAP ******/
    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
	   if(counts[b] > 0.0){
		MAP -= lgamma(counts[b] + Ps[b]);
		total += counts[b];
	   }
    }
    MAP += lgamma((double)total + npseudo);

    /**************** weight for column transfers *************/
    if(nblks > 1){
	weight +=  lnbico(on_col-nblks-1,nblks-1);
	if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }

    /**************** return MAP *************/
    free(observed);
#if 0
    indel_penalty=IndelPenaltyCMA(0,'S');
    return (MAP - weight + indel_penalty);
#elif 0
    indel_penalty=IndelPenaltySeqSet(data);
    return (MAP - indel_penalty);
#else
    // indel_penalty = JLH->IndelPenalty(0,mode,0);
    indel_penalty = NDL->IndelPenalty(0,mode,0);
    // fprintf(stderr,"Map=%.2f; weight = %.2f\n",MAP,weight);
    // fprintf(stderr," Final Map=%.2f\n",MAP-weight);
    return (MAP - weight + indel_penalty);
#endif
}

double  ssx_typ::CoreBildLLR(BooLean add_indels,char mode)
/**********************************************************************
 return the Bild Log-likelihood ratio for the aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 **********************************************************************/
{
    cma_typ	cma=this->CMA; 
    Int4        n,len,i,j,t,T,s,b,nblks=nBlksCMSA(cma);
    double      MAP,weight,D,d,indel_penalty=0.0;
    ss_type	data = DataCMSA(cma);

    if(this->ReCalcSMX) this->CalcSMX();
    /******** compute weight for number of blocks ********/
#if 0	// now modeled by indel_penalty
    for(weight=0.0, s=1; s<=N; s++){
	if(nBlksCMSA(cma) ==1 && nSites(1,s,SitesCMSA(cma)) == 0) continue;
	T = SqLenSeqSet(s,data) + nblks;
	for(t=1; t <= nblks; t++) T -= LengthCMSA(t,cma);
	if(T > nblks){
	   d=(double)Wt[s]/(double)wt_factor;
	   weight += d*lnbico(T, nblks);
	}
    }
#endif
    /******** end compute weight for number of blocks ********/

    for(MAP=0.0,t=1; t <= nblks; t++) {
	len = LengthCMSA(t,cma); 
	for(i=1; i <= len; i++) { MAP += DMS->bild(WtCnts[t][i]); }
    }

    /**************** return MAP *************/
    if(add_indels) indel_penalty = NDL->IndelPenalty(0,mode,Wt); else indel_penalty=0;
    // fprintf(stderr,"Map=%.2f; weight = %.2f; indel_penalty=%.2f\n",MAP,-weight,indel_penalty);
    // fprintf(stderr," Final Map=%.2f\n",MAP-weight);
    // return (MAP - weight + indel_penalty);
    return (MAP + indel_penalty);
}

#include "swaln.h" 

char	*ssx_typ::AlignP2P(ssx_typ *that,Int4 &start,Int4 &Score)
// return the P2P scoring object for doing an alignment.
{
	Int4	i,j,I,J,b,B,r,IntScore;
	double	scrI,scrJ,score;
	cma_typ cma=that->CMA;
	UInt4	*wtCntsI,**wtCntsJ,*wtCntsIJ; NEW(wtCntsIJ, nAlpha(AB) +3, UInt4);
	this->InitNDL(0);
	BooLean	debug=0;

	h_type	HG=0;
	e_type  CsqE=0,csqE=0;
	unsigned char *Csq=0,*csq=0;
if(debug){
	// HG=Histogram("BILD scores",-5000,5000,0.2);
	HG=Histogram("BILD scores",-5000,5000,50);
	CsqE=MkConsensusCMSA(CMA);
	csqE=MkConsensusCMSA(cma);
	Csq=SeqPtr(CsqE);
	csq=SeqPtr(csqE);
}
	Int4 D=that->WtN/(UInt4)wt_factor;
#if 0	// Bild profile-to-profile scoring.
	p2p_typ *p2p = new p2p_typ(TotalLenCMSA(cma), nBlksCMSA(CMA),LengthsCMSA(CMA),D);
	for(I=0,b=1; b <= nBlksCMSA(cma); b++){
	   for(i=1; i <= LengthCMSA(b,cma); i++){
	     I++; wtCntsI=that->WtCnts[b][i];
	     scrI=DMS->bild(wtCntsI);
	     for(J=0,B=1; B <= nBlksCMSA(CMA); B++){
	       wtCntsJ=this->WtCnts[B];
	       for(j=1; j <= LengthCMSA(B,CMA); j++){
	         J++; scrJ=DMS->bild(wtCntsJ[j]);
		 for(r=1; r <= nAlpha(AB); r++) wtCntsIJ[r]= wtCntsI[r] + wtCntsJ[j][r];
	         score=DMS->bild(wtCntsIJ);
// fprintf(stderr,"score[%d][%d][%d] = %g (%g; %g)\n",I,B,j,score,scrI,scrJ);
		 score = score - scrI - scrJ;
#if 0	// test...
		if(I == J) p2p->SetScore(I,B,j,1000);
		else p2p->SetScore(I,B,j,-1000);
#elif 0	// See if alignment is reasonable for 
IntScore=random_integer(2000) - 1010;
		 p2p->SetScore(I,B,j,IntScore);
#elif 0		// pairwise scores for the two sequences...
		 double dd=0.693*(double)valAlphaR(Csq[J],csq[I],AB)/2.0;  // half bits to nats.
		 if(HG) IncdHist(dd, HG);
		 IntScore = (Int4) floor(0.5+(double)pernats*dd);
		 p2p->SetScore(I,B,j,IntScore);
#else		// BILD scoring...
		 if(HG) IncdHist(score, HG);
		 IntScore = (Int4) floor(0.5+((double)pernats*score));
		 p2p->SetScore(I,B,j,IntScore);
#endif
	       }
	     }
	   }
	} free(wtCntsIJ);
#else	// alternative profile-to-profile scoring method...
	assert(nBlksCMSA(CMA) == 1);
	assert(nBlksCMSA(cma) == 1);
	p2p_typ *p2p = new p2p_typ(TotalLenCMSA(cma), nBlksCMSA(CMA),LengthsCMSA(CMA),D);
	for(I=0,b=1; b <= nBlksCMSA(cma); b++){
	   for(i=1; i <= LengthCMSA(b,cma); i++){
	     I++; wtCntsI=that->WtCnts[b][i];
	     // scrI=DMS->bild(wtCntsI);
	     for(J=0,B=1; B <= nBlksCMSA(CMA); B++){
	       wtCntsJ=this->WtCnts[B];
	       for(j=1; j <= LengthCMSA(B,CMA); j++){
	         J++; 
		 score=DMS->sdotfs(wtCntsJ[j],wtCntsI);
		 if(HG) IncdHist(score, HG);
		 IntScore = (Int4) floor(0.5+((double)pernats*score));
		 p2p->SetScore(I,B,j,IntScore);
	       }
	     }
	   }
	} free(wtCntsIJ);
#endif
if(HG){
   PutHist(stderr,60,HG); NilHist(HG);
   AlnSeqSW(11, 1, CsqE,csqE,AB);
   NilSeq(CsqE); NilSeq(csqE);
}
	// char *operation=JLH->GapAlnTraceP2P(p2p,start,Score);
	char *operation=NDL->GapAlnTraceP2P(p2p,start,Score);
fprintf(stderr,"Score = %d; start = %d\n",Score,start);
	delete p2p;
	return operation;
}

double  ssx_typ::RelEntropy(Int4 blk, Int4 col)
{
	UInt4 *wtCnt=WtCnts[blk][col];
	Int4 r;
        double RE,D,d,p=0,q=0,*BG=DMS->BackGrnd();
        for(D=d=0.0,r=1; r <= nAlpha(AB); r++){
                        D += (double) wtCnt[r]/(double) wt_factor;
                        d += BG[r];
        }
        for(RE=0,r=1; r <= nAlpha(AB); r++){
                        p = ((double) wtCnt[r]/(double) wt_factor)/D;
                        q = BG[r]/d;
                        if(p > 0) RE += p*log(p/q);
// fprintf(stderr,"%d: p=%g; q=%g; RE=%g\n",r,p,q,RE);
        } return RE;
}

double  ssx_typ::Entropy(Int4 blk, Int4 col)
{
	UInt4 *wtCnt=WtCnts[blk][col];
	Int4 r;
        double Ent,D,p=0;
        for(D=0.0,r=1; r <= nAlpha(AB); r++){ D += (double) wtCnt[r]/(double) wt_factor; }
        for(Ent=0,r=1; r <= nAlpha(AB); r++){
                        p = ((double) wtCnt[r]/(double) wt_factor)/D;
                        if(p > 0) Ent += p*log(p);
        } return Ent;
}

#if 0
   The program should calculate transition probabilities, so that M[i][j] is the probability of amino acid mutating to amino acid j.
Thus sum_j M[i][j] should indeed be 1.  I will describe below the reason I think this is not working out properly.  (A bug is also
possible.)

   The basic idea is as follows.  Given a substitution score matrix pam[i][j] of arbitrary scale, and a set of background frequencies back[j], the program first calculates the implicit scale lambda for the matrix.  Given this lambda, one should have:

exp(lambda*pam[i][j]) = q[i,j]/(p[i]*p[j]),

where q[i][j] is the probability of observing amino acids i and j aligned in related sequences.  Now, q[i,j] can be written either as p[i]*P(i->j) or as p[j]*P(j->i), where P(i->j) means the probability that amino acid i mutates to amino acid j.  (This is a simplifying assumption that Dayhoff makes that may not in fact be true.  In other words, she assumes that the number of times amino acid i mutates to amino acid j is balanced by the number of times amino acid j mutates to amino acid i.  Because the background frequencies are different, one must have asymmetrical mutation probabilities.)  Thus, in my program,

M[i][j]=back[j]*exp(lambda*pam[i][j])

can be interpreted as:

M[i][j] = p[j]*q[i,j]/(p[i]*p[j]) = p[j]*(p[i]*P[i->j])/(p[i]*p[j]) = P[i->j].

I think the problem is that the development above assumes consistency between the target and background frequencies, i.e. that the q[i,j] implicit in the substitution scores have the property:
sum_j q[i,j] = p[i].  To guarantee this is the case, one needs either to derive the p[i] from the substitution scores, or compositionally adjust the substitution scores to be consistent with the specified background frequencies; see [1].  Alternatively, if the scores are explicitly constructed to be consistent (like PAM or BLOSUM scores), one needs to use the background frequencies that these constructions used (see below).  Even this will not work perfectly, because of rounding to construct integral substitution scores.  One can also always take the calculated M[i][j], and normalize them so that sum_j M[i][j] = 1.
   I believe the PAM background frequencies are approximately:

0.08713, 0.04090, 0.04043, 0.04687, 0.03347, 0.03826, 0.04953, 0.08861, 0.03362, 0.03689, 0.08536, 0.08048, 0.01475, 0.03977, 0.05068, 0.06958, 0.05854, 0.01049, 0.02992, 0.06472

and the BLOSUM-62 background frequencies are approximately:

0.07422, 0.05161, 0.04465, 0.05363, 0.02469, 0.03426, 0.05431, 0.07415, 0.02621, 0.06792, 0.09891, 0.05816, 0.02499, 0.04742, 0.03854, 0.05723, 0.05089, 0.01303, 0.03228, 0.07292

where the amino acid order is:   ARNDCQEGHILKMFPSTWYV

[1] Yu, Y.-K., Wootton, J.C. & Altschul, S.F. (2003) "The compositional adjustment of amino acid substitution matrices,"
Proc. Natl. Acad. Sci. USA 100:15688-15693.

#endif

