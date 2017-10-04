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

#include "bpps_typ.h"

void bpps_typ::BinaryUpdateTheta_J(Int4 j, UInt4 **CntBG,UInt4 **CntFG)
/************************************************************************
  Quick recompute for jth position only.
 ************************************************************************/
{
        unsigned char   x,q_j;

	assert(j > 0 && j <= k);
        // for(Int4 j=1;j<=k;j++) 
	{
                assert(qst[j] != 0);
                q_j = query[j];
                if(q_j == 0) print_error("FATAL: X residues disallowed in query\n");
                double mat=PseudoCnts/wt_factor;                // add prior == 100.0
                double sum=PseudoCnts/wt_factor;                // add prior == 100.0
#if 1   // treat deletions like random background...
           if(del_as_random){
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(mat,mat,sum,sum,qst[j],CntBG[j][0],CntFG[j][0]);
           }
#endif
                for(x=startAB;x<=nAlpha(AB);x++) {
                   if(MemSset(x,qst[j])){ mat += CntBG[j][x] + CntFG[j][x]; }
                   else { sum += CntBG[j][x] + CntFG[j][x]; }
                }
                sum += mat;     // Failed to sum this in previous loop.
                assert(mat > SumXI[j]);
                mat -= SumXI[j];        // subtract out FG non-contamination, which doesn't count.
                sum -= SumXI[j];        // subtract out FG non-contamination, which doesn't count.
                theta[j][q_j] = mat/sum;	// not theta = 1 - theta for binomial model..
#if 0
                if(DebugCheck && j == 69){ fprintf(stderr,"theta[%d][%c] = %f; mat = %f; sum = %.2f\n",
                        j,AlphaChar(q_j,AB),theta[j][q_j],mat,sum); }
#endif
        }
}

void bpps_typ::BinaryUpdateSumXI_J(Int4 j, UInt4 **CountFG,UInt4 **CntBG)
/************************************************************************
  Quick 
 ************************************************************************/
{
        unsigned char   x;

	assert(j > 0 && j <= k);
        assert(qst[j] != 0); 
        alpha=a0/(a0+b0);
	{
           double matFG=0.0,matBG=0.0,sumBG=0.0;
           matBG=PseudoCnts/wt_factor;          // add prior == 100.0; m1
           sumBG=2.0*PseudoCnts/wt_factor;      // add prior *2; n2 + m2.
#if 0   // treat deletions like random background...
           if(del_as_random){
                dummy=Dummy=0.0;
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(matFG,matBG,dummy,Dummy,qst[j],CntBG[j][0],CountFG[j][0]);
                sumBG += CntBG[j][0];
           }
#endif
           for(x=startAB; x <= nAlpha(AB); x++){
                if(MemSset(x,qst[j])){ matFG += CountFG[j][x]; matBG += CntBG[j][x]; }
                sumBG += CntBG[j][x];
           } // theta_x is the fraction of background residues matching the pattern.
           double theta_x=matBG/sumBG;
           SumXI[j] = matFG*(alpha/(alpha+((1-alpha)*theta_x))); // from Equation [3].
           if(SumXI[j] != 0 && matFG != 0) assert(SumXI[j] < matFG);
	   // updateAlpha(CountFG);  // is a function of qst & SumXI[j]! (defined based on all j)
        }
}

double	bpps_typ::subNullMapBiNom_J(Int4 j, UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
// Considers the query residue set at each position to be null.
{
	Int4		t;
	UInt4	*CountBG_j,*CountFG_j;

        assert(j > 0 && j <= k);
   	assert(qst[j] != 0);
	unsigned char	x,q_j = query[j];
	double sum=PseudoCnts*wt_factor;
	double match=PseudoCnts*wt_factor;
#if 1   // treat deletions like random background...
           if(del_as_random){
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(match,match,sum,sum,qst[j],CntBG[j][0],CntFG[j][0]);
           }
#endif
        for(x=startAB;x<=nAlpha(AB);x++) {
	   if(MemSset(x,qst[j])) match += CntBG[j][x] + CntFG[j][x];
	   else sum += CntBG[j][x] + CntFG[j][x];
        } sum += match;
	NullTheta[j] = match/sum;

	q_j=query[j];
	double theta_x = NullTheta[j];
	CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
	double m1=0.0,m2=0.0;
	Int4 CardPttrnSet=0;
#if 1   // treat deletions like random background...
           if(del_as_random){
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(m1,m1,m2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
           }
#endif
	for(t=startAB; t<= nAlpha(AB); t++) {
	   if(MemSset(t,qst[j])){ CardPttrnSet++; m1+= CountBG_j[t] + CountFG_j[t]; }
	    else { m2+= CountBG_j[t] + CountFG_j[t]; }
	}
	double submap=0.0;
	m1 = m1*wt_factor; m2 = m2*wt_factor;
	submap += m1*log(theta_x) + m2*log((1.0-theta_x));
	assert(CardPttrnSet > 0);
	submap += logRho[j][0];	// log of one minus other Rho values..

	return submap;
}

void bpps_typ::FastUpdateAlpha(UInt4 **CountFG)
{
        register Int4   j;
        register char  	x;
        register double	N0=0.0,S=0.0;
	register UInt4	*CntFG_j;

        for(j=k; j > 0; j--){
           if(qst[j]){
		CntFG_j=CountFG[j];
		for(x=nAlpha(AB); x >= startAB;  x--){ if(MemSset(x,qst[j])) N0 += CntFG_j[x]; }
		S += SumXI[j];
	   }
        } alpha = (S+a0)/(a0+N0+b0);  // mean of beta = (S+a0)/(S + a0 + N0 - S + b0).
        assert(alpha <= 1.0 && alpha >= 0.0);
}

void	bpps_typ::UpdateNullBiNomTheta(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
   Considers the query residue set at each position to be null.
 ************************************************************************/
{
	Int4		j;
        for(j=1;j<=k;j++) {
                if(!qst[j]) continue;
                unsigned char x;
                double sum=PseudoCnts*wt_factor;
                double match=PseudoCnts*wt_factor;
#if 1   // treat deletions like random background...
           if(del_as_random){
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(match,match,sum,sum,qst[j],CntBG[j][0],CntFG[j][0]);
           }
#endif
                for(x=startAB;x<=nAlpha(AB);x++) {
                   if(MemSset(x,qst[j])) match += CntBG[j][x] + CntFG[j][x];
                   else sum += CntBG[j][x] + CntFG[j][x];
                } sum += match;
                NullTheta[j] = match/sum;
        }
}

void	bpps_typ::UpdateParameters(UInt4 **CntBG, UInt4 **CntFG)
{
	BinaryUpdatePseudoCnts(CntFG,CntBG);	// not a function of qst.
	alpha=a0/(a0+b0);
        BinaryUpdateSumXI(CntFG,CntBG);		// is a function of qst & alpha prior.
        FastUpdateAlpha(CntFG);			// is a function of qst & SumXI[j]! (defined based on all j)
	BinaryUpdateTheta(CntBG,CntFG);		// is a function of qst! (theta[j] update only)?
	UpdateNullBiNomTheta(CntBG,CntFG);
}

#if 1
double	bpps_typ::SeqLikelihoodFG(unsigned char *seq, UInt4 sqwt)
// Computes the effect on the LPR of adding seq to the FG partition.
{
	Int4		j;
	double		like=0.0,nlike=0.0;
	unsigned char	r,c;
	Int4    NumPttrnPos=0;          // NEW for PriorRi
	double	AveNumWtSeqFG=0.0;
	double	SqWt=sqwt*wt_factor; 
	for(j=k; j > 0; j--) {
           if(qst[j] == 0) continue;
           if(seq[j] == 0) continue;
	   else NumPttrnPos++;     
	   for(c=0,r=nAlpha(AB); r >= startAB; r--) if(MemSset(r,qst[j])) c++;
	   assert(c > 0 && Rho[j][c] > 0.0);
	   if(MemSset(seq[j],qst[j])){
                like += SqWt*log((1-alpha)*theta[j][query[j]] + alpha); 
	   	nlike += SqWt*log(NullTheta[j]); // this is for the NulModel.
	   } else {
                like += SqWt*log((1-alpha)*(1.0-theta[j][query[j]]));   // delta_j == 0.
	   	nlike += SqWt*log(1.0 - NullTheta[j]);		// this is for the NulModel.
	   }
	   like += logRho[j][c];
	   nlike += logRho[j][0];	// for Null model..
	   AveNumWtSeqFG +=SqWt;
	} 
	if(NumPttrnPos > 0 && LogRatioPriorRi != 0.0){
		AveNumWtSeqFG = AveNumWtSeqFG/NumPttrnPos;
		like += AveNumWtSeqFG * LogRatioPriorRi;	// 
		if(isnan(like)){
                        fprintf(stderr,"AveNumWtSeqFG=%g; NumPttrnPos=%g; LogRatioPriorRi=%g\n",
                                AveNumWtSeqFG,NumPttrnPos,LogRatioPriorRi);
                        abort();
                }
	} 
	return like;
	return like - nlike;
	return nlike;
}

double	bpps_typ::SeqLikelihoodBG(unsigned char *seq,UInt4 sqwt)
// Computes the effect on the LPR of adding seq to the BG partition.
//*************** (1-Rj) terms: *************** 
{
	Int4		j;
	double		like=0.0,nlike=0.0;
	unsigned char	r,c;
	double	SqWt=sqwt*wt_factor; 
	for(j=k; j > 0; j--) {
           if(qst[j] == 0) continue;
           if(seq[j] == 0) continue;
	   for(c=0,r=nAlpha(AB); r >= startAB; r--) if(MemSset(r,qst[j])) c++;
	   assert(c > 0 && Rho[j][c] > 0.0);
	   if(MemSset(seq[j],qst[j])){
		like += SqWt*log(theta[j][query[j]]);	
	   	nlike += SqWt*log(NullTheta[j]);		// this is for the NulModel.
	   } else {
                like += SqWt*log((1.0-theta[j][query[j]]));
	   	nlike += SqWt*log(1.0 - NullTheta[j]);		// this is for the NulModel.
	   }
	   like += logRho[j][c];
	   nlike += logRho[j][0];	// for Null model..
	} 
	// like += SqWt * log(1.0-PriorRi);
	// like += SqWt*LogRatioPriorRi;	//  not done for BG residues; see why?
	return like;
	return like - nlike;
	return nlike;
}
#endif

sst_typ	*bpps_typ::RtnOptPattern(UInt4 **CntBG, UInt4 **CntFG, sst_typ **LegalSST, unsigned char *seq,Int4 MaxCols,double &Lpr)
/************************************************************************
 ************************************************************************/
{
	const double NEG_INFINITY=-1.0e+100;
	double	submap,*theta_j,d_max=0.0;
	sst_typ	*opt_sst,best_sst=0;
	Int4		j,r,r_max;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	dh_type         dH=dheap(k+2,4);

	BinaryUpdatePseudoCnts(CntFG,CntBG);	// not a function of qst.
	alpha=a0/(a0+b0);
        BinaryUpdateSumXI(CntFG,CntBG);		// is a function of qst & alpha prior.
        FastUpdateAlpha(CntFG);			// is a function of qst & SumXI[j]! (defined based on all j)
	BinaryUpdateTheta(CntBG,CntFG);		// is a function of qst! (theta[j] update only)?
	for(j=1; j<=k; j++) {

	   if(seq[j]==0) continue;
	   sst_typ *sstJ=LegalSST[seq[j]];  
	   for(best_sst=0,d_max=NEG_INFINITY,r_max=0,r=1; sstJ[r] != 0; r++){    // check all legal patterns.
		this->AddColumn(j,sstJ[r]);  // add sstJ to position s; sets qst[j]=sstJ[r].

		BinaryUpdateSumXI_J(j,CntFG,CntBG);	// is a function of qst & alpha prior = a0/(a0+b0).
		FastUpdateAlpha(CntFG);			// is a function of qst & SumXI[j]! (defined based on all j)
		BinaryUpdateTheta_J(j,CntBG,CntFG);	// is a function of qst! (theta[j] update only)?
		
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
		theta_j = theta[j]; q_j=query[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		double m1,n1,m2,n2; m1=n1=m2=n2=0.0; 	
		Int4	CardPttrnSet=0;
#if 1   // treat deletions like random background...
           if(del_as_random){
                // CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
                CntDelete(n1,m1,n2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
           }
#endif
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			CardPttrnSet++; n1+= CountFG_j[t]; m1+=CountBG_j[t];
		   } else { n2+= CountFG_j[t]; m2+=CountBG_j[t]; }
		} submap=0.0;
		m1=m1*wt_factor; m2=m2*wt_factor;
		n1=n1*wt_factor; n2=n2*wt_factor;
		MatchFG[j] = n1; MisMatchFG[j] = n2;
		MatchBG[j] = m1; MisMatchBG[j] = m2;
		//*************** Rj terms: ***************
		submap += n1*log((1-alpha)*theta_j[q_j] + alpha); // delta_j == 1.
		submap += n2*log((1-alpha)*(1.0-theta_j[q_j]));	  // delta_j == 0.
		//*************** (1-Rj) terms: *************** 
		submap += m1*log(theta_j[q_j]);	
		submap += m2*log((1.0-theta_j[q_j]));
		assert(CardPttrnSet > 0);
		submap += logRho[j][CardPttrnSet];
		submap = submap - subNullMapBiNom_J(j,CntBG,CntFG);
                if(submap > 0.0 && submap > d_max){ d_max=submap; best_sst=sstJ[r]; }
	   } // end of legal sst loop.
           if(best_sst != 0){
		// fprintf(stderr,"d_max(%d) = %g\n",j,d_max);	// debug...
                this->AddColumn(j,best_sst);  // add sstJ to position j;
	   	BinaryUpdateSumXI_J(j,CntFG,CntBG);		// need to update SumXI[j] to match best.
		// No need to update BinaryUpdateTheta_J() because it isn't used (thought it likely is incorrect).
                insrtHeap(j,(keytyp)d_max,dH);  // save the best ...
           } else this->RmColumn(j);	// no need to update SumXI[j]; 
	} // end of k loop.
	while(ItemsInHeap(dH) > MaxCols){ j=delminHeap(dH); this->RmColumn(j);}

        // 5.Trim back patterns until optimum lpr is reached;
        double lpr=0,lpr_max=this->LPR(CntBG,CntFG);           // CntBG & CntFG initialized by InitCntsFG_BG( );
        sst_typ sst=0;
        do {
           if(ItemsInHeap(dH) == 0) break;
           j=delminHeap(dH); sst=qst[j]; this->RmColumn(j);
           lpr=this->LPR(CntBG,CntFG);           // set lpr here...
           if(lpr < lpr_max) { lpr=lpr_max; this->AddColumn(j,sst); break; } else lpr_max=lpr;
        } while(TRUE); Nildheap(dH);
        // if(fp) this->PutSubLPR(fp,0,CntBG,CntFG);
	NEW(opt_sst,k+3,sst_typ);
	for(j=1; j<=k; j++) { opt_sst[j]=qst[j]; }
	Lpr=lpr;
	return opt_sst;
}

