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

#include "pps_typ.h"
#include "blosum62.h"

#define	PPS_TYP_FIX_GAP_PROBLEM 1

double	pps_typ::GetGapsMatch(double CntX, sst_typ qst)
// return the total number of number of Pseudo-Matches for CntX
{
	register double rtn=0.0;
	register unsigned char x;
        for(x=1; x <= nAlpha(AB); x++){
	   if(MemSset(x,qst)) rtn += CntX*blosum62freq[x];
	} return rtn;
}

pps_typ::pps_typ(sst_typ *Qst,unsigned char *Query, Int4 Query_len, double A0, 
	double B0, double *B, a_type A, UInt4 WtFactor, char Model, double Sigma,
	UInt4 *card_set,double Rho,char alt_measure,BooLean UseGlobalSqWts)
{
	if(UseGlobalSqWts) StartAlpha=0; else StartAlpha=1;
	Init(Qst,Query,Query_len,A0,B0,B,A,WtFactor,Model,Sigma,card_set,Rho,alt_measure); }

void    pps_typ::Init(sst_typ *Qst,unsigned char *Query, Int4 Query_len,
		double A0, double B0, double *B, a_type A,UInt4 WtFactor,
		char Model,double Sigma,UInt4 *card_set,double Rho,char alt_measure)
{
	DebugCheck=FALSE;
	CardSet=card_set;	// pointer to array that gets changed from calling environment.
	query = Query; qst=Qst;
	AltMeasure= alt_measure;
	a0 = A0*WtFactor;	// = alpha prior.
	b0 = B0*WtFactor;	// = beta prior.
	alpha = a0/(a0+b0);
	if(Rho <= 0.0){
		print_error("pps_typ ( ): LnRho out of bounds");
	} else if(Rho <= 20.0){
		logRho=-Rho;		// penalty to add a column in nats 
		logOneMinusRho = log(1.0-exp(logRho));
	} else {
		logRho=-Rho;		// penalty to add a column in nats 
		logOneMinusRho = 0.0;
	}
	assert(Sigma <= 1.0 && Sigma >= 0.50);
	sigma=Sigma;	// stringency of constraints...
	// sigma=0.70;	// stringency of constraints...
	SetModel(Model);
	b = B; 
	AB = A;
#if 0	// get stdfreq
	double sum=0.0;
	Int4 x;
	for(x=0; x <= nAlpha(AB); x++){
		stdfreq[x] = blosum62freq[x]; sum+=stdfreq[x];
	}
	for(x=0; x <= nAlpha(AB); x++) stdfreq[x] = stdfreq[x]/sum;
#endif
	wt_factor=1.0/(double)WtFactor;
	k = Query_len;
	NEWP(theta,k+1,double);
	NEW(NullTheta,k+1,double);
	NEW(SumXI,k+1,double);
	NumCols=0;
	for(Int4 j=1;j<=k;j++){ 
		NEW(theta[j],nAlpha(AB)+2,double); 
		if(qst[j]){
		  NumCols++;
		  assert(MemSset(query[j],qst[j]));
		}
	}
        NEW(argum,nAlpha(AB)+3,double);
	NEW(SubNullMap,k+3,double); NEW(SubMap,k+3,double);
	NEW(subLPR,k+3,double);
	NEW(MatchFG,k+3,double); NEW(MisMatchFG,k+3,double);
	NEW(MatchBG,k+3,double); NEW(MisMatchBG,k+3,double);
	nullMap=(double)INT4_MIN;
	ComputeNull=TRUE;
}

double	pps_typ::ProbCols( )
// Crude estimate of the probability of the chance of observing the current
// number of  columns shared by two families.
{
	double	p=0.0;
	double	ncols=0.0;
	for(Int4 j=1;j<=k;j++){ 
		if(qst[j]){
			ncols += 1.0;
			p += blosum62freq[query[j]];
		}
	} p = p/ncols;	// average p...
	double P = CumBinomProb(ncols, (double) k, p);
	// fprintf(stdout,"CumBinomProb(%.0f,%d,%.3f) = %g\n",ncols,k,p,P);
	return P;
}

void pps_typ::Free()
{
	for(Int4 j = 1; j <= k; j++) free(theta[j]); 
	free(theta); free(NullTheta); free(SumXI); free(argum);
	free(SubNullMap); free(SubMap); free(subLPR);
	free(MatchFG); free(MisMatchFG); free(MatchBG); free(MisMatchBG);
}

//************************** Binary ************************************
void pps_typ::BinaryUpdateSumXI(UInt4 **CountFG,UInt4 **CntBG)
/************************************************************************
 Conditional on R, alpha and Theta, update the |Xi(.j)| for each column
 according to Equation 3.
 Greek letter Xi(i,j) is an indicator variable for position j of 
	sequence i belonging to the query cluster, where 
	Xi(i,j) == 1 --> x = seq[i][j] comes from foreground component.
	and Xi(i,j) == 0, the background component.
 The theta_x is the fraction of background residues matching the pattern.
 ************************************************************************/
{
        unsigned char	x;

        for(Int4 j=1;j<=k;j++){
	   if(!qst[j]) continue;	// null set here - ignore this position.
	   double matFG=0.0,matBG=0.0,sumBG=0.0;
	   // matFG=1.0/wt_factor;		// add prior == 100.0; n1
	   matBG=1.0/wt_factor;		// add prior == 100.0; m1
	   sumBG=2.0/wt_factor;		// add prior *2; n2 + m2.
	   for(x=1; x <= nAlpha(AB); x++){
		if(MemSset(x,qst[j])){
			matFG += CountFG[j][x];
			matBG += CntBG[j][x];
		} 
		sumBG += CntBG[j][x];
	   } 
	   // theta_x is the fraction of background residues matching the pattern.
	   double theta_x=matBG/sumBG;	
	   // could sample SumXI[j] from the binomial distribution with 
	   // alpha is the fraction of residues due to the background.
	   // Instead, taking SumXi[j] = mean of the Binomial distribution for alpha = n*p.
	   // where p = a/[a+(1-a)*theta_x]. 
	   SumXI[j] = matFG*(alpha/(alpha+((1-alpha)*theta_x))); // from Equation [3].
	   // SumXI[j] are the true counts (not due to BG matches by chance).
	   if(SumXI[j] != 0 && matFG != 0) assert(SumXI[j] < matFG);
	   // this must be true because require that some of the FG matches (matFG) are due 
	   // to // background contamination (== matFG - SumXI[j]).
        }
}

void pps_typ::BinaryUpdateTheta(UInt4 **CntBG,UInt4 **CntFG)
/************************************************************************
 Conditional on each |Xi(.j)| (i.e., the number of residues in the query set
  for position j. (??)), update alpha according to (1) and update theta_j
  according to (2).
 ************************************************************************/
{
	unsigned char	x,q_j;

        for(Int4 j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		if(q_j == 0) print_error("FATAL: X residues disallowed in query\n");
	        // argum[q_j]=1.0/wt_factor;		// add prior == 100.0
	        double mat =1.0/wt_factor;		// add prior == 100.0
	        double sum=1.0/wt_factor;		// add prior == 100.0
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])){
			mat += CntBG[j][x] + CntFG[j][x];
		   } else {
			sum += CntBG[j][x] + CntFG[j][x];
		   }
                }
		sum += mat;	// Failed to sum this in previous loop.
		// error here?! should not include in sum..
		assert(mat > SumXI[j]);
		mat -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
		sum -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
//              theta[j] = sample_from_dirichlet(argum);
		theta[j][q_j] = mat/sum;
#if 1
		if(DebugCheck && j == 69){ fprintf(stderr,"theta[%d][%c] = %f; mat = %f; sum = %.2f\n",
			j,AlphaChar(q_j,AB),theta[j][q_j],mat,sum); }
#endif
        } 
}

//************************** Binary ************************************

void pps_typ::updateSumXI(UInt4 **CountFG,UInt4 **CntBG)
/************************************************************************
 Conditional on R, alpha and Theta, update the |Xi(.j)| for each column
 according to Equation 3.
 Greek letter Xi(i,j) is an indicator variable for position j of 
	sequence i belonging to the query cluster, where 
	Xi(i,j) == 1 --> x = seq[i][j] comes from foreground component.
	and Xi(i,j) == 0, the background component.
 The theta_x is the fraction of background residues matching the pattern.
 ************************************************************************/
{
        unsigned char	x;

        for(Int4 j=1;j<=k;j++){
	   if(!qst[j]) continue;	// null set here - ignore this position.
	   double matFG=0.0,matBG=0.0,sumBG=0.0;
#if PPS_TYP_FIX_GAP_PROBLEM
	   if(StartAlpha == 0){
		matFG += GetGapsMatch(CountFG[j][0],qst[j]);
		matBG += GetGapsMatch(CntBG[j][0],qst[j]);
		sumBG += CntBG[j][0];	// pseudocounts added below.
	   }
#endif
	   for(x=1; x <= nAlpha(AB); x++){
		if(MemSset(x,qst[j])){
			matFG += CountFG[j][x];
			matBG += CntBG[j][x] + b[x];
		} 
		sumBG += CntBG[j][x] + b[x];
	   } 
	   // theta_x is the fraction of background residues matching the pattern.
	   double theta_x=matBG/sumBG;	
	   // could sample SumXI[j] from the binomial distribution with 
	   // alpha is the fraction of residues due to the background.
	   // Instead, taking SumXi[j] = mean of the Binomial distribution for alpha = n*p.
	   // where p = a/[a+(1-a)*theta_x]. 
	   SumXI[j] = matFG*(alpha/(alpha+((1-alpha)*theta_x))); // from Equation [3].
	   // SumXI[j] are the true counts (not due to BG matches by chance).
	   if(SumXI[j] != 0 && matFG != 0) assert(SumXI[j] < matFG);
	   // this must be true because require that some of the FG matches (matFG) are due 
	   // to // background contamination (== matFG - SumXI[j]).
        }
}

void pps_typ::updateAlpha(UInt4 **CountFG)
/************************************************************************
 Alpha is the parameter specifying the relative fraction of conserved 
 (query pattern set) and unconserved (non-query partition) residues in the 
 query partition.
 The prior for alpha is Beta(a0,b0) distributed where
	a0 specifies the prior pseudocounts for the query pattern set residues.
	b0 specifies the prior pseudocounts for unconserved residues.
 Don't adjust S and N0 for sequence weights as a0 and b0 are unadjusted.
	That is, a0 == 100 and wt_factor = 1/100 --> a0 == 1.
      S is the sum of Xi's over all aligned residues.
 Initialization function to get theta started.
 ************************************************************************/
{
        Int4 		j;
        unsigned char	x;

        double N0=0.0;
        for(j=1;j<=k;j++){
	   if(qst[j]){
#if PPS_TYP_FIX_GAP_PROBLEM
	      if(StartAlpha == 0){
		N0 += GetGapsMatch(CountFG[j][0],qst[j]);
	      }
#endif
	      for(x=1; x <= nAlpha(AB); x++){
		if(MemSset(x,qst[j])) N0 += CountFG[j][x]; 
		// N0 is the total number of residues in the query cluster
		// matching the query residue set at their specific positions.
	      }
	   }
	}
	double S=0.0;	// total foreground matches due to all
        for(j=1; j<=k; j++){ if(qst[j]){ S += SumXI[j]; } }
	alpha = (S+a0)/(a0+N0+b0);  // mean of beta = (S+a0)/(S + a0 + N0 - S + b0). 
	// S + a0 == number of matches due to true selective pressure.
	// N0 - S + b0 == number of matches due to background contamination.
	assert(alpha <= 1.0 && alpha >= 0.0);
	// fprintf(stderr,"alpha = %g; S = %g; N0 = %g\n",alpha,S,N0);
}

void pps_typ::updateTheta(UInt4 **CntBG,UInt4 **CntFG)
/************************************************************************
 Conditional on each |Xi(.j)| (i.e., the number of residues in the query set
  for position j. (??)), update alpha according to (1) and update theta_j
  according to (2).
 ************************************************************************/
{
        double  sum;
	unsigned char	x,q_j;

        for(Int4 j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		if(q_j == 0) print_error("FATAL: X residues disallowed in query\n");
		sum=0.0; argum[q_j]=0.0;
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])){
			if(x != q_j) argum[x]=0.0; // store collapsed residues in q_j
			argum[q_j] += CntBG[j][x] + CntFG[j][x] + b[x];
		   } else {
                        argum[x] = CntBG[j][x] + CntFG[j][x] + b[x];
			sum += argum[x];
		   }
                }
#if PPS_TYP_FIX_GAP_PROBLEM
		if(StartAlpha == 0){
                  for(x=1;x<=nAlpha(AB);x++) {
	   		if(MemSset(x,qst[j])){
			   argum[q_j] += CntBG[j][0]*blosum62freq[x];
			   argum[q_j] += CntFG[j][0]*blosum62freq[x];
			} else {
			   double tmp = CntBG[j][0]*blosum62freq[x];
			   tmp += CntFG[j][0]*blosum62freq[x];
			   argum[x] += tmp; sum += tmp;
			}
	   	  }
	   	}
#endif
		sum += argum[q_j];	// Failed to sum this in previous loop.
		// error here?! should not include in sum..
		assert(argum[q_j] > SumXI[j]);
		argum[q_j] -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
		sum -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
//              theta[j] = sample_from_dirichlet(argum);
		for(x=1;x<=nAlpha(AB);x++){ theta[j][x] = argum[x]/sum; }
        } 
}

#define MODIFIED_MAP	0

double	pps_typ::ProbRatio(unsigned char *seq, unsigned char **SqWt,
		UInt4 **CntBG, UInt4 **CntFG)
// see Gibbs sampling algorithm in section 4 of Jun's JASA paper.
// Step 1.
// Sample sequence i into either the query or the non-query partition...
{
	double		prob;
	Int4		j;
	unsigned char	q_j;
	
	alpha=a0/(a0+b0);
        updateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	updateTheta(CntBG,CntFG);
	for(prob=1.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		q_j=query[j];
		assert(theta[j][q_j] > 0.0 && theta[j][q_j] <= 1.0);
		if(MemSset(seq[j],qst[j])){
		   prob *= (1.0 + alpha*(1-theta[j][q_j])/theta[j][q_j]); 
		} else {
		   prob *= (1.0 - alpha);     
		}
	} return prob;
}

double pps_typ::MapMultiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 From equation 4 in Jun's draft.
 x[seq][col] == Jun's 20-dimensional vector (0,...,0,1,0,...,0) where
 the single 1 indicates the residue type and 
 
 CntBG[j][r] and CntFG[j][r] summarizes the information about the number 
	of main set residues in each column for the main set and query partition
	respectively.
 ************************************************************************/
{
	double		map,*theta_j;
	Int4		j,C;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;
	double		Qj_cntBG,Qj_cntFG;

        if(StartAlpha == 0) print_error("GapSqWts not allowed for Multinomial model");
	alpha=a0/(a0+b0);
	// updateTheta(CntBG);
        updateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	updateTheta(CntBG,CntFG);
	for(C=0,map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;	// skip over non-pattern positions...
		else C++;
		q_j=query[j];
		theta_j = theta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		Qj_cntBG=Qj_cntFG=0.0;
		// double	blsm_Qj=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			Qj_cntBG += CountBG_j[t]; Qj_cntFG += CountFG_j[t];
			// blsm_Qj+=blosum62P[q_j][t];
		   } else {
			map += (log(theta_j[t])*wt_factor*CountBG_j[t] 
				// + log((1-alpha)*blosum62P[q_j][t])*wt_factor*CountFG_j[t]);
				+ log((1-alpha)*theta_j[t])*wt_factor*CountFG_j[t]);
		   }
		}
		map += (log(theta_j[q_j])*wt_factor*Qj_cntBG 
				// + log((1-alpha)*blsm_Qj + alpha)*wt_factor*Qj_cntFG);
				+ log((1-alpha)*theta_j[q_j] + alpha)*wt_factor*Qj_cntFG);
	} 
	map += C*logRho; //  + (k-C)*log(1.0-rho);
	// assert(map != 0); // no non-null pattern positions.
	return map;
}

double pps_typ::NullMapMultiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
  Considers the query residue set at each position to be null.
 ************************************************************************/
{
	double		sum,map,*theta_j;
	Int4		j,t,C;
	UInt4	*CountBG_j,*CountFG_j,q_j;
	unsigned char	x;

     if(StartAlpha == 0) print_error("GapSqWts not allowed for Multinomial model");
     if(ComputeNull){
        for(j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
#if 0
                for(sum=0.0,t=1;t<=nAlpha(AB);t++) {
                        argum[t] = CntBG[j][t] + CntFG[j][t] + b[t];
			sum += argum[t];
                } for(t=1;t<=nAlpha(AB);t++){ 
			theta[j][t] = argum[t]/sum; 
		}
#else		// This computes frequencies based on all residues being in BG partition
		q_j = query[j];
		sum=0.0; argum[q_j]=0.0;
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])){
			if(x != q_j) argum[x]=0.0;
			argum[q_j] += CntBG[j][x] + CntFG[j][x] + b[x];
		   } else {
                        argum[x] = CntBG[j][x] + CntFG[j][x] + b[x];
			sum += argum[x];
		   }
                } sum += argum[q_j];
		// argum[q_j] -= SumXI[j]; // subtract nothing...
		for(x=1;x<=nAlpha(AB);x++){ 
		   if(!MemSset(x,qst[j])) theta[j][x] = argum[x]/sum; 
		   else theta[j][x] = 0.0;
		} theta[j][q_j] = argum[q_j]/sum;
#endif
        } 
	for(map=0.0,C=0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		else C++;
#if 0
		theta_j = theta[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		for(t=1; t<= nAlpha(AB); t++) {
			map += (log(theta_j[t])*wt_factor*(CountBG_j[t]+CountFG_j[t]));
		} 
#else
		q_j=query[j];
		theta_j = theta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double Qj_cntBG;
		Qj_cntBG=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			Qj_cntBG += CountBG_j[t] + CountFG_j[t];
		   } else {
			map += (log(theta_j[t])*wt_factor*(CountBG_j[t]+CountFG_j[t]));
		   }
		} map += (log(theta_j[q_j])*wt_factor*(Qj_cntBG));
#endif
	}
	// map += C*log(1.0-rho);	// C columns changed...
	map += C*logOneMinusRho;	// C columns changed...
	nullMap=map; ComputeNull=FALSE;
    } return nullMap;
}

double pps_typ::subMapMultiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
{
	double		submap,map,*theta_j;
	Int4		j;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	if(StartAlpha == 0) print_error("GapSqWts not allowed for Multinomial model");
	alpha=a0/(a0+b0);
	// updateTheta(CntBG);
        updateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	updateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
	   	if(!qst[j]){ SubMap[j]=0; continue; }
		theta_j = theta[j];
		q_j=query[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double Qj_cntBG,Qj_cntFG;
		Qj_cntBG=Qj_cntFG=0.0;
		submap=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			MatchBG[j]+=wt_factor*CountBG_j[t];
			MatchFG[j]+=wt_factor*CountFG_j[t];
			Qj_cntBG+=CountBG_j[t];
			Qj_cntFG+= CountFG_j[t];
		   } else {
			MisMatchBG[j]+=wt_factor*CountBG_j[t];
			MisMatchFG[j]+=wt_factor*CountFG_j[t];
			submap += (log(theta_j[t])*wt_factor*CountBG_j[t] 
				+ log((1-alpha)*theta_j[t])*wt_factor*CountFG_j[t]);
		   }
		}
		submap += (log(theta_j[q_j])*wt_factor*Qj_cntBG 
				+ log((1-alpha)*theta_j[q_j] + alpha)*wt_factor*Qj_cntFG);
		submap += logRho;
		SubMap[j] = submap;
		map += submap;
	} return map;
}

double pps_typ::subNullMapMultiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
// Considers the query residue set at each position to be null.
{
	double		sum,submap,map,*theta_j;
	Int4		j,t;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	q_j;

	if(StartAlpha == 0) print_error("GapSqWts not allowed for Multinomial model");
        for(j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		sum=0.0; argum[q_j]=0.0;
		Int4 x;
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])){
			if(x != q_j) argum[x]=0.0;
			argum[q_j] += CntBG[j][x] + CntFG[j][x] + b[x];
		   } else {
                        argum[x] = CntBG[j][x] + CntFG[j][x] + b[x];
			sum += argum[x];
		   }
                } sum += argum[q_j];
		for(x=1;x<=nAlpha(AB);x++){ 
		   if(!MemSset(x,qst[j])) theta[j][x] = argum[x]/sum; 
		} theta[j][q_j] = argum[q_j]/sum;
        } 
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		q_j=query[j];
		theta_j = theta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double Qj_cntBG;
		Qj_cntBG=0.0;
		submap=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			Qj_cntBG += CountBG_j[t] + CountFG_j[t];
		   } else {
			submap += (log(theta_j[t])*wt_factor*(CountBG_j[t]+CountFG_j[t]));
		   }
		}
		submap += (log(theta_j[q_j])*wt_factor*(Qj_cntBG));
		// submap += log(1.0-rho);
		submap += logOneMinusRho;
		SubNullMap[j]=submap;
		map += submap;
	} return map;
}

//************************* Collapsed MAPs Below *************************

double pps_typ::MapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 From equation 4 in Jun's draft.
 x[seq][col] == Jun's 20-dimensional vector (0,...,0,1,0,...,0) where
 the single 1 indicates the residue type and 
 
 CntBG[j][r] and CntFG[j][r] summarizes the information about the number 
	of main set residues in each column for the main set and query partition
	respectively.
 ************************************************************************/
{
	double		map,*theta_j;
	Int4		j;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	alpha=a0/(a0+b0);
        // updateSumXI(CntFG,CntBG);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	// updateTheta(CntBG,CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
	   	// if(!qst[j]){ map += log(1.0 - rho); continue; }
	   	if(!qst[j]){ continue; } // cancels out with null map...
		q_j=query[j];
		theta_j = theta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1,n1,m2,n2;
		m1=n1=m2=n2=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			m1 += CountBG_j[t]; n1 += CountFG_j[t];
		   } else {
			m2 += CountBG_j[t]; n2 += CountFG_j[t];
		   }
		}
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  double tmp = GetGapsMatch(CountFG_j[0],qst[j]);
		  n1 += tmp; 
		  n2 += CountFG_j[0] - tmp;

		  tmp = GetGapsMatch(CountBG_j[0],qst[j]);
		  m1 += tmp; 
		  m2 += CountBG_j[0] - tmp;
	   	}
#endif	// ######################################################
#if 0		// pre-Mathematica version
		map += (log(theta_j[q_j])*wt_factor*m1
				+ log((1-alpha)*theta_j[q_j] + alpha)*wt_factor*n1);
		map += (log(1.0-theta_j[q_j])*wt_factor*m2
				+ log((1-alpha)*(1.0 - theta_j[q_j]))*wt_factor*n2);
#else		// Post-Mathematica version
		n1 = wt_factor*n1; m1 = wt_factor*m1;
	 	n2 = wt_factor*n2; m2 = wt_factor*m2;
		map += n1*log((1-alpha)*theta_j[q_j] + alpha) + m1*log(theta_j[q_j]); 
		map += n2*log(1-alpha) + (n2+m2)*log(1.0-theta_j[q_j]);
		map += logRho;	// column penalty...
#endif
	} return map;
}

double pps_typ::NullMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
  Considers the query residue set at each position to be null.
 ************************************************************************/
{
	double		map;
	Int4		j,t;
	UInt4	*CountBG_j,*CountFG_j,q_j;
	unsigned char	x;

     if(ComputeNull){
        for(j=1;j<=k;j++) {	// compute theta_j's here...
	   	if(!qst[j]) continue;	// these columns cancel out for Net-LPR
		q_j = query[j];
		double sum=1.0/wt_factor;
		double match=1.0/wt_factor;
		// argum[q_j]=0.0;
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])) match += CntBG[j][x] + CntFG[j][x];
		   else sum += CntBG[j][x] + CntFG[j][x];
                } sum += match;
#if PPS_TYP_FIX_GAP_PROBLEM
		if(StartAlpha == 0){
		  print_error("PPS_TYP_FIX_GAP_PROBLEM not yet implemented");
                  for(x=1;x<=nAlpha(AB);x++) {
	   		if(MemSset(x,qst[j])){
			   argum[q_j] += CntBG[j][0]*blosum62freq[x];
			   argum[q_j] += CntFG[j][0]*blosum62freq[x];
			} else {
			   double tmp = CntBG[j][0]*blosum62freq[x];
			   tmp += CntFG[j][0]*blosum62freq[x];
			   argum[x] += tmp; sum += tmp;
			}
	   	  }
	   	}
#endif
		NullTheta[j] = match/sum;
        } 
	for(map=0.0,j=1; j<=k; j++) {
		// use same pattern but put all sequences into the background.
	   	if(!qst[j]) continue;	// these columns cancel out for Net-LPR
		q_j=query[j];
		double theta_x = NullTheta[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		double m1=0.0,m2=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){ m1+= CountBG_j[t] + CountFG_j[t]; }
		   else { m2 += CountBG_j[t] + CountFG_j[t]; }
		} 
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  print_error("PPS_TYP_FIX_GAP_PROBLEM not yet implemented");
		  double tmp = GetGapsMatch(CountFG_j[0],qst[j]);
		  tmp += GetGapsMatch(CountBG_j[0],qst[j]);
		  m1 += tmp; 
		  m2 += CountFG_j[0] + CountBG_j[0] - tmp;
	   	}
#endif	// ######################################################
		m1=wt_factor*m1; m2=wt_factor*m2;
		map += (m1*log(theta_x));
		map += (m2*log((1.0-theta_x)));
		// map += log(1.0-rho);	// columns out for null model.
		map += logOneMinusRho;	// columns out for null model.
	} nullMap=map; ComputeNull=FALSE;
    } return nullMap;
}

double pps_typ::subMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
{
	double		submap,map,*theta_j;
	Int4		j;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	alpha=a0/(a0+b0);
        // updateSumXI(CntFG,CntBG);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	// updateTheta(CntBG,CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
	   	if(!qst[j]){ SubMap[j]=0; continue; }
		theta_j = theta[j];
		q_j=query[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1,n1,m2,n2;
		m1=n1=m2=n2=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			// MatchBG[j]+=wt_factor*CountBG_j[t]; 
			// MatchFG[j]+=wt_factor*CountFG_j[t]; 
			n1+= CountFG_j[t]; m1+=CountBG_j[t];
		   } else {
			// MisMatchBG[j]+=wt_factor*CountBG_j[t]; 
			// MisMatchFG[j]+=wt_factor*CountFG_j[t]; 
			n2+= CountFG_j[t]; m2+=CountBG_j[t];
		   }
		}
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  print_error("PPS_TYP_FIX_GAP_PROBLEM not yet implemented");
		  double tmp = GetGapsMatch(CountBG_j[0],qst[j]);
		  MatchBG[j]+=wt_factor*tmp; m1 += tmp; 
		  MisMatchBG[j]+=wt_factor*(CountBG_j[0]-tmp); m2 += CountBG_j[0]-tmp;
		  
		  tmp = GetGapsMatch(CountFG_j[0],qst[j]);
		  MatchFG[j]+=wt_factor*tmp; n1 += tmp; 
		  MisMatchFG[j]+=wt_factor*(CountFG_j[0]-tmp); n2 += CountFG_j[0]-tmp;
	   	}
#endif	// ######################################################
		submap=0.0;
		// n1 = matching FG; m1 = matching BG
		// n2 = mis-match FG; m2 = mismatch BG.
		m1=m1*wt_factor; m2=m2*wt_factor;
		n1=n1*wt_factor; n2=n2*wt_factor;
		MatchFG[j] = n1; MisMatchFG[j] = n2;
		MatchBG[j] = m1; MisMatchBG[j] = m2;
#if 0		// Pre-Mathematica version
		submap += n1*log((1-alpha)*theta_j[q_j] + alpha);
		submap += n2*log((1-alpha)*(1.0-theta_j[q_j]));
		submap += m1*log(theta_j[q_j]);
		submap += m2*log((1.0-theta_j[q_j]));

			// + log((1-alpha)*(1.0-theta_j[q_j]) + alpha)*n2);
#else		// Post-Mathematica version
		//*************** Rj terms: ***************
		// matching FG residues; this is theta^alpha. 
		// (1-alpha)*theta_j[q_j] is the fraction of BG matching the query.
		// alpha*1.0 is the fraction of FG matching residues due to chance.
		submap += n1*log((1-alpha)*theta_j[q_j] + alpha); // delta_j == 1.
		//               FG 'contamination'       observed
		// Mismatching FG residues; assumes all are due to contamination.
		submap += n2*log((1-alpha)*(1.0-theta_j[q_j]));	  // delta_j == 0.
		// submap += n2*log(1-alpha);

		//*************** (1-Rj) terms: *************** 
		// matching BG residues; assumes none are due to non-contamination.
		submap += m1*log(theta_j[q_j]);	
		// Mismatching BG residues; assums all are due to theta
		submap += m2*log((1.0-theta_j[q_j]));
		// submap += (n2+m2)*log(1.0-theta_j[q_j]); moved part of term above below here...
#endif
		submap += logRho;
		SubMap[j] = submap;
		map += submap;
	} return map;
}

double pps_typ::subNullMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
// Considers the query residue set at each position to be null.
{
	double		submap,map;
	Int4		j,t;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	q_j;

        for(j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		double sum=1.0*wt_factor;
		double match=1.0*wt_factor;
		Int4 x;
                for(x=1;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])) match += CntBG[j][x] + CntFG[j][x];
		   else sum += CntBG[j][x] + CntFG[j][x];
                } sum += match;
#if PPS_TYP_FIX_GAP_PROBLEM
		if(StartAlpha == 0){
		  print_error("PPS_TYP_FIX_GAP_PROBLEM not yet implemented");
                  for(x=1;x<=nAlpha(AB);x++) {
	   		if(MemSset(x,qst[j])){
			   argum[q_j] += CntBG[j][0]*blosum62freq[x];
			   argum[q_j] += CntFG[j][0]*blosum62freq[x];
			} else {
			   double tmp = CntBG[j][0]*blosum62freq[x];
			   tmp += CntFG[j][0]*blosum62freq[x];
			   argum[x] += tmp; sum += tmp;
			}
	   	  }
	   	}
#endif
		NullTheta[j] = match/sum;
        } 
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		q_j=query[j];
		double theta_x = NullTheta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1=0.0,m2=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			m1+= CountBG_j[t] + CountFG_j[t];
		   } else {
			m2+= CountBG_j[t] + CountFG_j[t];
		   }
		}
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  print_error("PPS_TYP_FIX_GAP_PROBLEM not yet implemented");
		  double tmp = GetGapsMatch(CountFG_j[0],qst[j]);
		  tmp += GetGapsMatch(CountBG_j[0],qst[j]);
		  m1 += tmp; m2 += CountFG_j[0] + CountBG_j[0] - tmp;
	   	}
#endif	// ######################################################
		submap=0.0;
		m1 = m1*wt_factor; m2 = m2*wt_factor;
		submap += m1*log(theta_x) + m2*log((1.0-theta_x));
		// submap += log(1.0 - rho);
		submap += logOneMinusRho;
		map += submap;
		SubNullMap[j]=submap;
	} return map;
}

double pps_typ::MapBallInUrn(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 Ball in Urn model:
 ************************************************************************/
{
	double		map,m1,n1,m2,n2;
	UInt4		j,*countBG_j,*countFG_j;
	unsigned char	t;

	// fprintf(stderr,"************* MapBallInUrn( ) **************\n");
	// fprintf(stderr," ModelMode = %c; ModelSwitch = %c\n",ModelMode,ModelSwitch);
        updateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	if(alpha >= 1.0) alpha=0.999999999999999;
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		countBG_j = CntBG[j]; countFG_j = CntFG[j];
		m1=n1=m2=n2=0.0;
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			m1 += wt_factor*countBG_j[t]; n1 += wt_factor*countFG_j[t]; 
		   } else {
			m2 += wt_factor*countBG_j[t]; n2 += wt_factor*countFG_j[t]; 
		   }
		}
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  double tmp = GetGapsMatch(countBG_j[0],qst[j]);
		  m1 += wt_factor*tmp; 
		  m2 += wt_factor*(countBG_j[0]-tmp);

		  tmp = GetGapsMatch(countFG_j[0],qst[j]);
		  n1 += wt_factor*tmp; 
		  n2 += wt_factor*(countFG_j[0]-tmp);
	   	}
#endif	// ######################################################
		map+=submap_ball_in_urn(m1,m2,n1,n2);
	}
#if 0
	fprintf(stderr,"LPR = %g\n",map);
	map = floor(map+0.5);
	fprintf(stderr,"floored LPR = %g\n",map);
	return map;
#else
	// return floor(map+0.5);
	return map;
#endif
}

double pps_typ::subMapBallInUrn(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 Ball in Urn model:
 ************************************************************************/
{
	double		submap,map;
	UInt4		j,*countBG_j,*countFG_j;
	unsigned char	t;

        updateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	if(alpha >= 1.0) alpha=0.999999999999999;
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]) continue;
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
		countBG_j = CntBG[j]; countFG_j = CntFG[j];
		for(t=1; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			MatchBG[j]+=wt_factor*countBG_j[t];
			MatchFG[j]+=wt_factor*countFG_j[t];
		   } else {
			MisMatchBG[j]+=wt_factor*countBG_j[t];
			MisMatchFG[j]+=wt_factor*countFG_j[t];
		   }
		}
#if PPS_TYP_FIX_GAP_PROBLEM // ##############################################
	   	if(StartAlpha == 0){
		  double tmp = GetGapsMatch(countBG_j[0],qst[j]);
		  MatchBG[j]+=wt_factor*tmp; 
		  MisMatchBG[j]+=wt_factor*(countBG_j[0]-tmp); 
		  
		  tmp = GetGapsMatch(countFG_j[0],qst[j]);
		  MatchFG[j]+=wt_factor*tmp;
		  MisMatchFG[j]+=wt_factor*(countFG_j[0]-tmp);
	   	}
#endif	// ######################################################
		SubMap[j] = submap_ball_in_urn(MatchBG[j],MisMatchBG[j],
				MatchFG[j],MisMatchFG[j]);
		map+= SubMap[j]; SubNullMap[j]=0.0;
	}
	// return floor(map+0.5);
	return map;
}

const double TT=2.0;

double	pps_typ::BinomialWtdRE(double N, double n, double q, double p)
{
	return (BinomialDistribution(N,q,n)
			*(n*log(q) + (N-n)*log(1.0-q) - n*log(p) - (N-n)*log(1.0-p)));
}

double	pps_typ::CBKLD(double N, double n, double q, double p)
{

	double P=CumBinomProb(n,N,q)*(LnCBP(n,N,q)-LnCBP(n,N,p)); // selection for 
	P += CumBinomProb(N-n,N,1-q)*(LnCBP(N-n,N,1-q)-LnCBP(N-n,N,1-p)); // selection against 
	return P;
}

double pps_typ::submap_ball_in_urn(double m1, double m2, double n1, double n2)
{
	double		submap,p,q,pq,N,M,NN,MM;
	static long	calls=0;

	calls++;
	// n1 = wt_factor*n1; m1 = wt_factor*m1; n2 = wt_factor*n2; m2 = wt_factor*m2;
	N = n1+n2; M=m1+m2;
	p = ((1+m1)/(2+M)); q = ((1+n1)/(2+N));
	pq=((1+m1+n1)/(2+M+N));

	NN=(double) wt_factor*(double)CardSet[1];  // Wt number of FG seqs with & without deletions
	MM=(double) wt_factor*(double)CardSet[2];	// BG set.
	submap=0.0;
	double adjust=0.0;
   switch (ModelSwitch){
   case '0':	//  adjusted CBKLD model...
     {
	if(q != p) submap = 0.5*(CBKLD(N, n1,q,p) + CBKLD(M, m1,p,q));
	if(q < p) submap = -submap;
	adjust = 0.5*(CBKLD(NN, NN*sigma,sigma,1.0-sigma) 
				+ CBKLD(MM, MM*(1.0-sigma),1.0-sigma,sigma));
	if(sigma <= 0.5) adjust = -adjust;      // because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
     }
    break;
   case '1':	// adjusted CBLPR model ...
      {
	if(q > p) {		// positive selection for n1.
	   submap=0.5*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG selective pressure
	   submap+=0.5*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG selective pressure
	} else if(q < p){	// selection against n1.
	   submap=0.5*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG selective pressure
	   submap+=0.5*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // BG selective pressure
	   submap = -submap;	// negative pressure
	} else submap=0.0;
	adjust = 0.5*(LnCBP(NN*sigma,NN,sigma)-LnCBP(NN*sigma,NN,1-sigma));
	adjust += 0.5*(LnCBP(MM*sigma,MM,sigma)-LnCBP(MM*sigma,MM,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '2':	//  adjusted BKLD model...
     {
	if(q != p) submap = 0.5*(BKLD(N, n1,q,p) + BKLD(M, m1,p,q));
	if(q < p) submap = -submap;
	adjust = 0.5*(BKLD(NN, NN*sigma,sigma,1.0-sigma) 
				+ BKLD(MM, MM*(1.0-sigma),1.0-sigma,sigma));
	if(sigma <= 0.5) adjust = -adjust;      // because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
     }
    break;
   case '3':	//  Wally's model...
      {
	double OR,se;
	OR=log(n1) + log(m2) - log(n2) -log(m1);
	se=1.0/n1 + 1.0/n2 + 1.0/m1 + 1.0/m2;
	se=sqrt(se);
	submap = OR/se;
	// submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case 't':	//  adjusted BLPR model...
      {
	if(q != p) {
	   submap=0.5*(BKLDfactor(N,n1,q)-BKLDfactor(N,n1,p)); // FG selective pressure
	   submap+=0.5*(BKLDfactor(M,m2,1-p)-BKLDfactor(M,m2,1-q)); // BG selective pressure
	} if(q < p) submap = -submap;
	adjust = 0.5*(BKLDfactor(NN,NN*sigma,sigma)-BKLDfactor(NN,NN*sigma,1-sigma));
	adjust += 0.5*(BKLDfactor(MM,MM*sigma,sigma)-BKLDfactor(MM,MM*sigma,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '4':	// adjusted CBLPR model with additional factors...
      {
	if(q > p) {		// positive selection for n1.
	   submap=0.5*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG selective pressure
	   submap+=0.5*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG selective pressure
	   submap+=0.5*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // BG selective pressure
	   submap+=0.5*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG selective pressure
	} else if(q < p){	// selection against n1.
	   submap=0.5*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG selective pressure
	   submap+=0.5*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG selective pressure
	   submap+=0.5*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG selective pressure
	   submap+=0.5*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // BG selective pressure
	   submap = -submap;	// negative pressure
	} else submap=0.0;
	adjust = 0.5*(LnCBP(NN*sigma,NN,sigma)-LnCBP(NN*sigma,NN,1-sigma));
	adjust += 0.5*(LnCBP(MM*sigma,MM,sigma)-LnCBP(MM*sigma,MM,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '5':	//  adjusted BLPR model...
      {
	if(q != p) {
	   submap=0.5*(BKLDfactor(N,n1,q)-BKLDfactor(N,n1,p)); // FG selective pressure
	   submap+=0.5*(BKLDfactor(N,n2,1-q)-BKLDfactor(N,n2,1-p)); // FG selective pressure
	   submap+=0.5*(BKLDfactor(M,m2,1-p)-BKLDfactor(M,m2,1-q)); // BG selective pressure
	   submap+=0.5*(BKLDfactor(M,m1,p)-BKLDfactor(M,m1,q)); // BG selective pressure
	} if(q < p) submap = -submap;
	adjust = 0.5*(BKLDfactor(NN,NN*sigma,sigma)-BKLDfactor(NN,NN*sigma,1-sigma));
	adjust += 0.5*(BKLDfactor(MM,MM*sigma,sigma)-BKLDfactor(MM,MM*sigma,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '6':	//  adjusted/equally weighted CBKLD model...
     {
	// if(q != p) submap = WtdCumBnmlKLDivergence(N,n1,M,m1,1.0,1.0);
	if(q != p) submap = (M/(M+N))*CBKLD(N, n1,q,p) + (N/(M+N))*CBKLD(M, m1,p,q);
	else submap=0.0;
	if(q < p) submap = -submap;
	adjust = (MM/(MM+NN))*CBKLD(NN, NN*sigma,sigma,1.0-sigma) 
				+ (NN/(MM+NN))*CBKLD(MM, MM*(1.0-sigma),1.0-sigma,sigma);
	if(sigma <= 0.5) adjust = -adjust;      // because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
     }
    break;
   case '7':	// adjusted && equally-weighted CBLPR model ...
      {
	if(q > p) {		// positive selection for n1.
	   submap=(M/(M+N))*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG selective pressure
	   submap+=(N/(M+N))*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG selective pressure
	} else if(q < p){	// selection against n1.
	   submap=(M/(M+N))*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG selective pressure
	   submap+=(N/(M+N))*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // BG selective pressure
	   submap = -submap;	// negative pressure
	} else submap=0.0;
	adjust = (MM/(MM+NN))*(LnCBP(NN*sigma,NN,sigma)-LnCBP(NN*sigma,NN,1-sigma));
	adjust += (NN/(MM+NN))*(LnCBP(MM*sigma,MM,sigma)-LnCBP(MM*sigma,MM,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '8':	// adjusted && equally-weighted CBLPR model ...
      {
	if(q != p) {		// positive selection for n1.
	  submap=(M/(M+N))*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG selective pressure
	  submap+=(M/(M+N))*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG selective pressure
	  submap+=(N/(M+N))*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // BG selective pressure
	  submap+=(N/(M+N))*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG selective pressure
	  if(q < p) submap = -submap;	
	} else submap = 0.0;
	adjust = (MM/(MM+NN))*(LnCBP(NN*sigma,NN,sigma)-LnCBP(NN*sigma,NN,1-sigma));
	adjust += (NN/(MM+NN))*(LnCBP(MM*sigma,MM,sigma)-LnCBP(MM*sigma,MM,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case '9':	//  adjusted/equally weighted CBKLD model...try to eliminate biase against BG.
     {
	double P;
	if(q != p){
	 P = CumBinomProb(n1,N,q)*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // selection for 
	 P += CumBinomProb(n2,N,1-q)*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // selection against 
	 submap = (M/(M*N))*P;
	 P = CumBinomProb(m1,M,p)*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); // selection for 
	 P += CumBinomProb(m2,M,1-p)*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // selection against 
	 submap += (N/(M*N))*P;
	} else submap=0.0;
	if(q < p) submap = -submap;
	adjust = (MM/(MM*NN))*CBKLD(NN, NN*sigma,sigma,1.0-sigma) 
				+ (NN/(MM*NN))*CBKLD(MM, MM*(1.0-sigma),1.0-sigma,sigma);
	if(sigma <= 0.5) adjust = -adjust;      // because switches signs in this region
	submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
     }
    break;
   case 'D':	// equally-weighted cbp cross entropy... H(q||p) = sum q log 1/p
      {		// Made this symmetric. doesn't work?
	// if(q > p) {		// positive selection for n1.
	   submap=(M/(M+N))*CumBinomProb(n1,N,q)*(-LnCBP(n1,N,p));
	   submap+=(M/(M+N))*CumBinomProb(n2,N,1-q)*(-LnCBP(n2,N,1-p));
	   submap+=(N/(M+N))*CumBinomProb(m2,M,1-p)*(-LnCBP(m2,M,1-q)); // BG selective pressure
	   submap+=(N/(M+N))*CumBinomProb(m1,M,p)*(-LnCBP(m1,M,q)); // BG selective pressure
	// else if(q < p){	// selection against n1.
	if(q < p){	// selection against n1.
	   // submap=(M/(M+N))*CumBinomProb(n2,N,1-q)*(-LnCBP(n2,N,1-p));
	   // submap+=(N/(M+N))*CumBinomProb(m1,M,p)*(-LnCBP(m1,M,q)); // BG selective pressure
	   // submap = -submap;	// negative pressure
	} // else submap=0.0;
	// n1 = NN*sigma; n2 = NN*(1-sigma); m1 = MM*(1-sigma);
	adjust = (MM/(MM+NN))*(CumBinomProb(NN*sigma,NN,sigma)*(-LnCBP(NN*sigma,NN,1-sigma)));
	adjust += (MM/(MM+NN))*(CumBinomProb(NN*(1-sigma),NN,1-sigma)*(-LnCBP(NN*(1-sigma),NN,sigma)));
	adjust += (NN/(MM+NN))*(CumBinomProb(MM*sigma,MM,sigma)*(-LnCBP(MM*sigma,MM,1-sigma)));
	adjust += (NN/(MM+NN))*(CumBinomProb(MM*(1-sigma),MM,1-sigma)*(-LnCBP(MM*(1-sigma),MM,sigma)));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	// submap -= adjust;
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case 'V':	// Dual CBP model.  p = ((1+m1)/(2+M)); q = ((1+n1)/(2+N));
	if(q > p){
		submap += LnCBP(n1,N,q) - LnCBP(n1,N,p); // positive FG selective pressure
		submap += LnCBP(m2,M,(1-p)) - LnCBP(m2,M,(1-q)); // negative BG selective pressure
	} else if(q == p) submap += 0.0;	// no selective pressure
	else {
		submap += LnCBP(m1,M,p) - LnCBP(m1,M,q);	// negative selective pressure
		submap += LnCBP(n2,N,(1-q)) - LnCBP(n2,N,(1-p)); // positive selective pressure
	}
    break;
//*******************************************************************************************
   case 'Y':	// add on weighted mixture model (not tried yet)
	// MAY BE THE BEST...
	{
	   submap += (n1*log((1-alpha)*p + alpha) + m1*log(p));
	   submap += (n2*log(1-alpha) + (n2+m2)*log(1.0-p));
	   submap -= (n1+ m1)*log(pq);
	   submap -= (n2+m2)*log(1.0-pq);
	}
	if(q > p) submap = (1.0-(q-p))*submap;
	if(q > p) submap -= (q-p)*LnBinomialDistribution(N,p,n1);
    break;
   case 'Z':	//  not tried 
      {
	double adjust=0.0;
	if(q != p) {
	   submap = 0.5*(LnBinomialDistribution(N,q,n1) - LnBinomialDistribution(N,p,n1));
	   submap += 0.5*(LnBinomialDistribution(M,1-p,m2) - LnBinomialDistribution(M,1-q,m2));
	} if(q < p) submap = -submap;
	adjust -= 0.5*(LnBinomialDistribution(N*sigma,N,sigma)
				-LnBinomialDistribution(N*sigma,N,1-sigma));
	adjust -= 0.5*(LnBinomialDistribution(M*sigma,M,sigma)
				-LnBinomialDistribution(M*sigma,M,1-sigma));
	if(sigma <= 0.5) adjust = -adjust;	// because switches signs in this region
	submap -= log(k*1000.0);	// Bonferoni adjustment...
      }
    break;
   case 'W':	// one sided weighted mixture model (works fairly well).
	// BEST SO FAR...
	{	
	   submap += (n1*log((1-alpha)*p + alpha) + m1*log(p));
	   submap += (n2*log(1-alpha) + (n2+m2)*log(1.0-p));
	   submap -= (n1+ m1)*log(pq);
	   submap -= (n2+m2)*log(1.0-pq);
	} if(q > p) submap -= pow(q-p,2.0)*LnBinomialDistribution(N,p,n1);
    break;
   case 'M':	// one sided weighted mixture model (works fairly well).
	// BEST SO FAR...
	// seems to get trapped in infinite loops
	{	
	   submap += (n1*log((1-alpha)*p + alpha) + m1*log(p));
	   submap += (n2*log(1-alpha) + (n2+m2)*log(1.0-p));
	   submap -= (n1+ m1)*log(pq);
	   submap -= (n2+m2)*log(1.0-pq);
	}
	if(q > p) submap = (1.0-pow((q-p),2.0))*submap;
	if(q > p) submap -= pow(q-p,2.0)*LnBinomialDistribution(N,p,n1);
    break;
   case 'P':	// Dual CBP model.  p = ((1+m1)/(2+M)); q = ((1+n1)/(2+N));
	// APPEARS TO BE THE best yet.
	if(q > p){	// positive selective pressure on FG residue set n1 relative to BG
		submap += LnCBP(n1,N,q) - LnCBP(n1,N,p); // positive pressure on FG
		submap += LnCBP(m2,M,(1-p)) - LnCBP(m2,M,(1-q)); // + negative pressure on BG
		// == sum(0..m1) (N choose i) p^i * (1-p)^(N-i) == other tail prob...
	} else if(q == p) submap += 0.0;	// no relative selective pressure
	else {		// negative selective pressure on residue set n1 relative to BG
		submap += LnCBP(n2,N,(1-q)) - LnCBP(n2,N,(1-p)); // positive selective pressure
		// == sum(0..n1) CBP  (tail at other end...
		submap += LnCBP(m1,M,p) - LnCBP(m1,M,q);	// negative selective pressure
		// subtract these... --> (1-q) > (1-p)
	}
    break;
   case 'A':	// one sided weighted mixture model (works fairly well).
     {
	if(q != p) {
		submap = 0.5*(CBKLD(NN, n1,q,p) + CBKLD(MM, m1,p,q));
	} if(q < p) submap = -submap;
	submap -= 0.5*(CBKLD(NN, NN*sigma,sigma,1.0-sigma) 
				+ CBKLD(MM, MM*(1.0-sigma),1.0-sigma,sigma));
	submap -= log(k*1000.0);	// Bonferoni adjustment...
     }
    break;
   case 'B':	// BINOMIAL RELATIVE ENTROPY MEASURE.
    {		// BEST MODEL...but need to give equal weight to both?
	submap = BinomialWtdRE(N,n1,q,p) + BinomialWtdRE(M,m2,1-p,1-q);
	if(q < p) submap = -submap;
	// fprintf(stderr,"BNRE(%d,%d,%d,%d,%g,%g) = %g\n",
	// (Int4)N,(Int4)n1,(Int4)M,(Int4)m2,q,p,submap);
	submap -= log(k*1000);
     }
    break;
   case 'X':	// one sided weighted mixture model (works fairly well).
	// BEST SO FAR...
     {
	{	
	   submap += (n1*log((1-alpha)*p + alpha) + m1*log(p));
	   submap += (n2*log(1-alpha) + (n2+m2)*log(1.0-p));
	   submap -= (n1+ m1)*log(pq);
	   submap -= (n2+m2)*log(1.0-pq);
	}
	if(q > p) submap = (1.0-pow((q-p),2.0))*submap;
	double cbpmap = -LnBinomialDistribution(N,p,n1) + LnBinomialDistribution(M,(1-q),m2);
	// double cbpmap = LnCBP(n1,N,q) - LnCBP(n1,N,p); // positive pressure on FG
	// cbpmap += LnCBP(m2,M,(1-p)) - LnCBP(m2,M,(1-q)); // + negative pressure on BG
	if(q > p) submap -= pow(q-p,2.0)*cbpmap;
     }
    break;
   default: print_error("Urn model input error"); break;
   }
	// submap += log(rho/(1.0-rho));
	submap += logRho - logOneMinusRho;
	return submap;
}

double	*pps_typ::SubLPR(UInt4 **CntBG,UInt4 **CntFG)
{
	double	map0,map;
	map=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	subLPR[0]=map;
	map0=0.0;
	for(Int4 j=1; j <= k; j++){
	   if(qst[j]){
	     subLPR[j]=SubMap[j]-SubNullMap[j];
	   } else subLPR[j]=0.0;
	   map0+=subLPR[j];
	}
#if 0
	if(ModelMode=='U') map0=floor(map0+0.5);
	subLPR[0]=map0;	// TRY TO FIX ROUNDING ERROR!
	if(map != map0) fprintf(stderr,"map=%.2f;map0=%.2f\n",map,map0);
#endif
	return subLPR;
#if 0		// Original
	double	map0,map,sub;
	map=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	subLPR[0]=map;
	fprintf(stdout,
		"   Pattern:   MatFG  MisFG  % MatBG  MisBG  % subLPR\n");
		map0=0.0;
	for(Int4 j=1; j <= k; j++){
	   if(qst[j]){
	     sub=SubMap[j]-SubNullMap[j];
	     char tmp[30],*ptr; ptr=tmp;
	     for(Int4 r=StartAlpha; r <= nAlpha(AB); r++) 
		if(MemSset(r,qst[j])){ sprintf(ptr,"%c",AlphaChar(r,AB)); ptr++;}
	     sprintf(ptr,"%d",j); 
	     fprintf(stdout,
		"%10s: %5.0f %5.0f %5.0f %5.0f %5.1f (%.1f - %.1f)\n",
		tmp,MatchFG[j],MisMatchFG[j],
		MatchBG[j],MisMatchBG[j],
		sub,SubMap[j],SubNullMap[j]);
	     subLPR[j]=sub;
	   } else subLPR[j]=0.0;
	   map0+=subLPR[j];
	}
	fprintf(stdout," total map = %.1f (%.1f)(diff %g)\n",
			map,map0,map-map0);
	subLPR[0]=map0;	// TRY TO FIX ROUNDING ERROR!
	return subLPR;
#endif
}

void	pps_typ::PutSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG)
{
	 double	map=0.0,map0=0.0,sub;
	 // DebugCheck = TRUE;
	 map0=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	 // DebugCheck = FALSE;
	// this sets MatchFG, etc...
	 Int4 j,Min,Max,C;
	//********************* Determine range of histogram *************
	 double min=DBL_MAX,max=-9999999.0,Inc;
	 dh_type dH_LPR=dheap(k+2,4),dH=0;
	 dh_type dH_BN=dheap(k+2,4);
	 // char mode='E';
	 char mode=AltMeasure;
	 char *model_str=0;
	 char *units_str=0;
	 switch (mode){
		     case 'K': 		// KL-Divergence
			model_str=NewString("CBKLD");
			units_str=NewString("nats(rank)");
		     break;
		     case 'B': 		// LnBinomialDistribution...
			model_str=NewString("LnBinomial");
			units_str=NewString("ln(p)(rank)");
		     break;
		     case 'C': 	// CBLPR equally weighted binary model (option -method=7)
			model_str=NewString("CBLPR");
			units_str=NewString("ln(p)(rank)");
		     break;
		     case 'R': 	// weighted log-probability ratio...
			model_str=NewString("BallInUrnModel");
			units_str=NewString("ln(p)(rank)");
		     break;
		     case 'E': 	// binary relative entropy...
			model_str=NewString("Relative Entropy");
			units_str=NewString("0.001 nats(rank)");
		     break;
		     case 'U': 	// simply binary Urn model
			model_str=NewString("BinaryUrn");
			units_str=NewString("ln(p)(rank)");
		     break;
	 }
	 for(j=1; j <= k; j++){
		   if(qst[j]){
		     insrtHeap(j,(keytyp)-subLPR[j],dH_LPR);
		     double N,M,p,q,n1,n2,m1,m2,P;
		     n1=MatchFG[j]; n2=MisMatchFG[j];
		     N=n1+n2; q=(n1+1)/(N+2);

		     m1 = MatchBG[j]; m2 = MisMatchBG[j];
		     M=m1+m2; p=(m1+1)/(M+2);

		     switch (mode){
		     case 'K': 		// KL-Divergence
		     	P=0.5*(CBKLD(N,n1,q,p) + CBKLD(M,m1,p,q));
		        if(q < p) P = -P;
        	     	P -= 0.5*(CBKLD(N, N*sigma,sigma,1.0-sigma)
                        	+ CBKLD(M, M*(1.0-sigma),1.0-sigma,sigma));
        	     	// P -= log(k*1000.0);        // Bonferoni adjustment...
		      break;
		     case 'B': 		// P= LnBinomialDistribution(N,q,n1);
		     	P= LnBinomialDistribution(N,p,n1);
		      break;
		     case 'C': 	// CBLPR equally weighted binary model (option -method=7)
		       if(q > p){
          		P=(M/(M+N))*(LnCBP(n1,N,q)-LnCBP(n1,N,p)); // FG pressure 
          		P+=(N/(M+N))*(LnCBP(m2,M,1-p)-LnCBP(m2,M,1-q)); // BG pressure
		       } else if(q < p){
			P=(M/(M+N))*(LnCBP(n2,N,1-q)-LnCBP(n2,N,1-p)); // FG pressure
			P+=(N/(M+N))*(LnCBP(m1,M,p)-LnCBP(m1,M,q)); 
				// BG pressure: prob obs >= m2
			P = -P;       // negative selective pressure...
		       } else P=0.0;
		      break;
		     case 'R': 	// weighted log-probability ratio...
		     	P= (N/(M+N))*(LnCBP(n1,N,q) - LnCBP(n1,N,p));
		     	P+= (M/(M+N))*(LnCBP(m2,M,1-p) - LnCBP(m2,M,1-q));
		      break;
		     case 'E': 
			{
			  P=1000*(p*log(p/q) + (1-p)*log((1-p)/(1-q)));
			  P+=1000*(q*log(q/p) + (1-q)*log((1-q)/(1-p)));
			} break;
		     case 'U': P= -LnCBP(n1,N,p); break;
		     default: print_error("pps_typ input error"); break;
		     }
		     insrtHeap(j,(keytyp)-P,dH_BN);
		     min=MINIMUM(double,subLPR[j],min);
		     max=MAXIMUM(double,subLPR[j],max);
		   }
	 }
	 Int4	*rankLPR,*rankBNPD,rank;
	 double *keyBNPD;
	 keytyp	key;
	 NEW(rankLPR,k+3,Int4); 
	 NEW(rankBNPD,k+3,Int4);
	 NEW(keyBNPD,k+3,double);
	 dH=dheap(k+2,4);
	 for(rank=1;!emptyHeap(dH_LPR); rank++){
		   key=minkeyHeap(dH_LPR);
		   j=delminHeap(dH_LPR); rankLPR[j]=rank; 
		   insrtHeap(j,key,dH);
	 } Nildheap(dH_LPR); dH_LPR=dH;
	 dH=dheap(k+2,4);
	 for(rank=1;!emptyHeap(dH_BN); rank++){
		   key=minkeyHeap(dH_BN);
		   j=delminHeap(dH_BN); rankBNPD[j]=rank; keyBNPD[j]=-key;
		   insrtHeap(j,key,dH);
	 } Nildheap(dH_BN);  dH_BN=dH;
	 Min=(Int4) floor(min);
	 Max=(Int4) ceil(max);
	 Inc = floor((max-min)/30.0);
	 if(Inc < 1.0) Inc=0.5;
	 // fprintf(fp,"*************** Inc=%g ****************\n",Inc);
	//*******************************************************************
	 h_type HG=Histogram("pattern subLPRs",Min,Max,Inc);
	 // h_type HG=Histogram("pattern subLPRs",0,500,1.0);
	 Int4 round;
	 for(dH=dH_LPR,round=1; round <=2; dH=dH_BN,round++){
		   fprintf(fp,"              __Foreground___    __Background___     Information    %-15s WtNumSeq\n",model_str);
		   fprintf(fp,"    Pattern:  Match  Diverge     Match  Diverge       nats(rank)      %s\n",units_str);
		  Int4 *Order,Size=ItemsInHeap(dH);
		  NEW(Order,Size+3,Int4);
		  for(rank=1;!emptyHeap(dH); rank++){
		   j=delminHeap(dH); 
		   Order[rank]=j;
		   if(qst[j]){
		     char tmp[30],*ptr; ptr=tmp;
		     for(Int4 r=StartAlpha; r <= nAlpha(AB); r++){
			if(MemSset(r,qst[j])){ 
				sprintf(ptr,"%c",AlphaChar(r,AB)); ptr++;
			}
		     }
		     // if(fake2real) sprintf(ptr,"%d",fake2real[j]); 
		     if(fake2real) sprintf(ptr,"%d(%d)",fake2real[j],j); 
		     else sprintf(ptr,"%d",j);
		     fprintf(fp,"%12s: %5.0f %6.0f (%2.0f%c) %5.0f %6.0f (%2.0f%c)    %5.1f(%2d)",
			tmp,MatchFG[j],MisMatchFG[j],
			100.0*MatchFG[j]/(MatchFG[j]+MisMatchFG[j]),'%',
			MatchBG[j],MisMatchBG[j],
			100.0*MatchBG[j]/(MatchBG[j]+MisMatchBG[j]),'%',subLPR[j],rankLPR[j]);
		     if(round == 1) IncdHist(subLPR[j],HG);
		     fprintf(fp,"    %5.1f(%2d)", (double) keyBNPD[j], rankBNPD[j]);
#if 1	// DEBUG...
		     fprintf(fp,"      %5.1f", MatchFG[j]+MisMatchFG[j]+MatchBG[j]+MisMatchBG[j]);
#endif
		     {	// calculate binomial probability...
#if 0		// debug...
		        fprintf(fp," theta[%d]%c]=%f; SumXi=%.2f\n",
				j,AlphaChar(query[j],AB),
				theta[j][query[j]],SumXI[j]*wt_factor);
#else
			fprintf(fp,"\n");
#endif
		     }
		   } map+=subLPR[j];
		   if(rank % 10 == 0) fprintf(fp,"\n");
		  } fprintf(fp,"\n");
		  {
			fprintf(fp,"Pattern sets:\n");
		  	for(Int4 r=1; r <= Size; r++){
			   j=Order[r]; 
			   if(j!=0){
				PutSST(fp,qst[j],AB); 
				if(fake2real) fprintf(fp,"%d,",fake2real[j]);
				else fprintf(fp,"%d,",j);
			   }
			} fprintf(fp,"\n");
		  }
		  free(Order);
		  break;	// Skip round #2
	 }
	 fprintf(fp," total map = %.1f\n\n",map);
	 fprintf(fp," alpha=%.2f\n",alpha);
	 fprintf(fp," sigma=%.2f\n",sigma);
	 fprintf(fp," log(Rho) =%g; log(1-Rho)=%g\n",logRho,logOneMinusRho);
	 fprintf(fp," ModelMode = %c; ModelSwitch = %c\n",ModelMode,ModelSwitch);
	 fprintf(fp," Compare Mode = %c\n",mode);
	 PutHist(fp,60,HG);
	 NilHist(HG);
	 Nildheap(dH_LPR); Nildheap(dH_BN); 
	 free(rankLPR); free(rankBNPD);
	 free(model_str); free(units_str);
	 return ;
}


