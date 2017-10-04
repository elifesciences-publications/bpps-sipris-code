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
#include "blosum62.h"

bpps_typ::bpps_typ(sst_typ *Qst,unsigned char *Query, Int4 Query_len,
	double A0, double B0, a_type A, double **rho,double priorRi,BooLean UseGlobalSqWts)
// Original constructor...
{
	Int4	i,r;
	if(UseGlobalSqWts) StartAlpha=0; else StartAlpha=1;
	own_beta_query=FALSE; AB=A; RST=0; Rho=0;
	PriorRho=0; // this is used only in type-mode.
	NEW(beta,nAlpha(AB)+3,double);
	for(r=1; r <= nAlpha(AB); r++) beta[r]=WtFactor; // default = non-informed prior.
	NEWP(logRho,Query_len + 3, double); NEWP(Rho,Query_len + 3, double);
	for(i=1; i <= Query_len; i++){
	   NEW(logRho[i],nAlpha(AB) + 3, double); NEW(Rho[i],nAlpha(AB) + 3, double);
	} Init(Qst,Query,Query_len,A0,B0,rho,priorRi); 
	Type='G';
}

bpps_typ::bpps_typ(sst_typ *qsst, unsigned char *Qry, Int4 Query_len, char type, a_type A, 
	BooLean UseGlobalSqWts,double probRi)
{
	if(UseGlobalSqWts) StartAlpha=0; else StartAlpha=1;
	own_beta_query=TRUE; AB=A; 
	NEW(beta,nAlpha(AB)+3,double);
	for(Int4 r=1; r <= nAlpha(AB); r++) beta[r]=WtFactor; // default = non-informed prior.
	switch (type){
	    case 'r': case 'R': RST=new rst_typ('R',AB); if(probRi==0) probRi=0.5; break;
	    // case 'L': case 'l': RST=new rst_typ('L',AB); break;
	    case 'L': case 'l': RST=new rst_typ('G',AB); if(probRi==0) probRi=0.3; break;
	    case 'M': case 'm': RST=new rst_typ('G',AB); if(probRi==0) probRi=0.3; break;
	    default:  RST=new rst_typ('G',AB); if(probRi==0) probRi=0.3; break;
	}
#if 0	// DEBUG...
	assert(probRi <= 0.5 && probRi >= 0.3);
	probRi=0.5;
#endif
	switch (type){
	  case 'G': Initialize(qsst,90.0,10.0,probRi,0.00001,Query_len,Qry); break;
#if 1	// test parameters...       A0, B0,  Ri, rho
	  case 'L': Initialize(qsst,5.0,1.0,probRi,1e-4,Query_len,Qry); break;
	  case 'M': Initialize(qsst,4.0,1.0,probRi,2e-4,Query_len,Qry); break;
#if 1	// use for initial generation...need to allow for a lot of contamination.
	  case 'I': Initialize(qsst,6.0,3.0,probRi,0.001,Query_len,Qry); break;
#endif
	  // case 'M': Initialize(qsst,4.0,1.0,0.5,1e-2,Query_len,Qry); break;
	  // case 'R': Initialize(qsst,6.0,4.0,0.5,0.001,Query_len,Qry); break;
	  // case 'R': Initialize(qsst,6.0,4.0,0.5,0.001,Query_len,Qry); break;
	  case 'R': Initialize(qsst,6.0,4.0,probRi,0.001,Query_len,Qry); break;
	  case 'l': Initialize(qsst,48.0,2.0,probRi,0.0001,Query_len,Qry); break;
	  case 'm': Initialize(qsst,40.0,10.0,probRi,0.0001,Query_len,Qry); break;
	  case 'r': Initialize(qsst,42.0,8.0,probRi,0.0001,Query_len,Qry); break;
	  default: Initialize(qsst,48.0,2.0,probRi,0.0001,Query_len,Qry); break;
#else
	  case 'L': Initialize(qsst,48.0,2.0,0.3,0.003,Query_len,Qry); break;
	  case 'M': Initialize(qsst,40.0,10.0,0.3,0.003,Query_len,Qry); break;
	  case 'R': Initialize(qsst,60.0,40.0,0.5,0.0001,Query_len,Qry); break;
	  case 'l': Initialize(qsst,48.0,2.0,0.3,0.0001,Query_len,Qry); break;
	  case 'm': Initialize(qsst,40.0,10.0,0.3,0.0001,Query_len,Qry); break;
	  case 'r': Initialize(qsst,42.0,8.0,0.5,0.0001,Query_len,Qry); break;
	  default: Initialize(qsst,48.0,2.0,0.3,0.0001,Query_len,Qry); break;
#endif
	} 
	if(type=='R' || type == 'r') Type='R'; else Type='G';
}

void 	bpps_typ::Initialize(sst_typ *qsst, double A0, double B0,double priorRi, double rho,
					Int4 Query_len, unsigned char *Qry)
{
	sst_typ **LegalSST=RST->LegalResSets();
	unsigned char	*Query;
	Int4		i,x,s,r,Length=Query_len;

	PriorRho=rho;
	sst_typ **SST; NEWP(SST,Length+3,sst_typ);
	NEW(Query,Length+3,unsigned char);
	for(s=1; s <= Length; s++){
	   if(Qry){ 
	     r=Qry[s];
	     if(qsst[s]){
		if(!MemSset(r,qsst[s])){
		  char *str=GetPatternFromSST(qsst[s],AB); r=AlphaCode(str[0],AB); free(str);
		} 
	     } 
	   } else {
	     if(qsst[s]){
		char *str=GetPatternFromSST(qsst[s],AB); r=AlphaCode(str[0],AB); free(str);
	     } else r=0;
	   }
	   if(r==0) r=AlphaCode('A',AB);
	   Query[s]=r; SST[s]=LegalSST[Query[s]];
        }
	Rho=GetRhoCategoricalPriors(0, Length, rho, SST, AB); 
	NEWP(logRho,Query_len + 3, double); 
	for(i=1; i <= Length; i++) NEW(logRho[i],nAlpha(AB) + 3, double); 
	// for(s=1; s <= Length; s++){ free(SST[s]); } 
	free(SST);

	Init(qsst,Query,Length,A0,B0,0,priorRi);
	if(del_as_random) this->TreatDeletionAsRandom();
	own_beta_query=TRUE;
	// free(Query); // free pps after using...
}

void    bpps_typ::Init(sst_typ *Qst,unsigned char *Query, Int4 Query_len,
		double A0, double B0, double **rho,double priorRi)
{
	if(priorRi == 0.5){ PriorRi = 0.5; LogRatioPriorRi=0.0; }
	else {
		PriorRi=priorRi;
		assert(priorRi > 0.0 && priorRi < 1.0);
		LogRatioPriorRi = log(PriorRi/(1.0 - PriorRi));
		if(isnan(LogRatioPriorRi)) print_error("-Ri option input error");
		if(!isfinite(LogRatioPriorRi)) print_error("-Ri option input error");
		if(LogRatioPriorRi == HUGE_VAL){
		   // if(errno == EDOM || errno == ERANGE){ 
		   // fprintf(stderr,"LogRatioPriorRi = %lf;EDOM = %f; \n",LogRatioPriorRi);
		   // }
		   print_error("Fatal: PriorRi value is too extreme; exiting now...");
		}
	}
	del_as_random=0;	// 0 = ignored; 1 = random residue freqs; ? all non-functional.
	NullFreq=0; NullSST=0;
	PseudoCnts=0.0;
	DebugCheck=FALSE;
	query = Query; qst=Qst;
	AltMeasure= 'U';
	a0 = A0*WtFactor;	// = alpha prior.
	b0 = B0*WtFactor;	// = beta prior.
	alpha = a0/(a0+b0);
	b = beta; 
	k = Query_len;

	UpdateRho(rho);	// rho == 0 if called from initialize()

#if 0	// get stdfreq
	double sum=0.0;
	Int4 x;
	for(x=0; x <= nAlpha(AB); x++){
		stdfreq[x] = blosum62freq[x]; sum+=stdfreq[x];
	}
	for(x=0; x <= nAlpha(AB); x++) stdfreq[x] = stdfreq[x]/sum;
#endif
	wt_factor=1.0/(double)WtFactor;
	NEWP(theta,k+1,double); NEW(NullTheta,k+1,double); NEW(SumXI,k+1,double);
	NumCols=0;
	for(Int4 j=1;j<=k;j++){ 
		NEW(theta[j],nAlpha(AB)+2,double); 
		if(qst[j]){ NumCols++; assert(MemSset(query[j],qst[j])); }
	}
        NEW(argum,nAlpha(AB)+3,double);
	NEW(SubNullMap,k+3,double); NEW(SubMap,k+3,double);
	NEW(subLPR,k+3,double);
	NEW(subRawLPR,k+3,double);	// THIS APPEARS TO BE THE SAME AS SubMap array!!!
	NEW(MatchFG,k+3,double); NEW(MisMatchFG,k+3,double);
	NEW(MatchBG,k+3,double); NEW(MisMatchBG,k+3,double);
	nullMap=(double)INT4_MIN;
	ComputeNull=TRUE;
}

BooLean	bpps_typ::UpdateSST(sst_typ **&sst, e_type keyE,FILE *efp)
{
	Int4	j,r,ris;
	BooLean	Changed=FALSE;
	assert(LenSeq(keyE) == k && sst != 0);
        for(j=1; j <= k; j++){
	   if(sst[j]) free(sst[j]);
	   NEW(sst[j],k*(RST->MaxResSet()+3),sst_typ);;
           r = ResSeq(j,keyE);
	   for(ris=0; ris <= RST->NumResSets[r]; ris++){ sst[j][ris]=RST->ResidueSSet[r][ris];; }
	   sst[j][ris]=0;
           if(qst[j] && !MemSset(r,qst[j])){ // Fix incompatibility between qst and keyE 
		if(efp){ fprintf(stderr," Changed "); PutSST(stderr,qst[j],AB); }
		qst[j]=SsetLet(r);
		if(efp){ fprintf(stderr,"%d to ",j); PutSST(stderr,qst[j],AB); fprintf(stderr,"%d\n",j); }
		Changed=TRUE;   // qst is now different from the best.
           }
	} return Changed;
}

BooLean	bpps_typ::UpdateRho(sst_typ **&sst, e_type keyE)
{
	assert(PriorRho > 0.0);
	BooLean	Changed=this->UpdateSST(sst,keyE);
	for(Int4 j = 1; j <= k; j++){ free(Rho[j]); } free(Rho);
        Rho=GetRhoCategoricalPriors(0,k,PriorRho,sst,AB);
	this->UpdateRho( );
	// for(Int4 j=1; j <= k; j++) free(rho[j]); free(rho);
	return Changed;
}

void	bpps_typ::UpdateRho(double **rho)
{
	for(Int4 i=1; i <= k; i++){
	   double total=0.0;
	   // for(Int4 x=startAB; x <= nAlpha(AB); x++)
	   for(Int4 x=1; x <= nAlpha(AB); x++) {	// x is the size of the sequence set....
		if(rho) Rho[i][x] = rho[i][x];
		total += Rho[i][x];
		if(Rho[i][x] > 0.0) logRho[i][x] = log(Rho[i][x]); 
		else {
		   logRho[i][x] = -99999999.9;  // this should not be needed...sets shouldn't get this big...
		   // if it does, then it is likely that something happened to qst!
		}
	   }
	   if(total == 0.0){	// insertion here in query
	     if(rho)  assert(rho[i][0] == 1);
	     logRho[i][0] = 0.0; Rho[i][0] = 1.0;
	   } else {
	     assert(total > 0.0 && total < 1.0);
	     // logRho[i][0] = log(rho[i][0]); // this is not used..
	     logRho[i][0] = log(1.0-total); Rho[i][0] = 1.0-total; 	// these are used for null model...
	   }
	}
}

Int4    *bpps_typ::nBestPttrn(Int4 &N,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG)
{
         double map=0.0;
         map=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
         Int4 j,Min,Max,C;
         dh_type dH=dheap(k+2,4);
         for(j=1; j <= k; j++){
                   if(qst[j]){
		     double p = MatchFG[j]/(MatchFG[j]+MisMatchFG[j]);
		     double q = MatchBG[j]/(MatchBG[j]+MisMatchBG[j]);
		     if(p >= 0.85 && q <= 0.10){
                     	insrtHeap(j,(keytyp)-subLPR[j],dH);
		     }
                   }
         }
         Int4   rank;
        Int4 *Order,Size=ItemsInHeap(dH);
        NEW(Order,Size+3,Int4);
        for(rank=1;!emptyHeap(dH); rank++){
                   Order[rank]=delminHeap(dH);
        }
	// MatchFG[j]+MisMatchFG[j]+MatchBG[j]+MisMatchBG[j];
#if 1
        fprintf(stderr,"Pattern sets:\n");
        for(Int4 r=1; r <= Size; r++){
           j=Order[r];
           if(j!=0){
                PutSST(stderr,qst[j],AB);
                if(fake2real) fprintf(stderr,"%d,",fake2real[j]);
                else fprintf(stderr,"%d,",j);
           }
        } fprintf(stderr,"\n");
#endif
        Nildheap(dH);
        N=Size;
        return Order;
}

void	bpps_typ::UpdateRhoJ(Int4 j, double *rho)
{
	assert(j > 0 && j <= k);
	free(logRho[j]);free(Rho[j]); 
	NEW(logRho[j],nAlpha(AB) + 3, double); NEW(Rho[j],nAlpha(AB) + 3, double);
	double total=0.0;
	// for(Int4 x=startAB; x <= nAlpha(AB); x++)
	for(Int4 x=1; x <= nAlpha(AB); x++)
	{	// x is the size of the sequence set....
		total += rho[x];
		Rho[j][x] = rho[x];
		if(rho[x] > 0.0) logRho[j][x] = log(rho[x]); 
		else {
		   logRho[j][x] = -99999999.9;  // this should not be needed...sets shouldn't get this big...
		   // if it does, then it is likely that something happened to qst!
		}
	}
	if(total == 0.0){	// insertion here in query
	     assert(rho[0] == 1); logRho[j][0] = 0.0; Rho[j][0] = 1.0;
	} else {
	     assert(total > 0.0 && total < 1.0);
	     logRho[j][0] = log(1.0-total); 	// this is used for null model...
	     Rho[j][0] = 1.0-total; 	// this is used for null model...
	}
}

Int4	bpps_typ::Compare(FILE *fp,bpps_typ *that)
// compare this with xpps;
{
	Int4	i,x,n=0,ok=1;
	fprintf(fp,"Comparing two bpps_typ objects...");
	if(this->k != that->k){ fprintf(fp,"lengths differ (%d vs %d)\n",k,that->k); n++; ok=0; }
	if(this->a0 != that->a0){ fprintf(fp,"a0 differs (%.1f vs %.1f)\n",a0,that->a0); n++; }
	if(this->b0 != that->b0){ fprintf(fp,"b0 differs (%.1f vs %.1f)\n",b0,that->b0); n++; }
	if(this->PriorRi != that->PriorRi){
		 fprintf(fp,"PriorRi differs (%.1f vs %.1f)\n", PriorRi,that->PriorRi); n++;
	}
	if(ok){
	 for(Int4 i=1; i <= k; i++){
	   if(qst[i] != that->qst[i]){ 
		fprintf(fp,"\nqsst[%d] differ (",i);
		PutSST(fp,qst[i],AB); fprintf(fp," vs "); PutSST(fp,that->qst[i],AB); n++;
	   }
	   for(Int4 x=startAB; x <= nAlpha(AB); x++){	// x is the size of the sequence set....
		if(Rho[i][x] != that->Rho[i][x]){
		   fprintf(fp,"\nRho[%d] differ:" ,i); PutRho(fp,i); fprintf(fp," vs "); 
		   that->PutRho(fp,i);  n++; break;
		}
	   }
	 }
	}
	if(n > 0) fprintf(fp,"\nSummary: %d differences\n",n); else fprintf(fp,"checks out okay.\n");
	return n;
}

void	bpps_typ::PutRho(FILE *fp, Int4 i)
{
	assert(i > 0 && i <= k);
	fprintf(fp," ",i);
	// for(Int4 x=startAB; x <= nAlpha(AB); x++)
	for(Int4 x=1; x <= nAlpha(AB); x++)
	{	// x is the size of the sequence set....
		if(Rho[i][x] <= 0.0) break;
		fprintf(fp,"%d=%.3f " ,x,log(Rho[i][x]));
	} // fprintf(fp,"\n");
}

void	bpps_typ::PutRho(FILE *fp)
{
	for(Int4 i=1; i <= k; i++){
	   fprintf(fp,"%d: ",i);
	   // for(Int4 x=startAB; x <= nAlpha(AB); x++)
	   for(Int4 x=1; x <= nAlpha(AB); x++)
	   {	// x is the size of the sequence set....
		if(Rho[i][x] == 0.0) break;
		fprintf(fp,"%d %.3f" ,x,Rho[i][x]);
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
}

double	bpps_typ::ProbCols( )
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

void bpps_typ::Free()
{
	for(Int4 j = 1; j <= k; j++){ free(theta[j]); free(logRho[j]);free(Rho[j]); }
	free(logRho); free(theta); free(Rho);
	free(NullTheta); free(SumXI); free(argum);
	free(SubNullMap); free(SubMap); free(subLPR); free(subRawLPR);
	free(MatchFG); free(MisMatchFG); free(MatchBG); free(MisMatchBG);
	free(beta);
	if(RST) delete RST; 
	if(own_beta_query) { free(query); }   // otherwise  beta and Query not freed with pps.
}

//************************** Binary ************************************
#include "blosum62.h"

Int4    bpps_typ::CntDelete(double &matFG, double &matBG, double &misFG, double &misBG,
		sst_typ xsst, UInt4 DelBG,UInt4 DelFG)
{
	double	frqM=0.0,frqD=0.0;
	for(char x=1; x <= nAlpha(AB); x++){
		if(MemSset(x,xsst)) frqM += blosum62freq[x]; else frqD += blosum62freq[x];
	}
#if 1	// treat Root node different than other nodes.
	if(Type == 'R'){	// == Root node...
		matBG += frqM*DelBG; misBG += frqD*DelBG;
	} else {
		matFG += frqM*DelFG; misFG += frqD*DelFG;	
		matBG += frqM*DelBG; misBG += frqD*DelBG;
	}
#elif 1	// treat deletions as random background in both FG and BG...
		matFG += frqM*DelFG; misFG += frqD*DelFG;	
		matBG += frqM*DelBG; misBG += frqD*DelBG;
#elif 1 // count deletions as random background in BG only!
		matBG += frqM*DelBG; misBG += frqD*DelBG;
#else	// treat deletions as a special mismatch residue!
	misFG += DelFG;	
	misBG += DelBG;
#endif
}

void bpps_typ::BinaryUpdateSumXI(UInt4 **CountFG,UInt4 **CntBG)
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
#if 0
	   matBG=1.0/wt_factor;		// add prior == 100.0; m1
	   sumBG=2.0/wt_factor;		// add prior *2; n2 + m2.
#else
	   matBG=PseudoCnts/wt_factor;		// add prior == 100.0; m1
	   sumBG=2.0*PseudoCnts/wt_factor;	// add prior *2; n2 + m2.
#endif
#if 0	// treat deletions like random background...
	   if(del_as_random){
		dummy=Dummy=0.0;
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(matFG,matBG,dummy,Dummy,qst[j],CntBG[j][0],CountFG[j][0]);
		sumBG += CntBG[j][0];
	   }
#endif
	   for(x=startAB; x <= nAlpha(AB); x++){
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

#define UseJunsPatternWeighting 0
void bpps_typ::BinaryUpdateTheta(UInt4 **CntBG,UInt4 **CntFG)
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
#if 0
	        double mat =1.0/wt_factor;		// add prior == 100.0
	        double sum=1.0/wt_factor;		// add prior == 100.0
#else
	        double mat=PseudoCnts/wt_factor;		// add prior == 100.0
	        double sum=PseudoCnts/wt_factor;		// add prior == 100.0
#endif
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(mat,mat,sum,sum,qst[j],CntBG[j][0],CntFG[j][0]);
	   }
#endif
                for(x=startAB;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])){
			mat += CntBG[j][x] + CntFG[j][x];
		   } else {
			sum += CntBG[j][x] + CntFG[j][x];
#if UseJunsPatternWeighting
			theta[j][x] = CntBG[j][x] + CntFG[j][x] + b[x];
#endif
		   }
                }
		sum += mat;	// Failed to sum this in previous loop.
		// error here?! should not include in sum..
		assert(mat > SumXI[j]);
		mat -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
		sum -= SumXI[j];	// subtract out FG non-contamination, which doesn't count.
#if UseJunsPatternWeighting
                for(x=startAB;x<=nAlpha(AB);x++) {
		   if(!MemSset(x,qst[j])) theta[j][x] = theta[j][x]/sum;
		}
#endif
//              theta[j] = sample_from_dirichlet(argum);
		theta[j][q_j] = mat/sum;
#if 0
		if(DebugCheck && j == 69){ fprintf(stderr,"theta[%d][%c] = %f; mat = %f; sum = %.2f\n",
			j,AlphaChar(q_j,AB),theta[j][q_j],mat,sum); }
#endif
        } 
}

//************************** Binary ************************************

void bpps_typ::updateAlpha(UInt4 **CountFG)
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
	      for(x=startAB; x <= nAlpha(AB); x++){
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
	if(!(alpha <= 1.0 && alpha >= 0.0)){ 
		fprintf(stderr,"alpha = %g; S = %g; N0 = %g\n",alpha,S,N0);
		// fprintf(stderr,"alpha = %g\n",alpha);
		assert(alpha <= 1.0 && alpha >= 0.0);
	}
}

#if 0	// not used right now..
double	bpps_typ::ProbRatio(unsigned char *seq, unsigned char **SqWt,
		UInt4 **CntBG, UInt4 **CntFG)
// see Gibbs sampling algorithm in section 4 of Jun's JASA paper.
// Step 1.
// Sample sequence i into either the query or the non-query partition...
{
	double		prob;
	Int4		j;
	unsigned char	q_j;
	
	alpha=a0/(a0+b0);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
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
#endif

//************************* Collapsed MAPs Below *************************

#if 0	// old version.
void    bpps_typ::BinaryUpdatePseudoCnts(UInt4 **CntFG,UInt4 **CntBG)
{
	Int4	j,r;
	double	N=0.0;
	for(j=1; j<=k; j++) {
	   for(r=startAB; r<= nAlpha(AB); r++) { N += CntFG[j][r] + CntBG[j][r]; }
	} N = N/(double) k;	// average N;
	PseudoCnts=MAXIMUM(double,1.0,sqrt(N*wt_factor)/10.0);	// take root.
}
#elif 0 // new, faster version.

void    bpps_typ::BinaryUpdatePseudoCnts(UInt4 **CntFG,UInt4 **CntBG)
{
	register Int4	j,r;
	register UInt4	*cntFGj,*cntBGj;
	register double	N=0.0;
	for(j=k; j>=1; j--) {
	   cntFGj=CntFG[j]; cntBGj=CntBG[j];
	   for(r=nAlpha(AB); r >=startAB ; r--) { N += cntFGj[r] + cntBGj[r]; }
	} N = N/(double) k;	// average N;
	PseudoCnts=MAXIMUM(double,1.0,sqrt(N*wt_factor)/10.0);	// take root.
}
#else	// newer, even faster version.

void    bpps_typ::BinaryUpdatePseudoCnts(UInt4 **CntFG,UInt4 **CntBG)
{
	register Int4	r;
	register UInt4	*cntFGj,*cntBGj;
	register double	N=0.0;
	cntFGj=CntFG[1]; cntBGj=CntBG[1];
	for(r=nAlpha(AB); r >= 0; r--) { N += cntFGj[r] + cntBGj[r]; }
	PseudoCnts=MAXIMUM(double,1.0,sqrt(N*wt_factor)/10.0);	// take root.
}

#endif

double	bpps_typ::MapBiNom(UInt4 **CntBG, UInt4 **CntFG)
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
	double	AveNumWtSeqFG=0.0;	// NEW for PriorRi.
	Int4	NumPttrnPos=0;		// NEW for PriorRi.

	alpha=a0/(a0+b0);
	BinaryUpdatePseudoCnts(CntFG,CntBG);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
	   	// if(!qst[j]){ map += log(1.0 - rho); continue; }
	   	if(!qst[j]){ continue; } // cancels out with null map...
		else NumPttrnPos++;	// NEW for PriorRi.
		q_j=query[j];
		theta_j = theta[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1,n1,m2,n2;
		m1=n1=m2=n2=0.0;
#if 1	// Pattern residue set optimization.
		double	NumResFG[30];
		Int4	NumTypFG=0;
#endif
		Int4 CardPttrnSet=0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(n1,m1,n2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			CardPttrnSet++; m1 += CountBG_j[t]; n1 += CountFG_j[t];
			NumTypFG++; NumResFG[NumTypFG]=CountFG_j[t]*wt_factor; // Pattern set optimization.
		   } else {
			m2 += CountBG_j[t]; n2 += CountFG_j[t];
		   }
		}
		n1 = wt_factor*n1; m1 = wt_factor*m1;
	 	n2 = wt_factor*n2; m2 = wt_factor*m2;
		map += n1*log((1-alpha)*theta_j[q_j] + alpha) + m1*log(theta_j[q_j]); 
		map += n2*log(1-alpha) + (n2+m2)*log(1.0-theta_j[q_j]);
		assert(CardPttrnSet > 0);
		map += logRho[j][CardPttrnSet];	// column penalty...times # residues
		AveNumWtSeqFG += n1 + n2;	// New for PriorRi.
	}
	if(NumPttrnPos > 0 && LogRatioPriorRi != 0.0){	// apply priorRi in LPR computation.
		AveNumWtSeqFG = AveNumWtSeqFG/NumPttrnPos;
		map += AveNumWtSeqFG * LogRatioPriorRi; // LPR in FG versus NullLPR in BG 
	} return map;
}

double bpps_typ::NullMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
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
		double sum=PseudoCnts/wt_factor;
		double match=PseudoCnts/wt_factor;
#if 1	// treat deletions like random background...
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
	for(map=0.0,j=1; j<=k; j++) {
		// use same pattern but put all sequences into the background.
	   	if(!qst[j]) continue;	// these columns cancel out for Net-LPR
		q_j=query[j];
		double theta_x = NullTheta[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		double m1=0.0,m2=0.0;
		Int4 CardPttrnSet=0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		dummy=Dummy=0.0;
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(m1,m1,m2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){ m1+= CountBG_j[t] + CountFG_j[t]; CardPttrnSet++; }
		   else { m2 += CountBG_j[t] + CountFG_j[t]; }
		} 
		m1=wt_factor*m1; m2=wt_factor*m2;
		map += (m1*log(theta_x)) + (m2*log((1.0-theta_x)));
		assert(CardPttrnSet > 0);
		map += logRho[j][0];	// columns out for null model.
	} nullMap=map; ComputeNull=FALSE;
    } return nullMap;
}

#if UseJunsPatternWeighting
double	bpps_typ::JunsPatternWeighting(Int4 nAf, double theta_f, double *theta_j, double *Nb, double *Nf)
{ 
	double	logP=0.0;
	double	Nb_sum=0.0,Nf_sum=0.0;
	double	ave_theta=theta_f/(double)nAf;

	print_error("Not tested and not completed; this is bogus.");
        for(x=startAB;x<=nAlpha(AB);x++) {
	   logP += Nb[x] * log(theta[x]);
	   if(MemSset(x,qst[j])){ 		// x element of Af.
		logP += (Nf[x])*log((alpha*ave_theta/ ( (1.0-alpha)*(theta_f/(double)nAf) ) ) );
	   } else {				// x not element of Af.
		logP +=Nf[x]*log(1.0-alpha);
	   }
        }
	return LogP;
}
#endif

#if 0	// old code...
double bpps_typ::subMapBiNom2(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
{
	return subMapBiNom(CntBG, CntFG);	// testing...
	double		submap,map,*theta_j;
	Int4		j;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	double	AveNumWtSeqFG=0.0;	// NEW for PriorRi.
	Int4	NumPttrnPos=0;		// NEW for PriorRi.

	alpha=a0/(a0+b0);
	BinaryUpdatePseudoCnts(CntFG,CntBG);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
	   	if(!qst[j]){ SubMap[j]=0; continue; }
		else NumPttrnPos++;	// NEW for PriorRi.
		theta_j = theta[j];
		q_j=query[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1,n1,m2,n2;
		m1=n1=m2=n2=0.0;
#if 1	// Pattern residue set optimization.
		double	NumResFG[30];
		Int4	NumTypFG=0;
		double	NumResBG[30];
		Int4	NumTypBG=0;
#endif
		Int4	CardPttrnSet=0;
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			CardPttrnSet++;
			n1+= CountFG_j[t]; m1+=CountBG_j[t];
			NumTypFG++; NumResFG[NumTypFG]=CountFG_j[t]*wt_factor; // Pattern set optimization.
		   } else {
			n2+= CountFG_j[t]; m2+=CountBG_j[t];
		   }
		}
		submap=0.0;
		// n1 = matching FG; m1 = matching BG
		// n2 = mis-match FG; m2 = mismatch BG.
		m1=m1*wt_factor; m2=m2*wt_factor;
		n1=n1*wt_factor; n2=n2*wt_factor;
		MatchFG[j] = n1; MisMatchFG[j] = n2;
		MatchBG[j] = m1; MisMatchBG[j] = m2;
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
		assert(CardPttrnSet > 0);
		submap += logRho[j][CardPttrnSet];
		SubMap[j] = submap;
		map += submap;
		AveNumWtSeqFG += n1 + n2;	// New for PriorRi.

	}
	if(NumPttrnPos > 0 && LogRatioPriorRi != 0.0){	// apply priorRi in LPR computation.
		AveNumWtSeqFG = AveNumWtSeqFG/NumPttrnPos;
		map += AveNumWtSeqFG * LogRatioPriorRi; // LPR in FG versus NullLPR in BG 
		if(isnan(map)){
			fprintf(stderr,"AveNumWtSeqFG=%g; NumPttrnPos=%g; LogRatioPriorRi=%g\n",
				AveNumWtSeqFG,NumPttrnPos,LogRatioPriorRi);
			abort();
		}
	} return map;
}

double bpps_typ::subNullMapBiNom2(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
// Considers the query residue set at each position to be null.
{
	return subNullMapBiNom(CntBG, CntFG);
	double		submap,map;
	Int4		j,t;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	q_j,x;

        for(j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		double sum=PseudoCnts*wt_factor;
		double match=PseudoCnts*wt_factor;
                for(x=startAB;x<=nAlpha(AB);x++) {
		   if(MemSset(x,qst[j])) match += CntBG[j][x] + CntFG[j][x];
		   else sum += CntBG[j][x] + CntFG[j][x];
                } sum += match;
		NullTheta[j] = match/sum;
        } 
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]){ SubNullMap[j]=0; continue; }
		q_j=query[j];
		double theta_x = NullTheta[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		double m1=0.0,m2=0.0;
		Int4 CardPttrnSet=0;
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){ CardPttrnSet++; m1+= CountBG_j[t] + CountFG_j[t]; }
		   else { m2+= CountBG_j[t] + CountFG_j[t]; }
		}
		submap=0.0;
		m1 = m1*wt_factor; m2 = m2*wt_factor;
		submap += (m1*log(theta_x) + m2*log((1.0-theta_x)));
		assert(CardPttrnSet > 0);
		submap += logRho[j][0];	// log of one minus other Rho values..
		map += submap;
		SubNullMap[j]=submap;
	} return map;
}
#endif

double	*bpps_typ::SubRawLPR(UInt4 **CntBG,UInt4 **CntFG)
{
	double	map0,map;
	map=subMap(CntBG,CntFG);
	subRawLPR[0]=map;
	map0=0.0;
	for(Int4 j=1; j <= k; j++){
	   if(qst[j]){
	     subRawLPR[j]=SubMap[j];
	     map0+=subRawLPR[j];
	   } else subRawLPR[j]=0.0;
	} return subRawLPR;
}

double	*bpps_typ::SubLPR(UInt4 **CntBG,UInt4 **CntFG)
{
	double	map0,map;
	map=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	subLPR[0]=map;
	map0=0.0;
	for(Int4 j=1; j <= k; j++){
	   if(qst[j]){
	     subLPR[j]=SubMap[j]-SubNullMap[j];
	     map0+=subLPR[j];
	   } else subLPR[j]=0.0;
	}
	return subLPR;
}

void	bpps_typ::PutPattern(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG)
{
	 double	map=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	 Int4 j,Number=0;
	 dh_type dH=dheap(k+2,4);
	 for(j=1; j <= k; j++){ if(qst[j]) insrtHeap(j,(keytyp)-subLPR[j],dH); }
	 while(!emptyHeap(dH)){
	   j=delminHeap(dH); 
	   assert(qst[j] != 0);
	   char tmp[30],*ptr; ptr=tmp;
	   for(Int4 r=StartAlpha; r <= nAlpha(AB); r++){
		if(MemSset(r,qst[j])){ 
			sprintf(ptr,"%c",AlphaChar(r,AB)); ptr++;
		}
	   }
	   if(fake2real) sprintf(ptr,"%d(%d)",fake2real[j],j); 
	   else sprintf(ptr,"%d",j);
	   if(Number==0) fprintf(fp,"%s",tmp); else fprintf(fp,",%s",tmp);
	   Number++;
	} fprintf(fp,"\n");
	Nildheap(dH);
	return ;
}

void	bpps_typ::PutCDTreeSubLPR(FILE *fp,char *id, UInt4 **CntBG,UInt4 **CntFG)
#if 0	// Output for CDTree patterns.
cd00142: R135(89:0)=537.6; DE43(100:7)=525.6; YF116(94:5)=470.7; G147(94:4)=458.0; QK(81:1))35=455.9.
cd00891: T143=134.3; KR184=127.4; C136=125.7; G138=116.2; DE56=112.5.
cd00896: P177=59.4; P172=53.2; L148=50.6; D171=49.7; KR180=46.5.
#endif
{
	 SubLPR(CntBG,CntFG);
	 dh_type dH=dheap(k+2,4);
	 Int4 j,rank; 
	 double	map=0.0,TotalLPR=0.0;
	 for(j=1; j <= k; j++){
		if(qst[j]){
		   TotalLPR += subLPR[j];
		   insrtHeap(j,(keytyp)-subLPR[j],dH);
		}
	 } fprintf(fp,"%s: ",id);
	 for(rank=1;!emptyHeap(dH); rank++){
	   j=delminHeap(dH); 
	   if(qst[j]){
	     Int4 RunningPercent=(Int4)(floor(100.0*(subLPR[j]+map)/TotalLPR));
	     Int4 Percent=(Int4)(floor(100.0*subLPR[j]/TotalLPR));
	     Int4 Permille=(Int4)(floor(1000.0*subLPR[j]/TotalLPR));
	     for(Int4 r=StartAlpha; r <= nAlpha(AB); r++){
		if(MemSset(r,qst[j])){ fprintf(fp,"%c",AlphaChar(r,AB)); }
	     } fprintf(fp,"%d",j);
	     fprintf(fp,"(%.0f:%.0f)=%.1f",
		100.0*MatchFG[j]/(MatchFG[j]+MisMatchFG[j]),
		100.0*MatchBG[j]/(MatchBG[j]+MisMatchBG[j]),subLPR[j]);
	     // fprintf(fp,"(%d)",Percent);
	     fprintf(fp,"(%d)",Permille);
	     if(Percent < 1){ fprintf(fp,"."); break; }
	     if(RunningPercent >= 98){ fprintf(fp,"."); break; }
	     if(subLPR[j] < 3.0){ fprintf(fp,"."); break; }
	     if(emptyHeap(dH)) fprintf(fp,"."); else fprintf(fp,"; ");
	   } map+=subLPR[j];
	 } fprintf(fp,"\n\n");
	 Nildheap(dH); 
	 return ;
}


void	bpps_typ::PutSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG)
{
	 double	map=0.0,map0=0.0,sub,TotalLPR;
	 // DebugCheck = TRUE;
#if 0
	 map0=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
#else
	 SubLPR(CntBG,CntFG);
	 map0=subLPR[0];
#endif
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
		     case 'B': 		// LnBinomialDistribution...
			model_str=NewString("LnBinomial");
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
	 TotalLPR=0.0;
	 for(j=1; j <= k; j++){
		   if(qst[j]){
		     TotalLPR += subLPR[j];
		     insrtHeap(j,(keytyp)-subLPR[j],dH_LPR);
		     double N,M,p,q,n1,n2,m1,m2,P;
		     n1=MatchFG[j]; n2=MisMatchFG[j];
		     N=n1+n2; q=(n1+1)/(N+2);

		     m1 = MatchBG[j]; m2 = MisMatchBG[j];
		     M=m1+m2; p=(m1+1)/(M+2);

		     switch (mode){
		     case 'B': 		P= LnBinomialDistribution(N,p,n1); break;
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
		     default: print_error("bpps_typ input error"); break;
		     }
		     insrtHeap(j,(keytyp)-P,dH_BN);
		     min=MINIMUM(double,subLPR[j],min);
		     max=MAXIMUM(double,subLPR[j],max);
		   }
	 }
	 Int4	*rankLPR,*rankBNPD,rank;
	 double *keyBNPD;
	 keytyp	key;
	 NEW(rankLPR,k+3,Int4); NEW(rankBNPD,k+3,Int4); NEW(keyBNPD,k+3,double);
	 dH=dheap(k+2,4);
	 for(rank=1;!emptyHeap(dH_LPR); rank++){
		   key=minkeyHeap(dH_LPR); j=delminHeap(dH_LPR); 
		   rankLPR[j]=rank; insrtHeap(j,key,dH);
	 } Nildheap(dH_LPR); dH_LPR=dH;
	 dH=dheap(k+2,4);
	 for(rank=1;!emptyHeap(dH_BN); rank++){
		   key=minkeyHeap(dH_BN); j=delminHeap(dH_BN); 
		   rankBNPD[j]=rank; keyBNPD[j]=-key; insrtHeap(j,key,dH);
	 } Nildheap(dH_BN);  dH_BN=dH;
	 Min=(Int4) floor(min); Max=(Int4) ceil(max); Inc = ceil((max-min)/20.0);
	 if(Inc < 1.0) Inc=0.5;
	 if(Min == Max){ Max = Min + 1; }
	 // fprintf(fp,"*************** Inc=%g ****************\n",Inc);
	//*******************************************************************
	 if(Max <= Min) Max = Min + 4*Inc;
	 h_type HG=Histogram("pattern subLPRs",Min,Max,Inc);
	 // h_type HG=Histogram("pattern subLPRs",0,500,1.0);
	 Int4 round;
	 for(dH=dH_LPR,round=1; round <=2; dH=dH_BN,round++){
		   fprintf(fp,"              __Foreground___    __Background___     Information    %-15s %cinfo(sum )  WtNumSeq\n",model_str,'%');
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
#if 0	// output subMap values...
			fprintf(fp,"[%5.0f %6.0f] ", SubMap[j],SubNullMap[j]);
#endif
		     if(round == 1) IncdHist(subLPR[j],HG);
		     fprintf(fp,"    %5.1f(%2d)", (double) keyBNPD[j], rankBNPD[j]);
#if 0	// Add % of total map for each position.
		     Int4 Percent=(Int4)(floor(100.0*subLPR[j]/subLPR[0]));
		     Int4 RunningPercent=(Int4)(floor(100.0*(subLPR[j]+map)/subLPR[0]));
#else
		     Int4 Percent=(Int4)(floor(100.0*subLPR[j]/TotalLPR));
		     Int4 RunningPercent=(Int4)(floor(100.0*(subLPR[j]+map)/TotalLPR));
		     fprintf(fp,"    %2d%c(%3d%c)", Percent,'%',RunningPercent,'%');
#endif
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
		  	if(fake2real){	// print out column position as well...
			 for(Int4 r=1; r <= Size; r++){
			   j=Order[r]; 
			   if(j!=0){ PutSST(fp,qst[j],AB); fprintf(fp,"%d,",j); }
			 } fprintf(fp,"\n");
			}
		  }
		  free(Order);
		  break;	// Skip round #2
	 }
	 // fprintf(fp," total map = %.1f\n\n",map);
	 fprintf(fp," total LPR = %.1f\n\n",subLPR[0]);
	 fprintf(fp," alpha=%.5f\n",alpha);
	 fprintf(fp,"chn_see xxx -alpha=%.0f:%.0f=%.5f\n",a0*wt_factor,b0*wt_factor,alpha);
	 // fprintf(fp," log(Rho) =%g; log(1-Rho)=%g\n",logRho[j][1],logRho[j][0]); not valid...
	 fprintf(fp," Compare Mode = %c\n",mode);
	 PutHist(fp,60,HG);
	 NilHist(HG);
	 Nildheap(dH_LPR); Nildheap(dH_BN); 
	 free(rankLPR); free(rankBNPD);
	 free(keyBNPD);
#if 0
	 free(model_str); free(units_str);
#else	// NewString() uses new operation.
	 delete [] model_str; delete [] units_str;
#endif
	 return ;
}

void	bpps_typ::WriteSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG)
// for saving information on pattern residues and 
// read this in using csp->ReadBinomialBPPS(-cutoff,res_evals);
{
	 double	map=0.0,map0=0.0,sub;
	 // DebugCheck = TRUE;
	 map0=subMap(CntBG,CntFG) - subNullMap(CntBG,CntFG);
	 // DebugCheck = FALSE;
	// this sets MatchFG, etc...
	 Int4 j,C;
	//********************* Determine range of histogram *************
	 double min=DBL_MAX,max=-9999999.0,Inc;
	 dh_type dH_LPR=dheap(k+2,4),dH=0;
	 char mode=AltMeasure;
	 for(j=1; j <= k; j++){
		   if(qst[j]){ insrtHeap(j,(keytyp)-subLPR[j],dH_LPR); }
	 }
	 Int4	*rankLPR,rank;
	 keytyp	key;
	 NEW(rankLPR,k+3,Int4); 
	 dH=dheap(k+2,4);
	 for(rank=1;!emptyHeap(dH_LPR); rank++){
		   key=minkeyHeap(dH_LPR); j=delminHeap(dH_LPR); 
		   rankLPR[j]=rank; insrtHeap(j,key,dH);
	 } Nildheap(dH_LPR); 
	 fprintf(fp,"#             __Foreground___    __Background___     Information    %cinfo(sum )  WtNumSeq\n",'%');
	 fprintf(fp,"#   Pattern:  Match  Diverge     Match  Diverge       nats(rank)\n");
	 Int4 *Order,Size=ItemsInHeap(dH);
	 NEW(Order,Size+3,Int4);
	 for(rank=1;!emptyHeap(dH); rank++){
	   j=delminHeap(dH); Order[rank]=j;
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
#if 1	// Add % of total map for each position.
	     Int4 Percent=(Int4)(floor(100.0*subLPR[j]/subLPR[0]));
	     Int4 RunningPercent=(Int4)(floor(100.0*(subLPR[j]+map)/subLPR[0]));
	     fprintf(fp,"    %2d%c(%3d%c)", Percent,'%',RunningPercent,'%');
#endif
	     fprintf(fp,"      %5.1f\n", MatchFG[j]+MisMatchFG[j]+MatchBG[j]+MisMatchBG[j]);
	   } map+=subLPR[j];
	 } fprintf(fp,"\n");
	 free(Order);
	 fprintf(fp,"alpha=%.0f:%.0f=%.5f\n\n",a0*wt_factor,b0*wt_factor,alpha);
	 Nildheap(dH); 
	 free(rankLPR); 
	 return ;
}

double bpps_typ::subMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
{
	double		submap,map,*theta_j;
	Int4		j;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	t,q_j;

	double	AveNumWtSeqFG=0.0;	// NEW for PriorRi.
	Int4	NumPttrnPos=0;		// NEW for PriorRi.

	alpha=a0/(a0+b0);
	BinaryUpdatePseudoCnts(CntFG,CntBG);
        BinaryUpdateSumXI(CntFG,CntBG);
        updateAlpha(CntFG);
	BinaryUpdateTheta(CntBG,CntFG);
	for(map=0.0,j=1; j<=k; j++) {
		MatchBG[j]=MatchFG[j]=MisMatchBG[j]=MisMatchFG[j]=0.0;
	   	if(!qst[j]){ SubMap[j]=0; continue; }
		else NumPttrnPos++;	// NEW for PriorRi.
		theta_j = theta[j];
		q_j=query[j];
		CountBG_j = CntBG[j];
		CountFG_j = CntFG[j];
		double m1,n1,m2,n2;
		m1=n1=m2=n2=0.0;
		Int4	CardPttrnSet=0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(n1,m1,n2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			CardPttrnSet++;
			n1+= CountFG_j[t]; m1+=CountBG_j[t];
		   } else { n2+= CountFG_j[t]; m2+=CountBG_j[t]; }
		}
		submap=0.0;
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
		SubMap[j] = submap;
		map += submap;
		AveNumWtSeqFG += n1 + n2;	// New for PriorRi.

	}
	if(NumPttrnPos > 0 && LogRatioPriorRi != 0.0){
		AveNumWtSeqFG = AveNumWtSeqFG/NumPttrnPos;
		map += AveNumWtSeqFG * LogRatioPriorRi; // LPR in FG versus NullLPR in BG 
		if(isnan(map)){
			fprintf(stderr,"AveNumWtSeqFG=%g; NumPttrnPos=%g; LogRatioPriorRi=%g\n",
				AveNumWtSeqFG,NumPttrnPos,LogRatioPriorRi);
			abort();
		}
	}
	return map;
}

double bpps_typ::subNullMapBiNom(UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
// Considers the query residue set at each position to be null.
{
	double		submap,map;
	Int4		j,t;
	UInt4	*CountBG_j,*CountFG_j;
	unsigned char	q_j,x;

        for(j=1;j<=k;j++) {
	   	if(!qst[j]) continue;
		q_j = query[j];
		double sum=PseudoCnts*wt_factor;
		double match=PseudoCnts*wt_factor;
#if 1	// treat deletions like random background...
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
	for(map=0.0,j=1; j<=k; j++) {
	   	if(!qst[j]){ SubNullMap[j]=0; continue; }
		q_j=query[j];
		double theta_x = NullTheta[j];
		CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
		double m1=0.0,m2=0.0;
		Int4 CardPttrnSet=0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(m1,m1,m2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
		for(t=startAB; t<= nAlpha(AB); t++) {
		   if(MemSset(t,qst[j])){
			CardPttrnSet++; m1+= CountBG_j[t] + CountFG_j[t];
		   } else { m2+= CountBG_j[t] + CountFG_j[t]; }
		}
		submap=0.0;
		m1 = m1*wt_factor; m2 = m2*wt_factor;
		submap += m1*log(theta_x) + m2*log((1.0-theta_x));
		assert(CardPttrnSet > 0);
		submap += logRho[j][0];	// log of one minus other Rho values..
		map += submap;
		SubNullMap[j]=submap;
	} return map;
}

double bpps_typ::QuickLPRsubJ(Int4 j, UInt4 **CntBG, UInt4 **CntFG)
/************************************************************************
 ************************************************************************/
{
	double		sum,submap,*theta_j,m1,n1,m2,n2,theta_x;
	UInt4		*CountBG_j,*CountFG_j;
	unsigned char	t,q_j,x;

        assert(j > 0 && j<=k);
	if(!qst[j]) return 0.0;

	CountBG_j = CntBG[j]; CountFG_j = CntFG[j];
	alpha=a0/(a0+b0);
#if 1
	BinaryUpdatePseudoCnts(CntFG,CntBG); // doesn't change for column sampling
#else
	if(PseudoCnts == 0.0) BinaryUpdatePseudoCnts(CntFG,CntBG); // doesn't change for column sampling
#endif
#if 1
        BinaryUpdateSumXI(CntFG,CntBG); // *************
#else
        double matFG=0.0,matBG=0.0,sumBG=0.0;
        matBG=PseudoCnts/wt_factor;          // add prior == 100.0; m1
        sumBG=2.0*PseudoCnts/wt_factor;      // add prior *2; n2 + m2.
#if 1	// treat deletions like random background...
	   if(del_as_random){
		dummy=Dummy=0.0;
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(matFG,matBG,dummy,Dummy,qst[j],CountBG_j[0],CountFG_j[0]);
		sumBG += CountBG_j[0];
	   }
#endif
        for(x=startAB; x <= nAlpha(AB); x++){
                if(MemSset(x,qst[j])){
                        matFG += CountFG_j[x];
                        matBG += CountBG_j[x];
                } sumBG += CountBG_j[x];
        } theta_x=matBG/sumBG;
        SumXI[j] = matFG*(alpha/(alpha+((1-alpha)*theta_x))); // from Equation [3].
        if(SumXI[j] != 0 && matFG != 0) assert(SumXI[j] < matFG);
#endif
        updateAlpha(CntFG);
#if 1
	BinaryUpdateTheta(CntBG,CntFG); // ****************
#else 	//*************** BinaryUpdateTheta(CntBG,CntFG); // ****************
        double mat=PseudoCnts/wt_factor;                // add prior == 100.0
        sum=PseudoCnts/wt_factor;                // add prior == 100.0
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(mat,mat,sum,sum,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
        for(x=startAB;x<=nAlpha(AB);x++) {
             if(MemSset(x,qst[j])){ mat += CountBG_j[x] + CountFG_j[x]; }
             else { sum += CountBG_j[x] + CountFG_j[x]; }
        }
        sum += mat;     // Failed to sum this in previous loop.
        assert(mat > SumXI[j]);
        mat -= SumXI[j];        // subtract out FG non-contamination, which doesn't count.
        sum -= SumXI[j];        // subtract out FG non-contamination, which doesn't count.
        theta[j][q_j] = mat/sum;
#endif

	theta_j = theta[j];
	q_j=query[j];
	m1=n1=m2=n2=0.0;
	Int4	CardPttrnSet=0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(n1,m1,n2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
	for(t=startAB; t<= nAlpha(AB); t++) {
	   if(MemSset(t,qst[j])){
		CardPttrnSet++;
		n1+= CountFG_j[t]; m1+=CountBG_j[t];
	   } else { n2+= CountFG_j[t]; m2+=CountBG_j[t]; }
	}
	submap=0.0;
	m1=m1*wt_factor; m2=m2*wt_factor;
	n1=n1*wt_factor; n2=n2*wt_factor;
	//*************** Rj terms: ***************
	submap += n1*log((1-alpha)*theta_j[q_j] + alpha); // delta_j == 1.
	submap += n2*log((1-alpha)*(1.0-theta_j[q_j]));	  // delta_j == 0.
	//*************** (1-Rj) terms: *************** 
	submap += m1*log(theta_j[q_j]);	
	submap += m2*log((1.0-theta_j[q_j]));
	assert(CardPttrnSet > 0);
	submap += logRho[j][CardPttrnSet];

	//*************** Null component ******************
	sum=PseudoCnts*wt_factor;
	double match=PseudoCnts*wt_factor;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(match,match,sum,sum,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
        for(x=startAB;x<=nAlpha(AB);x++) {
	   if(MemSset(x,qst[j])) match += CountBG_j[x] + CountFG_j[x];
	   else sum += CountBG_j[x] + CountFG_j[x];
        } sum += match;
	theta_x=match/sum;
	m1=0.0,m2=0.0;
#if 1	// treat deletions like random background...
	   if(del_as_random){
		// CntDelete(matFG,matBG,misFG,misBG,qst[j],CntBG[j][0],CountFG[j][0]);
		CntDelete(m1,m1,m2,m2,qst[j],CountBG_j[0],CountFG_j[0]);
	   }
#endif
	for(t=startAB; t<= nAlpha(AB); t++) {
	   if(MemSset(t,qst[j])){ m1+= CountBG_j[t] + CountFG_j[t]; }
	   else { m2+= CountBG_j[t] + CountFG_j[t]; }
	} m1 = m1*wt_factor; m2 = m2*wt_factor;
	submap -= m1*log(theta_x) + m2*log((1.0-theta_x));
	submap -= logRho[j][0];	// log of one minus other Rho values..
	return submap;
}


