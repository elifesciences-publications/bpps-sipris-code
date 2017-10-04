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

#if !defined(_PPS_TYP_H_)
#define _PPS_TYP_H_
#include "alphabet.h"
#include "sset.h"
#include "probability.h"
#include "histogram.h"
#include "dheap.h"
#include "afnio.h"

class pps_typ {		// Partitioning with Pattern Selection.
public:
		// B = pseudocounts for SF; 
		pps_typ(sst_typ *Qst,unsigned char *Query, Int4 Query_len, double A0, 
			double B0, double *B, a_type A,UInt4 WtFactor,char Model,double Sigma,
			UInt4 *card_set,double Rho,char alt_measure,BooLean UseGlobalSqWts);
		~pps_typ() 	{ Free(); }
	double	LPR(UInt4 **CntBG,UInt4 **CntFG){
			double lpr = Map(CntBG,CntFG) - NullMap(CntBG,CntFG);
			// if(nullMap == (double) INT4_MIN) nullMap = NullMap(CntBG,CntFG);
			// fprintf(stderr,"lpr = %.2f\n",lpr);
			return lpr;
		}
	UInt4 NumColumns( ){ return NumCols; }
	double	ProbCols( );
	Int4	LenPattern( ){ return k; }
	BooLean	TheSamePattern(pps_typ *pps){
			Int4 k2 = pps->k;
			if(k != k2) return FALSE;
			for(Int4 j=1; j <= k; j++){
				if(qst[j] != pps->qst[j]) return FALSE;
			} return TRUE;	// same pattern.
		}
	void	AddColumn(Int4 j,sst_typ sstJ){
			if(j > 0 && j <= k){ 
				if(qst[j]==0) NumCols++;
				qst[j] = sstJ; ComputeNull=TRUE; 
			}
		}
	sst_typ	*RtnSST( ){ return qst; }
	sst_typ	*ModifySST(Int4 N,double min_nats, double min_freq){
		// Conceptually it changes the contrast of the alignment...
		// return top N most constrained residue positions.
		// pop off the N most strikingly conserved patterns.
		// WARNING: Assumes that subLPR[j] has been set!!!
		if(N==0) N=1000;
		dh_type dH=dheap(k+2,4);
		sst_typ *mqst;
		double	lpr;
		Int4	j,n;
		for(n=0,j=1;j<=k;j++){
		  if(qst[j] && subLPR[j] >= (double) min_nats){
		    double d=(double)MatchFG[j]/(double)(MatchFG[j]+MisMatchFG[j]);
		    if(d >= min_freq){ insrtHeap(j,(keytyp)-subLPR[j],dH); n++; }
		  }
		}
		if(n < N) N=n;
		NEW(mqst,k+5,sst_typ);
		for(n=0; !emptyHeap(dH) && n < N; n++){
		   lpr=-minkeyHeap(dH);
		   j=delminHeap(dH);
		   // fprintf(stderr,"j=%d; min_freq=%g; lpr=%g\n",j,min_freq,subLPR[j]);
		   mqst[j] = qst[j];
		}
		Nildheap(dH); 
		return mqst; 
	}
	sst_typ	*ModifiedSST(double min_nats,double min_freq){
		// WARNING: Assumes that subLPR[j] has been set!!!
			sst_typ *mqst;
			NEW(mqst,k+5,sst_typ);
			for(Int4 j=1;j<=k;j++){
			  if(qst[j] && subLPR[j] >= (double) min_nats){
			    double d=(double)MatchFG[j]/(double)(MatchFG[j]+MisMatchFG[j]);
			    // fprintf(stderr,"d=%g; min_freq=%g\n",d,min_freq);
			    if(d >= min_freq) mqst[j] = qst[j];
			  }
			}
			return mqst; 
		}
	sst_typ	RtnSST(Int4 j){
			if(j > 0 && j <= k){ return qst[j]; }
			else return 0;
		}
	void	RmColumn(Int4 j){
			if(j > 0 && j <= k){ 
				if(qst[j]) NumCols--;
				qst[j] = 0; ComputeNull=TRUE; 
			}
		}
#if 0
	double ColLPR(UInt4 *CntBG,UInt4 *CntFG){
			
		}
#endif
	double	*SubLPR(UInt4 **CntBG,UInt4 **CntFG);
	void	PutSubLPR(FILE *fp,UInt4 **CntBG,UInt4 **CntFG)
			{ PutSubLPR(fp, 0,CntBG,CntFG); }
	void	PutSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG);
	double	Alpha( ){ return alpha; }
	void	SetModel(char Model){
		  assert(Model == 'M' || Model == 'B' || Model == 'U' 
				|| (Model <= '9' && Model >= '0'));
	       	  if(Model <= '9' && Model >= '0'){
                	ModelSwitch = Model; ModelMode='U';
		  } else { ModelMode=Model; ModelSwitch = '?'; }
		}
private:
	void	Init(sst_typ *Qst,unsigned char *Query, Int4 Query_len, double A0, 
			double B0, double *B, a_type A,UInt4 WtFactor,
			char Model,double Sigma,UInt4 *card_set,double Rho, 
			char alt_measure);

	double  ProbRatio(unsigned char *, unsigned char **,
			UInt4 **CntBG, UInt4 **CntFG);
	double	BinomialDistribution(double N, double p, double n){
			return exp(LnBinomialDistribution(N,p,n)); 
		}
	double	LnBinomialDistribution(double N, double p, double n){
		 return (lngamma(N+1)-lngamma(n+1) - lngamma(N-n+1)
                                                + n*log(p) + (N-n)*log(1.0 - p)); 
		}
#if 0
	double  CBKLD(double N, double n, double q, double p) {
		   double P=CumBinomProb(n,N,q)*(LnCBP(n,N,q)-LnCBP(n,N,p)); 
        	   P+=CumBinomProb(N-n,N,1-q)*(LnCBP(N-n,N,1-q)-LnCBP(N-n,N,1-p));
        	   return P;
		}
#else	
	double  CBKLD(double N, double n, double q, double p);
#endif
	double  BKLDfactor(double N, double n, double p) {
			return (n*log(p) + (N-n)*log(1.0 - p));
		}
	double  BKLD(double N, double n, double q, double p) {
        	 double P=0.0,d;
		 d=LnBinomialProb(N,q,n);
		 P=exp(d)*(BKLDfactor(N,n,q)-BKLDfactor(N,n,p));
        	 d=LnBinomialProb(N,1-q,N-n);
        	 P+=exp(d)*(BKLDfactor(N,N-n,1-q)-BKLDfactor(N,N-n,1-p));
		 if(isfinite(P)) return P;
		 else {                  // Fix core dump...
                   fprintf(stderr,"NaN: N=%g; n=%g; q=%g; p=%g\n",N,n,q,p);
                   return 0.0;
        	 }
		}
	// Abstract model:
	char	ModelMode,ModelSwitch;
#if 0
	double	rho;	// prior for column being a pattern position (independent Bernoullis).
		// rho = 0.5 has no effect = non-informed prior...
		// put infinite weight on this prior (don't update).
#endif
	double	logRho; // penalties for adding a column.
	double	logOneMinusRho; // penalty for removing a column.

	double 	Map(UInt4 **CntBG,UInt4 **CntFG){
		  if(ModelMode == 'M') return MapMultiNom(CntBG,CntFG);
		  else if(ModelMode == 'U') return MapBallInUrn(CntBG,CntFG);
		  else return MapBiNom(CntBG,CntFG);
		}
	double	NullMap(UInt4 **CntBG, UInt4 **CntFG){
		  if(ModelMode == 'M') return NullMapMultiNom(CntBG,CntFG);
		  else if(ModelMode == 'U') return 0.0;
		  else return NullMapBiNom(CntBG,CntFG);
		}
	double 	subMap(UInt4 **CntBG,UInt4 **CntFG){
		  if(ModelMode == 'M') return subMapMultiNom(CntBG,CntFG);
		  else if(ModelMode == 'U') return subMapBallInUrn(CntBG,CntFG);
		  else return subMapBiNom(CntBG,CntFG);
		}
	double	subNullMap(UInt4 **CntBG, UInt4 **CntFG){
		  if(ModelMode == 'M') return subNullMapMultiNom(CntBG,CntFG);
		  else if(ModelMode == 'U') return 0.0;
		  else return subNullMapBiNom(CntBG,CntFG);
		}
	// Multinomial model:
	double  MapMultiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	NullMapMultiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapMultiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subNullMapMultiNom(UInt4 **CntBG, UInt4 **CntFG);
	// BallInUrn model:
	double	MapBallInUrn(UInt4 **CntBG, UInt4 **CntFG);
	double	MapBallInUrn2(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapBallInUrn(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapBallInUrn2(UInt4 **CntBG, UInt4 **CntFG);
	double 	submap_ball_in_urn(double m1, double m2, double n1, double n2);
	// Binomial model:
	double	MapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	NullMapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subNullMapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	void	BinaryUpdateSumXI(UInt4 **CountFG,UInt4 **CntBG);
	void	BinaryUpdateTheta(UInt4 **CntBG,UInt4 **CntFG);
	// End Binomial model:
	void 	Free();
	void	updateSumXI(UInt4 **,UInt4 **);
	// void	updateSumXI(UInt4 **);
	void	updateAlpha(UInt4 **CountFG);
	void	updateTheta(UInt4 **CntBG,UInt4 **CntFG);
	// void	updateTheta(UInt4 **CntBG);
	UInt4	NumCols;
	double	BinomialWtdRE(double N, double n, double p, double q);
	double  GetGapsMatch(double CntX, sst_typ qst);
#if 1	// NEW
	// unsigned char	**qf_set,**sf_set;	// residue sets.
	BooLean	ComputeNull;
#endif
	unsigned char *query;
	double	*SubNullMap,*SubMap,*MatchFG,*MisMatchFG,*MatchBG,*MisMatchBG;
	double	*subLPR;
	double	nullMap;
	sst_typ	*qst;	// query residue set (pattern) at each position.
	a_type	AB;	// alphabet of residue types.
	Int4 	k;	// k = length of the query (symbol used in Jun's paper).
	double 	alpha;	// == parameter specifying relative fraction of 
			 // 	conserved and unconserved residues.
	double	a0,b0;	 // prior pseudocounts for alpha (Beta distribution).
	double	**theta; // residue compositions for the non-query cluster.
	double	*NullTheta; // residue compositions for the binary Null model.
	double	*b;	 // priors for theta (product Dirichlets).
	double	*SumXI;	 // Sum of (greek letter) Xi, an indicator variable for query partition.
	double	wt_factor;	// integerized factor for sequence weights (using Henikoff method?).
				// e.g., wt_factor = 1/100 means that a count of 100 corresponds
				// to a single sequence.
	double	sigma;
	double	*argum;
	// double	stdfreq[30];
	char	AltMeasure;	// 
	UInt4	*CardSet;	// cardinality of sets.
	unsigned char	StartAlpha;
	BooLean	DebugCheck;
};

#endif
