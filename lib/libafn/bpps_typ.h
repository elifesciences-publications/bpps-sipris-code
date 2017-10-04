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

#if !defined (_BPPS_TYP_H_)
#define _BPPS_TYP_H_
#include "alphabet.h"
#include "sset.h"
#include "probability.h"
#include "histogram.h"
#include "dheap.h"
#include "afnio.h"
#include "sequence.h"
#include "rst_typ.h"

class bpps_typ {		// Binary Partitioning with Pattern Selection.
public:
		// B = pseudocounts for SF; 
		bpps_typ(sst_typ *Qst,unsigned char *Query, Int4 Query_len, double A0, 
			double B0, a_type A, double **Rho,double priorRi,BooLean UseGlobalSqWts);
private:
	void	Init(sst_typ *Qst,unsigned char *Query, Int4 Query_len, double A0, 
			double B0, double **Rho, double priorRi);
	static const UInt4  WtFactor=100;
	double	*beta;
	// static const Int4 startAB=0; //  treat deletions ('-'='x') as non-functional residues.
	// static const unsigned char startAB=0; // treat deletions ('-'='x') as non-functional residues.
	static const unsigned char startAB=1;  // treat seqs with deletions ('-'='x') as non-existent
	Int4    CntDelete(double &matFG, double &matBG, double &misFG, double &misBG,
					sst_typ xsst, UInt4 DelBG,UInt4 DelFG);
	unsigned char del_as_random;	// 0 == ignore; 1 == treat as random background.
	double	dummy,Dummy;
	char	Type;
public:
	bpps_typ(sst_typ *qsst, unsigned char *Qry, Int4 Query_len, char type, a_type A,
			BooLean UseGlobalSqWts, double probRi=0.5);
	rst_typ	*RtnRST(){ return RST; }
	BooLean	UpdateRho(sst_typ **&sst,e_type keyE);
	BooLean	UpdateSST(sst_typ **&sst,e_type keyE,FILE *efp=0);
private:
	void    Initialize(sst_typ *qsst, double A0, double B0,double priorRi, double rho,
                        Int4 Query_len, unsigned char *Qry);
	BooLean	own_beta_query;
	rst_typ	*RST;
	void	UpdateRho(double **rho=0);
public:
		~bpps_typ() 	{ Free(); }
// NEW for mcBPPS procedure:
	UInt4	RtnWtFactor(){ return WtFactor; }
	void	TreatDeletionAsRandom(){ del_as_random=1; }
	void	UpdateRhoJ(Int4 j,double *rho);
	void    PutRho(FILE *fp);
	void    PutRho(FILE *fp,Int4 i);
	Int4    Compare(FILE *fp,bpps_typ *that);
	double	**RtnRho(){ return Rho; }
	double	RawLPR(UInt4 **CntBG,UInt4 **CntFG){
		// returns the actual log probability; use this for optimization...
			double lpr = Map(CntBG,CntFG);
			return lpr;
		}
	double	LPR(UInt4 **CntBG,UInt4 **CntFG){
			double lpr = Map(CntBG,CntFG) - NullMap(CntBG,CntFG);
			// if(nullMap == (double) INT4_MIN) nullMap = NullMap(CntBG,CntFG);
			// fprintf(stderr,"lpr = %.2f\n",lpr);
			return lpr;
		}
	UInt4	NumColumns( ){ return NumCols; }
	Int4    *nBestPttrn(Int4 &N,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG);
	double	ProbCols( );
	Int4	LenPattern( ){ return k; }
	BooLean	TheSamePattern(bpps_typ *pps){
			Int4 k2 = pps->k;
			if(k != k2) return FALSE;
			for(Int4 j=1; j <= k; j++){
				if(qst[j] != pps->qst[j]) return FALSE;
			} return TRUE;	// same pattern.
		}
	void	AddColumn(Int4 j,sst_typ sstJ){
			if(j > 0 && j <= k){ 
				if(qst[j]==0 && sstJ) NumCols++;
				qst[j] = sstJ; ComputeNull=TRUE; 
			}
		}
	sst_typ	*RtnSST( ){ return qst; }
	void	ReplacePattern(sst_typ *new_qst)
		{
			NumCols=0;
			for(Int4 i=1; i <=k; i++){
			    qst[i]=new_qst[i];
			    if(qst[i]) NumCols++;
			} ComputeNull=TRUE;
		}
	sst_typ	*RtnCopySST( )
		{ 
			sst_typ *sst; NEW(sst,k+3,sst_typ);
			for(Int4 i=1; i <=k; i++) sst[i]=qst[i];
			return sst; 
		}
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
	unsigned char *RtnQuery(){ return query; }
	e_type	RtnQuerySeq(){ return MkSeq("bpps query seq", k, query); }
	double	*RtnBetaPriors(){ return b; }
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
	double	*SubLPR(UInt4 **CntBG,UInt4 **CntFG);
	double	*SubRawLPR(UInt4 **CntBG,UInt4 **CntFG);
	void    PutCDTreeSubLPR(FILE *fp,char *id,UInt4 **CntBG,UInt4 **CntFG);
	void	PutSubLPR(FILE *fp,UInt4 **CntBG,UInt4 **CntFG)
			{ PutSubLPR(fp, 0,CntBG,CntFG); }
	void	PutSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG);
	void    PutPattern(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG);
	double	PutBoltzmannLike(FILE *fp,double *StdFrq, Int4 j, UInt4 **CntBG,UInt4 **CntFG,Int4 n); // in omc_debug.cc
	void    WriteSubLPR(FILE *fp,Int4 *fake2real,UInt4 **CntBG,UInt4 **CntFG);
	double	Alpha( ){ return alpha; }
	double	Parameters(double &A0, double &B0){
			A0=a0*wt_factor; B0=b0*wt_factor; return alpha; 
		}
	//===== New public routines for tri-partitioned model. =====
	void    DefineNullModel(Int4 leng, sst_typ *null_sst, double **null_freq);
	sst_typ FindBestPattern(sst_typ *sstLegal, double *WtCnts);
	//===== End routines for tri-partitioned model. =====
	//===== Routine for pmcBPPS program. =====
	sst_typ *RtnOptPattern(UInt4 **, UInt4 **, sst_typ **, unsigned char *,Int4 ,double &);
	double	SeqLikelihoodFG(unsigned char * seq, UInt4 SqWt);
	double	SeqLikelihoodBG(unsigned char * seq, UInt4 SqWt);
	void	UpdateParameters(UInt4 **CntBG, UInt4 **CntFG);
	double	QuickLPRsubJ(Int4 j, UInt4 **CntBG, UInt4 **CntFG);
	void	PutParameters(FILE *fp){
			fprintf(fp,"A0=%.1f B0=%.1f Ri=%.3f \n",a0*wt_factor,b0*wt_factor,PriorRi);
		}
private:
	//===== Routines for pmcBPPS program. =====
	void	BinaryUpdateSumXI_J(Int4 j, UInt4 **CountFG,UInt4 **CntBG);
	void	BinaryUpdateTheta_J(Int4 j, UInt4 **CntBG,UInt4 **CntFG);
	double	subNullMapBiNom_J(Int4 j, UInt4 **CntBG, UInt4 **CntFG);
	void	UpdateNullBiNomTheta(UInt4 **CntBG, UInt4 **CntFG);
	double	subAlphaProb(double *n_f, double *n_n, double theta_f);
	double	subAlphaNullProb(double *n_f, double *n_n, double theta_f);
	//===== End of routines for pmcBPPS program. =====

	//===== New private routines for tri-partitioned model. =====
	//===== New private routines for testing model. =====
	double	NullMapBiNom2(UInt4 **CntBG, UInt4 **CntFG);
	double	subNullMapBiNom2(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapBiNom2(UInt4 **CntBG, UInt4 **CntFG);
	double	*NullFreq;
	sst_typ	*NullSST;	// optimum pattern at each position (passed in 
	//===== End routines for tri-partitioned model. =====

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
	double	rho;	// prior for column being a pattern position (independent Bernoullis).
		// rho = 0.5 has no effect = non-informed prior...
		// put infinite weight on this prior (don't update).
#elif 0
	double	logRho[32]; // penalties for adding a column.
	double	logOneMinusRho[32]; // penalty for removing a column.
#else	// categorical distribution passed in from calling environment...
	double	**logRho; // penalties for adding a column.
	double	**Rho; // penalties for adding a column.
#endif
	double	PriorRho;	// prior for column being a pattern position (independent Bernoullis).
	

	// implement priors for sampling sequences into foreground...
	double	PriorRi;	// prior probability of Ri == 1.
	double	LogRatioPriorRi;	// only computed in regular map!!!! not Null MAP.

	double 	Map(UInt4 **CntBG,UInt4 **CntFG){ return MapBiNom(CntBG,CntFG); }
	double	NullMap(UInt4 **CntBG, UInt4 **CntFG){ return NullMapBiNom(CntBG,CntFG); }
	double 	subMap(UInt4 **CntBG,UInt4 **CntFG){ return subMapBiNom(CntBG,CntFG); }
	double	subNullMap(UInt4 **CntBG, UInt4 **CntFG){ return subNullMapBiNom(CntBG,CntFG); }
	// BallInUrn model:
	double	MapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	NullMapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subMapBiNom(UInt4 **CntBG, UInt4 **CntFG);
	double	subNullMapBiNom(UInt4 **CntBG, UInt4 **CntFG);

	void	BinaryUpdatePseudoCnts(UInt4 **CountFG,UInt4 **CntBG);
	double	PseudoCnts;
	void	BinaryUpdateSumXI(UInt4 **CountFG,UInt4 **CntBG);
	void	updateAlpha(UInt4 **CountFG);
	void	FastUpdateAlpha(UInt4 **CountFG);
	void	BinaryUpdateTheta(UInt4 **CntBG,UInt4 **CntFG);
	// End Binomial model:
	void 	Free();
	UInt4	NumCols;
#if 1	// NEW
	// unsigned char	**qf_set,**sf_set;	// residue sets.
	BooLean	ComputeNull;
#endif
	unsigned char *query;
	double	*SubNullMap,*SubMap,*MatchFG,*MisMatchFG,*MatchBG,*MisMatchBG;
	double	*subLPR;
	double	*subRawLPR;
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
	double	*argum;
	// double	stdfreq[30];
	char	AltMeasure;	// 
	unsigned char	StartAlpha;
	BooLean	DebugCheck;
#if 0	// for categorical distribution of legal residue sets
	sst_typ	**LegalSST;	// legal sets.
	Int4	*NumSets;	//
	double	**prob_sets;	// should sum to 1.0 when the null set is included.
#endif
};

#endif
