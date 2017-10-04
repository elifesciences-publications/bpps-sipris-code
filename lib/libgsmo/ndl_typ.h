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

#if !defined (_NDL_TYP_)
#define	_NDL_TYP_
#include "cmsa.h"
#include "residues.h"
#include "histogram.h"
#include "gibbs.h"
#include "msaheap.h"
#include "sequence.h"
#include "smatrix.h"
#include "cluster.h"
#include "my_ncbi.h"
#include "blosum62.h"
#include "p2p_typ.h"
#include "lgm_typ.h"

/********************************************************************/
class ndl_typ {		// (position-specific) indel penalty type. 
public: 
		ndl_typ(Int4 pio,Int4 pdo,Int4 eie,Int4 ede,cma_typ c,lgm_typ *lgm,
							double pnts=1000,char wf=0)
			{  LGM=lgm; init(pio,pdo,eie,ede,pnts,c,wf); }
		~ndl_typ( ){ Free(); }

	//======================== ndl_typ.cc ============================
	BooLean	ChangePriors(Int4 pio,Int4 pdo,Int4 eie,Int4 ede)
		{ aa_per_ins_o=pio; aa_per_del_o=pdo;
        		exp_ins_ext=eie; exp_del_ext=ede; }

	//======================== ndl_hmm.cc ============================
	void	PutPenalties(FILE *fp, Int4 blk, Int4 I);
	void	PutPenalties(FILE *fp);
	double  IndelPenalty(FILE *fp,char mode,char *Wt=0);
	double  CalcIndelTransProb(FILE *fp,char mode,char *Wt);
	char    *GapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score,BooLean global=FALSE);
	void    AddIndels(Int4 sq,char Wt);
	void    RmIndels(Int4 sq,char Wt);
	void    SetPriorWt(double D){ assert(D >= 0.01 && D <= 100); prior_wt=D; }
	double	GetPriorWt( ){ return prior_wt; }
	void    InitIndels(char *Wt);	// recompute # transitions from scratch.
	double  GetParameters(Int4 &io,Int4 &d_o,Int4 &ie, Int4 &de,double &pw)
		{ io=aa_per_ins_o; d_o=aa_per_del_o; ie=exp_ins_ext; de=exp_del_ext; pw=prior_wt; }
	//======================== ndl_p2p.cc ============================
	char    *GapAlnTraceP2P(p2p_typ *P2P, Int4 &start,Int4 &score);
private:
	char    *gap_alignP2P(p2p_typ *P2P,Int4 *J,Int4 &start,Int4 &alnscore);
	//======================== ndl_p2p.cc ============================

	//======================== ndl_typ.cc ============================
	void	init(Int4, Int4, Int4, Int4, double,cma_typ,char);
	void	initTransitions();
	void	initTransProb();
	void	Free( ){ FreeTransitions( ); FreeTransProb( ); }
	void	FreeTransitions( );
	void	FreeTransProb( );
	//======================== ndl_typ.cc ============================
	cma_typ	CMA;
	UInt4   *Nmm,*Nmd,*Nmi,*Nm, *Nii,*Ndd,*Nim,*Ndm, *Nid,*Ndi,*Nsd,*Nsm;
	char	WtFact;
	BooLean	ReCalc;
	lgm_typ	*LGM;

	//======================== ndl_hmm.cc ============================
	char    *gap_align(Int4 n2,unsigned char *seq2,Int4 nmod,smx_typ *M,
				Int4 *J,Int4 *alnscore,Int4 *start,BooLean global=FALSE);
	void	RoundNearest(UInt4 &x,char wt){ x=(UInt4)floor(0.5+((double)x/(double)wt)); }
	UInt4	RoundNearest(double x,double wt){ return (UInt4)floor(0.5+(x/wt)); }
	//======================== ndl_hmm.cc ============================
	static const Int4 infinity = INT4_MAX/2;	// avoids over and underflow.
	// static const Int4 infinity = 999999999;
	// static const Int4 FloorI2I= 25;	// put a floor on i2i transition probability (25 == 0.025 nats).
	static const Int4 FloorI2I= 50;	// put a floor on i2i transition probability (25 == 0.025 nats).
	// static const Int4 FloorI2I= 70;	// put a floor on i2i transition probability (25 == 0.025 nats).
	static const Int4 FloorD2D= 50; // does mess up probability!!

	//======================== ndl_smpl.cc ============================
public:
	Int4	RtnFloorI2I(){ return FloorI2I; }
	Int4	RtnFloorD2D(){ return FloorD2D; }
	char    *SampleGapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score);
	void	GetTransProb(Int4 **mm,Int4 **mi,Int4 **md,Int4 **ii,Int4 **im,Int4 **dd,Int4 **dm,Int4 **bm)
		 { *mm=m2m; *mi=m2i; *md=m2d; *ii=i2i; *im=i2m; *dd=d2d; *dm=d2m; *bm=s2m; }
private:
	char    *sample_gap_align(Int4 sq_len,unsigned char *seq2,smx_typ smx, Int4 *J,
			Int4 *alnscore,Int4 *start);
	//======================== ndl_smpl.cc ============================
	
	Int4	*m2m,*m2i,*m2d,*m2s;
	Int4	*i2i,*i2m,*i2d;
	Int4	*d2d,*d2m,*d2i,*d2s;
	Int4	*s2m,*s2d;
	double	pernats;
        Int4    aa_per_ins_o;
        Int4    aa_per_del_o;
        Int4    exp_ins_ext;
        Int4    exp_del_ext;
        double	prior_wt;
#if 0
	static const double RatioPrior2Obs = 5.0;  // weight in favor of priors by this much.
	static const double TotalBeta = 5.0;  // Total number of counts as a function of temperature.
#endif

	//======================== ndl_debug.cc ============================
	void    PutDPMatrix(FILE *fp,Int4 **MAT,Int4 **TB,Int4 **DEL,Int4 **INS,
			Int4 sq_len,unsigned char *seq2,Int4 nmod,smx_typ *M,
			Int4 blk, Int4 startM, Int4 endM, Int4 startS, Int4 endS);
	void    put_dp_matrixSW(FILE *fp,Int4 **D, Int4** T, Int4 len, smx_typ M,
        		unsigned char *profile, unsigned char *seq, a_type A,
        		Int4 startM,Int4 RelstartM,Int4 start, Int4 end);
	void    put_dp_matrixSMX(FILE *fp, Int4 len, smx_typ M, unsigned char *profile,
        		unsigned char *seq, a_type A, Int4 startM,
			Int4 RelstartM,Int4 start,Int4 end);
	//======================== ndl_debug.cc ============================

	//======================== ndl_junk.cc ============================
	//======================== ndl_junk.cc ============================
	char	str[205];
};

#endif

