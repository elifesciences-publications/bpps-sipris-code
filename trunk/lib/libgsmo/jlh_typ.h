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

#if !defined (_JLH_TYP_)
#define	_JLH_TYP_
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

/********************************************************************/
class jlh_typ {		// Jun Liu HMM type
public: 
		jlh_typ(Int4 pio,Int4 pdo,Int4 eie,Int4 ede,cma_typ c,double pnts=1000,char wf=0)
			{ init(pio,pdo,eie,ede,pnts,c,wf); }
		~jlh_typ( ){ Free(); }

	//======================== jlh_typ.cc ============================
	BooLean	ChangePriors(Int4 pio,Int4 pdo,Int4 eie,Int4 ede)
		{ aa_per_ins_o=pio; aa_per_del_o=pdo;
        		exp_ins_ext=eie; exp_del_ext=ede; }

	//======================== jlh_hmm.cc ============================
	double  IndelPenalty(FILE *fp,char mode,char *Wt=0);
	char    *GapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score);
	void    AddIndels(Int4 sq,char Wt);
	void    RmIndels(Int4 sq,char Wt);
	void    InitIndels(char *Wt);	// recompute # transitions from scratch.

	//======================== jlh_p2p.cc ============================
	char    *GapAlnTraceP2P(p2p_typ *P2P, Int4 &start,Int4 &score);
private:
	char    *gap_alignP2P(p2p_typ *P2P,Int4 *J,Int4 &start,Int4 &alnscore);
	//======================== jlh_p2p.cc ============================

	//======================== jlh_typ.cc ============================
	void	init(Int4, Int4, Int4, Int4, double,cma_typ,char);
	void	Free( ){ }
	//======================== jlh_typ.cc ============================
	cma_typ	CMA;
	UInt8   Nmm,Nmd,Nmi,Nm, Nii,Ndd,Nim,Ndm, Nid,Ndi,Nsd,Nsm;
	char	WtFact;
	BooLean	ReCalc;

	//======================== jlh_hmm.cc ============================
	char    *gap_align(Int4 n2,unsigned char *seq2,Int4 nmod,smx_typ *M,
				Int4 *J,Int4 *alnscore,Int4 *start);
	void	RoundNearest(UInt8 &x,char wt){ x=(UInt8)floor(0.5+((double)x/(double)wt)); }
	UInt8	RoundNearest(UInt8 x,double wt){ return (UInt8)floor(0.5+((double)x/wt)); }
	//======================== jlh_hmm.cc ============================
	// static const Int4 infinity = INT4_MAX/2;	// avoids over and underflow.
	
	static const Int4 infinity = 999999999;
	Int4	m2m,m2i,m2d,m2s;
	Int4	i2i,i2m,i2d;
	Int4	d2d,d2m,d2i,d2s;
	Int4	s2m,s2d;
	double	pernats;
        Int4    aa_per_ins_o;
        Int4    aa_per_del_o;
        Int4    exp_ins_ext;
        Int4    exp_del_ext;
        Int4    prior_wt;
#if 0
	static const double RatioPrior2Obs = 5.0;  // weight in favor of priors by this much.
	static const double TotalBeta = 5.0;  // Total number of counts as a function of temperature.
#endif

	//======================== jlh_put.cc ============================
	void    put_dp_matrixSW(Int4 **D, Int4** T, Int4 len, smx_typ M,
        		unsigned char *profile, unsigned char *seq, a_type A,
        		Int4 startM,Int4 RelstartM,Int4 start, Int4 end);
	void    put_dp_matrixSMX(Int4 len, smx_typ M, unsigned char *profile,
        		unsigned char *seq, a_type A, Int4 startM,
			Int4 RelstartM,Int4 start,Int4 end);
	//======================== jlh_put.cc ============================

	//======================== jlh_debug.cc ============================

	//======================== jlh_junk.cc ============================
	//======================== jlh_junk.cc ============================
	char	str[205];
};

#endif

