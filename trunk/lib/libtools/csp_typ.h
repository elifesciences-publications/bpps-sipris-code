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

#if !defined(_CSP_TYP_)
#define _CSP_TYP_

#include "afnio.h"
#include "dheap.h"
#include "alphabet.h"
#include "sequence.h"
#include "residues.h"
#include "histogram.h"
#include "probability.h"
#include "random.h"
#include "sset.h"
#include "dsets.h"
#include "rst_typ.h"

#if 0	//************************************************************
	components: 
		1. Raw Pvalues (as -log10(p-values)). 
		2. LinearToLog values (for display).
		3. 
#endif	//************************************************************

class csp_typ {                 // class chain see p-value type.
public:
                csp_typ(double*,double*, double*,double*,a_type,double,rst_typ *);
	BooLean	ComputeBinomialBPPS(double, double *,double *);
	// BooLean	ReadBinomialBPPS(FILE *fp,double, double **);
	BooLean ComputeBinomialBPPS(double cutoff, double *pvalue,double *resEvals, 
								sst_typ psst, double lpr){
		// this routine not yet used; plan on writting LPR's from mcBPPS for chn_see.
			assert(cutoff <= 0.0); cutoff=-cutoff;
        		for(Int4 r=1; r <= nAlpha(AB); r++){
				if(MemSset(r,psst)) resEvals[r]=lpr; else resEvals[r]=0.0; 
			} *pvalue=lpr;
        		if(lpr < cutoff) return FALSE; else return TRUE;
		}
	BooLean	ComputeBinomial(double, double *,double *);
	void	SetModeBPPS(Int4 a0,Int4 b0, double a){
			assert(a > 0.0 && a <= 1.0);
			alpha = a; 
			assert(a0 > 0 && b0 > 0);
			A0=a0; B0=b0;
		}
                ~csp_typ( ){ Free(); }
private:
	// New BPPS version
	void    bnml_BPPS_dfs(Int4 rs, Int4 *ss, sst_typ sset, double n,
                	double m, double *min_Eval,double *min_p);
	void	Init(double*,double*,double*,double*, a_type,double,rst_typ *);
	double	CalcBiNomBPPS(double n, double m);
	double	CalcBiNomBPPS0(double n, double m); // can delete: old version
        void    Free( );                // free memory...

	// old ball-in-urn procedures here...
	BooLean ComputeBinomial0(double *, double , double *, double *);
	BooLean ComputeBinomial2(double *, double , double *, double *);
	void    cbp_dfs0(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
        	double *obs, double q, double *freq, double total,double *min_Eval,
        	double *min_p);
	BooLean cbp_dfs(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
        	double *obs, double q, double *freq, double total,double *min_Eval,
       		double *min_p);
	BooLean cbp_dfs2(Int4 rs, Int4 *ss, sst_typ sset, double Obs,
        	double *obs, double q, double *freq, double total,double *min_Eval,
       		double *min_p);

	BooLean	*IgnorePos;	// Started to try to improve -c<real> option...Finish later.
	a_type	AB;
	double	MinResFreqCutoff; // cutoff for reporting CBP p-values for conserved residue sets
	double	MinBlosum62Cutoff; // cutoff for defining related residues.
	double	AveBlosumScoreCutoff; // another cutoff for defining related residues.
	double	*ObsFG,*ObsBG;
	double	*FreqFG,*FreqBG;	// FreqFG CURRENTLY NOT USED...
	double	N,M;
	double	*ResEvals;
	ds_type	sets;
	double	alpha;		// contamination parameter for BPPS.
	Int4	A0,B0;
	rst_typ	*rst;
};

#endif

