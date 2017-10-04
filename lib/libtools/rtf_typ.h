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

#if !defined(_RTH_TYP_)
#define _RTH_TYP_

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
#include "csp_typ.h"

#if 0	//************************************************************
	components: 
		1. Raw Pvalues (as -log10(p-values)). 
		2. LinearToLog values (for display).
		3. 
#endif	//************************************************************

class rtf_typ {                 // class rich text histogram.
public:
                rtf_typ( ){ assert(!"Illegal constructor"); }
                rtf_typ(double,Int4,Int4,double, double**,double**,
			Int4,e_type,a_type,Int4,double,char);
                rtf_typ(double,Int4,Int4,double, double**,double**, 
			Int4,e_type,a_type,Int4,double,char,
			double**,double**,char,double,Int4,Int4);
#if 1		// for cmc_typ
		rtf_typ(double,Int4,Int4,double, double**,double**,
			Int4,e_type,a_type,Int4,double,char,
			double**,double**,double**);
#endif
                ~rtf_typ( ){ Free(); }
        void    PutLine(FILE *,Int4,Int4,Int4,Int4,char*,Int4,Int4,Int4);
        void    PutLine(FILE *,Int4,Int4,Int4,Int4,char*,Int4,Int4,Int4,Int4,Int4*);
        BooLean	Printable( ){ return (MaxPval > MinPval); }
        Int4	Length(){ return length; }
	// For omcBPPS cross conserved histogram: need rtf->SwapValue().
        double	*SwapValue(double *value,double *&pv){
		   static Int4 calls=0;
		   double *tmp=Value,*Tmp=PValue; PValue=pv; pv=Tmp; Value=value; 
		   if(0 && calls==0){ for(Int4 i=1; i <= length; i++)
		      fprintf(stderr,"Value[%d]=%.3f vs %.3f; PValue[%d]=%.3f vs %.3f\n",i,tmp[i],Value[i],i,Tmp[i],PValue[i]);
		      fprintf(stderr,"MinPval = %.3f\n",MinPval);
		   } calls++;
		   return tmp; 
		}
        double	*ResEvalues(Int4 i){
			if(i < 1 || i > length) return 0; else return res_evals[i]; }
	char    *Conserved(Int4 j);
	char    *Conserved(double *ResEvals);
	char    *Conserved(BooLean Ignore, double *ResEvals);
	char    ConservedState(char *,char,char);
	void	PutConservedState(FILE *,char ,char ,char ,char *);
	char    PvalueToChar(double factor, double pval);
	char    ShowInfo(FILE *, char , char , char *);
	void    PutLineBinomial(FILE *fptr,Int4 start,Int4 end, Int4 gstart,Int4 gend,
                	char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2);
	void    PutBarHeights(FILE *fptr,Int4 start,Int4 end, Int4 gstart,
                	Int4 gend, char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2);
	void    PutLineRelEntropy(FILE *fptr,Int4 start,Int4 end, Int4 gstart,
                Int4 gend, char *gnull,Int4 color,Int4 colorB,rtf_typ *rtf2);
	Int4	PutResEvals(FILE *fp,Int4 start,Int4 end, Int4 gstart,Int4 gend,char *gnull,
			double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,BooLean verbose,
			double WtNumSq,Int4 RawNumSeqs,BooLean PatternOnlyMode,
			Int4 MaxConcensusLines);
	Int4	PutResEvals(FILE *fp,Int4 start,Int4 end, Int4 gstart,Int4 gend,char *gnull,
			double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,BooLean verbose,
			double WtNumSq,Int4 RawNumSeqs,BooLean PatternOnlyMode,
			Int4 MaxConcensusLines, Int4 maxlen_gnull,Int4 *gnull_insrt_len);
	Int4	PutResEvals(FILE *fp,Int4 start,Int4 end, Int4 gstart,Int4 gend,char *gnull,
			double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,BooLean verbose,
			double WtNumSq,Int4 RawNumSeqs,BooLean PatternOnlyMode, Int4 MaxConcensusLines,
			Int4 maxlen_gnull,Int4 *gnull_insrt_len,BooLean UnderLine,BooLean IsFG);
#if 1	// for X-conserved; need a lot more work...
	Int4    PutResEvalsXC(FILE *fptr,Int4 start,Int4 end,Int4 gstart, Int4 gend,char *gnull,
			double **wtfreq,double *wtnsq,Int4 color,rtf_typ *rtf,
           		Int4 MaxConcensusLines,Int4 maxlen_gnull,Int4 *gnull_insrt_len, BooLean UnderLine);
#endif
	Int4	FontSize(){ return fontsize; }
	Int4	BarHeight(Int4 i){ 
			if(i < 1 || i > length) return 0;
			// return GetHistHeight(IgnorePos[i],Value[i]);
			return GetHistHeight(Value[i]);
		}
	Int4	ResBarHeight(Int4 i, Int4 r){ 
			if(i < 1 || i > length) return 0;
			if(r < 1 || r > nAlpha(AB)) return 0;
			// return GetHistHeight(IgnorePos[i],ResValue[i][r]); 
			return GetHistHeight(ResValue[i][r]); 
		}
	a_type	Alphabet(){ return AB; }
	double	GetValue(Int4 x){ if(x > 0 && x <= length) return Value[x]; else return 0.0; }
	double	GetPValue(Int4 x){ if(x > 0 && x <= length) return PValue[x]; else return 0.0; }
	double	GetLinearToLog( ){ return LinearToLog; }
	char    *DeriveStatus(char *statusSF,rtf_typ *rtfQ,char category);
	Int4	FixSetStatus(char *status,rtf_typ *rtfQ);
	void	ReInit(double LineToLog,double MinPvalue){
			for(Int4 s=1; s<=length; s++) free(ResValue[s]);
			free(ResValue); free(PValue); free(Value);
			// reinitialize histogram and highlighting parameters.
			init(LineToLog,binomial,MinPvalue);
		}
	// Frequency of conservation cutoffs 
	BooLean	SetMinResFreqCutoff(double value){
		// Set cutoff for reporting CBP p-values for conserved residue sets
			if(value >= 0.0 && value <= 1.00){
				MinResFreqCutoff=value;  return TRUE;
			} return FALSE;
		}
	BooLean	SetMinBlosum62Cutoff(double value){
		// Set cutoff for defining related residues.
			if(value >= -10000.0 && value <= 10.00){
			   MinBlosum62Cutoff=value; return TRUE;
			} return FALSE;
		}
	BooLean	SetAveBlosumScoreCutoff(double value){
		// Set second cutoff for defining related residues.
			if(value >= -10000.0 && value <= 1000.00){
				AveBlosumScoreCutoff=value; return TRUE;
			} return FALSE;
		}
	BooLean	IgnoreResPos(Int4 j){		// new ignore single position...
			if(IgnorePosition && Position==j) return TRUE; 
			else return FALSE;
		}
	void	SetExpPttrns(double x){ ExpPttrns=x; }
private:
	double	ExpPttrns; // as log10...
	void	Construct(double,Int4,Int4,double,double**,double**,
			double**,double**, Int4,e_type,a_type,Int4,
			double,char,char,double,Int4,Int4,double **);
	void    init(double LineToLog,double *PvalReal,double MinPvalue);
        double  MaxPvalue( ){ return MaxPval; }
	double	GetLineToLog(Int4 , double *, double ,double );
	BooLean	ComputeBinomial(double *, double , double *, double *,double *);
	double  ComputeRelEntropy(double *observed, double *freq);
	double  GetValue(double pvalue);
	char    IntegerToChar(Int4 x);
	char    FractionToChar(double fract);
	Int4    get_unit_height(Int4 fontsize);
	Int4	GetHistHeight(BooLean Ignore,double V){
		    Int4 x;
		    if(V < 0.0) return 0;
		    else if(Ignore) return 1;
		    x = (Int4) floor(V+0.5);
		    if(x > hist_height) return hist_height;
		    else if(x < 1) return 1;
		    else return x;
		}
	Int4	GetHistHeight(double V){
		    Int4 x;
		    if(V < 0.0) return 0;
		    x = (Int4) floor(V+0.5);
		    if(x > hist_height) return hist_height;
		    else if(x < 1) return 1;
		    else return x;
		}
        void    Free( );                // free memory...

	BooLean	*IgnorePos;	// Started to try to improve -c<real> option...Finish later.
        double  LinearToLog;
	Int4	bars_s_cut,minbars;
        Int4    fontsize,length;
        Int4    hist_height;            // max number of red bars == 30
        Int4    hist_up;
        double  *PValue,MaxPval,MinPval,TargetPval;
        double  *Value,MaxValue,MaxRatio;
	double	**res_evals,*binomial,**ResValue;
	double	*RelEntropy;
	a_type	AB;
	e_type	Query;	// If Query != NULL then use Query-centric mode...
        char    mode;
	Int4	Position;
	BooLean	IgnorePosition;	// new ignore single position...
	char	Residue;
	// experimental
	BooLean	cbp_dfs(Int4 rs, Int4 *ss, sst_typ sset, ds_type sets, double Obs,
		double *obs, double q, double *freq, double total,
		double *min_Eval,double *min_p,double *ResEvals);
	void	cbp_dfs0(Int4 rs, Int4 *ss, sst_typ sset, ds_type sets, double Obs,
		double *obs, double q, double *freq, double total,
		double *min_Eval,double *min_p,double *ResEvals);
	double	MinResFreqCutoff; // cutoff for reporting CBP p-values for conserved residue sets
	double	MinBlosum62Cutoff; // cutoff for defining related residues.
	double	AveBlosumScoreCutoff; // another cutoff for defining related residues.
	rst_typ	*rst;
};

Int4    get_font_color_code(char ColorCode);

#endif

