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

#if !defined(_RMA_TYP_)
#define _RMA_TYP_

#include "cmsa.h"
#include "HMM_typ.h"

class rma_typ {		// rich multiple alignment type
public:
                rma_typ( ){ assert(!"Illegal constructor"); }
                rma_typ(Int4,cma_typ);	// use default settings...
                rma_typ(Int4,char *,cma_typ);
//              rma_typ& operator=(const rma_typ&);     // assignment operator.
                ~rma_typ( ){ Free(); }
        void    Put(FILE *);
private:
        void    	init(Int4,char *pssm_arg,cma_typ);
        // void         copy(const rma_typ& idp);
        void    	Free();         // free memory...
	gsq_typ		csq;		// gapped consensus sequence.
	char		*operation;	// consensus sequence operation string.
	Int4		start;		// start for operation string.
	BooLean		own_cma;	// Is cma to be destroyed with rma_typ?
};

// char    FractionToCharRMA(double fract);
char    PvalueToCharRMA(double factor, double pval);
double *ComputRelEntropyRMA(Int4 length, char **conserved, char *info_show, Int4 t,
        double *observed, double ***Observed, double ***freq, sma_typ MA, a_type A,
        double cbp_cut,BooLean *use,char *pval_show,double *binomial_tail,
	double delta_p,double **ResEvals);
char    put_info_rma_typ(FILE *fptr, char info_show, char laststate, char *str);
char	get_conserved_state(char *conserved, char r,char c);
void	print_conserved_state(FILE *fptr, char state, char laststate,
                        char c, char *str);
Int4	get_font_color_code2(char ColorCode);
void	put_rasmol_rma(FILE *fp_ras, sma_typ MA, Int4 t,char su, a_type A, char *info_show,
        double *info_real);

double  ***SuperModelFreqsCMSA(Int4 rpts, char *pssm_arg, cma_typ cma,
	cma_typ main_cma,double ****Observed,Int4 JackCut);
float   **GetInfoCMSA(double purge_cut, double ***freq, sma_typ MA, 
	cma_typ cma,char mode, double ***Observed);
float   ***SubtractInfoCMSA(double purge_cut, Int4 rpts, char *pssm_arg, 
	cma_typ main_cma, cma_typ sub_cma, sma_typ subMA, char mode,Int4 JackCut);

void    NewCMA2RTF_HEADER(FILE *fptr,char PageSetUp);
void    NewCMA2RTF_HEADER(FILE *fptr,char PageSetUp,Int4 fontsize);
void    NewCMA2RTF_RETURNS(FILE *fptr,unsigned short numRtns);
void    NewCMA2RTF_TAIL(FILE *fptr);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, double ***freq, sma_typ MA,
        cma_typ cma,double ***Observed,char ColorCode);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, double ***freq, sma_typ MA,
        cma_typ cma,double ***Observed,FILE *fp_ras,char su,
	char ColorCode);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, double ***freq, sma_typ MA, cma_typ cma,
	double ***Observed,FILE *fp_ras,char su,double **MargProb,
	char *title,char ColorCode);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, double ***freq, sma_typ MA,
        cma_typ cma,FILE *fp_ras, char su,char ColorCode);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, cma_typ cma,double ***Observed,char ColorCode);
void    NewCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO,
        double infoHI, ss_type key_seq, cma_typ cma,char ColorCode);
void    SubtractInfoCMA2RTF(FILE *fptr, double cbp_cut, double purge_cut,
        double infoLO, double infoHI, ss_type key_seq, char PageSetUp,
        Int4 rpts, char *pssm_arg,cma_typ main_cma, cma_typ sub_cma,
        char mode,Int4 JackCut);
void    SubtractInfoCMA2RTF(FILE *fptr, char *filename, char subunit,
	double cbp_cut, double purge_cut, double infoLO, double infoHI, 
	ss_type key_seq, char PageSetUp, Int4 rpts, char *pssm_arg,
	cma_typ main_cma, cma_typ sub_cma, char mode,Int4 JackCut,
	double **MargProb);
double	WeightedCountsCMSA(cma_typ cma,Int4 *NoSeq,double *wNoSeq,double **wCnts);
double *HenikoffWeights2CMSA(cma_typ cma,double **wCnts);

#endif


