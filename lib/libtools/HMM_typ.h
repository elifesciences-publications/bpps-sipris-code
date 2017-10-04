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

#if !defined(_OHMM_TYP_)
#define _OHMM_TYP_

#include "cmsa.h"
#include "tpb_typ.h"
#include "new-dsc.h"
#include "psg.h"
#include "histogram.h"

class HMM_typ {		// position-specific-matrix type
public:
                HMM_typ( ){ assert(!"Illegal constructor"); }
                HMM_typ(Int4,cma_typ,double *);	// use default settings...
		HMM_typ(Int4,char *,cma_typ,Int4,double *);
//              HMM_typ& operator=(const HMM_typ&);     // assignment operator.
                ~HMM_typ( ){ Free(); }
        void	ReComputeSMX(cma_typ);
        void    Put(FILE *);
        char	*ArgTP( ){ return tpb->Arg( ); }
        smx_typ	*GetSMX(){ return smx; }
	void	PutAlign(FILE*, e_type,char *,Int4, Int4 , Int4 );
	void	PutAlign(FILE*, e_type,char *,Int4, Int4 , Int4 ,BooLean);
	void	SetRptPval(double p){ assert(p>0.0 && p<=1.0); target_pvalue=p; }
	Int4    GetScore(e_type,char *,Int4,Int4);
	char    *Align(e_type,Int4, Int4 *, Int4 *, Int4 *);
	double	*MarginalProb(e_type,char *,Int4,Int4);
	double	*MarginalProb0(e_type,char *,Int4,Int4);
	char    *SampleAlign(e_type,Int4,Int4 *,Int4 *,Int4 *,Int4);
	char    *Align(FILE*, e_type,Int4, Int4 *, Int4 *, Int4 *);
	double	*ThresholdScore(e_type, double, Int4 *, Int4 *);
	char	*FindBestRpts(e_type, Int4, Int4 *, Int4 *, Int4 *,Int4 *);
private:
	Int4		*InsEmitProb(unsigned char *, Int4);
	char		*Viterbi(e_type,Int4,Int4 *,Int4 *,Int4 *);
	void		GetProb(Int4,Int4,cma_typ);
        void    	init(Int4,char*,Int4,cma_typ,double *);
        // void         copy(const HMM_typ& tpb);
        void    	Free();         // free memory...
        unsigned short  MaxRpts, length, nblk,*lenblk;
	Int4		*start_block,*end_block;
	double		*prob_c,*prob_s,*prob_h;
	Int4		PerNats;
	tpb_typ		*tpb;	// insertion & deletion penalties.
	smx_typ		*smx;
	Int4		*ins_emit_prob,**mat_emit_prob;
	a_type		AB;
	char		*consensus;
	double		target_pvalue;  // for a single repeat...
	psg_typ		*psg;
	Int4		**MAT,**DEL,**INS;
	char		**traceM,**traceI,**traceD;
	void		inner_loop(Int4 j, Int4 *matj, Int4 *matjm1,
           		 Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           		 Int4 m2m, Int4 d2m, Int4 i2m, Int4  m2d, Int4  d2d, Int4  m2i, 
			 Int4 i2i, Int4 seq_len, unsigned char *seq, Int4 *scorej,Int4 *ins_emit);
	char    	*get_traceback(Int4, Int4, Int4, Int4 *, Int4 *);
};

#endif
