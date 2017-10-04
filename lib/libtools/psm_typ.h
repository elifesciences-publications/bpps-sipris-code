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

#if !defined(_PSM_TYP_)
#define _PSM_TYP_

#include "cmsa.h"
#include "new-dsc.h"
#include "psg.h"
// #include "gor_typ.h"
#include "histogram.h"

class psm_typ {		// position-specific-matrix type
public:
                psm_typ( ){ assert(!"Illegal constructor"); }
                psm_typ(Int4,cma_typ);	// use default settings...
                psm_typ(Int4,char *,cma_typ);
		psm_typ(Int4,char *,cma_typ,Int4);
//              psm_typ& operator=(const psm_typ&);     // assignment operator.
                ~psm_typ( ){ Free(); }
        void	ReComputeSMX(cma_typ);
        void    Put(FILE *);
        smx_typ	*GetSMX(){ return M; }
	double	RelMap(cma_typ);
	double	RelMap(Int4, cma_typ);
	void	PutAlign(FILE*, e_type,char *,Int4, Int4 , Int4 );
	void	SetRptPval(double p){ assert(p>0.0 && p<=1.0); target_pvalue=p; }
	char    *Align(e_type,Int4, Int4 *, Int4 *, Int4 *);
	char    *AlignHMM(e_type,Int4,Int4 *,Int4 *,Int4 *);
	char    *SampleAlign(e_type,Int4,Int4 *,Int4 *,Int4 *,Int4);
	char    *Align(FILE*, e_type,Int4, Int4 *, Int4 *, Int4 *);
	double	*ThresholdScore(e_type, double, Int4 *, Int4 *);
	char	*FindBestRpts(e_type, Int4, Int4 *, Int4 *, Int4 *,Int4 *);
	idp_typ	*IdpTyp( ){ return idp; }
	psg_typ	*PSG( ){ return psg; }
private:
	void		GetProb(Int4,Int4,cma_typ);
        void    	init(Int4,char*,char,Int4,cma_typ);
        // void         copy(const psm_typ& idp);
        void    	Free();         // free memory...
        unsigned short  MaxRpts, length, nblk,*lenblk;
	double		*prob_c,*prob_s,*prob_h;
	Int4		PerNats;
	idp_typ		*idp;	// insertion & deletion penalties.
	smx_typ		*M;
	Int4		**gpen;
	char		*consensus;
	double		target_pvalue;  // for a single repeat...
	psg_typ		*psg;
};

#endif
