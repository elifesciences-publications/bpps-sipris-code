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

#if !defined(_HMA_TYP_)
#define _HMA_TYP_

#include "HMM_typ.h"

class hma_typ {		// position-specific-matrix type
public:
                hma_typ( ){ assert(!"Illegal constructor"); }
                hma_typ(Int4,cma_typ,double *);	// use default settings...
		hma_typ(Int4,char *,cma_typ,Int4,double *);
		hma_typ(Int4,char *,cma_typ,Int4,double *,BooLean);
//              hma_typ& operator=(const hma_typ&);     // assignment operator.
                ~hma_typ( ){ Free(); }
	double  SampleRelMap(double);
	double	RelMap(){ return RelMap(0); }
	double	RelMap(Int4);
	void	Put(FILE*);
	char    *Align(Int4,Int4, Int4 *, Int4 *, Int4 *);
	char    *SampleAlign(Int4,Int4,Int4 *,Int4 *,Int4 *,Int4);
	char	*FindBestRpts(Int4, Int4, Int4 *, Int4 *, Int4 *,Int4 *);
	char    SampleGSeqGibbs(Int4,BooLean*,double*,double,double,UInt4);
	cma_typ	CMSA( ){ return cmsa; }
private:
	char    	*GetBestTransProbPairs( );
        void    	init(Int4,char*,Int4,cma_typ,double *);
        // void         copy(const hma_typ& tpb);
        void    	Free();         // free memory...
	char		**operation;
	Int4		*oper_len;
	Int4		*start;
	BooLean		*skipseq;
	HMM_typ		*hmm;
	cma_typ		cmsa;
	Int4		PerNats;
	BooLean		weightMAP;
};

#endif
