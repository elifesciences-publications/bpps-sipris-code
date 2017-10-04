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

#if !defined(_HSW_TYP_)
#define _HSW_TYP_

#include "cmsa.h"
#include "alphabet.h"
#include "set_typ.h"

typedef struct henikoff_wts_type {         // henikoff sequence weighting type
	double	**Weight;	// position specific sequence weights Weight[pos][sq];
	double	**WtCnts;	// Weighted counts at each position.
	double	**WtFreq;	
	double	*fract_seq_aln;	// fract_seq_aln[pos]
	Int4	Length;
	Int4	NWtSq;
	cma_typ	mcma;
	a_type	AB;
};

typedef henikoff_wts_type *hsw_typ;
hsw_typ GetSubHSW(hsw_typ HSW,set_typ Set, cma_typ cma, cma_typ mcma);
hsw_typ AddRandomHSW(hsw_typ HSW, cma_typ mcma,cma_typ rcma, cma_typ ccma);
void    FWriteHSW(FILE *fp,hsw_typ hsw);
hsw_typ FReadHSW(FILE *fp,a_type AB, cma_typ cma);
void    NilHSW(hsw_typ HSW);

#define	WeightsHSW(w)	((w)->Weight)
#define	WtCntsHSW(w)	((w)->WtCnts)
#define	WtFreqsHSW(w)	((w)->WtFreq)
#define	LengthHSW(w)	((w)->Length)

#if 0
class sqw_typ {
public:
	sqw_typ(){ assert(!"Illegal constructor"); }
	sqw_typ(double **weight, double **wtCnts, double **wtFreq; double *FractSeqAln,
		Int4 Length, Int4 NWtSq, cma_typ mcma, a_type AB){
		Weight;	// position specific sequence weights
	WtCnts;	// Weighted counts at each position.
	WtFreq;	
	fract_seq_aln;
	Length;
	NWtSq;
	mcma;
	AB;
	}
	hsw_typ *GetSubHSW(hsw_typ *HSW,set_typ Set, cma_typ cma);
	void    NilHSW(hsw_typ HSW)

private:
	double	**Weight;	// position specific sequence weights
	double	**WtCnts;	// Weighted counts at each position.
	double	**WtFreq;	
	double	*fract_seq_aln;
	Int4	Length;
	Int4	NWtSq;
	cma_typ	mcma;
	a_type	AB;
};
#endif

#endif

