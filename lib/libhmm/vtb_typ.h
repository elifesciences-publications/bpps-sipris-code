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

#if !defined(_VTB_TYP_)
#define _VTB_TYP_

#include "sequence.h"
#include "stdlib.h"
#include "shm_typ.h"
#include "hmm_math.h"

class vtb_typ {
public:
	vtb_typ() { assert(!"illegal constructor for vtb_typ()\n"); }
	vtb_typ(Int4 MaxProfLen, Int4 MaxSeqLen);
	~vtb_typ(){ Free(); }
	Int4    Viterbi(e_type sE,shm_typ *shm);
	char	*Viterbi(e_type sE,shm_typ *shm,Int4 *num_rpts,Int4 *total_score,Int4 *scores);
	char 	*Traceback(e_type sE,shm_typ *shm,Int4 *oper_length,Int4 *num_rpts,
        		Int4 *total_score, Int4 *scores);
private:
	Int4 	FastViterbi(e_type sE,Int4 prof_len,Int4 **mat_emit,Int4 **ins_emit,
			Int4 *m2m,Int4 *m2i,Int4 *m2d,Int4 *i2m,Int4 *i2i,Int4 *d2m,
			Int4 *d2d,Int4 *b2m,Int4 *m2e,Int4 n2b,Int4 n2n,Int4 e2c,
			Int4 e2j,Int4 c2t,Int4 c2c,Int4 j2b,Int4 j2j);
	Int4 	Viterbi(e_type sE,Int4 prof_len,Int4 **mat_emit,Int4 **ins_emit,
			Int4 *m2m,Int4 *m2i,Int4 *m2d,Int4 *i2m,Int4 *i2i,Int4 *d2m,
			Int4 *d2d,Int4 *b2m,Int4 *m2e,Int4 n2b,Int4 n2n,Int4 e2c,
			Int4 e2j,Int4 c2t,Int4 c2c,Int4 j2b,Int4 j2j);
	void	init(Int4 MaxProfLen, Int4 MaxSeqLen);
	void	Free();
	Int4	maxProfLen,maxSeqLen;
	Int4	**M,**D,**I;
	char 	**traceM,**traceD,**traceI;	
	char 	*traceJ,*traceB,*traceC,*traceN;	
	Int4	*N,*B,*E,*C,*J;
	Int4 	*jumpE;
Int4 *tempLong;
char *tempChar;
};

#endif
