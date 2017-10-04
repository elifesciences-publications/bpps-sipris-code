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

#if !defined(_SSEQ_)
#define _SSEQ_

/************* Sampling Block Alignment type ***************/
#include "prtn_model.h"
#include "smatrix.h"
#include "sequence.h"
#include "residues.h"
#include "alphabet.h"

class ssq_typ 
{
public:
		ssq_typ();
		ssq_typ(ptm_typ PM, double pernats,Int4 max_gap_len);
	void	make_ssq(ptm_typ PM, double pernats,Int4 max_gap_len);
	e_type 	sample_ssq(Int4 I);
	e_type  sample_ssq(Int4 I,Int4 block,Int4 column,Int4 inslen, double *left_relen, double *right_relen);
	void	Free();
		~ssq_typ();

private:
	double	**gpen, **match;
	Int4	*block_lengths, prof_len, nmbr_of_blocks,*start_prof, MaxGap, alph_len;
};

#endif

