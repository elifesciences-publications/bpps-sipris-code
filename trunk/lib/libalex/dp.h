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

#ifndef _DP_H_
#define _DP_H_
#include "cmsa.h"
#include "smatrix.h"
#include "routines.h"

class dp_typ{
public:
			dp_typ();
			~dp_typ();
			dp_typ(e_type E1, smx_typ M1, Int4 *inso, Int4 *inse, 
				Int4 *delo, Int4 *dele);
	char 		*TraceBackWRTpos(Int4 end_of_aln, Int4 *alignscore, 
				Int4 *start_of_aln, Int4 *oper_len);
	Int4 		*LastColumn();
	Int4		*MaxLast();
	Int4            *LocMax();
	void		PrintGaps();
private:
	Int4 		alph_len, seq_len, prof_len;
	Int4		*last_column, *max_last, *loc_max;
	char		*end_state;
	Int4		**MAT, **DEL, **INS;
	Int4 		**match;
	Int4		*io, *ie, *od, *de;
	e_type 		E;
	smx_typ 	M;

	void		Free();
	void		Mkdp_typ(e_type E1, smx_typ M1, Int4 *inso, Int4 *inse, 
				Int4 *delo, Int4 *dele);
};
	char 		**MaximizeSum(e_type E, smx_typ M1, Int4 *io, Int4 *ie, Int4 *od, Int4 *de, 
				Int4 start1, Int4 end1, Int4 start2, Int4 end2, Int4 flank, 
					Int4 *fstart1, Int4 *fstart2, Int4 *fend1, Int4 *fend2, 
						Int4 *foper_len1, Int4 *foper_len2, Int4 *falignscore1, 
							Int4 *falignscore2);
        char           **DoAlignm(e_type E, smx_typ M, Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 flank,
				Int4 cutoff, Int4 *starts, Int4 *ends, Int4 *alignscores, 
					Int4 *oper_lens, Int4 *n_rpts);

#endif
