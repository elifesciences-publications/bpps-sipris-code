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

/* pairaln.h - pairwise alignment methods . */
#if !defined(PAIRALN)
#define PAIRALN
#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "mheap.h"

/******************************** private **********************************/
BooLean relate_seq(register unsigned char *seq1, register unsigned char *seq2,
        register char **R, register Int4 n, register Int4 cutoff);
Int4     diagonal_score_seq(register unsigned char *seq1, 
	register unsigned char *seq2, register char **R, register Int4 n);
Int4     get_diagonal_ends_seq(unsigned char *seq1, unsigned char *seq2, 
	char **R, Int4 n, Int4 *begin, Int4 *end);
Int4     repeat_score_seq(register unsigned char *seq, register char **R,
                register Int4 o, register Int4 n);
/******************************** PUBLIC **********************************/
BooLean RelatedSeqsDiag(Int4 cutoff, unsigned char *seq1, unsigned char *seq2, 
	Int4 len, a_type A);
BooLean RelatedSeqs(Int4 cutoff, e_type E1, e_type E2, a_type A);
Int4	AlignSeqFastp(e_type E1, e_type E2, a_type A);
BooLean RelateSeqFastp(e_type E1, e_type E2, a_type A,Int4 score);
Int4	RepeatScoreSeq(e_type E, a_type A, Int4 *offset);
BooLean RelatedSeqsDiag(Int4 cutoff,unsigned char *seq1,unsigned char *seq2,
	register Int4 len, a_type A);
Int4    PutRepeatsSeq(FILE *fptr,e_type E, a_type A, Int4 cutoff);

/******************************** MACROS **********************************/


#endif
