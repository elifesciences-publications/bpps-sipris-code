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
#if !defined(SWALN)
#define SWALN
#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "blosum62.h"

/******************************** private **********************************/

/******************************** PUBLIC **********************************/
Int4    AlnSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A);
Int4    AlnSeqSW(FILE *fptr,Int4 a, Int4 b, e_type E1, e_type E2, a_type A);
Int4    FastAlnSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A);
float   RealAlnSeqSW(float a, float b, e_type E1, e_type E2);
Int4	SubSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A,
	Int4 *start, Int4 *end);
Int4	SubSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A,
	Int4 *start, Int4 *end,Int4 *startS);

#endif
