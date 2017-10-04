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

#include "sequence.h"
#include "alphabet.h"
#include "residues.h"
#include "smatrix.h"
#include "prtn_model.h"

typedef struct {
	Int4	***table;
	Int4	number,J;
	Int4	*counts,*nsize;
	a_type	A, dnaA;
	Int4	min_len,max_len;
	char	code;
	char	*filename;
	FILE	*fptr;
	Int4	num[7];
	e_type	*E[7],dnaE;
	Int4	max_seq, min_seq;
} dna2aa_type;
typedef dna2aa_type *n2a_typ;

#define	NumRFsN2A(n2a)		6
#define	nSqRFsN2A(rf, N)	(((rf) > 0) && ((rf) <= 6) ? ((N)->num[(rf)]):0)

/**************************private****************************/
Int4 ***create_table_n2a(Int4 code);
void	reset_n2a(n2a_typ n2a);
void reverse_compl_n2a(n2a_typ n2a);
BooLean	translate_n2a(BooLean rev, n2a_typ n2a);
/**************************PUBLIC*****************************/
n2a_typ	MkN2A(char *filename, Int4 code, Int4 max_seq, Int4 min_seq, 
		Int4 max_len, Int4 min_len, a_type A, a_type dnaA);
void	NilN2A(n2a_typ n2a);
e_type	SeqN2A(Int4 rf,Int4 i, n2a_typ n2a);
BooLean	ReadSeqN2A(n2a_typ n2a, BooLean rev);
