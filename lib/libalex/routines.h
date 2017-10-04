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

#include "residues.h"
#include "alphabet.h"
#include "sequence.h"
#include "prtn_model.h"
#include "smatrix.h"

Int4 *ShuffleArray(Int4 len_a);
e_type ConsensusSeq(smx_typ *M, Int4 nbr, Int4 rpts);
e_type ConsensusSeq(smx_typ M);
void PutAlign(e_type E1, smx_typ *M, Int4 nbl, char *operation, Int4 operation_length,
                        Int4 start, Int4 rpts);
void PutAlign(e_type E1, smx_typ M, char *operation, Int4 operation_length, Int4 start, Int4 score);
char *ReverseOperArray(char *operation1, Int4 oper_len1);
Int4 *ReversePen(Int4 *gap, Int4 len);
char *ConcatenateOperArrays(char **opers, Int4 n_opers, Int4 *oper_lens, 
Int4 *starts, Int4 *ends, Int4 seq_len, Int4 *operation_length);
Int4 CompareOperArr(char *oper1, char *oper2, Int4 *pos, 
                Int4 num_rpts1, Int4 num_rpts2, Int4 rpt_len);
Int4 CompareRpts(Int4 *p, Int4 *q, Int4 rpt_len, Int4 *pos);
char *ExtendOperArray(char *operation, Int4 operation_length, Int4 start);
