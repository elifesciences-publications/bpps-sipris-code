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

#if !defined (_PSIALN)
#define	_PSIALN
#include "stdinc.h"
#include "afnio.h"
/*************************** ADT ALIGNMENT ***************************
	
**********************************************************************/

/*************************** ALIGN type **************************/
typedef struct{
    char *seq;
    Int4 *res;
    Int4 start;
    Int4 end;
    Int4 pos;
}alnseq;

typedef struct{
    Int4	algn_length;
    Int4	seq_nmbr; 
    alnseq	*align;
    Int4	master;
} psialn_type;
typedef psialn_type *psi_typ;

/******************************* PRIVATE *****************************/

/******************************* PUBLIC ******************************/
psi_typ MkPsiAln(char *filename);
void    PutPsiAln(FILE *fp, psi_typ A);
void    NilPsiAln(psi_typ A);

Int4	LengMasterPsiAln(psi_typ A);
Int4	MasterIDPsiAln(psi_typ A);

BooLean InSeqPsiAln(Int4 s, Int4 i, psi_typ A);
// returns TRUE if position i in seq s is in align 'a'

Int4    MapSeqPsiAln(Int4 pos1, Int4 s2, psi_typ A);
// returns the residue position in s2 that is aligned with pos1 of 
// the master sequence of align 'a' in PSI-ALIGN A.
//  returns NULL if pos1 is not aligned with sequence s2


/**************************** macro operations **********************/

#endif

