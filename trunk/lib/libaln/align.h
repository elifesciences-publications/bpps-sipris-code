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

#if !defined (_ALIGNMENT_)
#define	_ALIGNMENT_
#include "afnio.h"
#include "stdinc.h"
#include "histogram.h"
#include "probability.h"
#include "mheap.h"
#include "msites.h"
#include "seqset.h"
#include "segment.h"
#include "pattern.h"
#include "profile.h"
/*************************** ADT ALIGNMENT ***************************
	
**********************************************************************/

/*************************** ALIGN type **************************/
typedef struct {
	ss_type		P;		/* aligned segment population */
	Int4		leng;		/* segment length */
        s_type          *S;       	/* pointer to array for segments */
        s_type          *maybe;       	/* array for possible segments */
        s_type          *tmp;       	/* temporary array for segments */
        s_type          *nmaybe;       	/* array for new possible segments */
        s_type          *omaybe;       	/* array for old possible segments */
        Int4             n;           	/* number of segments in alignment */
        Int4             nseg;           /* number of segments in population */
	Int4		N;
        Int4             *scores;        /* scores for refined alignment */
	mh_type		H;		/* min-max dheap */
	pfl_typ		M;		/* profile of alignment */
	BooLean		refined;
	a_type		A;		/* alphabet */
} align_type;
typedef align_type *aln_typ;

/******************************* PRIVATE *****************************/
BooLean XSegCharAlign(Int4 j, char *r, Int4 i, aln_typ Aln);
Int4     NumSeqSegLAlign(aln_typ Aln);
Int4     AddProfileAlign(s_type *fS,double cut,aln_typ Aln,BooLean fit);
Int4     RmProfileAlign(s_type *fS,s_type *mS, double cut, aln_typ Aln);
Int4	ProbProfileAlign(s_type *mS, aln_typ Aln, double cut);
void    PutPatternAlign(FILE *fptr, ptn_typ Q, aln_typ Aln);
Int4     put_align2(FILE *fptr, aln_typ Aln, BooLean verbose);
Int4     put_align(FILE *fptr, aln_typ Aln, BooLean verbose, ptn_typ Q);
Int4     FinalProbAlign(aln_typ Aln);
void	align_error(const char *s);

/******************************* PUBLIC ******************************/
aln_typ CreateAlign(ss_type P, Int4 length, s_type *S);
aln_typ *MkAlignsCombine(Int4 *numGrps, ptn_typ *motif, s_type **list,
        double hg_cut, ss_type P, double rcut, Int4 N);
BooLean SeqInAlign(short I, aln_typ Aln);
Int4     PutAlign(FILE *fptr, aln_typ Aln, BooLean verbose, ptn_typ Q);
Int4	RefineAlign(aln_typ Aln, double cut);
Int4     ScanfileAlign2(aln_typ Aln, BooLean verbose, FILE *snfptr);
Int4     ScanfileAlign(aln_typ Aln, BooLean verbose, FILE *snfptr);
void	NilAlign(aln_typ Aln);

/**************************** macro operations **********************/
#define	SegsAlign(A)		((A)->S)
#define	EmptyAlign(A)		((A)->S[0] == NULL)

#endif

