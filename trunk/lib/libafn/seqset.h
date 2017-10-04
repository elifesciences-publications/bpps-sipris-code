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

/***************** seqset.h - sequence set data type. ********************/
/*************************** Sequence Set ADT ******************************

	Define: P = <S,A>
	  where
		S is a set of sequences
		A is the sequence alphabet 

**************************************************************************/
#if !defined (SEQSET)
#define SEQSET
#include "stdinc.h"
#include "afnio.h"
#include "pseg.h"
#include "sequence.h"
#include "alphabet.h"
#include "random.h"
#include "histogram.h"
/******************************* PRIVATE ***********************************/
typedef struct seq_link {
	e_type	E;
	struct seq_link *next;
} sqlnk_type,*sqlnk_typ;

/****************************** seqset type ****************************/
typedef struct {
	char		*name;		/* input filename for entities */
	e_type		*entity;	/* array of sequence entities */
	Int4		nent;		/* number of input entities */
	Int4		max_leng;	/* sequence entity maximumlength */
	Int4		max_num_seqs;	/* maximum number of sequences */
	Int4		min_leng;	/* sequence entity minimum length */
	BooLean		xnu;		/* has data been x'ed? */
	Int4		*counts;	/* number of each residues */
	Int4		total;		/* total # residues */
	double		*tfreq;		/* residue total freqs for seqs */ 
	a_type		A;		/* sequence alphabet */
	double  	indel_penalty;  // log-odds penalty in nats (default=0.0)
} seqset_type;
typedef seqset_type	*ss_type;

/******************************* private *************************************/
void	seqset_error(const char *s);
ss_type  calcseqsetfreq(ss_type P);
ss_type  xnu_seqset(ss_type P);
Int4	count_seqset_entities(FILE *fptr, ss_type P,Int4 nsize[]);
FILE    *OpenSeqSetFile(ss_type P);
ss_type  seqset(char *filename,a_type A);
ss_type  fptr_seqset(FILE *fptr,a_type A);
ss_type  SeqSet_fptr(FILE *fptr,a_type A);

/******************************* PUBLIC *************************************/
/************************** seqset operations ***************************/
ss_type	Array2SeqSet(e_type *EList, Int4 N, char *name,a_type A);
ss_type MakeSeqSet(char *name,Int4 max_len, a_type A);
ss_type MkSeqSet(char *name,a_type A);
ss_type PSegSeqSet(ss_type data);
ss_type SeqSet(char *name,a_type A);
ss_type SeqSet1(char *filename,e_type E,a_type A);
ss_type MkXSeqSet(char *filename,a_type A);
ss_type CopySeqSet(ss_type P);
void	ReNumberSeqSet(ss_type P);
void	ReverseSeqSet(ss_type P);
ss_type MergeSeqSets(ss_type data1, ss_type data2);
e_type	ReplaceSeqSet(Int4 n, e_type E, ss_type P);
ss_type RmSeqSet(e_type E, ss_type P);
ss_type RemoveSeqSet(Int4 n, ss_type P);
ss_type AddSeqSet(e_type E, ss_type P);
void    RenameSeqSet(char *name, ss_type P);
ss_type NilSeqSet(ss_type P);
e_type  *NilSeqSetRtnSeqs(ss_type P);
double  SeqSetEntropy(ss_type P);
Int4    *LengthsSeqSet(ss_type P);
ss_type	PutSeqSet(FILE *fptr,ss_type P);		/* show seqset data */
ss_type PutSeqSetPIDs(FILE *fptr, ss_type P);
ss_type	PutSeqSetEs(FILE *fptr,ss_type P);
void	PutLengthsSeqSet(FILE *fptr,ss_type P);
void    PutTrimmedSeqSet(FILE *fptr, Int4 trim, ss_type P);
ss_type	PutSeqSetE(FILE *fptr,Int4 i, ss_type P);
ss_type	PutSeqSettFreqs(FILE *fptr,ss_type P);
ss_type	ShuffleSeqSet(ss_type P); 		/* shuffle seqset */
ss_type ShuffleSeqSet2(ss_type P);
double  LogL0SeqSet(ss_type P);
void    SetIndelPenaltySeqSet(double penalty, ss_type P);
BooLean *AreInSeqSet(ss_type data, ss_type keydata);
BooLean IsInSeqSet(ss_type data, char *seqid);
BooLean IsInSeqSet(e_type E, ss_type data);
Int4    AveSeqSeqSet(ss_type data);
/************************** sequence set defines ***************************/
#define MAX_NUMBER_SEQS		20000000
#define MAX_NUMBER_RESIDUES	2000000000
#define IndelPenaltySeqSet(P)	((P)->indel_penalty)
#define SeqSetA(P)		((P)->A)
#define nLetSeqSet(P)		nAlpha((P)->A)
#define CountsSeqSet(r,P)	((P)->counts[(r)])
#define TotalSeqSet(P)		((P)->total)
#define	tFreqSeqSet(P)		((P)->tfreq)
#define SeqP(n,i,P)		ResSeq(i,(P)->entity[(n)])
#define NSeqsSeqSet(P)		((P)->nent)
#define SqLenSeqSet(n,P)	LenSeq((P)->entity[(n)])
#define	SeqSeqSet(n,P)		SeqPtr((P)->entity[(n)])
#define	XSeqSeqSet(n,P)		XSeqPtr((P)->entity[(n)])
#define SeqSetE(n,P)		((P)->entity[(n)])
#define NameSeqSet(P)		((P)->name)
#define MinSeqSeqSet(P)		((P)->min_leng)
#define MaxSeqSeqSet(P)		((P)->max_leng)
#define	CntsSeqSet(P)		((P)->counts)
#define	XSeqSet(P)		xnu_seqset(P)
#define	XedSeqSet(P)		((P)->xnu)

#endif

