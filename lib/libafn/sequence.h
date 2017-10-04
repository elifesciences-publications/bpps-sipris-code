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

/* entity.h - entity data type. */
#if !defined(SEQUENCE)
#define SEQUENCE
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "dheap.h"
#include "random.h"
#include "histogram.h"
#include "rff_typ.h"
/**************************** Sequence ADT ***************************
		E == <I,S>	2-tuple
		I == sequence identification key
		S == biological sequence 

  Subseq = {|offset[extend]|}
  default: {|0(0)|}
  TaxSeq = {<phylum(kingdom)>}
  default: {<0(0)>}
  Both   = {|34(98)|<Chordata(M)>}
  Both   = {|34(98)|<Chordata(M11)>}	// M11 = Metazoa and genetic code 11.
**********************************************************************/
/*************************** Sequence type *****************************/
typedef struct seq_struct_type {
	UInt4	I;	/* identifier for sequence */
	Int4		n;	/* length of sequence */
	BooLean		xed;	/* sequences x'ed */
	unsigned char	*S;	/* sequence */
	unsigned char	*X;	/* if !xed X == S; else X'ed seq */
	char		*info;	/* description of entity */
	unsigned char	overlap;// maximum 'xxx' overlap on ends
	UInt4	offset;	/* offset from N-terminus */
	UInt4	extend;	/* extension beyond C-terminus */
	char		kingdom;// loosely defined.
	unsigned char	genetic_code;	// == 0 for protein.
	char		*phylum;// loosely defined.
	rff_typ		*rff;	// rich fasta format information...
	struct seq_struct_type *next;
} sequence_type;
typedef sequence_type	*e_type;

/******************************** private **********************************/
char    *copy_seq_info(e_type E);
Int4	get_diagonal_ends_seq(unsigned char *seq1, unsigned char *seq2, 
	char **R, Int4 n, Int4 *begin, Int4 *end);
e_type  EmptySeq(Int4 I, Int4 length);
e_type  AllocSeq(UInt4 I, Int4 length);
void	AllocXSeq(e_type E);
char    *get_info_seq(e_type E, char *info);
void	seq_error(const char *s);
/******************************** PUBLIC **********************************/
e_type  MkEmptySeq(Int4 I, char *info, Int4 length);
e_type  MakeSeq(char *id, char *descript, Int4 offset, Int4 extend, Int4 len,
	unsigned char *seq);
e_type  MakeSeq(char *id, char *descript, Int4 offset, Int4 extend, Int4 len,
        unsigned char *seq,char *phylum, char kingdom);
e_type  MkSeq(char *defline, Int4 len, unsigned char *seq);
e_type  MkSubSeq(Int4 start, Int4 end, e_type E);
e_type  CopySeq(register e_type E);
e_type  CopySeq(register e_type E, char *rffstr);
e_type  CopySeqPlus(e_type E, UInt4 p, UInt4 N);
void	NilSeq(e_type E);			// destructor.

// routine to check seuqence format and return the number of sequences and max length.
long    IsFastaFormatSeq(FILE *fp,long &MaxLen); 

// routings for reading input files.
e_type  ReadSeq(FILE *fptr, Int4 I, Int4 size, a_type A);
e_type  *ReadSeqFileFA(char *infile, a_type A, Int4 max_in_seq);
e_type  ReadSeqFA(char *infile, Int4 I, a_type A);
e_type  ReadGuideSeqFA(char *infile, Int4 I, BooLean **ignore, a_type A);
Int4    ReadSeqStr(FILE *fptr, register char *defline, register unsigned char *seq,
        Int4 max, a_type A);
e_type	ReadSeqStrII(FILE *fptr, UInt4 I, register char *defline, 
	register unsigned char *seq, Int4 max, a_type A);
void    ReadRichFASeq(FILE *fptr, char mode, e_type E);

// basic operations on sequences.
e_type  MergeSeqs(e_type *E);
Int4    CenterSeqHSP(Int4 offset, e_type E1, e_type E2, a_type A);
Int4    FindMaxWordSeq(Int4 q_start, Int4 s_start, Int4 N, e_type qE,
        e_type sE, Int4 *word_score, a_type A);
void    AdjustOffsetSeq(Int4 diff, e_type E);
void    AdjustCtermExtendSeq(Int4 diff, e_type E);
Int4    ConcatSeqList(e_type headE, e_type tailE);
void    SeverNextSeqList(e_type E);
void    ShortenSeqInfo(e_type E);

// Sequence comparison routines.
char	IsSubSeq(e_type E1, e_type E2);
char	IsSubSeq(e_type E1, e_type E2,Int4 *Start,BooLean IgnoreX);
char    IsSameSeq0(e_type E1, e_type E2,Int4 *Start,BooLean IgnoreX);
char    IsSameSeq(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap,BooLean IgnoreX);
char    IsSameSeqFast(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap);
char    IsSameSeqFast(e_type E1, e_type E2,Int4 *Start,Int4 *RtnNumX, Int4 MinOverlap);
char    IsSameSeqFastX(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap);
char    IsSameSeqFastX(e_type E1, e_type E2,Int4 *Start,Int4 *RtnNumX, Int4 MinOverlap);
// BooLean OverlappingSeq(Int4 *offset, Int4 MinOverlap, e_type E1, e_type E2);
// BooLean OverlappingSeq(Int4 *offset,Int4 MinOverlap,Int4 MaxMisMatch,e_type E1,e_type E2);
BooLean OverlappingSeq(Int4 &offset, e_type E1, e_type E2, Int4 MinOverlap, Int4 MaxMisMatch=0);
BooLean IdentSeqs(e_type E1, e_type E2);
BooLean FastIdentSeqs(register e_type E1, register e_type E2);
BooLean NonNullIdentSeqs(register e_type E1, register e_type E2);
BooLean SwissSeq(e_type E);
BooLean PdbSeq(e_type E);
BooLean EST_Seq(e_type E);

// Masking and unmasking routines.
void    MaskSeq(Int4 start, Int4 end, e_type E);
void    UnXSeq(e_type E);
void	AddXArraySeq(register e_type E);

// Get information routines.
Int4    GetFastaInfo(char *DBS_NAME, Int4 max, Int4 **pcounts,
        unsigned short **psize, a_type A);
Int4    GetLongFastaInfo(char *DBS_NAME,Int4 max,Int4 **pcounts,Int4 **psize,a_type A);
Int4    *NumResSeq(e_type E,a_type A);
void    NumResSeq(register e_type E,register UInt4 *num, a_type A);
Int4    NonXResSeq(e_type E,a_type A);
Int4    MaxSegSeq(e_type E);
double  *FreqResSeq(e_type E, double *freq, a_type A);
BooLean AddCountsSubSeq(UInt4 *counts,Int4 start,Int4 end,e_type E);
BooLean StringInSeqInfo(char *phrase, e_type E);
void	ChangeInfoSeq(char *new_info, e_type E);

// Output routines.
void    PutSeqID(FILE *fptr,e_type E);
void    PutShortSeqID(FILE *fptr,e_type E);
void    PutSeqID2(FILE *fptr,e_type E, Int4 len);
void    PutSeqInfo(FILE *fptr,e_type E);
void    PutSubSeqInfo(FILE *fptr,Int4 start,Int4 end, e_type E);
void    PutSeqInfo2(FILE *fptr,e_type E);
void    PutSeq(FILE *fptr,e_type E,a_type A);
void    PutSeq(FILE *fptr,Int4 subID, e_type E,a_type A);
void    PutSeq(FILE *fptr,Int4 subID, e_type E,a_type A,BooLean Rtns);
void    PutSubSeq(FILE *fptr, Int4 start, Int4 end, e_type E, a_type A);
void    PutDeleteSeq(FILE *fptr,Int4 start_del, Int4 end_del, e_type E,a_type A);
BooLean PutSuperSeq(FILE *fp, Int4 left, Int4 right, e_type sE, e_type oE,
        a_type A);
void    PutXSeq(FILE *fptr,e_type E,a_type A);
void    PutSeqMtfMask(FILE *fptr,Int4 M, Int4 *len, Int4 *p0, Int4 *p1,
        e_type E1,a_type A);
void    PutSeqRegion(FILE *fptr,Int4 start, Int4 length, e_type E, a_type A);
void    PutSeqRegion2(FILE *fptr,Int4 start, Int4 length, e_type E,
        Int4 flank, a_type A);
void    PutSeqRegionFormatMSA(FILE *fptr,Int4 start, Int4 length, e_type E,
        double p, a_type A);
Int4	PutDiagonalSeq(FILE *fptr, Int4 offset, e_type E1, e_type E2, a_type A);
sst_typ *SST_FromSeq(e_type keyE);

// Input and output of sequence information as character strings.
Int4	NumXSeq(e_type E);
void    StrSeqInfo(char *str,e_type E);
void    StrSeqID(char *str, Int4 n, e_type E);
BooLean	IsSameSeqID(e_type E1, e_type E2);
void    StrSeqDescript(char *str,Int4 n, e_type E);
Int4    SeqToString(register char *buffer, e_type E, register a_type A);
e_type  StringToSeq2(unsigned char *seq, Int4 length, char *info, UInt4 I);
e_type  StringToSeq(char *buffer, char *info, UInt4 I, a_type A);

// Randomization routines.
e_type  ReverseSeq(e_type E);
e_type	ShuffleSeq(e_type E);
e_type  ShuffleSubSeq(Int4 start, Int4 end, e_type E);
e_type  RandomSeq(Int4 length, Int4 I, double *freq, a_type A);
char    RandomResSeq(register e_type E);
unsigned char *ShuffleSeqArray(Int4 length, unsigned char *seq);

void    TaxAssignSeq(char *phylum, char kingdom, e_type E);
void    AddRFF2Seq(register e_type E, rff_typ *rff);
void	PutSeqRFF(FILE *fp,e_type E);
void	UnLabelSeq(e_type E);
void	LabelSeq(e_type E);
/**************************** Macros ********************************/
#define MAX_SEQ_DEFLINE_LENG		500
#define DEFAULT_INFO_LINE_LEN_SEQ	200
#define	GeneticCodeSeq(E)	((E)->genetic_code)
#define	KingdomSeq(E)		(toupper((E)->kingdom))
#define	kingdomSeq(E)		((E)->kingdom)
#define	PhylumSeq(E)		((E)->phylum)
#define	LabeledSeq(E)		(islower((E)->kingdom)?TRUE:FALSE)
#define	LenSeq(E)		((E)->n)
#define	NextSeq(E)		((E)->next)
#define	SetNextSeq(E,nE)	((E)->next=(nE))
#define OffSetSeq(E)		((E)->offset)
#define SetOffSetSeq(o,E)	((E)->offset = (o))
#define	CtermExtendSeq(E)	((E)->extend)
#define	SetCtermExtendSeq(o,E)	((E)->extend = (o))
#define	ResSeq(i,E)		((E)->S[(i)])
#define	IsStructSeq(E)		((E)->rff != 0 && (E)->rff->IsStruct())
#define	SeqIsRFF(E)		((E)->rff != 0)
#define	StructSeq(r,E)		(((E)->rff)?(E)->rff->Struct((r)):'u')
#define	XSeq(i,E)		((E)->X[(i)])
#define EqSeq(i,x,E)		((E)->S[(i)]=(unsigned char) (x))
#define EqXSeq(i,x,E)		((E)->X[(i)]=(unsigned char) (x))
#define RmSeq(E)		((E)->I=NULL)
#define EqSeqI(i,E)		((E)->I=(i))
#define SeqI(E)			((E)->I)
#define SeqKey(E)		((E)->info)
#define SeqPtr(E)		((E)->S)
#define XSeqPtr(E)		((E)->X)
#define XedSeq(E)		((E)->X != (E)->S)
#define SeqOverlap(E)		((E)->overlap)
#endif

