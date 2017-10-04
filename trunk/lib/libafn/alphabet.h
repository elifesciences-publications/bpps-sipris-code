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

/****************** alphabet.h - alphabet abstract data type.***************/
/* ident: Alpha  field: A  File: alphabet.h  FileId: Alpha/ type a_type
 *************************** alphabet datastructure **********************


  char: X  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  code:	0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20

*************************************************************************/
#if !defined(ALPHA)
#define ALPHA
#include "stdinc.h"
#include "sset.h"
#include "license.h"
#include "afnio.h"
/************************** alphabet datatype *****************************/
typedef struct {
	Int4     n;			/* number of LETTERS */
	Int4     N;			/* number of symbols (includes gap)*/
	char    *alphabet;		/* ALPHABET */
	char    *code2let;		/* CODE2LETTER */
	char    *code2lower;		/* CODE2LETTER lower case */
	char    *let2code;		/* LETTER2CODE */
	char    **R;			/* relatedness scoring matrix */
	Int4	loR;			/* lowest value in R */
	Int4	hiR;			/* highest value in R */
	Int4	**matrix;		/* blast-type matrix */
	char	*prs;			/* pairs string */
	char	**pairs;		/* residue pairs */
	BooLean	*paired;		/* is residue r paired? */
	Int4	npairs;			/* number of pairs */
	float	*surface;		// surface probability (0..1)
} alphabet_type;
typedef	alphabet_type *a_type;

/******************************** PRIVATE **********************************/
void	alpha_error(const char *s,a_type A);

/******************************** PUBLIC ***********************************/
/**************************** alphabet operations **************************/
a_type	MkAlpha(const char *letters,const char *score_matrix);	/* define alphabet */
a_type	DefaultMkAlpha(const char *letters,const char *score_matrix);
a_type  MkAlpha(const char *letters,const char *score_matrix,float *surface);
a_type  MkAlphabet(const char *map_s,const char *R,const char *prs);
a_type  CopyAlphabet(a_type A0);
a_type	DefAlphaR(const char *R,a_type A);	/* redefine R */
a_type  DefAlphaP(const char *prs,a_type A);
a_type	PutAlpha(FILE *fptr,a_type  A); /* output alphabet A to file */
void	PutAlphaR(FILE *fptr,a_type  A);/* output alphabet Pairs to file */
void    PutAlphaPairs(FILE *fptr, a_type A);
void    PutAlphaSurface(FILE *fptr,a_type A);
void	NilAlpha(a_type A);		/* make alphabet A undefined */
Int4    Str2SeqAlpha(char *str, unsigned char *seq, a_type A);

char    *GetPatternFromSST(sst_typ sst, a_type AB);
char    *SST2ArgStrAlpha(sst_typ *sst, Int4 len,a_type AB);
Int4	CardSstAlpha(sst_typ sst, a_type AB);

Int4    **GetBlastMatrixAlpha(a_type A);
BooLean MemberAlpha(char c, a_type A);
sst_typ *StringToSmallSetsAlpha(char *residue_str,Int4 max_pattern_length,
                const char *Usage,Int4 *pattern_length,a_type AB);
void    PutSST2Alpha(FILE *fp,sst_typ sst,a_type AB);
void    PutSST(FILE *fp,sst_typ sst,a_type AB);
void    PutSST(char *str,sst_typ sst,a_type AB);
Int4    CardSST(sst_typ sst,a_type AB);
double  **GetRhoCategoricalPriors(FILE *fp, Int4 Length, double rho, sst_typ **sst, a_type AB);
double  *GetJthRhoCategoricalPriors(FILE *fp, double rho, sst_typ *sstJ, a_type AB);
sst_typ *StringToSmallSet(char *residue_str,Int4 max_pattern_length,
                const char *Usage,Int4 *pattern_length,a_type AB);
Int4	ParseResidueSets(char *str, Int4 *pos, sst_typ *sst,a_type AB,const char *msg);
/**************************** alphabet defines ****************************/
#define nAlpha(A)		((A)->n)
#define NAlpha(A)		((A)->N)
#define UndefAlpha(A)		0
#define AlphaR(A)		((A)->R)
#define valAlphaR(c,d,A)	((A)->R[(c)][(d)])
#define lowAlphaR(A)		((A)->loR)
#define highAlphaR(A)		((A)->hiR)
#define valAlphaP(c,d,A)	((A)->pairs[(c)][(d)])
#define	AlphaChar(i,A)		((A)->code2let[(i)])
#define	AlphaCharLow(i,A)	((A)->code2lower[(i)])
#define	AlphaCode(c,A)		((A)->let2code[(c)])
#define	MapAlphaCode(A)		((A)->let2code)
#define	MapCodeAlpha(A)		((A)->code2let)
#define MemAlpha(c,A)		((Int4) (A)->let2code[(c)])
#define NPairsAlpha(A)		((A)->npairs)
#define PairedAlpha(c,A)	((A)->paired[(c)])
#define	AlphaSurface(A)		((A)->surface)

/************ CONSTANTS *****************/
#define	ALPHA_NUM_SYMBOLS		127

#endif

