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

/****************** strctprof.h ***************/
#if !defined(STRCTPROF)
#define STRCTPROF
#include "sequence.h"
#include "afnio.h"
#include "alphabet.h"
#include "histogram.h"
#include "probability.h"
#include "random.h"

/********************** structural profile datatype *************************/
typedef struct {
	char		*filename;
	e_type		E;	/** sequence of structure **/
	Int4		N;	/** length of profile **/
	a_type		A;	/** alphabet **/
	Int4		**S;	/** profile score matrix **/
	Int4		**I;	/** insertion matrix **/
	Int4		**D;	/** deletion matrix **/
	Int4		**score;/** profile scores **/
	Int4		*is;	/** insertion start penalties **/
	Int4		*ie;	/** insertion extend penalties **/
	Int4		*ds;	/** deletion start penalties **/
	Int4		*de;	/** deletion extend penalties **/
	float		**rS;	/** real profile score matrix **/
	float		**rI;	/** real insertion matrix **/
	float		**rD;	/** real deletion matrix **/
	BooLean		*null;	/** don't use these sites **/
} strct_pf_type;

typedef	strct_pf_type	*spf_typ;

/******************************** private **********************************/
void    strctprof_error(char *s);

/******************************** PUBLIC ***********************************/
spf_typ MkSProfile(double **S, e_type E, a_type A);
void	SetGapsSProfile(Int4 n, Int4 is, Int4 ie, Int4 ds, Int4 de, 
        spf_typ P);
Int4    PutAlnSProfile(e_type E, spf_typ P);
Int4    AlnSProfile(e_type E, spf_typ P);
void    PutSProfile(FILE *fptr, spf_typ P);
void    NilSProfile(spf_typ P);
spf_typ ReadSProfile(char *infile, e_type E, a_type A);
void    WriteSProfile(char *outfile, spf_typ P);
void	Gaps2ndarySProfile(Int4 is[3], Int4 ie[3], Int4 ds[3], Int4 de[3],
        char *ss, spf_typ P);
double  RealAlnSProfile(e_type E, spf_typ P);

/******************************** MACROS ***********************************/
#define SProfileA(P)		((P)->A)
#define LenSProfile(P)		((P)->N)

#endif
