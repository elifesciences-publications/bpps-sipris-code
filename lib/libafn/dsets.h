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

#if !defined(DSETS)
#define DSETS
#include "stdinc.h"

/************************ disjoint set data type **********************/
typedef struct {
	UInt4	*rank;
	UInt4	*p;
	UInt4	*max;
	UInt4	*leng;
	UInt4	N;
	UInt4	NSets;
} dsets_type;

typedef	dsets_type	*ds_type;

/*************************** PRIVATE *********************************/
void	make_dsets(ds_type D);

/*************************** PUBLIC *********************************/
ds_type	DSets(UInt4 n);
UInt4	findDSets(UInt4 x,ds_type S);
UInt4	linkDSets(UInt4 x,UInt4 y,ds_type S);
void	PutDSets(FILE *fptr,ds_type S);
void	PutDSet(FILE *fptr, UInt4 i,ds_type S);
void    NilDSets(ds_type S);
Int4    *AssignDSets(ds_type sets, Int4 **Cardinality, Int4 *NumSets);
Int4    *RtnOneDSet(ds_type sets, Int4 Element, Int4 &Cardinality);

/*************************** MACROS *********************************/
#define	pDSets(i,S)		((S)->p[(i)])
#define	NotUnitDSets(i,S)	((S)->rank[(i)])
#define maxDSets(i,S)		((S)->max[(i)])
#define NumDSets(S)		((S)->N)
#define NSetsDSets(S)		((S)->NSets)
#define EqlengDSets(L,i,S)	((S)->leng[(i)]=(L))
#define CanonicalDSet(i,S)	((maxDSets(findDSets(i,S),S))==(i))

#endif

