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

/************ mlist.h - multiple linked list abstract data type.***********/
#if !defined(MLIST)
#define MLIST
#include "stdinc.h"
/************************** alphabet datatype *****************************/
typedef struct {
	Int4     M;			/* number of elements for lists */
	Int4     N;			/* number of lists */
	Int4	*first;			/* pointers to beginning of lists */
	Int4	*item;			/* items on lists */
	Int4	*next;			/* pointer to next element */
	Int4	free;			/* pinter to free list */
} multi_list_type;
typedef	multi_list_type *ml_type;

/******************************** PRIVATE **********************************/
void	mlist_error(char *s);

/******************************** PUBLIC ***********************************/
/**************************** mlist operations **************************/
ml_type	MkMList(Int4 M, Int4 N);
void	Add2MList(register Int4 i, register Int4 n, register ml_type L);
Int4	GetListMList(register Int4 *list, register Int4 n, 
        	register ml_type L);
void	NilMList(ml_type L);

/***************************************************************************/
#define EmptyMList(n,L)	((L)->first[(n)]==0)
#define SizeMList(L)	((L)->M)

#endif

