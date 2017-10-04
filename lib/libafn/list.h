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

#if !defined (LIST)
#define LIST
#include "stdinc.h"
/************************************ list.h **********************************
 Header file for data structure representing list of small positive integers.	

list = 1 3 7 9 -

         0  1  2  3  4  5  6  7  8  9
       |-1| 3|-1| 7|-1|-1|-1| 9|-1| 0|		(-1 = not on list.)
            ^                       ^
            first                   last

 ******************************************************************************/

typedef struct {
	Int4	N;		/* list defined on ints in {1,...,N} */
	Int4	first;		/* beginning of list */
	Int4	last;		/* last element of the list */
	Int4	len;		/* number of elements in list */
	Int4	*next;		/* next[i] is successor of i in list */
} list_type;

typedef list_type *l_type;

void list_error(char *s);

/****************************** PUBLIC ***************************/
l_type	List(Int4 N);
l_type	MkList(Int4 N);
l_type	NilList(l_type L);
Int4	ItemList(Int4 i, l_type L);	/* access item i in THETA(i) time */
void	AppendList(Int4 i, l_type L);
void	RemoveList(Int4 i, l_type L);
BooLean RmList(register Int4 i, register l_type L);
void	AssignList(l_type L1, l_type L2);
void	ClearList(l_type L);		/* remove everything */
void	ResetList(Int4 N, l_type L);	/* re-initialize */
void	PutList(FILE *fptr, l_type L);

/********************** MACROS **********************************/
#define MbrList(i,L)		((L)->next[(i)] != -1)
#define	FirstItemList(L)	((L)->first)
#define	LastItemList(L)		((L)->last)
#define LengthList(L)		((L)->len)
#define SucList(i,L)		(((L)->next[(i)] == -1)?\
				(list_error("item not on list"),0):((L)->next[(i)]))

#endif

