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

#if !defined (OLIST)
#define OLIST
#include "stdinc.h"
/**************************  ADT olist ********************************
 Header file for data structure representing ordered list of small 
positive integers.

	list: x  x  x   ... x   where x < x < x < ... < x
               1  2  3       n         1   2   3         n

	and 0 < n <= N.	(note: next[x ] = x   .)
				     i     i+1
**********************************************************************/

typedef struct {
	Int4	N;		/* list defined on ints in {1,...,N} */
	Int4	n;		/* number of elements in list */
	Int4	*next;		/* next[i] is successor of i in list */
} orderlist_type;

typedef orderlist_type *ol_type;

void olist_error(char *s);
/****************************** PUBLIC ***************************/
ol_type	Olist(Int4 N);
Int4	RmOlist(Int4 i, ol_type L);
Int4	InsertOlist(Int4 i, ol_type L); /* access item i in THETA(i) time */
void	PutOlist(FILE *fptr, ol_type L);/* print item i on list */
void    GetOlist(register Int4 *array, register ol_type L);
void	NilOlist(ol_type L);

/********************** MACROS **********************************/

/* Return number of elements in list */
#define LengthOlist(L)	((L)->n)
#define ClearOlist(L)	((L)->next[0]=0)

#endif
