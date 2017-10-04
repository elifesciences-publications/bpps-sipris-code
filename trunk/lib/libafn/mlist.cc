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

#include "mlist.h"

ml_type MkMList(Int4 M, Int4 N)
/******************************************************************
  N = number of lists; M = number of elements.
  sets all N lists to empty and creates the list of free cells.
 ******************************************************************/
{
	ml_type	L;
	Int4	i;

	MEW(L,1,multi_list_type);
	L->M = M; L->N = N;
	NEW(L->first,N+2,Int4);	/* all set to zero = empty lists */
	MEW(L->item,M+2,Int4); MEW(L->next,M+2,Int4);
	for(L->free=1, i=1; i < M; i++){
		L->next[i]=i+1;
	}
	L->next[i] = 0;
	return L;
}

void	Add2MList(register Int4 i, register Int4 n, register ml_type L)
/******************************************************************
 add item i to the nth list.
 ******************************************************************/
 {
	register Int4	New,t;

	if(n > L->N || n <= 0) mlist_error("list number out of range");
	if((New = L->free) == 0) mlist_error("out of memory");
	L->free = L->next[New]; 
	t = L->first[n];		/* first[n] = t[?:?]->  */
	L->item[New] = i;
	L->next[New] = t;		/** New[item:t]-> t[?:?]-> **/
	L->first[n] = New; 		/** first[n] = New[i:t]-> t[?:?]-> */
}

Int4	GetListMList(register Int4 *list, register Int4 n, 
	register ml_type L)
/******************************************************************
 convert list to contain the nth list in L; returns list length.
 WARNING: assumes list[ ] is Int4 enough.
 ******************************************************************/
{
	register Int4	c;
	
	if((c=L->first[n])){
	   n=0; do{ list[n++]=L->item[c]; } while((c=L->next[c]));
	   return n;
	} return 0;
}

void    NilMList(ml_type L)
{
	free(L->first);
	free(L->item);
	free(L->next);
	free(L);
}

void    mlist_error(char *s){ fprintf(stderr,"mlist: %s\n",s); exit(1); }
