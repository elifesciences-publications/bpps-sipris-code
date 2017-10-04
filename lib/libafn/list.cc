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

#include "list.h"

l_type List(Int4 N) { return MkList(N); }

l_type MkList(Int4 N) {
	l_type L;
	Int4	i;

	NEW(L,1,list_type);
	L->N = N;
	NEW(L->next,N+2,Int4);
	L->first = L->last = 0;
	L->len = 0;
	for(i = 1; i <= N; i++) L->next[i] = -1;
	L->next[0] = 0;
	return L;
}

l_type NilList(l_type L) { free(L->next); free(L); return NULL; }

/* Change size of list. Discard old value. */
void ResetList(Int4 N, l_type L)
{
	Int4 i;
	L->N = N;
	free(L->next);
	NEW(L->next,N+1,Int4);
	L->first = L->last = 0;
	L->len = 0;
	for(i = 1; i <= N; i++) L->next[i] = -1;
	L->next[0] = 0;
}
	
/* Remove all elements from list. */
void ClearList(l_type L)
{
	Int4 i;
	while (L->first != 0) {
		i = L->first; L->first = L->next[i]; L->next[i] = -1;
	}
	L->last = 0; L->len = 0;
}

/* Return the i-th element, where the first is 1. */
Int4 ItemList(Int4 i, l_type L)
{
	Int4 j;
	if (i == 1) return L->first;
	if (i == L->len) return L->last;
	for (j = L->first; j != 0 && --i; j = L->next[j]) {}
	return j;
}

/* Add i to the end of the list. */
void AppendList(Int4 i, l_type L)
{
	if(L->next[i] != -1) list_error("AppendList: item already in list");
	if(L->first == 0) L->first = i;
	else L->next[L->last] = i;
	L->next[i] = 0; L->last = i; L->len++;
}

BooLean	RmList(register Int4 i, register l_type L)
/************************************************************************
 Remove element i from list. 
 ************************************************************************/
{
	register Int4	prev,nxt;

	if(L->next[i] == -1) return FALSE;
	if(i == L->first){
		L->first = L->next[L->first];
		if(L->first == 0) L->last= 0;
	} else {
	  for(prev=L->first; TRUE; prev=nxt) {
	    nxt = L->next[prev];
	    if(nxt==i){
		L->next[prev]= L->next[i];
		if(i == L->last) L->last = prev; 
		break;
	    }
	  }
	}
	L->next[i] = -1; L->len--;
	return TRUE;
}

void RemoveList(Int4 i, l_type L)
/* Remove the first i items. */
{
	Int4 f;
	while(L->first != 0 && i--) {
		f=L->first; L->first=L->next[f]; L->next[f]=-1; L->len--;
	}
	if(L->first == 0) L->last = 0;
}

/* Assign value of right operand to left operand. */
void AssignList(l_type L1, l_type L2)
{
	Int4 i;

	if (L1->N < L2->N) {
		L1->N = L2->N;
		free(L1->next);
		NEW(L1->next,L2->N + 1, Int4);
		L1->first = L1->last = 0;
		for(i = 1; i <= L1->N; i++) L1->next[i] = -1;
		L1->next[0] = 0;
	} else ClearList(L1);
	for(i = ItemList(1,L2); i != 0; i = SucList(i,L2))
		AppendList(i,L1);
}

void PutList(FILE *fptr, l_type L)
{
	Int4 i;
	for(i=L->first; i != 0; i=L->next[i])fprintf(fptr,"%d ",i);
}

void list_error(char *s){fprintf(stderr,"list: %s\n",s); exit(1); }
