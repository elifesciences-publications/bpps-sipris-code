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

#include "olist.h"

ol_type	Olist(Int4 N) 
/* Create an ordered list of zero items. */
{
	ol_type L;
	NEW(L,1,orderlist_type);
	L->N = N; 
	L->n = 0;
	assert(N > 0);
	// fprintf(stderr,"*** N=%d\n",N);
	NEW(L->next,N+1,Int4);
	L->next[0] = 0;		/* 0 = null */
	return L;
}

void	NilOlist(ol_type L) { free(L->next); free(L); }

Int4	RmOlist(Int4 i, ol_type L)
/* Remove item i from list L.  */
{
	Int4	x=0;
	do{
		if( L->next[x] == i){
			L->next[x] = L->next[L->next[x]];
			L->n--;
			return i;
		} else x =  L->next[x];
	} while( x != 0);
	return 0;
}

Int4	InsertOlist(Int4 i, ol_type L)
/* Add item i to list L in theta(i) time. */
{
	Int4 x,t;
	if(i > L->N) olist_error("item larger than N");
	for(x=0; L->next[x] != 0; x =  L->next[x]) {
	   if(L->next[x] >= i){
		if(L->next[x] == i) return 0;
		break;
	   }
	}
	t = L->next[x]; L->next[x] = i; L->next[i] = t;
	L->n++;
	return i;
}

void	GetOlist(register Int4 *array, register ol_type L)
/*  put the items into the array (starting from i=1) as ordered on list L. */
{
	register Int4 x,i;
	for(x=0,i=1; L->next[x] != 0; x=L->next[x],i++)
			array[i]=L->next[x];
	array[i] = 0;
}

void PutOlist(FILE *fptr, ol_type L)
/* Print the list L. */
{
	Int4 i;
	for(i=0; L->next[i] != 0; i=L->next[i])
		fprintf(fptr,"%d ", L->next[i]);
}

void olist_error(char *s){fprintf(stderr,"olist: %s\n",s); exit(1); }


