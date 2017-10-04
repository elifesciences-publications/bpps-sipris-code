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

#include "dsets.h"

ds_type	DSets(UInt4 n)
{
	ds_type	S;
	NEW(S,1,dsets_type);
	S->N = n;
	S->NSets = n;
	NEW(S->rank,(n+1),UInt4);
	NEW(S->max,(n+1),UInt4);
	NEW(S->leng,(n+1),UInt4);
	NEW(S->p,(n+1),UInt4);
	make_dsets(S);
	return (S);
}

void	make_dsets(ds_type D)
/* Create a new set containing the single element x, for each element in D */
/* put each element within it's own disjoint set */
{
	UInt4 i;

	for(i=0;i<=D->N;i++){
		D->p[i]=i;D->max[i]=i;D->rank[i]=0;
	}
}

UInt4	findDSets(UInt4 x,ds_type S)
/* Return the canonical element of the set containing element x */
/* disjoint set find() operation with path halving */
{ 
	UInt4	*p;
	p=S->p;
	while(p[p[x]]!=p[x]){p[x]=p[p[x]];x=p[x];} 
	return (p[x]); 
}

UInt4	linkDSets(UInt4 x,UInt4 y,ds_type S)
/* link(x,y): Form a new set that is the union of the two sets whose canonical
elements are x and y, destroying the two old sets.  Select and return
a canonical element for the new set.  This operation assumes that x != y. */
{
	UInt4	temp;

	if(x != y) S->NSets--;  // merging these results in one less set...
	if(S->rank[x] > S->rank[y]) { temp=x; x=y; y=temp; }
	else if(S->rank[x] == S->rank[y]) { S->rank[y]=S->rank[y] +1; }

	/* store maximum sequence in set with canonical element y */

	/* if max items are of equal length choose smaller item number */
	if(S->leng[S->max[x]]==S->leng[S->max[y]]) 
	    S->max[y]=S->max[x]=MINIMUM(UInt4,S->max[y],S->max[x]);   
	/* if not equal choose longer one */
	else if(S->leng[S->max[x]]>S->leng[S->max[y]])S->max[y]=S->max[x];
	else S->max[x] = S->max[y];

	S->p[x] = y;
	return(y);
}
	
void	PutDSets(FILE *fptr,ds_type S)
{	
	UInt4 i,j,s,set;
	for(set=i=1;i<=S->N;i++){
	   s = findDSets(i,S);   
	   if(s == i){
	      fprintf(fptr,"Set %d:\n",set); set++;
	      fprintf(fptr,"no\trank\tclass\n");
	      for(j=1;j<=S->N;j++){
	         if(s == findDSets(j,S)) PutDSet(fptr,j,S);
	      } fprintf(fptr,"\n");
	   }
	} fprintf(fptr,"  %d sets total\n",S->NSets);
}

void	PutDSet(FILE *fptr, UInt4 i,ds_type S) 
{ fprintf(fptr,"%d\t%d\t%d\n",i,S->rank[i],S->max[findDSets(i,S)]); }

void	NilDSets(ds_type S)
{ free(S->rank); free(S->max); free(S->leng); free(S->p); free(S); }

Int4	*AssignDSets(ds_type sets, Int4 **Cardinality, Int4 *NumSets)
// assign the N elements in sets to integers from 1 to NSets.
{
        Int4    i,j,s,set,N,NSets,n;
        Int4    *Set,*Card;

        NSets = sets->NSets; N = sets->N;
        NEW(Set, N+2, Int4);
        NEW(Card, NSets+2, Int4);
        for(set=0,i=1;i<=N;i++){
           s = findDSets(i,sets);
           if(s == i){
              set++;
              for(n=0,j=1;j<=N;j++) if(s==findDSets(j,sets)) { Set[j]=set; n++; }
	      Card[set]=n;
           }
        } *NumSets = set; *Cardinality = Card;
	assert(set == NSets);
        return Set;
}

Int4	*RtnOneDSet(ds_type sets, Int4 Element, Int4 &Cardinality)
// return the set containing Element
{
        Int4    i,j,s,set,N,NSets,n;
        Int4    *RtnSet,Card;

        NSets = sets->NSets; N = sets->N;
	if(Element < 1 || Element > N){ Cardinality=0; return 0; }
        NEW(RtnSet, N+2, Int4);	// set is no larger than N; can later determine exact size for space efficiency.
        s = findDSets(Element,sets);
        for(n=0,j=1;j<=N;j++) if(s==findDSets(j,sets)) { n++; RtnSet[n]=j; }
	Cardinality=n;
        return RtnSet;
}

