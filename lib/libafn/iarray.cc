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

#include "iarray.h"

void	BubbleSortIArray(Int4 *L)
{
	Int4 i,j,t,n;

	for(n=i=0; L[i]!= -1; i++) n++;
	for(i=1; i < n; i++){
		for(j=i;j > 0 && (L[j-1] > L[j]);j--){
			t=L[j]; L[j]=L[j-1]; L[j-1]=t;
		}
	}
}

Int4	MaxIntIArray(Int4 *I, Int4 *L)
/** sets *I to integer repeated most often; returns number of repeats **/
{
	Int4	i,N,maxI,n,x;

	BubbleSortIArray(L);
	i=0; maxI=-1; N = 0;
	while(L[i] != -1){
	   for(x=L[i], n=0; L[i] == x; ){
		n++; i++;
	   }
	   if(n > N) { N= n; maxI = x; }
	} 
	*I = maxI;
	return N;
}

Int4	RmCommonIArray(Int4 *LC,Int4 *L1)
/*  removes numbers from L1 that are present on LC */
{
	Int4 i,j,d1;

	for(i=j=d1=0; L1[i]!= -1 && LC[j]!= -1; ){
		if(L1[i] < LC[j]){
			L1[d1++] = L1[i++];
		}else if(L1[i] > LC[j]) j++;
		else { i++; j++; }
	}
	while(L1[i]!= -1) L1[d1++]=L1[i++];
	L1[d1]= -1; 
	return d1;
}

Int4	CommonIArray(Int4 *LC,Int4 *L1,Int4 *L2)
/*  removes numbers from L1 and L2 that they have in common and puts 
    them on LC */
{
	Int4 i,j,c,d1,d2;

	for(i=j=c=d1=d2=0; L1[i]!= -1 && L2[j]!= -1; ){
		if(L1[i] < L2[j]){
			L1[d1++] = L1[i++];
		}else if(L1[i] > L2[j]){
			L2[d2++] = L2[j++];
		} else {
			LC[c++] = L1[i]; i++; j++;
		}
	}
	while(L1[i]!= -1) L1[d1++]=L1[i++];
	while(L2[j]!= -1) L2[d2++]=L2[j++];
	LC[c]= -1; L1[d1]= -1; L2[d2]= -1;
	return c;
}

Int4	OffsetIArray(Int4 os,Int4 *L1,Int4 *L2)
/***************************************************************
  set L2 to be L1 with offset os.

	offsetIArray(2,[2,7,9,12,-1],L2) -> L2 = [4,9,11,14,-1]
  or
	offsetIArray(-2,[2,7,9,12,-1],L2) -> L2 = [0,5,7,10,-1]
	
  but   offsetIArray(-3,[2,7,9,12,-1],L2) ->
		iarray_error("attempt to shift numbers below 0")

***************************************************************/
{
	Int4	i;

	if((L1[0]+os)<0)iarray_error("attempt to shift number below 0");
	for(i=0; L1[i]!= -1;i++){ L2[i] = L1[i] + os;}
	L2[i] = -1;
	return i;
}

Int4	OffsetIArray2(Int4 os,Int4 *L)
{
	Int4	i;

	if((L[0]+os)<0)iarray_error("attempt to shift number below 0");
	for(i=0; L[i]!= -1;i++){ L[i] += os;}
	return i;
}

Int4	CopyIArray(Int4 *LC,Int4 *L)
{
	Int4	i;
	for(i=0; L[i]!= -1;i++){ LC[i] = L[i];}
	LC[i] = -1;
	return i;
}

Int4	IntersectIArrays(Int4 *LC,Int4 *L1,Int4 *L2)
/* set LC to be the intersection of L1 and L2; returns cardinality of LC */
{
	Int4 i,j,c;

	for(i=j=c=0; L1[i]!= -1 && L2[j]!= -1; ){
		if(L1[i] < L2[j]) i++; 
		else if(L1[i] > L2[j]) j++; 
		else { LC[c++] = L1[i]; i++; j++; }
	}
	LC[c]= -1; 
	return c;
}

Int4     UnionIArrays(Int4 *LC,Int4 *L1,Int4 *L2)
/* set LC to be the union of L1 and L2; returns cardinality of LC */
{
	Int4 i,j,c;

	for(i=j=c=0; L1[i]!= -1 && L2[j]!= -1; ){
	   if(L1[i] < L2[j]) LC[c++] = L1[i++];
	   else if(L1[i] > L2[j]) LC[c++] = L2[j++];
	   else { LC[c++] = L1[i]; i++; j++; }
	}
	while(L1[i]!= -1) LC[c++]=L1[i++];
	while(L2[j]!= -1) LC[c++]=L2[j++];
	LC[c]= -1; 
	return c;
}

void	PutIArray(FILE *fptr, Int4 *L)
{
	Int4	i;
	fprintf(fptr,"\n\t[");
	for(i=0; L[i]!= -1;i++){ 
		fprintf(fptr," %d",L[i]);
		if(i%10==9 && L[i+1] != -1) fprintf(fptr,"\n\t");
	}
	fprintf(fptr," ]\n");
}

void    iarray_error(char *s)
{ fprintf(stderr,"\tIArray error: %s\n\n",s); exit(1); }

