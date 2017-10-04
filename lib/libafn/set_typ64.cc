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

#include "set_typ.h"
UInt4 *SET_BIT_SETS=NULL;
unsigned char *CARDINALITY_SETS=NULL;
unsigned char *BIT_NUMBER_SETS=NULL;

void	initialize_sets(void)
/* initialize lookup tables used for bitwise operations */
{
	Int4 i,b,n;
	unsigned char bit;

	NEW(SET_BIT_SETS,32,UInt4);
	for(i=0,b=31;i<32;b--,i++) { SET_BIT_SETS[i] = 1 << b; }
	NEW(CARDINALITY_SETS,256,unsigned char);
	for(i=0;i<256;i++){ 
		for(b=n=0;b<8;b++){ if(i & (1<<b)) n++; }
		CARDINALITY_SETS[i] = n;
	}
	NEW(BIT_NUMBER_SETS,9,unsigned char);
	for(bit=1,b=7;b>=0;bit=(bit<<1),b--) BIT_NUMBER_SETS[b]=bit;
}

/*************** Create and Copy operations for sets ********************/
set_typ  MakeSet(Int4 N) /* create and return the null set B = { } */
{
	set_typ	B;

	if(N > MAX_SET_SIZE) seterror("Too many segments.");
	if(CARDINALITY_SETS == NULL) initialize_sets();
	MEW(B,1,set_type);
	B->N = N; 
	B->nint = (Int4) ceil((double)N/(double)32.0);
	B->nbyte= 4*B->nint;
#if 0	// quick fix for 64 bit processing
	NEW(B->i,B->nint,unsigned int);
#else	// 
	NEW(B->i,B->nint,UInt4);
#endif
	B->b = (unsigned char*) B->i;
	return (B);
}

void	CopySet(set_typ B1,set_typ B2)
/* copy set B2 into set B1 - set B2 is modified */
{
	Int4 j; 
	if(B1->N != B2->N) { fprintf(stderr,"N1=%d N2=%d\t",B1->N,B2->N); 
		seterror("incompatible sets!"); }
	for(j=0;j<B1->nint;j++) B1->i[j] = B2->i[j];
}

/**************** Modify operations for sets ********************/
void	ClearSet(set_typ B)
/* set set B equal to the empty set */
#if 0	// quick fix...
{ Int4 j; for(j=0;j<B->nint;j++) B->i[j] = (unsigned int) 0; }
#else
{ Int4 j; for(j=0;j<B->nint;j++) B->i[j] = (UInt4) 0; }
#endif

void	FillSet(set_typ B)
/* set set B equal to the set of integers from 0 to N-1 */
{	Int4 j; 
#if 1	// 0xff = 255 = binary 1111 1111.
	for(j=0;j<B->nbyte;j++) B->b[j] = (unsigned char) 0xff; 
#else	// 
	for(j=0;j<B->nbyte;j++) B->b[j] = (unsigned char) 0xff; 
#endif
	for(j=B->N;j<8*B->nbyte;j++) DeleteSet(j,B);
}

void	DeleteSet(Int4 element,set_typ B)
/* delete an element from set B */
{
	Int4	nbit,nint;
	if(element >= B->nbyte*8) seterror("input numbers too large");
#if 1					// 0x1f is hexadecimal 31 = binary 11111.
	nbit = element & 0x1f;		/* element mod 32 */
#else	// 64 bit.			// 0x3f is hexadecimal 63 = binary 111111.
	nbit = element & 0x3f;		/* element mod 64 */
#endif
	nint = element >> 5;		/* element / 32 */
	if(MemberSet(element,B)) B->i[nint] ^= SET_BIT_SETS[nbit];
}

void	AddSet(Int4 element,set_typ B)
/* add an element to set B */
{
	Int4	nbit,nint;
	if(element >= B->N) seterror("input numbers too large");
	nbit = element & 0x1f;		/* element mod 32 */
	nint = element >> 5;		/* element / 32 */
	B->i[nint] |= SET_BIT_SETS[nbit];
}

UInt4	MemberSet(register Int4 element, register set_typ B)
/* return 0 if element is not a member of B; else return non-zero */
{
	register Int4	nbit,nint;

	if(element >= B->nbyte*8) seterror("input numbers too large");
	nbit = element & 0x1f;
	nint = element >> 5;
	return (B->i[nint] & SET_BIT_SETS[nbit]);
}

/**************** Union operations for sets ********************/
void	UnionSet3(set_typ B1, set_typ B2, set_typ UB)
/* modifies B1 to equal B1 U B2; B1 is returned */
{
	Int4	j;
	for(j=0;j<B1->nint;j++) 
		UB->i[j] = B1->i[j] | B2->i[j];
}

set_typ	UnionSet(set_typ B1, set_typ B2)
/* modifies B1 to equal B1 U B2; B1 is returned */
{
	Int4	j;
	if(B1->N != B2->N) {
		fprintf(stderr,"N1=%d, N2=%d\n",B1->N,B2->N);
		seterror("N1!=N2 for Union operation");
	}
	for(j=0;j<B1->nint;j++) B1->i[j] = B1->i[j] | B2->i[j];
	return B1;
}

/*************** Intersect operations **********************/
void    IntersectNotSet(set_typ B1, set_typ B2, set_typ notIB)
/* modifies notIB to equal B1 intersect not B2 */
{
        Int4    j;
        if(B1->N != B2->N || B1->N != notIB->N)
                print_error("N1!=N2 for Intersection operation");
        for(j=0;j<B1->nint;j++) notIB->i[j] = B1->i[j] & ~(B2->i[j]);
}

void	IntersectNotSet(set_typ B1, set_typ B2)
/* modifies B1 to equal B1 intersect not B2 */
{
	Int4	j;
	if(B1->N != B2->N) seterror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) B1->i[j] = B1->i[j] & ~(B2->i[j]);
}

void	IntersectSet1(set_typ B1, set_typ B2,set_typ IB)
/* modifies IB to equal B1 intersect B2 */
{
	Int4	j;
	if(B1->N != B2->N || B1->N != IB->N) 
		seterror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) {
		IB->i[j] = B1->i[j] & B2->i[j];
	}
}

void	IntersectSet3(set_typ B1, set_typ B2)
/* modifies B1 to equal B1 intersect B2 */
{
	Int4	j;
	if(B1->N != B2->N) seterror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) B1->i[j] = B1->i[j] & B2->i[j];
}

/*************** Cardinality operations **********************/
Int4     CardInterSet(register set_typ B1,register set_typ B2)
/* return cardinality of the intersection of B1 and B2 */
{
	register Int4	i,b,byte,n=0;
	if(B1->N != B2->N) {
		fprintf(stderr,"N1=%d; N2=%d\n",B1->N,B2->N);
		seterror("N1!=N2 for CardIntersect operation");
	}
	for(i=0,byte=0; i < B1->nint; i++,byte+=4) {
	   if(B1->i[i] & B2->i[i])
	     for(b=0;b<4;b++)  {
		n+=CARDINALITY_SETS[(B1->b[byte+b])&(B2->b[byte+b])];
	     }
	} return n;
}

Int4     CardInterSetINotJ(register set_typ SetI,register set_typ SetJ)
/* return cardinality of the intersection of SetI and NotSetJ */
{
        register Int4    i,b,byte,n=0;
        if(SetI->N != SetJ->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",SetI->N,SetJ->N);
                seterror("N1!=N2 for CardIntersect operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           if(SetI->i[i] & SetJ->i[i])
             for(b=0;b<4;b++){
                n+=CARDINALITY_SETS[(SetI->b[byte+b]) & ~(SetJ->b[byte+b])];
             }
        } return n;
}

Int4     CardInterSetNotINotJ(register set_typ SetI,register set_typ SetJ)
/* return cardinality of the intersection of NotSetI and NotSetJ */
{
        register Int4    i,b,byte,n=0;
        if(SetI->N != SetJ->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",SetI->N,SetJ->N);
                seterror("N1!=N2 for CardIntersect operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           if(SetI->i[i] & SetJ->i[i])
             for(b=0;b<4;b++){
                n+=CARDINALITY_SETS[~(SetI->b[byte+b]) & ~(SetJ->b[byte+b])];
             }
        } return n;
}

Int4     CardInterSetNotIJ(register set_typ SetI,register set_typ SetJ)
/* return cardinality of the intersection of NotSetI and SetJ */
{
        register Int4    i,b,byte,n=0;
        if(SetI->N != SetJ->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",SetI->N,SetJ->N);
                seterror("N1!=N2 for CardIntersect operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           if(SetI->i[i] & SetJ->i[i])
             for(b=0;b<4;b++){
                n+=CARDINALITY_SETS[~(SetI->b[byte+b]) & (SetJ->b[byte+b])];
             }
        } return n;
}


Int4     CardSet(register set_typ B)
/* return cardinality of set B */
{
	register Int4	i,b,byte,n=0;

	for(i=0,byte=0; i < B->nint; i++,byte+=4) {
	   if(B->i[i])
		for(b=0;b<4;b++) 
			n+=CARDINALITY_SETS[(B->b[byte+b])];
	}
	return n;
}

/******************* List operations **********************/
Int4	*ListSet(set_typ B)
/* return list of elements in set B; list is  terminated by -1 */
{
	Int4	C,i,j=0,*L;
	
	C=CardSet(B);
	NEW(L,C+1,Int4);
	for(i=0;i<B->N;i++) if(MemberSet(i,B)) { L[j]=i; j++; }
	L[j]= -1;
	return L;
}

/******************* Output operations **********************/
set_typ	PutSet(FILE *fptr,set_typ B)
/* print set B to file pointed to by fptr */
{
	Int4	i,j=0,*m;
	
	m = ListSet(B);
	fprintf(fptr," {");
	for(i=0;m[i] != -1; i++) {
		if(m[i+1] != -1) fprintf(fptr,"%5d,",m[i]);
		else fprintf(fptr,"%5d",m[i]);
		if(i%10==9 && m[i+1] != -1) fprintf(fptr,"\n  ");
	}
	if(j>0) fprintf(fptr,"%5d }\n",m[i]);
	else fprintf(fptr," }\n");
	free(m);
	return B;
}

set_typ  NilSet(set_typ B)
/* destroy set B */
{
	if(B != NULL) { free(B->i); free(B); }
	return (set_typ) NULL;
}

void	seterror(const char *s) { fprintf(stderr,"Set Error: %s\n",s); exit(1); }


