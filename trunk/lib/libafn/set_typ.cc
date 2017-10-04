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
#if 1	// quick 64 bit fix:
unsigned int	*SET_BIT_SETS=NULL;
#else
UInt4 *SET_BIT_SETS=NULL;
#endif

unsigned char *CARDINALITY_SETS=NULL;
unsigned char *BIT_NUMBER_SETS=NULL;

void	initialize_sets(void)
/* initialize lookup tables used for bitwise operations */
{
	Int4 i,b,n;
	unsigned char bit;

#if 1	// quick 64 bit fix:
	NEW(SET_BIT_SETS,32,unsigned int);
#else
	NEW(SET_BIT_SETS,32,UInt4);
#endif
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
	NEW(B,1,set_type);
	B->N = N; 
	B->nint = (Int4) ceil((double)N/(double)32.0);
	B->nbyte= 4*B->nint;
#if 1	// quick fix for 64 bit processing
	NEW(B->i,B->nint+2,unsigned int);
#else	// 
	NEW(B->i,B->nint+2,UInt4);
#endif
	B->b = (unsigned char*) B->i;
	return (B);
}

set_typ	CopySet(set_typ B1)
{
	set_typ B2=MakeSet(B1->N);
	CopySet(B2,B1); // sets B2 to be the same as B1 (see below).
	return B2;
}

void	CopySet(set_typ B1,set_typ B2)
/* copy set B2 into set B1 - set B1 is modified */
{
	Int4 j; 
	if(B1->N != B2->N) { fprintf(stderr,"N1=%d N2=%d\t",B1->N,B2->N); 
		seterror("incompatible sets!"); }
	for(j=0;j<B1->nint;j++) B1->i[j] = B2->i[j];
}

/**************** Modify operations for sets ********************/
void	ClearSet(set_typ B)
/* set set B equal to the empty set */
#if 1	// quick 64 bit fix 
{ Int4 j; for(j=0;j<B->nint;j++) B->i[j] = (unsigned int) 0; }
#else
{ Int4 j; for(j=0;j<B->nint;j++) B->i[j] = (UInt4) 0; }
#endif

void	FillSet(set_typ B)
/* set set B equal to the set of integers from 0 to N-1 */
{	Int4 j; 
	for(j=0;j<B->nbyte;j++) B->b[j] = (unsigned char) 0xff; 
	for(j=B->N;j<8*B->nbyte;j++) DeleteSet(j,B);
}

void	DeleteSet(Int4 element,set_typ B)
/* delete an element from set B */
{
	Int4	nbit,nint;
	if(element >= B->nbyte*8){
	   fprintf(stderr,"DeleteSet( ): element = %d; B->nbyte = %d\n",element,B->nbyte);
	   fprintf(stderr,"DeleteSet( ): B->N = %d; B->nint = %d\n",B->N,B->nint);
	   assert(element < B->nbyte*8);
	   seterror(" DeleteSet() input numbers too large");
	}
#if 0	// for 64 bits
	nbit = element & 0x3f;	/* element mod 64 --> 64 bits in a Int4 */
	nint = element >> 6;	/* element / 64 --> number of longs to jump over */
#else
	nbit = element & 0x1f;	/* element mod 32 --> 32 bits in an int */
	nint = element >> 5;	/* element / 32 --> number of ints to jump over */
#endif
	if(MemberSet(element,B)) B->i[nint] ^= SET_BIT_SETS[nbit];
}

void	AddSet(Int4 element,set_typ B)
/* add an element to set B */
{
	Int4	nbit,nint;
	if(element >= B->N){
	   fprintf(stderr,"AddSet( ): element = %d; B->N = %d\n",element,B->N);
// assert(element <  B->N);
	   assert(element < B->N);
	   seterror("AddSet( ) input numbers too large");
	}
	nbit = element & 0x1f;		/* element mod 32 */
	nint = element >> 5;		/* element / 32 */
	B->i[nint] |= SET_BIT_SETS[nbit];
}

UInt4	MemberSet(register Int4 element, register set_typ B)
/* return 0 if element is not a member of B; else return non-zero */
{
	register Int4	nbit,nint;

	if(element >= B->nbyte*8){
	   fprintf(stderr,"MemberSet( ): element = %d; B->nbyte = %d\n",element,B->nbyte);
	   fprintf(stderr,"MemberSet( ): B->N = %d; B->nint = %d\n",B->N,B->nint);
	   // assert(element < B->nbyte*8);
	   seterror("input numbers too large");
	}
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

Int4     CardUnionSet(register set_typ B1,register set_typ B2)
/* return cardinality of the Union of B1 and B2 */
{
        register Int4   i,b,byte,n=0;
        if(B1->N != B2->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",B1->N,B2->N);
                seterror("N1!=N2 for CardUnionSet() operation");
        }
        for(i=0,byte=0; i < B1->nint; i++,byte+=4) {
           if(B1->i[i] | B2->i[i])
             for(b=0;b<4;b++)  {
                n+=CARDINALITY_SETS[(B1->b[byte+b]) | (B2->b[byte+b])];
             }
        } return n;
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
		seterror("N1!=N2 for CardInterSet() operation");
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
                seterror("N1!=N2 for CardInterSetINotJ() operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           // if(SetI->i[i] & (unsigned char)~SetJ->i[i])	//  very rare that complement of SetJ is 0.
           if(SetI->i[i])
             for(b=0;b<4;b++){
                n+=CARDINALITY_SETS[(SetI->b[byte+b] & (unsigned char)~SetJ->b[byte+b])];
             }
        } return n;
}

Int4     CardInterSetNotINotJ(register set_typ SetI,register set_typ SetJ)
/* return cardinality of the intersection of NotSetI and NotSetJ */
{
        register Int4    i,b,byte,n=0;
        if(SetI->N != SetJ->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",SetI->N,SetJ->N);
                seterror("N1!=N2 for CardInterSetNotINotJ() operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           // if(~SetI->i[i] & ~SetJ->i[i])
             for(b=0;b<4;b++){
#if 1		// undo automatic conversion to int type.
		n+=CARDINALITY_SETS[((unsigned char)~SetI->b[byte+b]) & ((unsigned char)~SetJ->b[byte+b])];
#else
		unsigned char b0 = (unsigned char)(~(SetI->b[byte+b])) & (~(SetJ->b[byte+b]));
		n+=CARDINALITY_SETS[b0];
                // n+=CARDINALITY_SETS[(~(SetI->b[byte+b])) & (~(SetJ->b[byte+b]))];
#endif
             }
        } return n - (SetI->nint*32 - SetI->N);	// // adjusts for highend overflow of complementary set.
	// if not using zero then subtract one from calling environment...
}

Int4     CardInterSetNotIJ(register set_typ SetI,register set_typ SetJ)
/* return cardinality of the intersection of NotSetI and SetJ */
{
        register Int4    i,b,byte,n=0;
        if(SetI->N != SetJ->N) {
                fprintf(stderr,"N1=%d; N2=%d\n",SetI->N,SetJ->N);
                seterror("N1!=N2 for CardInterSetNotIJ() operation");
        }
        for(i=0,byte=0; i < SetI->nint; i++,byte+=4) {
           // if(~SetI->i[i] & SetJ->i[i])
           if(SetJ->i[i])
             for(b=0;b<4;b++){
                n+=CARDINALITY_SETS[((unsigned char)~SetI->b[byte+b]) & (SetJ->b[byte+b])];
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

set_typ	ReadSet(FILE *fp)
{
	set_typ set=0;
	Int4	rtn;
	UInt4	i;

	if(CARDINALITY_SETS == NULL) initialize_sets();
#if 1	// read byte indicating null or non-null set within array.
	char	okay;
	rtn=fread(&okay,sizeof(char),1,fp);
	if(rtn != 1) seterror("ReadSet() input error 0");
	if(okay==0) return 0;
#endif
	NEW(set,1,set_type);
	rtn=fread(set,sizeof(set_type),1,fp);
	if(rtn != 1 || set->N >= MAX_SET_SIZE) seterror("ReadSet() input error 1");
	i=(UInt4) ceil((double)set->N/(double)32.0);
	if(i != set->nint){
		fprintf(stderr,"i=%d; set->nint=%d; set->N=%d\n",i,set->nint,set->N);
		seterror("ReadSet() input error 2");
	}
	NEW(set->i,set->nint+3,unsigned int);
	rtn=fread(set->i,sizeof(unsigned int),set->nint,fp);
	if(rtn != set->nint) seterror("ReadSet() input error 3");
	set->b = (unsigned char*) set->i; 
	return set;
}

set_typ	*ReadSets(FILE *fp,Int4 &Num)
{
	Int4 rtn,i,Number;
	set_typ *set=0;

	rtn=fread(&Number,sizeof(Int4),1,fp);
	if(rtn != 1) return 0;	// signal to calling enviroment if multiple files done being read.
	// fprintf(stderr,"Number = %d\n",Number);
	// if(rtn != 1) seterror("ReadSets() input error 1");
	Num=Number; NEW(set,Number+3,set_typ);
	for(i=1; i <= Number; i++){
	   set[i]=ReadSet(fp); 
	   // fprintf(stderr,"set %d/%d = %d\n",i,Number,CardSet(set[i]));
	} return set;
}

void	WriteSet(FILE *fp,set_typ set)
{
	Int4	rtn;

#if 1	// write byte indicating null or non-null set within array.
	char	okay;
	if(set==0){
	  okay=0; rtn=fwrite(&okay,sizeof(char),1,fp);
	  if(rtn != 1) seterror("WriteSet() input error 0");
	  return;
	} else {
	  okay=1; rtn=fwrite(&okay,sizeof(char),1,fp);
	  if(rtn != 1) seterror("WriteSet() input error 1");
	}
#endif

	rtn=fwrite(set,sizeof(set_type),1,fp);
	if(rtn != 1) seterror("WriteSet() input error 2");
	rtn=fwrite(set->i,sizeof(unsigned int),set->nint,fp);
	if(rtn != set->nint) seterror("WriteSet() input error 3");
}

void	WriteSets(FILE *fp,Int4 Number, set_typ *set)
{
	Int4	i,rtn;
	assert(Number > 0);
	rtn=fwrite(&Number,sizeof(Int4),1,fp);
	if(rtn != 1) seterror("WriteSets() input error 1");
	for(i=1; i <= Number; i++) WriteSet(fp,set[i]);
}

set_typ	ParseSet(char *str, Int4 &Low, Int4 &High)
// input string: "3,5-7,9,11-17"
// returns a set of size N=17+3 (0...N-1) // size larger than highest value..
{
        Int4    low,high,x,n,i,j,v,w;
	char	*str0=str;
	set_typ	RtnSet=0;
	const char msg[]="ParseSet() syntax error";

	// find the range of input and check syntax...
	if(!isdigit(str[0])) print_error(msg);
        for(x=0,low=-1,high=0,n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%d", &v) != 1){ print_error(msg); }
                else {
		   if(low==-1){ low = v; x = v; } else if(v <= x) print_error(msg);
		   else high=x=v; 
                   while(isdigit(str[0])) str++;
                   if(str[0] == '-'){
                        str++; if(!isdigit(str[0])) print_error(msg);
                        if(sscanf(str,"%d", &w) != 1) print_error(msg);
                        if(w <= v) print_error(msg);
                        for(i=v+1; i <= w; i++){ n++; }
			high=x=w;
                        while(isdigit(str[0])) str++;
                   }
                }
           } else print_error(msg);
        } str=str0; Low=low; High=high; 
	// fprintf(stderr,"Low=%d; High=%d\n",low,high);
	if(high<=0) print_error(msg);
	// Create a set 
	RtnSet=MakeSet(high+3);
	if(!isdigit(str[0])) print_error(msg);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%d", &v) != 1){ print_error(msg); }
                else {
		   AddSet(v,RtnSet);
                   while(isdigit(str[0])) str++;
                   if(str[0] == '-'){
                        str++; if(!isdigit(str[0])) print_error(msg);
                        if(sscanf(str,"%d", &w) != 1) print_error(msg);
                        if(w <= v) print_error(msg);
                        for(i=v+1; i <= w; i++){ n++; AddSet(i,RtnSet); }
                        while(isdigit(str[0])) str++;
                   }
                }
           } else print_error(msg);
        } return RtnSet;
}

char	*RtnStrSet(set_typ Set, Int4 &Low, Int4 &High)
{
	Int4	i,s,e,S,E,n=CardSet(Set);
	char	*str0,*str; NEW(str0,n+3,char); 
	for(S=-1,E=0,i=1; i < SetN(Set); i++){ if(MemberSet(i,Set)){ if(S < 0) S=i; else E=i; } }
	Low=S; High=E;
	for(str=str0,i=1; i < SetN(Set); i++){
	     if(MemberSet(i,Set)){
		sprintf(str,"%d",i); while(str[0] != 0) str++;
		for(s=e=i; MemberSet(e+1,Set); e++) ;
		if(e > s){ i=e; sprintf(str,"-%d",e); while(str[0] != 0) str++; }
		if(e < E){ sprintf(str,","); while(str[0] != 0) str++; }
	     }
	} return str0;
}

