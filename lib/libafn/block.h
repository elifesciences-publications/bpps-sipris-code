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

#if !defined(BLOCKS)
#define BLOCKS
#include "stdinc.h"

/*************************** ADT Block *******************************

    B = subset of { 0,...,n-1 } where n is an element of Postive integers.

**********************************************************************/
#if 1	// quick 64 bit fix
typedef	unsigned int	block_list_type;
#else
typedef	UInt4	block_list_type;
#endif
// typedef	unsigned short block_list_type;

typedef struct {
	block_list_type	*list;		/* list of non-zero words */
	block_list_type *list2;		/* 2nd list of non-zero words */
	unsigned char	*b;	/* as bytes = 8 bits */
#if 1	// quick 64 bit fix
	unsigned int	*i;		/* as Int4 ints = 32 bits */
#else
	UInt4	*i;		/* as Int4 ints = 32 bits */
#endif
	Int4		N;		/* maximum number in set */
	Int4		nbyte;		/* number of bytes in array */
	Int4		nint;		/* number of ints in array */
} block_type;

typedef	block_type	*b_type;

/*************************** Private *******************************/
extern unsigned char *CARDINALITY_BLOCKS;
#if 1	// quick 64 bit fix
extern unsigned int *SET_BIT_BLOCKS;
#else
extern UInt4 *SET_BIT_BLOCKS;
#endif
extern block_list_type *B1_list,*B2_list,*IB_list;
extern unsigned char *B1_byte,*B2_byte,*IB_byte;

void    initialize_blocks(void);
void    add_list_block(b_type B);
void    update_list_block(b_type B);
void    blockerror(char *s);
/*************************** Public ********************************/
b_type	Block(Int4 N);			/* create Null block */
b_type	BlockL(Int4 N);			/* create Null block with list */
void    CopyBlock(b_type B1,b_type B2);
void    CopyBlockL(b_type B1,b_type B2);
void    CopyBlockL2(b_type B1, b_type B2);
void    CopyBlockLL(register b_type B1,register b_type B2);
void    ClearBlock(b_type B);			/* zero out a set in B */
void    FillBlock(b_type B);
void	AddBlock(Int4 element,b_type B); 	/* add an element to set B */
void	DeleteBlock(Int4 element,b_type B);
UInt4    MemberBlock(register Int4 element, register b_type B);
b_type  UnionBlock(b_type B1, b_type B2);
void    UnionBlock3(b_type B1, b_type B2, b_type UB);
Int4	UnionBlockCL(b_type B1, b_type B2);
void	UnionBlockL(b_type B1, b_type B2);
void    UnionBlockLF(b_type B1, b_type B2, b_type UB);
void    IntersectNotBlock(b_type B1, b_type B2);
void    IntersectBlock3(b_type B1, b_type B2);
void	IntersectBlock1(b_type B1, b_type B2,b_type IB);
Int4     IntersectBlockCF(b_type B1,b_type B2,b_type IB);
void	IntersectBlockF(b_type B1,b_type B2,b_type IB);
void	IntersectBlockL(b_type B1,b_type B2);
Int4     IntersectBlockList(Int4  *L, b_type B1, b_type B2);
Int4     IntersectBlockLCXOR(b_type eB,b_type B, b_type IB);
Int4	CardBlock(b_type B);	/* return cardinality of a set in B */
Int4     CardBlockL(b_type B);
Int4	CardInterBlock(b_type B1,b_type B2);
Int4	CardInterBlockLF(b_type B1,b_type B2);
void	CardInterBlockLMult(register Int4 *C, b_type *B1,b_type B2);
Int4     *ListBlock(b_type B);	/* return list of members in B */
Int4     *ListBlockL(b_type B);	/* return list of members in B */
b_type	PutBlock(FILE *fptr,b_type B);	/* print block B */
b_type	NilBlock(b_type B);		/* destroy ISets */
/*************************** Macros **********************/
// #define END_BLOCK_LIST	65530		// for unsigned short
// #define MAX_BLOCK_SIZE	524000		/* = list of 65500 */
#define END_BLOCK_LIST	1000000
#define MAX_BLOCK_SIZE	((32*END_BLOCK_LIST) - 100) /* = list of 65500 */
#define BlockN(B)	((B)->N)	/* return max cardinality of sets */
#define EmptyBlockL(B)	((B)->list[0]==END_BLOCK_LIST)
#define ClearBlockL(B)	((B)->list[0]=END_BLOCK_LIST)

#endif
