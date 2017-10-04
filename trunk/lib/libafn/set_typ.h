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

#if !defined(SET)
#define SET
#include "stdinc.h"

/*************************** ADT Set *******************************

    B = subset of { 0,...,n-1 } where n is an element of Postive integers.

**********************************************************************/
typedef struct {
	unsigned char	*b;	/* as bytes = 8 bits */
#if 1	// quick fix for 64 bit compilations
	unsigned int	*i;		/* as int = 32 bits */
#else
	UInt4	*i;		/* as Int4 ints = 32 bits */
#endif
	Int4		N;		/* maximum number in set */
	Int4		nbyte;		/* number of bytes in array */
	Int4		nint;		/* number of ints in array */
} set_type;

typedef	set_type	*set_typ;

/*************************** Private *******************************/
extern unsigned char *CARDINALITY_SETS;
#if 1	// quick 64 bit fix
extern unsigned int *SET_BIT_SETS;
#else
extern UInt4 *SET_BIT_SETS;
#endif

void    initialize_sets(void);
void    seterror(const char *s);
/*************************** Public ********************************/
set_typ	MakeSet(Int4 N);			/* create Null set */
void    CopySet(set_typ B1,set_typ B2);
set_typ	CopySet(set_typ B1);
void    ClearSet(set_typ B);			/* zero out a set in B */
void    FillSet(set_typ B);
void	AddSet(Int4 element,set_typ B); 	/* add an element to set B */
void	DeleteSet(Int4 element,set_typ B);
UInt4    MemberSet(register Int4 element, register set_typ B);
set_typ  UnionSet(set_typ B1, set_typ B2);
void    UnionSet3(set_typ B1, set_typ B2, set_typ UB);
Int4     CardUnionSet(register set_typ B1,register set_typ B2);
void    IntersectNotSet(set_typ B1, set_typ B2, set_typ notIB);
void    IntersectNotSet(set_typ B1, set_typ B2);
void    IntersectSet3(set_typ B1, set_typ B2);
void	IntersectSet1(set_typ B1, set_typ B2,set_typ IB);
Int4	CardSet(set_typ B);	/* return cardinality of a set in B */
Int4    CardInterSet(register set_typ B1,register set_typ B2);
Int4    CardInterSetINotJ(register set_typ SetI,register set_typ SetJ);
Int4    CardInterSetNotIJ(register set_typ SetI,register set_typ SetJ);
Int4    CardInterSetNotINotJ(register set_typ SetI,register set_typ SetJ);
Int4    *ListSet(set_typ B);	/* return list of members in B */
set_typ	PutSet(FILE *fptr,set_typ B);	/* print set B */
set_typ	NilSet(set_typ B);		/* destroy ISets */

set_typ ReadSet(FILE *fp);	// binary files
set_typ *ReadSets(FILE *fp,Int4 &Num);
void    WriteSet(FILE *fp,set_typ set);
void    WriteSets(FILE *fp,Int4 Number, set_typ *set);
set_typ ParseSet(char *str, Int4 &Low, Int4 &High);
char    *RtnStrSet(set_typ Set,Int4 &Low, Int4 &High);

/*************************** Macros **********************/
#if 1	// 64 bit
#define MAX_SET_SIZE	50000000 /* arbitrary limit: can be increased to 2^63 */
#else
#define MAX_SET_SIZE	50000000 /* arbitrary limit: can be increased to 2^31 */
#endif

#define SetN(B)	((B)->N)	/* return max cardinality of sets */

#endif
