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

#if !defined(PATTERN)
#define PATTERN
#include "stdinc.h"
#include "alphabet.h"
#include "dheap.h"
#include "sset.h"
/**************************** pattern ADT ****************************
	A pattern is a sequence of k sets: G = <s1,s2,...,sk>
**********************************************************************/
typedef struct {
	sst_typ		*m;	/* pattern array of sets */
	unsigned char	i;	/* pointer to next position */
	unsigned char	k;	/* length of pattern */
	unsigned char	nlet;	/* number of letters in G */
} pattern_type;
typedef	pattern_type	*ptn_typ;
/******************************** PRIVATE **********************************/
Int4	AlignScorePatterns(ptn_typ G1,ptn_typ G2, Int4 *os, Int4 minid, 
	Int4 max, Int4 wt);
void	pattern_error(char *s);
#define MAX_PATTERN_LENGTH	200
/******************************** PUBLIC **********************************/
/************************* Pattern Operations **************************/
ptn_typ	Pattern(Int4 k,Int4 nlet);
ptn_typ	CopyPattern(ptn_typ E);
Int4     DepthPattern(ptn_typ G);
Int4     LengFormatPattern(ptn_typ G);
void	ShiftLPattern(ptn_typ P);
void	ShiftRPattern(ptn_typ P);
BooLean CombinablePatterns(ptn_typ G1,ptn_typ G2, Int4 *os);
ptn_typ	CombinePatterns(ptn_typ G1,ptn_typ G2, Int4 *os, Int4 minid, Int4 max,
	Int4 wt);
BooLean	SubPattern(ptn_typ G2, ptn_typ G1);
ptn_typ	PutPattern(FILE *fptr,ptn_typ G,a_type A);
Int4     ResPatterns(Int4 i,char *L, ptn_typ G, a_type A);
Int4     GetPattern(char **motif, ptn_typ G, a_type A);
ptn_typ MergePatterns(ptn_typ G1, ptn_typ G2, Int4 os);
ptn_typ	NilPattern(ptn_typ G);
Int4     LengthPattern(ptn_typ G);
BooLean IdenticalPatterns(ptn_typ G1,ptn_typ G2);
ptn_typ MultiVar2Pattern(Int4 n, ptn_typ G);
ptn_typ Var2Pattern(ptn_typ G);
ptn_typ MultiUSet2Pattern(Int4 n, ptn_typ G);
ptn_typ USet2Pattern(ptn_typ G);

/*************************** Defines ******************************/
#define iPattern(G)			((G)->i)
#define kPattern(G)			((G)->k)
#define setPattern(c,i,G)		((G)->m[(i)]=(1 << (c)))
#define AddPattern(c,i,G)		((G)->m[(i)]|=(1 << (c)))
#define RmPattern(c,i,G)		(((G)->m[(i)]&(1<<(c)))?\
						((G)->m[(i)]^=(1<<(c))):NULL)
#define setVPattern(i,G)		((G)->m[(i)]=0)
#define EqPattern(i,n,G)		((G)->m[(i)]==(n))
#define VarPattern(i,G)		((G)->m[(i)]==0)
#define NotVarPattern(i,G)		((G)->m[(i)]!=0)
#define UPattern(G)			UnivSset((G)->nlet)
#define EqUPattern(i,G)		((G)->m[(i)]==(UPattern(G)))
#define SetUPattern(i,G)		((G)->m[(i)]=(UPattern(G)))
#define NotEqUPattern(i,G)		((G)->m[(i)]!=(UPattern(G)))
#define NotEmptyPattern(i,G)		((G)->m[(i)]!=0)
#define MemPattern(i,c,G)		((G)->m[(i)] & (1 << c))
#endif
