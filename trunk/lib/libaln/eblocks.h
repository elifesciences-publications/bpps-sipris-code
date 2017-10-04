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

#if !defined (EBLOCKS)
#define	EBLOCKS
#include "afnio.h"
#include "stdinc.h"
#include "mheap.h"
#include "iarray.h"
#include "dheap.h"
#include "block.h"
#include "seqset.h"
#include "segment.h"
/*********************** ADT Elementary Blocks ************************
	
**********************************************************************/

/*************************** EBLOCKS type **************************/
typedef struct {
	ss_type		P;		/* aligned segment population */
	a_type		A;		/* alphabet */
	Int4		k_max;		/* segment size */
	Int4		w;		/* 2 * segment size (2k) */
	Int4		offset;		/* offset of elementary blocks */
        Int4             nseg;           /* number of segments in population */
        s_type          *segment;       /* array for segments */
	Int4		*N;		/* number of segments of length k */
	b_type		**eb,**b;	/* population elementary blocks */
	b_type		***eb2,***b2;	/* population elementary blocks */
	b_type		Bu;		/* universal block */
        b_type          dummy;          /* dummy block */
	char		*res;		/* ordered residues */
	char		**rel;		/* rel[c][0-n]: related res. w/ fq<c */
	Int4		mode;		/* determines # related residues */
} eblocks_type;
typedef eblocks_type *ebs_typ;

/******************************* PRIVATE *****************************/
/******************************** general *****************************/
char    SegValXEBlocks(Int4 i, s_type S, ebs_typ D);
void	NilEBlocks2(ebs_typ D);
void	eblocks_error(char *s);

/******************************* PUBLIC ******************************/
ebs_typ MkEBlocks(ss_type P, Int4 k_max);
void	AddEBlocks2(ebs_typ D, Int4 mode);
Int4     PutEBlocks(FILE *fptr,ebs_typ Ebs);
void    NilEBlocks(ebs_typ D);
void	ShiftREBLocks(ebs_typ D);
void	ShiftLEBLocks(ebs_typ D);
s_type  *List2SegsEBlocks(Int4 *L, ebs_typ D);
Int4     NumEBlocksPL(Int4 **List, b_type B, ebs_typ D);
Int4     NumEBlocksP(b_type B,ebs_typ D);
#if 0
Int4     NonOverlapEBlocks(Int4 **L, Int4 minseq, Int4 length, b_type B,
	ebs_typ D, Int4 *N);
#endif

/**************************** macro operations **********************/
#define	EBlocksN_k(k,D)		((D)->N[(k)])
#define	NSegsEBlocks(D)		((D)->nseg)
#define	EBlocksN(D)		((D)->N)
#define	EBlocksRes(i,D)		((D)->res[(i)])
#define	EBlocksBu(D)		((D)->Bu)
#define	EndFlagEBlocks(D)	-1
#define	RelatedEBlocks(c,D)	((D)->rel[(c)])
#define UnionEBlock(i,c,UB,D)	\
			UnionBlock((UB),(D)->b[(i)][(c)])
#define IntersectEBlockBu(i,c,B,D)	\
			IntersectBlock1((D)->b[(i)][(c)],(D)->Bu,(B))
#define IntersectEBlockCF(i,c,IB,IIB,D) \
			IntersectBlockCF((D)->b[(i)][(c)],(IB),(IIB))
#define IntersectEBlockLCXOR(i,c,EOB,IIB,D) \
			IntersectBlockLCXOR((D)->b[(i)][(c)],(EOB),(IIB))
#define IntersectEBlockF2(i,c,c2,IB,IIB,D) \
			IntersectBlockF((D)->b2[(i)][(c)][(c2)],(IB),(IIB))



#endif

