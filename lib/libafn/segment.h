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

#if !defined(SEGMENT)
#define SEGMENT
#include "stdinc.h"

/**************************** Segment ADT ****************************
		S == <I,o>	2-tuple

		o == offset from start of entity sequence
		I == sequence entity Identifier (E.I)

   [.............Int4.............] [.............Int4...............]
    <-------- Identifier -------->   <----------- offset ----------->

**********************************************************************/
#if 0
typedef struct {
	short	id;
	unsigned short	offset;
	char	sign;
} segment_type;
#endif

typedef UInt8	s_type;

/******************************** PRIVATE *********************************/
void	seg_error(char *s);
/******************************** PUBLIC *********************************/
/************************* Segment Operations **************************/
s_type	Segment(Int4 I, Int4 offset);
Int4     CopySegments(register s_type *SC, register s_type *S);
Int4     OffsetSegments(register Int4 os,register s_type *S1,
	register s_type *S2);
Int4     UnionSegments(s_type *SC,s_type *S1,s_type *S2);
Int4     IntersectSegments(register s_type *SC, register s_type *S1,
        register s_type *S2);
BooLean MemberSegments(s_type *L, s_type S);
Int4     DiffSegments(s_type *SD,s_type *S1,s_type *S2);
void	PutSegment(s_type S);
Int4	SegmentI(s_type S);
Int4	SegmentStart(s_type S);
Int4	BubbleSortSegments(s_type *S);
Int4     IntersectOSegments(Int4 *C, Int4 r, Int4 *os, Int4 k1, Int4 k2, 
	register s_type *S1, register s_type *S2);
/********************* Segment Macro Operations **************************/
#define NilSegment(S)		((S)=NULL)
#define NilSegmentList(S)	free(S)
#define CpSegment2List(i,L,S)		(L[(i)]=(S))
#endif

