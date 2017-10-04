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

#if !defined (IARRAY)
#define IARRAY
#include "stdinc.h"
/******************************* Integer Array ***************************
Header file for data structure representing array of non-negative integers.

	int *x;   - simple operations on integer lists.

	

 *************************************************************************/
/****************************** private ***************************/
void	iarray_error(char *s);
/****************************** PUBLIC ***************************/
Int4    MaxIntIArray(Int4 *I, Int4 *L);
Int4     OffsetIArray(Int4 os,Int4 *L1,Int4 *L2);
Int4     OffsetIArray2(Int4 os,Int4 *L);
void	BubbleSortIArray(Int4 *L);
Int4     RmCommonIArray(Int4 *LC,Int4 *L1);
Int4     CommonIArray(Int4 *LC,Int4 *L1,Int4 *L2);
Int4     IntersectIArrays(Int4 *LC,Int4 *L1,Int4 *L2);
Int4     UnionIArrays(Int4 *LC,Int4 *L1,Int4 *L2);
Int4     CopyIArray(Int4 *LC,Int4 *L);
void    PutIArray(FILE *fptr, Int4 *L);

/********************** MACROS **********************************/

#endif
