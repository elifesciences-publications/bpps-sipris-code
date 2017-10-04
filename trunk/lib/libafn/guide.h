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

#if !defined(GUIDE)
#define GUIDE
#include "stdinc.h"
#include "afnio.h"
#include "sequence.h"

/********************* guide sequence data type **********************/
typedef struct {
	e_type		E;
	Int4		len;
	BooLean		*ignore;
} guide_type;

typedef	guide_type	*gd_type;

/*************************** private *********************************/
void    guide_error(char *s);

/**************************** PUBLIC *********************************/
gd_type	MkGuide(char *infile, a_type A);
gd_type	CopyGuide(gd_type G);
void	PutGuide(FILE *fptr,gd_type G, a_type A);
void    NilGuide(gd_type G);

/**************************** MACROS *********************************/
#define GuideE(G)		((G)->E)
#define GuideIgnore(s,G)	(((s) > 0 && (s) <= (G)->len)?\
				  (G)->ignore[(s)]: FALSE)

#endif

