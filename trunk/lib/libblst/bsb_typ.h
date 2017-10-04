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

#ifndef __BSB_TYPE__
#define __BSB_TYPE__
// Taken from blastdef.h by Tom Madden
/* Revision 6.8  1999/03/21 19:40:30  madden */

#include "stdinc.h"
#include "my_ncbi.h"
#include "my_ncbimath.h"
#include "my_blastkar.h"
#include "gpxdrop.h"
#include "hspheap.h"
#include "brh_typ.h"

typedef struct gblast_search_block {
	Boolean		positionBased;	// psiblast??
	Boolean		posConverged;
	Int4		query_length;	// length of sequence. 
	unsigned char	*query_sequence;	// Actual (perhaps transl.) sequence.
	brh_typ		*results;
        Int4    	hitlist_count;	// number of results found.
	Int4		hitlist_max;	// maximum number of results.
	sbp_typ 	sbp;		// info on scoring. 
	double		ethresh; 	// for psiblast.
	Int4		pseudoCountConst;	// for psiblast.
	Int4		maxNumPasses;		// for psiblast.
} bsb_type,*bsb_typ;	// GBlastSearchBlk, *GBlastSearchBlkPtr;

bsb_typ GBlastSearchBlkNew(Int4 results_size);
void	ResetGBlastSearchBlk(bsb_typ bsb);
bsb_typ GBlastSearchBlkDelete(bsb_typ bsb);

#endif /* __BSB_TYPE__ */

