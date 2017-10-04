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

#ifndef __BRH_TYPE__
#define __BRH_TYPE__
// Adapted from blast code by Tom Madden
/* Revision 6.8  1999/03/21 19:40:30  madden */

#include "my_ncbi.h"
#include "my_ncbimath.h"
#include "my_blastkar.h"
#include "gpxdrop.h"
#include "hspheap.h"
#include "stdinc.h"

typedef struct _gblast_result_hitlist {
	hsp_typ	*hsp_array;	// An array holding the HSP's.
	Int4	array_max;	// length of array.
	double	best_evalue;	// best evalue in all the HSP's. 
	Int4	high_score; 	// HSP with highest score. 
	Int4	hspcnt;		// Number of HSP's. 
	unsigned int subject_id;	// ID of the subject.
	e_type	sE;		// AFN: subject sequence pointer
	Int4    subject_length; // length of the database sequence.
	// sap_typ seqalign; 	// alignment, if this a gapped calculation.
} brh_type,*brh_typ; // =  GBLASTResultHitlist, *GBLASTResultHitlistPtr;

brh_typ MakeGBLASTResultHitlist(Int4 size, e_type sE);
brh_typ NilBRH(brh_typ result);
void	AddHspBRH(hsp_typ hsp, brh_typ result);
sap_typ    ExtractAlnBRH(brh_typ result, e_type qE);
double	BestEvalBRH(brh_typ result);

typedef struct {
        mh_type         mH;             // minmax heap for results
        brh_typ		*brh_array;	// storage...
        Int4            size;           // heap size.
} brhheap_type,*rhp_typ;

rhp_typ MakeBRHHeap(Int4 hpsz);
void    NilBRHHeap(rhp_typ H);
brh_typ DelMinBRHHeap(double *key, rhp_typ H);
brh_typ SeeBRHHeap(Int4 item, rhp_typ H);
Int4    InsertBRHHeap(brh_typ result, rhp_typ H);


#endif /* __BRH_TYPE__ */

