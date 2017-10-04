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

/* purge.h - codes and constants for purge program. */
#if !defined (PURGE)
#define PURGE
#include "afnio.h"
#include "block.h"
#include "residues.h"
#include "scaninfo.h"
#include "seqset.h"
#include "pairaln.h"
#include "dheap.h"
#include "gblast.h"
/*** FOR PURGE SCAN INFO... ***/
#include "scanheap.h"
#include "swaln.h"
#include "karlin.h"
#include "residues.h"
#include "dsets.h"
#include "histogram.h"


BooLean	*RmHomologsList(Int4 cutoff, char method, Int4 minimum, Int4 maximum,
        BooLean query, Int4 T, Int4 N, e_type *List, a_type A);
Int4	RmHomologs(Int4 cutoff, char method, Int4 minimum, Int4 maximum,
	BooLean query, ss_type P, char *outfile, Int4 T);
e_type	*PurgeScanInfoRpts(Int4 cutoff, Int4 Nt_flank, Int4 Ct_flank,
	Int4 nblk, Int4 *blklen, snh_typ sH, a_type A);
e_type *PurgeScanInfo(Int4 cutoff, Int4 Nt_flank, Int4 Ct_flank, snh_typ sH, 
	a_type A);
Int4    PurgeFiles(char *infile, char *outfile, Int4 cutoff, Int4 inc,
        Int4 minseq, Int4 maxseq, a_type A);
Int4	RmIdent(e_type *List);

#if 0
        /********* N   A   C   T   G *******/
#define DNA_MTRX "-4  -4  -4  -4  -4 \
                  -4   5  -4  -4  -4 \
                  -4  -4   5  -4  -4 \
                  -4  -4  -4   5  -4 \
                  -4  -4  -4  -4   5 "
#endif

#endif

