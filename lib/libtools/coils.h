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

/* coils.h - codes and constants for ordered scan program. */
#if !defined (COILS)
#define COILS
#include "stdinc.h"
#include "afnio.h"
#include "probability.h"
#include "histogram.h"
#include "sequence.h"
#include "smatrix.h"
#include "wmodel.h"
#include "mheap.h"
#include "sites.h"
#include "pseg.h"
#include "seqset.h"
#include "guide.h"
#include "residues.h"
#include "fsm.h"

#include <gapxdrop.h>
#include <txalign.h>
#include <tofasta.h>
#include <simutil.h>
#include <sequtil.h>

#include <objres.h>
#include <objseq.h>
#include <seqport.h>

typedef struct {
	a_type  A;
	char	*snfile;
	BooLean	shuffle;
	Int4	N;		/** number of motif models **/
	double	total;
	double	*freq;
	wm_type  *M;
	char	method;		/** method used for model: g,c,r **/
	FILE	*dfp;		/** new database file pointer **/
	FILE	*recomb;	/** file pointer for recombinants **/
} coils_type;
typedef coils_type *csn_typ;

/********************************* private ********************************/
#define MAXSCN_BLOCK_SIZE  5000
#define MAX_BLOCK_LENGTH 100
#define MAX_NUM_MODELS	50
#define	ID_LENGTH	200

void	ReadCoils(csn_typ F, double pseudo,char *snfile, Int4 *counts);
/********************************* PUBLIC ********************************/
csn_typ MkCoils(double pseudo, char *snfile, Int4 *counts, a_type A, char method);
void    NilCoils(csn_typ F);
Int4    MotifInfoCoils(char *snfile, Int4 *counts, a_type A);
SeqEntryPtr     MySeqToSEP(unsigned char *seq, Int4 len_seq, char *info,
                Int4 id, a_type A);
SeqEntryPtr     MySeqToSEP(unsigned char *seq, Int4 len_seq, char *info,
                Int4 id, a_type A);
/********************************* MACROS ********************************/
#define ShuffleSegsCoils(F)	((F)->shuffle = TRUE)
#define SetMethodCoils(x,F)	((F)->method = (char)(x))
#define NumModelsCoils(F)    	((F)->N)
#define SMatrixCoils(m,F)	(((m) <= (F)->N && (m) > 0)?\
					GetSMatrixWModel((F)->M[(m)]): NULL)

#endif

