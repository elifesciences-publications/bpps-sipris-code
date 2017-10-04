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

#if !defined(_SBA_)
#define _SBA_

/************* Sampling Block Alignment type ***************/
#include "prtn_model.h"
#include "smatrix.h"
#include "sequence.h"
#include "residues.h"
#include "alphabet.h"
#include "histogram.h"
#include "dheap.h"

class sba_typ {
  public:
		sba_typ();
        	sba_typ(Int4 sequence_len, unsigned char *sequence,
                        Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 nblocks, 
				smx_typ *M, Int4 **gapfnct, double pernats);
		~sba_typ( );
        void    initialize(Int4 sequence_len, unsigned char *sequence,
                        Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 nblocks, 
				smx_typ *M, Int4 **gapfnct, double pernats);
	char	*TraceBack(Int4 *,Int4 *,Int4 *,double);
	char    *TraceBack(Int4 *,Int4 *, Int4 *, double ,double);  // tm==0.0 -> EM
	char    *TraceBack(Int4 *, Int4 *, Int4 *, double , double *RE);
	char    *EMTraceBack(Int4 *alignscore, Int4 *trace_length, 
			Int4 *start, double pernats);
	char    *TraceBack(Int4 *,Int4 *,Int4, double);
	char    *TraceBack(Int4 *,Int4 *,Int4, double,double);
	Int4    CountIndels(char *, Int4 , Int4 *, Int4 *);
  private:
	void		init();
	void		alloc_weights( );
	double		calc_indel_weight_sba(Int4 , Int4 );
	double		indel_weight_sba(Int4, Int4, Int4, double);
	double		indel_weight_sba(Int4, Int4, double);
	void	 	Free();
	BooLean 	Init();
	double		**MAT, **INSO, **INSE, **DELO, **DELE;
        double          **rMAT, **rINSO, **rINSE, **rDELO, **rDELE;
	double		*inso, *inse, *delo, *dele, **gpen, **match;
	double		totalLike;
	Int4		*block_lengths, prof_len, nmbr_of_blocks,*start_prof;
	Int4		alph_len,seq_len,num_stored,*maxgap;
	Int4		*store_i;
	char		*storeState;
	double		*storeLike;
	unsigned char	*seq;
	double		**weight;
};

#endif

