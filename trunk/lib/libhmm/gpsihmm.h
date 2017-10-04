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

#if !defined (_GPSIHMM_)
#define _GPSIHMM_
#include "afnio.h"
#include "random.h"
#include "residues.h"
#include "alphabet.h"
#include "dheap.h"
#include "mheap.h"
#include "sequence.h"
#include "probability.h"
#include "wdigraph.h"
#include "hmm_mem.h"

/*********************** Finite State Machine ************************
 *********************************************************************/

/*************************** generic gblast type **************************/
typedef struct {
	a_type  A;		/** alphabet **/
	Int4	**matrix;	// position-specific scoring matrix.
	e_type  E;              /** query sequence **/
	Int4	x_dropoff; 	// amount of drop to stop extension.
	/******************************************************/
	Int4	cutoff;		/** HSP score cutoff **/
	Int4	*hit_s;		/** hit start **/
	Int4	*hit_e;		/** hit end **/
	Int4	*hit_d;		/** hit diagonal **/
	Int4	*hit_s2;	/** hit start **/
	Int4	*hit_e2;	/** hit end **/
	Int4	*hit_d2;	/** hit diagonal **/
	BooLean	update;		/** have the hits been updated? **/
	double	*score;		/** key storage **/
	Int4	nhits;		/** number of hits **/
	mh_type	mH;		/** heap for best hits **/
	/******************************************************/
	Int4	zero;		/** diagonal zero value: e.g. -25,000 **/
	Int4	*diag0,*diag;	/** diagonal list [1,4,0,0,2] **/
	Int4	*ed0,*extdiag;	/** limit to which diagonal was extended **/
} gpsihmm_type;
typedef gpsihmm_type *gph_typ;
/*********************************************************************/
#define MAX_SEQ_LENG_GPB        30000
/******************************* private *******************************/
Int4    ExtendGPsiHMMStr(Int4 lenMaster, Int4 i1, Int4 len2, unsigned char *seq2,
        Int4 i2, Int4 *left, Int4 *right, Int4 **matrix,register Int4 x_dropoff, Int4 *m2m);
void	UpdateHitsGPsiHMM(gph_typ B);
/******************************* PUBLIC *******************************/
gph_typ MakeGPsiHMM(Int4 hpsz, e_type E, a_type A, Int4 max_proflen);
void    NilGPsiHMM(gph_typ B);
Int4  MatcherGPsiBlstStr(hit_typ head,Int4 query_length,Int4 uthreshold,
        Int4 ugpxdrop,Int4 **mtx,Int4 subj_len,register unsigned char *subj_ptr,
        register gph_typ B,Int4 *nwordhts);
/******************************** MACROS *****************************/
#define nHitsGPsiHMM(B)		(((B)->update)? UpdateHitsGPsiHMM(B),\
					(B)->nhits: (B)->nhits)
#define StartGPsiHMM(r,B)	(((B)->update)? UpdateHitsGPsiHMM(B),\
				  (B)->hit_s2[(r)]: (B)->hit_s2[(r)])
#define EndGPsiHMM(r,B)		(((B)->update)? UpdateHitsGPsiHMM(B),\
				  (B)->hit_e2[(r)]: (B)->hit_e2[(r)])
#define DiagGPsiHMM(r,B)		(((B)->update)? UpdateHitsGPsiHMM(B),\
				  (B)->hit_d2[(r)]: (B)->hit_d2[(r)])

#endif
