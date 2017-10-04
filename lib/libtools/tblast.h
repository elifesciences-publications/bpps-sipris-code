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

/* tblast.h - transitive blast program. */
#if !defined (TBLAST)
#define TBLAST
#include "gblast.h"
#include "swaln.h"
#include "histogram.h"
#include "purge.h"
#include "seqset.h"
#include "karlin.h"
#include "residues.h"
#include "dheap.h"
#include "mheap.h"
#include "sqheap.h"

/********************* Transitive Blast Type ***********************

/*********************************************************************/

/*************************** generic gblast type **************************/
typedef struct {
	a_type	A;		/** alphabet **/
	Int4	maxdepth;	/** search depth **/
	/************ database structures ****************/
	char	*name;		/** name of query file **/
	char	*dbs;		/** name of dbs to search **/
	double	*dbsfreq;	/** residue frequecy of database **/
	Int4	*counts;	/** numbers of residues in database **/
	Int4	total;		/** total residues in database **/
	Int4	number;		/** number of sequences in database **/
	double	avelen;		/** average sequence length **/
	char	*L;		/** list of dbs hits **/
				/** 1 = forget it; 0 = not hit; **/
				/** -1 = hit but try again. **/
	/************* sequence min/max heap *******************/
	Int4	maxhit;		/** maximum number of hits **/
	sh_typ	SH;		/** sequence heap for matches **/
	/************ statistics ************************/
	Int4	T;
	double	K,lambda,H;	/** Karlin-Altschul **/
	double	expect;		/** blast parameters **/
	double	qval;		/** max seq/seq p-value to querry **/
	h_type	HG;		/** histogram for testing **/
} tblast_type;
typedef tblast_type *tb_typ;
/*********************************************************************/
/* CONSTANTS */
#define MAX_HITS_TB     2000
#define MAX_IN_SEQS     500000
double expm1(double x);

/******************************* private *******************************/
Int4    run_tblast(e_type Q,e_type fullQ,double expect,BooLean lowH,
	Int4 depth,tb_typ T);

/******************************* PUBLIC *******************************/
tb_typ  MakeTBlast(char *name, char *dbs, Int4 maxhit, double q_value, 
	a_type A);
void    NilTBlast(tb_typ B);
Int4    RunTBlast(e_type Q,double expect,BooLean lowH,Int4 depth,tb_typ T);

/******************* EXPERIMENTAL *******************/
Int4    RunMultiTBlast(e_type *Q, double expect, BooLean lowH, Int4 depth, tb_typ T);


/*********************************************************************/
#define HistTBlast(B)	((B)->HG)

/* CODES */


#endif

