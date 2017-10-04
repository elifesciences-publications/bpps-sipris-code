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

/**************************** genetic.h - ********************************
  Genetic algorithm for propagation
 *************************************************************************/
#if !defined(GENETIC)
#define GENETIC
#include "gibbs.h"
#include "mheap.h"
#include "msaheap.h"
#include "binomial.h"
#include "histogram.h"

/************************ Genetic Data Structure  ************************/
typedef struct {
	ss_type		data;	/** input data set **/
	/*************** population *****************/
	mah_typ		maH;	/** population of multiple alignments **/
	Int4		size;	/** population size **/
	/*************** blocks & columns **********************/
	bn_type		Bc,Bb;	/** binomials for columns & blocks **/
	Int4		totcols;	/** total columns setting **/
	Int4		*ncols;		/** temp. store # columns **/
} genetic_type;

typedef genetic_type *ga_typ;

/******************************** private ********************************/

/********************************* PUBLIC ********************************/
ga_typ  MakeGenetic(Int4 size, Int4 aveblk, Int4 avecol, Int4 minlenseq);
void    PutGenetic(FILE *fp, ga_typ G);
cma_typ NilGenetic(ga_typ G);

void    SetGenetic(cma_typ msa, ga_typ GA);
Int4    GeneticSampB(Int4 m, Int4 max, ga_typ GA);
Int4    GeneticSampC(Int4 b, Int4 m, Int4 max, ga_typ GA);

Int4    GeneticSampC2(Int4 t, Int4 m,Int4 max,ga_typ GA);
Int4    BestBlksGenetic(ga_typ GA);
Int4    BestColsGenetic(Int4 t, ga_typ GA);

/********************************* MACROS ********************************/
#define	AveColsGenetic(G)	AveBinomial((G)->Bc)
#define	AveBlksGenetic(G)	AveBinomial((G)->Bb)
#define	HeapGenetic(G)		((G)->maH)

#endif

