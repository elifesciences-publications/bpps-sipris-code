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

/**************************** gibbs.h - *******************************
   Gibbs Sampling algorithms for local multiple alignment.
*************************************************************************/
#if !defined(GIBBS)
#define GIBBS
#include "cmsa.h"
#include "dheap.h"
#include "probability.h"
#include "mem_typ.h"

/********************************* PRIVATE ********************************/
typedef struct {
	cma_typ		cmsa; 			/** colinear msa **/
	BooLean		test[10];		/** test items -T option **/
	Int4		*pos;
	/****** sampling parameters ********/
        Int4		nruns,nconverge;
	BooLean		verbose,hot_gibbs;
	/******** propagation items *********/
	char		stage;		//sampling stage = { 0, 1, 2 }
	Int4		limit;
	Int4		maxlen;		// maximum length used for matrix.
	double		**matrix;	// alignment prob matrix[m][i]
	Int4		**strt;		// start[model][seq] 
	Int4		**end;		// end[model][seq]
	Int4		*oldsite;	// oldsite[model]
	/******** simulated annealing *********/
	char		mode;		// sampling mode.
	BooLean		mod_temp;	// modify the sampling temperature
	double		temp;		// inverse of sampling temperature
	double		temp0;		// sampling temperature in 'Kelvins'
	/******** gapped version *********/
	BooLean		gapped_gibbs;
	double		stage1_minmap;
} gibbs_sampler_type;

typedef gibbs_sampler_type *gs_type;

/******************************** private ********************************/
void    update_propagate_gibbs(gs_type G);
Int4	propagate_gibbs(Int4 n, double *L, BooLean *moved, gs_type G);
double	CondProbPropagate(Int4 t, Int4 n, gs_type G);
Int4    BestSitePropagate(register Int4 end, Int4 t, Int4 n, gs_type G);
Int4    ChooseSitePropagate(register Int4 s, Int4 t, Int4 n, gs_type G);
gs_type blank_gibbs( );

BooLean OptionsGibbs(Int4 argc, char *argv[], gs_type G);
double  GibbsSampler(BooLean (*Update)(char, gs_type),
        Int4 (*Sample)(Int4, double *, gs_type), gs_type G);
double  TotLikeGibbs(gs_type G);
double  LogLikePropagate(gs_type G);
BooLean EnlargeGibbs(gs_type G);

/********************************* PUBLIC ********************************/
BooLean SetTemperatureGibbs(double T, gs_type G);
double  RunGibbs(char *options, cma_typ *M);
double  RunGibbs2(Int4 nopt, char *options[], cma_typ *M);
double  RunMinMapGibbs(char *options,double minmap, cma_typ *M);
double  HotGibbs(char *options,cma_typ *M, double temp);
double  SimulatedAnnealingGibbs(Int4 nopt, char *options[], cma_typ *M, 
	char mode, double temp);
double  CoreGibbs(Int4 nopt, char *options[], cma_typ *M, 
	char mode, double temp, double minmap);
double  SimAnnealGibbs(char *options, cma_typ *M, char mode, double temp);
void    InitMAPGibbs(gs_type G);
cma_typ NilGibbs(gs_type G);
cma_typ GibbsCMSA(gs_type G);
BooLean UpdatePropagate(char c, gs_type G);
double	RunPropagate(FILE *fptr, gs_type G);

/********************************* MACROS ********************************/

#endif

