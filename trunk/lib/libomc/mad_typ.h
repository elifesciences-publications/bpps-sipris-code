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

#if !defined (_MAD_TYP_)
#define _MAD_TYP_
#include "sequence.h"
#include "set_typ.h"
#include "cmsa.h"
#include "residues.h"
#include "probability.h"
#include "table.h"
#include "sset.h"
#include "wdigraph.h"
#include "cth_typ.h"
#include "dheap.h"
#include "rst_typ.h"
#include "ctn_typ.h"
#include "swt_typ.h"
#include "swaln.h"


double  ChiSqTable(double *prob,t_type T);
// double  CnTable(t_type T,double *chsq, double *df);

class mad_typ { // multiple alignment dependence type...
public:
                mad_typ( ){ assert(!"Illegal constructor"); }
		mad_typ(swt_typ *inswt, set_typ inset, cma_typ CMA,double MinFreq,
				double max_gap_frq,char mode)
		{ init(inswt,inset,CMA,MinFreq,max_gap_frq,mode); }
                ~mad_typ( ){ Free(); }
	double	SearchSpace( );
	ctn_typ	*CliqueClusters(Int4,Int4*,double,double,sst_typ **&);
	char	**CliqueStrings(Int4,Int4*,double,double,sst_typ **&);
	sst_typ	**SortedResSets( );
	void	Verbose( ){  verbose=TRUE; }
	e_type	BstSeq(){ return bstE; }
	e_type	CnsSeq( ){ return csqE; }
private:
	// begin new stuff.
	BooLean	*skip;
	void    SubGrpCsqBsq(set_typ inset);
	e_type	bstE,csqE;
	swt_typ	*swt;
	// end new stuff.
        BooLean	TableItem(unsigned short i,unsigned short j,double Fisher_pcut);
	Int4	MaxNumEdges( );
	Int4	TotalNodes( );
	BooLean	verbose;
        void    Free();
        void    init(swt_typ*,set_typ, cma_typ , double,double,char);
	void	init_res_sets(){ init_res_sets('S'); }
	void	init_res_sets(char mode );
	void	Init_Res_Sets(const char *ResidueSets[21][18],const char NumberResSets_L[21]);
	Int4	max_residue_set;
	ctn_typ	*ctn;
	char    **sst_str;
	Int4	NumSsetStr;
	Int4	NumClust;
	Int4	NumCnTab;
	cti_typ	**CnTabList;
	Int4	*cell[3];
	Int4	expected( );
	double	E[3][3];
	// store CardSets for these...
        set_typ **Set;	// Set[i][rsets]
	set_typ	SetI,SetJ,USet,NotSetI,NotSetJ;
        Int4    Length;
        Int4    NumSeqs;
	e_type	keyE;
	Int4	**key_obs;
	double	**keyFreq;	// keyFreq[i][r];
	double	MinKeyFrq;	// minimum frequency of residues in cma...
	double	MaxGapFrq;	// maximum number of gap residues (x) allowed.
	Int4	*NumGaps;
	cma_typ	cma;
	a_type	AB;
	t_type	**table[21][21];	// Tables for computing...
	t_type	Tab;
	sst_typ **ResidueSSet;
	unsigned char	*NumResSets;
	cth_typ *cth;		// contingency table heap.
	h_type	HG;
	rst_typ	*rst;
};

#endif

