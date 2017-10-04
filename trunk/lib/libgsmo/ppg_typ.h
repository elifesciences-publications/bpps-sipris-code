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

/**************************** ppg_typ.h - *******************************
   Gibbs Sampling propagation algorithm for local multiple alignment.
*************************************************************************/
#if !defined(_PPG_TYP_)
#define _PPG_TYP_
#include "gibbs.h"
#include "cma_gmb.h"

class ppg_typ {
public:
	ppg_typ(gs_type gs){ G=gs; }
	~ppg_typ( ){ }
	BooLean Update(char c);
	Int4	Propagate(Int4 sq, double *L, BooLean *moved);
	Int4    PropagateClusters(set_typ Set, Int4 size, double *L);
	Int4    PropagateRepSets(Int4 percent_ident,double MaxFrctn, double &LLR);
	Int4    PropagateRandom(set_typ Set, Int4 size, double *L);
	Int4	MultiPropagate(set_typ Set, double *L,dh_type dH=0);
private:
	gs_type G;
	void    update_propagate_gibbs();
	double	CondProbPropagate(Int4 t, Int4 n);
	Int4    BestSitePropagate(register Int4 end, Int4 t, Int4 n);
	Int4    ChooseSitePropagate(register Int4 s, Int4 t, Int4 n);
	double	LogLikePropagate();
};

#endif

