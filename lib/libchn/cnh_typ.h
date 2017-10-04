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

#if !defined (_CNH_TYP_)
#define _CNH_TYP_

#include "c2h_typ.h"
#include "wdigraph.h"
#include "hpt_typ.h"
#include "dsets.h"
#include "tree2hpt.h"
#include "alphabet.h"
#include "residues.h"
#include "clique.h"

class cnh_typ {	// compare two hierarchies type.
public:
	cnh_typ( ) { assert(!"Illegal constructor"); }
	cnh_typ(Int4 argc, char *argv[]){ Init(argc,argv); }
	wdg_typ	GetConsTree( );
	~cnh_typ( ) { Free(); }

private:
	wdg_typ PurgeTree(Int4 Root, Int4 max_wt, Int4 &NumRm, wdg_typ Tree);
	wdg_typ RtnConsensusTree(Int4 Root, Int4 max_wt,wdg_typ Grph);
	void	PutConSetOverlap(FILE *fp,set_typ tstSet=0);
	set_typ	GetMemberGraph(wdg_typ G);	// return set of nodes with 'in'-arcs.
	void	Free();
	void	Init(Int4 argc, char *argv[]);
	char	**Name;
	wdg_typ **tree,grph,Grph;
	ds_type sets;
	hpt_typ **Hpt;
	Int4    NumTrees,***ident,*NumNodes,TotalNodes,Nnds;
	Int4	*Node2Tree,*Node2id,**Node;
	set_typ **NodeSet,**MiscSet;
	set_typ	*cns_set;
	FILE	*ofp;
	char	str[500],*infile;
	cma_typ	MainCMA;
	a_type	AB;
};

#endif

