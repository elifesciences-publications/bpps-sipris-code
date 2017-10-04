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

#if !defined (_CTN_TYP_)
#define _CTN_TYP_
#include "set_typ.h"
#include "clique.h"
#include "cth_typ.h"

// ctn_typ... // node item type
class ctn_typ {		// contingency table grf type.
		// [.....|.....|.....|.....]  
public:
		ctn_typ( ) { assert(!"Illegal constructor"); }
		ctn_typ(Int4 len, Int4 ms, Int4 TrueNumNodes){ Init(len,ms,TrueNumNodes); }
		~ctn_typ( ){ Free(); }
	unsigned short	NumNodes(){ return N; }
	unsigned short	Node(unsigned short site,unsigned char set)
   	// Note that: (Node - 1) % N + 1 == Site.
		{
			if(site > Len || set > max_set) return 0;
			// return ((Len*(set-1)) + site);  // Residue set clustering
			return ((max_set*(site-1)) + set);  // Residue site clustering
		}
	unsigned short	Site(unsigned short n){
			if(n > N) return 0;
			// return (((n-1)%Len)+1); // Residue set clustering
			return (((n-1)/max_set)+1);	// Residue site clustering
		}
	unsigned short	ResSet(unsigned short n){
			if(n > N) return 0;
			// return (((n-1)/Len)+1); // Residue set clustering
			return (((n-1)%max_set)+1); // Residue site clustering
		}
	BooLean	AddEdge(cti_typ *cti0);
	BooLean	AddEdge(unsigned short n1,unsigned short n2){
			if(n1 > N || n2 > N) return FALSE;
			grf->AddEdge(n1,n2);
			return TRUE;
		}
	vst_typ **RtnCliques( ) { return clique; }
	vst_typ **RtnConcise(Int4 *NmC) { *NmC=NumConcise; return concise; }
	vst_typ **CreateCliques(Int4 limit,Int4 hpsz, Int4 *NumClust, double pcut){
			return CreateCliques(0,limit,hpsz, NumClust, pcut); }
	vst_typ **CreateCliques(FILE *fp,Int4 limit,Int4 hpsz, Int4 *NumClust, double pcut);
	Int4	NumClusters(){ return NumCluster; }
	Int4	NumberConcise(){ return NumConcise; }
	Int4	PutClique(FILE *fp){ return 0; }
	Int4	SizeClique(Int4 c){ return 0; }
	Int4	VertexClique(Int4 n){ return 0; }
	unsigned short	MaxNumResSets( ){ return max_set; }
	void	ConciseCliques(Int4,Int4,Int4*);
private:
	void	Init(Int4 len, Int4 ms, Int4 TrueNumNodes);
	void	Free();
	grf_typ		*grf;
	cti_typ		**cti;
	set_typ		*NodeSet;
	unsigned short	*rank;
	unsigned short	*edges;
        // Node:  Node = (set * N) + site); // range 1...(set_num +1) *N
	unsigned short	Len;	// Total number of sites (sequence length).
	unsigned short	N;	// Total number of nodes.
	unsigned short	TrueN;	// True total number of nodes.
	unsigned short	max_set;	// max number of residue sets.
	// Int4	Rank;
	Int4	NumEdges;
	vst_typ	**clique;
	vst_typ	**concise;
	Int4	NumCluster,NumConcise;
};

#endif

