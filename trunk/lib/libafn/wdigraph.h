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

#if !defined (WDGRAPH)
#define WDGRAPH
#include "stdinc.h"
#include "list.h"
#include "random.h"
#include "set_typ.h"
/* Header file for data structure representing directed graph with
   weighted edges. */

// typedef Int4 node_typ;
// typedef Int4 edge_typ;

typedef struct {
		Int4 	h,t;	/* head and tail of the edge */
		Int4	wt;	/* edge weight */
		Int4	hnext;	/* link to next edge incident to head */
		Int4	tnext;	/* link to next edge incident to tail */
} wdgedge;
	
typedef struct {
	Int4	N, M;		/* max number of vertices and edges */
	Int4	n,m;		/* number of vertices and edges */
	Int4	*fout;		/* fout[v] is first edge out of v */
	Int4	*fin;		/* fin[v] is first edge into v */
	wdgedge	*edges;
} wdigraph;

typedef wdigraph* wdg_typ;

/************************* private *****************************/
void    bldadj_wdigraph(wdg_typ G);
void	esortWdgraph(wdg_typ G);
Int4	find_cycle_vertex(Int4 n, Int4 path[]);
void	wdigraph_error(char *s);
void	TreeDFS_WDG(FILE *fp, Int4 v, set_typ S, wdg_typ T); // for printing Newick format tree.

/************************* PUBLIC *****************************/
wdg_typ MakeWdgraph(Int4 N, Int4 M);
void    AddVerticesWdgraph(Int4 x, wdg_typ G);
wdg_typ MkWdgraph(Int4 N, Int4 M);
void	PrintNewickTreeWDG(FILE *fp,Int4 root, wdg_typ T);
void    ClearWdgraph(wdg_typ G);
Int4	JoinWdgraph(Int4 u, Int4 v, Int4 W, wdg_typ G);
void    PutWdgraph(FILE* f, wdg_typ G);
void    ToposortWdigraph(wdg_typ G, Int4* topo);
wdg_typ NilWdgraph(wdg_typ G);
Int4	ShortestPathWdgraph(wdg_typ G, Int4 s, Int4 *path, Int4 *dist);
void    TopoScanWdigraph(wdg_typ G, Int4 s, Int4 *path, Int4 *dist);
wdg_typ RtnMinSpanningTree(FILE *fp, Int4 Root, wdg_typ Grph);

//************************ New Routines *************************
Int4    NumEdgesOut(Int4 v, wdg_typ G);
Int4    NumEdgesIn(Int4 v, wdg_typ G);
BooLean IsJoined(Int4 i, Int4 j, wdg_typ G);
BooLean IsNodeOkay(Int4 i,wdg_typ G);

/************************* MACROS *****************************/
#define FirstOutWdgraph(v,G)	((G)->fout[(v)])
#define FirstInWdgraph(v,G)	((G)->fin[(v)])

/* Return next edge out of tail(e), following e. */
#define NextOutWdgraph(e,G)	((G)->edges[(e)].tnext)

/* Return next edge into head(e), following e. */
#define NextInWdgraph(e,G)	((G)->edges[(e)].hnext)

/* Return head of e. */
#define HeadWdgraph(e,G)	((G)->edges[(e)].h)

/* Return tail of e. */
#define	TailWdgraph(e,G)	((G)->edges[(e)].t)

/* Return weight of e.*/
#define WeightWdgraph(e,G)	((G)->edges[(e)].wt)
#define SetWeightWdgraph(w,e,G)	((G)->edges[(e)].wt=(w))
#define	WdgraphN(G)		((G)->N)
#define	WdgraphM(G)		((G)->M)
#define nWdgraph(G)		((G)->n)
#define mWdgraph(G)		((G)->m)

#endif
