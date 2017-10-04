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

#if !defined (DGRAPH)
#define DGRAPH
#include "stdinc.h"
#include "list.h"
#include "random.h"
/* Header file for data structure representing directed graph. */

typedef Int4 vertex;
typedef Int4 edge;
typedef Int4 weight;

typedef struct {
		vertex 	h,t;	/* head and tail of the edge */
		edge	hnext;	/* link to next edge incident to head */
		edge	tnext;	/* link to next edge incident to tail */
} dgedge;
	
typedef struct {
	Int4	N, M;		/* max number of vertices and edges */
	edge	*fout;		/* fout[v] is first edge out of v */
	edge	*fin;		/* fin[v] is first edge into v */
	dgedge	*edges;
	Int4	n,m;			/* number of vertices and edges */
} digraph;

typedef digraph* dg_type;
/************************* private *****************************/

Int4     JoinDgraph(vertex u, vertex v, dg_type G);
void    bldadj_digraph(dg_type G);
void    digraph_error(char *s);

/************************* PUBLIC *****************************/
BooLean GetEdgeDgraph(edge e, FILE* f, dg_type G);
dg_type Dgraph(Int4 N, Int4 M);
dg_type NilDgraph(dg_type G);
static char cflush(char c, FILE *f);
void    PutDgraph(FILE* f, dg_type G);
dg_type RandDgraph(Int4 n, double pd);
void    ToposortDigraph(dg_type G, Int4* topo);
Int4     Digraph_error(char *s);
BooLean GetDgraph(FILE* f,dg_type G);

/************************* MACROS *****************************/
/* Return first edge out of v. */
#define FirstOutDgraph(v,G)	((G)->fout[(v)])

/* Return first edge into v. */
#define FirstInDgraph(v,G)	((G)->fin[(v)])

/* Return next edge out of tail(e), following e. */
#define NextOutDgraph(e,G)	((G)->edges[(e)].tnext)

/* Return next edge into head(e), following e. */
#define NextInDgraph(e,G)	((G)->edges[(e)].hnext)

/* Return head of e. */
#define HeadDgraph(e,G)	((G)->edges[(e)].h)

/* Return tail of e. */
#define	TailDgraph(e,G)	((G)->edges[(e)].t)

/* Return weight of e.*/
#define	DgraphN(G)		((G)->N)
#define nDgraph(G)		((G)->n)

#endif
