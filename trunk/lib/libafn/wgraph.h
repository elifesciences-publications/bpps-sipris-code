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

// Header file for data structure representing undirected graph with
// weighted edges.

typedef int vertex;
typedef int edge;
typedef int weight;

class wgraph {
	int	N,M;			// max number of vertices and edges
	edge	*firstedge;		// firstedge[v] is first edge at v
	struct wgedge {
		vertex 	l,r;		// endpoints of the edge
		weight	wt;		// edge weight
		edge	lnext;		// link to next edge incident to l
		edge	rnext;		// link to next edge incident to r
	} *edges;

	void	bldadj();		// build adjacency lists
	bit	getedge(edge,FILE*);	// get edge from file
public:		wgraph(int=100,int=1000);
		~wgraph();
	int	n,m;			// number of vertices and edges
	edge	first(vertex);		// return first edge incident to v
	edge	next(vertex,edge);	// return next edge of v, following e
	vertex	left(edge);		// return left endpoint of e
	vertex	right(edge);		// return right endpoint of e
	vertex	mate(vertex,edge); 	// return other endpoint of e
	weight	w(edge);		// return weight of e
	int	join(vertex,vertex,weight); // join two vertices with an edge
	bit	get(FILE*);		// get graph from file
	void	put(FILE*);		// put graph to file
	void	rgraph(int,double,int);	// generate random graph
	void	esort();		// sort edges by weight
};

// Return first edge incident to x.
inline edge wgraph::first(vertex v) { return firstedge[v]; };

// Return next edge incident to v, following e.
inline edge wgraph::next(vertex v,edge e)
{ return edges[e].l == v ? edges[e].lnext : edges[e].rnext; }

// Return left endpoint of e.
inline vertex wgraph::left(edge e) { return edges[e].l; }

// Return right endpoint of e.
inline vertex wgraph::right(edge e) { return edges[e].r; }

// Return other endpoint of e.
inline vertex wgraph::mate(vertex v, edge e)
{ return edges[e].l == v ? edges[e].r : edges[e].l ; }

// Return weight of e.
inline weight wgraph::w(edge e) { return edges[e].wt; }

