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

#include "digraph.h"

/* Initialize to graph with N vertices and no edges. */
dg_type Dgraph(Int4 N, Int4 M)
{
	dg_type G;
	vertex	u;

	NEW(G,1,digraph);
	G->N = N; G->M = M+1;
	NEW(G->fin,N+1,edge);
	NEW(G->fout,N+1,edge);
	NEW(G->edges,M+1,dgedge);
	G->n = N; G->m = 0;
	for(u = 1; u <= N; u++) G->fout[u] = G->fin[u] = 0;
	return G;
}

dg_type NilDgraph(dg_type G)
{free(G->fin);free(G->fout);free(G->edges);free(G);return (dg_type) NULL;}

/* Read up to first occurrence of c or EOF. */
static char cflush(char c, FILE *f)
{
	char c1;
	while ((c1 = getc(f)) != EOF && c1 != c) {} return c1;
}

/* Read one edge from *f into edges[e]. Return TRUE on success, else FALSE. */
BooLean	GetEdgeDgraph(edge e, FILE* f, dg_type G)
{
	if (cflush('(',f) == EOF) return FALSE;
	fscanf(f,"%d",&G->edges[e].t);
	if (cflush(',',f) == EOF) return FALSE;
	fscanf(f,"%d",&G->edges[e].h);
	if (cflush(',',f) == EOF) return FALSE;
	return TRUE;
}

/* Join u and v with edge of given weight. Return edge number. */
Int4	JoinDgraph(vertex u, vertex v, dg_type G)
{
	Int4 m;
	G->m++; m = G->m;
	G->edges[m].t = u; 
	G->edges[m].h = v; 
	G->edges[m].tnext = G->fout[u]; G->fout[u] = m;
	G->edges[m].hnext = G->fin[v]; G->fin[v] = m;
	return m;
}

/* Build adjacency lists. */
void	bldadj_digraph(dg_type G)
{
	vertex	u;
	edge 	e;
	for(u = 1; u <= G->n; u++) G->fout[u] = G->fin[u] = 0;
	for(e = G->m; e >= 1; e--) {
		G->edges[e].tnext = G->fout[G->edges[e].t];
		G->fout[G->edges[e].t] = e;
		G->edges[e].hnext = G->fin[G->edges[e].h];
		G->fin[G->edges[e].h] = e;
	}
}

/* Get graph from f. */
BooLean GetDgraph(FILE* f, dg_type G)
{
	if(fscanf(f,"%d",&G->n) == EOF) return FALSE;
	if(G->n > G->N) digraph_error("input graph too big");
	for(G->m = 1; GetEdgeDgraph(G->m,f,G); G->m++)
	G->m--; bldadj_digraph(G);
	return TRUE;
}

/* Put graph out on f. */
void 	PutDgraph(FILE* f, dg_type G)
{
	edge e;
	fprintf(f,"%d\n",G->n);
	for(e = 1; e <= G->m; e++) {
		fprintf(f,"(%2d,%2d)  ", TailDgraph(e,G),
			HeadDgraph(e,G));
		if (e == G->m || e%5 == 0) putc('\n',f);
	}
	putc('\n',f);
}

dg_type RandDgraph(Int4 n, double pd)
{
/* Generate a random graph on n1 vertices with edge probability pd. */
	dg_type G = Dgraph(n,n*n);
	double p;
	vertex	u,v;

	G->m = 1;
	p = pd * INT4_MAX;
	for(u = 1; u <= n; u++) {
		for(v = 1; v <= n; v++) {
			if (u != v && Random() <= p) {
			   G->edges[G->m].t = u; 
			   G->edges[G->m].h = v;
			}
		}
	}
	G->m--; bldadj_digraph(G);
	return G;
}

/* sort the edges topologically. */
void	ToposortDigraph(dg_type G, Int4* topo)
{
        Int4 i; vertex u,v; edge e;
	l_type q = List(G->n+1);
        Int4* nin;
	NEW(nin,G->n+1,Int4);
        for(u = 1; u <= G->n; u++) {
                nin[u] = 0;
                for(e=FirstInDgraph(u,G);e!=0;e=NextInDgraph(e,G))
			nin[u]++;
                if(nin[u] == 0) AppendList(u,q);
        }
        i = 0;
        while((u=ItemList(1,q)) != 0) {
                RemoveList(1,q); topo[++i] = u;
                for(e=FirstOutDgraph(u,G);e!=0;e=NextOutDgraph(e,G)){
                        v = HeadDgraph(e,G);
                        if ((--(nin[v])) == 0) AppendList(v,q);
                }
        }
        if (i < G->n) digraph_error("toposort... graph has cycle");
	else {
	  for(u=1;u<= G->n; u++) fprintf(stderr,"%d ",topo[u]);
	  fprintf(stderr,"\n");
	}
}

void	digraph_error(char *s){fprintf(stderr,"digraph: %s\n",s);exit(1);}

