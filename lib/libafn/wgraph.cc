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

#include "stdinc.h"
#include "wgraph.h"

wgraph::wgraph(int N1, int M1) {
	N = N1; M = M1 + 1;
	firstedge = new edge[N+1];
	edges = new wgedge[M+1];
	n = N; m = 0;
	for (vertex u = 1; u <= n; u++) firstedge[u] = Null;
}

wgraph::~wgraph() { delete firstedge; delete edges; }

char cflush(char c, FILE *f) {
// Read up to first occurrence of c or EOF.
	char c1; while ((c1 = getc(f)) != EOF && c1 != c) {} return c1;
}

bit wgraph::getedge(edge e, FILE* f) {
// Read one edge from *f into edges[e]. Return true on success, else false.
	if (cflush('(',f) == EOF) return false;
	fscanf(f,"%d",&edges[e].l);
	if (cflush(',',f) == EOF) return false;
	fscanf(f,"%d",&edges[e].r);
	if (cflush(',',f) == EOF) return false;
	fscanf(f,"%d",&edges[e].wt);
	if (cflush(')',f) == EOF) return false;
	return true;
}

int wgraph::join(vertex u, vertex v, weight W) {
// Join u and v with edge of given weight. Return edge number.
	edges[++m].l = u; edges[m].r = v; edges[m].wt = W;
	edges[m].lnext = firstedge[u]; firstedge[u] = m;
	edges[m].rnext = firstedge[v]; firstedge[v] = m;
	return m;
}

void wgraph::bldadj() {
// Build adjacency lists.
	for (vertex u = 1; u <= n; u++) firstedge[u] = Null;
	for (edge e = m; e >= 1; e--) {
		edges[e].lnext = firstedge[edges[e].l];
		firstedge[edges[e].l] = e;
		edges[e].rnext = firstedge[edges[e].r];
		firstedge[edges[e].r] = e;
	}
}

bit wgraph::get(FILE* f) {
// Get graph from f.
	if (fscanf(f,"%d",&n) == EOF) return false;
	if (n > N) fatal("wgraph::get: input graph too big");
	for (m = 1; getedge(m,f); m++)
		if (m > M) fatal("wgraph::get: input graph too big");
	m--; bldadj();
	return true;
}

void wgraph::put(FILE* f) {
// Put graph out on f.
	fprintf(f,"%d\n",n);
	for (edge e = 1; e <= m; e++) {
		fprintf(f,"(%2d,%2d,%2d)  ",left(e),right(e),w(e));
		if (e == m || e%5 == 0) putc('\n',f);
	}
	putc('\n',f);
}

void wgraph::rgraph(int n1, double pd, int maxcost) {
// Generate a random graph on n1 vertices with edge probability pd.
// Generate random edge costs uniformly in the interval [1,maxcost].
	n = n1; m = 1;
	if (n > N) fatal("wgraph::rgraph: graph too big");
	double p = pd * BIGINT;
	for (vertex u = 1; u <= n; u++) {
		for (vertex v = u+1; v <= n; v++) {
			if (random() <= p) {
				if (m > M) fatal("wgraph::rgraph: graph too big");
				edges[m].l = u; edges[m].r = v;
				edges[m++].wt = randint(1,maxcost);
			}
		}
	}
	m--; bldadj();
}

void wgraph::esort() {
// Sort the edges by cost.
	for (int i = 2; i <= m; i++) {
		int j = i;
		while (j > 1 && edges[j].wt < edges[j-1].wt) {
			edges[0] = edges[j];
			edges[j] = edges[j-1];
			edges[j-1] = edges[0];
			j--;
		}
	}
	bldadj();
}
