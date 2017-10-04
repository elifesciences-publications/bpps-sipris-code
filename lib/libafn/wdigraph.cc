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

#include "wdigraph.h"

/**************** reuseable graph - not yet used *****************/
wdg_typ MakeWdgraph(Int4 N, Int4 M)
/* Initialize to graph with zero vertices and no edges. */
{
	wdg_typ G;

	G = MkWdgraph(N, M);
	G->n = 0;
	return G;
}

void	AddVerticesWdgraph(Int4 x, wdg_typ G)
{
	G->n += x;
	if(G->n > G->N) wdigraph_error("maximum # vertices exceeded");
}

void	ClearWdgraph(wdg_typ G)
/* Reinitialize graph to no edges. */
{
	Int4	u;

	for(u = 1; u <= G->n; u++) G->fout[u] = G->fin[u] = 0;
	G->n = 0; G->m = 0;
}

/**************** end reuseable graph - not yet used *****************/

wdg_typ MkWdgraph(Int4 N, Int4 M)
/* Initialize to graph with N vertices and no edges. */
{
	wdg_typ G;

	assert(N > 0 && N < INT4_MAX);
	assert(M > 0 && M < INT4_MAX);
	NEW(G,1,wdigraph);
	G->N = N; G->M = M+1;
	NEW(G->fin,N+1,Int4);
	NEW(G->fout,N+1,Int4);
	NEW(G->edges,M+1,wdgedge);
	G->n = N; G->m = 0;
	return G;
}

wdg_typ NilWdgraph(wdg_typ G)
{free(G->fin);free(G->fout);free(G->edges);free(G);return (wdg_typ) NULL;}

Int4	JoinWdgraph(Int4 u, Int4 v, Int4 W, wdg_typ G)
/* Join u and v with Int4 of given weight W. Return Int4 number. */
{
	Int4 m;

	assert(G->m < G->M);
	G->m++; m = G->m;
	G->edges[m].t = u; 
	G->edges[m].h = v; 
	G->edges[m].wt = W;
	G->edges[m].tnext = G->fout[u]; G->fout[u] = m;
	G->edges[m].hnext = G->fin[v]; G->fin[v] = m;
	return m;
}

void 	PutWdgraph(FILE* f, wdg_typ G)
/* Put graph out on f. */
{
	Int4 e;
	fprintf(f,"%d\n",G->n);
	for(e = 1; e <= G->m; e++) {
		fprintf(f,"(%2d,%2d,%2d)  ", TailWdgraph(e,G),
			HeadWdgraph(e,G),WeightWdgraph(e,G));
		if (e == G->m || e%5 == 0) putc('\n',f);
	}
	putc('\n',f);
}

/************************ PrintNewickTree ************************/
void    PrintNewickTreeWDG(FILE *fp,Int4 root, wdg_typ T)
{ // print tree in Newick format.
        // fprintf(fp,"("); TreeDFS(fp, root, T); fprintf(fp,")%d;\n",root);
        // PutWdgraph(stderr,T);
	set_typ S=MakeSet(T->N+1); ClearSet(S);
        TreeDFS_WDG(fp, root, S, T); fprintf(fp,";\n");
	NilSet(S);
}

void    TreeDFS_WDG(FILE *fp, Int4 v, set_typ S, wdg_typ T)
{
        // start from root, recursively print tree in Newick format.
        Int4    w,e,i,n=NumEdgesOut(v,T);
	if(MemberSet(v,S)) print_error("TreeDFS_WDG( ): input graph not a tree (cycle found).");
	AddSet(v,S);
        if(n > 0){
           fprintf(fp,"(");
           for(i=0,e=FirstOutWdgraph(v,T);e!=0;e=NextOutWdgraph(e,T)){
                i++; w = HeadWdgraph(e,T);      // a child of v.
                TreeDFS_WDG(fp,w,S,T);
                if(i < n) fprintf(fp,",");
           } fprintf(fp,")");
        } fprintf(fp,"%d_Set%d",v,v);      // output for HyperPartition.
}

/************************ Shortest Path Algorithms ************************
  bfs_shortest paths. see Tarjan.
/**************************************************************************/

wdg_typ RtnMinSpanningTree(FILE *fp, Int4 Root, wdg_typ Grph)
{
    assert(Grph != 0); assert(Root > 0);
    Int4  wt,i,r,N=WdgraphN(Grph),E=WdgraphM(Grph);
    Int4  *path,*dist; NEW(path,N+3,Int4); NEW(dist,N+3,Int4);
    if(ShortestPathWdgraph(Grph,Root,path,dist)) print_error("cycle in graph");
    wdg_typ Tree=MkWdgraph(N,E);
    for(i=1; i <= N; i++){
        // if(dist[i] > 0.0) continue;
        if(dist[i] < 0.0) continue;
        // if(NumEdgesOut(i,Grph) != 0) continue;
        if(fp) fprintf(fp," node %d: ",i);
        for(r=i; r != Root && r != 0; r = path[r]){
             if(fp) fprintf(fp,"%d ->",r);
             wt = dist[r]-dist[path[r]];
             if(!IsJoined(path[r],r,Tree)) JoinWdgraph(path[r],r,wt,Tree); // should end when path[r] == Root.
        } if(fp) fprintf(fp," %d (dist = %d)\n",Root,dist[i]);
    } if(fp) fprintf(fp,"\n\n");
    free(dist); free(path);
    return Tree;
}

Int4	ShortestPathWdgraph(wdg_typ G, Int4 s, Int4 *path, Int4 *dist)
/************************************************************************
// Bread First Scanning shortest paths algorithm [O(mn) time]:
//
// Find a shortest path tree rooted at s in graph G.
// Return the parent pointers of the tree in array p,
// and the distances from s to every vertex in array dist.
// If there is no path from s to v, then dist[v]=BIGINT and p[v]=Null.
// If there is a negative cycle reachable from s, then the parent
// pointers contain a cycle, and one of the vertices on a cycle
// is returned.
// ("What?  You put bugs in your programs?  What for?" :-)

   WARNING: NEED 'r != NULL' for path checking?  DEBUG!!!

        if(ShortestPathWdgraph(G, s, path, dist)){
                print_error("cycle in graph");
        }
        fprintf(stderr,"\nstart ");
        for(r = path[e]; r != s && r != NULL; r = path[r]){
                fprintf(stderr,"-> %d",r);
        } fprintf(stderr," (dist = %d)\n\n",-dist[e]);

/***********************************************************************/
{
        Int4	pass,v, w, last,e;
	l_type	q = List(G->n+1);

        for(v=1; v<=G->n; v++) { dist[v]=INT4_MAX; path[v]=0; }
        dist[s]=0; AppendList(s,q);
        pass = 0; last = s;
        while((v=ItemList(1,q)) != 0) {
           RemoveList(1,q); 
           for(e=FirstOutWdgraph(v,G);e!=0;e=NextOutWdgraph(e,G)){
		w = HeadWdgraph(e,G);
		// w = TailWdgraph(e,G);
		// if(dist[w] == INT4_MAX || 
		if((dist[v]+WeightWdgraph(e,G)) < dist[w]) {
                                path[w] = v;
                                dist[w] = dist[v] + WeightWdgraph(e,G);
                                if(!MbrList(w,q)) AppendList(w,q);
		}
           }
           if(v==last && ItemList(1,q) != 0) {
                        pass++; last = LastItemList(q);
           }
           if(pass == G->n) {
		/******************************************************
                    There is a negative cycle.  It can be found by
                    looking for a cycle among the parent pointers.
                    Return in negative_cycle_Int4 the number of
                    a Int4 which is in a cycle of parent pointers. 
		/******************************************************/
		NilList(q);
                return find_cycle_vertex(G->n, path);
           }
        } NilList(q);
	return 0;
}

Int4	find_cycle_vertex(Int4 n, Int4 path[])
/*********************************************************************
// Given parent pointers path defined on vertices {1, .., n},
// return a Int4 which is in a cycle of parent pointers.
// This routine will go into an infinite loop if there is no cycle.
/*********************************************************************/
{
        Int4	*mark,v,start = 1;

 	NEW(mark,n+1,Int4);
        for(v = 1; v <= n; v++) mark[v] = 0;

        while(TRUE) {
                /** Find an unmarked starting point to look for a cycle. **/
                for(; start <= n && mark[start] != 0; start++) {}
                if(start > n) { free(mark); return INT4_MAX; }

                mark[start] = start;
                for(v=path[start]; v != 0&& mark[v] == 0; v = path[v]) {
                        mark[v] = start;
                }
                /** Return start if we found a cycle. **/
                if (v != 0 && mark[v] == start) { free(mark); return v; }
        }
}

void	TopoScanWdigraph(wdg_typ G, Int4 s, Int4 *path, Int4 *dist)
/*******************************************************************
  Find shortest path in O(m) time using topological scanning;
  (see Tarjan, Data Structures & Network Algorithms, p 88-89).
  WARNING: G must be an acyclic graph.
/*******************************************************************/
{
        Int4	i,v,w,e,*topo;

	NEW(topo,G->n+2,Int4);
	ToposortWdigraph(G, topo);
	if(topo[1] != s) wdigraph_error("TopoScanWdigraph(): input error");
        for(v=1; v<=G->n; v++) { dist[v]=INT4_MAX; path[v]=0; }
	dist[s]=0;

	for(i=1; i <= G->n; i++){
	   v = topo[i];
           for(e=FirstOutWdgraph(v,G);e!=0;e=NextOutWdgraph(e,G)){
		w = HeadWdgraph(e,G);
		if((dist[v]+WeightWdgraph(e,G)) < dist[w]) {
                                path[w] = v;
                                dist[w] = dist[v] + WeightWdgraph(e,G);
		}
           }
        } 
	free(topo);
}

void	ToposortWdigraph(wdg_typ G, Int4 *topo)
/* sort the vertices topologically. */
{
        Int4	i,u,v,e;
	l_type	q = List(G->n+1);

        Int4* nin;
	NEW(nin,G->n+1,Int4);
        for(u = 1; u <= G->n; u++) {
                nin[u] = 0;
                for(e=FirstInWdgraph(u,G);e!=0;e=NextInWdgraph(e,G))
			nin[u]++;
                if(nin[u] == 0) AppendList(u,q);
        }
        i = 0;
        while((u=ItemList(1,q)) != 0) {
                RemoveList(1,q); topo[++i] = u;
                for(e=FirstOutWdgraph(u,G);e!=0;e=NextOutWdgraph(e,G)){
                        v = HeadWdgraph(e,G);
                        if ((--(nin[v])) == 0) AppendList(v,q);
                }
        }
        if (i < G->n) wdigraph_error("toposort... graph has cycle");
	else {
	  /******
	  for(u=1;u<= G->n; u++) fprintf(stderr,"%d ",topo[u]);
	  fprintf(stderr,"\n"); /******/
	}
	free(nin); NilList(q);
}

/******************************** private ********************************/
void	bldadj_wdigraph(wdg_typ G)
/* Build adjacency lists. */
{
	Int4	u;
	Int4 	e;

	for(u = 1; u <= G->n; u++) G->fout[u] = G->fin[u] = 0;
	for(e = G->m; e >= 1; e--) {
		G->edges[e].tnext = G->fout[G->edges[e].t];
		G->fout[G->edges[e].t] = e;
		G->edges[e].hnext = G->fin[G->edges[e].h];
		G->fin[G->edges[e].h] = e;
	}
}

void esortWdgraph(wdg_typ G)
/* Sort the edges by cost.*/
{
	Int4 i,j;
	for(i = 2; i <= G->m; i++) {
		j = i;
		while (j > 1 && G->edges[j].wt < G->edges[j-1].wt) {
			G->edges[0] = G->edges[j];
			G->edges[j] = G->edges[j-1];
			G->edges[j-1] = G->edges[0];
			j--;
		}
	}
	bldadj_wdigraph(G);
}

void	wdigraph_error(char *s){fprintf(stderr,"wdigraph: %s\n",s);exit(1);}

//********************************************************************************

Int4    NumEdgesIn(Int4 v, wdg_typ G)
{ Int4 e,n=0; for(e=FirstInWdgraph(v,G);e!=0;e=NextInWdgraph(e,G)){ n++; } return n; }

Int4    NumEdgesOut(Int4 v, wdg_typ G)
{ Int4 e,n=0; for(e=FirstOutWdgraph(v,G);e!=0;e=NextOutWdgraph(e,G)){ n++; } return n; }

BooLean IsJoined(Int4 i, Int4 j, wdg_typ G)
// i == head; j == tail.
{
       Int4    h,t,e;
       for(e=FirstOutWdgraph(i,G);e!=0;e=NextOutWdgraph(e,G)){
            h = HeadWdgraph(e,G); t = TailWdgraph(e,G); if(t==i && h==j) return TRUE;
        } return FALSE;
}

BooLean IsNodeOkay(Int4 i,wdg_typ G) { assert(i > 0 && i <= nWdgraph(G)); return TRUE; }


