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

/**************************************************************
Direct implementation of the Bron-Kerbosch procedure for 
finding all maximum cliques of an undirected graph.  This is
Algorithm 457 of the CACM Collected Algorithms.
/**************************************************************/

#include "clique.h"

//************************* Vertex set type *******************
void	vst_typ::Copy(vst_typ *to)
{
  Int4 loc = size;
  const Int4 *v1 = vertex;
  Int4 *v2 = to->vertex;
  to->size = size;
  while(loc-- > 0) { *v2++=*v1++; }
}

void	vst_typ::Add(Int4 v)
{ assert(size < N && size >= 0); vertex[size++] = v; }

Int4	vst_typ::Put(FILE *fp)
{
  Int4 i;
  fprintf(fp,"Size = %2d found\n", size);
  for (i=0; i<size; i++) { fprintf(fp,"%d ", vertex[i]); }
  fprintf(fp,"\n");
  return 1;
}

void	vst_typ::Put(FILE *fp,e_type E, Int4 ClusterSize,a_type AB)
{
    if(!(LenSeq(E) < N)) return;
    fprintf(fp,"#Vertex cluster:");
    Int4 last=0;
    for(Int4 i=0; i<size; i++) {
        Int4 v=vertex[i];
#if 1
	assert(v > 0);
	v = ((v-1)/ClusterSize) + 1;
#else
	if(v > LenSeq(E)){ v = (v)%LenSeq(E); if(v == 0) v = LenSeq(E); }
#endif
        Int4 s=v+OffSetSeq(E);
	assert(v <= LenSeq(E) && v < N);
        if(v != last) fprintf(fp," %c%d",AlphaChar(ResSeq(v,E),AB),s);
        // else fprintf(fp," %d", vertex[i]);
	last = v;
    } fprintf(fp,"\n");
}
//********************* end of vst_typ ***********************

//********************* begin of vsh_typ ***********************
void	vsh_typ::init(Int4 n,Int4 hs)
{
	hpsz=hs;
	N=n;
	vst = new vst_typ*[hpsz+2];
	for(Int4 r=1; r <=hpsz; r++) vst[r] = new vst_typ(N);
	mH = Mheap(hpsz,3);
}

void	vsh_typ::Put(FILE *fp)
{
	Int4	item;
	for(Int4 n=1; !EmptyMheap(mH); n++){
		double key = MinKeyMheap(mH);
		item = DelMinMheap(mH);
		vst[item]->Put(fp);
	}
}

void	vsh_typ::Free()
{
	for(Int4 r=1; r <=hpsz; r++) delete vst[r];
	delete [] vst;
	NilMheap(mH);
}

void	vsh_typ::Insert(vst_typ *vst0)
{
	double	key=(double)vst0->Size();
	Int4	item=InsertMheap(-key,mH);
	if(item){
		vst0->Copy(vst[item]);
	}
}

set_typ	*vsh_typ::ReturnSets(FILE *fp)
// return array of vertex sets for merged cliques.
{
	Int4	n,m,item;
	set_typ *set;

	NEW(set, hpsz+2, set_typ);
	for(n=1; !EmptyMheap(mH); n++){
		double key = MinKeyMheap(mH);
		item = DelMinMheap(mH);
		if(fp) vst[item]->Put(fp);
		set[n]=MakeSet(N+1);
		for(m=0; m < vst[item]->Size(); m++){
			AddSet(vst[item]->Vertex(m),set[n]);
		}
	} return set;
}

//********************* end of vsh_typ ***********************

//********************** simple graph type ***********************
void	grf_typ::init(Int4 n)
{
	Int4	r,s;
	N=n; connected=0;
	NEWP(connected,N+3,char);
	NEW(NumEdges,N+3,UInt4);
        for(r=0; r <= N; r++){ NEW(connected[r],N+3,char); }
	for(r=0; r <= N; r++) connected[r][r]=1;
}

void	grf_typ::Free( )
{
        for(Int4 r=0; r <= N; r++) free(connected[r]);
	free(connected); free(NumEdges); 
}

void	grf_typ::AddEdge(Int4 v1,Int4 v2) { AddEdge(v1,v2,1,1); }

void	grf_typ::AddEdge(Int4 v1,Int4 v2,char wt1,char wt2)
{
	assert(wt1 > 0 && wt1 <= 100); assert(wt2 > 0 && wt2 <= 100);
	if(v1 < 0 || v2 < 0 || v1 >= N || v2 >= N || v1==v2) return;
	if(connected[v1][v2]==0){ NumEdges[v1]++; NumEdges[v2]++; }
	connected[v1][v2]=wt1;
	connected[v2][v1]=wt2;
}

void	grf_typ::Put(FILE *fp,set_typ set) { Put(fp,set,FALSE); } 
void	grf_typ::PutWeighted(FILE *fp,set_typ set) { Put(fp,set,TRUE); } 

void	grf_typ::Put(FILE *fp,set_typ set,BooLean UseWts)
// ' ' = 0%;  '.' = 1-10%;  'o' = 10-25%;  'x' = 25-50%;  '*' = 50-100%
{
	Int4	i,j;

	if(set) assert(SetN(set) > N);
	fprintf(fp,"Graph connectivity:\n");
	for(i=0; i <= N; i++){
	   if(NumEdges[i] < 1) continue;	// skip these rows...
	   if(set && !MemberSet(i,set)) continue; // skipt these rows
	   Int4 n=0;
	   for(j=0; j <= N; j++){
	       if(NumEdges[j] < 1) continue; // skip these columns....
	       if(UseWts){
	          if(i == j) fprintf(fp,"O");
	          else if(connected[i][j]==0){ fprintf(fp," "); }
	          else if(connected[i][j] < 10){ n++; fprintf(fp,"."); }
	          else if(connected[i][j] < 25){ n++; fprintf(fp,"o"); }
	          else if(connected[i][j] < 50){ n++; fprintf(fp,"x"); }
	          else { n++; fprintf(fp,"*"); }
	       } else {
	          if(i == j) fprintf(fp,"o");
	          else if(connected[i][j]){ n++; fprintf(fp,"*"); }
	          else fprintf(fp,".");
	       }
	   } fprintf(fp," %d(%d).\n",i,n);
	} // fprintf(fp,"\n");
}

void	grf_typ::RmEdge(Int4 v1,Int4 v2)
{
	if(v1 < 0 || v2 < 0 || v1 >= N || v2 >= N || v1 == v2) return;
	if(connected[v1][v2]==1){
		connected[v1][v2]=connected[v2][v1]=0;
		NumEdges[v1]--; NumEdges[v2]--;
	}
}

//********************* end of simple graph ***********************

static Int4 extend_version2(char **connected, Int4 *old, Int4 ne, Int4 ce,
			vst_typ *compsub, vst_typ *best,vsh_typ *vsh,Int4 limit)
{
  Int4	newne,newce,i,j,count,pos=0,p,s,sel,minnod,nod,fixp;
  Int4	*New; NEW(New,ce+2,Int4); 

  // Determine each counter value and look for minimum:
  minnod = ce; nod=0;
  for(i=0; i<ce && minnod != 0; i++) {
    p = old[i]; count = 0;
    // Count disconnections:
    for(j=ne; j<ce && count < minnod; j++) {
      if(!connected[p][old[j]]) {
	count++;
	pos = j; // Save position of potential candidate.
      }
    }
    // Test new minimum:
    if(count < minnod) {
      fixp = p; minnod = count;
      if(i<ne) { s = pos; }	// it appears that, before setting pos=0, pos was sometimes uninitialized here...
      else { s = i; /* PreIncr: */ nod = 1; }
    }
  }
  // assert(s != 0);	// added to ensure that pos=0 is not used because I am not sure what effect this would have.
  // If fixed point initially chosen from candidates then 
  // number of diconnections will be preincreased by one.
  // BackTrackCycle:
  for(nod=minnod+nod; nod>=1; nod--) {
    /* Interchange */
    p = old[s];	old[s] = old[ne]; sel = old[ne] = p;

    /* Fill new set "not" */
    newne = 0;
    for(i=0; i<ne; i++) {
      if(connected[sel][old[i]]) { New[newne++] = old[i]; }
    }

    /* Fill new set "cand" */
    newce = newne;
    for(i=ne+1; i<ce; i++){
      if(connected[sel][old[i]]) { New[newce++] = old[i]; }
    }
    compsub->Add(sel);

    if(newce == 0) {
      if(best->Size() < compsub->Size()) {
	// found a max clique:
	compsub->Copy(best);
      }
      if(vsh && compsub->Size() >= limit) vsh->Insert(compsub);
      // if(compsub->Size() > 2) { compsub->Put(stdout); }
    } else {
      if(newne < newce) {
	extend_version2(connected,New,newne,newce,compsub,best,vsh,limit);
      }
    }
    // Remove from compsub:
    compsub->RmLast( ); // compsub->size--;
    // Add to "nod":
    ne++;
    if(nod > 1) {
      // Select a candidate disconnected to the fixed point.
      for(s=ne; connected[fixp][old[s]]; s++) ;

    } /* end Selection */
  } /* end BackTrackCycle */
  free(New); 
  return 1;
}

#if 0	//****************************************************
N total nodes with N1 red nodes and N-N1 black nodes.
Choose n elements at random.  The probability that the
group so chosen will contain x or more red elements is
given by:

        p=CumHypGeomProb(N1,N2,n,x)

If this probability is less than a specific cutoff c, x
is considered significant, in which case we merge the
cliques into a single set. 
(Adjustment for the number of sets?)
#endif	//****************************************************
set_typ	MergeSets(set_typ *set,double pcut)
{
	Int4	s1,s2;
	Int4	N1,N2,n,x;
	set_typ	USet;
	double	p;

	if(set[1]==0) return 0;
	USet=MakeSet(SetN(set[1]));
	CopySet(USet,set[1]);
	for(s1=1; set[s1] != NULL; s1++){
	   N1=CardSet(set[s1]);		// red nodes.
	   N2=SetN(set[s1])-N1;		// black nodes
#if 0
	   if(s1 == 1){
	     fprintf(stderr,"Set %d:\n",s1);
	     PutSet(stderr,set[s1]);
	   }
#endif
// if(s1 > 1) break;
	   for(s2=s1+1; set[s2] != NULL; s2++){
	      n=CardSet(set[s2]);	// n nodes chosen at random
	      x=CardInterSet(set[s1],set[s2]);	// x red out of the n nodes
	      if(s1 == 1 && x > 0){
	        p=CumHypGeomProb(N1,N2,n,x);	// NOT WORKING FOR x==0!!!
#if 0
	        fprintf(stderr,"  %d CumHypGeomProb(%d,%d,%d,%d) = %g\n",
			s2,N1,N2,n,x,p);
#endif
	        if(s1 == 1 && p <= pcut){
	           UnionSet(USet,set[s2]); // n=CardSet(USet);
	        }
	      }
	   }
	} return USet;
}

#include "dsets.h"

set_typ	*DontClusterSets(set_typ *set)
// return a copy of the same set array...
{
	Int4	N,n;

	for(N=0,n=1; set[n] != NULL; n++) N++;
	set_typ	*cluster; NEW(cluster,N+3,set_typ);
	for(n=1; n <= N; n++){ cluster[n] = set[n]; }
	return cluster;
}

set_typ	*ClusterSets(set_typ *set, double pval_cutoff,unsigned short TrueN)
// Merge all related cliques into clusters 
{
	Int4	c,N,n,m,s1,s2,N1,N2,r,x;
	double	p;

	// 1. Set up for algorithm...
	for(N=0,s1=1; set[s1] != NULL; s1++) N++;
	ds_type dst=DSets(N);	// rename ds_type to dst_typ
	// PutDSets(stdout,dst);

	// 2. Find related clusters...
	for(n=1; n < N; n++){
	   N1=CardSet(set[n]);		// red nodes.
	   // N2=SetN(set[n])-N1;		// black nodes (too large?)
	   assert(TrueN >= N1);
	   N2=TrueN-N1;		// black nodes 
	   // fprintf(stderr,"Set %d:\n",n);
	   // PutSet(stderr,set[n]);
	   s1 = findDSets(n,dst);
	   for(m=n+1; m <= N; m++){
	     r=CardSet(set[m]);			// r nodes chosen at random
	     x=CardInterSet(set[n],set[m]);	// x out of the r nodes are red
	     if(x > 0){		// NOT WORKING FOR x==0!!!
		if(r == x){
	           s2=findDSets(m,dst); // then merge these sets...
		   if(s1 != s2){ s1 = linkDSets(s1,s2,dst); }
		}
		// else if((p=CumHypGeomProb(N1,N2,r,x)) <= pval_cutoff){  
		else if(abs(r-x) < 2 && (p=CumHypGeomProb(N1,N2,r,x)) <= pval_cutoff){  
#if 0
	        fprintf(stdout,"  %d CumHypGeomProb(%d,%d,%d,%d) = %g\n",
			m,N1,N2,r,x,p);
#endif
	           s2=findDSets(m,dst);
		   if(s1 != s2){ s1 = linkDSets(s1,s2,dst); }
		}
	     }
	   }
	}

	// 2. Create related clusters...
	set_typ	*cluster;
	NEW(cluster,N+3,set_typ);
	for(c=0,n=1; n <= N; n++){
	   s1 = findDSets(n,dst);
	   if(s1 != n){ 
		UnionSet(set[s1],set[n]); // sets set[s1] == set[s1] U set[n];
	   } else { c++; cluster[c] = set[s1]; }
	} NilDSets(dst);
	return cluster;
}

#if 0	// Not used...
vst_typ	*BronKerbosch(Int4 N, char **connected,Int4 hpsz,Int4 limit)
// The input graph is expected 
{
  Int4 c,*ALL; 
  vsh_typ	*vsh=0;
  double	pcut=0.00001;

  if(hpsz > 0){ vsh = new vsh_typ(N,hpsz); }
  vst_typ *best=new vst_typ(N);
  vst_typ *compsub=new vst_typ(N);
  MEW(ALL,N+2,Int4);
  for(c=0; c<N; c++) { ALL[c] = c; }
  extend_version2(connected, ALL, 0, N, compsub,best,vsh,limit);
  // if(vsh){ vsh->Put(stderr); delete vsh; }
  if(vsh){
	// set_typ *set=vsh->ReturnSets(stderr); delete vsh; 
	set_typ *set=vsh->ReturnSets(0); delete vsh; 
	set_typ USet = MergeSets(set,pcut);
	fprintf(stderr,"Merged main set (%d nodes):\n",CardSet(USet));
	while(best->Size() > 0) best->RmLast();
	for(c=0; c<N; c++) if(MemberSet(c,USet)) best->Add(c);
	// PutSet(stderr,USet); 
	NilSet(USet);
    for(c=1; set[c] != NULL; c++){
	// PutSet(stderr,set[c]);
	NilSet(set[c]);
    } free(set);
  }
  delete compsub;
  free(ALL);
  return best;
}
#endif

vst_typ	**grf_typ::Bron_Kerbosch(FILE *fp,Int4 limit,Int4 hpsz, Int4 *NumClust,
	double pcut,set_typ *NodeSet,unsigned short TrueN, BooLean WithClustering)
// The input graph is expected 
// return multiple clique clusters
{
  Int4		s,C,c,*ALL; 
  vsh_typ	*vsh=0;
  // double	pcut=0.01; 

  if(hpsz > 0){ vsh = new vsh_typ(N,hpsz); }
  vst_typ **cluster=0;
  vst_typ *compsub=new vst_typ(N);
  vst_typ *best=new vst_typ(N);
  NEW(ALL,N+2,Int4);
  for(c=0; c<N; c++) { ALL[c] = c; }
  extend_version2(connected, ALL, 0, N, compsub,best,vsh,limit);
  // if(vsh){ vsh->Put(stderr); delete vsh; }
  // if(vsh){ vsh->Put(stderr); }	// this empties the heap so nothing is left...
  if(vsh){
    // set_typ *set=vsh->ReturnSets(stderr); delete vsh; 
    set_typ *set=vsh->ReturnSets(0); delete vsh; 
    for(c=0,s=1; set[s] != NULL; s++) c++;
    if(c > 1) pcut=pcut/(double)((c*(c-1))/2.0);
    set_typ *USet,*thisset;
    Int4 lastC=N; thisset=set;
    do {	// keep clustering sets until can't anymore...
      if(WithClustering) USet=ClusterSets(thisset, pcut,TrueN);
      else { USet=DontClusterSets(thisset); break; }
      for(C=1; USet[C] != NULL; ) C++; C--;
      if(C == lastC){ 
	if(thisset != set) free(thisset); 
	break; 
      } else {
	if(thisset != set) free(thisset);
	lastC=C; thisset = USet;
      } 
    } while(TRUE);
#if 1	// Remove small cliques (i.e., those whose size is <= limit).
    Int4	n,MinClique=limit;
    for(n=0,C=1; USet[C] != NULL; C++){
	   // PutSet(stderr,USet[C]);
	   if(CardSet(USet[C]) < MinClique){
		// fprintf(stderr,"CardSet(%d) = %d <= %d\n",C,CardSet(USet[C]),MinClique);
		USet[C]=0; 	// USet[C] is stored in set[c] array; no need to free.
	   } else n++;
    } C--; 
    if(n < C){ // some deleted;
	 NEW(thisset,n+3,set_typ);
	 for(n=1,c=1; c <= C; c++){
	   if(USet[c]){ thisset[n]=USet[c];  n++; }
	 } free(USet); USet=thisset;
    }
#endif
    Int4 Nsets=0;
    for(C=1; USet[C] != NULL; C++) Nsets++;
    if(USet[1]) NEWP(cluster,Nsets+3,vst_typ);
    // if(USet[1]) NEWP(cluster,N+3,vst_typ);
    for(C=1; USet[C] != NULL; C++){
	if(fp) fprintf(fp,"Clique Cluster %d (%d nodes):\n",C,CardSet(USet[C]));
	if(fp) PutSet(fp,USet[C]);
#if 1	// bug fix: 12-15-10 afn.
  	cluster[C]=new vst_typ(N);
	for(c=0; c < N; c++) if(MemberSet(c,USet[C])) cluster[C]->Add(c);
#endif
	// cluster[C]->Put(stdout,keyE,AB);
    } *NumClust=C-1;
    for(c=1; set[c] != NULL; c++){
	// if(fp) PutSet(fp,set[c]);
	NilSet(set[c]);
    } free(set); free(USet);
  } else *NumClust=0;
  delete compsub;
  delete best;
  free(ALL);
  return cluster;
}

