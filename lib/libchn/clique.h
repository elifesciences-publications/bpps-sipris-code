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

#ifndef _CLIQUE_H_
#define _CLIQUE_H_

#include "stdinc.h"
#include "afnio.h"
#include "mheap.h"
#include "set_typ.h"
#include "alphabet.h"
#include "sequence.h"
#include "probability.h"	// for CumHypGeomProb(N1,N2,n,x);

class vst_typ {	// Vertex set type
  public:	
	vst_typ( ){ assert(!"illegal constructor"); }
	vst_typ(Int4 n) { init(n); }
	~vst_typ( ){ Free(); }
	void	Copy(vst_typ *);
	void	Add(Int4 v);
	void	RmLast( ){ assert(size > 0); size--; }
	Int4	Put(FILE *);
	void    Put(FILE *fp,e_type E, Int4 ClusterSize,a_type AB);
	// void	Put(FILE *,e_type E, a_type AB);
	Int4	Size() { return size; } 
	Int4	Vertex(Int4 i){ assert(i >= 0 && i < N); return vertex[i]; }
	void	BubbleSort( ) {
		 Int4 i,j,t;
        	 for(i=1; i < size; i++){
                   for(j=i;j > 0 && (vertex[j-1] > vertex[j]);j--){
                        t=vertex[j]; vertex[j]=vertex[j-1]; vertex[j-1]=t;
                   }
        	 }
		}
  private:
  	Int4	N;		// maximum number of nodes.
  	Int4	size;		// number of selected nodes.
  	Int4	*vertex;	// selected nodes.
	void	init(Int4 n){ N=n; size=0; MEW(vertex,N+2,Int4); }
	void	Free(){ free(vertex); }
};

class vsh_typ {	// Vertex set heap type
  public:	
	vsh_typ( ){ assert(!"illegal constructor"); }
	vsh_typ(Int4 n,Int4 hs) { init(n,hs); }
	~vsh_typ( ){ Free(); }
	void    Insert(vst_typ *);
	void    Put(FILE *);
	set_typ *ReturnSets(FILE *);
  private:
	vst_typ	**vst;
	mh_type	mH;
	Int4	N,hpsz;
	void	init(Int4 n,Int4 hs);
	void	Free();
};

class grf_typ {		// simple graph type.
  public:	
	grf_typ( ){ assert(!"illegal constructor"); }
	grf_typ(Int4 n){ init(n); }
	~grf_typ(){ Free(); }
	void	AddEdge(Int4,Int4);
	void	AddEdge(Int4,Int4,char,char);
	BooLean	IsEdge(Int4 v1,Int4 v2){
			assert(v1 < N && v2 < N && v1 >= 0 && v2 >= 0);
			return (BooLean) connected[v1][v2];
		}
	void	RmEdge(Int4,Int4);
	void	Put(FILE *fp){ Put(fp,0); }
	void	Put(FILE *fp,set_typ set);
	void	PutWeighted(FILE *fp,set_typ set);
	void	Put(FILE *fp,set_typ set,BooLean UseWts);
	vst_typ	**Bron_Kerbosch(Int4 a,Int4 b,Int4 *c,double d,set_typ *e,unsigned short f){
			return Bron_Kerbosch(0,a,b,c,d,e,f,FALSE);
		}
	vst_typ	**Bron_Kerbosch_cluster(Int4 a,Int4 b,Int4 *c,double d,set_typ *e,unsigned short f)
			{ return Bron_Kerbosch(0,a,b,c,d,e,f,TRUE); }
	vst_typ	**Bron_Kerbosch_cluster(FILE *fp,Int4 a,Int4 b,Int4 *c,double d,set_typ *e,unsigned short f)
			{ return Bron_Kerbosch(fp,a,b,c,d,e,f,TRUE); }
  private:
	vst_typ	**Bron_Kerbosch(FILE *fp,Int4,Int4,Int4 *,double,set_typ*,unsigned short,BooLean);
	Int4	N;
	char	**connected;
	UInt4	*NumEdges;
	void	init(Int4 n);
	void	Free();
};

// vst_typ	*BronKerbosch(Int4 N, char **connected,Int4 hpsz,Int4 limit);

#endif

