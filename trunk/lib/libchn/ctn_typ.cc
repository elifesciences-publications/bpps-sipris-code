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

#include "ctn_typ.h"

void	ctn_typ::Init(Int4 len, Int4 ms, Int4 TrueNumNodes)
{
	Len=len; max_set=ms; N = Len*max_set;
	TrueN=TrueNumNodes;
	NEWP(cti,N+5,cti_typ);
	NEW(edges,N+5,unsigned short);
	NEW(rank,N+5,unsigned short);
	NEW(NodeSet,N+5,set_typ);
	grf = new grf_typ(N+2); // vertex numbering is from 0.
	concise=clique=0; NumConcise=NumCluster=0;
}

void	ctn_typ::Free( )
{
	Int4 i;
	for(i=0; i <= N; i++){
		// if(cti[i]) delete cti[i];
	} free(cti);
	free(edges); free(rank); free(NodeSet);
	delete grf;
	if(clique){
	for(i=1; i <= NumCluster; i++){ delete clique[i]; }
	free(clique);
	}
	if(concise){
		for(i=1; i <= NumConcise; i++){ delete concise[i]; }
		free(concise);
	}
}

BooLean	ctn_typ::AddEdge(cti_typ *cti0){
        Int4 n1=this->Node(cti0->SiteI,cti0->rsidI);
        Int4 n2=this->Node(cti0->SiteJ,cti0->rsidJ);
	// Save set information about this node...
	if(NodeSet[n1] == 0) NodeSet[n1] = cti0->SetI();
	if(NodeSet[n2] == 0) NodeSet[n2] = cti0->SetJ();
        // edge[n1][n2]=edge[n2][n1]=cti0;
        // if(sstList[n1]==0) sstList[n1]=cti0->sstI;
        // if(sstList[n2]==0) sstList[n2]=cti0->sstJ;
        // CnTabList[rank]=cti0;
        return this->AddEdge(n1,n2);
}

vst_typ **ctn_typ::CreateCliques(FILE *fp,Int4 limit,Int4 hpsz, Int4 *NumClust, double pcut)
{
	clique=grf->Bron_Kerbosch_cluster(fp,limit,hpsz,NumClust,pcut,NodeSet,TrueN); 
	// if(fp) grf->Put(fp);
	if(clique) assert(NumClust > 0); 
	NumCluster=*NumClust;
	return clique;
}

void    ctn_typ::ConciseCliques(Int4 MinClique,Int4 MaxRank,Int4 *Rank)
// Return clique of sites only...
{
  if(!clique) return;
  else {
    Int4        res,j,n,m,node;
    BooLean     found;
    NEWP(concise,NumCluster+3,vst_typ);
    for(m=0,n=1; n <= NumCluster; n++){
     found=FALSE;
     if(clique[n]->Size() >= MinClique){
        for(j=0; j < clique[n]->Size(); j++){
          node=clique[n]->Vertex(j);
          // res=this->Site(node);
          if(Rank[node] > MaxRank) continue;
          else {
                if(!found){
                        m++; found=TRUE;
                        concise[m]=new vst_typ(N+1);
                } concise[m]->Add(node);
          }
        }
     }
    } NumConcise=m;
  }
}

