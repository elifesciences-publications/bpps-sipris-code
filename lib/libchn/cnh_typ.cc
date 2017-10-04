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

#include "cnh_typ.h"

wdg_typ	cnh_typ::PurgeTree(Int4 Root, Int4 max_wt, Int4 &NumRm, wdg_typ Tre)
// 
{
    FILE	*fp=stderr;
    Int4	u,v,wt,nd,e,f,N=WdgraphN(Tre),E=WdgraphM(Tre),n=nWdgraph(Tre);
    double	MaxWt,MinWt,Wt,d;
    wdg_typ	tre=MkWdgraph(N,E);
    MaxWt=(double) max_wt;
    for(NumRm=0,nd=1; nd <= N; nd++){
	if(nd == Root) continue;
	e=FirstInWdgraph(nd,Tre); 
	if(e == 0) continue;
	f=NextInWdgraph(e,Tre); assert(f == 0);
	wt = WeightWdgraph(e,Tre); assert(wt >= 0);
	u = HeadWdgraph(e,Tre); v = TailWdgraph(e,Tre);
	Wt=(double) wt; d = Wt/MaxWt;
	if(d < 0.3334 && FirstOutWdgraph(u,Tre) == 0){ NumRm++; continue; }
	else JoinWdgraph(v,u,wt,tre);
    } return tre;
}

wdg_typ	cnh_typ::RtnConsensusTree(Int4 Root, Int4 max_wt,wdg_typ G)
// remove all but the best 
{
    FILE	*fp=stderr;
    Int4	wt,nd,e,N=WdgraphN(G),E=WdgraphM(G),n=nWdgraph(G);
    Int4	u,v,bst_e,bst_wt;
    wdg_typ	Tre=MkWdgraph(N,E);
    for(nd=1; nd <= N; nd++){
	if(nd == Root) continue;
	bst_e=bst_wt=0;
	for(e=FirstInWdgraph(nd,G);e!=0;e=NextInWdgraph(e,G)){
	    wt = WeightWdgraph(e,G);
	    assert(wt >= 0);
	    if(wt > bst_wt){ bst_e=e; bst_wt=wt; }
	} 
	if(bst_e == 0) continue;
	e = bst_e; // u = HeadWdgraph(e,G); v = TailWdgraph(e,G);
	JoinWdgraph(TailWdgraph(e,G),HeadWdgraph(e,G),bst_wt,Tre);
    } return Tre;
}

const char USAGE[]="usage: cnh_run infile [options]\n\
      NOTE: Main input file = infile.cnh\n\
	  infile.mma corresponds to the main sequence alignment\n\
          infile.cnh should contain > 1 & <= 100 prefix file names corresponding to:\n\
            <prefix>_new.mma \n\
            <prefix>_new.hpt: hierarchy for this prefix file\n\
            <prefix>_new.sets: sets corresponding to \n\
      options:\n\
        -first    run the first tree only\n\
        -x        dummy\n\n";


void	cnh_typ::Init(int argc, char *argv[])
// Run a comparison between two input hierarchies.
{
	char	*Argv[5]; 
	Int4	arg=1,a,i,j,k,n,*vctr,*vctrMap;
	BooLean	first_only=FALSE;
	FILE	*fp=0;

	if(argc < 2) print_error(USAGE);
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE);
           switch(argv[arg][1]) {
             case 'f':
		if(strcmp("-first",argv[arg]) == 0) first_only=TRUE;
                else print_error(USAGE);
	     case 'x': break;
             default: print_error(USAGE);
           }
	}
	ofp=stderr; fprintf(ofp,"argc = %d\n",argc);
	AB=MkAlphabet(AMINO_ACIDS,GBLAST_BLOSUM62,SREL26_BLSM62);
	infile=AllocString(argv[1]);
	fp=open_file(argv[1],".cnh","r");
	NEWP(Name,105,char);
	Int4	num_files=0;
	while(fscanf(fp,"%s\n",str) == 1){
		num_files++;
		Name[num_files]=AllocString(str);
		fprintf(stderr,"%d. %s\n",num_files,Name[num_files]);
		if(num_files > 100) print_error(USAGE);
	} fclose(fp);
	// 1. Gather information regarding similarities between all trees.
	TotalNodes=0;
	cns_set=0;
	NEWP(tree,num_files+3, wdg_typ);
	NEWP(ident,num_files+3, Int4 *);
	NEWP(NodeSet,num_files+3, set_typ);
	NEWP(MiscSet,num_files+3, set_typ);
	NEW(NumNodes,num_files+3, Int4);
	NEWP(Hpt,num_files+3,hpt_typ);
	for(NumTrees=0,arg=1; arg <= num_files ; arg++){
	   // 1a. For this tree find the corresponding trees 
	   NEW(tree[arg],num_files+3,wdg_typ);
	   NEWP(ident[arg],num_files+3,Int4);
	   Argv[2]=Name[arg];	// this is the main set...
	   NumTrees++;
	   for(i=0,a=1; a <= num_files; a++){
		Argv[1]=Name[a];
		c2h_typ c2h(2,Argv); c2h.NoTree(); c2h.Run(); i++; 
		if(a == arg){
			NodeSet[a]=c2h.RtnSetsI(n);
			MiscSet[a]=c2h.RtnMiscSetsI(NumNodes[a]);
			assert(n == NumNodes[a]);
			Hpt[a]=c2h.CopyOfHptI();
		}

		tree[arg][a]=c2h.RtnTreeI();	// sets TreeI=0; so need to free up.
		vctr=c2h.EqualArray[1];
		TotalNodes += vctr[0];
		for(n=0,k=1; k <= vctr[0]; k++){
			if(vctr[k] > 0){ n++; } fprintf(ofp,"%d",vctr[k]);
		} fprintf(ofp," (%d) %s\n",n,Name[arg]);

		vctr=c2h.EqualArray[0];
		vctrMap=c2h.NodeMap;
		for(n=0,k=1; k <= vctr[0]; k++){
			if(vctr[k] > 0){ n++; } fprintf(ofp,"%d",vctr[k]);
		} fprintf(ofp," (%d) %s\n",n,Name[a]);
		for(k=1; k <= vctrMap[0]; k++){
			fprintf(ofp,"%d ",vctrMap[k]);
		} fprintf(ofp,"\n");

		ident[arg][a]=c2h.RtnMatchesI( );
	   }
	   // if(first_only) break;
	}
#if 0
	if(first_only){
	   i=1;
	   fprintf(ofp,"================ %s ===============\n",Name[i]);
	   Int4 *Nid; NEW(Nid,ident[i][i][0]+3,Int4);
	   // for(j=1; j <= NumTrees; j++)
	   for(j=1; j <= num_files; j++){
		for(n=0,k=1; k <= ident[i][j][0]; k++){
			if(ident[i][j][k] > 0){ n++; Nid[k]++; } fprintf(ofp,"%d",ident[i][j][k]);
		} fprintf(ofp," (%d/%d) %s vs %s\n",n,ident[i][j][0],Name[i],Name[j]);
		// PutWdgraph(ofp,tree[i][j]);
	   } fprintf(ofp,"!");
	   for(k=2; k <= ident[i][i][0]; k++){ 
		if(Nid[k] >= 10) fprintf(ofp,"!"); else fprintf(ofp,"%d",Nid[k]);
	   } fprintf(ofp,"\n\n");
	   free(Nid);
	   exit(1);
	}
#endif
	// 2. Renumber all nodes in order to identify equivalence classes.
	NEW(Node2Tree,TotalNodes +3, Int4);
	NEW(Node2id,TotalNodes +3, Int4);
	NEWP(Node,NumTrees+3,Int4);
	Nnds=1;	// The universal root node == 1.
	for(i=1; i <= NumTrees; i++){
	   fprintf(ofp,"================ %s ===============\n",Name[i]);
	   NEW(Node[i],ident[i][i][0] +3, Int4);
	   Node[i][1]=1; Node2id[1]=1; 
	   Int4 *Nid; NEW(Nid,ident[i][i][0]+3,Int4);
	   for(k=2; k <= ident[i][i][0]; k++){
		Nnds++; Node2Tree[Nnds]=i; Node2id[Nnds] = k;  Node[i][k]=Nnds;
	   }
	   for(j=1; j <= NumTrees; j++){
		for(n=0,k=1; k <= ident[i][j][0]; k++){
			if(ident[i][j][k] > 0){ n++; Nid[k]++; } fprintf(ofp,"%d",ident[i][j][k]);
		} fprintf(ofp," (%d/%d) %s vs %s\n",n,ident[i][j][0],Name[i],Name[j]);
		// PutWdgraph(ofp,tree[i][j]);
	   } fprintf(ofp,"!");
	   for(k=2; k <= ident[i][i][0]; k++){ 
		if(Nid[k] >= 10) fprintf(ofp,"!"); else fprintf(ofp,"%d",Nid[k]);
	   } fprintf(ofp,"\n\n");
	   free(Nid);
	}
	if(first_only) exit(1);
	fp=open_file(argv[1],".mma","r"); MainCMA=ReadCMSA(fp,AB); fclose(fp);
	sets=DSets(Nnds+3);
	grph=MkWdgraph(Nnds+3, Nnds+3);
	Grph=MkWdgraph(Nnds+3, Nnds+3);
}

void	cnh_typ::Free()
{
	Int4	i,j;
	NilWdgraph(grph); NilWdgraph(Grph);
	for(i=1; i <= NumTrees; i++){
	   if(Hpt[i]) delete Hpt[i];
	   for(j=1; j <= NumTrees; j++){
	        if(tree[i][j]) NilWdgraph(tree[i][j]);
		if(ident[i][j]) free(ident[i][j]);
	   } free(ident[i]); free(tree[i]);
	   free(Name[i]); free(Node[i]);
	   for(j=1; j <= NumNodes[i]; j++){
		if(NodeSet[i][j]) NilSet(NodeSet[i][j]); 
		if(MiscSet[i][j]) NilSet(MiscSet[i][j]); 
	   } free(NodeSet[i]); free(MiscSet[i]);
	} free(tree); free(Hpt); free(ident); free(NodeSet);
	if(cns_set){
	   for(i=1; i <= Nnds; i++){
		if(cns_set[i]) NilSet(cns_set[i]);
	   } free(cns_set);
	}
	free(Node2Tree); free(Node2id); free(NumNodes);
	TotalNilCMSA(MainCMA);
	NilDSets(sets);
}

set_typ cnh_typ::GetMemberGraph(wdg_typ G)
{
	Int4	n;
	set_typ	rtn=MakeSet(WdgraphN(G) +9); ClearSet(rtn);
	for(n=1; n <=nWdgraph(G); n++){
		if(FirstInWdgraph(n,G) != 0) AddSet(n,rtn);
	} return rtn;
}

void    cnh_typ::PutConSetOverlap(FILE *fp, set_typ tstSet)
{
	Int4	V,U,i,j,v,n,card,maxI,minI,x;
	grf_typ	*grf = new grf_typ(Nnds+2); // vertex numbering is from 0.
	for(V=1; V <= Nnds; V++){
	   if(tstSet && !MemberSet(V,tstSet)) continue;
	   if(cns_set[V]){
	    set_typ same_set=MakeSet(Nnds+3); ClearSet(same_set); AddSet(V,same_set);
	    set_typ tmpI=CopySet(cns_set[V]),tmpU=CopySet(cns_set[V]); ClearSet(tmpU);
	    maxI=0; minI=NumSeqsCMSA(MainCMA);
	    for(U=1; U <= Nnds; U++){
	   	if(tstSet && !MemberSet(U,tstSet)) continue;
		if(cns_set[U] && V != U){
		  IntersectSet1(cns_set[V],cns_set[U],tmpI); 	// tmpI = V & U.
		  card = CardSet(tmpI);
		  if(card > 0){
		  	if(maxI < card) maxI=card;
		  	if(minI > card) minI=card;
			AddSet(U,same_set); 
#if 0
			Int4 C=CardSet(cns_set[V]),K=CardSet(cns_set[U]); 
			double d=(double)card/(double) MINIMUM(Int4,C,K);
		        if(d >= 0.50) grf->AddEdge(V,U);
#else
		        grf->AddEdge(V,U);
#endif
		  	// IntersectSet3(tmpI,cns_set[U]); // tmpI = tmpI && cns_set[U].
		  	UnionSet(tmpU,tmpI); 	// tmpU = tmpU U tmpI.
		  }
		}
	    } if(minI > maxI) minI=0;
	    i=Node2Tree[V]; v=Node2id[V];
	    set_typ S=0;
	    if(MiscSet[i][v]) S=MiscSet[i][v]; else S=NodeSet[i][v];
	    for(n=0,x=1; x <= Nnds; x++) if(findDSets(x,sets) == V) n++;
	    fprintf(stderr,
		"%d: %d.%d.%s (%d)(consensus of %d = %d; overlap union %d (range %d..%d); # nodes = %d)\n",
			V,i,v,Hpt[i]->SetName(v),CardSet(S),n,CardSet(cns_set[V]),
				CardSet(tmpU),minI,maxI,CardSet(same_set));
	    if(CardSet(same_set) > 1){
		fprintf(stderr,"   "); PutSet(stderr,same_set);
	    } NilSet(tmpU); NilSet(tmpI); 
	    NilSet(same_set);
	   }
	}
	Int4	NumClust,hpsz=50,limit=3;
	double	pcut=0.01;
	unsigned short TrueN=Nnds;
	set_typ	*nd_set=0;	// never used by grf_typ !!
	
	// grf->Put(stderr);
	vst_typ **clique=grf->Bron_Kerbosch_cluster(0,limit,hpsz,&NumClust,pcut,nd_set,TrueN);
	if(clique){
         for(n=1; n <= NumClust; n++){
	  fprintf(stderr,"Clique %d:",n);
	  for(j=0; j < clique[n]->Size(); j++){
		Int4 nd=clique[n]->Vertex(j);
	  	fprintf(stderr," %d",nd);
          } fprintf(stderr,"\n");
         } fprintf(stderr,"\n");
         for(n=1; n <= NumClust; n++){
		delete clique[n];
	 } free(clique);
	} delete grf;
}

wdg_typ	cnh_typ::GetConsTree()
{
	Int4	c,i,j,k,n,u,v,w,U,V,W,e,f,E,N,m,EqNd,x,y,z,UU,VV,X,Y,YY;
	double	d,dd;
	wdg_typ	G,H;
	//==============  1. Partition nodes into disjoint equivalence sets. ===============
	for(i=1; i <= NumTrees; i++){
	   for(j=1; j <= NumTrees; j++){
	      G=tree[i][j]; 
	      for(n=1; n <=nWdgraph(G); n++){
		e=FirstInWdgraph(n,G);
		if(e == 0) continue;	// should happen for root node.
		u=TailWdgraph(e,G); v=HeadWdgraph(e,G); w=WeightWdgraph(e,G);
		assert(n == v);	//  u -> v == w.
		U=Node[i][u]; V=Node[i][v]; W=Node[j][w];	// alternative numbering.
		if(W > 0){	// --> node i.v is equilvalent to node j.w
			if(i != j && NodeSet[i][v] && NodeSet[j][w]){
			   x=CardSet(NodeSet[i][v]); y=CardSet(NodeSet[j][w]);
			   z=CardInterSet(NodeSet[i][v],NodeSet[j][w]);
			   d=100*(double)z/(double)MINIMUM(Int4,x,y);
			   if(d < 98.0) fprintf(stderr,
				"%d:%d(%d) & %d:%d(%d) = %d (%.1f)\n",i,v,x,j,w,y,z,d);
			   if(MiscSet[i][v] && MiscSet[j][w]){
			      x=CardSet(MiscSet[i][v]); y=CardSet(MiscSet[j][w]);
			      z= CardInterSet(MiscSet[i][v],MiscSet[j][w]);
			      d=100*(double)z/(double)MINIMUM(Int4,x,y);
			      if(d < 98.0) fprintf(stderr,
				" <subtree> %d:%d(%d) & %d:%d(%d) = %d (%.1f)\n",i,v,x,j,w,y,z,d);
			   }
			}
			v=findDSets(V,sets); u = findDSets(W,sets);
			linkDSets(v,u,sets);
		}
	      }
	   } Hpt[i]->PutClean(stderr);
	}
#if 1
	//==============  1. Put other nodes into disjoint equivalence sets. ===============
	set_typ	Si,Sj;
	for(i=1; i < NumTrees; i++){
	   G=tree[i][i]; 
	   for(j=i+1; j <= NumTrees; j++){
	     H=tree[j][j]; 
	     for(n=1; n <=nWdgraph(G); n++){
		e=FirstInWdgraph(n,G);
		if(e == 0) continue;	// should happen for root node.
		u=TailWdgraph(e,G); v=HeadWdgraph(e,G); 
		assert(n == v);	//  u -> v == w.
		if(MiscSet[i][v]) Si=MiscSet[i][v]; else Si=NodeSet[i][v];
		U=Node[i][u]; V=Node[i][v];  // alternative numbering.
		VV=findDSets(V,sets); 
		for(m=1; m <=nWdgraph(H); m++){
		   f=FirstInWdgraph(m,H);
		   if(f == 0) continue;	
		   x=TailWdgraph(f,H); y=HeadWdgraph(f,H); 
		   assert(m == y);	//  u -> v == w.
		   X=Node[j][x]; Y=Node[j][y];  // alternative numbering.
		   YY=findDSets(Y,sets); 
		   if(VV == YY) continue;
		   if(MiscSet[j][y]) Sj=MiscSet[j][y]; else Sj=NodeSet[j][y];
		   c=CardSet(Si); k=CardSet(Sj);
		   z=CardInterSet(Si,Sj);
		   d=100*(double)z/(double)MINIMUM(Int4,c,k);
		   dd=100*(double)z/(double)MAXIMUM(Int4,c,k);
		   if(d >= 90.0 && dd >= 90.0){
			fprintf(stderr,
			   " <subtree> %d:%d(%d) & %d:%d(%d) = %d (%.1f)\n",i,v,c,j,y,k,z,d);
			linkDSets(findDSets(VV,sets),findDSets(YY,sets),sets);
		   }
	        }
	     }
	   }
	}
#endif
	//================== 2. Get size of canonical nodes. =======================
	NEW(cns_set,Nnds+3,set_typ);
	for(i=1; i <= NumTrees; i++){
	   for(j=1; j <= NumTrees; j++){
	     G=tree[i][j]; 
	     for(n=1; n <=nWdgraph(G); n++){
		e=FirstInWdgraph(n,G);
		if(e == 0) continue;	// should happen for root node.
		V=Node[i][n]; VV=findDSets(V,sets); 
		if(cns_set[VV] == 0) cns_set[VV]=CopySet(NodeSet[i][n]);
		else IntersectSet3(cns_set[VV],NodeSet[i][n]); // cns_set[VV]=cns_set[VV] && NodeSet[i][n]
	     }
	  }
	}
#if 1	// debug...
	PutConSetOverlap(stderr);
#endif

	//=================  3. Create a univeral graph. ========================
	for(i=1; i <= NumTrees; i++){
	   G=tree[i][i]; 
	   for(n=1; n <=nWdgraph(G); n++){
	      e=FirstInWdgraph(n,G);
	      if(e == 0) continue;	// should happen for root only.
	      u=TailWdgraph(e,G); v=HeadWdgraph(e,G); w=WeightWdgraph(e,G);
	      U=Node[i][u]; V=Node[i][v];  // alternative (global) numbering.
	      if(u == 1){
		assert(U == 1);
	        JoinWdgraph(1,V,1,grph); // add these to universal graph.
		if(V == 1){ 
			fprintf(stderr,"v = %d; u = %d; V= %d; U = %d\n",v,u,V,U);
			PutWdgraph(stderr,G); exit(1); 
		}
		  // do nothing....
	      } else {
	        JoinWdgraph(U,V,1,grph); // add these to universal graph.
	      }
	   }
	}
	fprintf(stderr,"=========== grph ============\n");
	PutWdgraph(stderr,grph);

	//================== 4. Create a consensus Graph. =======================
	for(n=1; n <= nWdgraph(grph); n++){
		e=FirstInWdgraph(n,grph);
		if(e == 0) continue;	// should happen for root only.
		u=TailWdgraph(e,grph); v=HeadWdgraph(e,grph); w=WeightWdgraph(e,grph);
		U=findDSets(u,sets); V=findDSets(v,sets); 
		if(V == v){
		    E = JoinWdgraph(U,V,1,Grph); // add these to consensus digraph.
		}
	}
	// fprintf(stderr,"=========== Grph ============\n");
	// PutWdgraph(stderr,Grph);

	//==================== 5. Add weights to consensus Graph. =====================
	for(n=1; n <= nWdgraph(grph); n++){
		e=FirstInWdgraph(n,grph);
		if(e == 0) continue;	// should happen for root only.
		u=TailWdgraph(e,grph); v=HeadWdgraph(e,grph); 
		U=findDSets(u,sets); V=findDSets(v,sets); 
		if(V != v){
		    for(E=FirstInWdgraph(V,Grph); E!=0; E=NextInWdgraph(E,Grph)){
			u = TailWdgraph(E,Grph);	// consensus parent node.
			if(U == u){
			   w = WeightWdgraph(E,Grph);
			   SetWeightWdgraph(w+1,E,Grph);
			}
		    }
		}
	}
	fprintf(stderr,"=========== Grph ============\n");
	PutWdgraph(stderr,Grph);

    	wdg_typ Tree=RtnConsensusTree(1,NumTrees,Grph);
	fprintf(stderr,"=========== Consensus Tree ============\n");
	PutWdgraph(stderr,Tree);
	PrintNewickTreeWDG(stderr,1,Tree);

	Int4	NumRm;
	wdg_typ TREE=0;
	do { TREE=PurgeTree(1,NumTrees,NumRm,Tree); NilWdgraph(Tree); Tree=TREE; } while(NumRm > 0);

	fprintf(stderr,"=========== Purged Consensus Tree ============\n");
	PutWdgraph(stderr,Tree);
	PrintNewickTreeWDG(stderr,1,Tree);

	set_typ	tmp=GetMemberGraph(Tree); PutConSetOverlap(stderr,tmp); NilSet(tmp);

	//============== 6. Renumber and rename the tree and the corresponding sets. ===============

	//============== 7. Get the hpt for the tree. ===============
#if 1
        int argcnt=0;
        char *argval[9];
        argval[0]=AllocString("tree2hpt"); argcnt++;
        argval[1]=AllocString("tmpfile"); argcnt++;	// this is not used..
        argval[2]=AllocString("-seed=0"); argcnt++; // don't reseed random generator.
        FILE *hpt_fp=tmpfile();
        FILE *ifp=tmpfile();
	PrintNewickTreeWDG(ifp,1,Tree); rewind(ifp);
        // 4. Create a new Hyperpartition based on processing analysis.
        tree2hpt(hpt_fp,ifp,argcnt,argval,Hpt[1]->NumberRandom());
	fclose(ifp); rewind(hpt_fp);
	hpt_typ *tHpt= new hpt_typ(hpt_fp); fclose(hpt_fp); 
	for(Int4 x=0; x < argcnt; x++) free(argval[x]);
	tHpt->Put(stderr);
#endif

	//============== 8. Get the seed cma files (*.sma) for the hpt. ===============
	//============== 9. output files (*.hpt && *.sma) for the hpt. ===============
	// need to create misc_set csq for internal nodes.
	// need to read in the main input alignment.
	// set_typ *st_set; NEW(st_set,tHpt->NumSets()+2,set_typ);
	// cma_typ *sma,cma; NEW(sma,tHpt->NumSets()+2,cma_typ);
	FILE *fp=open_file(infile,"_out.sma","w");
	RenameCMSA("Set1",MainCMA); PutConsensusCMSA(fp,MainCMA);
	///// sma[1]=MakeConsensusCMSA(MainCMA); 
	for(i=2; i < tHpt->NumSets(); i++){
		m=tHpt->ItoSetID(i); assert(m > 0 && m <= Nnds);
		set_typ Tmp=tHpt->MkSubTreeSet(i);
		set_typ TmpSq=CopySet(cns_set[m]); sprintf(str,"Set%d",m);
		for(n=0,j=1; j < tHpt->NumSets(); j++){
		   if(j == i) continue;
		   if(MemberSet(j,Tmp)){
			n++; m=tHpt->ItoSetID(j); assert(m <= Nnds);
		  	UnionSet(TmpSq,cns_set[m]); 	// TmpSq = TmpSq U cns_set[m].
		   }
		}
		FILE *tfp=tmpfile();
		// PutMergedCMSA(tfp,n,st_set,MainCMA,0);
		PutInSetCMSA(tfp,TmpSq,MainCMA); rewind(tfp); 
		cma_typ cma = ReadCMSA(tfp,AB); fclose(tfp);
		RenameCMSA(str,cma);
		PutConsensusCMSA(fp,cma);
		// sma[i]=MakeConsensusCMSA(cma); 
		// RenameCMSA(str,sma[i]);
		TotalNilCMSA(cma);
		NilSet(Tmp); NilSet(TmpSq);
	} fclose(fp);
	fp=open_file(infile,"_out.hpt","w"); tHpt->Put(fp); fclose(fp);
	
		
	// PutDSets(stderr,sets);

	return Tree; 
}


