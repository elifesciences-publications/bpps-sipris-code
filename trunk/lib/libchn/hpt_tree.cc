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

#include "hpt_typ.h"

BooLean	hpt_typ::IsTree(Int4 *&Parent)
// Does the Hpt correspond to a (properly formated) tree?
{
	Int4	r,c,R,C,P,G,NRows=NumElementarySets, NCols=NumberBPPS;
	Int4	FG,BG,OG; 	// 
	char	**cell=HyperPartition; // cell[Row][Column];
	Parent=0;

	if(NRows != NumberBPPS+1) return FALSE;	// 1 category (column) for each subtree (row).
	for(R=1; R < NRows; R++){
		if(cell[R][1] != '+') return FALSE;	// first column (root) all cells == '+' .
		if(cell[R][R] != '+') return FALSE;	// diagonal all '+'.
		for(c=R+1; c <= NCols; c++) if(cell[R][c] == '+') return FALSE; // upper right half != '+'.
	} if(cell[NRows][1] != '-') return FALSE;
	NEW(Parent,NRows+3,Int4); // free up from calling environment.
	for(C=2; C <= NCols; C++){	// C = child & column.
#if 0
		for(FG=BG=OG=0,r=1; r < NRows; r++){	// count the # FG & BG nodes in column R.
		   switch(cell[r][C]){
			case '+': FG++; break;
			case '-': BG++; break;
			case 'o': OG++; break;
			default: print_error("IsTree( ) input error");
		   }
		}
#endif
		if(cell[NRows][C] != 'o') return FALSE;	// Reject row must be 'o'.
		for(P=C-1; P > 0; P--) if(cell[C][P]== '+') break;  assert(P > 0); // found Parent.

		if(P == 1){		// Parent == root node. 
		   for(r=2; r < NRows; r++){	// row 1 can be either '-' or 'o' in this case.
			if(r == C) continue;	// already checked above.
			else if(cell[r][C] == 'o') return FALSE;  // with Root as parent, 'o' is disallowed.
			// can be either '+' or '-' at these positions.
		   }

		} else {	// an internal node is the parent; find the grandparent.
		   if(cell[1][C] != 'o') return FALSE;	// child must be 'o' in root row for this case.
		   // for(G=P-1; G > 0; G--) if(cell[C][G]== '+') break; assert(G > 0); // found Grandparent.
		   for(r=2; r < NRows; r++){	// row 1 == '-' or 'o'.
			if(r==C) continue;	// checked above.
// fprintf(stderr,"r=%d; P = %d; C= %d\n",r,P,C);
			if(cell[r][P] == 'o' && cell[r][C] != 'o') return FALSE; // both must be 'o'.
			if(cell[r][P] == '-' && cell[r][C] != 'o') return FALSE; // latter must be 'o'.
			if(cell[r][P] == '+' && cell[r][C] == 'o') return FALSE; // latter must be '+' or '-'.
			if(cell[r][P] != '+' && cell[r][C] == '+') return FALSE; // former must be '+'.
		   }
		} Parent[C]=P;	// use this to create the tree for Newick printout.
	} return TRUE;
}

BooLean hpt_typ::IsTree(Int4 *&Parent, Int4 *&children)
{
        BooLean rtn=this->IsTree(Parent);
        if(rtn){
           NEW(children,this->NumBPPS()+3,Int4);
           for(Int4 i=1; i<= this->NumBPPS(); i++){
              for(Int4 j=1; j <= this->NumBPPS(); j++){ if(Parent[j]==i) children[i]++; }
           }
        } else { children=0; }
	return rtn;
}

set_typ hpt_typ::DescendantSet(Int4 n)
// Return the set of all descendants of node n.

{
	set_typ set=0;
	Int4	i,j,k,*Parent=0;
        BooLean rtn=this->IsTree(Parent);
        if(rtn){
	   set=MakeSet(this->NumSets() +3);
	   do {		// for all descendants 'j' of 'n'; set Parent[j] = n.
             for(k=0,i=1; i<= this->NumSets(); i++){
	       if(Parent[i] == n){
                for(j=1; j<= this->NumSets(); j++){ if(Parent[j]==i){ Parent[j]=n; k++; } }
	       }
	     }
	   } while(k > 0);
           for(i=1; i<= this->NumSets(); i++){
		if(Parent[i] == n) AddSet(i,set); 
	   }
	} if(Parent) free(Parent);
	return set;
}

BooLean	hpt_typ::PutDFTree(FILE *fp)
// Warning: This assumes that the output file is for a GISMO hierarchical alignment MAPGAPS file.
// Thus it includes the main set.
{
	Int4	i,*rtn;
	if(this->RtnDFTree(rtn) == FALSE){ free(rtn); return FALSE; }
	fprintf(fp,"0");
	for(i=1; i < this->NumSets(); i++){ fprintf(fp,",%d",rtn[i]+1); } fprintf(fp,".\n");
	free(rtn);
	return TRUE;
}

BooLean	hpt_typ::RtnDFTree(Int4 *&rtn)
#if 0
	HyperParTition:                 output ==> [0,1,2,2,2,1,2,1,0]
	+-ooo-o- 1.Set1?		output ==> "0,1,2,3,3,3,2,3,2." from PutDFTree!!
	++----o- 2.Set8?
	+++---o- 3.Set7!
	++-+--o- 4.Set5!
	++--+-o- 5.Set6!
	+-ooo+-- 6.Set2?
	+-ooo++- 7.Set4!
	+-ooo-o+ 8.Set3!  		(used by mapgaps program)
#endif
{
	Int4	c,p,*Parent;
	NEW(rtn,this->NumSets() + 5, Int4);
	if(this->IsTree(Parent)){
		wdg_typ Tree=MkWdgraph(NumberBPPS+9,NumberBPPS+9);  // MkWdgraph(nodes,edges);
		for(c = NumberBPPS; c > 1; c--) {
		   p=Parent[c]; assert(p > 0);
		   JoinWdgraph(p,c,1,Tree);  // p ---> c.
		} free(Parent); 
        	set_typ S=MakeSet(WdgraphN(Tree)+1); ClearSet(S);
        	DFTreeSrch(rtn, 0, 1, S, Tree); // fprintf(fp,".\n");
		NilWdgraph(Tree); NilSet(S);
		return TRUE;
	} else {
		return FALSE; // eventually just add '0's to the array.
	}
}

// void    hpt_typ::DFTreeSrch(FILE *fp, Int4 depth, Int4 v, set_typ S, wdg_typ T)
void    hpt_typ::DFTreeSrch(Int4 *rtn, Int4 depth, Int4 v, set_typ S, wdg_typ T)
{
        // start from root, recursively print tree in Newick format.
        Int4    w,e,i,n=NumEdgesOut(v,T);
        if(MemberSet(v,S)) print_error("hpt_typ::DFTreeSrch( ): input graph not a tree (cycle found).");

#if 0
	if(depth > 0) fprintf(fp,","); 
	fprintf(fp,"%d",depth);
	// fprintf(fp,"%s\n",NameElementarySet[v]);
#else
        AddSet(v,S); i=CardSet(S); rtn[i]=depth;
#endif
        if(n > 0){
           for(i=0,e=FirstOutWdgraph(v,T);e!=0;e=NextOutWdgraph(e,T)){
                i++; w = HeadWdgraph(e,T);      // a child of v.
                DFTreeSrch(rtn,depth+1,w,S,T);
           } 
        } 
}

/************************ PrintNewickTree ************************/
void    hpt_typ::PrintSmartArtTree(FILE *fp,Int4 root, wdg_typ T)
{ // print tree in Newick format.
        // fprintf(fp,"("); TreeDFS(fp, root, T); fprintf(fp,")%d;\n",root);
        // PutWdgraph(stderr,T);
        set_typ S=MakeSet(WdgraphN(T)+1); ClearSet(S);
        SmartArtDFS(fp, 0, root, S, T); fprintf(fp,"\n");
        NilSet(S);
}

void    hpt_typ::SmartArtDFS(FILE *fp, Int4 depth, Int4 v, set_typ S, wdg_typ T)
{
        // start from root, recursively print tree in Newick format.
        Int4    w,e,i,n=NumEdgesOut(v,T);
        if(MemberSet(v,S)) print_error("hpt_typ::SmartArtDFS( ): input graph not a tree (cycle found).");
	for(i=1; i <= depth; i++) fprintf(fp,"\t"); fprintf(fp,"%s\n",NameElementarySet[v]);
        AddSet(v,S);
        if(n > 0){
           // fprintf(fp,"(");
           for(i=0,e=FirstOutWdgraph(v,T);e!=0;e=NextOutWdgraph(e,T)){
                i++; w = HeadWdgraph(e,T);      // a child of v.
                SmartArtDFS(fp,depth+1,w,S,T);
                // if(i < n) fprintf(fp,",");
           } // fprintf(fp,")");
        } 
}

BooLean	hpt_typ::PutAsSmartArt(FILE *fp)
{
	Int4	c,p,*Parent;
	if(this->IsTree(Parent)){
		wdg_typ Tree=MkWdgraph(NumberBPPS+9,NumberBPPS+9);  // MkWdgraph(nodes,edges);
		// for(c=2; c <= NumberBPPS; c++)
		for(c = NumberBPPS; c > 1; c--) {
		   p=Parent[c]; assert(p > 0);
		   JoinWdgraph(p,c,1,Tree);  // p ---> c.
		} PrintSmartArtTree(fp,1,Tree); NilWdgraph(Tree);
		free(Parent); return TRUE;
	} else print_error("hpt_typ::PutAsSmartArt() requires Tree formatting");
}

wdg_typ	hpt_typ::RtnAsTree( )
{
	Int4	c,p,*Parent;
	if(this->IsTree(Parent)){
		wdg_typ Tree=MakeWdgraph(NumberBPPS+9,NumberBPPS+9);  // MkWdgraph(nodes,edges);
		AddVerticesWdgraph(NumberBPPS,Tree);
		for(c=2; c <= NumberBPPS; c++){
		   p=Parent[c]; assert(p > 0);
		   JoinWdgraph(p,c,0,Tree);  // p ---> c.
		} return Tree; 
	} else print_error("hpt_typ::RtnAsTree() requires Tree formatting");
	return 0;
}

/************************ PrintNewickTree ************************/

BooLean	hpt_typ::PutAsTree(FILE *fp)
{
	Int4	c,p,*Parent;
	if(this->IsTree(Parent)){
		wdg_typ Tree=MkWdgraph(NumberBPPS+9,NumberBPPS+9);  // MkWdgraph(nodes,edges);
		// for(c=2; c <= NumberBPPS; c++)
		for(c = NumberBPPS; c > 1; c--) {
		   p=Parent[c]; assert(p > 0);
		   JoinWdgraph(p,c,1,Tree);  // p ---> c.
		} PrintNewickTree(fp,1,Tree); NilWdgraph(Tree);
		free(Parent); return TRUE;
	} else {
	   assert(Parent == 0);
	   if(this->IsScrambledTree(Parent)){
		wdg_typ Tree=MkWdgraph(NumberBPPS+9,NumberBPPS+9);  // MkWdgraph(nodes,edges);
                for(c = NumberBPPS; c > 1; c--) {
                   p=Parent[c]; assert(p > 0);
                   JoinWdgraph(p,c,1,Tree);  // p ---> c.
                } PrintNewickTree(fp,1,Tree); NilWdgraph(Tree);
		free(Parent); return TRUE;
	   } else if(Parent!=0){
		// for(c=2; c <= NumberBPPS; c++)
		for(c = NumberBPPS; c > 1; c--){
		  if(Parent[c]){
			fprintf(stderr,"parent %d --> child %d; error on %d\n",
				Parent[c],c,c+1); break; 
		  }
		} free(Parent);
		fprintf(stderr,"Hpt does not correspond to a tree\n");
		return FALSE;
	  } else print_error("PutAsTree(): this should not happen");
	}
}

void    hpt_typ::TreeDFS(FILE *fp, Int4 v, wdg_typ T)
// start from root, recursively print tree in Newick format.
{
      Int4    w,e,i,n=NumEdgesOut(v,T);
      if(n > 0){
	   fprintf(fp,"(");
	   for(i=0,e=FirstOutWdgraph(v,T);e!=0;e=NextOutWdgraph(e,T)){
		i++; w = HeadWdgraph(e,T);      // a child of v.
		TreeDFS(fp,w,T);
		if(i < n) fprintf(fp,",");
	   } fprintf(fp,")");
	  } 
	  assert(v > 0 && v <= NumElementarySets);
	  fprintf(fp,"%s",NameElementarySet[v]);
	  // fprintf(fp,"%d_%s",v, NameElementarySet[v]);
	  // fprintf(fp,"%d_Set%d",v,v);      // output for HyperPartition.
}

hpt_typ *hpt_typ::MkLineageTree(Int4 node, Int4 &new_node)
// collapse hpt_typ into a lineage to a given node's subtree.
{
        hpt_typ *hpt,*Hpt=this->Copy();
        Int4    n,m,*P,N=0;
	if(!Hpt->IsTree(P)) print_error("hpt_typ::MkLineageTree() input error");
	if(node < 2 || node >= this->NumSets()){
		print_error("hpt_typ::MkLineageTree(): input node out of range");
	}
	char *name,*Name=AllocString(Hpt->SetName(node));
	do {
           for(n=1; n < Hpt->NumSets(); n++){ if(strcmp(Name,Hpt->SetName(n))==0) break; }
	   assert(n < Hpt->NumSets()); node=n; N=Hpt->NumSets();
	   //fprintf(stderr,"node=%d; Name=%s\n",n,Name);
           set_typ line=Hpt->MkLineageSet(node);
           set_typ subtree=Hpt->MkSubTreeSet(node);
           set_typ leaves=Hpt->MkLeafSet();
           for(n=1; n < Hpt->NumSets(); n++){
              if(!MemberSet(n,subtree) && !MemberSet(n,line) && MemberSet(n,leaves)){
		hpt=Hpt->Delete(n); delete Hpt; Hpt=hpt; break; // need to redo the rest...
              }
           } NilSet(line); NilSet(subtree); NilSet(leaves);
	   // Hpt->Put(stderr,FALSE,FALSE,FALSE);
	} while(N > Hpt->NumSets()); free(Name); new_node=node;
        for(n=1; n < Hpt->NumSets(); n++){
	   name=Hpt->SetName(n);
           for(m=1; m < this->NumSets(); m++){
	      Name=this->SetName(m);
	      if(strcmp(Name,name)==0){
		   char *Arg=this->ArgmntStr[m];
		   // fprintf(stderr,"Arg[%d]=%s\n",m,Arg);
		   Hpt->SetArgStr(n,Arg); 
		   Hpt->ReSetArgv(n,this->nArg(m),this->Argv(m));
		   Arg=Hpt->ArgmntStr[n];
		   // fprintf(stderr,"***Arg[%d]=%s\n",n,Arg);
		   break;
	      }
	   }
	} free(P);
        // Hpt->Put(stderr,TRUE,FALSE,FALSE);
	return Hpt;
}

