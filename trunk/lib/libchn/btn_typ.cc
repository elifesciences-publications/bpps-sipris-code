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
#include "btn_typ.h"

btn_typ::btn_typ(Int4 id, char *Name )
{
	leftchild=0;
	rightchild=0;
	parent=0;
	if(Name) name = AllocString(Name);
	else {
		char str[50]; 
		sprintf(str,"Set%d",id);
		name = AllocString(str);
	}
	ID=id;
}

btn_typ::~btn_typ()
{
	if(leftchild != 0) delete leftchild; // free up all child nodes in tree.
	free(name);
}

btn_typ	*btn_typ::CopyTree()
{
	btn_typ *X=new btn_typ(this->ID,this->name);
	btn_typ *R=rightchild,*L=leftchild,*P=parent;
	if(R){ X->rightchild=R->CopyTree(); X->rightchild->parent = X; }
	if(L){ X->leftchild=L->CopyTree(); X->leftchild->parent = X; }
	return X;
}

Int4	btn_typ::MaxID()
// see what the maximum in subtree is.
{
	Int4	x,y,z;
	if(leftchild && rightchild){	// is an internal node
		x = leftchild->MaxID();
		y = rightchild->MaxID();
		z = MAXIMUM(Int4,x,y);
		return MAXIMUM(Int4,z,ID);
	} else return ID;
}

BooLean	btn_typ::IsOkay()
{
	if((leftchild && !rightchild) || (!leftchild && rightchild)){
		return FALSE;
	} else return TRUE;
}

char	**btn_typ::GetHPT(Int4 &N)
// this is assumed to be the root node...
{
	N=CountNodes(0);
	char **hpt;
	NEWP(hpt,N+3,char);
	for(Int4 i=0; i <= N; i++) NEW(hpt[i],N+3,char);
	return hpt;
}

hpt_typ	*btn_typ::ReturnHPT(Int4 depth, Int4 NumRandom)
{
	FILE *fp=tmpfile(); printHPT(fp,depth,NumRandom); rewind(fp);
	hpt_typ *hpt = new hpt_typ(fp); fclose(fp); return hpt;
}

Int4	btn_typ::printHPT(FILE *fp, Int4 depth, Int4 NumRandom)
/************************************************************
In the hyperpartition, each column corresponds to one node in the tree.

The ‘+’ rows in that column correspond to that node’s subtree,
which serves as the foreground.

The ‘-‘ rows in that column correspond to the rest of the parent node’s subtree,
which serves as the background

The remaining (non-participating) nodes in that column are labeled with an ‘o’.

(For the root node a set of random sequences serves as the background.) 
*************************************************************/
{
	// 1. Check to make sure that each node has a unique name.

	Int4 N=CountNodes(0),i,j,rtnN;

	// 2. Renumber nodes phylogenetically (using dfs) so that the hpt is easy to interpret.
	btn_typ **nodes = OrderAsDFS(rtnN,depth);
	assert(depth > 0 || rtnN == N);
	N=rtnN;
	if(0) fprintf(stderr,"\n");
	for(Int4 i=1; i <= N; i++){
		if(0){ fprintf(stderr,"%d: ",i); nodes[i]->Put(stderr); }
		nodes[i]->ID=i; 
	} if(0) fprintf(stderr,"\n");

	// 3. Allocate a matrix for the hyperpartition:
	char **hpt;
	NEWP(hpt,N+3,char);
	for(i=0; i <= N; i++) NEW(hpt[i],N+3,char);

	// 4. For each node (column in hpt) starting from the root...
	for(i=1; i <= N; i++){
         if(nodes[i]->IsRoot()) for(j=1; j <= N; j++) hpt[i][j]='+'; 
	 else {
   	   // 4.1. Set all nodes on the ith list to 'o'.		// non-participating subgroups.
	   for(j=1; j <= N; j++) hpt[i][j] = 'o';
   	   // 4.2. Reset all nodes in the parent's subtree to '-'.   // background nodes. (Caution: convert from binary to n-ary.
	   btn_typ *pnode=nodes[i]->FindNaryParent( );
	   if(pnode != 0){  		// don't do this for the root which has no parent;
	      assert(!pnode->IsLeafNode( ));	// can't be a leaf if is a parent node...
	      // Int4 *list=pnode->leftchild->GetSubTreeIDs(); // 
	      Int4 *list=pnode->GetSubTreeIDs(); // 
	      for(j=1; j <= N; j++){
		if(list[j] == 0) break;
		assert(list[j] > 0 && list[j] <= N);
		hpt[i][list[j]] = '-';
	      } free(list);
	   }

   	   // 4.3. Reset all nodes in the current node's subtree to '+'.    // foreground nodes.
	   assert(nodes[i]->ID > 0 && nodes[i]->ID <= N);
	   if(nodes[i]->IsLeafNode( )){ hpt[i][nodes[i]->ID]='+'; }
	   else {
	     Int4 *list=nodes[i]->GetSubTreeIDs(); // 
	     for(j=1; j <= N; j++){
		if(list[j] == 0) break;
		assert(list[j] > 0 && list[j] <= N);
		hpt[i][list[j]] = '+';
	     } free(list);
	   } 
	 }
	}

	// 5. Print the hpt.
	// fprintf(fp,"\n");
	fprintf(fp,"HyperParTition:\n");
	for(i=1; i <= N; i++) fprintf(fp,"!"); fprintf(fp,"\n");
	for(i=1; i <= N; i++){
	      for(j=1; j <= N; j++){ 
		 // fprintf(fp," %c ",hpt[j][i]);
		 fprintf(fp,"%c",hpt[j][i]);
	      } 
	      if(nodes[i]->IsLeafNode()) fprintf(fp," %d.%s!\n",i,nodes[i]->name);
	      else fprintf(fp," %d.%s?\n",i,nodes[i]->name);
	} fprintf(fp,"-");
	// fprintf(fp," - ");
	// for(j=2; j <= N; j++) fprintf(fp," o ");
	for(j=2; j <= N; j++) fprintf(fp,"o");
	fprintf(fp," %d.Random=%d.\n",N+1,NumRandom);
	// 6. Free the hpt.
	for(i=0; i <= N; i++) free(hpt[i]); free(hpt);
	free(nodes);
	return N;
}

Int4    *btn_typ::GetSubTreeIDs()
{
	Int4 N=CountNodes(0);
	Int4 *list; NEW(list, N+3,Int4);
	if(IsLeafNode( )){ list[1]=ID; }
	else {
	   Int4 Index=1; list[1]=ID; 
	   // Index++; list[Index]=leftchild->ID;
	   leftchild->SubTreeIDs(Index,N,list); 
	} return list;
}

Int4    btn_typ::SubTreeIDs(Int4 &Index, Int4 N,Int4 *list)
{
	Index++; 
	if(Index > N) fprintf(stderr,"Index = %d\n",Index);
	assert(Index <= N); list[Index] = ID;
	if(leftchild != 0){ leftchild->SubTreeIDs(Index,N,list); }
	if(rightchild != 0){ rightchild->SubTreeIDs(Index,N,list); }
	return Index;
}

btn_typ	*btn_typ::LastSibling()
// Get the last sibling to the current node.
{
	if(rightchild){
		if(rightchild->rightchild == 0) return rightchild;
		else return rightchild->LastSibling();
	} else return 0; // this == end of Sibling list.
}

BooLean	btn_typ::AddParent(btn_typ *Node)
// Add a parent node (Node) to 'this'.
{
	btn_typ *R,*L,*P,*tmp;
	P=parent; L = leftchild; R = rightchild;
	if(P==0){	// this is the root node.
		print_error("AddParent(): Root as target node disallowed");
	}
	if(this == P->leftchild){	// target is head of sibling list.
		P->leftchild=Node; Node->parent=P; 
		Node->leftchild=this; this->parent=Node; 
		Node->rightchild=R; // now siblings of new node.
		if(R) R->parent=Node; 
		this->rightchild=0; 
		// this->leftchild is unchanged.
	} else if(this == P->rightchild){	// target 
		P->rightchild=Node; Node->parent=P; 
		Node->leftchild=this; this->parent=Node; 
		if(R) R->parent=Node; Node->rightchild=R; 
		this->rightchild=0; 
		// this->leftchild is unchanged.
	} else {
		L = leftchild; R = rightchild;
		Node->Put(stderr);
		this->Put(stderr);
		// fprintf(stderr,"parent: "); P->Put(stderr); 
		// if(L){ fprintf(stderr,"leftchild: "); L->Put(stderr); } 
		// if(R){ fprintf(stderr,"rightchild: "); R->Put(stderr); }
		print_error("AddParent(): this should not happen");
	}  return TRUE;
}

btn_typ	*btn_typ::SplitOffParent(btn_typ *Node, Int4 nodeID)
// split off a parent node.
{
	btn_typ *R,*L,*P,*tmp;
	if(ID == nodeID){	// found target node.
		P=parent; L = leftchild; R = rightchild;
		if(P==0){	// this is the root node.
			print_error("SplitOffParent(): Root as target node disallowed");
			Node->parent = 0; Node->leftchild=this;
			Node->rightchild=0; this->parent=Node;
			return this;
		}
		if(this == P->leftchild){	// target is head of sibling list.
			P->leftchild=Node; Node->parent=P; 
			Node->leftchild=this; this->parent=Node; 
			Node->rightchild=R; // now siblings of new node.
			if(R) R->parent=Node; 
			this->rightchild=0; 
			// this->leftchild is unchanged.
		} else if(this == P->rightchild){	// target 
			P->rightchild=Node; Node->parent=P; 
			Node->leftchild=this; this->parent=Node; 
			if(R) R->parent=Node; Node->rightchild=R; 
			this->rightchild=0; 
			// this->leftchild is unchanged.
		} else {
			L = leftchild; R = rightchild;
			Node->Put(stderr);
			this->Put(stderr);
			// fprintf(stderr,"parent: "); P->Put(stderr); 
			// if(L){ fprintf(stderr,"leftchild: "); L->Put(stderr); } 
			// if(R){ fprintf(stderr,"rightchild: "); R->Put(stderr); }
			print_error("SplitOffParent(): this should not happen");
		}
		return this;
	} else {
		R=rightchild; L=leftchild;
		if(R){ tmp=R->SplitOffParent(Node,nodeID); if(tmp) return tmp; }
		if(L){ tmp=L->SplitOffParent(Node,nodeID); if(tmp) return tmp; }
		return 0;
	}
}

BooLean	btn_typ::AddLeaf(btn_typ *Node)
{
	btn_typ *R=this->rightchild,*L=this->leftchild,*tmp;
	if(L){		// this node already has some child nodes.
	   Node->rightchild = L; L->parent=Node;  // 
	   Node->leftchild=0;
	   Node->parent=this; this->leftchild=Node;
	} else {		// no subtree.
	   Node->parent=this; this->leftchild=Node;
	   Node->leftchild=0; Node->rightchild=0; 
	} return TRUE;
}

btn_typ	*btn_typ::SplitOffLeaf(btn_typ *Node, Int4 nodeID)
// add a child == Node; remember this is a binary tree!! So rightchild == siblings.
{
	btn_typ *R=this->rightchild,*L=this->leftchild,*tmp;
	if(this->ID == nodeID){	// found parent node to which node should be attached.
	    if(L){		// this node already has some child nodes.
		Node->rightchild = L; L->parent=Node;  // 
		Node->leftchild=0;
		Node->parent=this; this->leftchild=Node;
	    } else {		// no subtree.
		Node->parent=this; this->leftchild=Node;
		Node->leftchild=0; Node->rightchild=0; 
	    } return this;
	} else {
		if(R){ tmp=R->SplitOffLeaf(Node,nodeID); if(tmp) return tmp; }
		if(L){ tmp=L->SplitOffLeaf(Node,nodeID); if(tmp) return tmp; }
		return 0;
	}
}

BooLean	btn_typ::RemoveNode( )
// remove 'this' node and link its childrent to parent node.
{
	btn_typ *R=rightchild,*L=leftchild,*P=parent,*tmp;
	if(P){
	   if(this == P->rightchild){	// 'this' is within Sibling list.
		if(L){		// 'this' has a subtree attached.
			P->rightchild=L; L->parent = P;
			tmp=L->LastSibling(); 
			if(tmp){	// other siblings of L.
				tmp->rightchild=R; 
				if(R) R->parent=tmp;
			} else {	// L has no siblings.
				L->rightchild=R;
				if(R) R->parent=L;
			}
		} else {		// no subtree for 'this'.
			P->rightchild = R;  // simply delete 'this'
			if(R) R->parent=P;
		}
	   } else if(this == P->leftchild){  // 'this' = start of Sibling list.
		if(L){		// 'this' has a subtree.
			P->leftchild=L; L->parent=P;
			tmp=L->LastSibling();
			if(tmp){
				tmp->rightchild=R; 
				if(R) R->parent=tmp;
			} else {
				L->rightchild=R;
				if(R) R->parent=L;
			}
		} else {
			P->leftchild=R; 
			if(R) R->parent=P;
		} 
	   } else print_error("RemoveNode( ): this should not happen");
	   this->rightchild=this->leftchild=this->parent=0; return TRUE;
	} 
	print_error("RemoveNode(): Cannot be applied to the root node");
	return FALSE;
}

btn_typ	*btn_typ::Fuse(Int4 nodeID)
// Fuse node with its parent.
{
	btn_typ *R=rightchild,*L=leftchild,*P=parent,*tmp;
	if(ID == nodeID){	// found node to be fused.
		// btn_typ *P=this->FindNaryParent(),*tmp;
		if(P){
		   if(this == P->rightchild){	// 'this' is within Sibling list.
			if(L){		// 'this' has a subtree attached.
				P->rightchild=L; L->parent = P;
				tmp=L->LastSibling(); 
				if(tmp){	// other siblings of L.
					tmp->rightchild=R; 
					if(R) R->parent=tmp;
				} else {	// L has no siblings.
					L->rightchild=R;
					if(R) R->parent=L;
				}
			} else{		// no subtree for 'this'.
				P->rightchild = R;  // simply delete 'this'
				if(R) R->parent=P;
			}
		   } else if(this == P->leftchild){  // 'this' = start of Sibling list.
			if(L){		// 'this' has a subtree.
				P->leftchild=L; L->parent=P;
				tmp=L->LastSibling();
				if(tmp){
					tmp->rightchild=R; 
					if(R) R->parent=tmp;
				} else {
					L->rightchild=R;
					if(R) R->parent=L;
				}
			} else {
				P->leftchild=R; 
				if(R) R->parent=P;
			} 
		   } else print_error("Fuse( ): this should not happen");
		   this->rightchild=this->leftchild=this->parent=0;
		   return this;
		} else print_error("Cannot fuse a root node with parent");
	} else {
		if(R){ tmp=R->Fuse(nodeID); if(tmp) return tmp; }
		if(L){ tmp=L->Fuse(nodeID); if(tmp) return tmp; }
		return 0;
	}
}

btn_typ *btn_typ::FindNaryParent( )
// Find the parent node treating the tree as an N-ary tree.
// To find the N-ary parent node follow parent nodes until one is reached that is a leftchild of it's parent.
// The N-ary parent is this leftchild's parent.
{
	if(this->parent == 0) return 0;						     // root node.
	if(this->parent->leftchild == this) return this->parent;		     // true parent...
	if(this->parent->rightchild == this) return this->parent->FindNaryParent( ); // binary parent is an N-ary sibling.
	assert(!"this should not happen");
}


/************************************************************
3. For each node (column in hpt) starting from the root:
   3.0. Allocate a list of length == # nodes in tree.
   3.3. Reset all nodes in the current node's subtree to '+'.    // foreground nodes.
   
4. Print out the matrix, labeling each row with the nodes number and name.
*************************************************************/

btn_typ	**btn_typ::OrderAsDFS(Int4 &rtnN, Int4 max_depth)
// order the nodes on a list using a DFS convention.
{
	Int4	Index=0;
	btn_typ **nodes;
	Int4 N=CountNodes(0);
	NEWP(nodes,N+3,btn_typ);

	Index++; nodes[Index] = this;
	if(max_depth <= 0) OrderDFS(Index, N, nodes);
	else OrderLimitedDFS(Index, N, nodes,0,max_depth);
	rtnN=Index;
	return nodes;
}

Int4    btn_typ::OrderLimitedDFS(Int4 &Index, Int4 N, btn_typ **nodes,
		Int4 depth, Int4 max_depth)
{
	if(depth < max_depth && leftchild != 0){   // then look at child nodes
		Index++; 
		if(Index > N) fprintf(stderr,"Index = %d\n",Index);
		assert(Index <= N); nodes[Index] = leftchild;
		leftchild->OrderLimitedDFS(Index,N,nodes,depth+1,max_depth); 
	} else if(leftchild){ delete leftchild; leftchild=0; } // delete subtree...
	if(rightchild != 0){			   // these are sibling nodes in N-ary tree.
		Index++;
		if(Index > N) fprintf(stderr,"Index = %d\n",Index);
		assert(Index <= N); nodes[Index] = rightchild;
		rightchild->OrderLimitedDFS(Index,N,nodes,depth,max_depth); 
	}
	return Index;
}

Int4    btn_typ::OrderDFS(Int4 &Index, Int4 N, btn_typ **nodes)
{
	if(leftchild != 0){
		Index++; 
		if(Index > N) fprintf(stderr,"Index = %d\n",Index);
		assert(Index <= N); nodes[Index] = leftchild; leftchild->OrderDFS(Index,N,nodes); 
	}
	if(rightchild != 0){
		Index++;
		if(Index > N) fprintf(stderr,"Index = %d\n",Index);
		assert(Index <= N); nodes[Index] = rightchild; rightchild->OrderDFS(Index,N,nodes); 
	}
	return Index;
}

Int4	btn_typ::Print(FILE *fp,Int4 depth)
// print as an N-ary tree
// Encoding n-ary trees as binary trees
// 
// There is a one-to-one mapping between general ordered trees and binary trees, which 
// in particular is used by Lisp to represent general ordered trees as binary trees. 
// Each node N in the ordered tree corresponds to a node N' in the binary tree; the left 
// child of N' is the node corresponding to the first child of N, and the right child 
// of N' is the node corresponding to N 's next sibling --- that is, the next node in 
// order among the children of the parent of N. This binary tree representation of a 
// general order tree, is sometimes also referred to as a First-Child/Next-Sibling 
// binary tree, or a Doubly-Chained Tree, or a Filial-Heir chain.
// 
// One way of thinking about this is that each node's children are in a linked list, chained 
// together with their right fields, and the node only has a pointer to the beginning or 
// head of this list, through its left field.
// 
// The binary tree can be thought of as the original tree tilted sideways, with the 
// black left edges representing first child and the blue right edges representing 
// next sibling. 
{
	// if this is the parent node's rightchild then print comma 
	Int4	N=1,n;

	// case 1: leaf as last child of parent.
	if(leftchild==0 && rightchild==0){	
	  fprintf(fp,"%d",ID);
	  if(name) fprintf(fp,"_%s",name);
//fprintf(fp,"{%d}",N);
	} // case 2: internal node as last child of parent.
	else if(leftchild!=0 && rightchild==0){	
		fprintf(fp,"("); N+=leftchild->Print(fp,depth+1); fprintf(fp,"):%d",ID);
	  	if(name) fprintf(fp,"_%s",name);
//fprintf(fp,"{%d}",N);
	} // case 3: leaf node as middle child of parent.
	else if(leftchild==0 && rightchild!=0){
		fprintf(fp,"%d",ID);
		if(name) fprintf(fp,"_%s",name);
//fprintf(fp,"{%d}",N);
		fprintf(fp,","); N += rightchild->Print(fp,depth+1);
	} // case 4: internal node and middle child of parent.
	else if(leftchild!=0 && rightchild!=0){
		fprintf(fp,"("); N+=leftchild->Print(fp,depth+1); fprintf(fp,"):%d",ID);
		if(name) fprintf(fp,"_%s",name);
//fprintf(fp,"{%d}",N);
		fprintf(fp,","); N += rightchild->Print(fp,depth+1);
//fprintf(fp,"|%d|\n",N);
	}
// fprintf(fp,"[%d(%d)]\n",N,depth);
	return N;
}

Int4	btn_typ::CountNodes(Int4 depth)
{
	Int4	N=1;

	// case 2: internal node as last child of parent.
	if(leftchild!=0 && rightchild==0) N+=leftchild->CountNodes(depth+1);
	// case 3: leaf node as middle child of parent.
	else if(leftchild==0 && rightchild!=0) N += rightchild->CountNodes(depth+1);
	// case 4: internal node and middle child of parent.
	else if(leftchild!=0 && rightchild!=0){
		N+=leftchild->CountNodes(depth+1); N+=rightchild->CountNodes(depth+1);
	}
	// if(leftchild==0 && rightchild==0) ; // case 1: leaf as last child of parent. 
	return N;
}

BooLean btn_typ::IsRightChild()
{
	if(parent != 0 && this == parent->rightchild) return TRUE; else return FALSE;
}

BooLean btn_typ::IsLeftChild()
{
	if(parent != 0 && this == parent->leftchild) return TRUE; else return FALSE;
}

Int4	btn_typ::NumRightChildren()
{
	if(rightchild){ return (1 + rightchild->NumRightChildren()); } else  return 0;
}

void	btn_typ::PrintNewick(FILE *fp)
// Print as a Newick tree.
{
	btn_typ *R=rightchild,*L=leftchild;
	if(L){			// attached to a 'subtree'. 
		fprintf(fp,"(");
		L->PrintNewick(fp);
        	// fprintf(fp,"):");
        	fprintf(fp,")");
	} 
	fprintf(fp,"%d",ID);
	if(name) fprintf(fp,"_%s",name);
	// fprintf(fp,"\n\n");
	if(R){			// has siblings. 
		fprintf(fp,",");
		R->PrintNewick(fp);
	} 
}

void	btn_typ::PrintBinary(FILE *fp)
// Print as a binary tree.
{
	if(leftchild && rightchild){	// is an internal node
	    fprintf(fp,"("); 
	    leftchild->PrintBinary(fp); 
	    fprintf(fp,","); 
	    rightchild->PrintBinary(fp); 
	    fprintf(fp,"):");
	} else if(leftchild){	// an internal node
	    fprintf(fp,"("); leftchild->PrintBinary(fp); fprintf(fp,"):");
	} else if(rightchild){	// an internal node
	    fprintf(fp,"("); rightchild->PrintBinary(fp); fprintf(fp,"):");
	} 
	fprintf(fp,"%d",ID);
	if(name) fprintf(fp,"_%s",name);
	// fprintf(fp,"\n\n");
}

btn_typ *btn_typ::ReturnNode(char *node_name)
// Find the node with name == node_name.
{
	btn_typ *left,*right;
	if(strcmp(name,node_name) == 0) return this;
	if(leftchild != 0){
		left=leftchild->ReturnNode(node_name);
		if(left != 0) return left;
	}
	if(rightchild){
		right=rightchild->ReturnNode(node_name);
		if(right != 0) return right;
	} return 0;
}


Int4	btn_typ::printSubHPT(FILE *fp, char *root_node_name, Int4 depth)
/************************************************************
In the hyperpartition, each column corresponds to one node in the tree.

The ‘+’ rows in that column correspond to that node’s subtree,
which serves as the foreground.

The ‘-‘ rows in that column correspond to the rest of the parent node’s subtree,
which serves as the background

The remaining (non-participating) nodes in that column are labeled with an ‘o’.

(For the root node a set of random sequences serves as the background.) 
*************************************************************/
{
	// 1. Check to make sure that each node has a unique name.
	btn_typ *subroot = ReturnNode(root_node_name),*pnode,*rnode;
	if(subroot != 0){
	    pnode=subroot->parent; subroot->parent=0;
	    rnode=subroot->rightchild; subroot->rightchild=0;
	    Int4 rtn=subroot->printHPT(fp,depth,20000);
	    subroot->parent=pnode; subroot->rightchild=rnode;
	    return rtn;
	} else {
		fprintf(stderr,"Error: Subtree node not found!\n");
		return 0;
	}
}



