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

// Header file for a binary tree node data structure.
#if !defined (BTN_TYP)
#define BTN_TYP

#include "stdinc.h"
#include "afnio.h"
#include "hpt_typ.h"

class btn_typ  {
public:		
		btn_typ(Int4 id, char *Name);
		~btn_typ();
	char	*Name() { return name; }
	Int4	Ident() { return ID; }
	void	Isolate() { leftchild=rightchild=parent=0; }
	Int4	MaxID();		// return the maximum id in subtree
	void	AddLeftChild(btn_typ *child) { leftchild=child; child->parent=this; }
	void	AddRightChild(btn_typ *child) { rightchild=child; child->parent=this; }
	void	SetParent(btn_typ *node) { parent=node; }
	void	PrintBinary(FILE *fp);	// print as a binary tree in standard format
	void	PrintNewick(FILE *fp);	// print as a tree in Newick format
	Int4	Print(FILE *fp,Int4 depth);	// print as an N-ary tree in standard format
	hpt_typ *ReturnHPT(Int4 depth, Int4 NumRandom);
	Int4	PrintHPT(FILE *fp,Int4 depth,Int4 NumRandom){
			btn_typ *tmp=this->CopyTree();
			tmp->printHPT(fp,depth,NumRandom); 
			delete tmp;
		}
	Int4    PrintSubHPT(FILE *fp, char *root_node_name, Int4 depth){
			btn_typ *tmp=this->CopyTree();
			tmp->printSubHPT(fp, root_node_name, depth);
			delete tmp;
		}
	BooLean	IsOkay();		// check for correct properties.
	BooLean	IsRightChild();		
	BooLean	IsLeftChild();		
	BooLean	DeletedNode( ){ return (parent == 0 && leftchild==0 && rightchild==0); }
	BooLean	IsRoot(){ return (parent == 0); }
	BooLean	IsLeafNode(){ return (leftchild==0); }
	btn_typ	*LeftChild() { return leftchild; }
	btn_typ	*RightChild() { return rightchild; }
	btn_typ	*Parent(){ return parent; }
	Int4    NumRightChildren();
	Int4    CountNodes(Int4 depth);
	void	Put(FILE *fp){
		   btn_typ *R=rightchild,*L=leftchild,*P=parent;
		   if(this != 0) fprintf(fp,"%d: %s",ID,name);
		   if(P) fprintf(fp," parent=%d",P->ID);
		   if(R) fprintf(fp," rightchild=%d",R->ID);
		   if(L) fprintf(fp," leftchild=%d",L->ID);
		   fprintf(fp,"\n");
		}
	btn_typ	*Fuse(Int4 nodeID);	// return the node to be deleted; or null if not found.
	btn_typ	*LastSibling();
	btn_typ *CopyTree( );		// return a copy of 'this' tree.
	btn_typ *SplitOffLeaf(btn_typ *Node, Int4 nodeID);
	btn_typ *SplitOffParent(btn_typ *Node, Int4 nodeID);

	BooLean	AddLeaf(btn_typ *Node);
	BooLean AddParent(btn_typ *Node);
	BooLean	RemoveNode( );

	btn_typ	**ReturnAsList(Int4 &rtnN){ return OrderAsDFS(rtnN,0); }
	Int4	*GetSubTree() { return GetSubTreeIDs(); } 
	btn_typ *RtnNaryParent( ){ return this->FindNaryParent( ); }
private:		
	Int4	printHPT(FILE *fp,Int4 depth,Int4 NumRandom);	// print as a hyperpartition
	Int4    printSubHPT(FILE *fp, char *root_node_name, Int4 depth);
	btn_typ *ReturnNode(char *node_name);
	btn_typ *FindNaryParent( );
	Int4	SubTreeIDs(Int4 &Index, Int4 N,Int4 *list);
	Int4	*GetSubTreeIDs();
	btn_typ **OrderAsDFS(Int4 &rtnN, Int4 max_depth);
	Int4    OrderDFS(Int4 &Index, Int4 N, btn_typ **nodes);
	Int4    OrderLimitedDFS(Int4 &Index, Int4 N, btn_typ **nodes,
			Int4 depth, Int4 max_depth);
	char    **GetHPT(Int4 &N);
	Int4	ID;			// integer identifier
	char	*name;			// name of protein family
	btn_typ	*leftchild;		// = subtree of node
	btn_typ	*rightchild;		// == siblings of node.
	btn_typ	*parent;		
};



#endif

