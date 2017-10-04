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

// Header file for a binary search tree data structure.
// MaInt4ains a subset of items {1,...,n}, where each item has a key.
// Code for these routines was modeled on that in Chapter 13
// (Binary Search Trees) of Cormen, Leiserson, and Rivest's
// "Introduction to Algorithms".

typedef Int4 keytyp;
typedef Int4 item;

class bstree {
	Int4	N;			// max number of items in tree
	Int4	n;			// number of items in tree
	item	root;			// item at root of tree
	item	*leftc, *rightc, *p;	// left child, right child, and
					//    parent of each item
	keytyp	*kvec;			// kvec[i] is key of item i
	item	minimum(item);		// return smallest item in subtree
	item	maximum(item);		// return largest item in subtree
	void	clear2(item);		// recursive func. called by clear
public:		bstree(Int4=100);
		~bstree();
	keytyp	key(item);		// return the key of item
	char	member(item);		// return true if item in tree
	char	empty();		// return true if tree is empty
	void	insert(item,keytyp);	// insert item with specified key
	void	remove(item);		// remove item from tree
	item	search(keytyp);		// return an item with a given key
	item	findmin();		// return an item with smallest key
	item	findmax();		// return an item with largest key
	item	suc(item);		// return successor
	item	pred(item);		// return predecessor
	void	clear();		// remove everything
	void	print();		// print the tree	
	void	toNary(FILE *fp);	// afn: for my code.
	void	fromNary(char *str);	// afn: for my code.
};

// Return key of i.
inline keytyp bstree::key(item i) { return kvec[i]; }

// Return true if i in tree, else false.
inline char bstree::member(item i) { return (i <= N && i > 0 && p[i] != -1); }

// Return true if heap is empty, else false.
inline char bstree::empty() { return n == 0; }
