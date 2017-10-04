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
#include "bstree.h"

bstree::bstree(Int4 N1) {
	N = N1; n = 0;
	leftc = new item[N+1]; rightc = new item[N+1]; p = new item[N+1];
	kvec = new keytyp[N+1];
	root = 0;
	for (item i = 1; i <= N; i++) {
		leftc[i] = rightc[i] = 0;
		p[i] = -1;
	}
	leftc[0] = rightc[0] = 0;
	p[0] = -1;
}

bstree::~bstree() { delete leftc; delete rightc; delete p; delete kvec; }

void bstree::clear2(item i) {
// Traverse the elements of the tree, deleting them in postorder
	if (i == 0) return;
	clear2(leftc[i]);
	clear2(rightc[i]);
	leftc[i] = rightc[i] = 0;
	p[i] = -1;
}

void bstree::clear() {
	clear2(root); root = 0; n = 0;
}

void bstree::insert(item i, keytyp k) {
// Add item i with key k to tree.
	n++;
	kvec[i] = k; leftc[i] = 0; rightc[i] = 0;
	item y = 0;
	item x = root;
	while (x != 0) {
		y = x;
		if (kvec[i] < kvec[x]) x = leftc[x];
		else x = rightc[x];
	}
	p[i] = y;
	if (y == 0) { root = i; }
	else {
		if (kvec[i] < kvec[y]) leftc[y] = i;
		else rightc[y] = i;
	}
}

void bstree::remove(item i) {
// Remove item i from tree.
	item x,y;
	n--;
	if (leftc[i] == 0 || rightc[i] == 0) {
		y = i;
	} else {
		y = suc(i);
	}
	if (leftc[y] != 0) {
		x = leftc[y];
	} else {
		x = rightc[y];
	}
	if (x != 0) p[x] = p[y];
	if (p[y] == 0) {
		root = x;
	} else {
		if (y == leftc[p[y]]) {
			leftc[p[y]] = x;
		} else {
			rightc[p[y]] = x;
		}
	}
	if (y != i) {
	// Instead of copying key from y to i as in Sec 13.3 of CLR,
	// copy the tree poInt4ers from i to y.
	// For the reason why, see Exercise 13.3-5.
		p[y] = p[i];
		if (p[i] == 0) {
			root = y;
		} else {
			if (i == leftc[p[i]]) {
				leftc[p[i]] = y;
			} else {
				rightc[p[i]] = y;
			}
		}
		leftc[y] = leftc[i];
		if (leftc[y] != 0) p[leftc[y]] = y;
		rightc[y] = rightc[i];
		if (rightc[y] != 0) p[rightc[y]] = y;
	}
	p[i] = -1;
	leftc[i] = rightc[i] = 0;
}

item bstree::search(keytyp k) {
// Return an item with key k, or 0 if none exists.
	item x = root;
	while (x != 0 && k != kvec[x]) {
		if (k < kvec[x])
			x = leftc[x];
		else
			x = rightc[x];
	}
	return x;
}

item bstree::minimum(item i)
// Return an item with smallest key in the subtree rooted at item i,
// or 0 of i is not in the tree.
{
	item	x;
	if (i == 0 || p[i] == -1) return 0;
	for (x = i; leftc[x] != 0; x = leftc[x]) {}
	return x;
}

item bstree::findmin() {
// Return an item with smallest key value, or 0 if tree is empty.
	return minimum(root);
}

item bstree::maximum(item i)
// Return an item with largest key in the subtree rooted at item i,
// or 0 of i is not in the tree.
{
	item	x;
	if (i == 0 || p[i] == -1) return 0;
	for (x = i; rightc[x] != 0; x = rightc[x]) {}
	return x;
}

item bstree::findmax() {
// Return an item with largest key value, or 0 if tree is empty.
	return maximum(root);
}

item bstree::suc(item i) {
// Return the successor of item i in the tree,
// or 0 if i is the last item.
	if (i == 0 || p[i] == -1) return 0;
	if (rightc[i] != 0)
		return minimum(rightc[i]);
	item j = p[i];
	while (j != 0 && i == rightc[j]) {
		i = j; j = p[j];
	}
	return j;
}

item bstree::pred(item i) {
// Return the predecessor of item i in the tree,
// or 0 if i is the first item.
	if (i == 0 || p[i] == -1) return 0;
	if (leftc[i] != 0)
		return maximum(leftc[i]);
	item j = p[i];
	while (j != 0 && i == leftc[j]) {
		i = j; j = p[j];
	}
	return j;
}

void bstree::print() {
// Print the contents of the tree.
	Int4 i;
	printf("n=%2d root=%2d\nitem:", n, root);
	for (i = 1; i <= N; i++) printf(" %2d", i);
	printf("\nkvec:");
	for (i = 1; i <= N; i++) printf(" %2d", kvec[i]);
	printf("\n   p:");
	for (i = 1; i <= N; i++) printf(" %2d", p[i]);
	printf("\n  lc:");
	for (i = 1; i <= N; i++) printf(" %2d", leftc[i]);
	printf("\n  rc:");
	for (i = 1; i <= N; i++) printf(" %2d", rightc[i]);
#if 0	// N-ary tree conversion
	printf("\n  Level:");
	for (i = 1; i <= N; i++) printf(" %2d", rightc[i]);
#endif
	putchar('\n');
}


