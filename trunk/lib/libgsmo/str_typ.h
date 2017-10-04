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

#if !defined (STR_TYP)
#define STR_TYP
#include "stdinc.h"

// Header file for data structure representing a collection of lists that
// partition a set of small positive integers. The lists are circularly
// linked in both directions. There is no concept of a first or last
// node of a list. We can refer to a list by referring to any of its
// elements.

class str_typ {
	char	label;			// list defined on 'm' versus 'd'.
	str_typ *next;			// next[i] is successor of i in list
	str_typ	*prev;			// prev[i] is predecessor of i in list
public:		str_typ(char);
		~str_typ();
	str_typ	*Append(char);		// create a new link and insert after 'this'.
	str_typ	*Next() { return next; }	// return successor
	str_typ	*Prev() { return prev; }	// return predecessor
	char    *Return();
	void	Print(FILE *fp=stderr);		// print the list starting from 'this'
};

#endif

