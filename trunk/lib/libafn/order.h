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

#if !defined (ORDER)
#define	ORDER
#include "stdinc.h"
/************************ ADT Order ***********************************
Defines a model for the orderings of ntyp different types of linear 
elements having n (or npos) total positions in each of nseq sequences 
(of elements).

  pos:           0     1     2     3     4     n-1    n
  order:        ---[1]---[2]---[3]---[4]---....---[n]---

  seq[1]:	    A     B     A     C            B
  seq[2]:	    C     B     A     B            A
    :
    :
  seq[nseq]:	    B     C     A     B            A

 model[pos]:     ---[1]---[2]---[3]---....---[n]---
   A		  #A[1] #A[2] #A[3]  ....  #A[n] 
   B		  #B[1] #B[2] #B[3]  ....  #B[n] 
   C		  #C[1] #C[2] #C[3]  ....  #C[n] 

  pseudo:	   A*N0  ...etc...  (all positions have N0 of each element)
		   B*N0
		   C*N0 
N0 defines the prior probability N0 (in pseudocounts) that any 
element type will occur at any position in the ordering.
**********************************************************************/
/*************************** Order type **************************/
typedef struct {
	Int4	ntyp;		/* number of types of elements */
	Int4	npos;		/* total number of postions for elements */
	Int4	nseq;		/* number of sequences in the model */
	double	N0;		/* number of pseudo elements at each pos */
	double	**model; 	/* model[npos][ntyp] = # each type @ pos */
} order_type;
typedef order_type *o_type;

/******************************* private *****************************/
void    order_error(char *s);

/******************************* Public ******************************/
/***************************** operations ****************************/
o_type  Order(Int4 ntyp, Int4 npos, double npseudo);
void	PutOrder(FILE *fptr, o_type R);
Int4    *ConcensusOrder(o_type R);
void    InitOrder(o_type R);
void    Add2Order(Int4 *order, o_type R);
void    RmOrder(Int4 *order, o_type R);
double	RelProbOrder(Int4 *order, Int4 typ, Int4 pos, o_type R);
o_type  NilOrder(o_type R);

/**************************** macro operations **********************/
#define MAX_NUMBER_MOTIFS	26
#define nMotifsOrder(R)		((R)->ntyp)
#define nPosOrder(R)		((R)->npos)

#endif

