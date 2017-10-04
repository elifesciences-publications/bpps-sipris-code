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

#include "block.h"
unsigned char *CARDINALITY_BLOCKS=NULL;

/*************** operations **********************/
int	IntersectBlockLCXOR(b_type eB,b_type B, b_type IB, b_type EOB)
/* Sets IIB = eB intersect B and sets EOB = B xor IB */
{
	register unsigned char *eb=eB->b,*b=B->b,*ib=IB->b,*eob=EOB->b;
	register unsigned short *list=B->list,*eolist=EOB->list;
	register unsigned short *newlist=IB->list;
	register unsigned char *card = CARDINALITY_BLOCKS;
	register int k,ik,n=0;

	for(n=ik=k=0; (*list) != END_BLOCK_LIST; list++){
	   if(ib[*list]=eb[*list] & b[*list]) {
		newlist[ik++]=*list; n+=card[(ib[*list])];
		if((eob[*list]=b[*list]^ib[*list])){ eolist[k++] = *list; }
	   } else { eolist[k++] = *list; }
	}
	eolist[k]= END_BLOCK_LIST;
	newlist[ik]= END_BLOCK_LIST;
	return n;
}

