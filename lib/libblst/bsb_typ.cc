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

#include "bsb_typ.h"

bsb_typ	GBlastSearchBlkNew(Int4 results_size)
{
        bsb_typ	New;
        NEW(New,1,bsb_type);
        NEW(New->results,results_size+2,brh_typ);
        New->hitlist_max = results_size; New->hitlist_count = 0;
        return New;
}

void	ResetGBlastSearchBlk(bsb_typ bsb) // Free results for another round.
{

        if(bsb != NULL && bsb->results != NULL){
           brh_typ *results = bsb->results;
           for(Int4 i=0; i < bsb->hitlist_max; i++){
                if(results[i]) results[i] = NilBRH(results[i]);
           } bsb->hitlist_count = 0;
        } 
}

bsb_typ	GBlastSearchBlkDelete(bsb_typ bsb)
{
        if(bsb == NULL) return NULL;
	ResetGBlastSearchBlk(bsb);
        if(bsb->results != NULL) GMemFree(bsb->results); free(bsb);
	return NULL;
}

