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

#include "brh_typ.h"

brh_typ	NilBRH(brh_typ result)
{
	hsp_typ		hsp;
	register Int4	i;
        
        if (result == NULL) return NULL;
        for(i=0; i < result->hspcnt; i++) {
          hsp = result->hsp_array[i]; NilHSPTyp(hsp);
        }
        if(result->hsp_array != 0) free(result->hsp_array);
        free(result);
        return NULL;
}

brh_typ	MakeGBLASTResultHitlist(Int4 size, e_type sE)
{
        brh_typ	result;

        NEW(result,1,brh_type);
	NEW(result->hsp_array,size+1,hsp_typ);
	result->best_evalue = DBL_MAX;
	result->high_score= INT4_MIN;
	result->hspcnt=0;
	result->sE=sE;
	result->subject_id=SeqI(sE);
	result->array_max=size;
        return result;
}

sap_typ	ExtractAlnBRH(brh_typ result, e_type qE)
{
	sap_typ	sap,sap2,head=0;
	hsp_typ		hsp;

	for(Int4 h=0; h < result->hspcnt; h++){
	  hsp = result->hsp_array[h];
	  assert(hsp->gap_info != NULL);
	  sap=GXEBToGSeqAlign(hsp->gap_info,result->sE, qE);
	  GXEBDelete(hsp->gap_info); hsp->gap_info=0;
	  sap->evalue = hsp->e_value;
	  sap->score = hsp->score;
	  sap->bit_score = hsp->bit_score;
	  // sap->segs->subject_id = result->subject_id;
// fprintf(stderr,"sap->segs->subject_id = %d; result->subject_id = %d; %d\n",
//		sap->segs->subject_id,result->subject_id,SeqI(result->sE));
	  assert(sap->segs->subject_id == result->subject_id);
	  if(head == 0) { head = sap; }
          else { for(sap2=head;sap2->next;) sap2=sap2->next; sap2->next=sap; }
	}
	return head;
}

void	AddHspBRH(hsp_typ hsp, brh_typ result)
{
	assert(result->hspcnt < result->array_max);
	result->hsp_array[result->hspcnt] = hsp;
	result->hspcnt++;
	if(result->best_evalue > hsp->e_value) 
		result->best_evalue=hsp->e_value;
	if(result->high_score < hsp->score) 
		result->high_score=hsp->score;
// std::cerr << log10(BestEvalBRH(result)); std::cerr << std::endl;
// std::cerr << hsp->score; std::cerr << std::endl;
}

double  BestEvalBRH(brh_typ result){ return result->best_evalue; }

//======================= BRHHeap( ) ==========================

rhp_typ MakeBRHHeap(Int4 hpsz)
{
	rhp_typ	H;

	NEW(H,1,brhheap_type);
	NEW(H->brh_array,hpsz+2,brh_typ);
	H->mH=Mheap(hpsz, 3);
	H->size = hpsz;
	return H;
}

void    NilBRHHeap(rhp_typ H)
{
	brh_typ	brh; double	key;

	while((brh=DelMinBRHHeap(&key,H))) NilBRH(brh); 
        NilMheap(H->mH); free(H->brh_array);
	free(H);
}

Int4    InsertBRHHeap(brh_typ result, rhp_typ H)
{
	// double key = result->best_evalue;
	double key = result->high_score;

        // Int4 item=InsertMheap(key, H->mH);
        Int4 item=InsertMheap(-key, H->mH);
        if(item==NULL) return NULL;
        if(H->brh_array[item]!= NULL) NilBRH(H->brh_array[item]); 
        H->brh_array[item]=result; 
        return item;
}

brh_typ SeeBRHHeap(Int4 item, rhp_typ H) { return H->brh_array[item]; }

brh_typ DelMinBRHHeap(double *key, rhp_typ H)
{
        *key= MinKeyMheap(H->mH); 
        // *key= -MinKeyMheap(H->mH); 
	Int4 item = DelMinMheap(H->mH);
	brh_typ M = H->brh_array[item]; H->brh_array[item]=NULL;
        return M; 
}


