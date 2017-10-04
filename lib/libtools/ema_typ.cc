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

#include "ema_typ.h"

ema_typ	MakeEMA(Int4 fsq, Int4 ssq, Int4 start, Int4 end)
{ return initEMA(fsq,ssq,FALSE,start,end,0,0.0); }

ema_typ	MakeEMA(Int4 fsq, Int4 ssq, Int4 start, Int4 end, Int4 score,
	double Evalue)
{ return initEMA(fsq,ssq,TRUE,start,end,score,Evalue); }

ema_typ	initEMA(Int4 fsq, Int4 ssq, BooLean new_seq, Int4 start, 
	Int4 end, Int4 score, double Evalue)
{
	ema_typ	ema; NEW(ema,1,editcma_type_type);
	ema->next=0;
	ema->fsq=fsq; ema->ssq=ssq;
	ema->new_seq=new_seq;
	ema->start=start; ema->end=end;
	ema->score=score; ema->Evalue=Evalue;
	return ema;
}

#if 0
static Int4	overlap_ema(ema_typ ema1, ema_typ ema2)
{
	assert(ema2 && ema1);
	if(ema1->fsq != ema2->fsq) return 0;
	Int4 ol;
	
	if(ema1->start!= ema2->ssq) return 0;
}
#endif

ema_typ	PurgeEMA(ema_typ ema)
{
	ema_typ	tmp,nxt;
	for(tmp=ema; tmp; ){
	   if(tmp->next){
		nxt = tmp->next;
		if(tmp->fsq == nxt->fsq && tmp->start == nxt->start){
		   tmp->next = nxt->next; nxt->next=0; 
		   PutEMA(stderr,nxt); NilEMA(nxt);
		} else tmp= tmp->next;
	   } else break;
	} return ema;
}

ema_typ	AppendEMA(ema_typ ema, ema_typ ema2)
// add ema2 to list in such a way as to maintain correct ordering.
// delete overlaps as they occur...
{
	assert(ema2 && ema2->next == 0);
	if(ema == 0) return ema2;
	if(ema2->fsq < ema->fsq){ 
	    ema2->next = ema; ema = ema2;
	} else if(ema2->fsq == ema->fsq){   // then check subseq, etc...
	    if(ema2->ssq < ema->ssq){
		   ema2->next = ema; ema = ema2; 
	    } else if(ema->ssq == ema2->ssq){ // Only case where 'new_seq' matters.
		if(!ema->new_seq && ema2->new_seq) ema->next=AppendEMA(ema->next,ema2);
		else if(ema->new_seq && ema2->new_seq){
		   if(ema->start > ema2->start){ ema2->next=ema; ema=ema2; }
		   else if(ema->start == ema2->start){
			NilEMA(ema2); // ema2 is redundant; destroy it.
		   } else {    // ema->start < ema2->start
			ema->next = AppendEMA(ema->next, ema2);
		   }  
		} else if(ema->new_seq && !ema2->new_seq)
		   print_error("AppendEMA( ): can't add old seq after new_seqs");
		else // if(!ema->new_seq && !ema2->new_seq)
		   print_error("AppendEMA( ) old seqs need unique identifiers");
	    } else if(ema->next == 0) ema->next=ema2; 
	    else { // append further down the list...
		   ema->next = AppendEMA(ema->next, ema2);
	    }
	} else if(ema->next == 0){ ema->next=ema2; }
	else {	// ema2->fsq > ema->fsq...
	    ema->next = AppendEMA(ema->next, ema2);
	} return ema;
}

void    PutEMA(FILE *fp, ema_typ ema)
{
	Int4	N;
	ema_typ	tmp;
	for(N=0,tmp=ema; tmp; tmp= tmp->next){
	   if(!tmp->new_seq){
		fprintf(fp,"sq %d(%d): %d-%d\n",
		    tmp->fsq,tmp->ssq,tmp->start,tmp->end);
	   } else {
		fprintf(fp,"sq %d(%d..%d): %d-%d (%d -> %.2g)\n",
		    tmp->fsq,tmp->ssq,tmp->ssq+1,tmp->start,tmp->end,
			tmp->score,tmp->Evalue);
	   }
	   N++;
	}
}

Int4	LengthEMA(ema_typ ema)
{
	Int4	N;
	ema_typ	tmp;

	for(N=0,tmp=ema; tmp; tmp= tmp->next) N++;
	return N;
}

void	NilEMA(ema_typ ema)
{
	if(ema->next) NilEMA(ema->next);
	free(ema);
}

