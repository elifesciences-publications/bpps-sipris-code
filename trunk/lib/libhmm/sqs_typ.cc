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

#include "sqs_typ.h"

#define MAX_IN_SEQS	1000000

sqs_typ::sqs_typ(char *DBS_NAME,a_type A)
{
	AB=A;
	num_seqs=0;
	Init(DBS_NAME,MAX_IN_SEQS);
}

//******************** Add these later to sequence.h *******************
const unsigned char default_overlap = 10;	// Already in sequence.c

BooLean	AllocEconoSeq(e_type E,UInt4 I,unsigned char *S, char *info,unsigned short leng)
{
	Int4	i;
	NEW(E,1,sequence_type);
        E->S = S;
	for(i=-default_overlap; i<=0; i++) E->S[i]=0;
        E->overlap=default_overlap; E->offset=E->extend=0;
        E->info=info;  E->n = leng;
        E->I = I; E->X = E->S; E->xed = FALSE; E->rff=0;
        return TRUE;
}

e_type  *MkEconoSeq(Int4 num_seqs, unsigned char **S,char **info, unsigned short *length)
{
	e_type *E;
	MEW(E,num_seqs,e_type);
//	for(Int4 i=1;i<=num_seqs;i++) MEW(E[i],1,sequence_type);
	for(Int4 i=1; i <= num_seqs; i++){
		AllocEconoSeq(E[i],i,S[i],info[i],length[i]);
	} return E;
}

Int4	MaxSeqDeflineLen( ){  return 100; }

Int4	DefaultSeqOverlap(){ return default_overlap; }

void    NilEconoSeq(e_type E)
// Do the inverse of AllocEconoSeq;
{
        if(E!=NULL){
                if(E->xed) { E->X-= default_overlap; free(E->X); }
		E->S -= default_overlap; free(E->S);
                if(E->info!=NULL)free(E->info);
                if(E->rff) delete E->rff; free(E);
        }

}

//******************** Add these later to sequence.h *******************

void	sqs_typ::Init(char *DBS_NAME,Int4 max)
{
	Int4	i;
	num_seqs=GetFastaInfo(DBS_NAME,max,&counts,&seq_len,AB);// half of the time is spent here
	total_residues=0;
	for(i=0; i <= nAlpha(AB); i++){ total_residues+=counts[i]; }
	MEW(HugeSeqMem,total_residues+num_seqs*DefaultSeqOverlap()*3,unsigned char);
	MEW(HugeInfoMem,MaxSeqDeflineLen()*(num_seqs+3),char);
	NEWP(seq_ptr,num_seqs+3,unsigned char);
	NEWP(seq_info,num_seqs+3,char);
	unsigned char *tmp_seq_ptr=HugeSeqMem+DefaultSeqOverlap();
	char *tmp_seq_info=HugeInfoMem;
	FILE *fp=open_file(DBS_NAME,"","r");
	for(i=1; i <= num_seqs; i++){ // half of the time is spent here
		seq_ptr[i]=tmp_seq_ptr;
		seq_info[i]=tmp_seq_info;
		ReadSeqStr(fp, seq_info[i],seq_ptr[i],seq_len[i]+2,AB);
		tmp_seq_ptr += seq_len[i]+DefaultSeqOverlap()*2;
		tmp_seq_info += MaxSeqDeflineLen();
	}
	fclose(fp);
	E=MkEconoSeq(num_seqs,seq_ptr,seq_info,seq_len); // little time is spent here
}
