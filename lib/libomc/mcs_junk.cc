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

#include "mcs_typ.h"

void    mcs_typ::SetTheCSQ(FILE *fp,Int4 n,e_type csqE)
{
	Int4	j,ri,ris,Length;
	e_type  keyE=0;

	if(fp) fprintf(fp,"=========== Optimizing evolving consensus sequence ===========\n");
	assert(n > 0 && n <= Hpt->NumBPPS()); 
	// PutSeq(stderr,che[n]->KeyE( ),AB);
	// AlnSeqSW(11, 1, che[n]->KeyE( ),csqE,AB);
	if(csqE == 0) che[n]->UpdateConsensus(); else che[n]->SetConsensus(csqE);
	// PutSeq(stderr,che[n]->KeyE( ),AB); PutSeq(stderr,csqE,AB); exit(1);
	// AlnSeqSW(11, 1, che[n]->KeyE( ),csqE,AB);
	keyE=che[n]->KeyE( );
	bpps_typ *pps=che[n]->BPPS();  pps->UpdateRho(sst[n],keyE); 
#if 0
	Int4 mc=che[n]->RtnMaxNumColumns(); che[n]->SetMaxNumColumns(Length);
	for(j=1; j <= Length; j++){
		che[n]->RemoveColumn(j); che[n]->AddColumn(j,sst[n][j][1]); 
	} che[n]->SetMaxNumColumns(mc);
#endif
}

double  mcs_typ::SampleLeafParent(Int4 n)
{
        assert(n > 0 && n < Hpt->NumSets());
        // double d=che[n]->Sample(FALSE,sst[n]);
        double d=che[n]->Sample(TRUE,sst[n]);
        che[n]->PutInfo(stderr);
        che[n]->PutSubLPRs(stderr); return d;
}

