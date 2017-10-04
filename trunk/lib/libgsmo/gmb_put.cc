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

#include "gmb_typ.h"

void    gmb_typ::AlignSeq(FILE *fp,char mode, e_type qE,double temp)
{ SSX->AlignSeq(fp,mode,qE,temp); }

cma_typ	gmb_typ::CreateAlign(ss_type data)
// Align a set of input sequences to the existing cma file and return the new alignment.
{
	Int4	start,sq,score,*newpos,t;
	double	Temp=0.0;

	cma_typ	cma=SSX->RtnCMA( );
	gss_typ	*gss=gssCMSA(cma);
	cma_typ	rcma=EmptyCMSA(nBlksCMSA(cma), LengthsCMSA(cma), data,
			gss->GapOpen(),gss->GapExtend(), PerNatsCMSA(cma), 
			gss->LeftFlank(), gss->RightFlank());
	SSX->InitNDL(Temp); // initialize JHL parameters.
	gss_typ	*rgss=gssCMSA(rcma);
	for(sq=1; sq <= NSeqsSeqSet(data); sq++){
		e_type  sbjE = rgss->TrueSeq(sq);
		char *operation=SSX->GapAlnTrace(sbjE,Temp,start,score);
		// if(sq == 1){ this->PutGappedAln(stderr,sbjE,operation,start,cma); }
		Int4 trace_length=strlen(operation);
		gsq_typ *gsq = new gsq_typ[1];
		NEW(newpos,nBlksCMSA(rcma)+3,Int4);
		gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation,trace_length,start,sbjE,newpos);
		// gsq_typ *gsq=rgss->GetGSQ(sq);
		gsq_typ *ogsq=SwapGsqCMSA(sq,gsq,rcma); delete []ogsq;
		for(t=1; t<=nBlksCMSA(rcma); t++) AddSiteCMSA(t,sq,newpos[t],rcma);
		free(newpos); free(operation); 
	} return rcma;
}

cma_typ gmb_typ::CreateFullCMSA(ss_type Data, set_typ Set)
// Use 
{
        cma_typ cma=SSX->RtnCMA();
	ss_type data=TrueDataCMSA(cma);
        gsq_typ *gsq=0,*ogsq=0;
        gss_typ *gss=gssCMSA(cma);
	Int4    b,start,sq,score,*sites,J,trace_length;
	e_type	SqE,sqE;
	double	Temp=0.0;
	char	*operation;
        cma_typ CMA=EmptyCMSA(nBlksCMSA(cma),LengthsCMSA(cma),Data,gss->GapOpen(),
                        gss->GapExtend(),PerNatsCMSA(cma),gss->LeftFlank(),gss->RightFlank());
	SSX->InitNDL(Temp); // initialize JHL parameters.
        for(J=0,sq=1; sq <= NSeqsSeqSet(Data); sq++){
	    SqE=SeqSetE(sq,Data); 
            if(MemberSet(sq,Set)){
                J++; sqE=SeqSetE(J,data); assert(sqE==SqE);
		ogsq=gsqCMSA(J,cma); sites=GetPosSitesCMSA(J,cma);
                gsq = new gsq_typ[1]; gsq[0] = ogsq[0];
            } else {	// sample the sequence into the alignment...
		NEW(sites,nBlksCMSA(CMA)+3,Int4);
                operation=SSX->GapAlnTrace(SqE,Temp,start,score);
	    	trace_length=strlen(operation);
                gsq = new gsq_typ[1]; 
                gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation,trace_length,start,SqE,sites);
		free(operation); 
	    }
            ogsq=SwapGsqCMSA(sq,gsq,CMA); if(ogsq) delete []ogsq;
	    for(b=1; b<=nBlksCMSA(CMA); b++) AddSiteCMSA(b,sq,sites[b],CMA);
	    free(sites); 
        } return CMA;
}

cma_typ gmb_typ::CreatePurgedCMSA(Int4 percent_ident,e_type *Seqs,set_typ &Set, ss_type &data)
// remove similar sequences from the alignment, but retain the original alignments...
{
        cma_typ cma=SSX->RtnCMA();
	ss_type Data=TrueDataCMSA(cma);
        gsq_typ *gsq=0,*ogsq=0;
        gss_typ *gss=gssCMSA(cma);
	Int4    Nsq,Nset,i,j,b,sq,*sites,J;
	e_type	SqE,sqE;
	double	Temp=0.0;
	set_typ	InSet;
	char	*operation;

        Nsq=NumSeqsCMSA(cma); InSet=MakeSet(NumSeqsCMSA(cma)+4); FillSet(InSet);
        Set=RtnFastRepSetCMSA(stderr, percent_ident,InSet,cma); Nset=CardSet(Set);
        fprintf(stderr,"...file (\"%s\"): (%d/%d removed; %d remain).\n",
                                                NameCMSA(cma),Nsq-Nset,Nsq,Nset);
	NilSet(InSet);

	e_type	*RepSeqs; NEW(RepSeqs,Nsq+3,e_type); Nset=CardSet(Set);
        for(i=0,J=1; J <= Nsq; J++) if(MemberSet(J,Set)){ i++; RepSeqs[i]=Seqs[J]; }
	assert(i==Nset); data=Array2SeqSet(RepSeqs,Nset,NameCMSA(cma),AB);	// this frees up RepSeqs.

        cma_typ CMA=EmptyCMSA(nBlksCMSA(cma),LengthsCMSA(cma),data,gss->GapOpen(),
                        gss->GapExtend(),PerNatsCMSA(cma),gss->LeftFlank(),gss->RightFlank());
	SSX->InitNDL(Temp); // initialize JHL parameters.
        for(J=0,sq=1; sq <= NSeqsSeqSet(Data); sq++){
	    SqE=SeqSetE(sq,Data); 
            if(MemberSet(sq,Set)){
                J++; sqE=SeqSetE(J,data); assert(sqE==SqE);
		ogsq=gsqCMSA(sq,cma); sites=GetPosSitesCMSA(sq,cma);
                gsq = new gsq_typ[1]; gsq[0] = ogsq[0];
            	ogsq=SwapGsqCMSA(J,gsq,CMA); if(ogsq) delete []ogsq;
	        for(b=1; b<=nBlksCMSA(CMA); b++) AddSiteCMSA(b,J,sites[b],CMA);
	    	free(sites); 
	    }
        } return CMA;
}


