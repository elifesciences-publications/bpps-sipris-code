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

// set_typ	gmb_typ::RtnTestSet(set_typ Set)
Int4	*gmb_typ::RtnTestSet(set_typ Set)
// Return a set containing sequences not in Set but with similar scores to the Consensus for Set.
{
	FILE	*fp=stderr;
        cma_typ cma=SSX->RtnCMA();
        double  dd,d,DD,dm,*oldpos;
        Int4    sq,i,r,r0,scr,min,Nadd,N=NumSeqsCMSA(cma);
	// set_typ	tSet=CopySet(Set); ClearSet(tSet);
	Int4	*tSet; NEW(tSet,N+3,Int4);
        h_type  HG=Histogram("in-set sequence to consensus scores",0,500,10);
        e_type  CsqE=MKCsqSubCMSA(Set,cma);
        for(min=INT4_MAX,sq=1; sq <= N; sq++){
            assert(SeqIsInCMSA(sq,cma));
            if(MemberSet(sq,Set)){
                for(scr=0,i=1; i <= LengthCMSA(1,cma); i++){
                   r0=ResSeq(i,CsqE); r=ResidueCMSA(1,sq,i,cma);
                   scr += valAlphaR(r,r0,AB);
                } IncdHist((double)scr,HG);
                if(scr < min) min=scr;
            }
        } dm=(double)min;
        h_type  xHG=Histogram("out-of-set sequence to consensus scores",-200,200,10);
	if(min < 0) min=0; min = (Int4) floor(0.90*(double)min);
	dh_type dH=dheap(N+2,4);
	Int4	II=0;
        for(sq=1; sq <= N; sq++){
            if(!MemberSet(sq,Set)){
                for(scr=0,i=1; i <= LengthCMSA(1,cma); i++){
                   r0=ResSeq(i,CsqE); r=ResidueCMSA(1,sq,i,cma);
                   scr += valAlphaR(r,r0,AB);
                } IncdHist((double)scr,xHG);
		// if(scr >= min) AddSet(sq,tSet);
		if(scr >= min){ II++; tSet[II]=sq; }
		else insrtHeap(sq,(keytyp)-scr,dH);
            }
        } 
	// Nadd=(CardSet(Set)/2) -CardSet(tSet);
	dd=MeanHist(xHG) + 2.0*sqrt(VarianceHist(xHG));	// 2 stdev above mean...
	while(!emptyHeap(dH)){	// Put at least half as many in tSet as in Set.
		d=-(double)minkeyHeap(dH);
		if(d < dd && d < dm) break;
		else { sq=delminHeap(dH); II++; tSet[II]=sq; }
		// else { sq=delminHeap(dH); AddSet(sq,tSet); }
	} PutHist(fp,50,HG); NilHist(HG); PutHist(fp,50,xHG); NilHist(xHG); 
	NilSeq(CsqE); Nildheap(dH);
	return tSet;
}

// void	gmb_typ::PutSampleOutSqIn(FILE *fp,double Temp,set_typ Set, set_typ tstSet)
void	gmb_typ::PutSampleOutSqIn(FILE *fp,double Temp,set_typ Set, Int4 *tstSet)
{
	cma_typ	cma=SSX->RtnCMA();
	double	dd,d,DD,dm;
	Int4	i,j,sq,*oldpos;
	h_type  xHG=Histogram("out-of-set sequence LLR contributions",-200,+200,10);
#if 0
	Int4	i,r,r0,scr,min;
	h_type  HG=Histogram("in-set sequence LLR contributions",-200,+200,10);
	e_type  CsqE=MKCsqSubCMSA(Set,cma);
	for(min=INT4_MAX,sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    // if(!SeqIsInCMSA(sq,cma))
	    if(MemberSet(sq,Set)){
		for(scr=0,i=1; i <= LengthCMSA(1,cma); i++){ 
		   // if(AddOp[i] != 'd') continue;
		   r0=ResSeq(i,CsqE); r=ResidueCMSA(1,sq,i,cma);
		   scr += valAlphaR(r,r0,AB);
		} IncdHist((double)scr,HG);
		if(scr < min) min=scr;
	    }
	} PutHist(fp,50,HG); NilHist(HG); dm=(double)min;
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    if(!MemberSet(sq,Set)){
		d=this->SampleOutSqIn(sq,Temp,dd,CsqE,min);
		if(d >= dm){ AddSet(sq,Set); }
		else {
	   	   Int4	*oldpos=SSX->RemoveFromAlign(sq); free(oldpos);
		} IncdHist((double)d,xHG);
	    }
	} PutHist(fp,50,xHG); NilHist(xHG); 
#elif 0
	dm=0.0; dd=SSX->GapMap(Temp);
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    if(!MemberSet(sq,Set) && MemberSet(sq,tstSet)){
		d=this->SampleOutSqIn(sq,Temp,dd); DD=d-dd;
		if(DD >= dm){ AddSet(sq,Set); dd=d; }
		else { oldpos=SSX->RemoveFromAlign(sq); free(oldpos); }
		IncdHist((double)DD,xHG);
	    }
	} PutHist(fp,50,xHG); NilHist(xHG); 
#elif 0
	dd=SSX->GapMap(Temp);
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    if(!MemberSet(sq,Set) && MemberSet(sq,tstSet)){
		DD=dd;
		if(this->SampleOutSqIn(sq,Temp,dd)){
			AddSet(sq,Set); IncdHist((dd-DD),xHG); DD=dd;
		}
	    }
	} PutHist(fp,50,xHG); NilHist(xHG); 
#else	// Use a list instead...
	dd=SSX->GapMap(Temp);
	for(j=0,i=1; (sq=tstSet[i]) != 0; i++){
	    assert(sq <= NumSeqsCMSA(cma));
	    if(!MemberSet(sq,Set)){
		DD=dd;
		if(this->SampleOutSqIn(sq,Temp,dd)){
			j=0; AddSet(sq,Set); IncdHist((dd-DD),xHG); DD=dd;
		} else j++;
		if(j > 7) break;
	    }
	} PutHist(fp,50,xHG); NilHist(xHG); 
#endif
}

BooLean	gmb_typ::SampleOutSqIn(Int4 sq, double Temp, double &MAP, e_type CsqE,Int4 min)
// Sample sequence sq (which is not in the alignment), into the alignment.
{
	Int4    i,start,trace_length,score,*newpos;
        cma_typ CMA=SSX->RtnCMA();

        assert(!SeqIsInCMSA(sq,CMA));
        gss_typ *gss=gssCMSA(CMA);
        e_type  sbjE = gss->TrueSeq(sq);
        SSX->InitNDL(Temp);     // initialize JHL parameters?
        char    *operation=SSX->GapAlnTrace(sbjE,Temp,start,score);
        trace_length=strlen(operation);
        //========== 2. Create a gapped sequence using 'operation'. ===========
        gsq_typ *ogsq,*gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(CMA)+3,Int4);
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation,trace_length,start,sbjE,newpos);
        free(operation);
SSX->SetTempSqWt(sq);
        ogsq=SwapGsqCMSA(sq,gsq,CMA); SSX->AddToAlign(sq,newpos);
	free(newpos); delete []ogsq; 
	double dd=SSX->GapMap(Temp);
	if(dd-MAP > 0.0){	// Need to recalculate Sequence weights and indel penalties!
// SSX->UnSetTempSqWt(sq); // reinitialized, so no need to unset.
#if 0
		MAP=dd; delete SSX;
		SSX = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pernats,CMA,dms_mode);
#else
		MAP=dd; SSX->ReInitCnts();
#endif
		return TRUE;
	} else {
		Int4 *oldpos=SSX->RemoveFromAlign(sq); free(oldpos); 
SSX->UnSetTempSqWt(sq);
		return FALSE;
	} 
#if 0
	Int4	r,r0,scr=0;
	for(i=1; i <= LengthCMSA(1,CMA); i++){ 
		// if(AddOp[i] != 'd') continue;
		r0=ResSeq(i,CsqE); r=ResidueCMSA(1,sq,i,CMA);
		scr += valAlphaR(r,r0,AB);
	}
	if(scr < min){
	   Int4	*oldpos=SSX->RemoveFromAlign(sq); free(oldpos);
	} return scr;
#endif
}

