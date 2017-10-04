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

void    gmb_typ::init(Int4 r_per_io,Int4 r_per_do,Int4 x_ie, Int4 x_de, cma_typ cma,
				char dms_md, double prr_wt,ssx_typ *twin)
{
	if(PerNatsCMSA(cma) != 1000){
		fprintf(stderr,"PerNats(%.2f) reset to 1000\n",PerNatsCMSA(cma));
		SetPerNatsCMSA(1000,cma); SetPenaltyCMSA(10000,2000,cma);
	} pernats=1000; assert(x_de > 0);
	aa_per_io=r_per_io; aa_per_do=r_per_do; exp_ie=x_ie; exp_de=x_de; 
	dms_mode=dms_md;
	SetEx=0;
	// CMA=cma; 
	OPS=0;
	AddOp=0;
	MaxIter=INT4_MAX;
	AB=AlphabetCMSA(cma);
	NEWP(CsqCnts,nAlpha(AB) + 5, Int4);
	NEW(CsqQuery,NumSeqsCMSA(cma) + 5, e_type);
	// SaveBestCMSA(cma);	// save current as best
if(1)	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(SeqIsInCMSA(sq,cma)) ExtendFakeToRealCMSA(sq,cma);
	} 
// exit(1);
	SSX = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pernats,cma,dms_mode,twin);
	if(prr_wt > 0) SSX->SetPriorWt(prr_wt);
}

void    gmb_typ::Free( )
{
	cma_typ cma=SSX->RtnCMA();
	if(AddOp) free(AddOp);
	if(OPS) delete OPS;
	if(CsqQuery){
	   for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(CsqQuery[sq]) NilSeq(CsqQuery[sq]);
	   } free(CsqQuery);
	} free(CsqCnts); 
	delete SSX; 
}

void    gmb_typ::DeleteSSX(ssx_typ *ssx)
// WARNING: need to free cma AFTER deleting ssx!!!
{
	cma_typ cma=ssx->RtnCMA(); delete ssx; 
	gss_typ *gss=gssCMSA(cma); gss->~gss_typ(); NilCMSA(cma);
}


