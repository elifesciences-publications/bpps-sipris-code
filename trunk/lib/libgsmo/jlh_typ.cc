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

#include "jlh_typ.h"

void    jlh_typ::init(Int4 aa_per_io,Int4 aa_per_do,Int4 exp_ie, Int4 exp_de,
						double pnts, cma_typ cma,char wf)
{
	pernats=pnts;
	assert(pernats==1000);
	aa_per_ins_o=aa_per_io;
	aa_per_del_o=aa_per_do;
	exp_ins_ext=exp_ie;
	exp_del_ext=exp_de;
	prior_wt=1;
	CMA=cma;
	Nmm=Nmi=Nmd=Nii=Nim=Ndd=Ndm=Nid=Ndi=Nsd=Nsm=0;
	WtFact=wf;
	ReCalc=TRUE;
}

void	jlh_typ::RmIndels(Int4 sq,char Wt)
{
	UInt4 nmm=0,nmd=0,nmi=0,nm=0,nii=0,ndd=0,nim=0,ndm=0,nid=0,ndi=0,nsd=0,nsm=0;
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	gss_typ	*gss=gssCMSA(CMA); // assert(gss->Gapped());

	assert(sq > 0 && sq <= NumSeqsCMSA(CMA));
        NEW(pos,nBlksCMSA(CMA)+3,Int4);
	nblks=PosSites(sq,pos,SitesCMSA(CMA));
	if(nblks == 0){ free(pos);  return; } // sequence not in alignment; skip.
	assert(nblks == nBlksCMSA(CMA));
	gsq_typ *gsq=gss->GetGSQ(sq);
	gsq->FindIndels(nblks,pos,lens,nmm,nmi,nmd,nii,nim,ndd,ndm,ndi,nid,nsd,nsm);
	if(Wt != 0){
		Nmm-=nmm*Wt; Nmi-=nmi*Wt; Nmd-=nmd*Wt; Nii-=nii*Wt; Nim-=nim*Wt; Nid-=nid*Wt; 
		Ndd-=ndd*Wt; Ndm-=ndm*Wt; Ndi-=ndi*Wt; // Nsd-=nsd*Wt[sq]; Nsm-=nsm*Wt[sq];
	} else {
		Nmm-=nmm; Nmi-=nmi; Nmd-=nmd; Nii-=nii; Nim-=nim; 
		Nid-=nid; Ndd-=ndd; Ndm-=ndm; Ndi-=ndi;
	} ReCalc=TRUE; free(pos);
}

void	jlh_typ::AddIndels(Int4 sq,char Wt)
{
	UInt4	nmm=0,nmd=0,nmi=0,nm=0,nii=0,ndd=0,nim=0,ndm=0,nid=0,ndi=0,nsd=0,nsm=0;
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	gss_typ	*gss=gssCMSA(CMA); // assert(gss->Gapped());
	
	assert(sq > 0 && sq <= NumSeqsCMSA(CMA));
        NEW(pos,nBlksCMSA(CMA)+3,Int4);
	nblks=PosSites(sq,pos,SitesCMSA(CMA));
	if(nblks == 0){ free(pos);  return; } // sequence not in alignment; skip.
	assert(nblks == nBlksCMSA(CMA));
	gsq_typ *gsq=gss->GetGSQ(sq);
	gsq->FindIndels(nblks,pos,lens,nmm,nmi,nmd,nii,nim,ndd,ndm,ndi,nid,nsd,nsm);
	if(Wt != 0){
		Nmm+=nmm*Wt; Nmi+=nmi*Wt; Nmd+=nmd*Wt; Nii+=nii*Wt; Nim+=nim*Wt; Nid+=nid*Wt; 
		Ndd+=ndd*Wt; Ndm+=ndm*Wt; Ndi+=ndi*Wt; // Nsd+=nsd*Wt[sq]; Nsm+=nsm*Wt[sq];
	} else {
		Nmm+=nmm; Nmi+=nmi; Nmd+=nmd; Nii+=nii; Nim+=nim; 
		Nid+=nid; Ndd+=ndd; Ndm+=ndm; Ndi+=ndi;
	} ReCalc=TRUE; free(pos);
}

#define debug_JLH_penalty 0

void	jlh_typ::InitIndels(char *Wt)
// Jun Liu's statistical model for GISMO...
{
	gss_typ	*gss=gssCMSA(CMA); // assert(gss->Gapped());
	
	//========== 1. Get # observed indel and match states. ==============
	Nmm=Nmi=Nmd=Nii=Nim=Ndd=Ndm=Nid=Ndi=Nsd=Nsm=0;
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	NEW(pos,nBlksCMSA(CMA)+3,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(CMA); sq++){
		nblks=PosSites(sq,pos,SitesCMSA(CMA));
		if(nblks == 0) continue;	// sequence not in alignment; skip.
	   	assert(nblks == nBlksCMSA(CMA));
		gsq_typ *gsq=gss->GetGSQ(sq);
		UInt4 nmm=0,nmd=0,nmi=0,nm=0,nii=0,ndd=0,nim=0,ndm=0,nid=0,ndi=0,nsd=0,nsm=0;
		gsq->FindIndels(nblks,pos,lens,nmm,nmi,nmd,nii,nim,ndd,ndm,ndi,nid,nsd,nsm);
		if(Wt != 0){
			Nmm+=nmm*Wt[sq]; Nmi+=nmi*Wt[sq]; Nmd+=nmd*Wt[sq]; 
			Nii+=nii*Wt[sq]; Nim+=nim*Wt[sq]; Nid+=nid*Wt[sq]; 
			Ndd+=ndd*Wt[sq]; Ndm+=ndm*Wt[sq]; Ndi+=ndi*Wt[sq];
			// Nmd+=nsd*Wt[sq]; 
			// Nsd+=nsd*Wt[sq]; 
			// Nsm+=nsm*Wt[sq];
		} else {
			Nmm+=nmm; Nmi+=nmi; Nmd+=nmd; Nii+=nii; Nim+=nim; Nid+=nid; 
			Ndd+=ndd; Ndm+=ndm; Ndi+=ndi;
		}
#if 0	// DEBUG!!!
		if(sq == 469){
			fprintf(stderr,"Wt[sq]=%d; nmm=%d; nmi=%d; nmd=%d; nii=%d; nim=%d\n",
				Wt[sq],nmm,nmi,nmd,nii,nim);
			fprintf(stderr,"ndd=%d; ndm=%d; ndi=%d; nid=%d; nsd=%d; nsm=%d\n",
				ndd,ndm,ndi,nid,nsd,nsm);
			fprintf(stderr," Nmm=%d; Nmi=%d; Nmd=%d; Nii=%d; Nim=%d\n",
				nmm*Wt[sq],nmi*Wt[sq],nmd*Wt[sq],nii*Wt[sq],nim*Wt[sq]);
			fprintf(stderr,"Ndd=%d; Ndm=%d; Ndi=%d; Nid=%d; Nsd=%d; Nsm=%d\n",
				ndd*Wt[sq],ndm*Wt[sq],ndi*Wt[sq],nid*Wt[sq],nsd*Wt[sq],nsm*Wt[sq]);
		}
#endif
	} free(pos);
#if debug_JLH_penalty
        fprintf(stderr,"Nmm=%d; Nmi=%d; Nmd=%d; Nii=%d; Nim=%d\n",Nmm,Nmi,Nmd,Nii,Nim);
        fprintf(stderr,"  Ndd=%d; Ndm=%d; Ndi=%d; Nid=%d; Nsd=%d; Nsm=%d\n",Ndd,Ndm,Ndi,Nid,Nsd,Nsm);
#endif
}


#if 0
if(0){
        printf("m2m=%d; m2i=%d; m2d=%d; i2i=%d; i2m=%d; d2d=%d; d2m=%d; s2m=%d.\n",
                        m2m,m2i,m2d,i2i,i2m,d2d,d2m,s2m);
        printf("map = %g; penalty = %g; total map = %g\n",oldmap,penalty,oldmap+penalty);
}
#endif



