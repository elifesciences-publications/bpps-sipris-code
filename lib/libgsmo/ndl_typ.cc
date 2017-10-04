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

#include "ndl_typ.h"

void    ndl_typ::init(Int4 aa_per_io,Int4 aa_per_do,Int4 exp_ie, Int4 exp_de,
						double pnts, cma_typ cma,char wf)
{
	pernats=pnts;
	assert(pernats==1000 || pernats==693 || pernats == 1443);
	aa_per_ins_o=aa_per_io;
	aa_per_del_o=aa_per_do;
	exp_ins_ext=exp_ie;
	exp_del_ext=exp_de;
	// prior_wt=0.50;
	prior_wt=1.0;
	CMA=cma;
	WtFact=wf;
	initTransProb( );
	initTransitions( );
	ReCalc=TRUE;
}

void	ndl_typ::FreeTransProb( )
{
	free(m2m); free(m2i); free(m2d);
	free(i2i); free(i2m); free(i2d);
	free(d2d); free(d2m); free(d2i);
	free(s2d); free(s2m); 
}
void	ndl_typ::initTransProb( )
{
	Int4 len=TotalLenCMSA(CMA);
	NEW(m2m,len+9,Int4); NEW(m2i,len+9,Int4); NEW(m2d,len+9,Int4);
	NEW(i2m,len+9,Int4); NEW(i2i,len+9,Int4); NEW(i2d,len+9,Int4);
	NEW(d2m,len+9,Int4); NEW(d2d,len+9,Int4); NEW(d2i,len+9,Int4);
	NEW(s2m,len+9,Int4); NEW(s2d,len+9,Int4);
}

void	ndl_typ::initTransitions( )
{
	Int4 len=TotalLenCMSA(CMA);
	NEW(Nmm,len+9,UInt4); NEW(Nmi,len+9,UInt4); NEW(Nmd,len+9,UInt4);
	NEW(Nii,len+9,UInt4); NEW(Nim,len+9,UInt4); NEW(Nid,len+9,UInt4);
	NEW(Ndd,len+9,UInt4); NEW(Ndm,len+9,UInt4); NEW(Ndi,len+9,UInt4);
	NEW(Nsd,len+9,UInt4); NEW(Nsm,len+9,UInt4); // NEW(Nsi,len+9,UInt4);
}

void	ndl_typ::FreeTransitions( )
{
	free(Nmm); free(Nmi); free(Nmd);
	free(Nii); free(Nim); free(Nid);
	free(Ndd); free(Ndm); free(Ndi);
	free(Nsd); free(Nsm); // free(Nsi);
}

void	ndl_typ::RmIndels(Int4 sq,char Wt)
{
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	gss_typ	*gss=gssCMSA(CMA); // assert(gss->Gapped());

	assert(sq > 0 && sq <= NumSeqsCMSA(CMA));
        NEW(pos,nBlksCMSA(CMA)+3,Int4);
	nblks=PosSites(sq,pos,SitesCMSA(CMA));
	if(nblks == 0){ free(pos);  return; } // sequence not in alignment; skip.
	assert(nblks == nBlksCMSA(CMA));
	gsq_typ *gsq=gss->GetGSQ(sq);
	assert(Wt != 0);
	gsq->AddTransitions(nblks,pos,lens,Nmm,Nmi,Nmd,Nii,Nim,Ndd,Ndm,Ndi,Nid,Nsd,Nsm,-Wt);
		// Nsd-=nsd*Wt[sq]; Nsm-=nsm*Wt[sq];
	ReCalc=TRUE; free(pos);
}

void	ndl_typ::AddIndels(Int4 sq,char Wt)
{
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	gss_typ	*gss=gssCMSA(CMA); // assert(gss->Gapped());
	
	assert(sq > 0 && sq <= NumSeqsCMSA(CMA));
        NEW(pos,nBlksCMSA(CMA)+3,Int4);
	nblks=PosSites(sq,pos,SitesCMSA(CMA));
	if(nblks == 0){ free(pos);  return; } // sequence not in alignment; skip.
	assert(nblks == nBlksCMSA(CMA));
	gsq_typ *gsq=gss->GetGSQ(sq);
	assert(Wt != 0);
	gsq->AddTransitions(nblks,pos,lens,Nmm,Nmi,Nmd,Nii,Nim,Ndd,Ndm,Ndi,Nid,Nsd,Nsm,Wt);
	// Nsd+=nsd*Wt[sq]; Nsm+=nsm*Wt[sq];
	ReCalc=TRUE; free(pos);
}

#define debug_NDL_penalty 0

void	ndl_typ::InitIndels(char *Wt)
// Jun Liu's statistical model for GISMO...
{
	//========== 1. Get # observed indel and match states. ==============
	this->FreeTransitions( );
	this->initTransitions();
	Int4	nblks,*pos,*lens = LengthsCMSA(CMA);
	NEW(pos,nBlksCMSA(CMA)+3,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(CMA); sq++){
		nblks=PosSites(sq,pos,SitesCMSA(CMA));
		if(nblks == 0) continue;	// sequence not in alignment; skip.
	   	assert(nblks == nBlksCMSA(CMA));
		this->AddIndels(sq,Wt[sq]);
	} free(pos);
#if debug_NDL_penalty
        // fprintf(stderr,"Nmm=%d; Nmi=%d; Nmd=%d; Nii=%d; Nim=%d\n",Nmm,Nmi,Nmd,Nii,Nim);
        // fprintf(stderr,"  Ndd=%d; Ndm=%d; Ndi=%d; Nid=%d; Nsd=%d; Nsm=%d\n",Ndd,Ndm,Ndi,Nid,Nsd,Nsm);
#endif
}

void	gsq_typ::AddTransitions(Int4 nblks, Int4 *start, Int4 *len, 
	UInt4 *Nmm, UInt4 *Nmi, UInt4 *Nmd, UInt4 *Nii,
	UInt4 *Nim, UInt4 *Ndd, UInt4 *Ndm,
	UInt4 *Ndi,UInt4 *Nid,UInt4 *Nsd, UInt4 *Nsm, char SqWt)
// Adds the weighted observed numbers of transitions to passed in variables.
{
   assert(fakeE != NULL); 
   Int4	n,end,bk,i,j,J;
   char	state;
   for(J=0,bk = 1 ; bk <= nblks; bk++){
	end=start[bk] + len[bk] -1;
	assert(end <= LenSeq(fakeE));
   	state='S';
	for(i=start[bk]; i <= end; i++){
	   J++;	// position in arrays above...
           if(!deletion(i)){	// next state == 'M';
		switch(state){
		  case 'S': Nsm[J] += SqWt; break;
		  case 'M': Nmm[J] += SqWt; break;
		  case 'I': Nim[J] += SqWt; break;
		  case 'D': Ndm[J] += SqWt; break;
		  default: 
		     fprintf(stderr,"state = %c; sq=%s\n",state,SeqKey(realE));
		     print_error("gsq_typ::FindIndels 'M' input error"); 
		  break;
		} state='M';
	   } else {		// next state == 'D'
		  switch(state){
		   case 'S': Nsd[J] += SqWt; break;
		   case 'M': Nmd[J] += SqWt; break;
		   case 'D': Ndd[J] += SqWt; break;
		   case 'I': Nid[J] += SqWt; break;	// I->D transitions sometimes show up.
		   default: 
		     fprintf(stderr,"state = %c; sq=%s\n",state,SeqKey(realE));
		     print_error("gsq_typ::FindIndels 'D' input error"); 
		   break;
		  } state='D';
	   }
	   if(i != end && (n=insert(i))){	// is there a following insertion?
	   // Don't worry about insertions at ends...
		switch(state){
		   // e.g., "TaFfsv..." = previous insert followed by a match and more inserts.
		   // case 'I' should not occur...
		   // case 'S': break; // insert or match...
		   case 'M': Nmi[J] += SqWt; Nii[J]+=(n-1)*SqWt; break;
		   case 'D': Ndi[J] += SqWt; Nii[J]+=(n-1)*SqWt; break;
		   default: 
		     fprintf(stderr,"state = %c; n=%d; bk = %d; i = %d; sq=%s\n",
			state,n,bk,i,SeqKey(realE));
		     for(Int4 j=start[bk]; j <= end; j++){
			 fprintf(stderr,"indels[%d]=%d\n",j,indels[j]);
		     }
		     print_error("gsq_typ::FindIndels 'I' input error"); 
		   break;
		} state='I';
	   }
	}
   } return;
}

