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

#include "scaninfo.h"

sni_typ	ReMakeScanInfo(sni_typ I,float score,Int4 nblk, Int4 rpts, 
	float *blkscore, Int4 *site)
{
	Int4	i;

	I->S = score;
	if(nblk != I->nblk) {
		I->rpts = rpts;
		I->nblk = nblk;
		free(I->s); free(I->site);
		NEW(I->s,nblk+2,float);
		NEW(I->site,nblk+2,Int4);
	}
	for(i=1; i<=nblk; i++){
		I->s[i] = blkscore[i]; I->site[i] = site[i];
	}
	return I;
}

sni_typ	MakeRptScanInfo(Int4 r, Int4 Nt_flank, Int4 Ct_flank, sni_typ I)
/*** create a sni_typ data structure for the rth repeat ***/
{
	sni_typ	I2;
	Int4	i,j,n,start,end,offset;

	if(r < 1 || r > I->rpts) print_error("MakeRptScanInfo( ) input error");
	NEW(I2,1,scaninfo_type);
	I2->S = I->S;
	I2->method = I->method;
	I2->id = I->id;
	I2->rpts = 1;
	I2->next = NULL;
	I2->nblk = n = I->nblk/I->rpts; 
	NEW(I2->s,n+2,float);
	NEW(I2->site,n+2,Int4);
	start = 1 + (r-1)*n;
	Nt_flank = MINIMUM(Int4, Nt_flank, I->site[start] - 1);
	offset = I->site[start] - Nt_flank - 1;
/*****
fprintf(stderr,"offset = %d; Ntflank = %d; I->site[%d] = %d\n",
		offset , Nt_flank, start, I->site[start]);
/*****/
	for(j=start, i=1; i<=n; i++,j++){
		I2->site[i] = I->site[j] - offset;
		I2->s[i] = I->s[j]; 
	}
	start = I2->site[1] - Nt_flank;
	end = MINIMUM(Int4, I->site[j-1] + Ct_flank, LenSeq(I->E));
	I2->E = MkSubSeq(start,end,I->E);
	return I2;
}

sni_typ	MakeScanInfo(e_type E,float score,Int4 nblk, Int4 rpts,
	float *blkscore, Int4 *site, Int4 seqid,char method)
{
	sni_typ	I;
	Int4	i;

	NEW(I,1,scaninfo_type);
	/*** I->E = CopySeq(E); /** FOR TESTING **/
	/**/ I->E = E; /****/
	I->S = score;
	I->method = method;
	I->rpts = rpts;
	I->id = seqid;
	I->nblk = nblk;
	I->next = NULL;
	NEW(I->s,nblk+2,float);
	NEW(I->site,nblk+2,Int4);
	for(i=1; i<=nblk; i++){
		I->s[i] = blkscore[i]; I->site[i] = site[i];
	}
	return I;
}

sni_typ	AppendScanInfo(sni_typ I, sni_typ List)
{
	sni_typ	last;

	if(List == NULL) return I;
	for(last=List; last->next != NULL; ) last = last->next;
	last->next = I;
	return List;
}

void	PutScanInfo(FILE *fptr, sni_typ I, Int4 *len, BooLean msaformat, a_type A)
{
	e_type	E;
	Int4	s,i,m,m0,N0,n,nblks,gap,s2,rpts,gapped;
	float	score,p;
	char	method;

	method = MethodScanInfo(I);
	score = ScoreScanInfo(I);
	E = SeqScanInfo(I);
	if(!msaformat) fprintf(fptr,"[%3.2f] ", -score);
	else fprintf(fptr,"     ");
	PutSeqInfo2(fptr,E);

	N0=nblks = nBlksScanInfo(I);
	rpts = RptsScanInfo(I);
	if(rpts > 1) fprintf(fptr,"  (%d repeats)\n",rpts);

	for(m=1; m<=nblks; m++){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
	        if(!msaformat) fprintf(fptr," %4.1f: ",p);
		m0 = ((m-1) % N0) + 1;
		if(msaformat) {
		   PutSeqRegionFormatMSA(fptr,s, len[m0], E, p, A);
		} else {
		   PutSeqRegion2(fptr,s,len[m0],E,0,A);
		   /***** put insert lengths here *****/
		  if(m < nblks){
		   s2=SiteScanInfo(m+1,I);
		   gap = s2 - (s + len[m0]);
		   fprintf(fptr," (%d)",gap);	
		  } else {
		   s2=LenSeq(E)+CtermExtendSeq(E)+1;
		   gap = s2 - (s + len[m0]);
		   if(gap < 0) gap = 0;
		   fprintf(fptr," (%d)",gap);	
		  } 
		  fprintf(fptr,"\n");
		} 
		if((m%N0) == 0 && m < nblks) fprintf(fptr,"\n");
	}
        if(islower(method)) {
                fprintf(fptr,"  (gapped E-value used)");
        }
	fprintf(fptr,"\n");
}

sni_typ	RmLastScanInfo(sni_typ List)
{
	sni_typ	last,nxt;

	if(List == NULL) return NULL;
	if(List->next == NULL) return List;
	for(nxt=List; nxt->next->next != NULL; ) { nxt = nxt->next; }
	last = nxt->next; nxt->next = NULL;
	return last;
}

Int4	LengScanInfo(sni_typ List)
{
	sni_typ	last;
	Int4	n;

	if(List == NULL) return 0;
	for(n=1,last=List; last->next != NULL; n++) last = last->next;
	return n;
}

e_type	NilScanInfoRtnE(sni_typ I)
{
	Int4	n;
	e_type	E;
	
	n = LengScanInfo(I);
	if(n > 1) print_error("RtnNilScanInfo( ): repeats not yet implemented");
	else if(I == NULL) print_error("RtnNilScanInfo( ) input error");
	E = I->E;
	free(I->s); free(I->site); free(I);
	return E;
}

void    NilScanInfo(sni_typ I)
{
	sni_typ	tmp,nxt=I;
	
	while(nxt!= NULL){
		tmp = nxt->next; 
		NilSeq(nxt->E); free(nxt->s); free(nxt->site); free(nxt);
		nxt = tmp; 
	}
}



