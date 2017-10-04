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

void	gmb_typ::TestAddRemove()
{
	cma_typ	cma=SSX->RtnCMA();
	Int4    sq,i,r,r0,scr,min,Nadd,N=NumSeqsCMSA(cma);
	assert(nBlksCMSA(cma) == 1);
	char	**operation; NEWP(operation, N+4, char);
	Int4	**sites; NEWP(sites, N+4, Int4);
	Int4	score,*start; NEW(start, N+4, Int4);
#if 1
	for(sq=1; sq <= N; sq++){
	   assert(SeqIsInCMSA(sq,cma));
	   gss_typ *gss=gssCMSA(cma);
	   sites[sq]=GetPosSitesCMSA(sq,cma);
	   gsq_typ *gsq=gsqCMSA(sq,cma);
	   e_type sbjE=gss->TrueSeq(sq);
	   operation[sq]=SSX->GapAlnTrace(sbjE,0.0,start[sq],score);
	}
	for(sq=1; sq <= N; sq++){
	   Int4 *oldpos=SSX->RemoveFromAlign(sq);
	   assert(oldpos != 0);
	   assert(oldpos[1] == sites[sq][1]); free(oldpos);
	   assert(!SeqIsInCMSA(sq,cma));
	}
	fprintf(stderr,"LLR = %.3f\n",SSX->GapMap());
	for(sq=1; sq <= N; sq++){
	   assert(!SeqIsInCMSA(sq,cma));
	   gss_typ *gss=gssCMSA(cma);
	   gsq_typ *ogsq,*gsq = new gsq_typ[1]; 
	   e_type sbjE=gss->TrueSeq(sq);
	   Int4 *newpos; NEW(newpos,nBlksCMSA(cma)+3,Int4);
           gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                                operation[sq],strlen(operation[sq]),start[sq],sbjE,newpos);
	   ogsq=SwapGsqCMSA(sq,gsq,cma); SSX->AddToAlign(sq,newpos); delete []ogsq;
	   assert(SeqIsInCMSA(sq,cma));
	   free(newpos);
	   if(sites[sq]) free(sites[sq]); free(operation[sq]); 
	} free(sites); free(operation); free(start); 
#else
	for(sq=1; sq <= N; sq++){
	   gsq_typ *gsq=gsqCMSA(sq,cma);
	   sites[sq]=GetPosSitesCMSA(sq,cma);
	   operation[sq]=gsq->Operation(nBlksCMSA(cma),sites[sq],LengthsCMSA(cma));
	   Int4 *oldpos=SSX->RemoveFromAlign(sq);
	   assert(oldpos[1] == sites[sq][1]); free(oldpos);
	}
	fprintf(stderr,"LLR = %.3f\n",SSX->GapMap());
	for(sq=1; sq <= N; sq++){
	   SSX->AddToAlign(sq,sites[sq]); 
	   if(sites[sq]) free(sites[sq]); free(operation[sq]); 
	} free(sites); free(operation); free(start); 
#endif
	fprintf(stderr,"LLR = %.3f\n",SSX->GapMap());
}

void    gmb_typ::PutDebugClocl(Int4 strt,Int4 os,Int4 sq,char *op, e_type sbjE)
{
	fprintf(stderr,"%d: os = %d; strt=%d\n",sq,os,strt);
        if(strt < 1) this->PutGappedAln(stderr,sbjE,op,1);
	else this->PutGappedAln(stderr,sbjE,op,strt);
        fprintf(stderr,"op = %s\n",op); PutSeq(stderr,sbjE,AB);
}

void	gmb_typ::PutClusterIDs(FILE *fp,Int4 n, Int4 *cluster)
{
	fprintf(fp,"%d(%d): ",n,cluster[0]);
        for(Int4 x=1; cluster[x]; x++) fprintf(fp,"%d ",cluster[x]);
        fprintf(fp,"\n");
}

void    gmb_typ::PutCMA_SqIndels(FILE *fp, Int4 sq, Int4 *pos, gsq_typ *gsq)
{
	cma_typ cma=SSX->RtnCMA();
	Int4    *lens = LengthsCMSA(cma); 
	NEW(pos,nBlksCMSA(cma)+3,Int4);
	
	UInt4 nmm=0,nmd=0,nmi=0,nm=0,nii=0,ndd=0,nim=0,ndm=0,nid=0,ndi=0,nsd=0,nsm=0;
        gsq->FindIndels(nBlksCMSA(cma),pos,lens,nmm,nmi,nmd,nii,nim,ndd,ndm,ndi,nid,nsd,nsm);
	fprintf(fp,"sq=%d: nmm=%d; nmi=%d; nmd=%d; nii=%d; nim=%d\n",sq,nmm,nmi,nmd,nii,nim);
        fprintf(fp,"ndd=%d; ndm=%d; ndi=%d; nid=%d; nsd=%d; nsm=%d\n\n", ndd,ndm,ndi,nid,nsd,nsm);
}

void    gmb_typ::PutCMA_SqIndels(FILE *fp, Int4 sq)
{
	cma_typ cma=SSX->RtnCMA();
	gss_typ *gss=gssCMSA(cma);
	Int4    *pos,*lens = LengthsCMSA(cma); 
	NEW(pos,nBlksCMSA(cma)+3,Int4);
	
	if(PosSites(sq,pos,SitesCMSA(cma)) == 0) return;      // not in alignment; skip.
        assert(PosSites(sq,pos,SitesCMSA(cma)) == nBlksCMSA(cma));
	UInt4 nmm=0,nmd=0,nmi=0,nm=0,nii=0,ndd=0,nim=0,ndm=0,nid=0,ndi=0,nsd=0,nsm=0;
	gsq_typ *gsq=gss->GetGSQ(sq);
        gsq->FindIndels(nBlksCMSA(cma),pos,lens,nmm,nmi,nmd,nii,nim,ndd,ndm,ndi,nid,nsd,nsm);
	fprintf(fp,"sq=%d: nmm=%d; nmi=%d; nmd=%d; nii=%d; nim=%d\n",sq,nmm,nmi,nmd,nii,nim);
        fprintf(fp,"ndd=%d; ndm=%d; ndi=%d; nid=%d; nsd=%d; nsm=%d\n\n", ndd,ndm,ndi,nid,nsd,nsm);
	PutCMA_SqAln(fp,sq);
	free(pos);
}

void    gmb_typ::PutCMA_SqAln(FILE *fp, Int4 sq)
{
	cma_typ CMA=SSX->RtnCMA();
	BooLean skip[NumSeqsCMSA(CMA)+5];
	for(Int4 i=NumSeqsCMSA(CMA); i > 0; i--) skip[i]=TRUE; skip[sq]=FALSE;
	PutSelectCMSA(fp,skip,CMA);
}

Int4    gmb_typ::MaxSegSMatrix(unsigned char *maxseg, smx_typ M)
{
        Int4    s,i,r,K=M->K,max;

        for(i=1; i<=K; i++){
           for(max=INT4_MIN,r = 0; r <= M->nlet; r++){
                s = M->score[i][r];
                if(s > max){ maxseg[i]=r; max=s; }
           }
        } return i;
}

Int4	gmb_typ::PutGappedSeqAlnSMatrix(FILE *fp, char *operation, Int4 offset, Int4 n2,
        		unsigned char *seq2, Int4 nmod, smx_typ *M)
{
        Int4    o,m,i,j,k,r1,s,t,v,n1,*pos,p,p0,score;
	Int4	j2,s1,s0,g,gs,g_opt,J=strlen(operation);
	unsigned char	*seq1,*seq;
	char	*out[3];
	short	*mtf;
	a_type	A = SMatrixA(M[1]);
	Int4	fake_len=0,true_len=0;
	char	*BigStr=0;
	Int4	BigLen=0;

	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); NEW(mtf, n1+3, short); 
	for(s=0, m=1; m <= nmod; m++){
	    this->MaxSegSMatrix(seq1+s, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	} fake_len=s;
	/*** 3. Trace back step. ***/
	NEW(seq, J+3, unsigned char); 
	MEW(out[0],J+3,char); MEW(out[1],J+3,char); MEW(out[2],J+3,char); NEW(pos,J+3,Int4); 

	for(score=0,m=1,k=1,p=p0=1,o=j=i=1; operation[o] != 'E'; o++){ 
		switch(operation[o]){
		    case 'M': case 'm':
			v = ValSMatrix(k,seq2[j],M[m]); pos[p] = j; score+=v;
			out[1][p] = AlphaChar(seq1[i],A); out[2][p] = AlphaChar(seq2[j],A); 
			seq[p0]=seq2[j]; p0++;
			if(seq1[i]==seq2[j]) out[0][p] = ':';
			else if(v > 0) out[0][p] = '.'; else out[0][p] = ' ';
			p++; i++; j++; k++;
			break;
		    case 'D': 
		    case 'd': // Gap ('-') is in sequence
			score+=(Int4)floor(ExpScoreSMatrix(k,M[m])+0.5);
			out[1][p] = AlphaChar(seq1[i],A); out[0][p] = ' '; out[2][p] = '-';
			seq[p0]=0; p0++; pos[p]=j; p++; i++;  k++;
			break;
		    case 'I':   // Gap ('-') is in profile
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = tolower(AlphaChar(seq2[j],A)); pos[p] = j; p++; j++; break;
		    case 'i': j++; break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			assert(!"gmb_typ::PutGappedSeqAlnSMatrix(): this should not happen");
			print_error("gmb_typ::PutGappedSeqAlnSMatrix(): this should not happen"); 
			break;
		}
		if(k > LenSMatrix(M[m])){ // 4. Print out alignment 
		   if(m > 1)  fprintf(fp,"(%d)\n",pos[1]-pos[0]-1);
	   	   fprintf(fp,"\nsmx %c: ",(char)(m + 'A' - 1));
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[1][s]); 
	   	   fprintf(fp,"\n       ");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[0][s]); 
		   fprintf(fp,"\nseq  : ");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[2][s]);
	   	   fprintf(fp," %d-%d [%d]",pos[1]+offset,j-1+offset,score);
		   pos[0]=j-1;
		   if(m < nmod) { m++; k=1; p=p0=1; score=0; } else break;
		}
	} fprintf(fp,"(%d)\n",n2-pos[0]);
	/*** 5. free allocated memory ***/
	free(seq1); free(mtf); free(pos);  free(seq);
	free(out[0]); free(out[1]); free(out[2]);
	return o;
}

void    gmb_typ::PutGappedAln(FILE *fp,e_type sbjE,char *operation, Int4 start)
// This is slow because PutGappedSeqAlnSMatrix() calls smatrix_prob(); 
{
	cma_typ CMA=SSX->RtnCMA();
	PutSeqInfo2(fp,sbjE);
#if 0
        PutGappedSeqAlnSMatrix(fp,operation,start-1,
                LenSeq(sbjE)-start+1, SeqPtr(sbjE)+start-1, nBlksCMSA(CMA), SSX->RtnSMX());
#else
        this->PutGappedSeqAlnSMatrix(fp,operation,start-1,
                LenSeq(sbjE)-start+1, SeqPtr(sbjE)+start-1, nBlksCMSA(CMA), SSX->RtnSMX());
#endif
}

void	gmb_typ::PutOperations(FILE *fp,Int4 ii, Int4 s, char *operation, char *op)
{
	fprintf(fp,"%d(%d): \n",ii,s);
	fprintf(fp,"Operation(%d)=%s\n",strlen(operation),operation);
	fprintf(fp,"op(%d)=%s\n\n",strlen(op),op);
}

void	gmb_typ::PutMasterSlaveAln(Int4 s, e_type qE, Int4 start, char *qOp, 
				e_type E, Int4 strt, char *op, Int4 os, cma_typ cma)
// smx is globally defined within gmb_typ.
{
       Int4	trace_length=strlen(op);
       fprintf(stderr,"\n\n==== "); PutSeqInfo(stderr,E); fprintf(stderr," ====\n");
       fprintf(stderr,"\n%d(%d): Operation=%s\n",s,trace_length,qOp);
       PutGappedSeqAlnSMatrix(stderr,qOp,start-1, LenSeq(qE)-start+1, 
			SeqPtr(qE)+start-1, nBlksCMSA(cma),SSX->RtnSMX());
       fprintf(stderr,"\n  %d(%d): operation=%s\n",s,trace_length,op);
       PutGappedSeqAlnSMatrix(stderr,op,strt-1, LenSeq(E)-strt+1,
                SeqPtr(E)+strt-1, nBlksCMSA(cma),SSX->RtnSMX());
       fprintf(stderr,"\n  strt=%d; offset=%d\n",strt,os);
       // this->PutDiagonalSeq(stderr,os,qE,E); // DEBUG
       // exit(1);
}

void	gmb_typ::PutClusterAln(char *name,Int4 *cluster, cma_typ cma)
{
	Int4	ii,s;
	BooLean *skip; NEW(skip,NumSeqsCMSA(cma)+5,BooLean);
	for(ii=1; ii <= NumSeqsCMSA(cma); ii++) skip[ii]=TRUE;
	for(ii=1; (s=cluster[ii]) != 0; ii++){ skip[s]=FALSE; }
        FILE *tfp=open_file("junk_debug",".cma","w");
        PutSelectCMSA0(tfp,skip,cma);
        free(skip); fclose(tfp);
}

sap_typ	gmb_typ::MultiAlign(Int4 *ListE,Int4 *offset,e_type qE)
{
	sap_typ	sap=0,head=0,tail=0;
        Int4	i,j,os,x,s,sq=ListE[0],Qst,Qend,Sst,Send,len;
        e_type	sE;

	// 1. Align the rest of the sequences to the consensus 'cE'.
        for(j=0; j <= sq; j++){
		if(j==0){ os=0; sE=qE; }
		else { os=offset[j]; s=ListE[j]; sE=TrueSeqCMSA(s,SSX->RtnCMA()); }

		// this->PutDiagonalSeq(stderr,os,qE,sE); // DEBUG
		// fprintf(stderr,"os=%d\n",os);
		len=GetDiagonalEnds(os,qE,sE,Qst,Sst,Qend,Send);
		// WARNING: sap_typ start at zero while cma_typ starts at one!!!
        	char	*op; NEW(op,len+9,char); op[0]='E'; op[1]='M';
		for(x=2; x <= len; x++) op[x]='m'; op[x]='E';
		sap = ToGSeqAlign(strlen(op),op,qE,sE,Qst-1,Sst-1);
                if(head==0){ head=tail=sap; } else { tail->next=sap; tail=sap; }
		free(op);
	} 
	return head;
}

Int4	gmb_typ::PutMultiAlign(FILE *fp, Int4 *ListE,Int4 *offset,e_type qE)
// 2. Generate and print out the alignment.
{
        Int4	j,sq=ListE[0];
	
	sap_typ	head=this->MultiAlign(ListE,offset,qE);
        fprintf(fp,"cluster (%d seqs):",sq);
        for(j=1; j <= sq && j <= ListE[0]; j++) fprintf(fp," %d",ListE[j]);
        fprintf(stderr,"\n"); PutSeqInfo2(fp,qE); fprintf(fp,"\n");
        PutMultiGSeqAlign(fp,head,60,AB); FreeGSeqAlignList(head);
	return sq;
}

Int4	*gmb_typ::PutAlnToLongest(FILE *fp, Int4 **ListE, Int4 s)
{
	sap_typ	sap=0,head=0,tail=0;
        Int4	i,j,ssq,qsq,qst,sst,os,end,x,max,J,sq=ListE[s][0],Qst,Qend,Sst,Send,len;
        e_type	sE,qE;
	cma_typ	CMA=SSX->RtnCMA();

	// 1. find the longest sequence.
	for(max=0,J=0,j=1; j <= sq; j++){
		sE=TrueSeqCMSA(j,CMA);;
		if(LenSeq(sE) > max){ max = LenSeq(sE); J=j; }
	} if(J != 1){ qsq=ListE[s][J]; ListE[s][J]=ListE[s][1]; ListE[s][1]=qsq; }
	
	// 2. Align the rest of the sequences to the longest.
        qE=TrueSeqCMSA(qsq,CMA);;
	qst=TruePosCMSA(qsq,1,1,CMA); if(IsDeletedCMSA(1,qsq,1,CMA)) qst++;
        for(j=1; j <= sq; j++){
                ssq=ListE[s][j]; sE=TrueSeqCMSA(ssq,CMA);
		sst=TruePosCMSA(ssq,1,1,CMA); if(IsDeletedCMSA(1,ssq,1,CMA)) sst++;
		os=qst-sst;
		// if(fp) this->PutDiagonalSeq(fp,os,qE,sE); // DEBUG
		len=GetDiagonalEnds(os,qE,sE,Qst,Sst,Qend,Send);
		// if(fp) fprintf(fp,"%d: %d vs %d Qst=%d; Sst=%d; os=%d\n",s,qsq,ssq,Qst,Sst,os);
        	char	*op; NEW(op,len+9,char); op[0]='E'; op[1]='M';
		for(x=2; x <= len; x++) op[x]='m'; op[x]='E';
		// WARNING: sap_typ start at zero while cma_typ starts at one!!!
		if(fp){
		  sap = ToGSeqAlign(strlen(op),op,qE,sE,Qst-1,Sst-1);
                  if(head==0){ head=tail=sap; } else { tail->next=sap; tail=sap; }
		}
	} 

	// 3. Print out the alignment.
	if(fp){
          fprintf(fp,"cluster %d(%d seqs):",s,sq);
          for(j=1; j <= sq && j <=5; j++) fprintf(fp," %d",ListE[s][j]);
          fprintf(stderr,"\n"); PutSeqInfo2(fp,qE); fprintf(fp,"\n");
          PutMultiGSeqAlign(fp,head,60,AB); FreeGSeqAlignList(head);
	} 
}

Int4	gmb_typ::GetDiagonalEnds(Int4 os,e_type E1,e_type E2,
		Int4 &Strt1,Int4 &Strt2,Int4 &End1, Int4 &End2)
/*******************************************************************
 print the MSP for diagonal at offset 'os' between seq E1 and seq E2.
 *******************************************************************/
{
	Int4	n,n1,n2,score=0,start,end;
	unsigned char	*seq1,*seq2;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	if(os < (1-n2) || os >= n1) seq_error("offset out of range");
	seq1 = SeqPtr(E1); seq2 = SeqPtr(E2);
	if(os > 0){
	   n = MINIMUM(Int4,n1-os,n2);
	   score=get_diagonal_ends_seq(seq1+os,seq2,AlphaR(AB),n,&start,&end);
	   Strt1=start+os; End1=end+os; Strt2=start; End2=end;
	} else {
	   n = MINIMUM(Int4,n2+os,n1);
	   score=get_diagonal_ends_seq(seq2-os,seq1,AlphaR(AB),n,&start,&end);
	   Strt1=start; End1=end; Strt2=start-os; End2=end-os;
	} return end-start+1;
}

void	gmb_typ::PutListAsAln(FILE *fp, Int4 **ListE, Int4 *Strt, Int4 *Lngth, Int4 s)
{
	Int4	sq=ListE[s][0],i,j,ssq,qsq=ListE[s][1];
	sap_typ	sap=0,head=0,tail=0;
	cma_typ	CMA=SSX->RtnCMA();
        e_type	sE,qE=TrueSeqCMSA(qsq,CMA);;
        char	*operation; NEW(operation,Lngth[qsq]+9,char);
        operation[0]='E'; operation[1]='M';
        for(j=2; j < Lngth[qsq]; j++) operation[j]='m';
        if(Strt[qsq] > 0){ operation[j]='m'; j++; } operation[j]='E';
        Int4 qst=0; // Strt[qsq];
        // WARNING: sap_typ code treats deleted starting positions differently than qsq_typ.
	qst=TruePosCMSA(qsq,1,1,CMA); if(!IsDeletedCMSA(1,qsq,1,CMA)) qst--;
	gss_typ *gss=gssCMSA(CMA);
	// gss->Put(stderr,qsq);
	// char *qop=gss->Operation(stderr,qsq); fprintf(stderr,"%d: %s\n",qsq,qop); 
#if 0
char tmpstr[5000];
Int4	m;
for(m=1; m <= nBlksCMSA(CMA); m++){
	LenStr=gss->Region(gsq,tmpstr,SitePos(m,J,1,S),lenM);
}
#endif
        for(j=1; j <= sq; j++){
                ssq=ListE[s][j]; sE=TrueSeqCMSA(ssq,CMA);
                // Int4 sst=Strt[ssq]; if(sst > 0) sst=qst-1;
		Int4 sst=TruePosCMSA(ssq,1,1,CMA); if(!IsDeletedCMSA(1,ssq,1,CMA)) sst--;
#if 1
                sap = ToGSeqAlign(strlen(operation),operation,qE,sE,qst,sst);
                if(head==0){ head=tail=sap; } else { tail->next=sap; tail=sap; }
#else
	    char *op=gss->Operation(ssq); fprintf(stderr,"%d: %s\n",ssq,op); 
	    if(strcmp(qop,op)==0){	// aligned identically...
                sap = ToGSeqAlign(strlen(op),op,qE,sE,qst,sst);
                if(head==0){ head=tail=sap; } else { tail->next=sap; tail=sap; }
	    } free(op);
#endif
        } // free(qop);
        fprintf(stderr,"cluster %d(%d seqs)(query start=%d):",s,sq,Strt[qsq]);
        for(j=1; j <= sq && j <=5; j++) fprintf(stderr," %d",ListE[s][j]);
        fprintf(stderr,"\n"); PutSeqInfo2(stderr,qE); fprintf(stderr,"\n");
        PutMultiGSeqAlign(stderr,head,60,AB); FreeGSeqAlignList(head);
	free(operation);
}

Int4	gmb_typ::PutDiagonalSeq(FILE *fptr, Int4 offset, e_type E1, e_type E2)
/*******************************************************************
 print the MSP for diagonal at offset between seq1 and seq2.

 *******************************************************************/
{
	Int4	v,i,j,n,n1,n2,score=0,b,e,start,end,w=40;
	unsigned char	*seq1,*seq2;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	v = offset;
	if(v < (1-n2) || v >= n1) seq_error("offset out of range");
	seq1 = SeqPtr(E1); seq2 = SeqPtr(E2);
	if(v > 0){
	   n = MINIMUM(Int4,n1-v,n2);
	   score=get_diagonal_ends_seq(seq1+v,seq2,AlphaR(AB),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = MINIMUM(Int4,e,b+w-1);
		fprintf(fptr,"%4d ",v+b+OffSetSeq(E1));
		for(i=v+b,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c",AlphaChar(seq1[i],AB));
	  	fprintf(fptr," %4d ",v+end+OffSetSeq(E1));
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=v+b,j=b; j <= end; i++,j++)
			if(seq1[i]==seq2[j]) fprintf(fptr,":");
			else if(valAlphaR(seq1[i],seq2[j],AB) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4d ",b+OffSetSeq(E2));
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c", AlphaChar(seq2[j],AB));
		fprintf(fptr," %4d ",end+OffSetSeq(E2));
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %d\n\n",score);
	} else {
	   n = MINIMUM(Int4,n2+v,n1);
	   score=get_diagonal_ends_seq(seq2-v,seq1,AlphaR(AB),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = MINIMUM(Int4,e,b+w-1);
		fprintf(fptr,"%4d ",b+OffSetSeq(E1));
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c",AlphaChar(seq1[j],AB));
		fprintf(fptr," %4d ",end+OffSetSeq(E1));
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=b-v,j=b; j <= end; i++,j++)
			if(seq2[i]==seq1[j]) fprintf(fptr,":");
			else if(valAlphaR(seq2[i],seq1[j],AB) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4d ",b-v+OffSetSeq(E2));
		for(i=b-v,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c", AlphaChar(seq2[i],AB));
	  	fprintf(fptr," %4d ",end-v+OffSetSeq(E2));
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %d\n\n",score);
	}
	return score;
}

