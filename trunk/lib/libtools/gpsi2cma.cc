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

#include "gpsi2cma.h"

Int4    StartSubjectGSAP(sap_typ sap) { return sap->segs->starts[1]+1; }
// adds one to blast position.

Int4    EndSubjectGSAP(sap_typ sap)
// adds one to blast position.
{
	dsp_typ dsp = sap->segs; 
	return dsp->starts[2*(dsp->numseg-1)+1]+dsp->lens[dsp->numseg-1]; 
}

Int4    NumHSPsListGSAP(sap_typ head, Int4 max_sq_id, Int4 *min_qs, Int4 *max_qe)
{ return NumHSPsListGSAP(head,0,max_sq_id,min_qs,max_qe); }

Int4    NumHSPsListGSAP(sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
	Int4 *min_qs, Int4 *max_qe)
// return the number of aligned segment pairs in list...
// min_qs = min
// max_qs = 
{
	dsp_typ	dsp;
        UInt4   n=0,id;
	dsp = head->segs; e_type qE = dsp->qE;
	Int4		qs,qe,min=LenSeq(qE),max=0;

        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  id = SeqI(SubjectSeqGSAP(gsap));
	  if(id < min_sq_id || id > max_sq_id) continue; // NEW...
          dsp=gsap->segs;
	  qs = dsp->starts[0]; 
	  if(min > qs) min=qs;
	  qe = dsp->starts[2*(dsp->numseg-1)]+dsp->lens[(dsp->numseg-1)]-1; 
	  if(max < qe) max=qe;
          n++;
        } *min_qs=min; *max_qe=max;
	return n;
}

Int4    NumSeqsListGSAP(sap_typ head, Int4 max_sq_id,Int4 *min_qs, Int4 *max_qe)
{ return NumSeqsListGSAP(head,0,max_sq_id,min_qs,max_qe); }

Int4    NumSeqsListGSAP(sap_typ head, Int4 min_sq_id, Int4 max_sq_id,
	Int4 *min_qs, Int4 *max_qe)
// return the number of sequences in list...
// min_qs = min
// max_qs = 
{
	dsp_typ	dsp;
        UInt4   n=0,last_id=UINT4_MAX,id;
	dsp = head->segs; e_type qE = dsp->qE;
	Int4		qs,qe,min=LenSeq(qE),max=0;

        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  id = SeqI(SubjectSeqGSAP(gsap));
	  if(id < min_sq_id || id > max_sq_id) continue; 
          dsp=gsap->segs;
          if(dsp->subject_id != last_id){
		last_id= dsp->subject_id;
		qs = dsp->starts[0]; 
		if(min > qs) min=qs;
		qe = dsp->starts[2*(dsp->numseg-1)]+dsp->lens[(dsp->numseg-1)]-1; 
		if(max < qe) max=qe;
                n++;
          }
        } *min_qs=min; *max_qe=max;
	return n;
}

Int4	InDelsSeqAlign(sap_typ sap, Int4 *ins, Int4 *del,Int4 min_qs, Int4 max_qe)
{
        Int4		s,q,qs,ss,length,index,total,m,I=0,D=0,M=0;
	dsp_typ	dsp;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
        for(total=index=0; index<dsp->numseg; index++){ total += dsp->lens[index]; }
	{ // Deletions at N-terminal end...
		qs = dsp->starts[0];
		for(q=min_qs; q < qs; q++) D++;
	}
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1 && ss != -1){	// match region...
		    M+=length;
		} else if(qs == -1){		// insertion in subject sequence.
		    I+=length;
		} else if(ss == -1){		// deletion in subject sequence.
		    D+=length;
		}
	} 
	{ // Deletions at C-terminal end...
		Int4 qe = dsp->starts[2*(dsp->numseg-1)]+length;
		for(q=qe; q<=max_qe; q++) D++;
	}
	*ins=I; *del=D;
	return M;
}

void    tPutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank, 
	Int4 min_qs, Int4 max_qe, a_type AB)
// print output to a file in cma format...
{
        Int4		r,i,s,q,qs,ss,length,index,total;
	UInt4	qid,sid;
	dsp_typ	dsp;
	e_type		qE,sE;
	Int4		qos,sos,m;
	char		*qstr,*sstr,c;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
	fprintf(fp,">"); PutSeqInfo(fp,dsp->sE); 
	qid = dsp->query_id; qE = dsp->qE;
	sid = dsp->subject_id; sE = dsp->sE;
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);
        for(total=index=0; index<dsp->numseg; index++){
		total += dsp->lens[index];
	}
	ss = dsp->starts[1];
	// fprintf(fp,"{(");
	for(r=1; r <= ss; r++){
		c = AlphaChar(ResSeq(r,sE),AB); fprintf(fp,"%c",c);
	}
	fprintf(fp,"{()");
	// fprintf(fp,")");
	{
		qs = dsp->starts[0];
		for(q=min_qs; q < qs; q++) fprintf(fp,"-");
	}
	NEW(qstr,total+3,char); NEW(sstr,total+3,char);
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1 && ss != -1){	// match region...
		    for(i=0,r=ss+1; i < length; r++,i++){
			c= AlphaChar(ResSeq(r,sE),AB);
			fprintf(fp,"%c",c);
		    } qs++;
		} else if(qs == -1){		// insertion in subject sequence.
		    if(index==dsp->numseg-1) fprintf(fp,"(");
		    for(i=0,r=ss+1; i < length; r++,i++){
			c= AlphaChar(ResSeq(r,sE),AB);
			c = tolower(c); fprintf(fp,"%c",c);
		    } if(index==0) fprintf(fp,")");
		} else if(ss == -1){		// deletion in subject sequence.
			for(i=0; i < length; i++) fprintf(fp,"-"); 
		} else print_error("PutSeqAlignToCMA( ): dsp error");
		// if(index==dsp->numseg-1 && qs != -1) fprintf(fp,"(");
        }
	{
		Int4 qe = dsp->starts[2*(dsp->numseg-1)]+length;
		for(q=qe; q<=max_qe; q++) fprintf(fp,"-");
	}
	fprintf(fp,"()}");
	// fprintf(fp,"(");
	for( ; r <= LenSeq(sE); r++){
		c= AlphaChar(ResSeq(r,sE),AB);
		fprintf(fp,"%c",c);
	}
	// fprintf(fp,")}");
	// fprintf(fp,"}*");
	fprintf(fp,"*");
	free(qstr); free(sstr); 
}

void    QueryFirstSeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank,
	Int4 rightflank, a_type AB)
{ QueryFirstSeqAlignToCMA(fp,head,UINT4_MAX,leftflank,rightflank,AB); }

void    QueryFirstSeqAlignToCMA(FILE *fp, sap_typ head, Int4 max_sq_id,
	Int4 leftflank, Int4 rightflank, a_type AB)
{
	if(head==0) return;
	UInt4	last_id=UINT4_MAX;
	Int4		i,n,ins,del,mat,min_qs,max_qe,round;
	sap_typ 	gsap;

	n = NumSeqsListGSAP(head, max_sq_id, &min_qs, &max_qe);
	fprintf(fp,"[0_(1)=psiblast(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",n);
	fprintf(fp,"(%d)",max_qe-min_qs+1);
	for(i=0; i < (max_qe-min_qs+1); i++) fprintf(fp,"*");
	fprintf(fp,"\n"); i=0;
	dsp_typ dsp;
	for(round = 1; round <= 2; round++){
	  for(gsap=head; gsap!=NULL; gsap=gsap->next){
	   if(SeqI(SubjectSeqGSAP(gsap)) > max_sq_id) continue; // NEW...
	   if(round==1 && IsQueryGSAP(gsap) || round==2 && !IsQueryGSAP(gsap)){
	     dsp=gsap->segs; 
	     if(dsp->subject_id != last_id){
		last_id= dsp->subject_id; i++;
		mat = InDelsSeqAlign(gsap,&ins,&del,min_qs,max_qe);
		Int4 len = LenSeq(dsp->sE);
#if 1	// Use actual sequence identifier.
		fprintf(fp,"\n$%d=%d(%d):\n",SeqI(SubjectSeqGSAP(gsap)),len,mat+del);
#else
		fprintf(fp,"\n$%d=%d(%d):\n",i,len,mat+del);
#endif
		// fprintf(fp,">"); PutSeqInfo(fp,dsp->sE); 
		// fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
	  	tPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,AB);
	     } fprintf(fp,"\n");
	   }
	  }
	} fprintf(fp,"\n_0].\n");
}

void    tSeqAlignToCMA(FILE *fp,sap_typ head,Int4 leftflank,Int4 rightflank,a_type AB)
{ tSeqAlignToCMA(fp,head,0,INT4_MAX,leftflank,rightflank,AB); }

void    tSeqAlignToCMA(FILE *fp, sap_typ head, Int4 min_sq_id, Int4 max_sq_id, 
	Int4 leftflank, Int4 rightflank, a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// Add leftflank & rightflank.
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
	if(head==0) return;
	UInt4	last_id=0,id;
	Int4		i,n,ins,del,mat,min_qs,max_qe;

	n = NumSeqsListGSAP(head, min_sq_id,max_sq_id, &min_qs, &max_qe);
#if 1
	fprintf(fp,"[0_(1)=psiblast(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",n);
#else
	fprintf(fp,"[0_(1)=psiblast(%d){go=19,gx=2,pn=5.0,lf=%d,rf=%d}:\n",
			n,leftflank,rightflank);
#endif
	fprintf(fp,"(%d)",max_qe-min_qs+1);
	for(i=0; i < (max_qe-min_qs+1); i++) fprintf(fp,"*");
	fprintf(fp,"\n"); i=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
          id = SeqI(SubjectSeqGSAP(gsap));
          if(id < min_sq_id || id > max_sq_id) continue; 
	  dsp_typ dsp=gsap->segs; 
	  if(dsp->subject_id != last_id){
		last_id= dsp->subject_id; i++;
		mat = InDelsSeqAlign(gsap,&ins,&del,min_qs,max_qe);
		Int4 len = LenSeq(dsp->sE);
#if 1	// Use actual sequence identifier.
		fprintf(fp,"\n$%d=%d(%d):\n",last_id,len,mat+del);
#else
		fprintf(fp,"\n$%d=%d(%d):\n",i,len,mat+del);
#endif
		// fprintf(fp,">"); PutSeqInfo(fp,dsp->sE); 
		// fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
	  	tPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,AB);
	  } fprintf(fp,"\n");
	} fprintf(fp,"\n_0].\n");
}

//*************************************************************************
#if 1	// NEW function for conversion from psiblast to cma
//*************************************************************************

void    xPutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank, 
	Int4 min_qs, Int4 max_qe, Int4 SeqID, a_type AB)
// print output to a file in cma format...
// query is used as a template to define length and dimensions of cma.
// WARNING: OFFSETS ARE NOT SET CORRECTLY!!!
{
        Int4	r,s,q,qs,ss,length,index,total;
	UInt4	qid,sid;
	dsp_typ	dsp;
	e_type	qE,sE;
	Int4	qos,sos,m;
	Int4	i,n,ins,del,mat;
	char	*qstr,*sstr,c;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
	qid = dsp->query_id; qE = dsp->qE;
	sid = dsp->subject_id; sE = dsp->sE;

	mat = InDelsSeqAlign(sap,&ins,&del,min_qs,max_qe);
	Int4 len = LenSeq(dsp->sE);
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);

	// delimit subject sequence and create subseq.
	Int4 SS=StartSubjectGSAP(sap);
	Int4 ES=EndSubjectGSAP(sap);
	// SS=1; ES=LenSeq(sE);
	SS = MAXIMUM(Int4,SS-leftflank,1);
	ES = MINIMUM(Int4,ES+rightflank,LenSeq(sE));

        for(total=index=0; index<dsp->numseg; index++){
		total += dsp->lens[index];
	}
	ss = dsp->starts[1];
// printf("%d: SS=%d; ES=%d; LenSeq = %d; $%d=%d(%d)\n",SeqID,SS,ES,LenSeq(sE),
//		SeqID,ES-SS+1,mat+del);
#if 0	// Use actual sequence identifier.
	fprintf(fp,"\n$%d=%d(%d):\n",sid,ES-SS+1,mat+del);
#else
	Int4 eess,ssss;
	eess = ES; ssss=SS;
#if 0	// Fix Bug in this routine...?? or is problem with FakeSearch???
	if(dsp->starts[0] > min_qs){
	       	ssss = ssss + (dsp->starts[0] - min_qs); 
	}
#endif
	fprintf(fp,"\n$%d=%d(%d):\n",SeqID,eess-ssss+1,mat+del);
#endif
	// fprintf(fp,">"); PutSeqInfo(fp,dsp->sE); 
	fprintf(fp,">"); PutSubSeqInfo(fp,SS,ES,dsp->sE); 
	for(r=SS; r <= ss; r++){
		fprintf(fp,"%c",AlphaChar(ResSeq(r,sE),AB));
	} fprintf(fp,"{()");
	{
		qs = dsp->starts[0];
		for(q=min_qs; q < qs; q++) fprintf(fp,"-");
	}
	NEW(qstr,total+3,char); NEW(sstr,total+3,char);
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1 && ss != -1){	// match region...
		    for(i=0,r=ss+1; i < length; r++,i++){
			c= AlphaChar(ResSeq(r,sE),AB);
			fprintf(fp,"%c",c);
		    } qs++;
		} else if(qs == -1){		// insertion in subject sequence.
		    if(index==dsp->numseg-1) fprintf(fp,"(");
		    for(i=0,r=ss+1; i < length; r++,i++){
			c= AlphaChar(ResSeq(r,sE),AB);
			c = tolower(c); fprintf(fp,"%c",c);
		    } if(index==0) fprintf(fp,")");
		} else if(ss == -1){		// deletion in subject sequence.
			for(i=0; i < length; i++) fprintf(fp,"-"); 
		} else print_error("PutSeqAlignToCMA( ): dsp error");
		// if(index==dsp->numseg-1 && qs != -1) fprintf(fp,"(");
        }
	{
		Int4 qe = dsp->starts[2*(dsp->numseg-1)]+length;
		for(q=qe; q<=max_qe; q++) fprintf(fp,"-");
	} fprintf(fp,"()}");
	// for( ; r <= LenSeq(sE); r++){ fprintf(fp,"%c",AlphaChar(ResSeq(r,sE),AB)); }
	for( ; r <= ES; r++){ fprintf(fp,"%c",AlphaChar(ResSeq(r,sE),AB)); }
	fprintf(fp,"*"); fprintf(fp,"\n");
	free(qstr); free(sstr); 
}

void    xSeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank, Int4 rightflank,
	a_type AB)
{ xSeqAlignToCMA(fp,head,INT4_MAX,leftflank,rightflank,AB); }

void    xSeqAlignToCMA(FILE *fp,sap_typ head,Int4 max_sq_id,Int4 leftflank,
	Int4 rightflank, a_type AB)
{ xSeqAlignToCMA(fp,head,0,max_sq_id,leftflank,rightflank,AB); }

void    xSeqAlignToCMA(FILE *fp,sap_typ head,Int4 min_sq_id, Int4 max_sq_id,
	Int4 leftflank, Int4 rightflank, a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// Add leftflank & rightflank.
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
	if(head==0) return;
	Int4		sq,i,n,ins,del,mat,min_qs,max_qe;
	UInt4	id;
#if 0
	fprintf(fp,"[0_(1)=gold_pcna.full(%d){go=19,gx=2,pn=5.0,lf=%d,rf=%d}:\n",
			NumSeqsListGSAP(head),leftflank,rightflank);
#endif
	// n = NumHSPsListGSAP(head, &min_qs, &max_qe);
	n = NumHSPsListGSAP(head,min_sq_id,max_sq_id, &min_qs, &max_qe);
#if 0
	fprintf(fp,"[0_(1)=psiblast(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",n);
#else
	fprintf(fp,"[0_(1)=psiblast(%d){go=19,gx=2,pn=5.0,lf=%d,rf=%d}:\n",
			n,leftflank,rightflank);
#endif
	fprintf(fp,"(%d)",max_qe-min_qs+1);
	for(i=0; i < (max_qe-min_qs+1); i++) fprintf(fp,"*");
	fprintf(fp,"\n"); i=0;
#if 1	// output using sq_id to order....
	sap_typ *SAP;
	NEW(SAP,n+3,sap_typ); 
	dh_type dH = dheap(n+3,4);
	sap_typ gsap;
	for(sq=1,gsap=head; gsap!=NULL; gsap=gsap->next,sq++){
          id = SeqI(SubjectSeqGSAP(gsap));
          if(id < min_sq_id || id > max_sq_id) continue; 
	  SAP[sq]=gsap; insrtHeap(sq,(keytyp)id,dH); 
	}
	while(!emptyHeap(dH)){
	  sq=delminHeap(dH); gsap=SAP[sq];
          id = SeqI(SubjectSeqGSAP(gsap));
	  xPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,id,AB);
	} fprintf(fp,"\n");
	Nildheap(dH);
#else
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
          id = SeqI(SubjectSeqGSAP(gsap));
          if(id < min_sq_id || id > max_sq_id) continue; 
	  dsp_typ dsp=gsap->segs; i++;
	  // fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
#if 1	// Use actual sequence id...
	  xPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,id,AB);
#else
	  xPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,i,AB);
#endif
	} fprintf(fp,"\n");
#endif
	fprintf(fp,"_0].\n");
}

#endif

cma_typ	PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[], 
	const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
	Int4 minfix,BooLean UseAllHSPs,Int4 MaxSqID,cma_typ mtfcma,
	BooLean use_checkin)
{
	Int4	left_flank=0,right_flank=0;
	cma_typ	cma;

	if(qE==0){ qE = SeqSetE(1,data); UseAllHSPs=FALSE; }
	*QS=0;
	sap_typ	sap = RtnSAP_PsiBLAST(qE,data,argc,argv,USAGE,maxrounds,A,minfix,
		IncludeQuery,mtfcma,use_checkin);	// IncludeQuery==TRUE then sampling allowed.
				// this is TRUE currently only for main set.
	if(sap){
	  *QS=QueryStartGSeqAlnList(sap);
#if 1
	  FILE *fp=tmpfile();
#else
	  FILE *fp=open_file("junk_see",".cma","w");
#endif
	  if(IncludeQuery){
	    Int4 *starts,*lens;
	    NEW(starts, 5,Int4); NEW(lens,5,Int4);
	    starts[0]=starts[1]=0; lens[0]=LenSeq(qE);
	    sap_typ head=MakeGSeqAlign(1,SeqI(qE),SeqI(qE),qE,qE,starts,lens); 
	    head->next = sap;
	    if(UseAllHSPs) xSeqAlignToCMA(fp, head, MaxSqID,left_flank, right_flank, A);
	    else QueryFirstSeqAlignToCMA(fp,head,MaxSqID,left_flank,right_flank,A);
	    head->next = NULL; GSeqAlignFree(head);
	  } else {
	    if(UseAllHSPs) xSeqAlignToCMA(fp,sap,MaxSqID,left_flank, right_flank, A);
	    else QueryFirstSeqAlignToCMA(fp,sap,MaxSqID,left_flank,right_flank,A);
	  } 
#if 1
	  rewind(fp); cma=ReadCMSA(fp,A); fclose(fp); 
#else
	  fclose(fp); fp=open_file("junk_see",".cma","r");
	  cma=ReadCMSA(fp,A); fclose(fp); 
#endif
	  FreeGSeqAlignList(sap);
	  return cma;
	} else return 0; 
}

const char USAGE_PsiBLAST[]="USAGE: PsiBLAST [options]\n\
   options:\n\
        -C<file> write checkpoint file\n\
        -c<file> write sequences detected to file\n\
        -D       Don't SEG query sequence\n\
        -e<real> E-value cutoff for showing scores and alignments\n\
        -h<real> E-value cutoff for inclusion in next profile\n\
        -H<int>  heapsize (default 20000)\n\
        -j<int>  maximum number iterations (default: infinite)\n\
                  Input value must be > 1.\n\
        -L       do low complexity masking of database sequences\n\
        -nochk   don't create the default checkpoint file\n\
        -iter_ticks print only the iterations to stderr\n\
        -M<int>  minfix score (default 0)\n\
        -m<int>  print option m=0..6 (default 4)\n\
        -S<int>  Gibbs sample alignments <int> times at each iteration (default 0)\n\
        -T<int>  blast word hit threshold (default 11)\n\
        -X<int>  X dropoff for gapxdrop\n\
\n";

cma_typ	PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[], 
	const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
	Int4 minfix,BooLean UseAllHSPs)
{ return PsiBLAST(QS,qE,data,argc,argv,USAGE,maxrounds,A,
	IncludeQuery,minfix,UseAllHSPs,INT4_MAX,0,FALSE);}

cma_typ	PsiBLAST(Int4 *QS, e_type qE, ss_type data, int argc, char *argv[], 
	const char *USAGE, Int4 maxrounds, a_type A, BooLean IncludeQuery,
	Int4 minfix,BooLean UseAllHSPs,cma_typ mtfcma)
// This is used for mcma only...
{ return PsiBLAST(QS,qE,data,argc,argv,USAGE,maxrounds,A,
	IncludeQuery,minfix,UseAllHSPs,INT4_MAX,mtfcma,TRUE);}

sap_typ	RtnSAP_PsiBLAST(e_type qE, ss_type data, int argc, char *argv[], 
	const char *USAGE, Int4 maxrounds, a_type A, Int4 minfix, 
	BooLean use_gibbs,cma_typ mtfcma,BooLean use_checkin)
{
	Int4    T=11;
	time_t	time1=time(NULL);
	Int4	left_flank=0,right_flank=0;
        double	x_parameter=25.0;
	// Int4	maxrounds=1;
	Int4	width=60;
	UInt4	hpsz=20000,printopt=0;
	double	Ethresh=0.001,Ecutoff=0.001;
	BooLean	verbose=FALSE;
	BooLean	seg_seq=FALSE,success=FALSE;
	char	mask_dbs=' ';
	FILE    *seqfile=NULL;
	BooLean	GetMatrix=TRUE,fake_cma_srch=FALSE,notchk=FALSE,iter_ticks=FALSE;
	char	*checkout=0,str[200];
	char	*checkin=0; // recover...
	double	misalncut=2.0;
	UInt4	sampling_size=0;

	if(argc < 2) print_error(USAGE);
        for(Int4 arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') continue; // filenames... print_error(USAGE);
           switch(argv[arg][1]) {
	     case 'C': checkout = argv[arg]+2; break;
	     case 'c': seqfile = open_file(argv[arg]+2,"","w"); break;
             case 'D': seg_seq=TRUE; break;
             case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,USAGE); break;
	     case 'F': fake_cma_srch=TRUE; break;
	     case 'H': hpsz=IntOption(argv[arg],'H',1,1000000,USAGE); break;
             case 'h': Ethresh=RealOption(argv[arg],'h',0.0,10000,USAGE); break;
	     case 'j': maxrounds=IntOption(argv[arg],'j',0,1000,USAGE); break;
             case 'i': 
	     	if(strcmp("-iter_ticks",argv[arg]) == 0){ iter_ticks=TRUE; }
		break;
             case 'L': mask_dbs='x'; break;
	     case 'm': printopt=IntOption(argv[arg],'m',0,6,USAGE); break;
	     case 'n': if(strcmp("-nochk",argv[arg]) == 0){ notchk=TRUE; }
		break;
	     case 'R': if(use_checkin) checkin = argv[arg]+2; break;
	     case 'S': sampling_size=IntOption(argv[arg],'S',0,100,USAGE); break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,USAGE); break;
	     case 't': misalncut=RealOption(argv[arg],'t',0.0,1.0,USAGE); break;
             case 'v': verbose=TRUE; break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,USAGE); break;
	     case ' ': break; // ignore
	     case 'r': break; // ignore
	     case 'u': break; // ignore
             default: print_error(USAGE); break;
           }
        }
	if(checkout==0 && notchk==FALSE){
		sprintf(str,"%s.chk",argv[1]);
		checkout = str; // make checkpoint file by default
	}
	if(!use_gibbs) sampling_size=0;
	cma_typ	cma=0;

	// if(qE==0){ qE = SeqSetE(1,data); UseAllHSPs=FALSE; } // taken care of from calling
	// UseAllHSPs=TRUE;
#if 1	//***************************** Fake Search *************************
	gpsi_type *fake_gpsi=0;
	ss_type data1=0;
        if(mtfcma){
            data1=TrueDataCMSA(mtfcma);
	    PutSeq(stderr,qE,A); PutSeq(stderr,SeqSetE(1,data1),A);
	    assert(IdentSeqs(qE,SeqSetE(1,data1)));
            // use correct qE that is first in same as before...            
        }
#endif	//***************************** Fake Search *************************
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
	gpsi.KeepQuite();	// always use this option with these programs...
	if(seg_seq) ProcessSeqPSeg(17,2.2,2.5,100,qE,A);
	Int4 iteration=1;
	Int4	**mx=0;
        do { 
	   sap_typ sap;
#if 0	// NOTE: Not working...yet.
	   if(iteration==1 && mtfcma){	// new...
	     sap=gpsi.SearchWithPrior(T,misalncut,mask_dbs,mtfcma);
	   } else 
#endif
	   {
	     if(misalncut < 1.0) sap = gpsi.TrimmedSearch(T,misalncut,mask_dbs);
	     else sap = gpsi.Search(T,mask_dbs);
	   }
	   // if(sap) DriverNumberOfGapsGSAP(stderr,sap);
	   if(sap){ 
		if(iter_ticks){
			fprintf(stderr,"iteration %d: %d sequences aligned\n",
						iteration,NumSeqsListGSAP(sap));
		}
#if 1	// Append the fake_saps to this gpsi...
	     //***************************** Fake Search *************************
		   if(mtfcma && gpsi.NotConverged()){	// Then append the results and 
			FILE *null_fp=0;
	    		fake_gpsi = new gpsi_type(qE,data1,Ethresh,Ecutoff,
			    		x_parameter,hpsz,maxrounds,null_fp);
			fake_gpsi->KeepQuite();	// always use this option with these programs...
	    		fake_gpsi->FakeSearchCMSA(T,mtfcma);
			gpsi.AppendFakeResults(fake_gpsi);
	    		// Don't free this now, as I am giving the sap and results 
			// to the real gpsi( );
		   }
	     //***************************** Fake Search *************************
#endif
		   success=TRUE; 
	   } else break;
	   if(verbose) gpsi.Test('M');
	   if(minfix > 0){	// not sure what this is (finding informative regions)?
		// RemoveOverlapsGSAP(sap,mx,A);
	   	// if(minfix > 0) FixAlnRunOverGSAP(sap,mx,A,minfix);
		gpsi.ComputeMatrix(checkout);
		// if(gpsi.NotConverged( ) && mx) {  // for MargProb...
		if(mx) { 
			ComputeQueryScoreDSP(stderr,qE,mx); // needs double **mx;
			for(Int4 s=0; s < LenSeq(qE); s++) if(mx[s]) free(mx[s]); 
			free(mx);
		} if(gpsi.NotConverged( )) mx=gpsi.CopyOFposMatrix(); 
		// NOTE: mx starts from zero!!
	   } else gpsi.ComputeMatrix(checkout);
#if 1   // Add simple sampling procedure? // WARNING: Not at all tested!!
           for(Int4 s=1; s <= sampling_size; s++){
             if(gpsi.posMatrix){ // test sampling routines...
                sap_typ sap = gpsi.SampleMatrix(checkout); // recomputes matrix...
                if(verbose) gpsi.Test('M');
                // if(verbose) gpsi.ShowAlign(stdout,width,printopt);
             }
           }
#endif
	   iteration++;
	} while(gpsi.NotConverged( ));
	gpsi.ShowAlign(stdout,width,printopt);
	if(seqfile){ gpsi.PutSeqs(seqfile); fclose(seqfile); }
	if(verbose) fprintf(stderr,"time = %0.1f seconds\n",difftime(time(NULL),time1));
	if(success){
	  sap_typ sap = gpsi.StealSAP();
	  return sap;
	} else return 0; 
}


