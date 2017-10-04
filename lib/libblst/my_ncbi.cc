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

#include "my_ncbi.h"

void     *GMemNew(size_t size)
{ void *x_ptr=0; x_ptr=malloc(size); memset(x_ptr,0,size); return x_ptr; }

void * GMemFree(void *ptr){ if (ptr != NULL) free(ptr); return NULL; }

sap_typ GSeqAlignNew(void)
{ return (sap_typ)GMemNew(sizeof(sap_type)); }

sap_typ MakeSelfAlignedGSAP(e_type E)
// use this to create a sap of the sequence E aligned with itself
{
        Int4 *starts,*lens;
        NEW(starts, 3,Int4); NEW(lens,3,Int4);
        starts[0]=starts[1]=0; lens[0]=LenSeq(E);
        return MakeGSeqAlign(1,SeqI(E),SeqI(E),E,E,starts,lens);
}

sap_typ MakeGSeqAlign(Int2 numseg, UInt4 query_id, 
	UInt4 subject_id, e_type qE, e_type sE, Int4 *starts,
	Int4 *lens)
// not used in psiblast code yet...
{
        sap_typ sap = GSeqAlignNew();
        sap->dim =2; // only two dimention alignment.
        dsp_typ dsp = GDenseSegNew(); // Denseg Object is only type for GSeqAlign.
        dsp->dim = 2; dsp->numseg = numseg;
        dsp->subject_id = SeqI(sE); dsp->sE=sE;
        dsp->query_id=0; dsp->qE=qE;
        dsp->starts = starts; dsp->lens = lens;
        sap->segs = dsp; sap->next = NULL;
        return sap;
}

BooLean	IsSplitGSeqAlign(register sap_typ sap)
// If alignment is split return true else return false.
{
	if(sap == NULL) return FALSE;
	register char	state=' ';
	register Int4	i;
	register Int4	n=sap->segs->numseg;
	register Int4	*starts=sap->segs->starts;
        for(i=0; i < n; i++){
	     if((starts[2*i]) == -1){
		if(state=='I'){ return TRUE; } state='I';
	     } else if((starts[2*i+1]) == -1){
		if(state=='D'){ return TRUE; } state='D';
	     } else {
		if(state=='M'){ return TRUE; } state='M';
	     }
	} return FALSE;
}

BooLean	FixSplitsGSeqAlign(sap_typ sap)
// merge segs that can be merged...
// This fixes a gapxdrop bug.  This is an AFN function...
{
	if(sap == NULL) return FALSE;
	register dsp_typ dsp= sap->segs;
	register char	state=' ';
	register BooLean	merge;
	register Int4	i,j;

        for(i=0; i < dsp->numseg; i++){
	     merge = FALSE;
	     if((dsp->starts[2*i]) == -1){
		if(state=='I'){ merge = TRUE; } state='I';
	     } else if((dsp->starts[2*i+1]) == -1){
		if(state=='D'){ merge = TRUE; } state='D';
	     } else {
		if(state=='M'){ merge = TRUE; } state='M';
	     }
	     if(merge){		// delete lens[i] and move all segs up by one.
		   dsp->lens[i-1] += dsp->lens[i];
		   dsp->numseg--;
        	   for(j=i; j < dsp->numseg; j++){
			dsp->starts[2*j] = dsp->starts[2*(j+1)];
			dsp->starts[(2*j)+1] = dsp->starts[(2*(j+1))+1];
		   	dsp->lens[j] = dsp->lens[j+1];
		   } i--; 
	     }
	} return TRUE;
}

#define DEBUG_PutMultiGSeqAlign 1
void    PutMultiGSeqAlign(FILE *fp, sap_typ head, Int4 width, a_type AB)
#if 0
	WARNING: This code was core dumping for ce2cma when subject sequences
	having deletions at beginning or end of alignment are input!!!
	Insure indicates array bound erros at:
		if(state != 'I') for(j=0; j < maxins[q]; j++) str[row][t++] = '.';
	for maxins[q] and str[row][t++]...!
	I have not fixed this yet....
#endif
{
	UInt4	qid,sid;
	dsp_typ	dsp;
	sap_typ	gsap;
	e_type		qE,sE;
        Int4		r,i,j,k,s,q,qs,ss,t,len,end,m,N,template_len,querySize,row;
        Int4    	*maxins,*cumlen,*start,*sqid,*Cterm,*Nterm,cum_ins_len;
	char		**str,c,state;

	if(head == NULL) return;

	// 0. get key information ...
	dsp = head->segs; qid = dsp->query_id; qE = dsp->qE;
	querySize=LenSeq(qE);
        NEW(maxins,querySize+3, Int4); // initialized to zero inserts

	// 1. compute template length...
        for(N=0,gsap=head; gsap!=NULL; gsap=gsap->next){
	   // PutGSeqAlign(fp, gsap, width, AB); // DEBUG
	   N++; dsp = gsap->segs;
           q = dsp->starts[0]; s = dsp->starts[1];
	   assert(q > -1); q=q-1;	// store max insert right after q. 
	   // Tolerate no deletions in query at begin!!!
	   // assert(q > 0); q=q-1;	// q is less than 0 but must get incremented...
#if 0	//**********************************************************************
	// gapxdrop sometimes returns two qs == -1's in a row, in which case 
	// we need to increment the previous length rather than take the current
	// length to avoid errors.
#endif	//**********************************************************************
	   cum_ins_len=0;
           for(i=0; i < dsp->numseg; i++){
                qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
                len = dsp->lens[i];
#if 1
                if(qs == -1){	 // == an insert relative to query after 'q'
			if(len > maxins[q]) maxins[q] = len;
		} else { 
			q+=len; 
		}
#else
                if(qs == -1){	 // == an insert relative to query after 'q'
		   cum_ins_len += len;
		   if(cum_ins_len > len){ 
			if(maxins[q] < cum_ins_len) maxins[q] = cum_ins_len;
		   } else if(len > maxins[q]) maxins[q] = len;
		} else { q+=len; cum_ins_len=0; }
#endif
           }
        }
        NEW(cumlen,querySize+2, Int4); // cumulative template length < q.
	for(template_len=0, q=0; q< querySize; ){
		template_len += 1 + maxins[q]; q++;
		cumlen[q]=template_len;
	}
	NEWP(str,N+2,char); for(i=0;i<=N;i++) NEW(str[i],template_len+2,char);
        NEW(start,N+2, Int4); // starting positions.
        NEW(sqid,N+2, Int4); // sequence id.
        NEW(Nterm,N+2, Int4); NEW(Cterm,N+2, Int4);
	// fprintf(stderr,"template_len = %d\n",template_len);

	// 2. store multiple alignment...
        for(row=1,gsap=head; gsap!=NULL; gsap=gsap->next){
	   dsp = gsap->segs; 
	   sid = dsp->subject_id; sE = dsp->sE; sqid[row]=sid;
	   q = dsp->starts[0]; s = dsp->starts[1];
	   start[row]=s+1; t=cumlen[q];	// jump to start in template
	   Nterm[row] = t; state=' ';
	   for(i=0; i<t; i++){ str[row][i]=' '; }
	   for(m=i=0; i < dsp->numseg; i++){
		qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
                len = dsp->lens[i]; end = len-1;
		if(qs != -1 && ss != -1){  // matches at this point.
		   q = qs-1; // first fill in previous insert...
		   if(q >= 0 && state=='D')
			for(j=0; j < maxins[q]; j++) str[row][t++] = '.';
		   for(q++,j=0,r=ss+1; j < end; r++,j++){
			str[row][t++] = AlphaChar(ResSeq(r,sE),AB); 
			for(k=0; k < maxins[q]; k++){ str[row][t++] = '.'; }
			q++;
		   } str[row][t++] = AlphaChar(ResSeq(r,sE),AB); 
#if DEBUG_PutMultiGSeqAlign
	   if(state == 'M'){
		GSeqAlignPut(fp,gsap); // DEBUG
	   	PutGSeqAlign(fp, gsap, width, AB); // DEBUG
		assert(state != 'M');
	   }
#endif
		   state = 'M';
		} else if(qs == -1){ 	   // insert at this point
		   for(j=0,r=ss+1; j < len; r++,j++){ 
				c = AlphaChar(ResSeq(r,sE),AB);
				str[row][t++] = tolower(c);
		   }
		   // note: q retained from last segment!!!
		   for( ; j < maxins[q]; j++){ str[row][t++] = '.'; }
#if DEBUG_PutMultiGSeqAlign
	   if(state == 'I'){
		GSeqAlignPut(fp,gsap); // DEBUG
	   	PutGSeqAlign(fp, gsap, width, AB); // DEBUG
		assert(state != 'I');
	   }
#endif
		   state = 'I';
		} else {		   // deletion at this point
		   q = qs-1; // if(state != 'I') then first fill in previous insert...
		   if(state != 'I') for(j=0; j < maxins[q]; j++) str[row][t++] = '.';
		   for(q++,k=0; k < end; k++){ 
			str[row][t++] = '-'; // one for query match
			for(j=0; j < maxins[q]; j++) str[row][t++] = '.';
			q++;
		   } str[row][t++] = '-'; // one for query match
#if DEBUG_PutMultiGSeqAlign
	   if(state == 'D'){
		GSeqAlignPut(fp,gsap); // DEBUG
	   	PutGSeqAlign(fp, gsap, width, AB); // DEBUG
		assert(state != 'D');
	   }
#endif
		   state = 'D';
		} 
	   }
	   if(t > template_len){
		GSeqAlignPut(fp,gsap); // DEBUG
	   	PutGSeqAlign(fp, gsap, width, AB); // DEBUG
		fprintf(fp,"t = %d; template_len = %d\n",t,template_len);
	   	assert(t <= template_len);
	   }
	   Cterm[row] = t-1;
	   for( ; t < template_len; t++) str[row][t]=' ';
	   row++;
	}

	// 3. output multiple alignment...
	for(i=0; i < template_len; i+=width){
	   end = MINIMUM(Int4,template_len,i+width);
	   for(row=1; row <= N; row++){
	    if(Nterm[row] >= end  || Cterm[row] < i) continue;
	    ss = start[row];
	    fprintf(fp,"%-6d %-4d ",sqid[row],ss);
	    for(j=i; j < end; j++){ 
		fprintf(fp,"%c",str[row][j]); 
		if(isalpha(str[row][j])) ss++;
	    } fprintf(fp," %d\n",ss-1); start[row] = ss;
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n\n");

	// 4. deallocate memory ...
	for(i=0;i<=N;i++) free(str[i]);
        free(maxins); free(str); free(cumlen); free(start); free(sqid);
	free(Nterm); free(Cterm);
}

Int4    LengthListGSAP(sap_typ head){ return NumHSPsListGSAP(head); }

Int4	NumHSPsListGSAP(sap_typ head)
// return the number of gapped aligned segments in list...
{
	UInt4	n=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next) n++;
	return n;
}

sap_typ SortBySeqIDGSAP(sap_typ sap)
// Contains private info!  Need to move this to sap_typ.c
// needs dheap.h file...
{
        sap_typ next,*array;
        Int4    n,i;

        if(sap == 0) return sap;

        n=NumHSPsListGSAP(sap);
        dh_type dH=dheap(n+3,4);
        NEW(array,n+3,sap_typ);

        for(n=1,next=sap; next != NULL; next=NextGSAP(next)){
                e_type  sE=SubjectSeqGSAP(next);
                Int4 id = SeqI(sE);
                Int4 score = ScoreGSAP(next);
                // assert(score >= 1);
                insrtHeap(n,((double)id + (double)1.0/score),dH);
                array[n] = next; n++;
        }
        for(next=0; (i=delminHeap(dH)) != 0; ){
           if(next){ next->next=array[i]; next = next->next; }
           else { sap=next=array[i]; }
        } next->next = NULL;
        Nildheap(dH); free(array);
        return sap;
}

Int4	NumSeqsListGSAP(sap_typ head)
// return the number of sequences in list...
{
	UInt4	n=0,last_id=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
		dsp_typ dsp=gsap->segs; 
		if(dsp->subject_id != last_id){
			last_id=dsp->subject_id; n++;
		}
	} return n;
}

#if 1	// output for vsi files...
void    PutSameSeqsListGSAP(FILE *fp, sap_typ head, Int4 MinLen, double fract_IDs, a_type AB)
// output sequences in list...
{
	UInt4	last_id=0;
	Int4	minlen,alnlen, Ident, Inserts;
	double	fraction;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	   minlen=IdentitiesGSAP(&alnlen,&Ident,&Inserts,gsap);
	   if(Inserts==0){
	     fraction=(double)Ident/(double)alnlen;
	     if(fraction >= fract_IDs && alnlen >= MinLen ){
		dsp_typ dsp=gsap->segs; 
		if(dsp->subject_id != last_id){
			last_id= dsp->subject_id;
			PutSeq(fp,dsp->sE,AB);
		}
	     }
	   }
	}
}

#endif

Int4	PutSeqsListGSAP(FILE *fp, sap_typ head, a_type AB)
// output sequences in list...
{
	UInt4	last_id=0;
	Int4	n=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
		dsp_typ dsp=gsap->segs; 
		if(dsp->subject_id != last_id){
			last_id= dsp->subject_id;
			PutSeq(fp,dsp->sE,AB);
			n++;
		}
	} return n;
}

void    PutSubSeqsListGSAP(FILE *fp, sap_typ head, Int4 left, Int4 right,
	a_type AB)
{
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  dsp_typ dsp=gsap->segs; 
	  PutGSubSeq(fp, gsap, left, right,AB);
	} 
}

Int4	InfoGSAP(Int4 *QS, Int4 *QE, Int4 *SS, Int4 *SE, sap_typ sap)
{
        Int4		s,q,qs,ss,length,index;
	e_type		qE,sE;
	Int4		qos,sos;

        assert(sap->segs != NULL);
	dsp_typ	dsp = sap->segs; 
	qE = dsp->qE; sE = dsp->sE;
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);
        for(s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1) q+=length;
		if(ss != -1) s+=length;
        }
	*QS = dsp->starts[0] +1 + qos; *SS = dsp->starts[1]+1 + sos;
	*QE = *QS+q-1; *SE = *SS+s-1;
	return s;
}

sap_typ ConcatenateGSAP(sap_typ head, sap_typ sap)
// Attach sap to tail of head.
{
        if(head == 0) return sap;
        else if(sap == 0) return head;
        for(sap_typ tail=head; tail != 0 ; tail=tail->next){
                if(tail->next== 0){ tail->next=sap; break; }
        } return head;
}

void    PutMergedGSAP(FILE *fp, sap_typ head, Int4 left, Int4 right,
        a_type AB)
/*********************************************************************
        if(s2 <= e1 || s1 <= e2)   then segments overlap and merge
   s1            e1                             s1       e1
   ================-------------   -------------===========-----
   -------------===========-----   ================-------------
                s2       e2        s2            e2
 *********************************************************************/
{
	ds_type 	sets;
	sap_typ		sap1,sap2;
	dsp_typ    dsp1,dsp2;
	Int4		N,s,s1,s2,i,j;
	Int4		qs,qe,ss1,se1,len;
	Int4		*start,*end;

	N = NumSeqsListGSAP(head);
	sets = DSets(N);
	NEW(start,N+3,Int4); NEW(end,N+3,Int4);
	for(i=1,sap1=head; sap1!=NULL; sap1=sap1->next,i++){
		len=InfoGSAP(&qs,&qe,&ss1,&se1,sap1);
		start[i]=ss1; end[i] = se1;
	}
	for(i=1,sap1=head; sap1!=NULL; sap1=sap1->next,i++){
	  s1 = findDSets(i,sets);
	  dsp1=sap1->segs; 
	  for(j=i+1,sap2=sap1->next; sap2!=NULL; sap2=sap2->next,j++){
		dsp2 = sap2->segs;
		if(dsp1->subject_id == dsp2->subject_id){
		   s2 = findDSets(j,sets);
		   if(s1 != s2){ 
		   	if((start[s2] <= end[s1] && end[s2] >= start[s1])
				|| (start[s1] <= end[s2] && end[s1] >= start[s2])){
				s = linkDSets(s1,s2,sets); 
				start[s] = MINIMUM(Int4,start[s1],start[s2]);
				end[s] = MAXIMUM(Int4,end[s1],end[s2]);
#if 0
	fprintf(stderr,"Fusing %d and %d: start = %d, end = %d\n",
		s1,s2,start[s],end[s]);
#endif
				s1 = s;
			}
		   }
		}
	  }
	} 
	for(s=1,sap1=head; sap1!=NULL; sap1=sap1->next,s++){
	   s1 = findDSets(s,sets);
	   if(s == s1){
	        dsp1=sap1->segs;
		PutSubSeq(fp, start[s]-left, end[s]+right, dsp1->sE,AB);
	   }
	} free(start); free(end); NilDSets(sets);
}

BooLean	IsQueryGSAP(sap_typ sap)
// if sap query sequence is identical to sap subject sequence
// then return TRUE; otherwise return false.
{ dsp_typ dsp=sap->segs; return (IdentSeqs(dsp->sE, dsp->qE)); } 

#if 0
BooLean	SubjectOrQuerySubSeqGSAP(sap_typ sap)
// if sap 
{
}
#endif

BooLean	OverlapGSAP(sap_typ sap1, sap_typ sap2)
// Do sap1 and sap2 beInt4 to the same sequence and do regions overlap?
// WARNING: NOT TESTED!!!
{
	Int4		qs,qe,ss1,ss2,se1,se2;
	Int4		len1,len2;
	dsp_typ    dsp1,dsp2;

	dsp1 = sap1->segs; dsp2 = sap2->segs;
	if(dsp1->subject_id != dsp2->subject_id) return FALSE;
	len1=InfoGSAP(&qs,&qe,&ss1,&se1,sap1);
	len2=InfoGSAP(&qs,&qe,&ss2,&se2,sap2);
	if(ss2 <= se1 || ss1 <= se2) return TRUE;	
	return FALSE;
}

void    PutGSubSeq(FILE *fp, sap_typ sap, Int4 left, Int4 right, a_type AB)
{
        Int4		ss,length,index,total;
	Int4		start,end;
	dsp_typ	dsp;
	e_type		sE;

        assert(sap->segs != NULL);
	dsp = sap->segs; sE = dsp->sE;
	start = dsp->starts[1]+1; // subject start site.
        for(total=index=0; index < dsp->numseg; index++){
		ss = dsp->starts[2*index+1]; // subject starts...
		if(ss != -1){
                     length = dsp->lens[index];
		     total+=length;
		} 
        } end = start + total;
	PutSubSeq(fp, start-left, end+right, dsp->sE,AB);
}

void    PutOneSeqAlign(FILE *fp, sap_typ sap, Int4 width, a_type AB)
// put only the first sequence's alignments
{
	UInt4   last_id=0;
	sap_typ	gsap=0;
        dsp_typ dsp=sap->segs; last_id=dsp->subject_id;
        fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE);
        fprintf(fp,"\t\tlength = %d\n\n",LenSeq(dsp->sE));
	for(gsap=sap; gsap!=NULL; gsap=gsap->next){
          dsp=gsap->segs;
          if(dsp->subject_id == last_id){
		PutGSeqAlign(fp,gsap,width,AB);
          } else break;
        } fprintf(fp,"\n");
}

sap_typ	MinEvalSAP(sap_typ head, double *MinEval)
// set MinEval to the minimum evalue for the first subject sequence in head.
// return the sap for the next sequence or null if end of list
{
	sap_typ sap=0;
	UInt4	id=0;
	double	eval,min=DBL_MAX;

	if(!head){ *MinEval=DBL_MAX; return 0; }
	dsp_typ dsp=head->segs; id = dsp->subject_id;
	for(sap=head; sap; sap=sap->next){
		dsp=sap->segs;
		if(dsp->subject_id == id){
			eval = EvalueGSAP(sap);
			if(eval < min) min = eval;
		} else { *MinEval = min; return sap; } // found another sap;
	}
	assert(sap == 0);
	*MinEval = min;	// reached end of list; sap == 0;
	return sap;
}

void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
{
	UInt4	last_id=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  dsp_typ dsp=gsap->segs; 
	  if(dsp->subject_id != last_id){
		last_id= dsp->subject_id;
		fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE); 
		fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
		// fprintf(fp,"\tlength = %d; ",LenSeq(dsp->sE));
		// fprintf(fp,"\tQuery: "); PutSeqInfo(fp,dsp->qE); 
		// fprintf(fp,"\n\n");
	  } PutGSeqAlign(fp,gsap,width,AB);
	} fprintf(fp,"\n");
}

void    PutGSeqAlign(FILE *fp, sap_typ sap, Int4 width, a_type AB)
{
        Int4		r,i,j,s,q,qs,ss,length,index,total;
	UInt4	qid,sid;
	dsp_typ	dsp;
	e_type		qE,sE;
	Int4		qos,sos;
	Int4		end,m;
	char		*qstr,*sstr;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
	qid = dsp->query_id; qE = dsp->qE;
	sid = dsp->subject_id; sE = dsp->sE;
        assert(qE && sE);
// PutSeq(stderr,sE,AB);
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);
        for(total=index=0; index<dsp->numseg; index++){
		total += dsp->lens[index];
	}
	NEW(qstr,total+3,char); NEW(sstr,total+3,char);
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1){
		   if(ss == -1){
			for(i=0,r=qs+1; i < length; r++,i++){
			  qstr[q++] = tolower(AlphaChar(ResSeq(r,qE),AB)); 
			}
		   } else {
			for(i=0,r=qs+1; i < length; r++,i++){
			  qstr[q++] = AlphaChar(ResSeq(r,qE),AB); 
			}
		   }
		} else {
			for(i=0; i < length; i++){ qstr[q++] = '-'; }
		}
		if(ss != -1){
		   if(qs == -1){
			for(i=0,r=ss+1; i < length; r++,i++){
			  sstr[s++] = tolower(AlphaChar(ResSeq(r,sE),AB));
			}
		   } else {
			for(i=0,r=ss+1; i < length; r++,i++){
			  sstr[s++] = AlphaChar(ResSeq(r,sE),AB);
			}
		   }
		} else {
			for(i=0; i < length; i++){ sstr[s++] = '-'; }
		} 
        } qs = dsp->starts[0] +1; ss = dsp->starts[1]+1; 
	Int4 aln_len;
	// fprintf(fp,"\tScore = %.1f bits (%d), Expect = %1.2g\n",
	double percent_ids=PercentIdentGSAP(&aln_len, sap);
	fprintf(fp,"\tScore = %.0f bits (%d), Expect = %1.2g\n",
		sap->bit_score,sap->score,sap->evalue);
	fprintf(fp,"\tIdentities = %d/%d(%.1f)\n",
		(Int4)(percent_ids*aln_len),aln_len,100.0*percent_ids);
	for(i=0; i < total; i+=width){
	    end = MINIMUM(Int4,total,i+width);
	    fprintf(fp,"\n QUERY: %-4d ",qs+qos);
	    for(j=i; j < end; j++){ 
		fprintf(fp,"%c",qstr[j]);  
		if(qstr[j] != '-') qs++;
	    }
	    fprintf(fp," %d\n             ",qs-1+qos);
	    for(j=i; j < end; j++){ 
		if(qstr[j] == '-' || sstr[j] == '-') fprintf(fp," ");
		else {
		  q=AlphaCode(qstr[j],AB); s=AlphaCode(sstr[j],AB);
		  if(q==s) fprintf(fp,"%c",sstr[j]); 
		  else if(valAlphaR(q,s,AB) > 0) fprintf(fp,"+"); 
		  else fprintf(fp," ");
		}
	    }
	    fprintf(fp,"\n SBJCT: %-4d ",ss+sos);
	    for(j=i; j < end; j++){ 
		fprintf(fp,"%c",sstr[j]); 
		if(sstr[j] != '-') ss++;
	    }
	    fprintf(fp," %d\n",ss-1+sos);
	}
	fprintf(fp,"\n");
	free(qstr); free(sstr); 
}

void	GSeqAlignPut(FILE *fp, sap_typ gsap)
{
	fprintf(fp,"dim = %d\n",gsap->dim);
	GDenseSegPut(fp, gsap->segs);
	fprintf(fp,"\n");
}

sap_typ GSeqAlignFree(sap_typ sap)
{
    if (sap == NULL) return (sap_typ)NULL;
    GDenseSegFree((dsp_typ)sap->segs);
    return (sap_typ)GMemFree(sap);
}

Int4	IdentitiesGSAP(Int4 *alnlen, Int4 *Ident,Int4 *Inserts, sap_typ sap)
{
        Int4		i,j,qs,ss,length,r_q,r_s;
	double		total,ident,insert;
	e_type		qE,sE;

        assert(sap != NULL && sap->segs != NULL);
	dsp_typ dsp = sap->segs; qE = dsp->qE; sE = dsp->sE;
        for(total=ident=insert=0.0,i=0; i <dsp->numseg; i++){
		qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
		if(qs != -1 && ss != -1){
                   length = dsp->lens[i]; r_q=qs+1; r_s=ss+1;
		   for(j=0; j < length; j++){
			if(ResSeq(r_q,qE)==ResSeq(r_s,sE)) ident++;
			total++; r_q++; r_s++;
		   }
		} else insert++;
        }
	*alnlen = total; *Ident = ident; *Inserts=insert;
	return MINIMUM(Int4,LenSeq(qE),LenSeq(sE));
}

double	PercentIdentGSAP(Int4 *net_len, sap_typ sap)
{
        Int4		i,j,qs,ss,length,r_q,r_s;
	double		total,ident;
	e_type		qE,sE;

        assert(sap != NULL && sap->segs != NULL);
	dsp_typ dsp = sap->segs; qE = dsp->qE; sE = dsp->sE;
        for(total=ident=0.0,i=0; i <dsp->numseg; i++){
		qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
		if(qs != -1 && ss != -1){
                   length = dsp->lens[i]; r_q=qs+1; r_s=ss+1;
		   for(j=0; j < length; j++){
			if(ResSeq(r_q,qE)==ResSeq(r_s,sE)) ident++;
			total++; r_q++; r_s++;
		   }
		}
        } *net_len = total;
	return ident/total;
}

//===================== From my old psi-blast code =============

void	FreeGSeqAlignList(sap_typ sap)
{
	while(sap != NULL){
		sap_typ tmp = sap->next; sap->next=NULL;
		GSeqAlignFree(sap); sap = tmp;
	}
}

// get the number of segments in the alignment.
// operation: 
// 'i' = insertion in sequence outside of motifs
// 'I' = insert in sequence within motif (need to delete this)
//  'M' = match to start of a motif block
//  'm' = match to other sites in motif
//  'D' = deletion of sequence within motif
//  'd' = deletion of sequence outside of motif (not used)

static Int4	NumSegsSeqAlign(char *operation)
{
	register Int4	n,i;
	register char	state=' ',o;

	for(n=0, i=1; operation[i] != 'E'; i++){
	    o = operation[i];
	    switch(o){
		case 'm': case 'M': if(state!='m'){ state='m'; n++; } break;
		case 'i': case 'I': if(state!='i'){ state='i'; if(n>0) n++; } break;
		case 'd': case 'D': if(state!='d'){ state='d'; if(n>0) n++; } break;
		default: print_error("This should not happen");
	    }
	}
	return n;
}

sap_typ	ToGSeqAlign(Int4 numopers, char *operation, e_type qE, e_type sE, 
	Int4 start1, Int4 start2)
// taken from my 'GetNCBISegsSeqAlign( )' routine in BLAST/seqalign.c.
// converts an alignment operation array from my code to my generic NCBI
// seq_align_pointer...
{
        sap_typ	sap;
	Int4	n,i,k,len,s1,s2;
	char	o,state;

#if 0 // prune end of operations...
	Int4 i=numopers;
	char o;

	do {
	    i--;
	    if(i < 0){
		fprintf(stderr,"operations = %s (%d)\n",operation,numopers);
		 print_error("ToNCBISeqAlign( ) input error");
	    }
	    o = operation[i]; 
	} while(o != 'm' && o != 'M');
	i++; operation[i] = 'E';
#endif
        sap = GSeqAlignNew();
	sap->next=NULL;
        sap->dim = 2;		// two sequences
        dsp_typ segs = GDenseSegNew();
        segs->dim = 2;
	n = NumSegsSeqAlign(operation);
        segs->numseg = n;
        segs->query_id = SeqI(qE); segs->qE = qE;
        segs->subject_id = SeqI(sE); segs->sE = sE;

        NEW(segs->lens,n+3,Int4);
        NEW(segs->starts,2*n+3,Int4);
	state = ' ';
	for(len=0,k=n=0, s1=start1,s2=start2, i=1; operation[i] != 'E'; i++){
	    o = operation[i];
	    switch(o){
		case 'm': case 'M': 
		  if(state!='m'){
		  	if(n > 0){ segs->lens[k] = len; k++; }
			segs->starts[n] = s1; n++; segs->starts[n] = s2; n++;
			state='m'; len=0;
		  } s1++; s2++; len++;
	          break;
		case 'i': case 'I': 
		  if(state!='i'){ 
		     if(n > 0){
			segs->lens[k] = len; k++;
			segs->starts[n]=-1; n++; segs->starts[n]=s2; n++;
		     } state='i'; len=0;
		  } s2++; len++;
		  break;
		case 'd': case 'D':
		  if(state!='d'){
		     if(n > 0){
			segs->lens[k] = len; k++;
			segs->starts[n]=s1; n++; segs->starts[n]=-1; n++;
		     } state='d'; len=0; 
		  } s1++; len++;
		  break;
		default: print_error("This should not happen");
	    }
	} segs->lens[k] = len; // get length of last segment
        sap->segs = segs;
        return sap;
}

void	PutRasMolGSAP(FILE *fp, sap_typ sap, Int4 N_Colors, const char **Colors)
// print out a rasmol script of sap
{
   assert(sap && sap->segs);
   Int4 num_segs = sap->segs->numseg;
   Int4 *seq_starts = sap->segs->starts;
   Int4 *segment_lens = sap->segs->lens;
   Int4  j,i=0, c=0,site;

  // fprintf(fp,"set ambient 60\n");
  // fprintf(fp,"background gray\n");
  fprintf(fp,"backbone 1\ncolor [60,60,60]\n");
  for(j=c=i=0; j < num_segs; j++, i=i+2) {
    if(seq_starts[i] < 0 || seq_starts[i+1] < 0) continue;
    site = seq_starts[i] + 1;
    fprintf(fp,"select %d-%dA\n color %s\n",
        site, site + segment_lens[j] - 1, Colors[c]);
    fprintf(fp,"backbone 100\n");
    if (++c >= N_Colors) c=0;
    site = seq_starts[i+1] + 1;
    fprintf(fp,"select %d-%dB\n color %s\n",
        site, site + segment_lens[j] - 1, Colors[c]);
    fprintf(fp,"backbone 50\n");
    if (++c >= N_Colors) c=0;
  }
}

//*************************** dsp_typ *******************
void	GDenseSegPut(FILE *fp, dsp_typ gdsp)
{
	fprintf(fp,"dim = %d; numseg = %d\n",gdsp->dim,gdsp->numseg);
	for(Int4 i=0; i < gdsp->numseg; i++){
	   fprintf(fp,"len = %d; starts = { %d",gdsp->lens[i],gdsp->starts[i*gdsp->dim]);
	   for(Int4 d=1; d < gdsp->dim; d++){
		fprintf(fp,", %d",gdsp->starts[i*gdsp->dim+d]);
	   }
	   fprintf(fp," }\n");
	} 
	fprintf(fp,"\n");
}

dsp_typ GDenseSegNew(void)
{ dsp_typ	dsp; NEW(dsp,1,GDenseSeg); return dsp; }

dsp_typ GDenseSegFree(dsp_typ dsp)
{
   if(dsp){ 
	if(dsp->starts) free(dsp->starts); 
	if(dsp->lens) free(dsp->lens); 
	free(dsp); 
   } return NULL; 
}

unsigned char *GetSequenceWithGDenseSeg(dsp_typ dsp, Boolean query, Int4Ptr start, 
	Int4Ptr length)
// This needs to use my sequence structure to return sequence array.
{
        Int4		index, offset;
	e_type		E;

        if (dsp == NULL) return NULL;
        if(query == TRUE) { offset=0; E = dsp->qE; } else { offset=1; E = dsp->sE; }
        *start = dsp->starts[offset];
        *length = 0;
        for (index=0; index<dsp->numseg; index++){
                if (dsp->starts[offset+2*index] != -1)
                        *length += dsp->lens[index];
        } return (SeqPtr(E)+(*start)+1);  // WARNING: NOT SURE THIS IS CORRECT!!??
}

// Routines to convert SeqAlign to CMA output file...

void    SeqAlignToCMA(FILE *fp, sap_typ head, Int4 leftflank, Int4 rightflank,
	a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// Add leftflank & rightflank.
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
	UInt4	last_id=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  dsp_typ dsp=gsap->segs; 
	  if(dsp->subject_id != last_id){
		last_id= dsp->subject_id;
		fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE); 
		fprintf(fp,"\t\tlength = %d\n\n",LenSeq(dsp->sE));
	  } PutSeqAlignToCMA(fp,gsap,leftflank,rightflank,AB);
	} fprintf(fp,"\n");
}

void    PutSeqAlignToCMA(FILE *fp, sap_typ sap, Int4 leftflank, Int4 rightflank, a_type AB)
// print output to a file in cma format...
{
        Int4		r,i,s,q,qs,ss,length,index,total,m;
	UInt4	qid,sid;
	dsp_typ	dsp;
	e_type		qE,sE;
	Int4		qos,sos;
	char		*qstr,*sstr,c;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
	qid = dsp->query_id; qE = dsp->qE;
	sid = dsp->subject_id; sE = dsp->sE;
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);
        for(total=index=0; index<dsp->numseg; index++){
		total += dsp->lens[index];
	}
	NEW(qstr,total+3,char); NEW(sstr,total+3,char);
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1){
			for(i=0,r=qs+1; i < length; r++,i++){
			  qstr[q++] = AlphaChar(ResSeq(r,qE),AB); 
			}
		} else {	// insertion in subject sequence.
			for(i=0; i < length; i++){ qstr[q++] = '-'; }
		}
		if(ss != -1){
			for(i=0,r=ss+1; i < length; r++,i++){
			  sstr[s++] =c= AlphaChar(ResSeq(r,sE),AB);
if(qs == -1) c = tolower(c);  // insertion in subject sequence.
fprintf(fp,"%c",c);
			}
		} else { 	// deletion in subject sequence.
			for(i=0; i < length; i++){ 
			   sstr[s++] = '-'; 
fprintf(fp,"-");
			}
		} 
        } qs = dsp->starts[0] +1; ss = dsp->starts[1]+1; 
#if 0
	Int4 aln_len;
	// fprintf(fp,"\tScore = %.1f bits (%d), Expect = %1.2g\n",
	double percent_ids=PercentIdentGSAP(&aln_len, sap);
	fprintf(fp,"\tScore = %.0f bits (%d), Expect = %1.2g\n",
		sap->bit_score,sap->score,sap->evalue);
	fprintf(fp,"\tIdentities = %d/%d(%.0f)\n",
		(Int4)(percent_ids*aln_len),aln_len,percent_ids);
	for(i=0; i < total; i+=60){
	    end = MINIMUM(Int4,total,i+60);
	    fprintf(fp,"\n QUERY: %-4d ",qs+qos);
	    for(j=i; j < end; j++){ 
		fprintf(fp,"%c",qstr[j]);  
		if(qstr[j] != '-') qs++;
	    }
	    fprintf(fp," %d\n             ",qs-1+qos);
	    for(j=i; j < end; j++){ 
		if(qstr[j] == '-' || sstr[j] == '-') fprintf(fp," ");
		else {
		  q=AlphaCode(qstr[j],AB); s=AlphaCode(sstr[j],AB);
		  if(q==s) fprintf(fp,"%c",sstr[j]); 
		  else if(valAlphaR(q,s,AB) > 0) fprintf(fp,"+"); 
		  else fprintf(fp," ");
		}
	    }
	    fprintf(fp,"\n SBJCT: %-4d ",ss+sos);
	    for(j=i; j < end; j++){ 
		fprintf(fp,"%c",sstr[j]); 
		if(sstr[j] != '-') ss++;
	    }
	    fprintf(fp," %d\n",ss-1+sos);
	}
	fprintf(fp,"\n");
#endif
	free(qstr); free(sstr); 
}

Int4	DelimitGSeqAlignList(sap_typ head,Int4 **SS,Int4 **SP,Int4 **ES,Int4 **EP)
{
	UInt4	last_id=0;
	Int4		ssq,esq,spf,epf,N;
	Int4		*ss,*sp,*es,*ep,sq;
	sap_typ		gsap;
	
	N = NumSeqsListGSAP(head);
	NEW(ss,N+3,Int4); NEW(es,N+3,Int4); NEW(sp,N+3,Int4); NEW(ep,N+3,Int4);
	ssq=spf=INT4_MAX; esq=epf=0;
	for(sq=0,gsap=head; gsap!=NULL; gsap=gsap->next){
	  dsp_typ dsp=gsap->segs; 
	  if(dsp->subject_id != last_id){
		if(last_id > 0){
			ss[sq]=ssq; sp[sq]=spf; es[sq]=esq; ep[sq]=epf;
		}
		last_id= dsp->subject_id; sq++;
		ssq=spf=INT4_MAX; esq=epf=0;
	  } DelimitGSeqAlign(gsap,&ssq,&spf,&esq,&epf);
	} if(last_id > 0){ ss[sq]=ssq; sp[sq]=spf; es[sq]=esq; ep[sq]=epf; }
	*SS=ss; *SP=sp; *ES=es; *EP=ep;
	return N;
}

Int4	QueryStartGSeqAlnList(sap_typ sap)
{
	assert(sap != NULL && sap->segs != NULL);
	Int4	qs,start=INT4_MAX;
	for(sap_typ gsap=sap; gsap!=NULL; gsap=gsap->next){
		dsp_typ	dsp = sap->segs;
		qs = dsp->starts[0];
		if(qs < start) start = qs;
	} return start;
}

void    DelimitGSeqAlign(sap_typ sap, Int4 *SS,Int4 *SP,Int4 *ES,Int4 *EP)
{
        Int4		qs,ss,length,index,ssq,esq,spf,epf,tmp;
	UInt4	qid,sid;
	dsp_typ	dsp;
	e_type		qE,sE;

        assert(sap->segs != NULL);
	dsp = sap->segs; 
	qid = dsp->query_id; qE = dsp->qE;
	sid = dsp->subject_id; sE = dsp->sE;
        assert(qE && sE);
	ssq=spf=INT4_MAX; esq=epf=0;
        for(index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1){
			if(spf > qs) spf=qs;
			tmp=qs+length-1;
			if(epf < tmp) epf = tmp;
		}
		if(ss != -1){
			if(ssq > ss) ssq=ss;
			tmp=ss+length-1;
			if(esq < tmp) esq = tmp;
		}
        }
	if(*SS > ssq) *SS=ssq+1;
	if(*ES < esq) *ES=esq+1;
	if(*SP > spf) *SP=spf;
	if(*EP < epf) *EP=epf;
}

Int4    *PutScwrlSeqGSAP(FILE *fp,sap_typ sap, e_type keyE, a_type A)
// make a sequence for Scwrl homology modeling
// subject sequence corresponds to pdb file; query to concensus 
// and keyE to sequence to be modeled.
{
        Int4            r,i,j,qs,ss,length,r_s;
        e_type          qE,sE;
        char            c;
        Int4            *map;// map from keyE to pdbE ;

        assert(sap != NULL && sap->segs != NULL);
        dsp_typ dsp = sap->segs; qE = dsp->qE; sE = dsp->sE;
        assert(LenSeq(keyE) == LenSeq(qE));
        NEW(map,LenSeq(keyE)+3,Int4);
        ss = dsp->starts[1];
        for(r_s=1; r_s <= ss; r_s++){
                if(r_s % 50 == 0) fprintf(fp,"\n");
                c=AlphaChar(ResSeq(r_s,sE),A);
                fprintf(fp,"%c",tolower(c));
        }
        for(i=0; i <dsp->numseg; i++){
                qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
                length = dsp->lens[i];
                if(qs != -1 && ss != -1){ // in aligned region use keyE.
                   r=qs+1; r_s=ss+1;
                   for(j=0; j < length; j++){
                        if(r_s % 50 == 0) fprintf(fp,"\n");
                        fprintf(fp,"%c",AlphaChar(ResSeq(r,keyE),A));
                        map[r]=r_s;
                        r++; r_s++;
                   }
                } else if(ss != -1){    // put subject pdb in gaps.
                   r_s=ss+1;
                   for(j=0; j < length; j++){
                     if(r_s % 50 == 0) fprintf(fp,"\n");
                     c=AlphaChar(ResSeq(r_s,sE),A);
                     fprintf(fp,"%c",tolower(c));
                     r_s++;
                   }
                }
        }
        for(   ;r_s <= LenSeq(sE); r_s++){
                if(r_s % 50 == 0) fprintf(fp,"\n");
                c=AlphaChar(ResSeq(r_s,sE),A);
                fprintf(fp,"%c",tolower(c));
        } fprintf(fp,"\n");
        return map;
}

Int4    *PutSCGenSeqGSAP(FILE *fp,sap_typ sap, e_type keyE, a_type A)
// make a cobbled sequence for SCGEN homology modeling
// subject sequence corresponds to pdb file; query to concensus 
// and keyE to sequence to be modeled.
{
        Int4            r,i,j,qs,ss,length,r_s;
        e_type          qE,sE;
        char            c;
        Int4            *map;// map from keyE to pdbE ;

        assert(sap != NULL && sap->segs != NULL);
        dsp_typ dsp = sap->segs; qE = dsp->qE; sE = dsp->sE;
        assert(LenSeq(keyE) == LenSeq(qE));
        NEW(map,LenSeq(keyE)+3,Int4);
        ss = dsp->starts[1];
	// 1. Output original sequence on first line...
        for(r_s=1; r_s <= LenSeq(sE); r_s++){
	     c=AlphaChar(ResSeq(r_s,sE),A); fprintf(fp,"%c",c);
	} fprintf(fp,"\n");
	// 2. Output cobbled sequence on second line...
        for(r_s=1; r_s <= ss; r_s++){
                c=AlphaChar(ResSeq(r_s,sE),A); fprintf(fp,"%c",c);
        }
        for(i=0; i <dsp->numseg; i++){
                qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
                length = dsp->lens[i];
                if(qs != -1 && ss != -1){ // in aligned region use keyE.
                   r=qs+1; r_s=ss+1;
                   for(j=0; j < length; j++){
                        fprintf(fp,"%c",AlphaChar(ResSeq(r,keyE),A));
                        map[r]=r_s; r++; r_s++;
                   }
                } else if(ss != -1){    // put subject pdb in gaps.
                   r_s=ss+1;
                   for(j=0; j < length; j++){
                     c=AlphaChar(ResSeq(r_s,sE),A);
                     fprintf(fp,"%c",c); r_s++;
                   }
                }
        }
        for(   ;r_s <= LenSeq(sE); r_s++){
                c=AlphaChar(ResSeq(r_s,sE),A); fprintf(fp,"%c",c);
        } fprintf(fp,"\n");
        return map;
}

static h_type gHG=0;

BooLean	NumberOfGapsGSAP(h_type HG, sap_typ sap)
// print output to a file in cma format...
{
        Int4            s,qs,ss,length,index;
        dsp_typ    	dsp;
        e_type          qE,sE;
	Int4		pos;

        assert(sap->segs != NULL);
        dsp = sap->segs;
        qE = dsp->qE; sE = dsp->sE;
        for(index=0; index<dsp->numseg; index++){
            qs = dsp->starts[2*index];
            ss = dsp->starts[2*index+1];
            length = dsp->lens[index];
            if(qs != -1 && ss != -1){   // matching region
		pos = qs + length - 1;
            } else if(qs == -1){        // insertion(s) in subject sequence.
		IncdMHist(pos,length,HG); 
            } else if(ss == -1){        // deletion(s) in subject sequence.
		for(pos=qs,s=0; s < length; s++){ IncdMHist(pos,1,HG); pos++; } 
            } else {        // insertion in query or subject sequence.
		print_error("This should not happen!");
            }
	} return TRUE;
}

Int4    DriverNumberOfGapsGSAP(FILE *fp,sap_typ head)
// return the number of sequences in list...
{
        UInt4   n=0,last_id=0;

	h_type HG=Histogram("number of gaps",0,LenSeq(head->segs->qE),1.0); 
        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
		dsp_typ dsp=gsap->segs;
                if(dsp->subject_id != last_id){
                        last_id=dsp->subject_id; n++;
                }
		NumberOfGapsGSAP(HG,gsap);
			//PutGSeqAlign(fp,gsap,60,A);
        } 
	PutHist(fp,60,HG); NilHist(HG);
	return n;
}

BooLean	DeleteGDenseSeg(BooLean first, dsp_typ dsp)
{
	Int4	i,index;

	if(dsp->numseg <= 2) return FALSE; // don't delete them all...
	if(first){	// then delete the first segment.
	   dsp->numseg = dsp->numseg - 2;
	 for(i=2,index=0; index<dsp->numseg; index++,i++){
	   dsp->starts[2*index] = dsp->starts[2*i];
           dsp->starts[2*index+1] = dsp->starts[2*i+1];
	   dsp->lens[index] = dsp->lens[i];
	 }  
	} else {
	   dsp->numseg = dsp->numseg - 2;	// don't need to reinitialize
	}
	return TRUE;
}

BooLean StartEndScoreGSAP(sap_typ sap, Int4 *Start, Int4 *End, Int4 **mx, a_type A)
// compute and return the starting and ending subsegment scores.
{
        Int4            I,J,r,i,j,s,qs,ss,length,index,subscore,end;
        e_type          qE,sE;
        Int4            ins_open=11,ins_extend=1; // fix this later!!
        Int4            del_seg[2];

        assert(sap->segs != NULL);
        dsp_typ dsp= sap->segs;
        if(dsp->numseg < 3) return FALSE;
        qE = dsp->qE; sE = dsp->sE;
        // 1. first check the score
        for(end=0,index=0; end <= 1; index = dsp->numseg-1,end++){
          for(subscore=0,I=index,J=0; J <= 1 ; I++,J++){
            qs = dsp->starts[2*I];
            ss = dsp->starts[2*I+1];
            length = dsp->lens[I];
            if(qs != -1 && ss != -1){   // matching region
                for(i=0,j=qs,s=ss+1; i < length; s++,j++,i++){
                  r=ResSeq(s,sE);
                  if(mx) subscore += mx[j][r];
                  else subscore += valAlphaR(r,ResSeq(j+1,qE),A);
                  if(mx) fprintf(stderr,"%d: %c %c = %d\n",
                                j,AlphaChar(ResSeq(j+1,qE),A),AlphaChar(r,A),subscore);
                }
                // fprintf(stderr,"Matching region score = %d\n",subscore);
            } else {        // insertion in query or subject sequence.
                subscore -= (ins_extend*length + ins_open);
            }
          } del_seg[end]=subscore;
        } *Start=del_seg[0]; *End=del_seg[1];
        return TRUE;
}

BooLean DeleteOverlapGSAP(sap_typ sap1,sap_typ sap2, Int4 **mx,a_type A)
// Assume that sap1 tail overlaps with sap2 head.
// NEW routine...
{
        Int4    dummy,s2_score,e1_score;
        BooLean multiple1,multiple2;
        multiple1=StartEndScoreGSAP(sap1,&dummy,&e1_score,mx,A);
        multiple2=StartEndScoreGSAP(sap2,&s2_score,&dummy,mx,A);
        if(multiple1 && multiple2){     // either end can be delete
                if(e1_score > s2_score) DeleteGDenseSeg(TRUE,sap2->segs);
                else DeleteGDenseSeg(FALSE,sap1->segs);
        } else if(multiple2){   // only head of sap2 can be deleted
                DeleteGDenseSeg(TRUE,sap2->segs);
        } else if(multiple1){   // only tail of sap1 can be deleted
                DeleteGDenseSeg(FALSE,sap1->segs);
        } else return FALSE;            // neither can be deleted.
        return TRUE;
}

sap_typ AlignToGSAP(e_type qE, e_type sE, char *operation, Int4 start_qE, Int4 start_sE)
// Written by Aleksandar.
{
        Int4 i=1,k,j,s,t,seg;
        Int4 *starts, *lens;
        short numseg = 0;
        while(operation[i] != 'E'){
                if(tolower(operation[i]) != tolower(operation[i-1])) numseg++;
                i++;
        }
        MEW(starts,2*numseg+1,Int4); MEW(lens,numseg+2,Int4);
	for(i=2,k=1,seg=0; operation[i] != 'E'; i++){
                if(tolower(operation[i]) != tolower(operation[i-1])) {
                        lens[seg++] = k;
                        k = 1;
                } else { k++; }
        }
        lens[seg] = k;
        for(k=1,s=1,t=1,j=0,i=0;i<numseg;i++){
                switch(operation[t]){
                        case 'd':
                        case 'D': starts[j++] = k+start_qE-2;
                                  starts[j++] = -1; k += lens[i];
                                  break;
                        case 'm':
                        case 'M': starts[j++] = k+start_qE-2;
                                  starts[j++] = s+start_sE-2;
                                  k += lens[i]; s += lens[i];
                                  break;
                        case 'I': starts[j++] = -1;
                                  starts[j++] = s+start_sE-2;
                                  s += lens[i];
                                  break;
                } t += lens[i];
        }
#if 0
	for(i=0;i<=2*numseg-1;i++) printf("starts[%d]=%d\n",i,starts[i]);
	for(i=0;i<numseg;i++) printf("lens[%d]=%d\n",i,lens[i]);
#endif
        return MakeGSeqAlign(numseg, SeqI(qE), SeqI(sE), qE, sE, starts, lens);
}

Int4	ComputeQueryScoreDSP(FILE *fp,e_type qE,Int4 **mx)
{
	Int4            r,j,s,s0,w=3,n,start,end;
        double		score,*vals;
	dh_type		dH;

	NEW(vals,2*w+4,double);
	h_type HG; Int4 pos=1;
	if(fp) HG=Histogram("avg. score",0,LenSeq(qE),1.0); 
	for(j=0,s=1; s <= LenSeq(qE); j++,s++){
           score = 0;
	   end=MINIMUM(Int4,s+w,LenSeq(qE));
	   start=MAXIMUM(Int4,s-w,1);
	   dH=dheap(2*w+3,4);
	   for(n=0,j=1,s0=start; s0 <= end; s0++,j++){
              r=ResSeq(s0,qE);
	      vals[j] = (mx[s0-1][r]);
              // score += mx[s0-1][r];
	      insrtHeap(j,-((keytyp)vals[j]),dH);
	      // insrtHeap(j,((keytyp)vals[j]),dH);
	      n++;
	   } 
           score = 0;
	   for(s0=0; s0 <= w; s0++){
		j = delminHeap(dH); score += vals[j];
	   } while(delminHeap(dH) != 0) ; 
	   Nildheap(dH);
	   if(fp){ IncdMHist(pos,(double)score*100.0/(double)(w+1),HG); pos++; }
	} free(vals);
	if(fp){ 
		PutHist(stderr,60,HG); NilHist(HG); 
        	fprintf(fp,"computed score = %d\n",(Int4)score);
	} return score;
}

Int4	ComputeScoreDSP(FILE *fp,dsp_typ dsp,Int4 **mx)
{
	Int4            r,i,j,s,qs,ss,length,index,total;
        e_type          qE,sE;
        Int4            score,subscore;
        Int4            ins_open=11,ins_extend=1; // fix this later!!

        assert(dsp != NULL);
        qE = dsp->qE; sE = dsp->sE;
        for(total=index=0; index<dsp->numseg; index++){
                total += dsp->lens[index];
        }
	h_type HG;
if(fp) HG=Histogram("running score",0,total,1.0); 
Int4 pos=1;
        score = 0;
        for(index=0; index<dsp->numseg; index++){
                qs = dsp->starts[2*index];
                ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
                subscore=0;
                if(qs != -1 && ss != -1){   // matching region
                        for(i=0,j=qs,s=ss+1; i < length; s++,j++,i++){
                          r=ResSeq(s,sE);
                          subscore += mx[j][r];
                          score += mx[j][r];
if(fp){ IncdMHist(pos,score,HG); pos++; }
// IncdMHist(pos,subscore,gHG);
                        }
                } else {        // insertion in query or subject sequence.
                        subscore = -(ins_extend*length + ins_open);
			score += subscore;
if(fp) for(i=0,j=qs; i < length; j++,i++) { IncdMHist(pos,score,HG); pos++;  }
                }
if(fp){
		fprintf(fp," subscore(%d) = %d\n",index+1,subscore);
                fprintf(fp," running score = %d\n",score);
}
        } 
if(fp){ PutHist(stderr,60,HG); NilHist(HG); }
        if(fp) fprintf(fp,"computed score = %d\n",score);
        return score;
}

Int4	ComputeScoreDSP(dsp_typ dsp,Int4 **mx){ return ComputeScoreDSP(0,dsp,mx); }

Int4 MaxSegScoreDSP(dsp_typ dsp, Int4 **mtx, e_type sE, Int4 *sp_nmbr)
// Derived from Aleksandar's HMM-BLAST routines...
// finds max scoring sp; carry over to my_ncbi.c
{
        Int4 i,j,temp,score=0,startQ,startS;
        unsigned char *ptr = SeqPtr(sE);

        Int4 *starts = dsp->starts;
        Int4 *lens = dsp->lens;
        *sp_nmbr = 1;
        for(i=1;i<=dsp->numseg;i++) {
                startQ = starts[2*i-2];
                if(startQ == -1) continue;
                startS = starts[2*i-1];
                if(startS == -1) continue;
                for(temp=0,j=1;j<=lens[i-1];j++) {
                        temp += mtx[startQ+j-1][ptr[startS+j-1]];
                }
                if(temp > score) { score = temp; *sp_nmbr = i; }
        }
        return score;
}

Int4 MaxWordScoreDSP(dsp_typ dsp, Int4 sp_nmbr, Int4 word_len, Int4 **mtx, 
		Int4 mtx_len, e_type sE)
// Derived from Aleksandar's HMM-BLAST routines...
{
        Int4 i,j,temp=0,score=0,first,last,shortWord;
	Int4	startQ,startS;
        unsigned char *ptr = SeqPtr(sE);
        Int4 *starts = dsp->starts;
        Int4 *lens = dsp->lens;

        startQ = starts[2*sp_nmbr-2];
        startS = starts[2*sp_nmbr-1];
        Int4 tmtx = mtx_len-startQ+1;
        Int4 tseq = LenSeq(sE)-startS+1;

        if(tmtx < word_len || tseq < word_len) return 0;
        for(j=1;j<=word_len;j++) {
                temp += mtx[startQ+j-1][ptr[startS+j-1]];
        }
        if(temp > score) { score = temp; }
        for(i=2;i<=lens[sp_nmbr-1]-word_len+1;i++) {
                first = mtx[startQ+i-2][ptr[startS+i-2]];
                last = mtx[startQ+i+word_len-2][ptr[startS+i+word_len-2]];
                temp += last - first;
                if(temp > score) { score = temp; }
        }
        return score;
}

BooLean	FixAlnRunOverGSAP0(FILE *fp, sap_typ sap, Int4 **mx, Int4 minsubscore,
	a_type A)
// remove subsegments in the sap below minsubscore.
{
        Int4            r,i,j,s,qs,ss,length,index;
        dsp_typ    	dsp;
        e_type          qE,sE;
        Int4            score,subscore,end;
	Int4		ins_open=11,ins_extend=1; // fix this later!!
	Int4		del_seg[2];

        assert(sap->segs != NULL);
        dsp = sap->segs;
	if(dsp->numseg < 3) return FALSE;
        qE = dsp->qE; sE = dsp->sE;
        score = 0;
	// 1. first check the score 
        for(end=0,index=0; end <= 1; index = dsp->numseg-2,end++){
	  subscore=0;
          for(Int4 I=index,J=0; J <= 1 ; I++,J++){
            qs = dsp->starts[2*I];
            ss = dsp->starts[2*I+1];
            length = dsp->lens[I];
            if(qs != -1 && ss != -1){   // matching region
                for(i=0,j=qs,s=ss+1; i < length; s++,j++,i++){
		  r=ResSeq(s,sE);
		  if(mx) subscore += mx[j][r];
		  else {
			Int4 qr=ResSeq(j+1,qE);
			subscore += valAlphaR(r,qr,A);
		  }
		}
            } else {        // insertion in query or subject sequence.
		subscore -= (ins_extend*length + ins_open);
            }
	  }
	  if(subscore < minsubscore) del_seg[end]=subscore;
	  else del_seg[end]=0;
	}
if(fp){
	if(del_seg[0]) fprintf(fp,"1st subscore = %d\n",del_seg[0]);
	if(del_seg[1]) fprintf(fp,"2nd subscore = %d\n",del_seg[1]);
if(mx){
  if(dsp->qE == dsp->sE) ComputeScoreDSP(stderr,dsp,mx);
	else ComputeScoreDSP(dsp,mx);
  }
}
	if(del_seg[0]) del_seg[0]=DeleteGDenseSeg(TRUE,dsp);
	if(del_seg[1]) del_seg[1]=DeleteGDenseSeg(FALSE,dsp);
	// fprintf(fp,"computed score = %d\n",score);
	return (del_seg[0] || del_seg[1]);
}

Int4    FixAlnRunOverGSAP(sap_typ head,Int4 **mx,a_type A,Int4 minfix)
// return the number of sequences in list...
{
	FILE	*fp=stderr;
        UInt4   n=0,last_id=0;

gHG=Histogram("subscores",0,500,2.0); 
        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
                dsp_typ dsp=gsap->segs;
                if(dsp->subject_id != last_id){
                        last_id=dsp->subject_id; n++;
if(dsp->subject_id == 120 || dsp->subject_id == 77 || dsp->subject_id == 20){
			fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE);
			fprintf(fp,"\t\tlength = %d\n\n",LenSeq(dsp->sE));
}
                }
			//PutGSeqAlign(fp,gsap,60,A);
if(dsp->subject_id == 120 || dsp->subject_id == 77 || dsp->subject_id == 20){
			PutGSeqAlign(fp,gsap,60,A);
		}
		while(FixAlnRunOverGSAP0(fp,gsap,mx,minfix,A)) ;
#if 1
if(dsp->subject_id == 120 || dsp->subject_id == 77 || dsp->subject_id == 20){
			PutGSeqAlign(fp,gsap,60,A);
		}
#endif
        } 
PutHist(stderr,60,gHG); NilHist(gHG);
	return n;
}

static Int4	CountNumSegs(register char *hmm_path)
{
        register Int4	i,numseg = 0;
	register char	c,state;

	if(hmm_path[0] != 'B') print_error("CountNumSegs() input error");
        for(state='B',i=1; (c=hmm_path[i]); i++){
             switch(c){
		case 'm': if(state!='m') numseg++; break;
		case 'd': if(state!='d') numseg++; break;
		case 'i': if(state!='i') numseg++; break;
		case 'E': return numseg; break;
		case 'n':	// N-terminal insertions...ignore.
		case 'c':	// C-terminal insertions...ignore.
		default: print_error("CountNumSegs() input error");
		   break;
	     } state=c;
		
	} print_error("CountNumSegs() input error");
}

sap_typ PathHMMToGSAP(e_type hmmE, e_type sE, char *hmm_path, double *Evalue,
	double *bit_score,Int4 score,Int4 *scores)
// Create a sap_typ from a path through an HMM.
// Uses Eddy's Plan7 architecture
// State types: 'S','n','B','m','d','i','E','c','j','T'.
// Uses a finite state machine to parse the states.
// Input string ends in 0.
{
        Int4	i,j,seg,Nrpts,rpt,numseg = 0;
	char	state,c;
	Int4    *starts, *lens;
	Int4	qs,ss,L;
	sap_typ	head=0,sap=0;

	if(hmm_path[0] != 'S') print_error("PathHMMToGSAP() parse error 1");
        for(Nrpts=0,i=1; (c=hmm_path[i]); i++){
		if(c=='B') Nrpts++;
	} if(Nrpts==0) return NULL;
	starts=lens=0;
	
        for(j=seg=qs=ss=rpt=0,state='S',i=1; (c=hmm_path[i]); i++){
             switch(c){
		case 'n':	// N-terminal insertions..
		   if(!strchr(" nS",state)) print_error("PathHMMToGSAP() parse error 2");
		   ss++; 
		   break;
		case 'B':	// Begin state; encountered a repeat...
		   if(!strchr(" njES",state)){
	     		fprintf(stderr,"state=%c; c = %c\npath=%s\n",state,c,hmm_path);
			print_error("PathHMMToGSAP() parse error 3");
		   }
		   rpt++;
		   numseg = CountNumSegs(hmm_path+i);
		   // fprintf(stderr,"rpt(%d) = %d segs\n",rpt,numseg);
		   if(numseg==0) print_error("PathHMMToGSAP() parse error 4");
		   NEW(starts,2*numseg+2,Int4);
		   NEW(lens,numseg+2,Int4);
		   qs=0; seg=j=L=0;
		   break;
		case 'm':	// match state; only reachable from 'B','m','i' or 'd'.
		   if(!strchr(" dimB",state)) print_error("PathHMMToGSAP() parse error 5");
		   if(state!='m'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = ss; seg++; starts[seg] = qs; seg++;
			L=0;
		   } qs++; ss++; L++;
		   break;
		case 'd':	// delete state; only reachable from 'd', 'm' or 'B'.
		   if(!strchr(" dmB",state)) { 
			fprintf(stderr,"current state = %c; last state = %c\n",
				c,state);
			print_error("PathHMMToGSAP() parse error 6");
		   }
		   if(state!='d'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = -1; seg++; starts[seg] = qs; seg++;
			L=0;
		   } qs++; L++;
		   break;
		case 'i':	// insert state; only reachable from 'i' or 'm'.
		   if(!strchr(" im",state)) print_error("PathHMMToGSAP() parse error 7");
		   if(state!='i'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = ss; seg++; starts[seg] = -1; seg++;
			L=0;
		   } ss++; L++;
		   break;
		case 'E':	// End state; only reachable from 'd',or 'm'.
		   if(qs != LenSeq(hmmE)) print_error("PathHMMToGSAP() seq length error 8");
		   if(!strchr(" dm",state)) print_error("PathHMMToGSAP() parse error 9");
		   if(seg > 0){ lens[j]=L; }
		   if(head){
			// sap->next= MakeGSeqAlign(numseg,SeqI(hmmE),SeqI(sE),hmmE,sE,starts,lens);
			sap->next= MakeGSeqAlign(numseg,SeqI(sE),SeqI(hmmE),sE,hmmE,starts,lens);
			sap=sap->next;
			sap->score = scores[rpt];
			sap->evalue = Evalue[rpt];
			sap->bit_score = bit_score[rpt];
		   } else {
			// head=sap=MakeGSeqAlign(numseg,SeqI(hmmE),SeqI(sE),hmmE,sE,starts,lens);
			head=sap=MakeGSeqAlign(numseg,SeqI(sE),SeqI(hmmE),sE,hmmE,starts,lens);
			sap->score = scores[rpt];
			sap->evalue = Evalue[rpt];
			sap->bit_score = bit_score[rpt];
		   } qs=0;
		   break;
		case 'j':	// repeat insert state; only reachable from 'j' or 'E'.
		   if(!strchr(" jE",state)) print_error("PathHMMToGSAP() parse error 10");
		   ss++;
		   break;
		case 'T':
		   if(!strchr(" cE",state)) print_error("PathHMMToGSAP() parse error 11");
        	   return head;
		   break;
		case 'c':	// Stop here and return head.
		   if(!strchr(" cE",state)) print_error("PathHMMToGSAP() parse error 12");
		   ss++;
        	   return head;
		   break;
		default: print_error("PathHMMToGSAP() parse error 13");
		   break;
	     }
	     if(ss > LenSeq(sE)) print_error("PathHMMToGSAP() subject seq error 14");
	     if(qs > LenSeq(hmmE)) print_error("PathHMMToGSAP() hmm seq error 15");
	     // fprintf(stderr,"state=%c; c = %c\n",state,c);
	     state=c;
        } print_error("PathHMMToGSAP() parse error 16");
}

Int4    GetStartEndGSAP(sap_typ sap, UInt4 *Start, UInt4 *End)
{
        Int4            r,i,j,s,length,index,total;
	UInt4	ss,es;

        assert(sap->segs != NULL);
        dsp_typ dsp = sap->segs;
        for(total=index=0; index<dsp->numseg; index++){
                total += dsp->lens[index];
        }
        ss = dsp->starts[1]+1;
        for(s=index=0; index<dsp->numseg; index++){
                if(dsp->starts[2*index+1] != -1) s += dsp->lens[index];
        } es = ss + s -1;
        *Start=ss; *End=es;
        return (es-ss+1);
}

void    PutDeleteBestHSPs(FILE *fp, sap_typ head, a_type AB)
// Output sequences hit with best HSP deleted...
{
        UInt4   last_id=0,Start,End;
        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	  // 1. Get start and end of hsp:
	  Int4    len = GetStartEndGSAP(gsap,&Start,&End);
          fprintf(stderr,"Start = %d; End = %d; len = %d\n",Start,End,len);
          dsp_typ dsp=gsap->segs;
          if(dsp->subject_id != last_id){	// delete only Int4est HSP
                last_id= dsp->subject_id;
                // fprintf(fp,"\n>"); PutSeqInfo(fp,dsp->sE);
                // fprintf(fp,"\t\tlength = %d\n\n",LenSeq(dsp->sE));
		PutDeleteSeq(fp,Start,End,dsp->sE,AB);
	  }
        } 
}

void	GetSeqAlignTable(unsigned short *TabP, unsigned char *TabR, sap_typ sap,
	a_type AB)
// print output to a file in cma format...
// query is used as a template to define length and dimensions of cma.
// WARNING: OFFSETS ARE NOT SET CORRECTLY!!!
{
        Int4		s,q,qs,ss,length,index;
	Int4		qos,sos,m;
	unsigned short	r;

        assert(sap->segs != NULL);
	dsp_typ	dsp = sap->segs; 
	e_type	qE = dsp->qE, sE = dsp->sE;

	Int4	i,n;
	Int4 len = LenSeq(dsp->sE);
	qos = OffSetSeq(qE); sos = OffSetSeq(sE);
        for(m=s=q=index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs != -1 && ss != -1){	// match region...
fprintf(stderr,"ss=%d; qs=%d; len=%d\n",ss,qs,length);
		    for(i=0,r=ss+1,s=qs+1; i < length; r++,s++,i++){
			TabP[s] = r+sos;
			TabR[s] = ResSeq(r,sE);
		    } 
fprintf(stderr,"TabP[qs+1]=%d\n",TabP[qs+1]);
		} else if(qs == -1){		// insertion in subject sequence.
		} else if(ss == -1){		// deletion in subject sequence.
		} else print_error("GetSeqAlignTable( ): dsp error");
        }
}


void    SeqAlignToTable(FILE *fp,sap_typ head, a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
        if(head==0) return;

        Int4	r,i,j,n = NumHSPsListGSAP(head);
	if(n > 30) print_error("Too many aligned sequences for tabular output");

	// Get sequence information and print out first line of table.
	unsigned short	**TabP;
	unsigned char	**TabR;
	sap_typ		gsap;
	e_type		qE,sE;

	qE = QuerySeqGSAP(head);
	NEWP(TabP,n+3,unsigned short);
	NEWP(TabR,n+3,unsigned char);
	PutShortSeqID(fp,qE); fprintf(fp,"\t");
        for(i=1,gsap=head; gsap!=NULL; gsap=gsap->next,i++){
		sE = SubjectSeqGSAP(gsap);
		PutShortSeqID(fp,sE); fprintf(fp,"\t");
		NEW(TabP[i],LenSeq(qE)+3,unsigned short);
		NEW(TabR[i],LenSeq(qE)+3,unsigned char);
          	GetSeqAlignTable(TabP[i],TabR[i],gsap,AB);
	} fprintf(fp,"\n");

	// Output table information.
	for(j=1; j <= LenSeq(qE); j++){
	  r = ResSeq(j,qE);
	  fprintf(fp,"%c%d\t",AlphaChar(r,AB),j);
	  for(i=1; i <= n; i++){
		if(TabP[i][j] != 0)
		  fprintf(fp,"%c%d\t",AlphaChar(TabR[i][j],AB),TabP[i][j]);
		else fprintf(fp,"\t");
	  } fprintf(fp,"\n");
	} fprintf(fp,"\n");

        for(i=1,gsap=head; gsap!=NULL; gsap=gsap->next,i++){
	  free(TabP[i]); free(TabR[i]); 
        } free(TabP); free(TabR); 
}

sap_typ TrimSAP(sap_typ sap, Int4 trimleft, Int4 trimright)
//trims trimleft residues in subj seq from sap at the begining and
//trimright at the end
{ 
	e_type qE = QuerySeqGSAP(sap);
	e_type sE = SubjectSeqGSAP(sap);
	dsp_typ dsp = sap->segs;
	Int4 numseg = dsp->numseg;
	Int4 *starts= dsp->starts;
	Int4 *lens = dsp->lens;
	Int4 len_s = starts[2*numseg-1]+lens[numseg-1]-1;
	assert(trimleft >=0 && trimright >= 0 && trimleft + trimright <= len_s);
	if(trimleft + trimright == len_s) return NULL;
	Int4 newNumseg, block_left, block_right, offset_left=0, offset_right=0, i, j;
	Int4 x,y;
	if(trimleft == 0) {
		block_left = 0; offset_left = 0;
	} else {	
	   for(i=0;i<numseg;i++) {
		if((x=starts[2*i+1]+lens[i]) >= (y=starts[1]+trimleft)) break;
	   }
	   if((x==y) && (i != numseg-1) && starts[2*i] != -1 && starts[2*i+1] != -1) { 
		i += 2; block_left = i; offset_left = 0; 
	   } else if (starts[2*i] == -1 || starts[2*i+1] == -1) {
			i++; block_left = i; offset_left = 0;
	   } else {
		block_left = i;
		for(j=1; j<=lens[i];j++) {
		  if(starts[2*i+1]+j >= starts[1]+trimleft) { offset_left = j; break; }
		}
	   }
	}
	if(trimright == 0) { block_right = numseg-1; offset_right = 0; }
	else {
	   for(i=0;i<=numseg-1;i++) {
		if((x=starts[2*i+1]+lens[i]) >= 
			(y=starts[2*numseg-1]+lens[numseg-1]-trimright)) break; 
	   }
	   if(starts[2*i] == -1 || starts[2*i+1] == -1){
			i--; block_right = i;  offset_right = 0;
	   } else {
		for(block_right=i,j=0;j<lens[i];j++) {
		   if(starts[2*i+1]+lens[i]-j <= 
			starts[2*numseg-1] + lens[numseg-1] - trimright) {
				offset_right = j; break;
				}
		}
	   }
	}
	newNumseg = block_right - block_left + 1;
	sap_typ newSAP = GSeqAlignNew();
        dsp_typ newDSP = GDenseSegNew();
        newDSP->dim = 2; newDSP->numseg = newNumseg;
        newDSP->subject_id = SeqI(sE); newDSP->sE=sE;
        newDSP->query_id=0; newDSP->qE=qE;
	Int4 *newLens,*newStarts;
	NEW(newLens,newNumseg,Int4);
	NEW(newStarts,2*newNumseg,Int4);
	if(block_left == block_right) { 
		newLens[0] = lens[block_left] - offset_left - offset_right;
	} else { 
		newLens[0] = lens[block_left] - offset_left; 
		newLens[newNumseg-1] = lens[block_right] - offset_right; 
		for(i=1;i<newNumseg-1;i++) newLens[i] = lens[block_left+i]; 
	}
	newStarts[0] = starts[2*block_left] + offset_left;
	newStarts[1] = starts[2*block_left+1] + offset_left;
	for(i=1;i<=newNumseg-1;i++) { 
		newStarts[2*i] = starts[2*(block_left+i)]; 
		newStarts[2*i+1] = starts[2*(block_left+i)+1];
	}
        newDSP->starts = newStarts; newDSP->lens = newLens;
        newSAP->segs = newDSP; newSAP->next = NULL;
        newSAP->dim =2;
	return newSAP;
}

void	MkGlobalSeqAlign(sap_typ sap,a_type AB)
// make the Subject globally aligned with query...
// WARNING: This assumes that sap is for a mapgaps alignment; will cause major problems for MCMC alignments.
{
	Int4		r,i,j,s,e,q,qs,ss,len,length,index,total,m,start,end;
	Int4		qs1,qs2,ss1,ss2,len1,len2;
	Int4		lenQ,lenS;
	dsp_typ		dsp = sap->segs; 
	e_type		qE=dsp->qE,sE=dsp->sE;
	
	// if(dsp->numseg < 2) return;

	// 1. first check last segment (easiest to fix).
	i=dsp->numseg - 1;	// last segment.
	qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1]; len = dsp->lens[i];
	// qs = 0...lenQ-1; ss = 0...lenS-1
	lenQ=LenSeq(qE)-qs; // length of full query seq.
	lenS=LenSeq(sE)-ss; // length of full query seq.
	if(lenQ > len && lenS > len)
	{
		// fprintf(stderr,"FOUND DELETION ON END\n");
		if(lenS >= lenQ) dsp->lens[i] = lenQ;
		else dsp->lens[i] = lenS;
		// PutGSeqAlignList(stderr, sap, 60, AB);
	}

	// 2. next check first segment.
	i=0;	// first segment.
	qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1]; len = dsp->lens[i];
	if(qs > 0 && ss > 0){
		// fprintf(stderr,"FOUND DELETION AT BEGIN\n");
		if(ss >= qs){	// e.g., if ss = 12, qs = 5 then lengthen by qs = 5.
			dsp->starts[2*i] = 0;
			dsp->starts[2*i+1] = ss - qs;	// e.g., ss = 12 - 5 = 7.
			dsp->lens[i] = len + qs;	// e.g., len = len + 5.
		} else {	// e.g., if ss = 5, qs = 12 then lengthen by ss = 5
			dsp->starts[2*i] = qs - ss;	// e.g., qs = 12 - 5 = 7.
			dsp->starts[2*i+1] = 0;  
			dsp->lens[i] = len + ss;
		}
		// PutGSeqAlignList(stderr, sap, 60, AB);
	}
#if 0
	i=0; qs1 = dsp->starts[2*i]; ss1 = dsp->starts[2*i+1]; len1 = dsp->lens[i];
	i=1; qs2 = dsp->starts[2*i]; ss2 = dsp->starts[2*i+1]; len2 = dsp->lens[i];
	if(ss1 == -1 && qs1 != -1 && ss2 != -1 && qs2 != -1){
fprintf(stderr,"FOUND DELETION ON END\n");
		// if ss is deleted at end and matches at next segment.
		// then fix this as much as possible.
	    start = ss2+1;  // start of subject seq. in alignment
	    if(len1 <= start){  
		// if the first segment can be filled with the subject seq.
		// then collapse first two segments into one.
		i=0; qs1 = dsp->starts[2*i]; ss1 = dsp->starts[2*i+1]; len1 = dsp->lens[2*i];
		j=1; dsp->starts[2*j]=dsp->starts[2*i];  // query start = beginning of alignment.
		     dsp->starts[2*j+1]=start-len1;    // back up subject start by len1.
		     dsp->lens[j] = len1+len2;
		
		// delete first segment.
        	for(i=0,j=1; j < dsp->numseg; i++,j++){
		  dsp->starts[2*i] = dsp->starts[2*j] ;
		  dsp->starts[2*i+1] = dsp->starts[2*j+1] ;
		  dsp->lens[i] = dsp->lens[j];
		} dsp->numseg--;
	    } else if(start > 1){
	    }
	} 	// else do nothing...
#endif
}

void	MakeGlobalSeqAlign(sap_typ head,a_type AB)
{
	Int4 i=0;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
		// fprintf(stderr,"Globalizing align %d\n",++i);
		MkGlobalSeqAlign(gsap,AB);
	}
}

