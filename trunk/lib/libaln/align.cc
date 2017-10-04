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

#include "align.h"

aln_typ	CreateAlign(ss_type P, Int4 length, s_type *S)
{
	aln_typ	Aln;
	Int4     end,i,j,I,s,n,e,N;
	unsigned char	*seq;
	e_type  E;

	NEW(Aln,1,align_type);
	Aln->P = P;
	Aln->A = SeqSetA(P);
	Aln->refined = FALSE;
	Aln->leng = length;
	for(n=0; S[n]!=NULL; n++) ;
	Aln->S = S;
	N = Aln->N = MAXIMUM(Int4,8*n,NSeqsSeqSet(Aln->P));
	Aln->H = Mheap(N,3); 
	Aln->scores = NULL; 
	Aln->n = n;
        for(Aln->nseg=0,I=1 ;I <= NSeqsSeqSet(Aln->P);I++) {
		Aln->nseg += MAXIMUM(Int4,0,SqLenSeqSet(I,Aln->P) - length + 1);
        }
	Aln->M = MkProfile(Aln->nseg,length,tFreqSeqSet(Aln->P),Aln->A);
	for(i=0; S[i]!=NULL; i++){
 	    I=SegmentI(S[i]); E = SeqSetE(I,Aln->P);
	    seq = XSeqPtr(E); s = SegmentStart(S[i]);
	    e = s+length-1; end = LenSeq(E);
	    if(s > 0 && e <= end){ Add2Profile(seq,s,Aln->M); }
	}
	MEW(Aln->maybe,Aln->nseg+2,s_type);
        for(j=0,I=1;I <= NSeqsSeqSet(Aln->P);I++) {
                E=SeqSetE(I,Aln->P);
		end = (Int4)LenSeq(E) - length + 1;
                for(s=1; s <= end; s++) {
                    Aln->maybe[j]= Segment(I,s);
                    j++;
                }
        }
	Aln->maybe[j]=NULL;
	return Aln;
}

aln_typ	*MkAlignsCombine(Int4 *numGrps, ptn_typ *motif, s_type **list, 
	double hg_cut, ss_type P, double rcut, Int4 N)
{
     Int4	os,group,g,i,j,k,k1,k2,numG;
     Int4	N1,N2,n,x,sC,maxnseg;
     s_type	*LC,*L2,*L1,*tL;
     ptn_typ	Q,Q1,Q2,sQ;
     aln_typ	*align;
     BooLean	verbose=FALSE;

     numG=*numGrps;
     if(numG < 1) align_error("CombineAlign( ) number = 0");
     NEW(align,numG+2,aln_typ);
     for(maxnseg=0,g=0; g<numG; g++) {
	k = LengthPattern(motif[g]); 
	align[g] = CreateAlign(P, k, list[g]); 
	if(verbose) PutAlign(stdout, align[g], TRUE, motif[g]);
	sC = RefineAlign(align[g], 1.0);
	maxnseg += sC;
     }
    if(numG > 1){
     NEW(L2,maxnseg+2,s_type); NEW(L1,maxnseg+2,s_type);
     NEW(LC,maxnseg+2,s_type);
     for(group=numG,j=0; group > 1 && j < numG; j++){
	 if((Q1=motif[j]) !=NULL){
           fprintf(stderr,".");
	   k1 = LengthPattern(Q1); 
	   N1 = CopySegments(L1,SegsAlign(align[j]));
	   BubbleSortSegments(L1);
	   for(i=0; i < numG; i++){	  
		if((Q2=motif[i]) != NULL && Q2 != Q1){
	   	   k2 = LengthPattern(Q2);
	  	   n = CopySegments(L2,SegsAlign(align[i]));
	   	   BubbleSortSegments(L2);
                   if(IntersectOSegments(&x,7,&os,k1,k2,L1,L2) == 1 && x>1){
			if(os > 0) {
				OffsetSegments(-os,L2,LC);
				tL = LC; LC = L2; L2 = tL;
			} else {
				OffsetSegments(os,L1,LC);
				tL = LC; LC = L1; L1 = tL;
			} N2 = N - N1;
			if(CumHypGeomProb(N1,N2,n,x) < hg_cut 
				|| x==N1 || x==n){
			   if(os > 0) { sQ = MergePatterns(Q2,Q1,-os); }
			   else { sQ = MergePatterns(Q1,Q2,os); }
			   sC = UnionSegments(LC,L2,L1);
			   if(verbose){
				fprintf(stdout,"\t*** %d %d\n",i,j);
				PutAlign(stdout,align[j],TRUE,Q1);
				fprintf(stdout,"\tplus...\n");
				PutAlign(stdout,align[i],TRUE,Q2);
				fprintf(stdout,"\t= merged alignment...\n");
				fprintf(stdout,"\tos = %d;\n", os);
			   }
			   NilPattern(Q1); NilPattern(Q2);
			   NilAlign(align[i]); NilAlign(align[j]);
			   motif[j] = motif[i] = NULL; 
			   align[j] = align[i] = NULL;
			   if(j > i) j = i; 
			   Q=motif[j]=Q1=sQ; 
			   k1 = LengthPattern(sQ);
			   NEW(tL,sC+3,s_type); CopySegments(tL,LC);
			   align[j]=CreateAlign(P, k1, tL); 
			   if(verbose) PutAlign(stdout,align[j],TRUE,Q);
			   sC = RefineAlign(align[j], rcut);
			   if(verbose) PutAlign(stdout,align[j],TRUE,Q);
			   if(sC > (N1 + n)){
			       maxnseg -= (N1 + n); maxnseg += sC;
			       NilSegmentList(L2); NEW(L2,maxnseg+2,s_type); 
			       NilSegmentList(L1); NEW(L1,maxnseg+2,s_type);
			       NilSegmentList(LC); NEW(LC,maxnseg+2,s_type);
			   }
	  	   	   N1 = CopySegments(L1,SegsAlign(align[j]));
	   		   BubbleSortSegments(L1);
			   group--; i=0;  			/* restart */
			   if(group < 2){ i=numG; j=-1; break; } /** quit **/
			} else {
			   CopySegments(L1,SegsAlign(align[j]));
	   		   BubbleSortSegments(L1);
			}
                   } 
		}
	   }
	 }
     }
     NilSegmentList(L2); NilSegmentList(L1); NilSegmentList(LC);
    }
     for(k=i=0; i < numG; i++){
	   if(motif[i] != NULL){
		motif[k] = motif[i]; align[k] = align[i]; k++;
	   }
     }
     motif[k] = NULL; align[k] = NULL;
     *numGrps=k;
     return align;
}

void	NilAlign(aln_typ Aln)
{ 
	NilMheap(Aln->H);
	if(Aln->M != NULL) NilProfile(Aln->M); 
	if(Aln->scores != NULL) free(Aln->scores); 
	NilSegmentList(Aln->S); 
	NilSegmentList(Aln->maybe); free(Aln); 
}

Int4	RefineAlign(aln_typ Aln, double cut)
/****************************************************************
 Adds related segments and removes unrelated segments from alignment.
 returns cardinality of new segment list.
 ****************************************************************/
{
	Int4	n,i,C;
	s_type	*S2,*S3,*S4,*S=Aln->S;
	BooLean	redo;

	if(Aln == NULL || S[0]==NULL)  return 0;
	if(Aln->refined) align_error("already refined");
	else Aln->refined = TRUE;
	C = Aln->n;
	if(C > 3){
	   NEW(S3,Aln->N+2,s_type); 
	   NEW(S4,Aln->N+2,s_type);
	   i=0; do{ S3[i]=S[i]; } while(S[i++] != NULL);
	   C=AddProfileAlign(S4,1.0, Aln,FALSE);
	   n=RmProfileAlign(S3,S4,1.0, Aln); redo = TRUE;
	   if(n>0){
	      C=n; n=AddProfileAlign(S4,1.0,Aln,TRUE);
	      if(C==n){
		for(redo=FALSE, i=0; S4[i] != NULL && S3[i]!=NULL; i++){
			if(S4[i] != S3[i]) { redo = TRUE; break; }
		}
		if(S4[i] != S3[i]) redo = TRUE; 
	      } else if(n==0) redo = FALSE;
	      if(redo){
		S2 = S3; S3 = S4; S4 = S2;
		n=AddProfileAlign(S4,1.0,Aln,TRUE);
	      } 
	      if(n>0) {
		n = RmProfileAlign(S3,S4,cut, Aln);
#if 0
	        /*** ProbProfileAlign(S3, Aln, cut); /***/
/** TEST ** Aln->S=S3; PutAlign(stdout,Aln,TRUE,NULL);/****/
#endif
	      } else return 0;
	   }
	   NilSegmentList(S);  NEW(S,n+3,s_type);
	   i=0; do{ S[i]=S3[i]; } while(S3[i++] != NULL);
	   NilSegmentList(S3);  NilSegmentList(S4);  Aln->S=S; 
	   return n;
	} else return C; 
}

Int4	AddProfileAlign(s_type *fS,double cut,aln_typ Aln,BooLean fit)
/* go through segment population and select all high scoring segments */
{
	Int4	I,I2,i,k,sdex,itm,s,e,end,score,len[3];
	Int4	N = Aln->N;
	s_type  *S2,S;
	mst_type	St;
	e_type	E;
	unsigned char	*seq;
	double	cutoff;
	mh_type	H;

	H = Aln->H;
	k = Aln->leng; 
	cutoff = cut/(double)Aln->nseg;
	NEW(S2,N+2,s_type); 
        for(I2=0,sdex=0; (S=Aln->maybe[sdex]) != NULL; sdex++) {
	    I=SegmentI(S); 
	    if(I2 != I){ I2 = I; E=SeqSetE(I,Aln->P); seq = XSeqPtr(E); }
	    s = SegmentStart(S);
	    if(s > 0 && (s+k-1) <= (Int4) LenSeq(E)){
		score = ScoreProfile(seq, s, Aln->M);
		if(ProfileProb(score,Aln->M) < cutoff){
		   if((itm=InsertMheap(-(keytyp)score,H)) != NULL){
			CpSegment2List(itm,S2,S);
		   } 
		} 
	    }
        }
	if(fit){ len[1]=k; St=MkMSites(1,len,Aln->P); }
	NilProfile(Aln->M); 
	Aln->M = MkProfile(Aln->nseg,k,tFreqSeqSet(Aln->P),Aln->A);
	for(i=0; (itm = DelMinMheap(H))!= 0; ){
	    if(fit){
		fS[i]=NULL;
		I = SegmentI(S2[itm]); s = SegmentStart(S2[itm]);
		if(!OccupiedMSite(1,I,s,St)){
			AddMSite(1,I,s,St); fS[i] = S2[itm]; 
		}
	    } else fS[i] = S2[itm];
	    if(fS[i]!=NULL){
		I=SegmentI(fS[i]); E = SeqSetE(I,Aln->P);
		seq = XSeqPtr(E); s = SegmentStart(fS[i]);
		e = s+k-1; end = LenSeq(E);
		if(s > 0 && e <= end){ Add2Profile(seq,s,Aln->M); }
		i++;
	    }
	}
	fS[i] = NULL; NilSegmentList(S2); 
	if(fit) NilMSites(St);
	return i; 
}
  
Int4	RmProfileAlign(s_type *fS,s_type *mS, double cut, aln_typ Aln)
/* go through segment population and select all high scoring segments 
   WARNING: assumes segments are ordered by decreasing profile score! */
{
	Int4	I,i,j,k,itm,s,score,min_s;
	s_type  *S2,S;
	e_type	E;
	unsigned char	*seq,*minseq;
	double	cutoff;
	mh_type	H;
	Int4	N = Aln->N;

	H = Aln->H; k = Aln->leng; i = Aln->nseg;
	cutoff = cut/(double)i; NEW(S2,N+2,s_type); 
	for(i=0; mS[i]!=NULL; i++){
	    fS[i] = mS[i]; 
 	    I=SegmentI(mS[i]); E = SeqSetE(I,Aln->P);
	    seq = XSeqPtr(E); s = SegmentStart(mS[i]);
	    if(s <= 0 || (s+k-1) > (Int4) LenSeq(E))
		print_error("Segment error in RmProfileAlign( )");
	    if(mS[i+1] == NULL) { 
		minseq = seq; min_s = s; RmProfile(seq,s,Aln->M); 
	    }
	}
	fS[i] = NULL;
	while(i>0){
		score = ScoreProfile(minseq, min_s, Aln->M);
		if(ProfileProb(score,Aln->M) > cutoff) { 
		  fS[i-1] = NULL;	/* wipe out min segment */
		  for(j=0; fS[j] != NULL; j++){
	             I = SegmentI(fS[j]); E=SeqSetE(I,Aln->P); 
		     seq = XSeqPtr(E); s = SegmentStart(fS[j]);
		     score = ScoreProfile(seq, s, Aln->M);
		     if((itm=InsertMheap(-(keytyp)score,H)) != 0){
			CpSegment2List(itm,S2,fS[j]);
		     } else print_error("this should not happen!");
		  }
		  for(S=NULL,i=0; (itm = DelMinMheap(H))!=0; i++){
			S = fS[i] = S2[itm]; 
		  }
		  if(S==NULL) break; /** no segments left **/
	          I = SegmentI(S); E=SeqSetE(I,Aln->P); 
		  minseq = XSeqPtr(E); min_s = SegmentStart(S);
		  RmProfile(minseq,min_s,Aln->M); 
		} else { Add2Profile(minseq,min_s,Aln->M); break; /*done*/}
		fS[i] = NULL;
        } 
	fS[i] = NULL; NilSegmentList(S2); 
	return i; 
}

Int4	ProbProfileAlign(s_type *mS, aln_typ Aln, double cut)
/* Get profile probabilities for segments on list mS */
{
	Int4	itm,n,N,I,i,k,s,score,*scores;
	s_type	*fS;
	e_type	E;
	unsigned char	*seq;
	double	prob;
	mh_type	H;

	k = Aln->leng; 
	for(n=0; mS[n]!=NULL; ) n++;
	NEW(Aln->scores,n+2,Int4); 
	NEW(scores,n+2,Int4); 
	H = Mheap(n,3); NEW(fS,n+2,s_type);
	N = Aln->nseg;;
	for(i=0; mS[i]!=NULL; i++){
 	    I=SegmentI(mS[i]); E = SeqSetE(I,Aln->P);
	    seq = XSeqPtr(E); s = SegmentStart(mS[i]);
#if 0
RmProfile(seq,s,Aln->M); 
#endif
	    score = ScoreProfile(seq, s, Aln->M);
	    prob = (double)N * ProfileProb(score,Aln->M);
	    score = (Int4) (-100.0 * log10(prob) - 0.5); 
	    if(prob < cut && (itm=InsertMheap(-(keytyp)score,H)) != 0){
#if 0
Add2Profile(seq,s,Aln->M);
#endif
			CpSegment2List(itm,fS,mS[i]);
			scores[itm] = score;
	    }
#if 0
	else {
		n--; mS[i] = mS[n]; mS[n] = NULL;
		while(DelMinMheap(H)!= 0) ;
		if(i==0) break;
		else i=0; 
	    } else print_error("this should not happen!"); 
#endif
	}
	for(i=0; (itm = DelMinMheap(H))!= 0; i++){
	    mS[i] = fS[itm]; 
	    Aln->scores[i] = scores[itm];
        }
	mS[i] = NULL; 
	NilSegmentList(fS); NilMheap(H); free(scores);
	return i; 
}

BooLean	XSegCharAlign(Int4 j, char *r, Int4 i, aln_typ Aln)
{
        Int4 x,I;
	unsigned char c;
        e_type E;
	s_type *S=Aln->S;

        x = j+SegmentStart(S[i]);
        I=SegmentI(S[i]); E=SeqSetE(I,Aln->P);
        if(x >= 1 && x <= (Int4) LenSeq(E)){
		c = ResSeq(x,E); *r = AlphaChar(c,Aln->A);
                return (c != XSeq(x,E));
        } else { *r = ' '; return FALSE; }
}

Int4     NumSeqSegLAlign(aln_typ Aln)
/* returns the number of sequence entities for the segments on list S */
{
	Int4     I,i,n,N;
	BooLean	*hit;
	s_type *S=Aln->S;

	N = NSeqsSeqSet(Aln->P);
	MEW(hit,N+2,BooLean);
	for(I=1;I<=N;I++) hit[I]=FALSE;
	for(i=n=0; S[i] != NULL; i++) {
		I=SegmentI(S[i]);
		if(!hit[I]) { n++; hit[I] = TRUE; } 
	}
	free(hit);
	return n;
}

Int4	ScanfileAlign2(aln_typ Aln, BooLean verbose, FILE *snfptr)
/* create a scan file from segments in S */
{
	Int4	c,i,j,I,m,n,nseq,N,numsegs,length,s,times=100,first,laln;
	ss_type	P=Aln->P;
	a_type	A=Aln->A; 
	s_type *S=Aln->S;
	e_type	E;
	BooLean	*rich;
	unsigned char	*str,d;
	char	*tstr;
	double	**seg_freq,*site_freq,info,maxinfo,summax,*tfreq;

	length = Aln->leng;
	tfreq = tFreqSeqSet(Aln->P);
	nseq = NumSeqSegLAlign(Aln); N = NSeqsSeqSet(Aln->P);
	if(snfptr == NULL) return 0;
	for(n=numsegs=0; S[n]!=NULL; n++) numsegs++;
	NEW(rich,length,BooLean);
	NEW(site_freq,nAlpha(A)+2,double); NEWP(seg_freq,length+2,double);
	for(i=0; i<length; i++) { NEW(seg_freq[i],nAlpha(A)+2,double); }

	fprintf(stderr,"\n");
	for(summax=0.0, n=0; n<times; n++) {
	  maxinfo = -999;
	  for(j=0; j<=10; j++) {
	   for(i=0; S[i]!= NULL; i++) {
		I=SegmentI(S[i]); E=SeqSetE(I,Aln->P);
		d = RandomResSeq(E); 
#if 0
fprintf(stderr,"%c",AlphaChar(d,A));
#endif
		site_freq[d] += 1.0;
	   }
	   for(info=0.0,c=0; c <= nAlpha(A); c++) { 
		site_freq[c] /= numsegs; 
		if(site_freq[c] > 0.0){
	   	    info += site_freq[c] * log(site_freq[c]/tfreq[c]);
		}
	   	site_freq[c] = 0.0;
	   }
#if 0
fprintf(stderr,"info = %g\n",info);
#endif
	   maxinfo = MAXIMUM(double,info,maxinfo);
          }
	  summax += maxinfo;
	}
	maxinfo = summax/times;
	for(i=m=0; S[i]!=NULL; i++){
	   I=SegmentI(S[i]); E=SeqSetE(I,Aln->P); s=SegmentStart(S[i]); 
	   str=SeqPtr(E);
	   for(j=0; j<length; j++) seg_freq[j][str[s++]]+=1.0; 
	}
	for(first=-1,j=0; j < length; j++) {
	      for(info=0.0,c=0; c <= nAlpha(A); c++) { 
		seg_freq[j][c] /= numsegs; 
		if(seg_freq[j][c] > 0.0){
	   	    info += seg_freq[j][c] * log(seg_freq[j][c]/tfreq[c]);
		}
	      }
	      if(verbose){
		fprintf(stderr,"maxinfo = %g\n",maxinfo);
		fprintf(stderr,"column %d: info = %g ",j+1,info);
		if(info >= maxinfo) fprintf(stderr,"*\n");
		else fprintf(stderr,".\n");
	      }
	      if(info >= maxinfo) {
		if(first<0) first = j;
		rich[j] = TRUE; laln = j;
	      } else rich[j] = FALSE;
	}
    if((laln-first+1) > 5){ 
	if(verbose) fprintf(snfptr,"AL   ");
	for(j=first; j <= laln; j++) {
		if(rich[j]) fprintf(snfptr,"*");
		else fprintf(snfptr,".");
	}
	fprintf(snfptr,"\n");
	for(i=0; S[i]!=NULL; i++){
	   if(verbose) fprintf(snfptr,"     ");
	   I=SegmentI(S[i]); E=SeqSetE(I,Aln->P); s=SegmentStart(S[i]); 
	   str=SeqPtr(E);
	   for(j=0; j < first; j++) s++;
	   for(; j <= laln; j++) {
		fprintf(snfptr,"%c", AlphaChar(str[s++],A));
	   }
	   if(verbose) {
		fprintf(snfptr," (%d ",SegmentStart(S[i])); 
		tstr = SeqKey(E);
		for(j=0;!isspace(tstr[j]) && tstr[j] != 0; j++) 
			fprintf(snfptr,"%c",tstr[j]);
	        fprintf(snfptr,")\n");
	   } else fprintf(snfptr,"\n");
	}
	fprintf(snfptr,"\n");
    }
	for(j=0; j < length; j++) free(seg_freq[j]);
	free(seg_freq); free(site_freq); free(rich);
	return numsegs;
}

Int4	PutAlign(FILE *fptr, aln_typ Aln, BooLean verbose, ptn_typ Q)
{
	if(Q != NULL) return put_align(fptr, Aln, verbose, Q);
	else return put_align2(fptr, Aln, verbose);
}

Int4	put_align2(FILE *fptr, aln_typ Aln, BooLean verbose)
{
	Int4	i,j,I,nseq,length,start;
	char	*tmp;
	s_type *S= Aln->S;
	e_type	E;

	length = Aln->leng;
	nseq = NumSeqSegLAlign(Aln);
	for(i=0; S[i]!=NULL; i++){
		I=SegmentI(S[i]); E=SeqSetE(I,Aln->P);
		start = SegmentStart(S[i]);
		fprintf(fptr,"%-4d",I);
		PutSeqRegion(fptr,start, length, E, Aln->A);
		if(verbose) { 
			fprintf(fptr,"  "); tmp = SeqKey(E);
			for(j=0;!isspace(tmp[j]) && tmp[j] != 0; j++) 
				fprintf(fptr,"%c",tmp[j]);
			fprintf(fptr,"\n");
		}
		else fprintf(fptr,"\n");
		
	}
	fprintf(fptr,"    (%d segments in %d sequences)\n\n", i, nseq);
	return nseq;
}

Int4	FinalProbAlign(aln_typ Aln)
/* Get profile probabilities for segments on list mS */
{
	Int4	itm,n,N,I,i,k,s,score,*scores;
	s_type	*fS,*mS=Aln->S;
	e_type	E;
	unsigned char	*seq;
	double	prob;
	mh_type	H;

	k = Aln->leng; 
	for(n=0; mS[n]!=NULL; ) n++;
	NEW(Aln->scores,n+2,Int4); 
	NEW(scores,n+2,Int4); 
	H = Mheap(n,3); NEW(fS,n+2,s_type);
	N = Aln->nseg;;
	for(i=0; mS[i]!=NULL; i++){
 	    I=SegmentI(mS[i]); E = SeqSetE(I,Aln->P);
	    seq = XSeqPtr(E); s = SegmentStart(mS[i]);
	    RmProfile(seq,s,Aln->M);
	    score = ScoreProfile(seq, s, Aln->M);
	    prob = (double)N * ProfileProb(score,Aln->M);
	    score = (Int4) (-100.0 * log10(prob) - 0.5); 
	    if(score > 0 && (itm=InsertMheap(-(keytyp)score,H)) != NULL){ 
			Add2Profile(seq,s,Aln->M);
			CpSegment2List(itm,fS,mS[i]);
			scores[itm] = score;
	    }
	}
	for(i=0; (itm = DelMinMheap(H))!=NULL; i++){
	    mS[i] = fS[itm]; 
	    Aln->scores[i] = scores[itm];
        }
	mS[i] = NULL; 
	NilSegmentList(fS); NilMheap(H); free(scores);
	return i; 
}

Int4	put_align(FILE *fptr, aln_typ Aln, BooLean verbose, ptn_typ Q)
{
	Int4	c,i,j,k,end,I,I2,n,nseq,length,start,*scores;
	BooLean	up;
	a_type	A; 
	s_type *S;
	e_type	E;
	char	*tstr,*str,str2[10],L2[33],h,d;

	FinalProbAlign(Aln);
	S= Aln->S;
	scores=Aln->scores;
	A=Aln->A;  
	length = Aln->leng;
	if(length != LengthPattern(Q)) align_error("put_align( ) length error");
	fprintf(fptr,"\n%-*s   %-3s %5s %7s\n", LengFormatPattern(Q),
		"SEG","SEQ","SCORE", "LOCATION");
	nseq = NumSeqSegLAlign(Aln);
	for(n=0; S[n]!=NULL; n++) ;
	NEW(str,(2*length+1),char);
	if(S[0]!=NULL) { I2= SegmentI(S[0]); sprintf(str2,"%-4d",I2); }
	for(I2=k=i=0; S[i]!=NULL; i++){
		I=SegmentI(S[i]); E=SeqSetE(I,Aln->P);
		for(j=k=0,h=' ',up=FALSE;j<length;j++) {
	   	   if(NotEqUPattern(j,Q) && NotEmptyPattern(j,Q)) {
		   	if(!up) { str[k++]=' '; up=TRUE;}
			if(XSegCharAlign(j,&d,i,Aln)) h = 'L';
			c = AlphaCode(d,A);
			if(!MemPattern(j,c,Q)) {
				if(isalpha(d)) d=tolower(d);
				n = ResPatterns(j,L2,Q,A);
			} 
			str[k++]=d;
		   } else {
		   	if(up) { str[k++]=' '; up=FALSE; }
			if(XSegCharAlign(j,&d,i,Aln)) h='L';
			if(isalpha(d)) d = tolower(d);
			str[k++]=d;
		   }
		}
		str[k]=0; 
		if(I2 != I) { sprintf(str2,"%d",I); I2 = I; }
		else sprintf(str2,"    ");
		start = MAXIMUM(Int4,1,SegmentStart(S[i]));
		end = MINIMUM(Int4,(SegmentStart(S[i])+length-1),LenSeq(E));
		if(scores != NULL) fprintf(fptr,
			"%s %4s  %-3d%c (%d-%d)", str,
			str2, scores[i], h, start, end);
		else fprintf(fptr,"%s %4s   - %c (%d-%d)",
				str, str2, h, start, end);
		if(verbose){
			fprintf(fptr,"  "); tstr = SeqKey(E);
			for(j=0;!isspace(tstr[j]) && tstr[j] != 0; j++) 
				fprintf(fptr,"%c",tstr[j]);
		}
		fprintf(fptr,"\n");
	}
	for(j=0;j<k;j++) fprintf(fptr,"_");
	fprintf(fptr,"    (%d segments in %d sequences)\n", i, nseq);
	if(verbose) PutPatternAlign(fptr, Q, Aln);
	else fprintf(fptr,"\n\n");
	free(str);
	return nseq;
}

void	PutPatternAlign(FILE *fptr, ptn_typ Q, aln_typ Aln)
{
	Int4	tmax,c,i,length;
	BooLean	up;
	char	**motif;

	length = LengthPattern(Q);
	NEWP(motif,nAlpha(Aln->A)+2,char);
	for(c=0;c<=nAlpha(Aln->A);c++) NEW(motif[c],length+2,char);
	tmax = GetPattern(motif,Q,Aln->A); 
	for(c=0;c<tmax;c++) {
	   if(c>0) fprintf(fptr,"\n"); 
	   for(i=0,up=FALSE;i<length;i++) {
	   	if(NotEqUPattern(i,Q) && NotEmptyPattern(i,Q)) {
		   if(!up) { fprintf(fptr," "); up=TRUE; }
		} else if(up) { fprintf(fptr," "); up=FALSE; }
		fprintf(fptr,"%c",motif[c][i]);
	   }
	}
	fprintf(fptr,"\n"); 
	for(c=0;c<=nAlpha(Aln->A);c++) free(motif[c]);
	free(motif);
	fprintf(fptr,"\n"); 
}

BooLean	SeqInAlign(short I, aln_typ Aln)
/** return TRUE if a segment from sequence I is in the alignment
    else return FALSE **/
{
	Int4	i;
	s_type *S=Aln->S;

	for(i=0; S[i]!=NULL; i++) if(SegmentI(S[i]) == I) return TRUE;
	return FALSE;
}

Int4	ScanfileAlign(aln_typ Aln, BooLean verbose, FILE *snfptr)
/* create a scan file from segments in S */
{
	Int4	i,j,I,n,numsegs,length,s;
	a_type	A=Aln->A; 
	s_type *S=Aln->S;
	e_type	E;
	char	*tstr;
	unsigned char *str;

	length = Aln->leng;
	if(snfptr == NULL) return 0;
	for(n=numsegs=0; S[n]!=NULL; n++) numsegs++;
    if(length >= 5){ 
	if(verbose) fprintf(snfptr,"AL   ");
	for(j=0; j < length; j++) fprintf(snfptr,"*");
	fprintf(snfptr,"\n");
	for(i=0; S[i]!=NULL; i++){
	   if(verbose) fprintf(snfptr,"     ");
	   I=SegmentI(S[i]); E=SeqSetE(I,Aln->P); s=SegmentStart(S[i]); 
	   str=SeqPtr(E);
	   for(j=0; j < length; j++) {
		fprintf(snfptr,"%c", AlphaChar(str[s++],A));
	   }
	   if(verbose) {
		fprintf(snfptr," (%d ",SegmentStart(S[i])); 
		tstr = SeqKey(E);
		for(j=0;!isspace(tstr[j]) && tstr[j] != 0; j++) 
			fprintf(snfptr,"%c",tstr[j]);
	        fprintf(snfptr,")\n");
	   } else fprintf(snfptr,"\n");
	}
	fprintf(snfptr,"\n");
    }
	return numsegs;
}

void	align_error(const char *s) { fprintf(stderr,"Asset: %s\n",s); exit(1); }

