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

#include "pairaln.h"

/****************************************************************
 first:  (seq2=1)
    from s1 = 1
               v
               DSIMVNPPY             s1 = s2 = 1;
               IATNPPYVGHIKG    	   s1++;
               ^
   to s1 = n1-3;
               v
         DSIMVNPPY                       s1 = n1 - 3 = 9 - 3 = 6
               IATNPPYVGHIKG             (s1...s1+end && s2...end) 
               ^				s2 = 1;
 second:  
              v                                   v
              IATNPPYVGHIKG    		IATNPPYVGHIKG
               DSIMVNPPY        to                DSIMVNPPY             
               ^                                  ^			
******************************************************************/
BooLean	relate_seq(register unsigned char *seq1, register unsigned char *seq2,
	register char **R, register Int4 n, register Int4 cutoff)
{
	register Int4 min=0,sum=0;

	while(n > 0){
        	if(min < (sum += R[seq1[n]][seq2[n]])){
                       if(sum-min >= cutoff) return TRUE;
        	} else min = sum;
		n--;
        };
	return FALSE;
}

BooLean	RelatedSeqsDiag(Int4 cutoff,unsigned char *seq1,unsigned char *seq2,
	register Int4 len, a_type A)
{ return relate_seq(seq1,seq2,AlphaR(A),len,cutoff); }


BooLean	RelatedSeqs(Int4 cutoff, e_type E1, e_type E2, a_type A)
/*********************************************************************
Compare all diagonals of E1 and E2 to see if the two have an MSP
with score > cutoff.
 *********************************************************************/
{
	Int4	i,s,n1,n2,end,w;
	unsigned char	*seq1,*seq2;

	w = cutoff/highAlphaR(A);
	n1 = LenSeq(E1); seq1 = XSeqPtr(E1); 
	n2 = LenSeq(E2); seq2 = XSeqPtr(E2);
        for(end = n1 - w, s=0; s < end; s++,i--) {
	   	i = MINIMUM(Int4,n1-s,n2);
		if(relate_seq(seq1+s, seq2, AlphaR(A),i,cutoff)) return TRUE;
        }
        for(end = n2 - w, s=1; s <= end; s++) {
		i = MINIMUM(Int4,n2-s,n1);
		if(relate_seq(seq2+s, seq1, AlphaR(A),i,cutoff)) return TRUE;
	}
	return FALSE;
}

Int4	repeat_score_seq(register unsigned char *seq, register char **R, 
		register Int4 o, register Int4 n)
/********************************************************************
             /ptr=seq1+s
seq1: 	MQNKSQKETGDILGISQMHVSRL
seq2:	     MPPLFVMNNEILMHLRALKKTKKDVS
	     |...... n .......|
 ********************************************************************/
{
	register Int4 min=0,sum=0,score=-9999;

	for(seq++; n > o; seq++){
        	if(min > (sum += R[seq[0]][seq[o]])){
			min = sum;
		} else if(score < (sum-min)) score = (sum-min);
		n--;
        };
	return score;
}

Int4	PutRepeatsSeq(FILE *fptr,e_type E, a_type A, Int4 cutoff)
/** look for internal repeats **/
{
	Int4	i,s,n,r;
	unsigned char	*seq;

	n = LenSeq(E); seq = XSeqPtr(E); 
        for(r=0,i=1; i < n; i++) {
		s = repeat_score_seq(seq,AlphaR(A),i,n);
		if(s >= cutoff) { 
			PutDiagonalSeq(fptr, i, E, E, A);
			r++;
		}
        }
	return r;
}

Int4	RepeatScoreSeq(e_type E, a_type A, Int4 *offset)
/** look for internal repeats **/
{
	Int4	i,best,s,n,score;
	unsigned char	*seq;

	n = LenSeq(E); seq = XSeqPtr(E); 
        for(score= -9999,i=1; i < n; i++) {
		s = repeat_score_seq(seq,AlphaR(A),i,n);
		if(s > score) { best = i; score = s; }
        }
	*offset = best;
	return score;
}

Int4	diagonal_score_seq(register unsigned char *seq1, 
	register unsigned char *seq2, register char **R, register Int4 n)
/********************************************************************
             /ptr=seq1+s
seq1: 	MQNKSQKETGDILGISQMHVSRL
seq2:	     MPPLFVMNNEILMHLRALKKTKKDVS
	     |...... n .......|
 ********************************************************************/
{
	register Int4 min=0,sum=0,score=-9999;

	while(n > 0){
        	if(min > (sum += R[seq1[n]][seq2[n]])){
			min = sum;
		} else if(score < (sum-min)) score = (sum-min);
		n--;
        };
	return score;
}

Int4	AlignSeqFastp(e_type E1, e_type E2, a_type A)
/*******************************************************************
 Align two sequences using the fastp algorithm.
 see W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448 and 
     Wilbur WJ; Lipman DJ (1983) Rapid similarity searches of nucleic acid 
     and protein data banks.  Proc Natl Acad Sci U S A 80: 726-730.
 *******************************************************************/
{
	Int4	v,i,x,s,n,r,n1,n2,score=0,item,*off,*off0;
	Int4	*D,**Q,*nQ,nsave=10,*lastpos,*lastpos0,maxspace=10;
	unsigned char	*seq1,*seq2;
	mh_type	H;
	BooLean	hits=FALSE;
	/** may want to remove spacer stuff - doesn't seem to help **/
	keytyp	*key,*key0;  /** freq. of 1to5-spaces between matches **/

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2);
	/** Construct a lookup table for the query sequence (n1) using k=1 
	    words with positions of all words in seq1. **/
	NEW(nQ,nAlpha(A)+3,Int4);
	MEW(Q,nAlpha(A)+3,Int4*);
	for(i=0; i<= nAlpha(A); i++) NEW(Q[i],n1+3,Int4);
	for(s=1; s <= n1; s++) {
		r = seq1[s];
		Q[r][nQ[r]] = s; nQ[r]++;
	}
        /** For each residue in a database sequence look up positions for 
	    all identical residues in the table & calculate the difference 
	    in position for all matches.  (Store all matches with the same 
            offset in a second table). ***/
	v = n1 + n2; 
	NEW(off0,v+4,Int4); off = off0 + n2 + 2;
	NEW(lastpos0,v+4,Int4); lastpos = lastpos0 + n2 + 2;
	NEW(key0,v+4,keytyp); key = key0 + n2 + 2;
	for(s=1; s <= n2; s++) {
		r = seq2[s];
		for(i=0; i< nQ[r]; i++){
			x = Q[r][i] - s; /** find offset **/
			off[x]++;
			if(lastpos[x]!=0){
				n = s - lastpos[x]; /** get spacer **/
				if(n <= maxspace) key[x]+=4.0/(keytyp)n;
			} else lastpos[x] = s;
		}
	}
        /** locate regions of similarity using offsets with high # matches **/
	H = Mheap(nsave,3);
	NEW(D,nsave+1,Int4);
	for(i = 1-n2; i < n1; i++){
		if(i < 0) v = n2 + i; /** get length of diagonal **/
		else v = n1 - i;
		key[i] += (keytyp)off[i] - (keytyp)v/(keytyp)nAlpha(A); /**/
		if(key[i] > 3.0) {
		   item=InsertMheap(-key[i],H);
		   D[item] = i;
		   /**** printf("key[%d]=%g (v=%d; percent = %g)\n",
			i,key[i],v, 100.0*(double)key[i]/(double)v); /****/
		   /**** printf("off[%d]=%d (O-E=%f; v=%d; ave=%g)\n",
                        i,off[i],(keytyp)v/(keytyp)nAlpha(A),
			v, (double)off[i]/(double)v); /****/
		}
	}
        /** compute the score for the ten offsets of highest similarity **/
	/*** printf("\nn1 = %d; n2 = %d\n",n1,n2);/****/
	while((item=DelMinMheap(H))!=NULL){
		i = D[item];
		/*** fprintf(stderr,"key[%d] = %g; s[i] = ",i,key[i]); /****/
		if(i > 0){	/** -> start at seq1[i] **/
		   n = MINIMUM(Int4,n1-i,n2);
		   s = diagonal_score_seq(seq1+i, seq2, AlphaR(A),n);
		} else {	/** -> start at seq2[i] **/
		   n = MINIMUM(Int4,n2+i,n1);
		   s= diagonal_score_seq(seq2-i, seq1, AlphaR(A),n);
		}
		if(s > score) { score = s; v = i; }
		/*** fprintf(stderr,"%d; score = %d\n",s,score); /***/
		/*** PutDiagonalSeq(stdout, i, E1, E2, A); /** TEST **/
		hits=TRUE;
	}
	/** output high scoring segment **/
	if(hits) PutDiagonalSeq(stdout, v, E1, E2, A);
	else printf("no match\n");
	/** deallocate memory **/
	NilMheap(H); free(D);
	for(i=0; i<= nAlpha(A); i++) free(Q[i]);
	free(Q); free(nQ); 
	free(lastpos0); free(key0); free(off0);
	return score;
}

BooLean	RelateSeqFastp2(e_type E1, e_type E2, a_type A,Int4 score)
/** use this to test fastp alignment, etc. **/
{
	if(score <= AlignSeqFastp(E1, E2, A)) return TRUE;
	else return FALSE;
}

BooLean	RelateSeqFastp(e_type E1, e_type E2, a_type A,Int4 score)
/*******************************************************************
 See if two sequences are related with at or above a cutoff  score 
     using the fastp algorithm.
 *******************************************************************/
{
	Int4	i,x,*D,**Q,*nQ,r,s,n,n1,n2,item,nsave=10;
	unsigned char	*seq1,*seq2;
	mh_type	H;
	Int4	*off,*off0,*lastpos,*lastpos0; 
	keytyp	*key,*key0;
	Int4	maxspace=10;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2);
	NEW(nQ,nAlpha(A)+3,Int4);
	MEW(Q,nAlpha(A)+3,Int4*);
	for(i=0; i<= nAlpha(A); i++) MEW(Q[i],n1+3,Int4);
	for(i=1; i <= n1; i++) {
		r = seq1[i]; Q[r][nQ[r]] = i; nQ[r]++;
	}
	n = n1 + n2; 
	NEW(off0,n+4,Int4); off = off0 + n2 + 2;
	NEW(lastpos0,n+4,Int4); lastpos = lastpos0 + n2 + 2;
	NEW(key0,n+4,keytyp); key = key0 + n2 + 2;
	for(i=1; i <= n2; i++) {
		r = seq2[i];
		for(n=0; n < nQ[r]; n++){
			x = Q[r][n] - i; /** find offset **/
			off[x]++;
			if(lastpos[x]!=0){
				s = i - lastpos[x]; /** get spacer **/
				if(s <= maxspace) key[x]+=4.0/(keytyp)s;
			} else lastpos[x] = i;
		}
	}
       /** locate regions of similarity using offsets with high # matches **/
	H = Mheap(nsave,3);
	MEW(D,nsave+2,Int4);
	for(i = 1-n2; i < n1; i++){
		if(i < 0) n = n2 + i; /** get length of diagonal **/
		else n = n1 - i;
		key[i] += (keytyp)off[i] - (keytyp)n/(keytyp)nAlpha(A);
		if(key[i] > 3.0) {
			item=InsertMheap(-key[i],H);
			D[item] = i;
		}
	}
        /** compute the score for the ten offsets of highest similarity **/
	for(i=0; i<= nAlpha(A); i++) free(Q[i]);
	free(Q); free(nQ); free(lastpos0); free(key0); free(off0);
	while((item=DelMinMheap(H))!=NULL){
		i = D[item];
		if(i > 0){	/** -> start at seq1[i] **/
		   n = MINIMUM(Int4,n1-i,n2); 
		   if(relate_seq(seq1+i,seq2,AlphaR(A),n,score)){
			NilMheap(H); free(D); return TRUE;
		   }
		} else {	/** -> start at seq2[i] **/
		   n = MINIMUM(Int4,n2+i,n1); 
		   if(relate_seq(seq2-i,seq1,AlphaR(A),n,score)){
			NilMheap(H); free(D); return TRUE;
		   }
		}
	}
	NilMheap(H); free(D); return FALSE;
}

