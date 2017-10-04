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

#include "swaln.h"

Int4	SubSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A, 
	Int4 *start, Int4 *end)
{ Int4 startS; return SubSeqSW(a,b,E1,E2,A,start,end,&startS); }

Int4	SubSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A, 
	Int4 *start, Int4 *end,Int4 *startS)
/*******************************************************************
 Perform a Smith-Waterman alignment on sequences E1 and E2 and 
 return the score and set 'start' and 'end' equal to the start and
 end of the subsequence of E1 that aligns with E2.
 Modification: set startS = start of subject sequence.
********************************************************************/
{
	Int4	i,j,r1,s,t,n1,n2,score,max_i,max_j;
	unsigned char	*seq1,*seq2;
	Int4	p,w1,**D,**T,**E,**F;
	e_type	subE;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	/** seq1 = SeqPtr(E1); seq2 = SeqPtr(E2);	/*** no seg ***/
	/**/ seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2); /*** seg ***/
	/*** 1. Allocate memory. ***/
	MEW(D,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(E,n1+3,Int4*); MEW(F,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		NEW(D[i],n2+3,Int4); NEW(T[i],n2+3,Int4); 
		NEW(E[i],n2+3,Int4); NEW(F[i],n2+3,Int4); 
	}
	/*** 2. Dynamic programming step. ***/
	for(w1 = a+b, score = 0, i=1; i<= n1; i++) {
	   for(r1=seq1[i], j=1; j<= n2; j++) {
		s = D[i-1][j-1] + valAlphaR(r1,seq2[j],A); t=0;
                E[i][j] = MAXIMUM(Int4,D[i-1][j]-w1,E[i-1][j]-b);
                if(s < E[i][j]) { s = E[i][j]; t=1; }
                F[i][j] = MAXIMUM(Int4,D[i][j-1]-w1,F[i][j-1]-b);
                if(s < F[i][j]) { s = F[i][j]; t=-1; }
		T[i][j] = t;
		if(s > 0){
		   D[i][j] = s; 
		   if(s > score){ score = s; max_i = i; max_j = j; }
		}
	   }
	}
	/*** 3. Trace back step. ***/
	for(i=max_i,j=max_j; i > 0 && j > 0; ){
		if(D[i][j] == 0) break;
		switch(T[i][j]){
		  case 0: i--; j--; break;
		  case 1: i--; break;
		  case -1: j--; break;
		  default: print_error("this should not happen"); break;
		}
	} *start = i+1; *startS = j+1; *end = max_i;
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(D[i]);free(T[i]);free(E[i]);free(F[i]);}
	free(D); free(T); free(E); free(F); 
	return score;
}

static void    put_dp_alnsw(Int4 **D, Int4** T, e_type E1, e_type E2, a_type A,
        Int4 startA,Int4 endA, Int4 startB, Int4 endB)
{
        Int4    r,m,i,j,lenA=LenSeq(E1),lenB=LenSeq(E2);
	unsigned char *seqB=SeqPtr(E2),*seqA=SeqPtr(E1);

	endA = MINIMUM(Int4,endA,lenA); endB = MINIMUM(Int4,endB,lenB);
	startA = MAXIMUM(Int4,startA,1); startB = MAXIMUM(Int4,startB,1);
        fprintf(stderr,"     |");
        for(i=startA; i<=endA; i++)
		fprintf(stderr,"  %c  |", AlphaChar(seqA[i],A));
        fprintf(stderr,"\n     |");
        for(i=startA; i<=endA; i++) fprintf(stderr,"%3d  |",i);
        fprintf(stderr,"\n");
        for(r=startB; r<=endB; r++) {
           fprintf(stderr,"%c%4d|",AlphaChar(seqB[r],A),r);
           for(i=startA; i<=endA; i++) {
                if(D[i][r] == SHRT_MIN) fprintf(stderr," -inf ");
                else fprintf(stderr,"%5d",D[i][r]);
                if(T[i][r] == -1) fprintf(stderr,"^");
                else if(T[i][r] == 0) fprintf(stderr,"\\");
                else fprintf(stderr,"<");
           }
           fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
}

Int4	AlnSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A)
{ return AlnSeqSW(stdout,a,b,E1,E2,A); }

Int4	AlnSeqSW(FILE *fptr,Int4 a, Int4 b, e_type E1, e_type E2, a_type A)
/*******************************************************************
 Perform a Smith-Waterman alignment on sequences E1 and E2 
 returns optimum subalignment score.  W(k) = a + bk.
 D[0][j]=D[i][0]=0;
 		MAX{ D[i-1][j-1] + S(r1,r2),
 D[i][j] = 		D[i-k][j] - W(k) for k=1..i-1,
			D[i][j-k]- W(k) for k=1..j-1}.

 Find maximum score and trace back to get alignment.
 see T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197.
********************************************************************/
{
	Int4	i,j,r1,s,t,v,n1,n2,score,max_i,max_j;
	unsigned char	*seq1,*seq2;
	char	*out[3];
	Int4	*pos[3],p,w1,**D,**T,**E,**F;
	UInt4	os1,os2;

Int4    gINS,*gDEL,s0,s1;
	os1=OffSetSeq(E1); os2=OffSetSeq(E2);
	n1 = LenSeq(E1); n2 = LenSeq(E2);
	seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2);
	/*** 1. Allocate memory. ***/
	MEW(D,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(E,n1+3,Int4*); MEW(F,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		NEW(D[i],n2+3,Int4); NEW(T[i],n2+3,Int4); 
		NEW(E[i],n2+3,Int4); NEW(F[i],n2+3,Int4); 
	}
	/*** 2. Dynamic programming step. ***/
#if 0
gINS=0; NEW(gDEL,n2+3,Int4);
#endif
	for(w1 = a+b, score = 0, i=1; i<= n1; i++) {
	   for(r1=seq1[i], j=1; j<= n2; j++) {
#if 0
		s = D[i-1][j-1] + valAlphaR(r1,seq2[j],A); t=0;

		if((s0=D[i-1][j]-w1) > (s1=E[i-1][j]-b)){
			gDEL[j]=1;  E[i][j] = s0; 
		} else { gDEL[j]++;  E[i][j] = s1; }
                if(s < E[i][j]){ s=E[i][j]; t=gDEL[j]; } 

		if((s0=D[i][j-1]-w1) > (s1=F[i][j-1]-b)){
			gINS=-1;  F[i][j] = s0; 
		} else { gINS--; F[i][j] = s1; }
                if(s < F[i][j]){ s=F[i][j]; t=gINS; } 

		T[i][j] = t;
		if(s > 0){
			D[i][j] = s; 
			if(s > score){ score = s; max_i = i; max_j = j; }
		}
	   }
#endif
#if 1
		s = D[i-1][j-1] + valAlphaR(r1,seq2[j],A);
                /*********************************************************
                O(m*n) modified Gotoh method for w(k)=u*k+v where u,v >= 0
                 *********************************************************/
                t=0;
                E[i][j] = MAXIMUM(Int4,D[i-1][j]-w1,E[i-1][j]-b);
                if(s < E[i][j]) { s = E[i][j]; t=1; }
                F[i][j] = MAXIMUM(Int4,D[i][j-1]-w1,F[i][j-1]-b);
                if(s < F[i][j]) { s = F[i][j]; t=-1; }
		T[i][j] = t;
		if(s > 0){
			D[i][j] = s; 
			if(s > score){ score = s; max_i = i; max_j = j; }
		}
	   }
#endif
	}
// put_dp_alnsw(D, T, E1, E2, A, 130, 150, 74,127);
	/*** 3. Trace back step. ***/
	MEW(out[0],n1+n2+3,char); MEW(out[1],n1+n2+3,char);
	MEW(out[2],n1+n2+3,char);
	NEW(pos[1],n1+n2+3,Int4); NEW(pos[2],n1+n2+3,Int4);
#if 0
Int4 store,t0;
#endif
	Int4 total_aln=0,total_match=0;
	for(p=1,i=max_i,j=max_j; i > 0 && j > 0; ){
		if(D[i][j] == 0) break;
#if 0
		t0=T[i][j];	
// fprintf(stderr,"i=%d; j=%d; t0 =%d\n",i,j,t0);
	   do {
		if(t0 > 0){ t=1; t0--; }
		else if(t0 < 0){ t=-1; t0++; }
		else { t=0; }
		switch(t){
		  case 0:
			v = valAlphaR(seq1[i],seq2[j],A);
			pos[1][p] = i; pos[2][p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A);
			if(seq1[i]==seq2[j]){ total_match++; out[0][p] = ':'; }
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i--; j--;
			break;
		  case 1:
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			pos[1][p] = i; p++; i--;
			break;
		  case -1:
			out[2][p] = AlphaChar(seq2[j],A);
			out[0][p] = ' '; out[1][p] = '-';
			pos[2][p] = j; p++; j--;
			break;
		  default: print_error("this should not happen"); break;
		} total_aln++;
	   } while(t0 != 0);
#endif
#if 1
		switch(T[i][j]){
		  case 0:
			v = valAlphaR(seq1[i],seq2[j],A);
			pos[1][p] = i; pos[2][p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A);
			if(seq1[i]==seq2[j]){ total_match++; out[0][p] = ':'; }
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i--; j--;
			break;
		  case 1:
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			pos[1][p] = i; p++; i--;
			break;
		  case -1:
			out[2][p] = AlphaChar(seq2[j],A);
			out[0][p] = ' '; out[1][p] = '-';
			pos[2][p] = j; p++; j--;
			break;
		  default: print_error("this should not happen"); break;
		} total_aln++;
#endif
	}
	/** 4. Print out alignment **/
	fprintf(fptr,"\nscore = %d; %d/%d = %.1f%% identity\n\n",score,
		total_match,total_aln,
		100.0*(double)total_match/(double)total_aln);
	for(i=p-1; i > 0 ; ){
	   v = MINIMUM(Int4,50,i); fprintf(fptr,"  ");
	   for(j=i; j > i-v ; j--) {
		if(pos[1][j] && (pos[1][j]+os1) % 10 == 0) {
			fprintf(fptr,"%4d",pos[1][j]+os1); j-=3;
		} else fprintf(fptr," ");
	   }
	   fprintf(fptr,"\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fptr,"%c",out[1][j]); fprintf(fptr,"\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fptr,"%c",out[0][j]); fprintf(fptr,"\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fptr,"%c",out[2][j]); fprintf(fptr,"\n  ");
	   for(j=i; j > i-v ; j--) {
		if(pos[2][j] && (pos[2][j]+os2) % 10 == 0) {
			fprintf(fptr,"%4d",pos[2][j]+os2); j-=3;
		} else fprintf(fptr," ");
	   }
	   fprintf(fptr,"\n\n"); i-=v;
	}
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(D[i]);free(T[i]);free(E[i]);free(F[i]);}
	free(D); free(T); free(E); free(F); 
	free(out[0]); free(out[1]); free(out[2]);
	free(pos[1]); free(pos[2]);
// free(gDEL);
	return score;
}

Int4	FastAlnSeqSW(Int4 a, Int4 b, e_type E1, e_type E2, a_type A)
/*******************************************************************
 Perform a Smith-Waterman alignment on sequences E1 and E2 
 but compute score only.
********************************************************************/
{
	Int4	i,j,s,n1,n2,score,w1,**D,**E,**F;
	unsigned char	*seq1=XSeqPtr(E1),*seq2=XSeqPtr(E2);

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	MEW(D,n1+3,Int4*); MEW(E,n1+3,Int4*); MEW(F,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		NEW(D[i],n2+3,Int4); 
		NEW(E[i],n2+3,Int4); NEW(F[i],n2+3,Int4); 
	}
	for(w1 = a+b, score = 0, i=1; i<= n1; i++) {
	   for(j=1; j<= n2; j++) {
		s = D[i-1][j-1] + valAlphaR(seq1[i],seq2[j],A);
                E[i][j] = MAXIMUM(Int4,D[i-1][j]-w1,E[i-1][j]-b);
                if(s < E[i][j]) { s = E[i][j]; }
                F[i][j] = MAXIMUM(Int4,D[i][j-1]-w1,F[i][j-1]-b);
                if(s < F[i][j]) { s = F[i][j]; }
		if(s > 0) { 
			D[i][j] = s; 
			if(s > score){ score = s; }
		}
	   }
	}
	for(i=0; i<=n1; i++) {free(D[i]);free(E[i]);free(F[i]);}
	free(D); free(E); free(F); 
	return score;
}

float	RealAlnSeqSW(float a, float b, e_type E1, e_type E2)
/*******************************************************************
 Perform a Smith-Waterman alignment on sequences E1 and E2 
 but compute (a real) score only.
********************************************************************/
{
	Int4	i,j,n1,n2,w1;
	float	score, s,**D,**E,**F;
	unsigned char	*seq1=XSeqPtr(E1),*seq2=XSeqPtr(E2);

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	MEW(D,n1+3,float*); MEW(E,n1+3,float*); MEW(F,n1+3,float*);
	for(i=0; i<= n1; i++) { 
		NEW(D[i],n2+3,float); 
		NEW(E[i],n2+3,float); NEW(F[i],n2+3,float); 
	}
	for(w1 = a+b, score = 0, i=1; i<= n1; i++) {
	   for(j=1; j<= n2; j++) {
		s = D[i-1][j-1] + blosum62[seq1[i]][seq2[j]];
                E[i][j] = MAXIMUM(float,D[i-1][j]-w1,E[i-1][j]-b);
                if(s < E[i][j]) { s = E[i][j]; }
                F[i][j] = MAXIMUM(float,D[i][j-1]-w1,F[i][j-1]-b);
                if(s < F[i][j]) { s = F[i][j]; }
		if(s > 0) { 
			D[i][j] = s; 
			if(s > score){ score = s; }
		}
	   }
	}
	for(i=0; i<=n1; i++) {free(D[i]);free(E[i]);free(F[i]);}
	free(D); free(E); free(F); 
	return score;
}

