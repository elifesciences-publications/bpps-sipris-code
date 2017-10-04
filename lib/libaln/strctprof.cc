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

#include "strctprof.h"

void	WriteSProfile(char *outfile, spf_typ P)
{
        Int4    n,d;
        a_type  A = P->A;
	FILE	*fptr;

	fptr = open_file(outfile,"","w");
        for(n=1; n <= P->N; n++){
                fprintf(fptr,"%2d:",n);
                for(d=0; d<=nAlpha(A); d++){
                        fprintf(fptr," %3d", P->score[n][d]);
                }
                fprintf(fptr,"\n");
	}
	fclose(fptr);
}

spf_typ	ReadSProfile(char *infile, e_type E, a_type A)
{
        Int4    N,n,d,**score,i;
	FILE	*fptr;
	spf_typ	P;

	N = LenSeq(E);
	NEW(P,1,strct_pf_type);
	P->N = N = LenSeq(E);
	P->A = A;
	P->E = E;
	MEW(P->S,N+2,Int4*); MEW(P->I,N+2,Int4*); MEW(P->D,N+2,Int4*);
	MEW(P->rS,N+2,float*); MEW(P->rI,N+2,float*); MEW(P->rD,N+2,float*);
	NEW(P->is,N+2,Int4); NEW(P->ie,N+2,Int4);
	NEW(P->ds,N+2,Int4); NEW(P->de,N+2,Int4);
	NEWP(score,N+2,Int4);
	P->score = score;
	fptr = open_file(infile,"","r");
        for(n=1; n <= P->N; n++){
	   NEW(score[n],nAlpha(A)+2,Int4);
           if(fscanf(fptr,"%d:",&i) != 1) print_error("read error");
	   if(i!= n){
                fprintf(stderr,"n=%d; i=%d\n",n,i);
		print_error("*.prof input error");
	   }
	   for(d=0; d<=nAlpha(A); d++){
		if(fscanf(fptr," %d", &score[n][d]) != 1)
			print_error("read error");
	   }
	   P->is[n] = 120; P->ds[n] = 240;
	   P->ie[n] = 40; P->de[n] = 80;
	}
	fclose(fptr);
	return P;
}

spf_typ MkSProfile(double **S, e_type E, a_type A)
{
	spf_typ	P;
	Int4	r,n,N;
	double	s;
	h_type	H;

	NEW(P,1,strct_pf_type);
	P->N = N = LenSeq(E);
	P->A = A;
	P->E = E;
	MEW(P->S,N+2,Int4*); MEW(P->I,N+2,Int4*); MEW(P->D,N+2,Int4*);
	MEW(P->rS,N+2,float*); MEW(P->rI,N+2,float*); MEW(P->rD,N+2,float*);
	NEW(P->is,N+2,Int4); NEW(P->ie,N+2,Int4);
	NEW(P->ds,N+2,Int4); NEW(P->de,N+2,Int4);
	NEWP(P->score,N+2,Int4);
	H = Histogram("profile scores",-200,+200,5);
	for(n=1; n<=N; n++){
	   NEW(P->score[n],nAlpha(A)+2,Int4);
	   for(r=0; r<=nAlpha(A); r++){
		s = 10*S[r][n];
		P->score[n][r] = (Int4) floor(s+0.5);
		IncdHist(s,H);
		if(s < -200) printf("n=%d; r=%d; score = %g\n",n,r,s);
		/****/
	   }
	   P->is[n] = 120;
	   P->ds[n] = 240;
	   P->ie[n] = 40;
	   P->de[n] = 80;
	}
	PutHist(stdout,60,H);
        NilHist(H);
	return P;
}

void	Gaps2ndarySProfile(Int4 is[3], Int4 ie[3], Int4 ds[3], Int4 de[3],
	char *ss, spf_typ P)
/************************************************************************
 ************************************************************************/
{
	Int4	i,r,h=0,s=1,c=2;

        for(i=1; i<=P->N; i++){
           if(ss[i] =='h') {       /** helix **/
                   if(i<P->N && ss[i+1]=='c'){
                        SetGapsSProfile(i,is[c],ie[c],ds[h],de[h],P);
                   } else SetGapsSProfile(i,is[h],ie[h],ds[h],de[h],P);
           } else if(ss[i] =='s'){ /** strands **/
                   if(i<P->N && ss[i+1]=='c'){
                        SetGapsSProfile(i,is[c],ie[c],ds[s],de[s],P);
                   } else SetGapsSProfile(i,is[s],ie[s],ds[s],de[s],P);
           } else {                   /** coil **/
		   SetGapsSProfile(i,is[c],ie[c],ds[c],de[c],P);
	/*******
		   SetGapsSProfile(i,0,0,0,0,P);
		if((i<P->N && ss[i+1]!='c') || (i > 1 && ss[i-1] != 'c')){
                   SetGapsSProfile(i,is[h],ie[h],ds[h],de[h],P);
		} else {
		}
	/*******/
	/******/
		for(r=0; r<=nAlpha(P->A); r++){
			P->score[i][r] *= 0.1;
		}
	/******/
	/******
		for(r=0; r<=nAlpha(P->A); r++){
			P->score[i][r] = 0;
		}
	/******/
           }
	}
}

void	SetGapsSProfile(Int4 n, Int4 is, Int4 ie, Int4 ds, Int4 de, 
	spf_typ P)
{
	if(n <= P->N && n > 0){
	   P->is[n] = is; P->ie[n] = ie;
	   P->ds[n] = ds; P->de[n] = de;
	}
}

void	NilSProfile(spf_typ P)
{
	free(P->S); free(P->I); free(P->D);
	free(P->rS); free(P->rI); free(P->rD);
	free(P->is); free(P->ie); free(P->ds); free(P->de);
	free(P);
}

void	PutSProfile(FILE *fptr, spf_typ P)
{
        Int4    n,d;
        a_type  A = P->A;

        fprintf(fptr,"\n    ");
        for(d=0; d<=nAlpha(A); d++){
                fprintf(fptr,"  %c ", AlphaChar(d,A));
        }
        fprintf(fptr,"\n");
        for(n=1; n <= P->N; n++){
                fprintf(fptr,"%2d: ",n);
                for(d=0; d<=nAlpha(A); d++){
                        fprintf(fptr,"%3d ", P->score[n][d]);
                }
                fprintf(fptr,"\n");
	}
	fprintf(fptr,"\npos:  is:  ds:  ie:  de:\n");
        for(n=1; n <= P->N; n++){
		fprintf(fptr,"%3d: %4d %4d %4d %4d\n",
			n,P->is[n],P->ds[n], P->ie[n], P->de[n]);
			
        }
        fprintf(fptr,"\n");
}

Int4	PutAlnSProfile(e_type E, spf_typ P)
/*******************************************************************
 Perform a Smith-Waterman profile alignment on sequence E and 
 return optimum subalignment score.  

 S(0,m) = S(n,0) = 0;
 S(n,m) = MAXIMUM(S(n-1,m-1)+s(seq(m),n), I(n,m),D(n,m))

 where  D(n,m) = MAXIMUM(S(n-1,m) - del(n,1), D(n-1,m) - de(n)),

        I(n,m) = MAXIMUM(S(n,m-1) - ins(n,1), I(n,m-1) - ie(n)),

        del(n,1) = ds(n) + de(n), and isn(n,1) = is(n) + ie(n).

 Used O(m*n) modified method like that of Gotoh.
 Find maximum score and trace back to get alignment.
 see T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197.
********************************************************************/
{
	Int4	s,t,v,n,m,N,M,score,max_n,max_m,p;
	Int4	*pos[3],**S,**T,**D,**I,ins_n,del_n,de_n,ie_n;
	unsigned char	*seq0,*seq1,*seq2;
	char	*out[3];
	e_type	E0;
	a_type	A = P->A;

	M = LenSeq(E); N = P->N;
	S = P->S; I = P->I;  D= P->D;
	seq1 = XSeqPtr(E); seq2 = XSeqPtr(P->E);
	E0 = CopySeq(P->E); seq0=XSeqPtr(E0);
	MEW(T,N+3,Int4*);
	for(n=0;n<=N; n++){
		seq0[n]=0;
		NEW(S[n],M+3,Int4); NEW(T[n],M+3,Int4); 
		NEW(D[n],M+3,Int4); NEW(I[n],M+3,Int4); 
	}
	for(score=0, n=1; n<= N; n++) {
	   ins_n = P->is[n] + P->ie[n];
	   del_n = P->ds[n] + P->de[n];
	   de_n = P->de[n]; ie_n = P->ie[n];
	   for(m=1; m<= M; m++) {
		s = S[n-1][m-1] + P->score[n][seq1[m]];
		t=0;
		D[n][m] = MAXIMUM(Int4,S[n-1][m]-del_n,D[n-1][m]-de_n);
		I[n][m] = MAXIMUM(Int4,S[n][m-1]-ins_n,I[n][m-1]-ie_n);
		if(s < D[n][m]) { s = D[n][m]; t=-1; }
		if(s < I[n][m]) { s = I[n][m]; t=1; }
		if(s > 0){
			S[n][m] = s; 
			if(s > score){ score = s; max_n = n; max_m = m; }
		}
		T[n][m] = t;
	   }
	}
	printf("max = score(%d,%d) = %d\n",max_n,max_m,score);
	/*** Trace back step ***/
	MEW(out[0],N+M+3,char); MEW(out[1],N+M+3,char);
	MEW(out[2],N+M+3,char);
	NEW(pos[1],N+M+3,Int4); NEW(pos[2],N+M+3,Int4);
	for(p=1,s=score,n=max_n,m=max_m; n > 0 && m > 0; ){
		if(S[n][m] == 0) break;
		switch(T[n][m]){
		  case 0:
			v = P->score[n][seq1[m]];
			pos[1][p] = m; pos[2][p] = n;
			out[1][p] = AlphaChar(seq1[m],A);
			out[2][p] = AlphaChar(seq2[n],A);
			/*** seq0 for threading ***/
			seq0[n] = seq1[m];
			/*** seq0 for threading ***/

			if(seq1[m]==seq2[n]) out[0][p] = ':';
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; n--; m--;
			break;
		  case 1:
			out[1][p] = AlphaChar(seq1[m],A);
			out[0][p] = ' '; out[2][p] = '-';
			pos[1][p] = m; p++; m--;
			break;
		  case -1:
			out[2][p] = AlphaChar(seq2[n],A);
			out[0][p] = ' '; out[1][p] = '-';
			pos[2][p] = n; p++; n--;
			break;
		  default: print_error("this should not happen"); break;
		}
	}
	/** print out alignment **/
	printf("\n\n");
	for(n=p-1; n > 0 ; ){
	   v = MINIMUM(Int4,50,n); printf("  ");
	   for(m=n; m > n-v ; m--) {
		if(pos[1][m] && pos[1][m] % 10 == 0) {
			printf("%4d",pos[1][m]); m-=3;
		} else printf(" ");
	   }
	   printf("\n     ");
	   for(m=n; m > n-v ; m--) printf("%c",out[1][m]); printf("\n     ");
	   for(m=n; m > n-v ; m--) printf("%c",out[0][m]); printf("\n     ");
	   for(m=n; m > n-v ; m--) printf("%c",out[2][m]); printf("\n  ");
	   for(m=n; m > n-v ; m--) {
		if(pos[2][m] && pos[2][m] % 10 == 0) {
			printf("%4d",pos[2][m]); m-=3;
		} else printf(" ");
	   }
	   printf("\n\n"); n-=v;
	}
	PutSeq(stdout,E0,A);
	NilSeq(E0);
	/*** free allocated memory ***/
	for(n=0; n<=N; n++){ free(T[n]);free(I[n]);free(D[n]);free(S[n]); }
	free(T); free(pos[1]); free(pos[2]);
	free(out[0]); free(out[1]); free(out[2]);
	return score;
}

Int4	AlnSProfile(e_type E, spf_typ P)
/**********************************************************************
 Perform a Smith-Waterman profile alignment on sequence E without 
 traceback and return the optimum subalignment score.  

    DELETION:				INSERTION:

     M Y P R O F I L E                   P R O F - - I L E
    (h h h h c c h h h)                 (h h c c     h h h )
     S E Q U - - E N C                   S E Q U E N C E X 
           0 1 2 3                             0 1 2 3
                                                      
   0. S[n-1][m-1]+score[n][m]          0. S[n-1][m-1]+score[n][m]
   1. D[n][m]=S[n-1][m]-del[n]         1. I[n][m]=S[n][m-1]-ins[n]
   2. S[n][m]=D[n-1][m]-de[n]          2. S[n][m]=I[n][m-1]-ie[n]
   3. S[n-1][m-1]+score[n][m]          3. S[n-1][m-1]+score[n][m]

   M  Y  P  R  O  F  I  L  E            M  Y  P  R  O  F  I  L  E
 S[*][ ][ ][ ][ ][ ][ ][ ][ ]         S[ ][ ][*][ ][ ][ ][ ][ ][ ] m++
 E[ ][*][ ][ ][ ][ ][ ][ ][ ]         E[ ][ ][ ][*][ ][ ][ ][ ][ ] | 
 Q[ ][ ][*][ ][ ][ ][ ][ ][ ]         Q[ ][ ][ ][ ][*][ ][ ][ ][ ] |
 U[ ][ ][ ][0][1][2][ ][ ][ ]         U[ ][ ][ ][ ][ ][0][ ][ ][ ] V
 E[ ][ ][ ][ ][ ][ ][3][ ][ ]         E[ ][ ][ ][ ][ ][1][ ][ ][ ]  
 N[ ][ ][ ][ ][ ][ ][ ][.][ ]         N[ ][ ][ ][ ][ ][2][ ][ ][ ]
 C[ ][ ][ ][ ][ ][ ][ ][ ][.]         C[ ][ ][ ][ ][ ][ ][3][ ][ ]  
 E[ ][ ][ ][ ][ ][ ][ ][ ][ ]         E[ ][ ][ ][ ][ ][ ][ ][.][ ]
 X[ ][ ][ ][ ][ ][ ][ ][ ][ ]         X[ ][ ][ ][ ][ ][ ][ ][ ][.]

 **********************************************************************/
{
	Int4	s,n,m,N,M,score,**S,**D,**I,ins,del,de,ie;
	unsigned char	*seq;

	M = LenSeq(E); N = P->N;
	S = P->S; I = P->I;  D= P->D;
	seq = XSeqPtr(E); 
	for(n=0;n<=N; n++){
	   NEW(S[n],M+3,Int4); NEW(D[n],M+3,Int4); NEW(I[n],M+3,Int4); 
	}
	for(score=0, n=1; n<= N; n++) {
	   ins = P->is[n] + P->ie[n]; ie = P->ie[n];
	   del = P->ds[n] + P->de[n]; de = P->de[n];
	   for(m=1; m<= M; m++) {
		s = S[n-1][m-1] + P->score[n][seq[m]];
		D[n][m] = MAXIMUM(Int4,S[n-1][m]-del,D[n-1][m]-de);
		I[n][m] = MAXIMUM(Int4,S[n][m-1]-ins,I[n][m-1]-ie);
		if(s < D[n][m]) s = D[n][m]; 
		if(s < I[n][m]) s = I[n][m]; 
		if(s > 0){
			S[n][m] = s; 
			if(s > score) score = s; 
		}
	   }
	}
	for(n=0; n<=N; n++){ free(I[n]);free(D[n]);free(S[n]); }
	return score;
}

double	RealAlnSProfile(e_type E, spf_typ P)
{
	Int4	s,n,m,N,M,ins,del,de,ie;
	float	score,**S,**D,**I;
	unsigned char	*seq;

	M = LenSeq(E); N = P->N;
	S = P->rS; I = P->rI;  D= P->rD;
	seq = SeqPtr(E); 
	for(n=0;n<=N; n++){
	   NEW(S[n],M+3,float); NEW(D[n],M+3,float); NEW(I[n],M+3,float); 
	}
	for(score=0, n=1; n<= N; n++) {
	   ins = P->is[n] + P->ie[n]; ie = P->ie[n];
	   del = P->ds[n] + P->de[n]; de = P->de[n];
	   for(m=1; m<= M; m++) {
		s = S[n-1][m-1] + P->score[n][seq[m]];
		D[n][m] = MAXIMUM(float,S[n-1][m]-del,D[n-1][m]-de);
		I[n][m] = MAXIMUM(float,S[n][m-1]-ins,I[n][m-1]-ie);
		if(s < D[n][m]) s = D[n][m]; 
		if(s < I[n][m]) s = I[n][m]; 
		if(s > 0){
			S[n][m] = s; 
			if(s > score) score = s; 
		}
	   }
	}
	for(n=0; n<=N; n++){ free(I[n]);free(D[n]);free(S[n]); }
	return (double) score;
}
