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

void	ndl_typ::PutDPMatrix(FILE *fptr,Int4 **MAT,Int4 **TB,Int4 **DEL,Int4 **INS,
		Int4 sq_len,unsigned char *seq2,Int4 nmod,smx_typ *M,
		Int4 blk, Int4 startM, Int4 endM, Int4 startS, Int4 endS)
//******************** debug using matrix output. *********************
{
	Int4	i,j,k,m,hmm_len;
	assert(blk <= nmod);
	cma_typ	cma=CMA;
	for(hmm_len=0, m = 1; m <= nmod; m++){ hmm_len += LenSMatrix(M[m]); }
	unsigned char *seq1; NEW(seq1, hmm_len+3, unsigned char);
	a_type	A=AlphabetCMSA(cma);
  	if(endM > LenSMatrix(M[blk])) endM=LenSMatrix(M[blk]);

#if 0	// for multiple blocks where start is 
	{ for(Int4 x=0, m = 1; m <= nmod; m++){ 
        	MinSegSMatrix(seq1+x,M[m]); // NEED ONLY FOR OUTPUT.
		x += LenSMatrix(M[m]);
	  } e_type E=MkSeq("min_seq for debugging",hmm_len,seq1);
	  PutSeq(stderr,E,A);
	}
	{ for(Int4 x=0, m = 1; m <= nmod; m++){ 
        	MaxSegSMatrix(seq1+x,M[m]); // NEED ONLY FOR OUTPUT.
		x += LenSMatrix(M[m]);
	  } e_type E=MkSeq("max_seq for debugging",hmm_len,seq1);
	  PutSeq(stderr,E,A);
	}
  { Int4 s;
        for(s=m=1; m <= nmod; m++){
            for(Int4 j=1; j<= LenSMatrix(M[m]); j++){ 
		if(m==blk && j==i){ startM=s; break; }  
		else s++;
	    } if(m==blk) break;
        }
  }
#endif
        /** get concensus sequence for smatrix **/
        MaxSegSMatrix(seq1,M[blk]); // NEED ONLY FOR OUTPUT.

	// taken from PutSMatrix():
        fprintf(fptr,"\n    ");
        for(k=0; k<=nAlpha(A); k++) fprintf(fptr,"   %c ",AlphaChar(k,A));
        fprintf(fptr,"\n");
        for(j=startM; j <= endM; j++){
           fprintf(fptr,"%2d: ",j);
           for(k=0; k<=nAlpha(A); k++){ fprintf(fptr,"%4d ",
		(Int4) floor(((double)ValSMatrix(j,k,M[blk])/100.0)+0.5));
           } fprintf(fptr,"\n");
        } fprintf(fptr,"\n");

        fprintf(stderr,"SMX[i][j]:\n");
        this->put_dp_matrixSMX(fptr,sq_len, M[blk],seq1,seq2,A,startM,endM,startS,endS);
        fprintf(stderr,"MAT[i][j]:\n");
        this->put_dp_matrixSW(fptr,MAT, TB, sq_len, M[blk],seq1,seq2,A,startM,endM,startS,endS);
        fprintf(stderr,"DEL[i][j]:\n");
        this->put_dp_matrixSW(fptr,DEL, TB, sq_len, M[blk],seq1,seq2,A,startM,endM,startS,endS);
        fprintf(stderr,"INS[i][j]:\n");
        this->put_dp_matrixSW(fptr,INS, TB, sq_len, M[blk],seq1,seq2,A,startM,endM,startS,endS);
        fprintf(stderr,"TB[i][j]:\n");
        this->put_dp_matrixSW(fptr,TB, TB, sq_len, M[blk],seq1,seq2,A,startM,endM,startS,endS);
}

void    ndl_typ::put_dp_matrixSW(FILE *fptr,Int4 **D, Int4** T, Int4 len, smx_typ M,
        unsigned char *profile, unsigned char *seq, a_type A,
        Int4 startM,Int4 endM,Int4 start, Int4 end)
{
        Int4    r,i,lenM,v;

        lenM=LenSMatrix(M);
	assert(startM > 0 && endM <= lenM && startM <= endM);
        // start = MAXIMUM(Int4,start,0); end = MINIMUM(Int4,len,end);
        // start = MINIMUM(Int4,start,end); end = MAXIMUM(Int4,start,end);
        fprintf(stderr,"     |");
        for(i=startM-1; i<=endM; i++) fprintf(stderr,"  %c  |",AlphaChar(profile[i],A));
        fprintf(stderr,"\n     |");
        for(i=startM-1; i <= endM; i++) fprintf(stderr,"%3d  |",i); fprintf(stderr,"\n");
        for(r=start-1; r<=end; r++) {
           fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
           for(i=startM-1; i<=endM; i++) {
		v=D[i][r];
                if(v <= -infinity) fprintf(stderr," -inf");
                else { if(D != T) v = (Int4) floor(((double)v/100.0)+0.5); 
		fprintf(stderr,"%5d",v); }
                if(T[i][r] < 0) fprintf(stderr,"^");
                else if(T[i][r] == 0) fprintf(stderr,"\\");
                else fprintf(stderr,"<");
           } fprintf(stderr,"\n");
        } fprintf(stderr,"\n");
}

void    ndl_typ::put_dp_matrixSMX(FILE *fptr, Int4 len, smx_typ M, unsigned char *profile,
        unsigned char *seq, a_type A, Int4 startM,Int4 endM,Int4 start, Int4 end)
{
        Int4    r,m,i,lenM;

        lenM=LenSMatrix(M);
	assert(startM > 0 && endM <= lenM && startM <= endM);
        // end = MINIMUM(Int4,len,end); start = MAXIMUM(Int4,start,0);
        // end = MAXIMUM(Int4,start,end); start = MINIMUM(Int4,start,end);
        fprintf(stderr,"     |");
        for(i=startM,m=startM; m<=endM; m++,i++) fprintf(stderr,"  %c  |",
                                AlphaChar(profile[i],A));
        fprintf(stderr,"\n     |");
        for(m=startM; m<=endM; m++) fprintf(stderr,"%3d  |",m); fprintf(stderr,"\n");
        for(r=start; r<=end; r++) {
           fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
           for(m=startM; m<=endM; m++) {
                fprintf(stderr,"%5d|",(Int4) floor(((double)ValSMatrix(m,seq[r],M)/100.0)+0.5));
           } fprintf(stderr,"\n");
        } fprintf(stderr,"\n");
}

