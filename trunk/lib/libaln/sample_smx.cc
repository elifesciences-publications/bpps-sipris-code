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

#include "smatrix.h"

////////////////////// SAMPLING METHODS //////////////////////////
char	*SampleOperationSeqAlnSMatrixSW(int a, int b, int n2, unsigned char *seq2, 
	UInt4 offset, int nmod, smx_typ *M, int **gapscore, int *J, Int4 *alnscore)
// sample an alignment returned as a set of match,insert, & delete operations.
{ return sample_gapped_aln_seq_smatrixSW(a,b,n2,seq2,nmod,M,gapscore,J,alnscore); }

Int4	PutSampledSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	UInt4 offset, Int4 nmod, smx_typ *M, Int4 **gapscore)
{
	char	*operation;
	Int4	J,alnscore;

	operation=sample_gapped_aln_seq_smatrixSW(a,b,n2,seq2,nmod,M,gapscore,&J,&alnscore);
	// fprintf(fp,"Put: %s\n\n",operation);
	put_seqaln_smatrixSW(fp, operation, n2, seq2, offset,  J, nmod, M);
	// fprintf(stderr,"operations = %s\n",operation);
	free(operation);
	return alnscore;
}

char	*SampleGapOperationsSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore)
{
	Int4    s,N; 
	return sample_gapped_aln_seq_smatrixSW(a,b,len,seq2,nmod,M,gapscore,&N,&s);
}

#include "random.h"
static void	put_dp_sampled_matrixSW(double **D, Int4 len, smx_typ M,
	unsigned char *profile, unsigned char *seq, a_type A,
	Int4 startM,Int4 start, Int4 end, Int4 startMtf, Int4 endMtf)
{
	Int4	r,m,i,lenM;
	
	end = MINIMUM(Int4,len,end);
	start = MAXIMUM(Int4,start,0);
	end = MAXIMUM(Int4,start,end);
	start = MINIMUM(Int4,start,end);
	lenM=LenSMatrix(M);
	fprintf(stderr,"     |");
	for(i=startM-1,m=0; m<=lenM; m++,i++){
	    if(m >= startMtf && m <= endMtf)
		fprintf(stderr,"  %c  |", AlphaChar(profile[i],A));
	}
	fprintf(stderr,"\n     |");
	for(m=0; m<=lenM; m++){
	    if(m >= startMtf && m <= endMtf) fprintf(stderr,"%3d  |",m);
	}
	fprintf(stderr,"\n");
	for(r=start-1; r<=end; r++) {	
	   fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
	   for(i=startM-1,m=0; m<=lenM; m++,i++) {
	     if(m >= startMtf && m <= endMtf){
		fprintf(stderr,"%5.1g",D[i][r]);
		fprintf(stderr," ");
	     }
	   }
	   fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
}

void	foo_sample(double *smx, Int4 n2, double *mat, double *ins,
	double *del,double *m1del,double *m1mat,double *m1ins, unsigned char *seq2,
	double L_io, double L_ie, double L_do, double L_de)
	
{
	register Int4	j,jm1;

	for(jm1=0,j=1; j<= n2; j++,jm1++) {  
		mat[j] = smx[seq2[j]]*(m1mat[jm1] + m1del[jm1] + m1ins[jm1]);
		del[j] = L_do*m1mat[j] + L_de*m1del[j];
		ins[j] = L_io*mat[jm1] + L_ie*ins[jm1];
	}
}

void	foo_sample_edge(double *smx, Int4 n2, double *mat, double *ins,
	double *del,double *m1del,double *m1mat,double *m1ins, unsigned char *seq2,
	double L_do, double L_de)
{
	register Int4	j,jm1;

	for(jm1=0,j=1; j<= n2; j++,jm1++) {  
		mat[j] = smx[seq2[j]]*(m1mat[jm1] + m1del[jm1] + m1ins[jm1]);
		del[j] = L_do*m1mat[j] + L_de*m1del[j];
		ins[j] = mat[jm1] + del[jm1] + ins[jm1];
	}
}

/////////////// Sampling version of SW-ALIGNMENT WITH GAP FUNCTION /////////////
//
//  Work in progress.  Not yet finished...
//
///////////////////////////////////////////////////////////////////
char	*sample_gapped_aln_seq_smatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *J, Int4 *alnscore)
{
	Int4	m,i,j,k,r1,t,v,n1,max_i,max_j,lenM;
	Int4	t0,j2,g,gs,g_opt,w=a+b,J0,**T,**TB,score;
	double	**MAT,**DEL,**INS,gINS,s,s0,s1,**SMX,*smx;
	double	**DELX, **INSX;
	double	*delx, *insx,*m1delx,*m1insx;
	double	sum,rand_num,*mat,*ins,*del,*m1del,*m1mat,*m1ins;
	double	L_io,L_do,L_ie,L_de;
	unsigned char	*seq1;
	short	*mtf;
	a_type	A = SMatrixA(M[1]);
	char	*operation;
	const double	adjust=1e-100;
	BooLean	flag=TRUE;
	Int4	j0,jm1;

	if(gapscore != NULL) gapscore--; // want gapscore[m-1] for block m.
	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); NEWP(TB, nmod+3, Int4); 
	for(j=0, m=1; m <= nmod; m++){
	    NEW(TB[m], n2+3, Int4); 
	    MaxSegSMatrix(seq1+j, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ j++; mtf[j] = m; }
	}
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,double*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,double*); MEW(INS,n1+3,double*);
	MEW(DELX,n1+3,double*); MEW(INSX,n1+3,double*);
	MEW(SMX,n1*3,double*);
	for(i=0; i<= n1+1; i++) { 
#if 0
		MEW(MAT[i],n2+3,double); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,double); MEW(INS[i],n2+3,double); 
		MEW(DELX[i],n2+3,double); MEW(INSX[i],n2+3,double); 
#endif
#if 1
		NEW(MAT[i],n2+3,double); NEW(T[i],n2+3,Int4); 
		NEW(DEL[i],n2+3,double); NEW(INS[i],n2+3,double); 
		NEW(DELX[i],n2+3,double); NEW(INSX[i],n2+3,double); 
#endif
		NEW(SMX[i],nAlpha(A)+3,double);
	}
#if 0   // old sampling temperature (1 bit)
	for(k=1,i=1,m=lenM=0; i<= n1; i++) {
	   if(k > lenM) { k=1; m++; lenM=LenSMatrix(M[m]); }
	   for(j=0; j<=nAlpha(A); j++){
		SMX[i][j] = pow(2.0,(double)ValSMatrix(k,j,M[m]));
	   } k++;
	}
	L_do = pow(2.0,(double)-w); L_de = pow(2.0,(double)-b);
	L_io = pow(2.0,(double)-w); L_ie = pow(2.0,(double)-b);
#endif 
#if 1   // new sampling temperature (0.5 bit)
	for(k=1,i=1,m=lenM=0; i<= n1; i++) {
	   if(k > lenM) { k=1; m++; lenM=LenSMatrix(M[m]); }
	   for(j=0; j<=nAlpha(A); j++){
		SMX[i][j] = pow(2.0,2.0*(double)ValSMatrix(k,j,M[m]));
	   } k++;
	}
	L_do = pow(2.0,2.0*(double)-w); L_de = pow(2.0,2.0*(double)-b);
	L_io = pow(2.0,2.0*(double)-w); L_ie = pow(2.0,2.0*(double)-b);
#endif
	// Make GLOBAL ALIGNMENT with respect to profile,
	MAT[0][0] = 1.0; DEL[0][0] = 1.0; INS[0][0] = 1.0;
	for(i=1; i<= n1; i++) {  
		DEL[i][0] = 0.0;
		DELX[i][0] = L_de*DELX[i-1][0];
		// INS[i][0] = L_de*INS[i-1][0];
		INS[i][0] = 0.0;
		INSX[i][0] = 1.0;
		MAT[i][0] = INSX[i-1][0];  // Local on ends of profile 
		T[i][0] = 1; // for full alignment
	}
	for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
		MAT[0][j] = 1.0; T[0][j]=-1; // traceback full alignment 
		DEL[0][j] = INS[0][j] = 1.0;
		DELX[0][j] = INSX[0][j] = 1.0;
	}
//////////// 2. Dynamic programming-like Forward step.  ////////////////
	for(k=i=1,m=lenM=0; i<= n1; i++) {
	   if(k > lenM) { k=1; m++; if(m <= nmod) lenM=LenSMatrix(M[m]); }
#if 0
	   if(k!=lenM) foo_sample(SMX[i], n2, MAT[i], INS[i], DEL[i],DEL[i-1],MAT[i-1],
		INS[i-1], seq2, L_io, L_ie, L_do, L_de);
	   else foo_sample_edge(SMX[i], n2, MAT[i], INS[i], DEL[i],DEL[i-1],MAT[i-1],
                INS[i-1], seq2, L_do, L_de);
#endif
#if 1
	   mat=MAT[i]; m1mat=MAT[i-1];
	   ins=INS[i]; m1ins=INS[i-1];
	   insx=INSX[i]; m1insx=INSX[i-1];
	   del=DEL[i]; m1del=DEL[i-1];
	   delx=DELX[i]; m1delx=DELX[i-1];
	   for(smx=SMX[i],jm1=0,j=1; j<= n2; j++,jm1++) {  
		mat[j] = smx[seq2[j]]
		   *(m1mat[jm1] + m1ins[jm1] + m1insx[jm1] 
			+ m1del[jm1] + m1delx[jm1]);
		if(k==lenM) {
			ins[j] = mat[jm1];
			insx[j] = ins[jm1] + insx[jm1];
		} else {
			ins[j] = L_io*mat[jm1];
			insx[j] = L_ie*(ins[jm1] + insx[jm1]);
		}
		del[j] = L_do*m1mat[j];
		delx[j] = L_de*(m1delx[j] + m1del[j]);
	   }
#endif
	   k++;
	}
#if 0
Int4 startS=55,endS=70,motif=3,startMtf=4,endMtf=13;
	for(i=1; i<=n1; i++) if(mtf[i]==motif) break; 
	fprintf(stderr,"DEL[i][j]:\n");
	put_dp_sampled_matrixSW(DEL, n2, M[motif],seq1,seq2,A,i,startS,endS,
		startMtf,endMtf);
	fprintf(stderr,"DELX[i][j]:\n");
	put_dp_sampled_matrixSW(DELX, n2, M[motif],seq1,seq2,A,i,startS,endS,
		startMtf,endMtf);
	fprintf(stderr,"INS[i][j]:\n");
	put_dp_sampled_matrixSW(INS, n2, M[motif],seq1,seq2,A,i,startS,endS,
		startMtf,endMtf);
	fprintf(stderr,"INSX[i][j]:\n");
	put_dp_sampled_matrixSW(INSX, n2, M[motif],seq1,seq2,A,i,startS,endS,
		startMtf,endMtf);
	fprintf(stderr,"MAT[i][j]:\n");
	put_dp_sampled_matrixSW(MAT, n2, M[motif],seq1,seq2,A,i,startS,endS,
		startMtf,endMtf);
	// fprintf(stderr,"SMX[i][j]:\n");
	// put_dp_matrixSMX(n2, M[motif],seq1,seq2,A,i,startS,endS);
#endif
//////////// 2b. Sample a global alignment with respect to profile. ///////////
	score=0; max_i=n1; max_j=n2;
	char state=' ';
	for(m=nmod+1,k=0,i=max_i,j=max_j; j > 0 && i > 0; ){
	   if(k==0){	// use block gap functions
	      m--; if(m > 0) k=LenSMatrix(M[m]); state = ' ';
	   }
	   // de <- do;  ie <- io only...
	   switch(state){
	     case 'd':   // deletion extention...
	       sum = DEL[i][j] + DELX[i][j]; 
	       rand_num = sum * SampleUniformProb();
	       if((rand_num -= DEL[i][j]) < 0.0) { 
	 	 score-=w; T[i][j]=1; i--; k--; state='D';
	       } else { score-=b; T[i][j]=1; i--; k--; state='d'; } 
	       break;
	     case 'i':   // insertion extention...
	       sum = INS[i][j] + INSX[i][j]; 
	       rand_num = sum * SampleUniformProb();
	       if((rand_num -= INS[i][j]) < 0.0) { 
		 if(k != LenSMatrix(M[m])) score-=w; T[i][j]=-1; j--; state='I';
	       } else {
		 if(k != LenSMatrix(M[m])) score-=b; T[i][j]=-1; j--; state='i';
	       }
	       break;
	     case 'D':   // deletion opening... from match state only
		   score += ValSMatrix(k,seq2[j],M[m]);
		   T[i][j]=0; i--; j--; k--; state='M';
	       break;
	     case 'I':   // insertion opening...from match state only
		   score += ValSMatrix(k,seq2[j],M[m]);
		   T[i][j]=0; i--; j--; k--; state='M';
	       break;
	     default:    // match 
	        sum = DEL[i][j] + MAT[i][j] + INS[i][j] + DELX[i][j] + INSX[i][j];
	        rand_num = sum * SampleUniformProb();
	        if((rand_num -= DEL[i][j]) < 0.0) {
			score-=w; T[i][j]=1; i--; k--; state='D';
		} else if((rand_num -= DELX[i][j]) < 0.0) { 
		    score-=b; T[i][j]=1; i--; k--; state='d';
	        } else if((rand_num -= MAT[i][j]) < 0.0) {
		   score += ValSMatrix(k,seq2[j],M[m]);
		   T[i][j]=0; i--; j--; k--; state='M';
	       } else if((rand_num -= INS[i][j]) < 0.0) { 
		   if(k != LenSMatrix(M[m])){ score-=w; }
		   T[i][j]=-1; j--; state='I';
	       } else { if(k != LenSMatrix(M[m])) score-=b; T[i][j]=-1; j--; state='i';}

	       break;
	   }
	}
//////////// 3. Trace back step. ///////////////////////////////
// So far unchanged below this line; use as common function???
	// operations: 
	// 'i' = insertion in sequence outside of motifs
	// 'I' = insert in sequence within motif (need to delete this)
	//  'M' = match to start of a motif block
	//  'm' = match to other sites in motif
	//  'D' = deletion of sequence within motif
	//  'd' = deletion of sequence outside of motif (not used)
	NEW(operation,n1+n2+3,char);
	m=nmod; k=LenSMatrix(M[m]);
	for(J0=0,j=n2; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
	for(i=n1; i > 0 || j > 0; ){  // full alignment...
	  t0 = T[i][j];
	  do {
	    if(t0 > 0){ t=1; t0--; }
            else if(t0 < 0){ t=-1; t0++; } else { t=0; }
	    switch(t){
		case 0:
		   i--; J0++;
		   if(k==1) operation[J0] = 'M'; else operation[J0] = 'm';
		   if(FALSE && gapscore != NULL && k==1 && m > 1){
			for(g=TB[m][j], j--; g < 0; g++){
			    J0++; operation[J0] = 'i'; j--;
			}
		   } else j--; k--; break;
		case 1:  // Gap ('-') is in sequence; add X's for gaps
			J0++; 
			if(k==1) operation[J0] = 'D'; else operation[J0] = 'd';
			i--;  k--; break;
		case -1:	// Gap ('-') is in profile
			if(m < 1 || k==LenSMatrix(M[m]))
				{ J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			j--; break;
		default: print_error("this should not happen"); break;
	     }
	     if(k==0){ m--; if(m > 0) k=LenSMatrix(M[m]); }
	  } while(t0 != 0);
	}
	while(j > 0) { J0++; operation[J0] = 'i'; j--; }
	/*** 4b. reverse operations ***/
        char *operation2;
        NEW(operation2,J0+3,char); operation2[0] = 'E';
        for(j=J0,i=1; i <= J0; j--,i++){
                operation2[i] = operation[j];
        } free(operation); operation2[i] = 'E'; i++; operation2[i] = 0;
	/*** 5. free allocated memory ***/
	for(m=1; m <= nmod; m++) free(TB[m]); free(TB);
	for(i=0; i<=n1+1; i++) {
		free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);
		free(DELX[i]);free(INSX[i]); free(SMX[i]);
	}
	free(MAT); free(T); free(DEL); free(INS); free(DELX); free(INSX); 
	free(seq1); free(mtf); free(SMX);
	*alnscore = (Int4) score; *J = J0;
// printf("%s\n",operation2);
// printf("gapscore = %d\n",(Int4)gapscore);
        return operation2;
}

