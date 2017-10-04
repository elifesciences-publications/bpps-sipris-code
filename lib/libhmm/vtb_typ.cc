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

#include "vtb_typ.h"

#define DEBUG	0
#define NEG_INF_VTB	-99999999

vtb_typ::vtb_typ(Int4 MaxProfLen, Int4 MaxSeqLen)
{ init(MaxProfLen,MaxSeqLen); }

void vtb_typ::init(Int4 MaxProfLen, Int4 MaxSeqLen)
{
	Int4 	i;
	maxProfLen = MaxProfLen;
	maxSeqLen = MaxSeqLen;

//	Int4 *tempLong;
	MEW(tempLong,3*(maxSeqLen+3)*(maxProfLen+7),Int4);

	NEWP(M,maxSeqLen+2,Int4);NEWP(D,maxSeqLen+2,Int4);
	NEWP(I,maxSeqLen+2,Int4);NEWP(traceM,maxSeqLen+2,char);
	NEWP(traceD,maxSeqLen+2,char);NEWP(traceI,maxSeqLen+2,char);

	M[0] = tempLong;
	for(i=1;i<=maxSeqLen+1;i++) M[i] = M[i-1] + maxProfLen+2;
	D[0] = M[maxSeqLen+1] + maxProfLen+2;
	for(i=1;i<=maxSeqLen+1;i++) D[i] = D[i-1] + maxProfLen+2;
	I[0] = D[maxSeqLen+1] + maxProfLen+2;
	for(i=1;i<=maxSeqLen+1;i++) I[i] = I[i-1] + maxProfLen+2;

	N = I[maxSeqLen+1] + maxProfLen+2;
	B = N + maxSeqLen+2;
	E = B + maxSeqLen+2;
	J = E + maxSeqLen+2;
	C = J + maxSeqLen+2;
	jumpE = C + maxSeqLen+2;

//	char *tempChar;
	MEW(tempChar,3*(maxSeqLen+3)*(maxProfLen+7),char);

	traceM[0] = tempChar;
	for(i=1;i<=maxSeqLen+1;i++) traceM[i] = traceM[i-1] + maxProfLen+2;
	traceD[0] = traceM[maxSeqLen+1] + maxProfLen+2;
	for(i=1;i<=maxSeqLen+1;i++) traceD[i] = traceD[i-1] + maxProfLen+2;;
	traceI[0] = traceD[maxSeqLen+1] + maxProfLen+2;
	for(i=1;i<=maxSeqLen+1;i++) traceI[i] = traceI[i-1] + maxProfLen+2;;

	traceB = traceI[maxSeqLen+1] + maxProfLen+2;
	traceC = traceB + maxSeqLen+2;
	traceJ = traceC + maxSeqLen+2;
	traceN = traceJ + maxSeqLen+2;
}

void vtb_typ::Free( )
{
free(tempLong);free(tempChar);
if(M!=NULL) free(M);
if(D!=NULL) free(D);
if(I!=NULL) free(I);
if(traceM!= NULL) free(traceM);
if(traceD!= NULL) free(traceD);
if(traceI!= NULL) free(traceI);

#if 0
	Int4 i;
	if(M!=NULL){ for(i=0;i<=maxSeqLen+1;i++) free(M[i]); free(M); }
	if(D!=NULL){ for(i=0;i<=maxSeqLen+1;i++) free(D[i]); free(D); }
	if(I!=NULL){ for(i=0;i<=maxSeqLen+1;i++) free(I[i]); free(I); }
	if(traceM!= NULL){ for(i=0;i<=maxSeqLen+1;i++) free(traceM[i]); free(traceM); }
	if(traceD!= NULL){ for(i=0;i<=maxSeqLen+1;i++) free(traceD[i]); free(traceD); }
	if(traceI!= NULL){ for(i=0;i<=maxSeqLen+1;i++) free(traceI[i]); free(traceI); }
	if(traceJ != NULL) { free(traceJ); }
	if(traceB != NULL) { free(traceB); }
	if(traceC != NULL) { free(traceC); }
	if(traceN != NULL) { free(traceN); }
	if(jumpE != NULL) { free(jumpE); }
#endif
}

char    *vtb_typ::Viterbi(e_type sE,shm_typ *shm,Int4 *num_rpts,Int4 *total_score, Int4 *scores)
{
	Int4	oper_length,Nrpts,tot_score;
	Int4 raw_score = Viterbi(sE,shm->GetLength(),shm->Matrix(),
		shm->InsMatrix(),shm->MatToMat(),shm->MatToIns(),
                shm->MatToDel(),shm->InsToMat(),shm->InsToIns(),
		shm->DelToMat(),shm->DelToDel(),shm->BegToMat(),
                shm->MatToEnd(),shm->NB(),shm->NN(),shm->EC(),
		shm->EJ(),shm->CT(),shm->CC(),shm->JB(),shm->JJ());
	char *operation=Traceback(sE,shm,&oper_length,&Nrpts,&tot_score,scores);
	*num_rpts=Nrpts; *total_score = tot_score;
	return operation;
}

Int4 vtb_typ::Viterbi(e_type sE,shm_typ *shm)
// Warning: need to use renormalized input hmm.
// Eventually modify this to eliminate traceback matrix.
{
	Int4 score = FastViterbi(sE,shm->GetLength(),shm->Matrix(),
		shm->InsMatrix(),shm->MatToMat(),shm->MatToIns(),
                shm->MatToDel(),shm->InsToMat(),shm->InsToIns(),
		shm->DelToMat(),shm->DelToDel(),shm->BegToMat(),
                shm->MatToEnd(),shm->NB(),shm->NN(),shm->EC(),
		shm->EJ(),shm->CT(),shm->CC(),shm->JB(),shm->JJ());
	return score;
}

Int4 vtb_typ::Viterbi(e_type sE,Int4 prof_len,Int4 **mat_emit,Int4 **ins_emit,
		Int4 *m2m,Int4 *m2i, Int4 *m2d,Int4 *i2m,Int4 *i2i,Int4 *d2m,
		Int4 *d2d,Int4 *b2m,Int4 *m2e,Int4 n2b,Int4 n2n,Int4 e2c,
		Int4 e2j,Int4 c2t,Int4 c2c,Int4 j2b,Int4 j2j)
{
	Int4 		score;
	Int4 		seq_len = LenSeq(sE);
	Int4 		i,j,im1,jm1,Temp,max;
	Int4 		*Mim1,*Iim1,*Dim1;
	Int4 		Jim1,Bim1,Cim1;
	Int4 		bb,cc,dd,ee,ii,jj,mm,nn;
	Int4 		best_j;
	unsigned char 	*seq = SeqPtr(sE),res;

	N[0] = 0;
	traceN[0] = 'S';
	for(i=1;i<=seq_len-1;i++) { N[i] = N[i-1] + n2n; traceN[i] = 'N'; }
	N[seq_len] = NEG_INF_VTB;
	B[0] = n2b;
	traceB[0] = 'S';
	J[0] = C[0] = NEG_INF_VTB;
	for(i=0;i<=seq_len;i++) { 
		M[i][0] = D[i][0] = D[i][1] = I[i][0] = I[i][prof_len] =  NEG_INF_VTB; 
	}
	for(j=1;j<=prof_len;j++){
	  M[0][j] = D[0][j] = I[0][j] = I[1][j] = I[seq_len][j] = NEG_INF_VTB; 
	}
#if DEBUG
for(j=1;j<=prof_len;j++) {
printf("b2m[%d]=%d  m2e[%d]=%d\n",j,b2m[j],j,m2e[j]);
}
#endif
	for(i=1,im1=0;i<=seq_len;i++,im1++) {
	  res = seq[i];
	  Mim1 = M[im1], Iim1 = I[im1], Dim1 = D[im1];
	  Jim1 = J[im1]; Bim1 = B[im1]; Cim1 = C[im1];
          for(j=1,jm1=0;j<=prof_len;j++,jm1++) {
		mm = Mim1[jm1] + m2m[jm1]; 
		ii = Iim1[jm1] + i2m[jm1];
		dd = Dim1[jm1] + d2m[jm1];
		bb = Bim1 + b2m[j];
		if(mm>=ii && mm>=dd && mm>=bb){
		  M[i][j] = mat_emit[j][res] + mm; traceM[i][j] = 'M'; 
		} else if (ii>=dd && ii>=bb) {
		  M[i][j] = mat_emit[j][res] + ii; traceM[i][j] = 'I'; 
		} 
		else if (dd>=bb) { M[i][j] = mat_emit[j][res] + dd; traceM[i][j] = 'D'; }
		else { M[i][j] = mat_emit[j][res] + bb; traceM[i][j] = 'B'; }
#if DEBUG
printf("mat_emit[%d][%d]=%d\n",j,res,mat_emit[j][res]);
#endif

#if DEBUG
printf("mm=%d ii=%d dd=%d bb=%d\n",mm,ii,dd,bb);
#endif
		mm = Mim1[j] + m2i[j];
		ii = Iim1[j] + i2i[j];
		if(mm>=ii) { I[i][j] = ins_emit[j][res] + mm; traceI[i][j] = 'M'; }
		else { I[i][j] = ins_emit[j][res] + ii; traceI[i][j] = 'I'; }
#if DEBUG
printf("mm=%d ii=%d\n",mm,ii);
#endif
		mm = M[i][jm1] + m2d[jm1];
		dd = D[i][jm1] + d2d[jm1];
		if(mm>=dd) { D[i][j] = mm; traceD[i][j] = 'M'; }
		else { D[i][j] = dd; traceD[i][j] = 'D'; }
#if DEBUG
printf("mm=%d dd=%d\n",mm,dd);
printf("traceM[%d][%d]=%c traceI[%d][%d]=%c traceD[%d][%d]=%c\n",
		i,j,traceM[i][j],i,j,traceI[i][j],i,j,traceD[i][j]);
printf("M[%d][%d]=%d I[%d][%d]=%d D[%d][%d]=%d\n",
		i,j,M[i][j],i,j,I[i][j],i,j,D[i][j]);
#endif
          }
	  for(max=INT4_MIN,j=1;j<=prof_len;j++) { 
	    if((Temp = M[i][j] + m2e[j]) > max) { max = Temp; best_j = j; }
	  }
          E[i] =  Temp;
#if DEBUG
printf("E[%d]=%d ",i,E[i]);
#endif
	  jumpE[i] = best_j;
	  ee = E[i] + e2j;
	  jj = Jim1 + j2j;
	  if(ee>=jj) { J[i] = ee; traceJ[i] = 'E'; }
	  else { J[i] = jj; traceJ[i] = 'J'; }
#if DEBUG
printf("J[%d]=%d ",i,J[i]);
#endif
	  nn = N[i] + n2b;
	  jj = J[i] + j2b;
	  if(nn>=jj) { B[i] = nn; traceB[i] = 'N'; }
	  else { B[i] = jj; traceB[i] = 'J'; }
#if DEBUG
printf("B[%d]=%d ",i,B[i]);
#endif
	  ee = E[i] + e2c;
	  cc = Cim1 + c2c;
	  if(ee>=cc) { C[i] = ee; traceC[i] = 'E'; }
	  else { C[i] = cc; traceC[i] = 'C'; }
#if DEBUG
printf("C[%d]=%d\n",i,C[i]);
#endif
	}
	score = C[seq_len] + c2t;
//printf("--- SCORE=%d ------\n",score);
	return score;
}
#if 1
Int4 vtb_typ::FastViterbi(e_type sE,Int4 prof_len,Int4 **mat_emit,Int4 **ins_emit,
		Int4 *m2m,Int4 *m2i, Int4 *m2d,Int4 *i2m,Int4 *i2i,Int4 *d2m,
		Int4 *d2d,Int4 *b2m,Int4 *m2e,Int4 n2b,Int4 n2n,Int4 e2c,
		Int4 e2j,Int4 c2t,Int4 c2c,Int4 j2b,Int4 j2j)
{
	Int4 		score;
	Int4 		seq_len = LenSeq(sE);
	Int4 		i,j,im1,jm1,Temp,max;
	Int4 		*Mim1,*Iim1,*Dim1;
	Int4 		Jim1,Bim1,Cim1;
	Int4 		bb,cc,dd,ee,ii,jj,mm,nn;
	Int4 		best_j;
	unsigned char 	*seq = SeqPtr(sE),res;

	N[0] = 0;
	traceN[0] = 'S';
	for(i=1;i<=seq_len-1;i++) N[i] = N[i-1] + n2n;
	N[seq_len] = NEG_INF_VTB;
	B[0] = n2b;
	J[0] = C[0] = NEG_INF_VTB;
	for(i=0;i<=seq_len;i++) { 
		M[i][0] = D[i][0] = D[i][1] = I[i][0] = I[i][prof_len] =  NEG_INF_VTB; 
	}
	for(j=1;j<=prof_len;j++){
	  M[0][j] = D[0][j] = I[0][j] = I[1][j] = I[seq_len][j] = NEG_INF_VTB; 
	}
#if DEBUG
for(j=1;j<=prof_len;j++) {
printf("b2m[%d]=%d  m2e[%d]=%d\n",j,b2m[j],j,m2e[j]);
}
#endif
	for(i=1,im1=0;i<=seq_len;i++,im1++) {
	  res = seq[i];
	  Mim1 = M[im1], Iim1 = I[im1], Dim1 = D[im1];
	  Jim1 = J[im1]; Bim1 = B[im1]; Cim1 = C[im1];
          for(j=1,jm1=0;j<=prof_len;j++,jm1++) {
		mm = Mim1[jm1] + m2m[jm1]; 
		ii = Iim1[jm1] + i2m[jm1];
		dd = Dim1[jm1] + d2m[jm1];
		bb = Bim1 + b2m[j];
		if(mm>=ii && mm>=dd && mm>=bb){
		  M[i][j] = mat_emit[j][res] + mm;
		} else if (ii>=dd && ii>=bb) {
		  M[i][j] = mat_emit[j][res] + ii; 
		} 
		else if (dd>=bb) M[i][j] = mat_emit[j][res] + dd;
		else M[i][j] = mat_emit[j][res] + bb;
#if DEBUG
printf("mat_emit[%d][%d]=%d\n",j,res,mat_emit[j][res]);
#endif

#if DEBUG
printf("mm=%d ii=%d dd=%d bb=%d\n",mm,ii,dd,bb);
#endif
		mm = Mim1[j] + m2i[j];
		ii = Iim1[j] + i2i[j];
		if(mm>=ii) I[i][j] = ins_emit[j][res] + mm; 
		else I[i][j] = ins_emit[j][res] + ii;
#if DEBUG
printf("mm=%d ii=%d\n",mm,ii);
#endif
		mm = M[i][jm1] + m2d[jm1];
		dd = D[i][jm1] + d2d[jm1];
		if(mm>=dd) D[i][j] = mm;
		else D[i][j] = dd;
#if DEBUG
printf("mm=%d dd=%d\n",mm,dd);
printf("M[%d][%d]=%d I[%d][%d]=%d D[%d][%d]=%d\n",
		i,j,M[i][j],i,j,I[i][j],i,j,D[i][j]);
#endif
          }
	  for(max=INT4_MIN,j=1;j<=prof_len;j++) { 
	    if((Temp = M[i][j] + m2e[j]) > max) { max = Temp; best_j = j; }
	  }
          E[i] =  Temp;
#if DEBUG
printf("E[%d]=%d ",i,E[i]);
#endif
	  jumpE[i] = best_j;
	  ee = E[i] + e2j;
	  jj = Jim1 + j2j;
	  if(ee>=jj) J[i] = ee;
	  else J[i] = jj; 
#if DEBUG
printf("J[%d]=%d ",i,J[i]);
#endif
	  nn = N[i] + n2b;
	  jj = J[i] + j2b;
	  if(nn>=jj) B[i] = nn; 
	  else B[i] = jj; 
#if DEBUG
printf("B[%d]=%d ",i,B[i]);
#endif
	  ee = E[i] + e2c;
	  cc = Cim1 + c2c;
	  if(ee>=cc) C[i] = ee; 
	  else C[i] = cc; 
#if DEBUG
printf("C[%d]=%d\n",i,C[i]);
#endif
	}
	score = C[seq_len] + c2t;
//printf("C[%d]=%d c2t=%d\n",seq_len,C[seq_len],c2t);
//printf("--- SCORE=%d ------\n",score);
	return score;
}
#endif
char *vtb_typ::Traceback(e_type sE,shm_typ *shm,Int4 *oper_length,Int4 *num_rpts,
	Int4 *total_score, Int4 *scores)
{
	char 		*operation,*back_operation;
	Int4 		*back_scores;
	Int4 		seq_len = LenSeq(sE);
	char 		state = traceC[seq_len];
	Int4 		prof_len = shm->GetLength();
	Int4 		i = seq_len, j = prof_len, l;
	Int4		k = 0, s, end_score = 0, m = 1;
	Int4		max_path,MaxRpts,rpt=1;
	Int4 		bits = shm->Bits();
	Int4 		**mat_emit = shm->Matrix();
	Int4 		**ins_emit = shm->InsMatrix();
	Int4 		nalpha = nAlpha(shm->GetAlphabet());
	double		*p,*tot_p,tp,sum;
	unsigned char 	*seq = SeqPtr(sE);
	Int4 		*nule = shm->NULE();
	double 		*nl;
	Int4 		*total_emiss_dist, **rep_emiss_dist;
	NEW(p,nalpha+2,double);
	NEW(tot_p,nalpha+2,double);
	NEW(nl,nalpha+2,double);
	for(l=1;l<=nalpha;l++) nl[l] = (1./nalpha)*exp2(nule[l]/((double) bits));
	*total_score = C[seq_len] + shm->CT();

#if 0
	Predicted maximum back_operation length is:
		seqlen [nmic]
		+ MaxRpts*(prof_len - 1) [md]
		+ 2*MaxRpts [BE] + 2 [ST].
	Predicted maximum number of repeats is MaxRpts = seq_len.
	Therefore the maximum back_operation length ==
		seqlen + seqlen*(prof_len - 1) + 2*seqlen + 2.
	for seqlen == 30,000 and hmmlen == 1000:
		= 30060002 bytes or about 30 Mbytes.
	But...
	  in practice assume that Nrpts = 2*seqlen/prof_len;
	Then maximum array for seqlen == 30,000 and hmmlen == 1000:
		MaxRpts = (2*seqlen/prof_len) = 60.
		max_path = seqlen + MaxRpts*(prof_len - 1) + 2*MaxRpts + 2.
			= 90062 or about 90 Kbytes.
#endif
	MaxRpts=(2*seq_len/prof_len)+5;
	if(MaxRpts < 2) MaxRpts=2;
	NEW(total_emiss_dist,nalpha+2,Int4);
	NEWP(rep_emiss_dist,MaxRpts+2,Int4);
	for(l=1;l<=MaxRpts;l++) NEW(rep_emiss_dist[l],nalpha+2,Int4);

	max_path=seq_len + MaxRpts*(prof_len) + 2*MaxRpts + 2 + 10000;// need fix A.P.
        MEW(back_operation,max_path+3,char);
	MEW(back_scores,seq_len+3,Int4);
	back_operation[k++] = 'T';
	*num_rpts = 0;
        while(state != 'S'){
           switch(state){
		case 'M': back_operation[k++] = 'm';
			  state = traceM[i][j];
			  for(l=1;l<=nalpha;l++) {
			  	p[l] += (tp = nl[l]*exp2(mat_emit[j][l]/((double) bits)));
			  	tot_p[l] += tp;
			  }
			  i--; j--; break;
		case 'D': back_operation[k++] = 'd';
			  state = traceD[i][j];
			  j--; break;
		case 'I': back_operation[k++] = 'i';
			  for(l=1;l<=nalpha;l++) {
			  	p[l] += (tp = nl[l]*exp2(ins_emit[j][l]/((double) bits)));
			  	tot_p[l] += tp;
			  }
			  state = traceI[i][j];
			  i--; break;
		case 'B': back_scores[m] = end_score - B[i] + shm->NB() + i*(shm->NN());
			  for(s=1;s<=j;s++) back_operation[k++]='d';
			  back_operation[k++] = 'B';
			  state = traceB[i];
			  for(sum = 0.,l=1;l<=nalpha;l++) {
				sum += p[l];
			  }
			  for(l=1;l<=nalpha;l++) {
			  	rep_emiss_dist[m][l] = (Int4) floor(0.5 + bits*log2(p[l]/(sum*nl[l])));
				p[l] = 0.;
			  } 
			  m++;
			  j=prof_len; break;
		case 'N': if((state = traceN[i]) == 'N') {
                                back_operation[k++] = 'n'; 
                                i--; 
			  } break;
		case 'J': if((state = traceJ[i]) == 'J') {
                                back_operation[k++] = 'j';
                                i--; 
			  } break;
		case 'E': *num_rpts += 1;
			  end_score = E[i] + shm->EC() + (seq_len-i)*(shm->CC()) + shm->CT();
			  back_operation[k++]='E';
			  j = jumpE[i];
			  for(s=1;s<=prof_len-j;s++) back_operation[k++]='d';
			  state = 'M';
			  break;
		case 'C': back_operation[k++] = 'c';
			  i--; state = traceC[i];
			  break; 
		default: 
		   fprintf(stderr,"state = '%c'(%d)\n",state,state);
		   print_error("Traceback 1: error in the traceback matrix!");
           }
        }
	for(sum = 0.,l=1;l<=nalpha;l++) {
		sum += tot_p[l];
	}

	for(l=1;l<=nalpha;l++) {
		total_emiss_dist[l] = (Int4) floor(0.5 + bits*log2(tot_p[l]/(sum*nl[l])));
	}
	back_operation[k] = 'S';
	*oper_length = k+1;
	MEW(operation,*oper_length+2,char);
	for(s=0;s <= k;s++) operation[s] = back_operation[k-s];
	operation[s] = '\0';
 	i = seq_len, j = prof_len;
	state = traceC[seq_len]; 
	m=1;
	Int4 total_adj_score = 0, back_adj_score = 0;
        while(state != 'S'){
           switch(state){
		case 'M': state = traceM[i][j];
			  total_adj_score += total_emiss_dist[seq[i]];
			  back_adj_score += rep_emiss_dist[m][seq[i]];
			  i--; j--; break;
		case 'D': state = traceD[i][j];
			  j--; break;
		case 'I': state = traceI[i][j];
			  total_adj_score += total_emiss_dist[seq[i]];
			  back_adj_score += rep_emiss_dist[m][seq[i]];
			  i--; break;
		case 'B': state = traceB[i];
			  back_adj_score -= 8*bits;
			  back_scores[m] -= ILogsum(0, back_adj_score, bits);
			  back_adj_score = 0;
			  m++;
			  j=prof_len; break;
		case 'N': if((state = traceN[i]) == 'N') {
                                i--; 
			  } break;
		case 'J': if((state = traceJ[i]) == 'J') {
                                i--; 
			  } break;
		case 'E': j = jumpE[i];
			  state = 'M';
			  break;
		case 'C': i--; state = traceC[i];
			  break; 
		default: 
		   fprintf(stderr,"state = '%c'(%d)\n",state,state);
		   print_error("Traceback 2: error in the traceback matrix!");
           }
        }
	total_adj_score -= 8*bits;
	*total_score -= ILogsum(0,total_adj_score,bits);
	for(s=1;s<=*num_rpts;s++) scores[s] = back_scores[*num_rpts-s+1]; 
	free(back_operation); free(back_scores);
	free(p); free(tot_p); free(nl); free(total_emiss_dist); 
	for(l=1;l<=MaxRpts;l++) free(rep_emiss_dist[l]);
	free(rep_emiss_dist);
        return operation;
}
