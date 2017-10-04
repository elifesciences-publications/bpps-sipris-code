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


#include "NWprof.h"

char *nwp_algn_smatrix(Int4 seq_len, unsigned char *seq, smx_typ SM, Int4 *io, 
Int4 *ie,Int4 *od, Int4 *de, Int4 *alignscore,Int4 *J, Int4 *start)
{
	int i, j, k, oper_len;
	char *operation, *back_operation;
	char **tMAT, **tINSRT, **tDLT;
	double **MAT, **INSRT, **DLT;
	double m, r, d;
	double **mtch;
	int prof_len = LenSMatrix(SM);
	a_type A = SMatrixA(SM);
	char state;
	double score=0, max;
	int i_max;

	NEWP(MAT,seq_len+3,double);NEWP(INSRT,seq_len+3,double);NEWP(DLT,seq_len+3,double);
        NEWP(tMAT,seq_len+3,char);NEWP(tINSRT,seq_len+3,char);NEWP(tDLT,seq_len+3,char);
	for(i=0;i<=seq_len+2;i++){
		NEW(MAT[i],prof_len+3,double);NEW(INSRT[i],prof_len+3,double);
		NEW(DLT[i],prof_len+3,double);
                NEW(tMAT[i],prof_len+3,char);NEW(tINSRT[i],prof_len+3,char);
                NEW(tDLT[i],prof_len+3,char);
	}
	NEWP(mtch,seq_len+3,double);
	for(i=0;i<=seq_len+2;i++) NEW(mtch[i],prof_len+3,double);
	NEW(back_operation,seq_len+prof_len+2,char);
	NEW(operation,seq_len+prof_len+4,char);

 	for(i=1;i<=nAlpha(A);i++)
		for(j=1;j<=prof_len;j++)
			mtch[i][j]=(double)ValSMatrix(j,i,SM);

        tMAT[0][0]=tINSRT[0][0]=tDLT[0][0]='d';
	tMAT[1][0]=tINSRT[1][0]=tDLT[1][0]='d';
	tMAT[0][1]=tINSRT[0][1]=tDLT[0][1]='d';
	tMAT[1][1]=tINSRT[1][1]=tDLT[1][1]='d';

	MAT[0][1]=INSRT[0][1]=DLT[0][1]=od[1]; 
	MAT[1][1]=mtch[seq[1]][1]; INSRT[1][1]=io[1]; DLT[1][1]=od[1];

	Int4 res;
	for(i=2;i<=seq_len;i++){
		res=seq[i];
		MAT[i][1]=mtch[res][1]; INSRT[i][1]=INSRT[i-1][1]+ie[1]; DLT[i][1]=od[1];
		tMAT[i][1]=tINSRT[i][1]=tDLT[i][1]='d';
	}

	for(j=2;j<=prof_len;j++){
		MAT[0][j]=INSRT[0][j]=DLT[0][j]=MAT[0][j-1]+de[j];
		tMAT[0][j]=tINSRT[0][j]=tDLT[0][j]='d';
	}

// forward step...

	for(i=1;i<=seq_len;i++){
		res=seq[i];
		for(j=2;j<=prof_len;j++){
			m=MAT[i-1][j-1]+mtch[res][j];
			r=INSRT[i-1][j-1]+mtch[res][j];
			d=DLT[i-1][j-1]+mtch[res][j];
			if((m>=r)&&(m>=d)) {MAT[i][j]=m;tMAT[i][j]='m';}
			else if(r>=d) {MAT[i][j]=r;tMAT[i][j]='i';}
			     else {MAT[i][j]=d;tMAT[i][j]='d';}

			m=MAT[i-1][j]+io[j]; 
			if (i==1) r=INSRT[i-1][j]+io[j]; 
			else r=INSRT[i-1][j]+ie[j];
			if(m>=r) {INSRT[i][j]=m; tINSRT[i][j]='m';}
			else {INSRT[i][j]=r; tINSRT[i][j]='i';}

			m=MAT[i][j-1]+od[j]; d=DLT[i][j-1]+de[j];
			if (m>=d) {DLT[i][j]=m; tDLT[i][j]='m';}
			else {DLT[i][j]=d; tDLT[i][j]='d';}
		}			
	}

// trace back step..

	max=DLT[seq_len][prof_len];
	back_operation[1]='d';
	state=tDLT[seq_len][prof_len];
	score=DLT[seq_len][prof_len];
	i_max=seq_len;

	for(i=1;i<=seq_len;i++){
		if (MAT[i][prof_len]>max) 
			{max=MAT[i][prof_len]; back_operation[1]='m'; state=tMAT[i][prof_len];
			score=MAT[i][prof_len]; i_max=i-1;}
		if (DLT[i][prof_len]>max) 
			{max=DLT[i][prof_len]; back_operation[1]='d'; state=tDLT[i][prof_len];
			score=DLT[i][prof_len]; i_max=i;}
	}

	k=2; i=i_max; j=prof_len-1;	
	if(i==0) state='d';

	while(j!=0){
	  start[1]=i; start[2]=j;
	   switch(state){
		case 'i': back_operation[k]='i'; state=tINSRT[i][j]; i--; k++; break;
		case 'd': back_operation[k]='d'; state=tDLT[i][j]; j--; k++; break;	
		case 'm': back_operation[k]='m'; state=tMAT[i][j]; i--; j--; k++; break;
		default : print_error("should not happen");
	   }
	}

	oper_len=k-1;
	operation[0]='E';
	for(i=1;i<=oper_len;i++) operation[i]=back_operation[oper_len+1-i];
	operation[i]='E';

	*alignscore=score;
	*J=oper_len;

        for(i=0;i<=seq_len+2;i++){
                free(MAT[i]);free(INSRT[i]); free(DLT[i]);
                free(tMAT[i]);free(tINSRT[i]); free(tDLT[i]);
	}
	free(MAT);free(INSRT);free(DLT);
	free(tMAT);free(tINSRT);free(tDLT);
        for(i=0;i<=seq_len+2;i++) free(mtch[i]);
	free(mtch);

        free(back_operation);

	return operation;
}
