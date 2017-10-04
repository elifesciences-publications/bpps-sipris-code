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

#include "s-w.h"

char *sw_loc_algn_smatrix(Int4 seq_len, unsigned char *seq, smx_typ SM, Int4 *io, 
Int4 *ie,Int4 *od, Int4 *de, Int4 *alignscore,Int4 *J, Int4 *start)
{
	int i, j, k, signal=0, oper_len;
	char *operation, *back_operation;
	char Tstate;
	double **MAT, **INSRT, **DLT;
	double m, r, d;
	double **mtch;
	int prof_len = LenSMatrix(SM);
	a_type A = SMatrixA(SM);
	double score=0, max=0;
	int max_i,max_j;

	NEWP(MAT,seq_len+3,double);NEWP(INSRT,seq_len+3,double);NEWP(DLT,seq_len+3,double);
	for(i=0;i<=seq_len+2;i++){
		NEW(MAT[i],prof_len+3,double);NEW(INSRT[i],prof_len+3,double);
		NEW(DLT[i],prof_len+3,double);
	}

	NEWP(mtch,seq_len+3,double);
	for(i=0;i<=seq_len+2;i++) NEW(mtch[i],prof_len+3,double);
	NEW(back_operation,seq_len+prof_len+2,char);
	NEW(operation,seq_len+prof_len+4,char);

 	for(i=1;i<=seq_len;i++)
		for(j=1;j<=prof_len;j++)
			mtch[i][j]=(double)ValSMatrix(j,seq[i],SM);

// forward step...

	for(i=1;i<=seq_len;i++){
		for(j=1;j<=prof_len;j++){
			m=MAT[i-1][j-1]+mtch[i][j];
			r=INSRT[i-1][j-1]+mtch[i][j];
			d=DLT[i-1][j-1]+mtch[i][j];
			if((m>=r)&&(m>=d)&&(m>0)) MAT[i][j]=m;
			else if((r>=d)&&(r>0)) MAT[i][j]=r;
			     else if (d>0) MAT[i][j]=d;
				  else MAT[i][j]=0;

			if (MAT[i][j] > max){ max=MAT[i][j]; max_i=i; max_j=j;}

			m=MAT[i-1][j]+io[j]; r=INSRT[i-1][j]+ie[j];
			if((m>=r)&&(m>0)) INSRT[i][j]=m;
			else if (r>0) INSRT[i][j]=r;
			     else INSRT[i][j]=0;

			m=MAT[i][j-1]+od[j]; d=DLT[i][j-1]+de[j];
			if ((m>=d)&&(m>0)) DLT[i][j]=m;
			else if (d>0) DLT[i][j]=d;
			     else DLT[i][j]=0;

						
		}			
	}


// trace back step..

	i=max_i;j=max_j;k=1;
	score=MAT[i][j];
	Tstate='m';

	while(signal==0){
	   switch(Tstate){
		case 'i': back_operation[k]='i'; i--; k++; 
			  if ((MAT[i][j]+io[j]>=INSRT[i][j]+ie[j])&&(MAT[i][j]>0)) Tstate='m';
			  else if (INSRT[i][j]>0) Tstate='i';
			  else Tstate='s';
			  break;
		case 'd': back_operation[k]='d'; j--; k++; 
			  if ((MAT[i][j]+od[j+1]>=DLT[i][j]+de[j+1])&&(MAT[i][j]>0)) Tstate='m';
                          else if (DLT[i][j]>0) Tstate='d';
                          else Tstate='s';
			  break;	
		case 'm': back_operation[k]='m'; i--; j--; k++;
			  if ((MAT[i][j]>=DLT[i][j])&&(MAT[i][j]>=INSRT[i][j])&&(MAT[i][j]>0)) Tstate='m';
                          else if ((INSRT[i][j]>=DLT[i][j])&&(INSRT[i][j]>0)) Tstate='i';
			  else if (DLT[i][j]>0) Tstate='d';
                          else Tstate='s';
			  break;
		case 's': oper_len=k-1; signal=1; start[1]=i+1; start[2]=j+1;
			  break;
		default : print_error("this should not happen");
	   }
	}	
printf("start=(%d,%d) end =(%d,%d)\n", start[1], start[2], max_i, max_j);
	operation[0]='E';
	for(i=1;i<=oper_len;i++) operation[i]=back_operation[oper_len+1-i];
	operation[i]='E';

	*alignscore=score;
	*J=oper_len;

        for(i=0;i<=seq_len+2;i++){
                free(MAT[i]);free(INSRT[i]); free(DLT[i]);
	}
	free(MAT);free(INSRT);free(DLT);
        for(i=0;i<=seq_len+2;i++) free(mtch[i]);
	free(mtch);

        free(back_operation);

	return operation;
}
