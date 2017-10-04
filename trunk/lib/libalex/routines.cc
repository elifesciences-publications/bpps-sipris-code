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

#include "routines.h"

Int4 *ShuffleArray(Int4 len_a)
{
	Int4 i;
        Int4    *sh;
        dh_type H;

        NEW(sh,len_a+1,Int4);

        H = dheap(len_a+1,4);
        for(i=1; i <= len_a; i++) insrtHeap(i,((keytyp)Random()),H);
        for(i=1; i <= len_a; i++) sh[i] = delminHeap(H);

        Nildheap(H);

        return sh;
}

e_type ConsensusSeq(smx_typ M)
{
        e_type          E;
        unsigned        char *cons;
        Int4            i,k,s,p=1,max_k,len;
        double          max_val;
        a_type 		A = SMatrixA(M);

        len = LenSMatrix(M);

        NEW(cons,len+3,unsigned char);

	for(s=1;s<=LenSMatrix(M);s++){
		max_k=1; max_val = ValSMatrix(s,1,M);
		for(k=2; k<=nAlpha(A); k++){
			if(ValSMatrix(s,k,M) > max_val) {max_val=ValSMatrix(s,k,M); max_k=k;}
		}
		cons[p++] = max_k;
	}   

        E = MkSeq("cons", len, cons);
        free(cons);
    
        return E;
}

e_type ConsensusSeq(smx_typ *M, Int4 nbr, Int4 rpts)
{
        e_type          E;
        unsigned        char *cons;
        Int4            i,k,s,j,n,p=1,max_k,totlen=0;
        double          max_val;
        a_type 		A = SMatrixA(M[1]);
        for(i=1;i<=nbr;i++) totlen += LenSMatrix(M[i]);

        NEW(cons,rpts*totlen+3,unsigned char);

        for(j=1;j<=nbr;j++){
                for(s=1;s<=LenSMatrix(M[j]);s++){
                        max_k=1; max_val = ValSMatrix(s,1,M[j]);
                        for(k=2; k<=nAlpha(A); k++){
                                if(ValSMatrix(s,k,M[j]) > max_val) {max_val=ValSMatrix(s,k,M[j]); max_k=k;}
                        }
                        for(n=0;n<rpts;n++){
                                cons[n*totlen+p] = max_k;
                        }
                        p++;
                }   
        }

        E = MkSeq("cons", rpts*totlen, cons);
        free(cons);
    
        return E;
}

void PutAlign(e_type E1, smx_typ M, char *operation, Int4 operation_length, Int4 start, Int4 score)
{
// E2 is profile consensus

        Int4    i, x;
        Int4    tmp = 0, e1 = start, e2 = start-1; 
        Int4    j = 1, s = 1, k = start;
        e_type  E2 = ConsensusSeq(M);
	a_type 	A = SMatrixA(M);

        unsigned char *ptr1 = SeqPtr(E1);
        unsigned char *ptr2 = SeqPtr(E2);

        char *seq1, *seq2, *beetw;

        NEW(seq1,operation_length+2,char);
        NEW(seq2,operation_length+2,char);
        NEW(beetw,operation_length+2,char);
        for(i=1; i<=operation_length; i++){
                if (operation[i] == 'M' || operation[i] == 'm'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = AlphaChar(ptr2[j],A);
                        if(ptr1[k] == ptr2[j]) beetw[i] = ':';
                        else if(ValSMatrix(s,ptr1[k],M) > 0) beetw[i] = '.';
                        else beetw[i] = ' ';
                        j++;k++;s++;
                }
                else if (operation[i] == 'D' || operation[i] == 'd'){
                        seq1[i] = '-'; seq2[i] = AlphaChar(ptr2[j],A); beetw[i] = ' ';
                        j++;s++;
                }
                else if (operation[i] == 'I' || operation[i] == 'i'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = '-'; beetw[i] = ' ';
                        k++;
                }
                else print_error("PutAlign(): error in operation array");
        }

        printf("\n\n");

	x=1;
        for(j=1; j<=operation_length; j++){
        	if (seq2[j] == '-') x++;
                else break;
        }
        for(j=x; j<=operation_length; j++)
                printf("%c",seq2[j]);
        printf(" score=%d\n",score);
        for(j=x; j<=operation_length; j++)
                printf("%c",beetw[j]);
        printf("\n");
        for(tmp=0,j=x; j<=operation_length; j++){
                printf("%c",seq1[j]);
                if (seq1[j] != '-') tmp++;
        }
        e2+=(x+tmp-1); e1=e2-tmp+1;
        printf(" %d-%d",e1,e2);
        printf("\n\n\n");

        free(seq1); free(seq2); free(beetw);
	NilSeq(E2);
}

void PutAlign(e_type E1, smx_typ *M, Int4 nbl, char *operation, Int4 operation_length, 
                        Int4 start, Int4 rpts)
{
        Int4    i, x;
        Int4    tmp = 0, e1 = start, e2 = start-1; 
        Int4    j = 1, s = 1, t = 1, k = start;
        e_type  E2 = ConsensusSeq(M, nbl, rpts);
        Int4    rpt_len = LenSeq(E2)/rpts;
        Int4    *lengths;
	a_type 	A = SMatrixA(M[1]);

        NEW(lengths,rpts + 1,Int4);

        unsigned char *ptr1 = SeqPtr(E1);
        unsigned char *ptr2 = SeqPtr(E2);

        char *seq1, *seq2, *beetw;

        NEW(seq1,operation_length+2,char);
        NEW(seq2,operation_length+2,char);
        NEW(beetw,operation_length+2,char);
        for(i=1; i<=operation_length; i++){
                if (operation[i] == 'M' || operation[i] == 'm'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = AlphaChar(ptr2[j],A);
                        if(ptr1[k] == ptr2[j]) beetw[i] = ':';
                        else if(ValSMatrix(s,ptr1[k],M[t]) > 0) beetw[i] = '.';
                        else beetw[i] = ' ';
                        j++;k++;
                        if(s>=LenSMatrix(M[t])) {s=1;t++;}
                        else s++;
                }
                else if (operation[i] == 'D' || operation[i] == 'd'){
                        seq1[i] = '-'; seq2[i] = AlphaChar(ptr2[j],A); beetw[i] = ' ';
                        j++;
                        if(s>=LenSMatrix(M[t])) {s=1;t++;}
                        else s++;
                }
                else if (operation[i] == 'I' || operation[i] == 'i'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = '-'; beetw[i] = ' ';
                        k++;
                }
                else print_error("PutAlign(): error in operation array");
        }

        printf("\n\n");

        Int4 br = 0; j=1; i=1;
        while (i <= operation_length){
                while(br < rpt_len){
                        if (seq2[i++] != '-') br++;
                        else lengths[j] += 1;
                }
                br = 0; j++;
        }

        Int4 total=0;
        for(i=1; i<=rpts; i++){
                x=1;
                for(j=1; j<=rpt_len + lengths[i]; j++){
                        if (seq2[total+j] == '-') x++;
                        else break;
                }
                for(j=x; j<=rpt_len + lengths[i]; j++)
                        printf("%c",seq2[total+j]);
                printf("\n");
                for(j=x; j<=rpt_len + lengths[i]; j++)
                        printf("%c",beetw[total+j]);
                printf("\n");
                for(tmp=0,j=x; j<=rpt_len + lengths[i]; j++){
                        printf("%c",seq1[total+j]);
                        if (seq1[total+j]!='-') tmp++;
                }
                e2+=(x+tmp-1);e1=e2-tmp+1;
                printf(" %d-%d",e1,e2);
                printf("\n\n\n");
                total += rpt_len + lengths[i];
        }

        free(seq1); free(seq2); free(beetw); free(lengths);
	NilSeq(E2);
}
#if 0
char *ReverseOperArray(char *operation1, Int4 oper_len1)
{
        Int4    i, k=1;
        char    *operation;
                           
        NEW(operation,oper_len1+1,char);
        
        operation[0] = operation[oper_len1+1] = 'E';
        
        for(i=oper_len1;i>=1;i--){
                switch(operation1[i]){
                        case 'M': if(operation1[i+1]=='i' || operation1[i+1]=='E') 
					operation[k++]='M';
                                  else operation[k++] = 'm';
                                  break;
                        case 'm': if(operation1[i+1]=='i' || operation1[i+1]=='E') 
					operation[k++]='M';
                                  else operation[k++] = 'm';
                                  break;
                        case 'i': operation[k++] = 'i';
                                  break;
                        case 'I': operation[k++] = 'I';
                                  break;
                        case 'd': if(operation1[i+1]=='i' || operation1[i+1]=='E') 
					operation[k++]='D';
                                  else operation[k++] = 'd';
                                  break;
                        case 'D': if(operation1[i+1]=='i' || operation1[i+1]=='E') 
					operation[k++]='D';
                                  else operation[k++] = 'd';
                                  break;
                        default : print_error("ReverseOperArray(): error in the input array");
                }
        }

        return operation;
}
#endif

char *ReverseOperArray(char *operation1, Int4 oper_len1)
{
        Int4 i;
        char *operation;
        
        NEW(operation,oper_len1+2,char);
	operation[0] = 'E'; operation[oper_len1+1]= 'E';
        
        for(i=1; i<=oper_len1; i++){
		operation[i] = operation1[oper_len1-i+1];
        }
	if (operation[1]=='m') operation[1]='M';
	else if(operation[1]=='d') operation[1]='D';
	else print_error("ReverseOperArray(): error in operation array");

        if (operation[oper_len1]=='M') operation[oper_len1]='m';
        else if(operation[oper_len1]=='D') operation[oper_len1]='d';
        else print_error("ReverseOperArray(): error in operation array");

        return operation;  
}

Int4 *ReversePen(Int4 *gap, Int4 len)
{
	Int4 i;
	Int4 *ngap;

	NEW(ngap,len+2,Int4);
	
	for(i=1; i<=len; i++) ngap[i] = gap[len-i+1];

	return ngap;
}

char *ConcatenateOperArrays(char **opers, Int4 n_opers, Int4 *oper_lens, 
	Int4 *starts, Int4 *ends, Int4 seq_len, Int4 *operation_length)
{
	char *operation;
	Int4 i,j,k,tot_len=0;

	for(i=1;i<=n_opers;i++)
		tot_len+=oper_lens[i];
	NEW(operation,tot_len+seq_len,char);

	operation[0]='E';
	i=1;
	for(j=1;j<=starts[1]-1;j++) operation[i++] = 'i';
	for(j=1;j<=oper_lens[1];j++) operation[i++]=opers[1][j]; 

	for(k=2;k<=n_opers;k++){
		for(j=ends[k-1]+1;j<=starts[k]-1;j++){
			operation[i++]='i';
		}
		for(j=1;j<=oper_lens[k];j++){
			operation[i++]=opers[k][j];
		}
	}

	for(j=ends[n_opers]+1;j<=seq_len;j++) operation[i++]='i';
	operation[i]='E';
	*operation_length = i-1;

	return operation;
}

Int4 CompareOperArr(char *oper1, char *oper2, Int4 *pos, 
                Int4 num_rpts1, Int4 num_rpts2, Int4 rpt_len)
{
        Int4    i,j,k,br,score=0,max,j_best,totscore=0;
        Int4    **p,**q;

        NEWP(p,num_rpts1+2,Int4);
        NEWP(q,num_rpts2+2,Int4);
        for(i=1;i<=num_rpts1+1;i++){
                NEW(p[i],rpt_len+1,Int4);
        }
        for(i=1;i<=num_rpts2+1;i++){
                NEW(q[i],rpt_len+1,Int4);
        }
        i=1;j=1;k=1;br=0;
        while(oper1[k]!='E'){
                switch(oper1[k]){
                        case 'i': br++; k++; break; 
                        case 'I': br++; k++; break;
                        case 'D': p[i][j]=-br; 
                                  if (j==rpt_len) {i++;j=1;} 
                                  else j++;
                                  k++; break;
                        case 'd': p[i][j]=-br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        case 'm': br++; p[i][j]=br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        case 'M': br++; p[i][j]=br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        default : print_error("CompareOperArray(): error in operation array1");
                }
        }
        i=1;j=1;k=1;br=0;
        while(oper2[k]!='E'){
                switch(oper2[k]){
                        case 'i': br++; k++; break;
                        case 'I': br++; k++; break;
                        case 'D': q[i][j]=-br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        case 'd': q[i][j]=-br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        case 'm': br++; q[i][j]=br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        case 'M': br++; q[i][j]=br; 
                                  if (j==rpt_len) {i++;j=1;}
                                  else j++;
                                  k++; break;
                        default : print_error("CompareOperArray(): error in operation array2");
                }
        }
        j=1; j_best=0;
                for(i=1;i<=num_rpts1;i++){
                        max=0;
                        while(j <= num_rpts2){
                                score = CompareRpts(p[i],q[j],rpt_len, pos);
                                if(score > max) {max = score; j_best = j;}
                                j++;
                        }
                        j = j_best+1;
                        totscore += max;
                }
        for(i=0;i<=num_rpts1+1;i++) free(p[i]); free(p);
        for(i=0;i<=num_rpts2+1;i++) free(q[i]); free(q);

        return totscore;
}

Int4 CompareRpts(Int4 *p, Int4 *q, Int4 rpt_len, Int4 *pos)
{
        Int4 score=0, k=1;

        for(Int4 i=1;i<=rpt_len;i++){
                if(i==pos[k]){
                        k++;
                        if(p[i]==q[i]) score+=1;
                }
        }


        return score;
}

char *ExtendOperArray(char *operation, Int4 operation_length, Int4 start)
{
        Int4    i;
        char    *oper;
        NEW(oper,operation_length+start+3,char);
        oper[0]='E';
        for(i=1;i<start;i++) oper[i]='i';
        for(i=0;i<operation_length;i++) oper[start+i]=operation[i+1];
        oper[start+operation_length]='E';
        return oper;
}

