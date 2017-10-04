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

#include "n2a.h"

Int4 ***create_table_n2a(Int4 code)
{
	const Int4 table1[5][5][5]=
		{
		{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
		{{0,0,0,0,0},{0,10,6,10,6},{5,5,5,5,5},{0,11,4,11,4},{0,17,17,19,17}},//A
		{{0,0,0,0,0},{0,9,12,9,12},{20,20,20,20,20},{11,11,11,11,11},{18,18,18,18,18}},//C
		{{0,0,0,0,0},{0,8,7,8,7},{3,3,3,3,3},{2,2,2,2,2},{16,16,16,16,16}},//G
		{{0,0,0,0,0},{0,21,14,21,14},{4,4,4,4,4},{0,21,1,13,1},{0,18,15,18,15}}//U
		};

        Int4 ***table,i,j,k;
        NEWPP(table,5,Int4);
        for(i=0;i<=4;i++) {
                NEWP(table[i],5,Int4);
                for(j=0;j<=4;j++)
                        NEW(table[i][j],5,Int4);
        }

        for(i=0;i<=4;i++) {
                for(j=0;j<=4;j++){
                        for(k=0;k<=4;k++){
                                table[i][j][k]=table1[i][j][k];
                        }
                }
        }
        switch (code){
                case  1:   break;
                case  2: table[1][3][1]=0;table[1][3][3]=0;table[1][4][1]=19;table[4][3][1]=13;
                           break;
                case  3: table[1][4][1]=19;table[2][4][4]=5;table[2][4][2]=5;table[2][4][1]=5;
                           table[2][4][3]=5;table[2][4][0]=5;table[4][3][1]=13;table[2][3][1]=0;
                           table[2][3][2]=0; 
                           break;
                case  4: table[4][3][1]=13;
                           break;
                case  5: table[1][3][1]=4;table[1][3][3]=4;table[1][3][0]=4;table[1][4][1]=19;
                           table[4][3][1]=13;
                           break;
                case  6: table[4][1][1]=9;table[4][1][3]=9;
                           break;
                case  9: table[1][1][1]=6;table[1][3][1]=4;table[1][3][3]=4;table[1][3][0]=4;
                           table[4][3][1]=13;
                           break;
                case 10: table[4][3][1]=1;
                           break;
                case 11:   break;
                case 12: table[2][4][3]=4;
                           break;
                case 13: table[1][3][1]=2;table[1][3][3]=2;table[1][4][1]=19;table[4][3][1]=13;
                           break;
                case 14: table[1][1][1]=6;table[1][3][1]=4;table[1][3][3]=4;table[1][3][0]=4;
                           table[4][1][1]=14;table[4][3][1]=13;
                           break;
                case 15: table[4][1][3]=9;
                           break;
                case 16:   break;
                case 21:   break;
                default: print_error("no such conversion table");
        }
        return table;
}

n2a_typ	MkN2A(char *filename, Int4 code, Int4 max_seq, Int4 min_seq, 
		Int4 max_len, Int4 min_len, a_type A, a_type dnaA)
{
	n2a_typ	n2a;
	const char aa[]="XCGASTNDEQKRHWYFVILMP";
	const char nt[]="NACGT";

	if(nAlpha(A) != 20 || nAlpha(dnaA) != 4) print_error("MkN2A( ): inconsistent alphabet");
	else {
		for(Int4 i=0;i<=20;i++){
			if(AlphaChar(i,A)!=aa[i]) print_error("MkN2A( ): inconsistent alphabet");
		}
		for(Int4 j=0;j<=4;j++){
			if(AlphaChar(j,dnaA)!=nt[j]) print_error("MkN2A( ): inconsistent alphabet");
		}
	}	
	NEW(n2a,1,dna2aa_type);
	for(Int4 rf=1; rf <=NumRFsN2A(n2a); rf++){
		n2a->num[rf]=0;
		n2a->E[rf]=NULL;
	}
	n2a->table=create_table_n2a(code);
	n2a->dnaE=NULL;
	n2a->A = A; 
	n2a->J = 0; 
	n2a->dnaA = dnaA; 
	n2a->fptr = open_file(filename,"","r");
	n2a->filename = AllocString(filename);
	n2a->max_seq=max_seq;
	n2a->min_seq=min_seq;
	n2a->max_len=max_len;
	n2a->min_len=min_len;
	Int4 *counts, *nsize;
	n2a->number = GetLongFastaInfo(filename, max_seq, &counts, &nsize, dnaA);
	n2a->counts=counts;
	n2a->nsize=nsize;
	return n2a;
}

void reset_n2a(n2a_typ n2a)
{
    if(n2a->dnaE != NULL){
	for(Int4 rf=1; rf <=NumRFsN2A(n2a); rf++){
	     for(Int4 n = 1; n <= n2a->num[rf]; n++){
		NilSeq(n2a->E[rf][n]);
	     }
	     free(n2a->E[rf]);
	}
	NilSeq(n2a->dnaE);
    }
    n2a->dnaE=NULL;   
}

void NilN2A(n2a_typ n2a)
{
	Int4	i,j;
	reset_n2a(n2a);
	for(i=0;i<=4;i++){
		for(j=0;j<=4;j++) free(n2a->table[i][j]); 
		free(n2a->table[i]);
	} free(n2a->table);
	free(n2a->filename);
	free(n2a->counts);
	free(n2a->nsize);
	fclose(n2a->fptr);
	free(n2a);
}

void reverse_compl_n2a(n2a_typ n2a)
{
        Int4    i, j;
	char a;

        ReverseSeq(n2a->dnaE);
	unsigned char *ptr = SeqPtr(n2a->dnaE);

        for(i=1;i<=LenSeq(n2a->dnaE);i++){
		a = AlphaChar(ptr[i],n2a->dnaA);
                if (a == 'A')  ptr[i] = AlphaCode('T',n2a->dnaA);
                else if (a == 'T') ptr[i] = AlphaCode('A',n2a->dnaA);
		else if (a == 'G') ptr[i] = AlphaCode('C',n2a->dnaA);
		else if (a == 'C') ptr[i] = AlphaCode('G',n2a->dnaA);		
        }
}

e_type SeqN2A(Int4 rf, Int4 i, n2a_typ n2a)
{
	if (rf > NumRFsN2A(n2a) || rf < 1) 
		print_error("SeqN2A( ): wrong number of reading frames");
	else if(i > n2a->num[rf] || i <1) 
		print_error("SeqN2A( ): wrong number of seqs in a reading frame");
	else
		return n2a->E[rf][i];
}

BooLean	ReadSeqN2A(n2a_typ n2a, BooLean rev)
{
	n2a->J++; 
	reset_n2a(n2a);
	if(n2a->J > n2a->number) return NULL;

	n2a->dnaE = ReadSeq(n2a->fptr,n2a->J,n2a->nsize[n2a->J],n2a->dnaA);
	if(n2a->dnaE == NULL) return NULL;
	else {
	   for(Int4 rf=1;rf<=NumRFsN2A(n2a);rf++) 
		NEW(n2a->E[rf],LenSeq(n2a->dnaE)/n2a->min_len+5,e_type);
	   translate_n2a(FALSE, n2a);
	   if (rev) translate_n2a(TRUE, n2a);
	   return TRUE;
	}
}

BooLean	translate_n2a(BooLean rev, n2a_typ n2a)
{
        Int4            start=1, end, rf, j, k, m, n, r, t, x, y, lng, lngth, len[4];
        unsigned char   a, *seq, *P;
        char            *st, str[500], temp1[50], *temp2, *p, c;

        NEW(st,500,char);
        StrSeqInfo(st,n2a->dnaE);
        for(j=0;st[j]!=' ';j++) temp1[j]=st[j]; temp1[j]='\0';
        temp2=st; temp2+=j;

        if (rev) {reverse_compl_n2a(n2a); c='-'; x=4; y=6;}
        else {c='+'; x=1; y=3;}

        lngth = LenSeq(n2a->dnaE);

        NEW(P,lngth/3+5,unsigned char);
        seq = SeqPtr(n2a->dnaE);
        r = lngth%3;
        if (r == 0) {len[1] = lngth/3; len[2] = lngth/3-1; len[3] = lngth/3-1;}
        else if (r == 1) {len[1] = lngth/3; len[2] = lngth/3; len[3] = lngth/3-1;}
        else {len[1] = lngth/3; len[2] = lngth/3; len[3] = lngth/3;}

        for(j=1,rf=x; rf<=y; j++, rf++){
                n=1;k=1;t=1;
                while(k<=len[j]){
                        m=3*(k-1)+j;
                        P[t]=n2a->table[seq[m]][seq[m+1]][seq[m+2]];
                        if(P[t]==21 || k==len[j]){
                                if(t>n2a->min_len){
                                        end=m-1;start=end-(3*t)+4;
                                        sprintf(str,"%s (nt %c%d:%d)%s",temp1,c,j,start,temp2);
					if(P[t]!=21 && k==len[j]) lng=t;
					else lng=t-1;
                                        n2a->E[rf][n]=MkSeq(str,lng,P);
                                        n++;
                                }
                                t=1;k++;
                        }
                        else {k++;t++;}
                }
                n2a->num[rf]=n-1;
        }
        free(st);
        free(P);
        return 1;
}
