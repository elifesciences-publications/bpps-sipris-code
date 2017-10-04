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

#include "psialn.h"

psi_typ MkPsiAln(char *filename)
{
        psi_typ	A;
        Int4	i, j, k, l, seq_id, s=-1, x, poscntr,start,pos;
        Int4	al_n=1, al_l=0, seq_n, spacer,last_spacer;
        char	*str, *t,c;
	Int4	master=-1;
	BooLean	flag;
	FILE	*fp;

	fp=fopen(filename,"r");
	NEW(str,600,char);
        NEW(A,1,psialn_type);
	// count up the number of sequences...
	A->algn_length=0; flag=TRUE; last_spacer=-1;
        for(seq_n=-1; fgets(str,550,fp)!=NULL; ){
                t=str; sscanf(t,"%d",&seq_id);
		if(master < 0) master = seq_id;
		if(seq_id == master) {  // then determine align length
			spacer=j=0; t=str;
        		while(!isspace(t[0])){ t++; spacer++; } // sequence number
        		while(isspace(t[0])){ t++; spacer++; } 
        		while(!isspace(t[0])){ t++; spacer++; } // sequence position
        		while(isspace(t[0])){ t++; spacer++; } 
        		while(!isspace(t[0])){t++;j++; } // alignment
			A->algn_length += j;
			if(last_spacer >= 0 && last_spacer != spacer) 
				print_error("MkPsiAln( ) input error");
			else last_spacer=spacer;
		} 
                if(seq_id > seq_n) seq_n=seq_id;
                if(str[0]=='\n') continue;
        } seq_n++;
        A->seq_nmbr=seq_n;
	A->master = master;
	// WARNING: this can't handle more than one repeat in the align!!!
        NEW(A->align,A->seq_nmbr+1,alnseq);
        for(j=0;j < A->seq_nmbr;j++) {
		NEW(A->align[j].seq,A->algn_length+3,char);
		NEW(A->align[j].res,A->algn_length+3,Int4);
	}
        rewind(fp);
#if 0
fprintf(stderr,"seq_nmbr=%d; spacer=%d; algn_length=%d\n",
		A->seq_nmbr,spacer,A->algn_length);
#endif

        for(i=1; (fgets(str,550,fp)!=NULL); ){
                if(!isdigit(str[0])) continue;
		t=str; 
		if(sscanf(t,"%d %d", &j, &start) != 2)
		    print_error("MkPsiAln( ) input error"); // jth sequence; pos = start
                if (A->align[j].pos < 1) {
			A->align[j].start = start-1;
			A->align[j].end = start-1;
			A->align[j].pos = 1;;
		}
		t = str + spacer;
// fprintf(stderr,"%s\n",t);
		for(k=0; TRUE; k++){ // go to end of line...
			c = t[k];
			pos = A->align[j].pos;
			if(isalpha(c)){
				A->align[j].end++; // find end of 
				A->align[j].seq[pos] = c;
			} else if(isdigit(c)) {    // at end of alignment
				A->align[j].pos--; break;
			}
			A->align[j].pos++;
			A->align[j].res[pos] = A->align[j].end;
                }
        }
	fclose(fp); free(str);
        return A;
}

void    PutPsiAln(FILE *fp, psi_typ A)
{
        int	i,j,k,s,signn,pos,start,end,res;
	char	c;

    start = 1;
    for(end= MINIMUM(Int4, 60,A->algn_length); end <= A->algn_length; ){
	   s = A->master;  // first print master sequence
	   res = A->align[s].start;
	   fprintf(fp,"%-2d %5d ", s, A->align[s].res[start]);
	   for(i = start; i <= end; i++){
		c = A->align[s].seq[i];
		if(isalpha(c)){ res++; fprintf(fp,"%c",c); }
		else if(c == NULL) fprintf(fp,"-");
		else print_error("error in PutPsiAln( )");
	   } fprintf(fp," %d\n",A->align[s].res[end]);

	   for(s=0; s < A->seq_nmbr; s++){
		if(s == A->master) continue;
		else if(start <= A->align[s].pos){
		     res = A->align[s].start;
		     fprintf(fp,"%-2d %5d ", s, A->align[s].res[start]);
	   	     for(i = start; i <= end; i++){
			c = A->align[s].seq[i];
			if(isalpha(c)){ res++; fprintf(fp,"%c",c); }
			else if(c == NULL) fprintf(fp,"-");
			else print_error("error in PutPsiAln( )");
		     } fprintf(fp," %d\n",A->align[s].res[end]);
		}
	   }
	   fprintf(fp,"\n");
	   if(end < A->algn_length){
		start = end + 1;
		end += 60; 
	   	end = MINIMUM(Int4, end,A->algn_length);
	   } else break;
	}
	fprintf(fp,"\n");
}

void    NilPsiAln(psi_typ A)
{
        Int4    i,s;

        if(A != NULL) {
            if(A->align != NULL){
		for(s=0; s < A->seq_nmbr; s++){
			  free(A->align[s].seq);
			  free(A->align[s].res);
		}
                free(A->align);
            }
            free(A);
        }
}

Int4    LengMasterPsiAln(psi_typ A)
{
	return A->align[A->master].end;
}

BooLean InSeqPsiAln(Int4 s, Int4 i, psi_typ A)
// returns TRUE if position i in seq s is in align 'a'
{
	if(s < 0 || s >= A->seq_nmbr) return FALSE;
	else if(i <= A->align[s].start || i > A->align[s].end) return FALSE;
	else return TRUE;
}

Int4    MasterIDPsiAln(psi_typ A) { return A->master; }

Int4    MapSeqPsiAln(Int4 pos1, Int4 s, psi_typ A)
// returns the residue position in s that is aligned with pos1 of 
// the master sequence PSI-ALIGN A.
//  returns NULL if pos1 is not aligned with sequence s
{
	Int4	i,j,k,res,m,r;

	if(s < 0 || s >= A->seq_nmbr) return NULL;
	m = A->master;  // master sequence
	if(!InSeqPsiAln(m,pos1,A)) return NULL;
	res = A->align[m].start;
	if(res >= pos1) return NULL;
	for(i = 1; i <= A->algn_length; i++){
	   if(isalpha(A->align[m].seq[i])) res++;
	   if(res == pos1){
// fprintf(stderr,"res = %d; s = %d; i = %d\n", res,s,i);
		if(!isalpha(A->align[s].seq[i])) return NULL;
		else return A->align[s].res[i];
	   }
	}
	return NULL;
}



