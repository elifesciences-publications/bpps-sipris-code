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

#include "gpsihmm.h"

gph_typ MakeGPsiHMM(Int4 hpsz,e_type E,a_type A,Int4 proflen)
{
	gph_typ B;
	Int4	i,z;
	Int4 len = LenSeq(E);
	NEW(B,1,gpsihmm_type);
	B->A = A; B->E = E; B->mH = NULL;
	B->diag0 = B->ed0 = NULL;
	B->hit_s = B->hit_e = B->hit_d = NULL;
	B->hit_s2 = B->hit_e2 = B->hit_d2 = NULL; B->score = NULL;
        MEW(B->diag0,len+proflen+5,Int4);
        B->diag = B->diag0 + len;
        NEW(B->ed0,len+proflen+5,Int4);
        B->extdiag = B->ed0 + len;
        z = B->zero = - len;  
        for(i = z; i <= proflen; i++) B->diag[i] = z;
	B->mH = Mheap(hpsz,3);
        NEW(B->hit_s,hpsz+5,Int4); NEW(B->hit_s2,hpsz+5,Int4);   
        NEW(B->hit_e,hpsz+5,Int4); NEW(B->hit_e2,hpsz+5,Int4);
        NEW(B->hit_d,hpsz+5,Int4); NEW(B->hit_d2,hpsz+5,Int4);
        B->nhits = 0;
        NEW(B->score,hpsz+3,double);
	return B;
}

void	NilGPsiHMM(gph_typ B)
{
	if(B->diag0 != NULL) free(B->diag0); 
	if(B->ed0 != NULL) free(B->ed0);
	if(B->mH != NULL) NilMheap(B->mH);
	if(B->hit_s != NULL) free(B->hit_s);
	if(B->hit_e != NULL) free(B->hit_e);
	if(B->hit_d != NULL) free(B->hit_d);
	if(B->hit_s2 != NULL) free(B->hit_s2);
	if(B->hit_e2 != NULL) free(B->hit_e2);
	if(B->hit_d2 != NULL) free(B->hit_d2);
	if(B->score!= NULL) free(B->score);
	free(B);
}

void	UpdateHitsGPsiHMM(gph_typ B)
{
      Int4	n,i,r;
      n = ItemsInMheap(B->mH); 
      for(r=1; ItemsInMheap(B->mH) > 0; r++){
	  B->score[r] = -(double) MinKeyMheap(B->mH); 
	  i = DelMinMheap(B->mH); 
	  B->hit_s2[r] = B->hit_s[i]; 
	  B->hit_e2[r] = B->hit_e[i]; 
	  B->hit_d2[r] = B->hit_d[i]; 
      } 
      B->update = FALSE;
      B->nhits=n;
}

Int4    ExtendGPsiHMM(Int4 lenMaster, Int4 i1, Int4 len2, unsigned char *seq2,
        Int4 i2, Int4 *left, Int4 *right, Int4 **matrix, register Int4 x_dropoff)
/****************************************************************
                    start (*left=i1+1 *right=i2+1)
                     <|------> (max = 40)
                word=SPH                word = 14
                    -278474154
        Query:   481 SPHSPLVRA 489
                      PHSPL+RA
        Sbjct:   488 IPHSPLLRA 496      done: *left=i1+1; *right=i2+7.
NOTE: there may be a slight speed up by setting x_dropoff to a constant (e.g. 22).
 ****************************************************************/
//left and right are positions in the subject sequence (which starts from 1) 
//where hsp starts and ends
{
        register Int4   i,s,max,end,**mtx;
        register unsigned char  *p2;
        
        mtx=matrix+i1; p2=seq2+i2;
        if((max=len2-i2) > (s=lenMaster-i1)) end=s; else end=max;
        *right = 1;
        s = mtx[1][p2[1]];
        for(max=s, i=2; i <= end; i++) {
                if((s += mtx[i][p2[i]]) > max) { *right=i; max=s; }
                else if((s <= (max-x_dropoff)) || s < 0) break;
        } *right += i2; *left = 1;
        if(i1<=i2) end=-i1; else end=-i2;
        for(s=max, i=0; i > end; i--){
                if((s += mtx[i][p2[i]]) > max) { *left=i; max=s; }
                else if((s <= (max-x_dropoff)) || s < 0) break;
        } *left += i2;
        return max;
}

Int4	ExtendGPsiHMMStr(Int4 lenMaster, Int4 i1, Int4 len2, unsigned char *seq2,
	Int4 i2, Int4 *left, Int4 *right, Int4 **matrix, register Int4 x_dropoff, Int4 *m2m)
/****************************************************************
                    start (*left=i1+1 *right=i2+1)
		     <|------> (max = 40)
		word=SPH		word = 14
	            -278474154	
	Query:   481 SPHSPLVRA 489	
	              PHSPL+RA
	Sbjct:   488 IPHSPLLRA 496	done: *left=i1+1; *right=i2+7.
NOTE: there may be a slight speed up by setting x_dropoff to a constant (e.g. 22).
 ****************************************************************/
{
	register Int4	i,s,max,end,**mtx,*tr;
	register unsigned char	*p2;

	mtx=matrix+i1; p2=seq2+i2; tr = m2m + i1;
	if((max=len2-i2) > (s=lenMaster-i1)) end=s; else end=max; 
	*right = 1;
	s = mtx[1][p2[1]];
	for(max=s, i=2; i <= end; i++) {
		if((s += (mtx[i][p2[i]] + tr[i-1])) > max) { *right=i; max=s; }
		else if((s <= (max-x_dropoff)) || s < 0) break; 
	} *right += i2; *left = 1;
	if(i1<=i2) end=-i1; else end=-i2;
	for(s=max, i=0; i > end; i--){
		if((s += (mtx[i][p2[i]] + tr[i])) > max) { *left=i; max=s; }
		else if((s <= (max-x_dropoff)) || s < 0) break; 
	} *left += i2;
	return max;
}

Int4  MatcherGPsiBlstStr(hit_typ head,Int4 query_length,Int4 uthreshold,
	Int4 ugpxdrop,Int4 **mtx,Int4 subj_len,register unsigned char *subj_ptr,
	register gph_typ B,Int4 *nwordhts)
{
	hit_typ tmp;
	tmp = head;
        register Int4   i,d,*ed;
        Int4            n,*next,nxt,score,s,e,max=0,item;
        while(DelMinMheap(B->mH) != NULL) ; B->update = TRUE;
        nxt=B->zero; next=B->diag; ed=B->extdiag;
	*nwordhts=0;
	while (tmp != 0) {
		*nwordhts+=1;
		i = tmp->pos_seq;
                d = tmp->pos_hmm - i;
                if(i > ed[d] + 1){
                  // i should point to the middle of the word hit A.P.
                  // in ExtendGPsiHMMStr both matrix and subj_ptr start from 1 A.P.
		  score=ExtendGPsiHMM(query_length,tmp->pos_hmm,subj_len,subj_ptr,
		     i,&s,&e,mtx, ugpxdrop);
                  if(score > max) max = score;
                  if(ed[d] == 0){ next[nxt] = d; nxt = d; }
                  ed[d] = e;
                  if(score >= uthreshold && (item=InsertMheap(-score,B->mH))!=NULL){
                        B->hit_s[item]=s; B->hit_e[item]=e; B->hit_d[item]=d;
                  }
                }
		tmp = tmp->next;
        }
        for(d=B->zero; next[d]!=B->zero; d=n){ // free extend list.
                n=next[d]; ed[n]=0; next[d]=B->zero;
        }
        return max;
}
