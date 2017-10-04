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

#include "gpsiblst.h"

static Int4 count_hits_gpsiblst(Int4 len, Int4 **matrix, Int4 T, a_type A)
// Calculate the maximum number of words that can possibly score >= T
// at an arbitrary sequence site. (Needed for allocating list.) 
{
	Int4		r0,r1,s,hits,i0,i1,i2,*m0,*m1;
	register Int4	r2,t,*m2;

	for(hits=0,i0=1,i1=2,i2=3; i2 <= len; i0++,i1++,i2++){
// fprintf(stderr,"pos(%d): ",i0);
	   m0=matrix[i0]; m1=matrix[i1]; m2=matrix[i2];
	   for(r0=0; r0 <= nAlpha(A); r0++){
		s = T - m0[r0];
		for(r1=0; r1 <= nAlpha(A); r1++){	
		   t = s - m1[r1];
		   for(r2=0; r2 <= nAlpha(A); r2++){	
			if(t <= m2[r2]){
				hits++;
#if 0 	// add to acceptor states.
fprintf(stderr,"%c%c%c(%d) ",AlphaChar(r0,A), AlphaChar(r1,A),
	AlphaChar(r2,A), matrix[i0][r0]+matrix[i1][r1]+matrix[i2][r2]);
#endif
			} 
		   }
		}
	   } 
// fprintf(stderr,"\n");
	} return hits;
}

gpb_typ MakeGPsiBlst(Int4 hpsz,Int4 cutoff,Int4 T,e_type E,a_type A,
	Int4 **matrix,Int4 x_dropoff)
{
	gpb_typ B;
	Int4	n,q,q0,q1,q2,c,d,e,x,y,z,hits,nlet,i;

	if(T >= 40) print_error("T too large; no hits possible.");
	Int4 len = LenSeq(E);
	NEW(B,1,gpsiblast_type);
	B->matrix=matrix;
	B->A = A; B->E = E; B->mH = NULL;
	B->diag0 = B->ed0 = NULL;
	B->hit_s = B->hit_e = B->hit_d = NULL;
	B->hit_s2 = B->hit_e2 = B->hit_d2 = NULL; B->score = NULL;
	B->x_dropoff = x_dropoff;
	NEW(B->tmp,len+2,Int4); 

	// Create Finite Automaton.
	// 1. Create transition function: d(r,q).
	n = nAlpha(A) +1; n = 1 + n + n*n;  
	NEWP(B->d,nAlpha(A)+2,Int4);
	for(c=0; c <= nAlpha(A); c++) NEW(B->d[c],n+2,Int4);
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q++; B->d[c][0] = q; q0=q;		// q0 = state "^c"
	    for(d=0; d <= nAlpha(A); d++){
	    	q++; B->d[d][q0] = q;        	// q1 = state "cd"
	    }
	}
	B->nQ = q;				// nQ = #states 
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q0 = B->d[c][0];			// q = state "^c"
	    for(d=0; d <= nAlpha(A); d++){
	    	q1 = B->d[d][q0];		// q = state "cd" 
	    	q2 = B->d[d][0];		// q2 = state "^d" 
		for(e=0; e <= nAlpha(A); e++){
			q = B->d[e][q2];	// q2 = state "de" 
			B->d[e][q1] = q;	// "cde" = "^de" = "de" 
		}
	    }
	}
	// 2. create acceptance states & lists[q][r][1..]
	hits=count_hits_gpsiblst(LenSeq(E), matrix, T, A);
// std::cerr << hits; std::cerr << std::endl;

	nlet = nAlpha(A) +1; // n = B->nQ; // ????
	B->mlist=MkMList(hits+5,MaxStateGPB(nlet));

	Int4			r0,r1,s,i0,i1,i2,*m0,*m1;
	register ml_type	mlist = B->mlist;
	register Int4 		r2,t,*m2;
	for(i0=1,i1=2,i2=3; i2 <= len; i0++,i1++,i2++){
	   // look through 3-words hits for "matches".
	   m0=matrix[i0]; m1=matrix[i1]; m2=matrix[i2];
	   for(r0=0; r0 <= nAlpha(A); r0++){
	        q0 = B->d[r0][0];
		s = T - m0[r0];
		for(r1=0; r1 <= nAlpha(A); r1++){	
	           q1 = B->d[r1][q0];
		   t = s - m1[r1];
		   for(r2=0; r2 <= nAlpha(A); r2++){	
			if(t <= m2[r2]) // then add to acceptor states.
				{ Add2MList(i0,StateGPB(q1,r2), mlist); }
		   }
		}
	   } 
	}
        MEW(B->diag0,len+MAX_SEQ_LENG_GPB+5,Int4);
        B->diag = B->diag0 + MAX_SEQ_LENG_GPB;
        NEW(B->ed0,len+MAX_SEQ_LENG_GPB+5,Int4);
        B->extdiag = B->ed0 + MAX_SEQ_LENG_GPB;
        z = B->zero = - MAX_SEQ_LENG_GPB;
        for(i = z; i <= len; i++) B->diag[i] = z;
        B->cutoff = cutoff; B->mH = Mheap(hpsz,3);
        NEW(B->hit_s,hpsz+5,Int4); NEW(B->hit_s2,hpsz+5,Int4);
	NEW(B->hit_e,hpsz+5,Int4); NEW(B->hit_e2,hpsz+5,Int4);
        NEW(B->hit_d,hpsz+5,Int4); NEW(B->hit_d2,hpsz+5,Int4);
        B->nhits = 0;
        NEW(B->score,hpsz+3,double);
// fprintf(stderr,"total hits = %d; list size = %d\n",hits,SizeMList(B->mlist));
	return B;
}

void	NilGPsiBlst(gpb_typ B)
{
	for(Int4 c=0; c <= nAlpha(B->A); c++) { free(B->d[c]); }
	NilMList(B->mlist);
	if(B->diag0 != NULL) free(B->diag0); 
	if(B->ed0 != NULL) free(B->ed0);
	free(B->tmp); free(B->d); 
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

Int4	MatcherGPsiBlst(FILE *fptr, e_type E, gpb_typ B)
{
	Int4	num,n,q,score,s,c,d,e,i,j,hits=0,max,*tmp=B->tmp,maxoff;
	a_type	A=B->A;
	Int4        lft,rt;

	assert(MAX_SEQ_LENG_GPB > LenSeq(E));
	n = LenSeq(E) - 1; c = XSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(max=0, i=2; i <= n; c=d, i++){
	   d = XSeq(i,E);
	   q = B->d[d][q];		// q = state "xd" 
	   e = XSeq(i+1,E);		// e = next token (mealy model)
					// then q + e signals acceptance.
	   if(!EmptyMList(StateGPB(q,e),B->mlist)){
		num=GetListMList(tmp, StateGPB(q,e), B->mlist);
		for(j=0;  j < num; j++){
			s = tmp[j];
			score=ExtendGPsiBlstStr(LenSeq(B->E),s,LenSeq(E), XSeqPtr(E),
					i-1, &lft, &rt, B->matrix,B->x_dropoff);
			if(score > max){ max = score; maxoff = s-i+1; }
		} hits+=j;
	   }
	}
	if(max > 0 && fptr != NULL) {
		fprintf(fptr,"hits = %d; max score = %d\n",hits,max); 
		fprintf(fptr,"score = %d; offset = %d\n", max,maxoff);
		PutDiagonalSeq(fptr, maxoff, B->E, E, A);
	}
	return max;
}

void	UpdateHitsGPsiBlst(gpb_typ B)
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

Int4	MatcherGPsiBlstStr(Int4 len, register unsigned char *seq, register gpb_typ B)
/************************************************************************
    diag = [z,..., z, z, z, z, z,...,z ]  next = last = z.
    pos  =  z,...,-2,-1, 0, 1, 2,...,len 
      ed = [0,..., 0, 0, 0, 0, 0,...,0 ]  

		d=1;  diag[z]=1;  ed[1]=33;  next=1;
		d=-2; diag[1]=-2; ed[-2]=25; next=-2;
		d=0;  diag[-2]=0; ed[0]=9;   next=0;

    diag = [1,..., 0, z, z,-2, z,...,z ]  (list(z) = 1 -> -2 -> 0 -> z)
    pos  =  z,...,-2,-1, 0, 1, 2,...,len  
      ed = [0,...,25, 0, 9,33, 0,...,0 ]  
/************************************************************************/
{
	register Int4	i,d,*ed,q,j;
	Int4		*next,nxt,n,score,s,e,max,*tmp=B->tmp,item;
	char **R=AlphaR(B->A);

	while(DelMinMheap(B->mH) != NULL) ; B->update = TRUE;
	// n = LenSeq(B->E);
	assert(LenSeq(B->E) < MAX_SEQ_LENG_GPB);
	nxt=B->zero; next=B->diag; ed=B->extdiag; 
	q = B->d[seq[1]][0];		/** q = state "^c" **/
	for(max=0, i=2; i < len; i++){
	   q = B->d[seq[i]][q];		/** q = state "xd" **/
				/** seq[i+1] = next token (mealy model) **/
	   if((n=GetListMList(tmp,StateGPB(q,seq[i+1]),B->mlist))){ 
					/** then q + t signals acceptance **/
	     for(j=0;  j < n; j++){
		d = tmp[j] - (i-1);	/** get diagonal **/
		if(i > ed[d]){
		  // It looks like i2 here points to beginning of word hit.
		  score=ExtendGPsiBlstStr(LenSeq(B->E),tmp[j],len,seq,i-1,&s,&e,
			B->matrix,B->x_dropoff);
		  if(score > max) max = score;
		  if(ed[d] == 0){ next[nxt] = d; nxt = d; }
		  ed[d] = e; 
		  if(score >= B->cutoff 
				&& (item=InsertMheap(-score,B->mH))!=NULL){
			B->hit_s[item]=s; B->hit_e[item]=e; B->hit_d[item]=d; 
		  }
		}
	     } 
	   }
	}
	for(d=B->zero; next[d]!=B->zero; d=n){	// free extend list.
		n=next[d]; ed[n]=0; next[d]=B->zero; 
	}
	return max;
}

Int4	ExtendGPsiBlstStr(Int4 lenMaster, Int4 i1, Int4 len2, unsigned char *seq2,
	Int4 i2, Int4 *left, Int4 *right, Int4 **matrix, register Int4 x_dropoff)
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
	register Int4	i,s,max,end,**mtx;
	register unsigned char	*p2;

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

