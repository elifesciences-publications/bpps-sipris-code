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

#include "gblast.h"

gb_typ MkGBlast2(Int4 hpsz,Int4 cutoff,Int4 T,e_type E,a_type A,Int4 xdrop)
{ gb_typ B=MkGBlast(hpsz,cutoff,T,E,A); B->x_dropoff=xdrop; return B; }

gb_typ MkGBlast(Int4 hpsz, Int4 cutoff, Int4 T, e_type E, a_type A)
{
	gb_typ	B;
	Int4	n,i,z;

	B = MakeGBlast(T, E, A);
	n = LenSeq(E);
	MEW(B->diag0,n+MAX_SEQ_LENG_GB+5,Int4);
	B->diag = B->diag0 + MAX_SEQ_LENG_GB;
	NEW(B->ed0,n+MAX_SEQ_LENG_GB+5,Int4);
	B->extdiag = B->ed0 + MAX_SEQ_LENG_GB;
	z = B->zero = - MAX_SEQ_LENG_GB;
	for(i = z; i <= n; i++) B->diag[i] = z;
	B->cutoff = cutoff;
	B->mH = Mheap(hpsz,3);
	NEW(B->hit_s,hpsz+5,Int4);
	NEW(B->hit_e,hpsz+5,Int4);
	NEW(B->hit_d,hpsz+5,Int4);

	NEW(B->hit_s2,hpsz+5,Int4);
	NEW(B->hit_e2,hpsz+5,Int4);
	NEW(B->hit_d2,hpsz+5,Int4);
	B->nhits = 0;

	NEW(B->score,hpsz+3,double);
	return B;
}

static void	calc_maxhits_gblast(a_type A)
/** calculate the maximum number of words that can possibly score >= T
    at an arbitrary sequence site. (Needed for allocating lists.) **/
{
	unsigned char	x1,x2,x3,y1,y2,y3;
	Int4	i,s,score[55],s1,s2,s3,max[55];

	for(i=0; i<=50;i++) max[i]=0;
	for(x1=0;x1<=nAlpha(A); x1++){
	  for(x2=0;x2<=nAlpha(A); x2++){
	    for(x3=0;x3<=nAlpha(A); x3++){
		for(i=0; i<=51;i++) score[i]=0;
		for(y1=0;y1<=nAlpha(A); y1++){
		  s=valAlphaR(x1,y1,A);
		  for(y2=0;y2<=nAlpha(A); y2++){
		    s+=valAlphaR(x2,y2,A);
		    for(y3=0;y3<=nAlpha(A); y3++){
		        s+=valAlphaR(x3,y3,A);
			if(s>=0){ score[s]++; }
		        s-=valAlphaR(x3,y3,A);
		    }
		    s-=valAlphaR(x2,y2,A);
		  }
		}
		for(i=50; i>=0;i--) score[i]+=score[i+1];
		for(i=0; i<=50;i++) if(score[i]>max[i]) max[i]=score[i];
	    }
	  }
	}
	fprintf(stdout,"Int4 maxhits[50] = { ");
	for(i=0; i<=50;i++){
		fprintf(stdout,"%d, ",max[i]);
		if(i%10 ==9) printf("\n");
	} printf("}; \n");
	print_error("calc max hits....");
}

gb_typ MakeGBlast(Int4 T, e_type E, a_type A)
{
	gb_typ B;
	Int4	n,q,q0,q1,q2,r0,r1,r2,s,t,c,d,e,i,x,y,z,hits,nlet;
	char	**best,*b0,*b1,*b2; // letters ordered by decreasing score.
	keytyp	key;
	dh_type	H;
	Int4        time1;
	// Eventually put this into alphabet data structure:
	Int4	blsm62max[40] = { 2418, 1826, 1540, 1264, 1130, 995, 
			759, 627, 458, 304, 222, 164, 105, 83, 70, 61, 
			58, 58, 58, 49, 28, 15, 10, 7, 4, 2, 1, 1, 1, 
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}; 

	if(T >= 34) print_error("MakeGBlast( ) T too large; no hits possible.");
	NEW(B,1,gblast_type);
	B->x_dropoff=22;
	B->A = A; B->E = E; B->mH = NULL;
	B->diag0 = B->ed0 = NULL;
	B->hit_s = B->hit_e = B->hit_d = NULL;
	B->hit_s2 = B->hit_e2 = B->hit_d2 = NULL; B->score = NULL;
	/***** order letters by decreasing score *****/
	NEWP(best,nAlpha(A)+2,char);
	H = dheap(nAlpha(A)+3,3);
	for(i=0,c=0; c <= nAlpha(A); c++){
	    NEW(best[c],nAlpha(A)+2,char);
	    for(d=0; d <= nAlpha(A); d++){
		key = (keytyp) -valAlphaR(c,d,A);
		insrtHeap(d+1,key, H);
	    } // fprintf(stderr,"%c: ", AlphaChar(c,A));
	    for(i=0; (d = delminHeap(H)) != NULL; i++){
		best[c][i] = d - 1;
	        // fprintf(stderr,"%c ", AlphaChar(d-1,A)); 
	    } // fprintf(stderr,"\n"); 
	}
	Nildheap(H);
	/***** create finite automaton *****/
	/*** create transition function d(r,q). ***/
	n = nAlpha(A) +1; n = 1 + n + n*n;  
	NEWP(B->d,nAlpha(A)+2,Int4);
	for(c=0; c <= nAlpha(A); c++) NEW(B->d[c],n+2,Int4);
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q++; B->d[c][0] = q; q0=q;		/** q0 = state "^c" **/
	    for(d=0; d <= nAlpha(A); d++){
	    	q++; B->d[d][q0] = q;        	/** q1 = state "cd" **/
	    }
	}
	B->nQ = q;				/** nQ = #states **/
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q0 = B->d[c][0];			/** q = state "^c" **/
	    for(d=0; d <= nAlpha(A); d++){
	    	q1 = B->d[d][q0];		/** q = state "cd" **/
	    	q2 = B->d[d][0];		/** q2 = state "^d" **/
		for(e=0; e <= nAlpha(A); e++){
			q = B->d[e][q2];	/** q2 = state "de" **/
			B->d[e][q1] = q;	/** "cde" = "^de" = "de" **/
		}
	    }
	} // create acceptance states & lists[q][r][1..].
	time1=time(NULL);
	nlet = nAlpha(A) +1;
	n = B->nQ;
	B->mlist=MkMList(blsm62max[T]*LenSeq(E),MaxStateGB(nlet));
        // fprintf(stderr,"\nallocation time: %d seconds\n", time(NULL)-time1);
	time1=time(NULL);

	n = LenSeq(E) - 1;
	NEW(B->tmp,n+3,Int4); 
	c = XSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(hits=0,i=2; i <= n; c=d, i++){
	   d = XSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XSeq(i+1,E);
#if 0
fprintf(stderr,"%c%c%c: (q=%d; d(%c,q)=%d)", 
		AlphaChar(c,A),AlphaChar(d,A),AlphaChar(e,A),
		q,AlphaChar(e,A),B->d[e][q]);
#endif
	   // look through related 3-words for "matches".
	   b0 = best[c]; b1 = best[d]; b2 = best[e];
	   for(x=0; x <= nAlpha(A); x++){
	        r0 = b0[x]; q0 = B->d[r0][0];
		s = T - (Int4) valAlphaR(c,r0,A);
		for(y=0; y <= nAlpha(A); y++){	
	           r1 = b1[y]; q1 = B->d[r1][q0];
		   t = s - (Int4) valAlphaR(d,r1,A);
		   for(z=0; z <= nAlpha(A); z++){	
			r2 = b2[z];
			if(t <= (Int4)valAlphaR(e,r2,A)){
				Add2MList(i-1,StateGB(q1,r2), B->mlist);
				hits++; // add to acceptor states
#if 0	
fprintf(stderr," %c%c%c(%d)", AlphaChar(r0,A), AlphaChar(r1,A),
					AlphaChar(r2,A),s+t+u);
#endif
			} else { break; } 
		   } if(z==0) break;
		}
		if(y==0) break;
	   } //fprintf(stderr,"\n");  
	}
        // fprintf(stderr,"\nneighborhood time: %d seconds\n", time(NULL)-time1);
	for(c=0; c <= nAlpha(A); c++) free(best[c]); free(best);
#if 0
	fprintf(stderr,"hits = %d\n",hits); 
        fprintf(stderr,"\tcreate time: %d seconds\n", time(NULL)-time1);
#endif
	return B;
}

void	NilGBlast(gb_typ B)
{
	Int4	c;

	for(c=0; c <= nAlpha(B->A); c++) {
		free(B->d[c]);
	}
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

Int4	MatcherGBlast(FILE *fptr, e_type E, gb_typ B)
{
	Int4	num,n,q,score,s,c,d,e,i,j,hits=0,max,*tmp=B->tmp,maxoff;
	a_type	A=B->A;
	Int4        time1,lft,rt;

	time1=time(NULL);
	n = LenSeq(E) - 1; c = XSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(max=0, i=2; i <= n; c=d, i++){
	   d = XSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XSeq(i+1,E);		/** e = next token (mealy model) **/
				/** then q + e signals acceptance **/
	   if(!EmptyMList(StateGB(q,e),B->mlist)){
		num=GetListMList(tmp, StateGB(q,e), B->mlist);
		for(j=0;  j < num; j++){
			s = tmp[j];
			score = ExtendGBlastStr(B->E,s,LenSeq(E), XSeqPtr(E),
					i-1, &lft, &rt, AlphaR(A),B->x_dropoff);
			if(score > max){ max = score; maxoff = s-i+1; }
		} hits+=j;
	   }
	}
	if(max > 0 && fptr != NULL) {
		fprintf(fptr,"hits = %d; max score = %d\n",hits,max); 
		fprintf(fptr,"score = %d; offset = %d\n", max,maxoff);
		PutDiagonalSeq(fptr, maxoff, B->E, E, A);
        	fprintf(fptr,"\tmatcher time: %d seconds\n", 
			time(NULL)-time1);
	}
	return max;
}

BooLean	FastMatcherGBlastStr(unsigned char *seq, Int4 length, gb_typ B, Int4 score)
/************************************************************************
  GBlast for coiled coils heuristic program.
/************************************************************************/
{
	Int4	num,n,q,s,c,d,e,i,j,*tmp=B->tmp;

	n = length - 1; 
	c = seq[1];
	q = B->d[c][0];		/** q = state "^c" **/
	for(i=2; i <= n; c=d, i++){
	   d = seq[i];
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = seq[i+1];		/** e = next token (mealy model) **/
					/** then q + e signals acceptance **/
	   if(!EmptyMList(StateGB(q,e),B->mlist)){ 
		num=GetListMList(tmp,StateGB(q,e), B->mlist);
		for(j=0;  j < num; j++){
		   s = tmp[j];
		   if(FastExtendGBlastStr(B->E,s,length, seq,i-1,
		   		AlphaR(B->A),score)) return TRUE;
		}
	   }
	}
	return FALSE;
}

BooLean	FastMatcherGBlast(e_type E, gb_typ B, Int4 score)
/************************************************************************
  GBlast for purge program.
/************************************************************************/
{
	Int4	num,n,q,s,c,d,e,i,j,*tmp=B->tmp;

	n = LenSeq(E) - 1; 
	c = XSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(i=2; i <= n; c=d, i++){
	   d = XSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XSeq(i+1,E);		/** e = next token (mealy model) **/
					/** then q + e signals acceptance **/
	   if(!EmptyMList(StateGB(q,e),B->mlist)){ 
		num=GetListMList(tmp,StateGB(q,e), B->mlist);
		for(j=0;  j < num; j++){
		   s = tmp[j];
		   if(score <= ExtendGBlast(B->E, s, LenSeq(E), XSeqPtr(E),i-1,
		   			AlphaR(B->A))) return TRUE;
#if 0
		   if(FastExtendGBlastStr(B->E,s,LenSeq(E), XSeqPtr(E),i-1,
		   		AlphaR(B->A),score)) return TRUE;
#endif
		}
	   }
	}
	return FALSE;
}

Int4	ExtendGBlast(e_type E1, Int4 i1, Int4 len2, unsigned char *seq2,
	Int4 i2, register char **R)
/****************************************************************

 ****************************************************************/
{
	register Int4	i,s,max,end;
	register unsigned char	*p1,*p2;

	p1 = XSeqPtr(E1)+i1; p2=seq2+i2; 
	if((max=len2-i2) > (s=LenSeq(E1)-i1)){ end=s; } else { end=max; }
	s = R[p1[1]][p2[1]];
	for(max=s, i=2; i <= end; i++) {
		if((s += R[p1[i]][p2[i]]) > max) max=s;
		else if((s <= (max-22))) break; 
	}
	if(i1<=i2) end=-i1; else end=-i2;
	for(s=max, i=0; i > end; i--){
		if((s += R[p1[i]][p2[i]]) > max) max=s; 
		else if((s <= (max-22)) || s < 0) return max; 
	}
	return max;
}

BooLean	FastExtendGBlastStr(e_type E1, Int4 i1, Int4 len2, unsigned char *seq2,
	Int4 i2, register char **R, Int4 score)
{
	register Int4	i,s,max,end;
	register unsigned char	*p1,*p2;

	p1 = XSeqPtr(E1)+i1; p2=seq2+i2; 
	if((max=len2-i2) > (s=LenSeq(E1)-i1)){ end=s; } else { end=max; }
	s = R[p1[1]][p2[1]];
	for(max=s, i=2; i <= end; i++) {
		if((s += R[p1[i]][p2[i]]) > max) { 
			if(s >= score) return TRUE; 
			max=s; 
		} else if((s <= (max-22))) break; 
	}
	if(i1<=i2) end=-i1; else end=-i2;
	for(s=max, i=0; i > end; i--){
		if((s += R[p1[i]][p2[i]]) > max) { 
			if(s >= score) return TRUE; 
			max=s;
		} else if((s <= (max-22)) || s < 0) return FALSE; 
	}
	return FALSE;
}

void	UpdateHitsGBlast(gb_typ B)
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

double	SumStatGBlastStr(double N, double M, double lambda,
	double H, double K, gb_typ B)
{
      Int4	i,r,*hsp;
      double	sum,d,adj;

      if(B->update) UpdateHitsGBlast(B);
      if(B->nhits == 0) return 1.0;
      NEW(hsp,B->nhits+3,Int4);
      r=ConsistentGBlast(N,M, hsp, lambda,H,K,B);
      for(sum=0.0,i=r; i>0; i--){ sum+=(double)B->score[hsp[i]]; }
      sum *= lambda;
      /** printf("sum = %g; r = %d; r! = %g\n",sum,r,factrl(r)); /****/
      /******************************************************
	Adjust for edge effects: compute the expected length 
	for HSPs & subtract from actual lengths.
      /******************************************************/
      if(r>1) d = (double)((r-1)*MAX_OVERLAP_GB)*H/2.0; 
      else d = 0.0;   /** if r > 1 then allow for some overlap. **/
      adj = (sum-d)/H;
      
      N = MAXIMUM(double,(N-adj),1.0);
      M = MAXIMUM(double,(M-adj),1.0);
      sum -= (double)r*log(K*N*M);
      if(r>1) sum+=lngamma(r+1.0); /** ln r! ordering adjustment **/
      free(hsp);
      /** use Phil Green's adjustment for multiple tests. **/
      return (pow(2.0,(double)r)*SumPStd(r,sum));
}

Int4	ConsistentGBlast(Int4 N, Int4 M, Int4 *segs, double lambda,
        double H, double K, gb_typ B)
/************************************************************************

                           cx1                     
            _________________________________
          s|               .                 |   Attempts to find an optimum
          e|           s1  .                 |   arrangement of r_opt 
          q|            \  .                 |   segments that maximizes
           |             \ .                 |   the significance.
          Y|              \.                 |                        
       cy1 |...............*   (len1 of r1)  |   The selected HSPs are given
           |                \                |   in seg[1..r_opt] and 
           |             s2  \               |   r_opt is returned.   
           |              \   \              |                        
           |               \   e1            |   A topological sort   
           |                \                |   algorithm is used to 
       cy2 |.................*               |   find optimal.      
           |                 .\              |                         
           |                 . \ (len2 of r2)|   overlapx=4 & overlapy=1.
           |                 .  \            |   and len1 = len2 = 7. 
           |                 .   e2          |   overlap = 4 (i.e., max). 
           |_________________________________|
             seq X = querry
                             cx2

/************************************************************************/
{
	Int4	score,r,r_opt,s,e;
	Int4	*dist,*path,maxoverlap=MAX_OVERLAP_GB,maxfrac=MAX_FRACTION_GB;
	Int4	overlap,overlapx,overlapy,len1,len2;
	Int4	r1,r2,e1,e2,s1,s2,d1,d2,cx1,cy1,cx2,cy2;
	double	adj;
	wdg_typ	G;

	if(B->update) UpdateHitsGBlast(B);
	if(B->nhits < 1) print_error("this should not happen");
	/** adjust for length using just the MSP score **/
	s=(Int4)floor(((double)B->score[1]*lambda/H)+0.5);
	N = MAXIMUM(Int4,(N-s),1);
	M = MAXIMUM(Int4,(M-s),1);
	adj = log(K*M*N);  /** adjust score for search space **/
	NEW(dist,B->nhits+5,Int4);
	NEW(path,B->nhits+5,Int4);
	e = B->nhits+2;  e *= e;
	G = MkWdgraph(B->nhits+2,e);
	s = B->nhits+1; e = B->nhits+2;
	/*** add starting & ending arcs to graph ****/
	for(r=1; r <= B->nhits; r++){
	   JoinWdgraph(s, r, 0, G);
	   score=(Int4)floor(10*((double)B->score[r]*lambda-adj)+0.5);
	   JoinWdgraph(r, e, -score, G);
	}
	/*** add legitimate arcs to graph ****/
	for(r1=1; r1 < B->nhits; r1++){
	  s1=B->hit_s2[r1]; e1=B->hit_e2[r1]; d1=B->hit_d2[r1];
	  len1=e1-s1;
	  cx1 = (s1+len1/2); cy1 = cx1+d1;
	  for(r2=r1+1; r2 <= B->nhits; r2++){
	    s2=B->hit_s2[r2]; e2=B->hit_e2[r2]; d2=B->hit_d2[r2];
	    len2 = e2-s2;
	    cx2 = (s2 + len2/2); cy2 = cx2+d2;
	    /*** fprintf(stderr,"(%d:%d) cx1=%4d; cx2=%4d; cy1=%4d; cy2=%4d\n",
			r1,r2,cx1,cx2,cy1,cy2); /******/
	    if(cx1 < cx2 && cy1 < cy2){  	/** i.e.,  r1 < r2 **/
	      overlap = MAXIMUM(Int4,(e1-s2),((e1+d1)-(s2+d2)));
	      if(overlap <= maxoverlap){
		if(overlap < 1 || 
			((len1/overlap)>=maxfrac && (len2/overlap)>=maxfrac)){
		  score=(Int4)
		    floor(10*((double)B->score[r1]*lambda-adj)+0.5);
		  JoinWdgraph(r1, r2, -score, G);
		}
	      }
	    } else if(cx1 > cx2 && cy1 > cy2){  /** i.e., r2 < r1 **/
	      overlap = MAXIMUM(Int4,(e2-s1),((e2+d2)-(s1+d1)));
	      if(overlap <= maxoverlap){
		if(overlap < 1 || 
			((len1/overlap)>=maxfrac && (len2/overlap)>=maxfrac)){
		  score=(Int4)
		    floor(10*((double)B->score[r2]*lambda-adj)+0.5);
		  JoinWdgraph(r2, r1, -score, G);
		}
	      }
	    }
	  }
	}
	// PutWdgraph(stdout, G); 
	TopoScanWdigraph(G, s, path, dist);
	// fprintf(stderr,"\nstart "); 
	for(r_opt=0,r=path[e]; r!=s; r=path[r]){ r_opt++; segs[r_opt]=r; }
	// *debug** { fprintf(stderr,"-> %d",r); } 
	// fprintf(stderr," (dist = %d)\n",-dist[e]); 
	NilWdgraph(G); free(dist); free(path);
	return r_opt;
}

#if 0	// NOT CURRENTLY WORKING...
/****************************************************************
 FOR DATABASE SEARCH ON STRINGS:
 ****************************************************************/
e_type	SubSeqGBlast(char *id, Int4 len, unsigned char *seq, Int4 flank, gb_typ B)
{
      Int4	r,s,e,d,max,min,min_d,max_d;
      Int4 	i,r_opt,*hsp;

      if(B->update) UpdateHitsGBlast(B);
      NEW(hsp,B->nhits+3,Int4);
      r_opt=ConsistentGBlast(LenSeq(B->E),len, hsp, B);

      for(max=0,min=len,i=r_opt; i > 0; i--){
	r = hsp[i];
	s = B->hit_s2[r]; e = B->hit_e2[r]; d = B->hit_d2[r]; 
	if(e > max) { max = e; max_d = d; }
	if(s < min) { min = s; min_d = d; }
      }
      max = MINIMUM(Int4,max+flank,len);
      min = MAXIMUM(Int4,min-flank,1);
      free(hsp);
      return MkSeq(id, (max-min+1), (seq+min-1));
}

Int4    ColBlcksGBlastStr(Int4 N, Int4 M, Int4 number, unsigned char *seq,
	Int4 *col, char **R, gb_typ B)
{
      Int4	i,j,s,e,d,r,r_opt,*segs,c;

      if(B->update) UpdateHitsGBlast(B);
      NEW(segs,B->nhits+3,Int4);
      r_opt=ConsistentGBlast(N,M, segs, B);

      for(c=0,i=r_opt; i > 0; i--){
	r = segs[i];
	s = B->hit_s2[r]; e = B->hit_e2[r]; d = B->hit_d2[r]; 
	for(j=0; s <= e; s++){
	   if(R[ResSeq(s+d,B->E)][seq[s]] > 0) c++;
	}
      } 
      free(segs);
      *col = c;
      return r_opt;
}

Int4    PutHitsGBlastStr(FILE *fp, Int4 N, Int4 M, Int4 dbslen, gb_typ B)
{
      Int4	i,s,e,x,d,r,*segs,adj,len=M;
      double	sum,pval,snorm;

      if(B->update) UpdateHitsGBlast(B);
      NEW(segs,B->nhits+3,Int4);
      r=ConsistentGBlast(N,M, segs, B);

      for(sum=0.0,i=r; i > 0; i--){ sum+=(double)B->score[segs[i]]; }
      sum *= B->lambda;
      /******************************************************
	Adjust for edge effects: compute the expected length 
	for HSPs & subtract from actual lengths.
        Note: need to convert score to bits to calc. edge effect 
	because H is given in bits.
      /******************************************************/
      adj = (Int4) floor((sum/B->H)+0.5);
      N = MAXIMUM(Int4,(N-adj),1);
      M = MAXIMUM(Int4,(M-adj),1);
      sum -= (double)r*log(B->K*N*M);
      if(r>1) sum+=lngamma(r+1.0); /** ln r! ordering adjustment? **/
      /** use Phil Green's adjustment for multiple tests. **/
      pval = (pow(2.0,(double)r)*SumPStd(r,sum));

      for(i=r; i > 0; i--){
	x = segs[i];
        snorm = B->lambda*B->score[x]-log(B->K*N*M);
	s = B->hit_s2[x]; e = B->hit_e2[x]; d = B->hit_d2[x]; 
        fprintf(fp," %5d(%5d): %5d %5d %5d\n",
			(Int4)B->score[x],(Int4)snorm,s,e,d);
        fprintf(fp,"        query: %5d %5d\n",s+d,e+d);
      } 
      fprintf(fp,"SumP(%d,%d) = %g; N = %d; M = %d\n",
		r,(Int4)sum,pval*(double)dbslen/(double)len,N,M);
      free(segs);
      return B->nhits;
}
#endif

Int4	MatcherGBlastStr(Int4 len, register unsigned char *seq, 
	register gb_typ B, char **R)
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

	assert(len < MAX_SEQ_LENG_GB);
	while(DelMinMheap(B->mH) != NULL) ; B->update = TRUE;
	n = LenSeq(B->E);
	nxt=B->zero; next=B->diag; ed=B->extdiag; 
	q = B->d[seq[1]][0];		/** q = state "^c" **/
	for(max=0, i=2; i < len; i++){
	   q = B->d[seq[i]][q];		/** q = state "xd" **/
				/** seq[i+1] = next token (mealy model) **/
	   if((n=GetListMList(tmp,StateGB(q,seq[i+1]),B->mlist))){ 
					/** then q + t signals acceptance **/
	     for(j=0;  j < n; j++){
		d = tmp[j] - (i-1);	/** get diagonal **/
		if(i > ed[d]){
		  score=ExtendGBlastStr(B->E,tmp[j],len,seq,i-1,&s,&e,R,B->x_dropoff);
		  if(score > max) max = score;
#if 0   // ************************* DEBUG ***************************
		  if(score >= B->cutoff) 
		    fprintf(stderr,"\n%d; s=%d; e=%d; q=%d; dbs=%d; ed[%d]=%d\n",
			score,s,e,tmp[j],i,d,ed[d]);
#endif	//************************* DEBUG ***************************
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

Int4	ExtendGBlastStr(e_type E1,Int4 i1,Int4 len2,unsigned char *seq2,
	Int4 i2,Int4 *left,Int4 *right,register char **R,register Int4 x_dropoff)
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
	register Int4	i,s,max,end;
	register unsigned char	*p1,*p2;

	p1 = XSeqPtr(E1)+i1; p2=seq2+i2; 
	if((max=len2-i2) > (s=LenSeq(E1)-i1)){ end=s; } else { end=max; }
	*right = 1;
	s = R[p1[1]][p2[1]];
	for(max=s, i=2; i <= end; i++) {
		if((s += R[p1[i]][p2[i]]) > max) { *right=i; max=s; }
		else if((s <= (max-x_dropoff)) || s < 0) break; 
	}
	*right += i2; *left = 1;
	if(i1<=i2) end=-i1; else end=-i2;
	for(s=max, i=0; i > end; i--){
		if((s += R[p1[i]][p2[i]]) > max) { *left=i; max=s; }
		else if((s <= (max-x_dropoff)) || s < 0) break; 
	} *left += i2;
	return max;
}

/******************************* private *******************************/

void	gblast_error(char *s) { fprintf(stderr,"gblast: %s\n",s); exit(1); }

