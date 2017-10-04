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

#include "ndl_typ.h"

#define debug_NDL_penalty 0

char	*ndl_typ::GapAlnTraceP2P(p2p_typ *P2P, Int4 &start,Int4 &score)
// return the trace, the start position of the trace, and the trace length.
{
	Int4    i,j,newlen,lenTrace;
	char	*operation,*tmp;

	tmp=this->gap_alignP2P(P2P,&lenTrace,start,score);
	// put_seqaln_smatrixSW(stdout,tmp,LenSeq(E),SeqPtr(E),OffSetSeq(E),lenTrace,nmod,M);
#if 0
std::cerr << "tmp traceback"; std::cerr << std::endl;
std::cerr << tmp; std::cerr << std::endl;
#endif
	for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
	NEW(operation,newlen+3,char); operation[0]='E';
	for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
	operation[i]='E'; i++;
	operation[i]=0; free(tmp);
#if 0
std::cerr << "real traceback"; std::cerr << std::endl;
std::cerr << operation; std::cerr << std::endl;
std::cerr << *start; std::cerr << " = start\n";
#endif
	return operation;
}

#define debugger_on 0		// These routines are not working!!!
#define use_juns_hmm 1
//======================= FAST gapped align for sampling ================
char	*ndl_typ::gap_alignP2P(p2p_typ *P2P,Int4 *J,Int4 &start,Int4 &alnscore)
{
	Int4	m,i,j,k,n1,n2,score,g,gINS;
	Int4	**MAT,**T,**DEL,**INS;
	Int4	s_mm,s_mi,s_md,s_dd,s_dm,s_ii,s_im;
	Int4	s_di,s_id;	// Jun's HMM.
	char	*operation;
	char	flag=0;

	/** get total length of profile **/
	n1=P2P->RtnTotLenP(); n2=P2P->RtnLenC();
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
	}
	// Make GLOBAL ALIGNMENT with respect to profile,
/********************************************************************************
 ********************************************************************************/
	MAT[0][0] = 0; 		// corresponds to S0 state.
	// DEL[0][0] = d2d[1]-m2d[1];  // --> DEL[1][0] = m2d. note m2d == s2d...
	// DEL[0][0] = -m2d[1];  // --> DEL[1][0] = -s2d - d2d. note m2d == s2d...
	DEL[0][0] = 0;
	INS[0][0] = 0;
	T[0][0] = 0;	
	for(i=1; i<= n1; i++) {  // j==0 implies deletion of model position route.
	    DEL[i][0] = DEL[i-1][0] - d2d[i];	// penalty for deletions at end (-m2d too?)
	    INS[i][0] = -infinity;	// disallow internal inserts from start of sequence.
					// requires d2i transitions.
	    // MAT[i][0] = 0;
	    MAT[i][0] = -infinity;  
	    T[i][0] = T[i-1][0] + 1; // deletion traceback route 
	}
	for(j=1; j<= n2; j++) { 
	    	MAT[0][j] = 0;  
		DEL[0][j] = 0;     // Make LOCAL with respect to sequence.
		// DEL[0][j] = -infinity;
		// DEL[0][j] = - s2d[1];
		INS[0][j] = 0;	// no penalty for inserts before HMM...
		T[0][j]= T[0][j-1] -1; // insert traceback route 
		// T[0][j]= -1; // traceback insert route 
	} DEL[0][0] = 0;  
	Int4 *gdel; NEW(gdel,n2+3,Int4);
	// 2. Dynamic programming step. 
	for(k=m=i=1; i<= n1; i++) {
	   Int4	*mat_i=MAT[i],*mat_im1=MAT[i-1];
	   Int4 *ins_i=INS[i],*ins_im1=INS[i-1];
	   Int4 *del_i=DEL[i],*del_im1=DEL[i-1];
	   register Int4 jj,jm1;
	   for(gINS=0,jm1=0,jj=1; jj<= n2; jj++,jm1++) {  
		//========= 1. Choose optimum between m->i & i->i transitions. ==========
	   	Int4	t,s_max;
	        if(k==P2P->RtnLenP(m)){	// NO PENALTIES for inserts at ends!
		  s_mi = mat_i[jm1];	// insertion opening score
		  s_ii = ins_i[jm1];	// insertion extension score
	        } else {
                  s_mi=mat_i[jm1]-m2i[i];	// insertion opening score
		  s_ii=ins_i[jm1]-i2i[i];	// insertion extension score
		}
		// save the maximum score and trace back...
                if(s_mi > s_ii){ t=gINS=-1; ins_i[jj] = s_mi; s_max=s_mi; }
		else { gINS--; t=gINS; ins_i[jj] = s_ii; s_max=s_ii; }
		// else { t=-1; ins_i[jj] = s_ii; s_max=s_ii; }
#if use_juns_hmm		// Jun's HMM...
                if(k==P2P->RtnLenP(m)){
		   s_di = del_i[jm1] - d2m[i];	// no penalty: d2s = 0;
		   if(s_di > s_mi && s_di > s_ii){
		     gINS=-1; ins_i[jj]=s_di; t=-1; s_max=s_di; 
		   }
		}
#endif
if(flag==0){
if(ins_i[jj] > 1000000000){ fprintf(stderr,"INS[%d][%d]=%d; m=%d; k=%d\n",i,jj,ins_i[jj],m,k); flag=1; }
}

		//========== 2. Choose optimum between m->d & d->d transitions. ===========
		s_md=mat_im1[jj]-m2d[i];	// note m2d == s2d 
		s_dd=del_im1[jj]-d2d[i];	// deletion extension score
                if(s_md > s_dd){
			gdel[jj]=1;  del_i[jj]=s_md; 
			if(s_md > s_max){ s_max=s_md; t=1; }
		} else {
		 	gdel[jj]++; del_i[jj] = s_dd; 
			if(s_dd > s_max){ s_max=s_dd; t=gdel[jj]; }
			// if(s_dd > s_max){ s_max=s_dd; t=1; }
		}
#if use_juns_hmm		// Jun's HMM...
                if(k==1){
		   s_id = ins_im1[jj]-m2d[i];	// - s2d?
		   if(s_id > s_md && s_id > s_dd){
			gdel[jj]=1; del_i[jj]=s_id;
			if(s_id > s_max){ s_max=s_id; t=1; }
		   }
		}
#endif
if(flag==0){
if(del_i[jj] > 1000000000){ fprintf(stderr,"DEL[%d][%d]=%d; m=%d; k=%d\n",i,jj,del_i[jj],m,k); flag=1; }
}

		// 3. Choose optimum among m->m, i->m or d->m transitions.
		// score = P2P->RtnScore(jj,m,k);  // match score at i,j cell.
		score = P2P->RtnScore(jj,m,k);  // match score at i,j cell or k in m.
#if 0
if(score > 1000000000){ fprintf(stderr,"score=%d; i=%d; jj=%d; m=%d; k=%d\n",score,i,jj,m,k); }
if(score < -2000000000){ fprintf(stderr,"score=%d; i=%d; jj=%d; m=%d; k=%d\n",score,i,jj,m,k); }
#endif
	        if(k==1) {	// s2m PENALTY for match at start of motif!!!!
		  s_mm=mat_im1[jm1] + score - s2m[i];	// should be m2s???
		  s_im=ins_im1[jm1] + score;	// no s->m penalty.
	        } else {
		  s_mm=mat_im1[jm1] + score - m2m[i];
		  s_im=ins_im1[jm1] + score - i2m[i];
		}
		s_dm=del_im1[jm1]+score-d2m[i];

		if(s_mm >= s_dm){
			if(s_mm >= s_im){
				mat_i[jj] = s_mm; 
				if(s_mm >= s_max){ s_max=s_mm; t=0; }
			} else {
				mat_i[jj]=s_im; 
				if(s_im >= s_max){ s_max=s_im; t=0; }
			}
		} else if(s_dm >= s_im){
			mat_i[jj]=s_dm; 
			if(s_dm >= s_max){ s_max=s_dm; t=0; }
		} else {
			mat_i[jj]=s_im; 
			if(s_im >= s_max){ s_max=s_im; t=0; }
		} T[i][jj]=t;	// record optimum route to i,j cell.
if(flag==0){
   if(mat_i[jj] > 1000000000){
	fprintf(stderr,"MAT[%d][%d]=%d; m=%d; k=%d; T = %d\n",i,jj,mat_i[jj],m,k,T[i][jj]); flag=1; 
	fprintf(stderr,"m2m=%d; m2i=%d; i2i=%d; i2m=%d; m2d=%d; d2d=%d; d2m=%d; s2m=%d.\n",
                        m2m[i],m2i[i],i2i[i],i2m[i],m2d[i],d2d[i],d2m[i],s2m[i]);
	fprintf(stderr,"s_mm=%d; s_im=%d; s_dm=%d.\n",s_mm,s_im,s_dm);
	fprintf(stderr,"s_mi=%d; s_ii=%d; s_md=%d.\n",s_mi,s_ii,s_md);
   }
}
	   } k++;
	   if(k > P2P->RtnLenP(m)) { k = 1; m++; }
	} free(gdel);
	// 2b. ********** Find optimum global score with respect to profile. *******
	// 3. *************** Trace back step. ***************/
	// operations: 
	// 'i' = insertion in sequence outside of motifs
	// 'I' = insert in sequence within motif (need to delete this)
	//  'M' = match to start of a motif block
	//  'm' = match to other sites in motif
	//  'D' = deletion of sequence within motif
	//  'd' = deletion of sequence outside of motif (not used)
#if 0	//********************************************************************

#endif	//********************************************************************
	Int4	s,t,t0,J0,max_i,max_j,end;
	max_i=n1;		// global with respect to profile.
	for(score=INT4_MIN, j=1; j<= n2; j++){
		if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
// if(score > 1000000000){ fprintf(stderr,"MAT: score=%d; s=%d; j=%d; max_j=%d; max_i=%d\n",score,s,j,max_j,max_i); }
	}
#if 1	// Determine whether optimum alignment ends in a match or a deletion.
	for(j=1; j<= n2; j++){
		if((s=DEL[max_i][j]) > score){ score=s; max_j=j; }
// if(score > 1000000000){ fprintf(stderr,"DEL: score=%d; s=%d; j=%d; max_j=%d; max_i=%d\n",score,s,j,max_j,max_i); }
	}
#endif	// NOTE: insertions at end can't improve the score...
	NEW(operation,n1+n2+3,char); operation[0]='E';
	m=P2P->RtnNumBlks(); k=P2P->RtnLenP(m);
	for(J0=0,j=n2; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
	// for(i=n1; i > 0 || j > 0; ){  // full alignment...
	for(i=max_i; i > 0; ){  // full alignment...Stop at end of profile
	  t0 = T[i][j];
	  // fprintf(stderr,"T[%d][%d]=%d\n",i,j,t0);
	  do {
	    if(t0 > 0){ t=1; t0--; }
            else if(t0 < 0){ t=-1; t0++; } else { t=0; }
	    switch(t){
		case 0:	  // Match state found...
		   i--; J0++;
		   if(k==1) operation[J0] = 'M'; else operation[J0] = 'm';
		   j--; k--; break;
		case 1:   // Deletion: gap ('-') is in sequence; add X's for gaps
			J0++; 
			if(k==1) operation[J0] = 'D'; else operation[J0] = 'd';
			i--;  k--; break;
		case -1:  // Insertion: gap ('-') is in profile
			if(m < 1 || k==P2P->RtnLenP(m))
				{ J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			j--; break;
		default: print_error("this should not happen"); break;
	     }
	     if(k==0){ m--; if(m > 0) k=P2P->RtnLenP(m); }
	  } while(t0 != 0);
	}
	start=j+1; J0++; operation[J0]='E';
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(T); free(DEL); free(INS); 
	alnscore = score; *J = J0;
        return operation;
}


