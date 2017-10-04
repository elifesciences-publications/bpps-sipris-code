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

///////////////// HMM-ALIGNMENT WITH TRANSITION PROBABILITIES //////////////////
//
// Perform an alignment on multiple smatrix *M and sequence seq2
// returns optimum subalignment score.  
//
// Find maximum score and trace back to get alignment.
// 			D = DEL; M = MATCH; I = INS;
//             .....
//             :DMI:    
//             .....    gap = TB[m][j]
//      .......m....... | ..........m+1.......... 
//     : D   R   F   R :v: G   L   V   L   I   S : (concensus)
//  ....................................                
//   , : ,,: ,,: ,,: ,,: : ,,: ,,: ,,:             0          
//  ...:...:...:...:...:.:...:...:...:...                'a' = affine penalty
// I,, :   :   :   :   : :   :   : I : I :         1
//  ...:...:...:...:...:.:...:...:.|.\.|.:               ' ' = undefined.
// N,, :   :   :   :   : :   :   : H :\H :         2
//  ...:...:...:...0...:.:...:...:.|.:.|.:               '*' = NEG_INF.
// S,, :   :   : ->1`?--+->? :   : a : a :         :
//  ...:...:...:...:...:|:...\...:.|.:.|.:.........      ',' = zero.
// A,, :   :   :   :   :|:   : ? : d_: d :   :   : j-2
//  ...:...:...:...:...:|:...:...:...\.|.:...:...:.     S
// V,, :   :   :   :   :|:   :   :   :\H :   :   : j-1  E
//  ...:...:...:...:...:|:...:...:...:.|.:...:...:.     Q
// L,, :   :   :   :   :|:   :   :   : e :   :   : j    U
//  ...:...:...:...:...:|....:...:...:.|.:...:...:.     E
// I   :   :   :   :   :+->? .   :   : k_:   :   : j+1  N
//             :...:...:.:...\...:...:...0...:...:.     C
// S               :   : :   : ? :   :   :\G_:   : j+2  E
//               ..:...:.:...:...:...:...:...\...:.
// P               :   : :   :   :   :   :   :   : j+3
//                .:...:.:...:...:...:...:...:...:.
// N               :   : :   :   :   :   :   :   : j+4 (n2)
//                 :...:.:...:...:...:...:...:...:.
//   0   1   2  ... i-2  i-1  i  i+1 i+2 i+3 i+4
//   k = 1   2  ...lenM   1   2   ...       lenM
//                         MODEL            (n1)
//
//        AAAA   BBBBBB
//        DRFR...GLVLIS model concensus
//          .    ..::::
//       XXINLAESAVLISKHI sequence
//
//                  :
//    INS: [////]---|   Don't use this for k == lenM!
//         [\\\\\\\]|
//
//                  :
//    DEL: [///////]|   Don't use this for k == 1
//         [\\\\]---|
//
//                  :
//    MAT: [///////]|
//         [\\\\\\\]|
//
///////////////////////////////////////////////////////////////////

char	*ndl_typ::SampleGapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score)
// return the trace, the start position of the trace, and the trace length.
{
	Int4    i,j,newlen,lenTrace;
	char	*operation,*tmp;
	assert(nmod == 1);

	tmp=this->sample_gap_align(LenSeq(E),SeqPtr(E),M[1],&lenTrace,score,start);
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

//======================= FAST gapped align for sampling ================
char	*ndl_typ::sample_gap_align(Int4 sq_len,unsigned char *seq2,smx_typ smx,
			Int4 *J,Int4 *alnscore,Int4 *start)
#if 0	// 2b. ********** Find optimum global score with respect to profile. *******
operations: 
	'i' = insertion in sequence outside of profile
	'I' = insert in sequence within motif (need to delete this)
	'M' = match to start of a motif block
	'm' = match to other sites in motif
	'D' = deletion of sequence within motif
	'd' = deletion of sequence outside of motif (not used)
#endif
{
	Int4	m,i,j,hmm_len,g,gINS;
	double	**MAT,**DEL,**INS;
	char	*operation;

	hmm_len=LenSMatrix(smx); 	// get length of profile.
#if 1
	double	*M2M,*M2I,*M2D;
	double	*I2M,*I2I;
	double	*D2M,*D2D;
	double	*S2M;
	MEW(M2M,hmm_len+3,double); MEW(M2D,hmm_len+3,double); MEW(M2I,hmm_len+3,double); 
	MEW(I2M,hmm_len+3,double); MEW(I2I,hmm_len+3,double);
	MEW(D2D,hmm_len+3,double); MEW(D2M,hmm_len+3,double); MEW(S2M,hmm_len+3,double); 
	for(i=1; i<= hmm_len; i++) { 
	   M2M[i]=(double)m2m[i]/(double)pernats;
	   M2I[i]=(double)m2i[i]/(double)pernats;
	   M2D[i]=(double)m2d[i]/(double)pernats;
	   I2I[i]=(double)i2i[i]/(double)pernats;
	   I2M[i]=(double)i2m[i]/(double)pernats;
	   D2M[i]=(double)d2m[i]/(double)pernats;
	   D2D[i]=(double)d2d[i]/(double)pernats;
	   S2M[i]=(double)s2m[i]/(double)pernats;
	}
#endif
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,hmm_len+3,double*); 
	MEW(DEL,hmm_len+3,double*); MEW(INS,hmm_len+3,double*);
	for(i=0; i<= hmm_len; i++) { 
		MEW(MAT[i],sq_len+3,double); 
		MEW(DEL[i],sq_len+3,double); MEW(INS[i],sq_len+3,double); 
	}
	// Make GLOBAL ALIGNMENT with respect to profile,
/********************************************************************************
  j == sequence. i == HMM.
       DEL[i][j]               INS[i][j]              MAT[i][j]
       j\i: 0 : 1 : 2 : 3 :    j\i: 0 : 1 : 2 : 3 :   j\i: 0 : 1 : 2 : 3 :
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        0 | 0 |-dd|-dd|   |     0 | 0 |inf|inf|inf|    0 | 0 |inf|inf|inf|
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        1 |inf|   |   | M |     1 | 0 |   |   |   |    1 | 0 | M |   |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        2 |inf|   |   |   |     2 | 0 |   |   |   |    2 | 0 |   | M |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        3 |inf|   |   |   |     3 |   | M |   |   |    3 | 0 |   |   | M |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
 HMM=i: MWYYWCFWFCFWWF	           MWYYWCFWFCFWWF        MWYYWCFWFCFWWF
          ::::::::::::             ::::::::::::::        ::::::::::::::
 seq=j: --YYWCFWFCFWWF         (yf)MWYYWCFWFCFWWF        MWYYWCFWFCFWWF

 ********************************************************************************/
	MAT[0][0] = 0; 		// corresponds to S0 state.
	// DEL[0][0] = D2D[1]-M2D[1];  // --> DEL[1][0] = M2D. note M2D == S2D...
	// DEL[0][0] = -M2D[1];  // --> DEL[1][0] = -S2D - D2D. note M2D == S2D...
	DEL[0][0] = 0;
	INS[0][0] = 0;
	for(i=1; i<= hmm_len; i++) {  // j==0 implies deletion of model position route.
	    DEL[i][0] = DEL[i-1][0] - D2D[i];	// penalty for deletions at end (-M2D too?)
	    INS[i][0] = -infinity;	// disallow internal inserts from start of sequence.
					// requires D2I transitions.
	    // MAT[i][0] = 0;
	    MAT[i][0] = -infinity;  
	}
	for(j=1; j<= sq_len; j++) { // Make LOCAL with respect to sequence.
	    	MAT[0][j] = 0;  
		DEL[0][j] = 0;
		// DEL[0][j] = -infinity;
		INS[0][j] = 0;	// no penalty for inserts before HMM...
	} DEL[0][0] = 0;  
	double	dd,DD,score,s_mm,s_mi,s_md,s_dd,s_dm,s_ii,s_im;

h_type	HG=0;
// HG=Histogram("gapped scores",-500,500,50);
	//=========== 2. Dynamic programming step.  ===================
	double *mat_i,*mat_im1;
	double *ins_i,*ins_im1;
	double *del_i,*del_im1;
	register Int4 jj,jm1;
	for(i=1; i<= hmm_len; i++) {
	   mat_i=MAT[i]; mat_im1=MAT[i-1];
	   ins_i=INS[i]; ins_im1=INS[i-1];
	   del_i=DEL[i]; del_im1=DEL[i-1];
	   for(gINS=0,jm1=0,jj=1; jj<= sq_len; jj++,jm1++) {  

		//========= 1. Compute probabilities for each transition. ===========
		//========= 1a. Transitions into insert state. ===========
	        if(i==hmm_len){		// NO PENALTIES for inserts at ends!
		  s_mi = mat_i[jm1];	// insertion opening score
		  s_ii = ins_i[jm1];	// insertion extension score
	        } else {
                  s_mi=mat_i[jm1] - M2I[i];	// insertion opening score
		  s_ii=ins_i[jm1] - I2I[i];	// insertion extension score
		}

		//========= 1b. Transitions into delete state. ===========
		s_md=mat_im1[jj] - M2D[i];	// note M2D == S2D 
		s_dd=del_im1[jj] - D2D[i];	// deletion extension score

		//========= 1c. Transitions into match state. ===========
		score = (double)ValSMatrix(i,seq2[jj],smx)/(double)pernats;  // match score at i,j cell.
	        if(i==1) {	// S2M PENALTY for match at start of hmm!!!!
		  s_mm=mat_im1[jm1] + score - S2M[i];	// should be M2S???
		  s_im=ins_im1[jm1] + score;		// no s->m penalty.
	        } else {
		  s_mm=mat_im1[jm1] + score - M2M[i];
		  s_im=ins_im1[jm1] + score - I2M[i];
		}
		s_dm=del_im1[jm1] + score - D2M[i];

		//========= 2. Sum the likelihoods for each target state. ===========
		// eventually create an array of exp(x) values...
		// store the best on a heap???
		
		//========= 2a. sum for transitions into insert state. ===========
		if(s_mi < -200) s_mi=-200.0; 
		if(s_ii < -200) s_ii=-200.0; 
		ins_i[jj] = log(exp(s_mi) + exp(s_ii)); 
		
		//========= 2b. sum for transitions into delete state. ===========
		if(s_md < -200) s_md=-200.0; 
		if(s_dd < -200) s_dd=-200.0; 
		del_i[jj] = log(exp(s_md) + exp(s_dd)); 
		
		//========= 2c. sum for transitions into match state. ===========
		if(s_mm < -200) s_mm=-200.0; 
		if(s_dm < -200) s_dm=-200.0; 
		if(s_im < -200) s_im=-200.0; 
		mat_i[jj] = log(exp(s_dm) + exp(s_mm) + exp(s_im)); 
if(HG) { IncdHist(ins_i[jj], HG); IncdHist(del_i[jj], HG); IncdHist(mat_i[jj], HG); }
	   }
	} 
if(HG){ PutHist(stderr,60,HG); NilHist(HG); }

	//=================== 3. Sample a traceback based on likelihoods. ==============
	Int4	s,t,t0,J0,max_i,max_j,end,im1;
	char	state=0;
	max_i=hmm_len;		// global with respect to profile.

	//============ 3a. Sample the starting point for the traceback...
	// for now just take the most likely starting point...
	// Determine whether optimum alignment ends in a match or a deletion.
	for(score=(double)INT4_MIN, j=1; j<= sq_len; j++){
		if((dd=MAT[max_i][j]) > score){ score=dd; max_j=j; }
	} state='M';
	for(j=1; j<= sq_len; j++){
		if((dd=DEL[max_i][j]) > score){ score=dd; max_j=j; state='D'; }
	} // NOTE: insertions at end can't improve the score...

	NEW(operation,hmm_len+sq_len+3,char); operation[0]='E';
	for(J0=0,j=sq_len; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }

	//============ 3b. Sample back to the beginning point ...
	// Start at position i=max_i & j = max_j;
	double	total,rand_no;
	// fprintf(stderr,"max_j = %d; max_i=%d\n",max_j,max_i);
	for(i=max_i,jj=max_j; i > 0 && jj > 0; ){  // full alignment...Stop at end of profile
	    //========= 3b-i. Sample the next state based on the current state ====.
	    switch(state){
		case 'M':	  // Match state found...
		   J0++;
		   if(i==1) operation[J0] = 'M'; else operation[J0] = 'm';
		   i--; jj--; break;
		case 'D':   // Deletion: gap ('-') is in sequence; add X's for gaps
			J0++; 
			if(i==1) operation[J0] = 'D'; else operation[J0] = 'd';
			i--;  break;
		case 'I':  // Insertion: gap ('-') is in profile
			if(i < 1 || i==LenSMatrix(smx)) { J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			jj--; break;
		default: print_error("this should not happen"); break;
	    } if(i <= 0 || j <= 0) break;
	    jm1=jj-1; im1=i-1;
	    mat_i=MAT[i]; mat_im1=MAT[im1];
	    ins_i=INS[i]; ins_im1=INS[im1];
	    del_i=DEL[i]; del_im1=DEL[im1];
	    rand_no = (double) Random()/(double) RANDOM_MAX;
	    switch(state){		// How did we get into this [i][jj] state?
		case 'M':
		   // score = (double) ValSMatrix(i,seq2[jj],smx);  // added to all, so can ignore.
	    	   s_dm=exp(del_i[jj]); s_im=exp(ins_i[jj]);
		   s_mm=exp(mat_i[jj]);
		   total=s_dm + s_im + s_mm;
		   dd=s_mm/total;
		   if((rand_no -= dd) <= 0.0){ state='M';  break; } 
		   dd=s_dm/total;
		   if((rand_no -= dd) <= 0.0){ state='D';  break; } 
		   dd=s_im/total;
	// WARNING: This is failing...occasionally...
		   assert((rand_no -= dd) <= 0.0); 
		   state='I';  break; 
		case 'D':	  
		   s_md=exp(mat_i[jj]); s_dd=exp(del_i[jj]);
		   total=s_dd + s_md;
		   dd=s_dd/total;
		   if((rand_no -= dd) <= 0.0){ state='D';  break; } 
		   dd=s_md/total;
		   assert((rand_no -= dd) <= 0.0); state='M';  break;  
		case 'I':	  
		   s_mi=exp(mat_i[jj]); s_ii=exp(ins_i[jj]);
		   total=s_ii + s_mi;
		   dd=s_ii/total;
		   if((rand_no -= dd) <= 0.0){ state='I';  break; } 
		   dd=s_mi/total;
		   assert((rand_no -= dd) <= 0.0); state='M';  break;  
		   break;
	    }
	}
	if(jj <= 0){
		while(i > 1){ J0++; operation[J0]='d'; i--; }
		if(i > 0){ J0++; operation[J0]='D'; i--; }
	}
	// fprintf(stderr,"jj = %d; i=%d\n",jj,i);
	*start=jj+1; J0++; operation[J0]='E';
	*alnscore = (Int4) ceil(((double)pernats*score)-0.5); *J = J0;

	//=============== 4. free allocated memory ======================
	for(i=0; i<=hmm_len; i++) {free(MAT[i]); free(DEL[i]);free(INS[i]);}
	free(MAT); free(DEL); free(INS); 
#if 1
	free(M2M); free(M2D); free(M2I); 
	free(I2M); free(I2I);
	free(D2D); free(D2M); free(S2M); 
#endif
        return operation;
}


