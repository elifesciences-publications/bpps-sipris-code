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

#include "jlh_typ.h"

#define debug_JLH_penalty 0

double	jlh_typ::IndelPenalty(FILE *fp,char mode,char *Wt)
// Test Jun Liu's HMM statistical model for GISMO...
// Wt == array of integer weights for each sequence; used merely as a BooLean indicator here.
{
	cma_typ cma=CMA;
	//========== 1. Compute effective # of observed HMM transitions. ==============
#if debug_JLH_penalty 
	// fprintf(stderr,"Nmm=%d; Nmi=%d; Nmd=%d; Nii=%d; Nim=%d\n",Nmm,Nmi,Nmd,Nii,Nim);
	// fprintf(stderr,"  Ndd=%d; Ndm=%d; Ndi=%d; Nid=%d; Nsd=%d; Nsm=%d\n\n",Ndd,Ndm,Ndi,Nid,Nsd,Nsm);
#endif
	double	D,Dmm,Dmd,Dmi,Dm, Dii,Ddd,Dim,Ddm, Did,Ddi,Dsd,Dsm;
	double	PriorDownWt=10;	// downweight placed on priors... 10 == 1/ten weight of data.
	double	DataDownWt=1;	// downweight placed on data... 1 == no downweight...probably best.
	double	Adjust=DataDownWt*(double)WtFact;
	if(Wt == 0) Adjust=DataDownWt; 
	Dmm=(double)Nmm/(double)Adjust; Dmd=(double)Nmd/(double)Adjust; Dmi=(double)Nmi/(double)Adjust;
	Ddd=(double)Ndd/(double)Adjust; Ddi=(double)Ndi/(double)Adjust; Ddm=(double)Ndm/(double)Adjust;
	Dim=(double)Nim/(double)Adjust; Dii=(double)Nii/(double)Adjust; Did=(double)Nid/(double)Adjust;
	Dsd=(double)Nsd/(double)Adjust; Dsm=(double)Nsm/(double)Adjust;
#if debug_JLH_penalty 
	fprintf(stderr," Nmm=%.1f; Nmi=%.1f; Nmd=%.1f\n Nii=%.1f; Nim=%.1f\n",Dmm,Dmi,Dmd,Dii,Dim);
	fprintf(stderr," Ndd=%.1f; Ndm=%.1f;\n (Ndi=%.1f; Nid=%.1f; Nsd=%.1f; Nsm=%.1f)\n",
					Ddd,Ddm,Ddi,Did,Dsd,Dsm);
	fprintf(stderr,"  aa_per_ins_o=%d;",aa_per_ins_o);
	fprintf(stderr,"  exp_ins_ext=%d\n",exp_ins_ext);
	fprintf(stderr,"  aa_per_del_o=%d;",aa_per_del_o);
	fprintf(stderr,"  exp_del_ext=%d\n",exp_del_ext);
	fprintf(stderr,"  mode='%c'\n",mode);
#endif

	//========== 2. Compute Observed HMM transitions. ==============
	if(Wt == 0) Adjust=PriorDownWt; else Adjust=PriorDownWt*(double)WtFact; 
	UInt8	ae,be,ao,bo,n_all,n_m,n_i,n_d;
	D = (double)(Nmm+Nmi+Nmd+Nii+Nim+Ndd+Ndm+Ndi+Nid+Nsd+Nsm)/(double)Adjust;
	n_all = (UInt8) floor(0.5+D);	
	ao=(UInt8) floor(0.5+(double)n_all/(double)aa_per_ins_o);  // expected total # of insert openings.
	if(ao < 1) ao=1;
	bo=(UInt8) floor(0.5+(double)n_all/(double)aa_per_del_o);  // expected total # of deletion openings.
	if(bo < 1) bo=1;

	ae=ao*exp_ins_ext;	// expected total # of insertion extensions.
	be=bo*exp_del_ext;	// expected total # of deletion extensions.
	n_i=ao + ae;		// total number of expected insertions.
	n_d=bo + be;		// expected total # of expected deletions.
	n_m = n_all-n_d-n_i;	// expected total # of matching residues.

	// assert(n_all < LONG_MAX/4);	// make sure don't get overflow.
// std::cerr << "DEBUG 1" << std::endl;
	//   RoundNearest(xx,WtFact);
#if debug_JLH_penalty
	fprintf(stderr," n_mm=%d; n_mi =%d; n_md=%d\n",n_m-ao-bo,ao,bo);
	fprintf(stderr," n_ii=%d; n_im=%d\n",ae,n_i-ae);
	fprintf(stderr," n_dd=%d; n_dm=%d\n",be,n_d-be);
	fprintf(stderr,"Prior: n_m=%d; n_all=%d; n_i=%d; n_d=%d\n",n_m,n_all,n_i,n_d);
#endif

	//========== 3. Compute posterior indel penalties. ==============
	double Bo,Be,Ao,Ae;
	double	penalty=0.0;
	double d_m=(double)n_m,d_i=(double)n_i,d_d=(double)n_d;
	Dm=Dmm+Dmi+Dmd;
#if 0	// ignore data!!
	Ao=(double)ao; Ae=(double)ae; Bo=(double)bo; Be=(double)be;
#elif 1	// MLE of prior parameters.
	Ao=(double)ao; Ae=(double)ae; Bo=(double)bo; Be=(double)be;
#if debug_JLH_penalty
	fprintf(stderr,"\nao = %g; ae = %g; ao = %g; ae = %g\n",Ao,Ae,Bo,Be);
#endif
	Bo=(Dmd+Bo)/(Dmm+Dmi+Dmd+d_m);
	Be=(Ddd+Be)/(Ddm+Ddd+d_d);
	Ao=(Dmi+Ao)/(Dmm+Dmi+Dmd+d_m);
	Ae=(Dii+Ae)/(Dim+Dii+d_i);
#if debug_JLH_penalty
	fprintf(stderr,"Ao(m->i) = %g; Ae(i->i) = %g; Bo(m->d) = %g; Be(d->d) = %g\n",Ao,Ae,Bo,Be);
	fprintf(stderr,"  1-Ao-Bo(m->m) = %g;\n\n",1-Ao-Bo);
#endif
#else
	Bo=(double)(Dmd+bo)/(double)(Dmm+Dmi+Dmd+n_m);
	Be=(double)(Ddd+be)/(double)(Ddm+Ddd+n_d);
	Ao=(double)(Dmi+ao)/(double)(Dmm+Dmi+Dmd+n_m);
	Ae=(double)(Dii+ae)/(double)(Dim+Dii+n_i);
#endif
	long double Penalty=0.0;
	// numerator... from match states.
	Penalty += lgammal(Dmi+Ao) + lgammal(Dmd+Bo);
	Penalty += lgammal(Dmm+d_m-Ao-Ao) + lgammal(d_m);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= lgammal((double)(Dm+d_m)) + lgammal(Ao) + lgammal(Bo);
	Penalty -= lgammal(d_m-Ao-Ao);

	// numerator...from insert states.
	Penalty += lgammal(Dii+Ae) + lgammal(Dim+d_i-Ae) + lgammal(d_i);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= lgammal(Dim+Dii+d_i) + lgammal(Ae) + lgammal(d_i-Ae);

	// numerator...from delete states.
	Penalty += lgammal(Ddd+Be) + lgammal(Ddm+d_d-Be) + lgammal(d_d);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= lgammal(Ddd+Ddm+d_d) + lgammal(Be) + lgammal(d_d-Be);
	penalty = (double) Penalty;

	// std::cerr << "DEBUG 2" << std::endl;
#if debug_JLH_penalty
	std::cerr << "gap penalty = " << penalty << std::endl;
#endif
	// std::cerr << "DEBUG 3" << std::endl;
	//========== 4. Compute indel penalties. ==============
	if(mode != 'S' && mode != 's'){ 	// MLE of parameters
	    //========== 4a. Use average penalties (don't sample). ==============
#if 1	// MLEs from paper:	
	    Bo=(double)(Dmd+bo)/(double)(Dmm+Dmi+Dmd+n_m);
	    Be=(double)(Ddd+be)/(double)(Ddm+Ddd+n_d);
	    Ao=(double)(Dmi+ao)/(double)(Dmm+Dmi+Dmd+n_m);
	    Ae=(double)(Dii+ae)/(double)(Dim+Dii+n_i);
#else 	// pseudocounts only...
	    Ao=(double)(ao)/(double)(n_m); Ae=(double)(ae)/(double)(n_i);
	    Bo=(double)(bo)/(double)(n_m); Be=(double)(be)/(double)(n_d);
#endif
#if debug_JLH_penalty
	    fprintf(stderr,"\nPrior only: Ao = %g; Ae = %g; Bo = %g; Be = %g\n", 
	    	   (double)(ao)/(double)(n_m),(double)(ae)/(double)(n_i),
	    	   (double)(bo)/(double)(n_m),(double)(be)/(double)(n_d));
	    fprintf(stderr,"Posterior: Ao = %g; Ae = %g; Bo = %g; Be = %g\n",Ao,Ae,Bo,Be);
	    fprintf(stderr,"In thousandths nats: Ao = %g; Ae = %g; Bo = %g; Be = %g\n\n",
			pernats*log(Ao), pernats*log(Ae), pernats*log(Bo), pernats*log(Be));
	    fflush(stderr);
#endif
	} else {	
	   //========== 5. Sample indel penalties. ==============
	   static Int4 Seed=0;
	   if(Seed==0) {
             Seed=-Random();  // need to initialize with a negative number.
             fprintf(stderr,"INITIALIZING JLH SEED (%d)\n",Seed);
	     fflush(stderr);
           } // Parameters based both on observed and prior counts.
#if 0
	// Add temperature factor!!!
	
#endif
	   UInt8 alpha,beta;
	   double  All,BetaDensity=1000.0;  // BetaDensity corresponds to temperature.
	   for(Int4 iter=1; iter < 10; iter++)
	   {
		// temperature == BetaDensity.
	     All=(double)(Dmd+Dmm+Dmi+n_m)/BetaDensity; 
	     alpha=RoundNearest(Dmd+bo,All);
	     beta=RoundNearest(Dmm+Dmi+n_m-bo,All);
	     if(alpha < 1) alpha=1; if(beta < 1) beta=1;
	     Bo=betadev(alpha,beta,&Seed);

	     All=(double)(Ddd+Ddm+n_d)/BetaDensity; 
	     alpha=RoundNearest(Ddd+be,All);
	     beta=RoundNearest(Ddm+n_d-be,All);
	     if(alpha < 1) alpha=1; if(beta < 1) beta=1;
	     Be=betadev(alpha,beta,&Seed);

	     All=(double)(Dmi+Dmm+n_m-bo)/BetaDensity; 
	     alpha=RoundNearest(Dmi+ao,All);
	     beta=RoundNearest(Dmm+n_m-ao-bo,All);
	     if(alpha < 1) alpha=1; if(beta < 1) beta=1;
	     Ao=(1-Bo)*betadev(alpha,beta,&Seed);

	     All=(double)(Dii+Dim+n_i)/BetaDensity; 
	     alpha=RoundNearest(Dii+ae,All);
	     beta=RoundNearest(Dim+n_i-ae,All);
	     if(alpha < 1) alpha=1; if(beta < 1) beta=1;
	     Ae=betadev(alpha,beta,&Seed);
	     break;	// one iteration only...
	     fprintf(stderr,"Ao = %g; Ae = %g; Bo = %g; Be = %g\n",Ao,Ae,Bo,Be);
	   }
	} // fprintf(stderr,"\nAo = %g; Ae = %g; Bo = %g; Be = %g\n\n",Ao,Ae,Bo,Be);
	m2m = -(Int4)floor((pernats*log(1.0 -Ao-Bo))+0.5);
	m2i = -(Int4)floor((pernats*log(Ao))+0.5);
	m2d = -(Int4)floor((pernats*log(Bo))+0.5);
	i2i = -(Int4)floor((pernats*log(Ae))+0.5);
	i2m = -(Int4)floor((pernats*log(1.0 -Ae))+0.5);
	d2d = -(Int4)floor((pernats*log(Be))+0.5);
	d2m = -(Int4)floor((pernats*log(1.0-Be))+0.5);
	s2m = -(Int4)floor((pernats*log(1.0 -Bo))+0.5);
#if debug_JLH_penalty
	fprintf(stderr,"m2m=%d; m2i=%d; i2i=%d; i2m=%d; m2d=%d; d2d=%d; d2m=%d; s2m=%d.\n",
                        m2m,m2i,i2i,i2m,m2d,d2d,d2m,s2m);
// exit(1);
#endif
	ReCalc=FALSE;
	return penalty;
}

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

char	*jlh_typ::GapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score)
// return the trace, the start position of the trace, and the trace length.
{
	Int4    i,j,newlen,lenTrace;
	char	*operation,*tmp;

	tmp=this->gap_align(LenSeq(E),SeqPtr(E),nmod,M,&lenTrace,score,start);
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
char	*jlh_typ::gap_align(Int4 n2,unsigned char *seq2,Int4 nmod,smx_typ *M,
			Int4 *J,Int4 *alnscore,Int4 *start)
{
	Int4	m,i,j,k,n1,score,g,gINS;
	Int4	**MAT,**T,**DEL,**INS;
	Int4	s_mm,s_mi,s_md,s_dd,s_dm,s_ii,s_im;
	Int4	s_di,s_id;	// Jun's HMM.
	char	*operation;

	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
#if debugger_on	// debug using matrix output.
                NEW(MAT[i],n2+3,Int4); NEW(T[i],n2+3,Int4);
                NEW(DEL[i],n2+3,Int4); NEW(INS[i],n2+3,Int4);
#else 	
		MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
#endif

	}
	// Make GLOBAL ALIGNMENT with respect to profile,
/********************************************************************************
  j == sequence. i == HMM.
       DEL[i][j]               INS[i][j]              MAT[i][j]
       j\i: 0 : 1 : 2 : 3 :    j\i: 0 : 1 : 2 : 3 :   j\i: 0 : 1 : 2 : 3 :
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        0 | 0 | md|+dd|+dd|     0 |inf|inf|inf|inf|    0 | 0 |inf|inf|inf|
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        1 |inf|   |   | M |     1 | 0 |   |   |   |    1 | 0 | M |   |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        2 |inf|   |   |   | M   2 | 0 |   |   |   |    2 | 0 |   | M |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        3 |inf|   |   |   |     3 | 0 | M |   |   |    3 | 0 |   |   | M |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
 mtf A: MWYYWCFWFCFWWF	           MWYYWCFWFCFWWF        MWYYWCFWFCFWWF
          ::::::::::::             ::::::::::::::        ::::::::::::::
        --YYWCFWFCFWWF         (yf)MWYYWCFWFCFWWF        MWYYWCFWFCFWWF

/********************************************************************************/
	MAT[0][0] = 0; 		// corresponds to S0 state.
	DEL[0][0] = d2d-m2d;  // --> DEL[1][0] = m2d. note m2d == s2d...
	INS[0][0] = 0;
	T[0][0] = 0;	
	for(i=1; i<= n1; i++) {  // j==0 implies deletion of model position route.
	    DEL[i][0] = DEL[i-1][0] - d2d;	// penalty for deletions at end (-m2d too?)
	    INS[i][0] = -infinity;	// disallow internal inserts from start of sequence.
					// requires d2i transitions.
	    MAT[i][0] = 0;
	    MAT[i][0] = -infinity;  
	    T[i][0] = T[i-1][0] + 1; // deletion traceback route 
	}
	for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
	    	MAT[0][j] = 0;  
		DEL[0][j] = -infinity;
		INS[0][j] = 0;	// no penalty for inserts before motif...
		T[0][j]= T[0][j-1] -1; // insert traceback route 
		// T[0][j]= -1; // traceback insert route 
	} DEL[0][0] = 0;  
	Int4 *gdel; NEW(gdel,n2+3,Int4);
	// 2. Dynamic programming step. 
	for(k=m=i=1; i<= n1; i++) {
	   gINS=0;
	   Int4 *mat_im1=MAT[i-1]; 
	   Int4 *mat_i=MAT[i];
	   Int4 *ins_i=INS[i];
	   Int4 *ins_im1=INS[i-1];
	   Int4 *del_i=DEL[i];
	   Int4 *del_im1=DEL[i-1];
	   register Int4 jj,jm1;
	   smx_typ  smx=M[m];
	   for(jm1=0,jj=1; jj<= n2; jj++,jm1++) {  
		// 1. Choose optimum among m->i vs i->i transitions.
	   	Int4	t,s_max;
	        if(k==LenSMatrix(smx)){	// NO PENALTIES for inserts at ends!
		  s_mi = mat_i[jm1];	// insertion opening score
		  s_ii = ins_i[jm1];	// insertion extension score
	        } else {
                  s_mi=mat_i[jm1]-m2i;	// insertion opening score
		  s_ii=ins_i[jm1]-i2i;	// insertion extension score
		}
                if(s_mi > s_ii){ t=gINS=-1; ins_i[jj] = s_mi; s_max=s_mi; }
		else { gINS--; t=gINS; ins_i[jj] = s_ii; s_max=s_ii; }
		// else { t=-1; ins_i[jj] = s_ii; s_max=s_ii; }
#if use_juns_hmm		// Jun's HMM...
                if(k==LenSMatrix(smx)){
		   s_di = del_i[jm1]-d2m;	// no penalty: d2s = 0;
		   if(s_di > s_mi && s_di > s_ii){
		     gINS=-1; ins_i[jj]=s_di; t=-1; s_max=s_di; 
		   }
		}
#endif
		// 2. Choose optimum among m->d vs d->d transitions.
		s_md=mat_im1[jj]-m2d;	// note m2d == s2d 
		s_dd=del_im1[jj]-d2d;	// deletion extension score
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
		   s_id = ins_im1[jj]-m2d;	// - s2d?
		   if(s_id > s_md && s_id > s_dd){
			gdel[jj]=1; del_i[jj]=s_id;
			if(s_id > s_max){ s_max=s_id; t=1; }
		   }
		}
#endif

		// 3. Choose optimum among m->m, i->m or d->m transitions.
		score = ValSMatrix(k,seq2[jj],smx);  // match score at i,j cell.
	        if(k==1) {	// s2m PENALTY for match at start of motif!!!!
		  s_mm=mat_im1[jm1] + score - s2m;	// should be m2s???
		  s_im=ins_im1[jm1] + score;	// no s->m penalty.
	        } else {
		  s_mm=mat_im1[jm1] + score - m2m;
		  s_im=ins_im1[jm1]+score-i2m;
		}
		s_dm=del_im1[jm1]+score-d2m;

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
#if 0		// debug...
	    if(m==2 && jj > 290 && jj < 316){
		Int4 tm1;
		if(t ==0){  // match: tm1 can be anything...
		} else if(t < 0){	// insertion...tm1 is I or M.
		   tm1=T[i][jj + t];	// points back one insert...
		   if(tm1 > 0){ fprintf(stderr,"s_max = %d\n",s_max); }
		   assert(tm1 <= 0);	// insertion or match...
		} else {	// t > 0 --> deletion...tm1 is D or M.
		   tm1=T[i-t][jj];	// point back t deletions...
		   if(tm1 < 0){ fprintf(stderr,"s_max = %d\n",s_max); }
		   assert(tm1 >= 0);	// deletion or match
		}
	    }
#endif
	   } k++;
	   if(k > LenSMatrix(M[m])) { k = 1; m++; }
	} free(gdel);
#if debugger_on   //******************** debug using matrix output. *********************
        /** get concensus sequence for smatrix **/
	a_type A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	Int4 startM,startS=0,endS=n2,motif=2;
	assert(motif <= nmod);
	unsigned char *seq1;
	startS=290; endS=316;
	// if(m==1){ i=0; startS=0; } else i=1;
        NEW(seq1, n1+3, unsigned char);

#if 1
	{ for(Int4 x=0, m = 1; m <= nmod; m++){ 
        	MinSegSMatrix(seq1+x,M[m]); // NEED ONLY FOR OUTPUT.
		x += LenSMatrix(M[m]);
	  } e_type E=MkSeq("min_seq for debugging",n1,seq1);
	  PutSeq(stderr,E,A);
	}
	{ for(Int4 x=0, m = 1; m <= nmod; m++){ 
        	MaxSegSMatrix(seq1+x,M[m]); // NEED ONLY FOR OUTPUT.
		x += LenSMatrix(M[m]);
	  } e_type E=MkSeq("max_seq for debugging",n1,seq1);
	  PutSeq(stderr,E,A);
	}
#endif
  startM=0;
  i=10;	// must be that 0 < i < lenM;
  { Int4 s;
        for(s=m=1; m <= nmod; m++){
            for(Int4 j=1; j<= LenSMatrix(M[m]); j++){ 
		if(m==motif && j==i){ startM=s; break; }  
		else s++;
	    } if(m==motif) break;
        }
  }
	m=motif;
        MaxSegSMatrix(seq1+(startM-i),M[motif]); // NEED ONLY FOR OUTPUT.
	PutSMatrix(stderr, M[motif]);
        fprintf(stderr,"DEL[i][j]:\n");
        this->put_dp_matrixSW(DEL, T, n2, M[motif],seq1,seq2,A,startM,i,startS,endS);
        fprintf(stderr,"INS[i][j]:\n");
        this->put_dp_matrixSW(INS, T, n2, M[motif],seq1,seq2,A,startM,i,startS,endS);
        fprintf(stderr,"MAT[i][j]:\n");
        this->put_dp_matrixSW(MAT, T, n2, M[motif],seq1,seq2,A,startM,i,startS,endS);
#if 1
        fprintf(stderr,"SMX[i][j]:\n");
        this->put_dp_matrixSMX(n2, M[motif],seq1,seq2,A,startM,i,startS,endS);
#endif
        fprintf(stderr,"T[i][j]:\n");
        this->put_dp_matrixSW(T, T, n2, M[motif],seq1,seq2,A,startM,i,startS,endS);
#endif  //*****************************************************************************
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
            G   P   T   M           G   H   W   A          P   H   W   T
          :   :   : i :   :       :   :   :i-1: i :      :i-3:   :   : i :
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
     W    |   |   |   |   |  W j-3|   |   | s |   | W    |   |   |   |   |
       ---+---+---+---+---+    ---+---+---+-^-+---+   ---+---+---+---+---+
     P    |   | s_|   |   |  P    |   |   | | |   | P  j | s<----------D |
       ---+---+---\---+---+    ---+---+---+-|-+---+   ---+---+---+---+---+
     S  j |   |   |`M |   |  S    |   |   | I |   | S    |   |   |   |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
     A    |   |   |   |   |  A  j |   |   |   |   | A    |   |   |   |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
          T[i][j] == 0            T[i][j] == -2          T[i][j] == 3
             s --> M                s --> I                 s --> D 
              
    model[i]   PT                     W--A                  PhwT
	       :.                     :  :                  :  :
    seq[j]     PS                     WpsA                  P--S

#endif	//********************************************************************
	Int4	s,t,t0,J0,max_i,max_j,end;
	max_i=n1;		// global with respect to profile.
	for(score=INT4_MIN, j=1; j<= n2; j++){
		if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
	}
#if 1	// Determine whether optimum alignment ends in a match or a deletion.
	for(j=1; j<= n2; j++){
		if((s=DEL[max_i][j]) > score){ score=s; max_j=j; }
	}
#endif	// NOTE: insertions at end can't improve the score...
	NEW(operation,n1+n2+3,char); operation[0]='E';
	m=nmod; k=LenSMatrix(M[m]);
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
			if(m < 1 || k==LenSMatrix(M[m]))
				{ J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			j--; break;
		default: print_error("this should not happen"); break;
	     }
	     if(k==0){ m--; if(m > 0) k=LenSMatrix(M[m]); }
	  } while(t0 != 0);
	}
	*start=j+1; J0++; operation[J0]='E';
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(T); free(DEL); free(INS); 
	*alnscore = score; *J = J0;
        return operation;
}


