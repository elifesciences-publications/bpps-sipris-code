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

void	ndl_typ::PutPenalties(FILE *fp)
{
	cma_typ cma=CMA;
	Int4	   b,i,j,len=TotalLenCMSA(CMA);
        h_type  m2iHG = Histogram("m2i penalties",1000,10000,250);
        h_type  m2dHG = Histogram("m2d penalties",1000,10000,250);
        h_type  i2iHG = Histogram("i2i penalties",0,1000,50);
        h_type  d2dHG = Histogram("d2d penalties",0,1000,50);
   	for(j=0,b=1; b <= nBlksCMSA(cma); b++){
   	   for(i=1; i <= LengthCMSA(b,cma); i++){
		j++;
		IncdHist(m2i[j], m2iHG); IncdHist(m2d[j], m2dHG);
		IncdHist(i2i[j], i2iHG); IncdHist(d2d[j], d2dHG);
	   }
	}
	PutHist(stderr,60,m2iHG); NilHist(m2iHG);
	PutHist(stderr,60,m2dHG); NilHist(m2dHG);
	PutHist(stderr,60,i2iHG); NilHist(i2iHG);
	PutHist(stderr,60,d2dHG); NilHist(d2dHG);
}

void	ndl_typ::PutPenalties(FILE *fp, Int4 blk, Int4 I)
{
	cma_typ cma=CMA;
	Int4	b,i,j;
   	for(j=0,b=1; b <= nBlksCMSA(cma); b++){
   	   for(i=1; i <= LengthCMSA(b,cma); i++){
		j++;
		if(b == blk && i == I){
		  fprintf(fp," *** %d: m2m=%d; m2i=%d; i2i=%d; i2m=%d; m2d=%d; d2d=%d; d2m=%d; s2m=%d ***\n",
                        	j,m2m[j],m2i[j],i2i[j],i2m[j],m2d[j],d2d[j],d2m[j],s2m[j]);
		  return;
		}
	   }
	} 
}

#define debug_NDL_penalty 0

double	ndl_typ::IndelPenalty(FILE *fp,char mode,char *Wt)
// Generalization of Jun Liu's HMM statistical model for GISMO...
// Wt == array of integer weights for each sequence; used merely as a BooLean indicator here.
{
   cma_typ cma=CMA;
   double  D,penalty=0.0;
   Int4	   i,j,len=TotalLenCMSA(CMA);
   assert(Wt != 0);
#if debug_NDL_penalty 
   fprintf(stderr,"  aa_per_ins_o=%d; exp_ins_ext=%d\n",aa_per_ins_o,exp_ins_ext);
   fprintf(stderr,"  aa_per_del_o=%d; exp_del_ext=%d\n",aa_per_del_o,exp_del_ext);
   fprintf(stderr,"  mode='%c'\n",mode);
#endif

   assert(exp_del_ext > 0);
   assert(exp_ins_ext > 0);
   for(i=1; i <= len; i++){
	//========== 1. Compute the total Observed HMM transitions times wt_factor. ==============
	UInt4	ae,be,ao,bo,n_all,n_m,n_i,n_d;
	D = Nmm[i]+Nmi[i]+Nmd[i]+Nii[i]+Nim[i]+Ndd[i]+Ndm[i]+Ndi[i]+Nid[i]+Nsd[i]+Nsm[i];
	// D = Nmm[i]+Nmd[i]+Ndd[i]+Ndm[i]+Nim[i];
	n_all = (UInt4) floor(0.5+D*prior_wt);	// true D = D/(double)WtFact;

	// obtain the pseudo counts for computations (i.e., multiplied x WtFact).
	ao=(UInt4) floor(0.5+(double)n_all/(double)aa_per_ins_o);  // expected total # of insert openings.
	if(ao < 1) ao=1;
	bo=(UInt4) floor(0.5+(double)n_all/(double)aa_per_del_o);  // expected total # of deletion openings.
	if(bo < 1) bo=1;
	ae = ao*exp_ins_ext;	// expected total # of insertion extensions.
	be = bo*exp_del_ext;	// expected total # of deletion extensions.
	n_i = ao + ae;		// total number of expected insertions.
	n_d = bo + be;		// expected total # of expected deletions.
	n_m = n_all-n_d-n_i;	// expected total # of matching residues.

	// assert(n_all < LONG_MAX/2);	// make sure don't get overflow.
	//   RoundNearest(xx,WtFact);
#if debug_NDL_penalty
	fprintf(stderr,"%d.Data: N_mm=%d; N_mi =%d; N_md=%d;",i,Nmm[i],Nmi[i],Nmd[i]);
	 fprintf(stderr," N_ii=%d; N_im=%d;",Nii[i],Nim[i]);
	 fprintf(stderr," N_dd=%d; N_dm=%d\n",Ndd[i],Ndm[i]);
	fprintf(stderr,"%d.Prior: n_mm=%d; n_mi =%d; n_md=%d;",i,n_m-ao-bo,ao,bo);
	 fprintf(stderr," n_ii=%d; n_im=%d;",ae,n_i-ae);
	 fprintf(stderr," n_dd=%d; n_dm=%d\n",be,n_d-be);
	fprintf(stderr,"   n_m=%d; n_all=%d; n_i=%d; n_d=%d\n",n_m,n_all,n_i,n_d);
#endif
	//========== 2. Compute posterior indel probability. ==============
	double Penalty=0.0;
	// numerator... from match states.
	UInt4 Nm=Nmi[i]+Nmm[i]+Nmd[i];
	Penalty += LGM->LGamma(Nmi[i]+ao) + LGM->LGamma(Nmd[i]+bo);
	Penalty += LGM->LGamma(Nmm[i]+n_m-ao-bo) + LGM->LGamma(n_m);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= LGM->LGamma(Nm+n_m) + LGM->LGamma(ao) + LGM->LGamma(bo);
	Penalty -= LGM->LGamma(n_m-ao-bo);

	// numerator...from insert states.
	Penalty += LGM->LGamma(Nii[i]+ae) + LGM->LGamma(Nim[i]+n_i-ae) + LGM->LGamma(n_i);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= LGM->LGamma(Nim[i]+Nii[i]+n_i) + LGM->LGamma(ae) + LGM->LGamma(n_i-ae);

	// numerator...from delete states.
	Penalty += LGM->LGamma(Ndd[i]+be) + LGM->LGamma(Ndm[i]+n_d-be) + LGM->LGamma(n_d);
	assert(!isinf(Penalty));
	// denominator...
	Penalty -= LGM->LGamma(Ndd[i]+Ndm[i]+n_d) + LGM->LGamma(be) + LGM->LGamma(n_d-be);
	penalty += (double) Penalty;

	// std::cerr << "DEBUG 2" << std::endl;
#if debug_NDL_penalty
	std::cerr << "  gap penalty = " << penalty << std::endl;
#endif
	//========== 3. Compute indel penalties. ==============
	// Call as a routine?: void ndl_typ::CalcIndelTransProb(FILE *fp,char mode,char *Wt)
	double Bo,Be,Ao,Ae;	// Dirichlet parameters alpha_o, alpha_e, beta_o, beta_e
	if(TRUE || mode != 'S' && mode != 's'){ 	// MLE of parameters
	    //========== 4a. Use average penalties (don't sample). ==============
#if 1	// MLEs from paper:	
	    Ao=(double)(Nmi[i]+ao)/(double)(Nmm[i]+Nmi[i]+Nmd[i]+n_m);
	    Ae=(double)(Nii[i]+ae)/(double)(Nim[i]+Nii[i]+n_i);
	    Bo=(double)(Nmd[i]+bo)/(double)(Nmm[i]+Nmi[i]+Nmd[i]+n_m);
	    Be=(double)(Ndd[i]+be)/(double)(Ndm[i]+Ndd[i]+n_d);
#else 	// pseudocounts only...
	    Ao=(double)(ao)/(double)(n_m); Ae=(double)(ae)/(double)(n_i);
	    Bo=(double)(bo)/(double)(n_m); Be=(double)(be)/(double)(n_d);
#endif
#if debug_NDL_penalty
	    fprintf(stderr,"Prior only: Ao = %g; Ae = %g; Bo = %g; Be = %g\n", 
	    	   (double)(ao)/(double)(n_m),(double)(ae)/(double)(n_i),
	    	   (double)(bo)/(double)(n_m),(double)(be)/(double)(n_d));
	    fprintf(stderr,"Posterior: Ao = %g; Ae = %g; Bo = %g; Be = %g\n",Ao,Ae,Bo,Be);
	    fprintf(stderr,"In thousandths nats: Ao = %g; Ae = %g; Bo = %g; Be = %g\n",
			pernats*log(Ao), pernats*log(Ae), pernats*log(Bo), pernats*log(Be));
	    fflush(stderr);
#endif
	} else {	
	 //========== 5. Sample indel penalties. ==============
	 double	Dmm,Dmd,Dmi,Dm, Dii,Ddd,Dim,Ddm, Did,Ddi,Dsd,Dsm,wt=(double)WtFact;
	 Dmm=(double)Nmm[i]/wt; Dmd=(double)Nmd[i]/wt; Dmi=(double)Nmi[i]/wt;
	 Ddd=(double)Ndd[i]/wt; Ddi=(double)Ndi[i]/wt; Ddm=(double)Ndm[i]/wt;
	 Dim=(double)Nim[i]/wt; Dii=(double)Nii[i]/wt; Did=(double)Nid[i]/wt;
	 Dsd=(double)Nsd[i]/wt; Dsm=(double)Nsm[i]/wt;
	 Dm=Dmm+Dmi+Dmd;
#if debug_NDL_penalty 
	 fprintf(stderr,"****************  column =%d ****************\n",i);
	 fprintf(stderr," Nmm=%.1f; Nmi=%.1f; Nmd=%.1f\n Nii=%.1f; Nim=%.1f\n",Dmm,Dmi,Dmd,Dii,Dim);
	 fprintf(stderr," Ndd=%.1f; Ndm=%.1f;\n (Ndi=%.1f; Nid=%.1f; Nsd=%.1f; Nsm=%.1f)\n",
					Ddd,Ddm,Ddi,Did,Dsd,Dsm);
#endif
	   static Int4 Seed=0;
	   if(Seed==0) {
             Seed=-Random();  // need to initialize with a negative number.
             fprintf(stderr,"INITIALIZING NDL SEED (%d)\n",Seed);
	     fflush(stderr);
           } // Parameters based both on observed and prior counts.
#if 0
	// Add temperature factor!!!
	
#endif
	   UInt4 alpha,beta;
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

	m2m[i] = -(Int4)floor((pernats*log(1.0 -Ao-Bo))+0.5);
	m2i[i] = -(Int4)floor((pernats*log(Ao))+0.5);
	m2d[i] = -(Int4)floor((pernats*log(Bo))+0.5);
	i2i[i] = -(Int4)floor((pernats*log(Ae))+0.5);
	i2m[i] = -(Int4)floor((pernats*log(1.0 -Ae))+0.5);
	d2d[i] = -(Int4)floor((pernats*log(Be))+0.5);
	d2m[i] = -(Int4)floor((pernats*log(1.0-Be))+0.5);
	s2m[i] = m2m[i];
#if 1	// compensate for low i2i transitions penalty; does screw up probability!!
	i2i[i]+=FloorI2I; 
	d2d[i]+=FloorD2D; 
#elif 1
	const Int4 FloorTP= 50; 	// put a floor on how low the penalty can go (0.050 nats).
	// const Int4 FloorTP= 25; 	// put a floor on how low the penalty can go (0.025 nats).
	// const Int4 FloorTP= 10; 	// put a floor on how low the penalty can go (0.010 nats).
	// const Int4 FloorTP= 5; 	// put a floor on how low the penalty can go (0.010 nats).
	// m2i[i]+=FloorTP; 
	// i2m[i]+=FloorTP; 
	i2i[i]+=FloorTP; 
	d2d[i]+=FloorTP; // m2d[i]+=FloorTP; 
#else
	const Int4 FloorTP= 50; 	// put a floor on how low the penalty can go (0.05 nats).
        m2m[i]+=FloorTP; m2i[i]+=FloorTP; m2d[i]+=FloorTP; i2i[i]+=FloorTP; 
	i2m[i]+=FloorTP; d2d[i]+=FloorTP; d2m[i]+=FloorTP; s2m[i]+=FloorTP;
#endif

#if debug_NDL_penalty
	fprintf(stderr,"m2m=%d; m2i=%d; i2i=%d; i2m=%d; m2d=%d; d2d=%d; d2m=%d; s2m=%d.\n\n",
                        m2m[i],m2i[i],i2i[i],i2m[i],m2d[i],d2d[i],d2m[i],s2m[i]);
// exit(1);
#endif
    } ReCalc=FALSE;
    // d2d[1]=m2d[1];	// disfavor deletions at the N-terminal end...
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

char	*ndl_typ::GapAlnTrace(e_type E,Int4 nmod, smx_typ *M, Int4 *start,Int4 *score,BooLean global)
// return the trace, the start position of the trace, and the trace length.
{
	Int4    i,j,newlen,lenTrace;
	char	*operation,*tmp;

#if 1
	tmp=this->gap_align(LenSeq(E),SeqPtr(E),nmod,M,&lenTrace,score,start,global);
#else	// DEBUG: sequences with multiple matches to model.
	if(LenSeq(E)==471){
		a_type AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
		PutSeq(stderr,E,AB);
#if 0	// mulitple domain search..
		// Need to modify (extend) m2m, m2i,...etc arrays.
		smx_typ MM[5]; MM[1]=M[1]; MM[2]=M[1];
		unsigned char *ptr=SeqPtr(E);
		Int4	lensq=LenSeq(E),os=OffSetSeq(E);
		tmp=this->gap_align(lensq,ptr,2,MM,&lenTrace,score,start,global);
		for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
		NEW(operation,newlen+3,char); operation[0]='E';
		for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
		operation[i]='E'; i++;
		operation[i]=0; free(tmp);
		put_seqaln_smatrixSW(stderr,operation,lensq,ptr,os,lenTrace,2,MM);

#else	// limit scope of search: region of seq searched.
		assert(nmod==1);
		unsigned char *ptr=SeqPtr(E);
		Int4	lensq=LenSeq(E),os=OffSetSeq(E);
		Int4	shift,Len=LenSMatrix(M[1]);
		double d=(double)lensq/(double)Len;
		if(d > 2.0){
		  shift=(2*Len/3) + 1; 
		  lensq -= shift;

		  tmp=this->gap_align(lensq,ptr,nmod,M,&lenTrace,score,start,global);
		  for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
		  NEW(operation,newlen+3,char); operation[0]='E';
		  for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
		  operation[i]='E'; i++;
		  operation[i]=0; free(tmp);

		  lenTrace=strlen(operation);
	   Int4 newpos[5],lengths[5]; newpos[1]=0; newpos[2]=0; lengths[1]=Len;
           gsq_typ *gsq = new gsq_typ[1]; // NEEDS TO BE AN ARRAY FOR gss_typ!!!
           gsq->initialize(10,10,operation,lenTrace,*start,E,newpos);
	   gsq->Put_cma_format(stderr,1,1,newpos,lengths,AB);
	   delete [] gsq;

		  put_seqaln_smatrixSW(stderr,operation,lensq,ptr+(*start-1),(*start-1),lenTrace,nmod,M);

		  ptr+= shift; os += shift;
		  tmp=this->gap_align(lensq,ptr,nmod,M,&lenTrace,score,start,global);
		  for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
		  NEW(operation,newlen+3,char); operation[0]='E';
		  for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
		  operation[i]='E'; i++;
		  operation[i]=0; free(tmp);
		  lenTrace=strlen(operation);
		  put_seqaln_smatrixSW(stderr,operation,lensq,ptr+(*start-1),shift+(*start-1),lenTrace,nmod,M);

		}
#endif
		fprintf(stderr,"\nNeed to turn off debugging section of ndl_typ::GapAlnTrace() in ndl_hmm.cc!\n\n");
		exit(1);
	} else tmp=this->gap_align(LenSeq(E),SeqPtr(E),nmod,M,&lenTrace,score,start,global);
#endif
	// put_seqaln_smatrixSW(stderr,tmp,LenSeq(E),SeqPtr(E),OffSetSeq(E),lenTrace,nmod,M);
#if 0
std::cerr << "tmp traceback"; std::cerr << std::endl;
std::cerr << tmp; std::cerr << std::endl;
#endif
	NEW(operation,lenTrace+3,char); operation[0]='E';
	for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
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

#define debugger_on 1		// these debugging routines are not working!!!
#define test_insert2delete 0	// allow insertion to deletion transitions; not fully implemented.
//======================= FAST gapped align for sampling ================
char	*ndl_typ::gap_align(Int4 sq_len,unsigned char *seq2,Int4 nmod,smx_typ *M,
			Int4 *J,Int4 *alnscore,Int4 *start,BooLean global)
{
	Int4	m,i,j,k,hmm_len,score,g,gINS;
	Int4	**MAT,**TB,**DEL,**INS;
	Int4	s_mm,s_mi,s_md,s_dd,s_dm,s_ii,s_im;
	Int4	s_di,s_id;	// Jun's HMM.
	char	*operation;

	/** get total length of profile **/
	for(hmm_len=0, m = 1; m <= nmod; m++){ hmm_len += LenSMatrix(M[m]); }

	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,hmm_len+3,Int4*); MEW(TB,hmm_len+3,Int4*);
	MEW(DEL,hmm_len+3,Int4*); MEW(INS,hmm_len+3,Int4*);
	for(i=0; i<= hmm_len; i++) { 
#if debugger_on	// debug using matrix output.
                NEW(MAT[i],sq_len+3,Int4); NEW(TB[i],sq_len+3,Int4);
                NEW(DEL[i],sq_len+3,Int4); NEW(INS[i],sq_len+3,Int4);
#else 	
		MEW(MAT[i],sq_len+3,Int4); MEW(TB[i],sq_len+3,Int4); 
		MEW(DEL[i],sq_len+3,Int4); MEW(INS[i],sq_len+3,Int4); 
#endif

	}
	// Make GLOBAL ALIGNMENT with respect to profile,
/********************************************************************************
  j == sequence. i == HMM.
       DEL[i][j]               INS[i][j]              MAT[i][j]
       j\i: 0 : 1 : 2 : 3 :    j\i: 0 : 1 : 2 : 3 :   j\i: 0 : 1 : 2 : 3 :
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        0 | 0 |-dd|-dd|   |     0 | 0 |inf|inf|inf|    0 | 0 |inf|inf|inf|
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        1 | 0 |   |   | M |     1 | 0 |   |   |   |    1 | 0 | M |   |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        2 | 0 |   |   |   |     2 | 0 |   |   |   |    2 | 0 |   | M |   |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
        3 | 0 |   |   |   |     3 |   | M |   |   |    3 | 0 |   |   | M |
       ---+---+---+---+---+    ---+---+---+---+---+   ---+---+---+---+---+
 HMM=i: MWYYWCFWFCFWWF             MWYYWCFWFCFWWF        MWYYWCFWFCFWWF
          ::::::::::::             ::::::::::::::        ::::::::::::::
 seq=j: --YYWCFWFCFWWF         (lf)MWYYWCFWFCFWWF        MWYYWCFWFCFWWF

 ********************************************************************************/
	MAT[0][0] = 0; 		// corresponds to S0 state.
	// DEL[0][0] = d2d[1]-m2d[1];  // --> DEL[1][0] = m2d. note m2d == s2d...
	// DEL[0][0] = -m2d[1];  // --> DEL[1][0] = -s2d - d2d. note m2d == s2d...
	DEL[0][0]=0;
	INS[0][0]=0;
	TB[0][0]=0;	
	for(i=1; i<= hmm_len; i++) {  // j==0 implies deletion of model position route.
	    DEL[i][0] = DEL[i-1][0] - d2d[i];	// penalty for deletions at start.
	    INS[i][0] = -infinity;	// disallow internal inserts from start of sequence.
					// requires d2i transitions.
	    // MAT[i][0] = 0;
	    MAT[i][0] = -infinity;  
	    TB[i][0] = TB[i-1][0] + 1; // deletion traceback route 
	}
	for(j=1; j<= sq_len; j++) { 
	    	MAT[0][j] = 0;  
if(global) DEL[0][j] = -infinity;
else 
		DEL[0][j] = 0;     // Make LOCAL with respect to sequence.
		// DEL[0][j] = - s2d[1];
		INS[0][j] = 0;	// no penalty for inserts before HMM...
		TB[0][j]= TB[0][j-1] -1; // insert traceback route 
		// TB[0][j]= -1; // traceback insert route 
	} DEL[0][0] = 0;  
	Int4 *gdel; NEW(gdel,sq_len+3,Int4);
	BooLean	*Internal; NEW(Internal,sq_len+3,BooLean);
	// 2. Dynamic programming step. 
	for(k=m=i=1; i<= hmm_len; i++) {
	   Int4	*mat_i=MAT[i],*mat_im1=MAT[i-1];
	   Int4 *ins_i=INS[i],*ins_im1=INS[i-1];
	   Int4 *del_i=DEL[i],*del_im1=DEL[i-1];
	   smx_typ  smx=M[m];
	   Int4	len_smx=LenSMatrix(smx);
	   register Int4 jj,jm1;
	   for(gINS=0,jm1=0,jj=1; jj<= sq_len; jj++,jm1++) {  
	   	Int4	t,s_max;
		//========= 1. Choose optimum between m->i & i->i transitions. ==========
#if 0	// this might be messing things up; Turned off by afn 11-7-2014.
	        if(k==len_smx){	// NO PENALTIES for inserts at ends!
		  s_mi = mat_i[jm1];	// insertion opening score
		  s_ii = ins_i[jm1];	// insertion extension score
	        } else 
#endif
		{
                  s_mi=mat_i[jm1]-m2i[i];	// insertion opening score
		  s_ii=ins_i[jm1]-i2i[i];	// insertion extension score
		}
		// save the maximum score and trace back...
                if(s_mi > s_ii){ t=gINS=-1; ins_i[jj] = s_mi; s_max=s_mi; }
		else { gINS--; t=gINS; ins_i[jj] = s_ii; s_max=s_ii; }
		// else { t=-1; ins_i[jj] = s_ii; s_max=s_ii; }
#if test_insert2delete		// Jun's HMM...
                if(k==len_smx){
		   s_di = del_i[jm1] - d2m[i];	// no penalty: d2s = 0;
		   if(s_di > s_mi && s_di > s_ii){
		     gINS=-1; ins_i[jj]=s_di; t=-1; s_max=s_di; 
		   }
		}
#endif

		//========== 2. Choose optimum between m->d & d->d transitions. ===========
#if 0	// bug fix?: afn 11-7-2014.	(This appeared to work somewhat...but not the best...)
		if(Internal[j]){	// 0 to i consists of deletions to the N-terminal end.
		   s_md=mat_im1[jj]-m2d[i];
		   s_dd=del_im1[jj]-d2d[i];
		} else {		// no deletion penalties on the N-terminal end.
		   s_md=mat_im1[jj];	// note m2d == s2d 
		   s_dd=del_im1[jj];	// deletion extension score
		}
#else	// original code
		s_md=mat_im1[jj]-m2d[i];	// note m2d == s2d 
		s_dd=del_im1[jj]-d2d[i];	// deletion extension score
#endif
                if(s_md > s_dd){
			gdel[jj]=1;  del_i[jj]=s_md; if(s_md > s_max){ s_max=s_md; t=1; }
		} else {
		 	gdel[jj]++; del_i[jj]=s_dd; if(s_dd > s_max){ s_max=s_dd; t=gdel[jj]; }
			// if(s_dd > s_max){ s_max=s_dd; t=1; }
		}
#if test_insert2delete		// Jun's HMM...
                if(k==1){
		   s_id = ins_im1[jj]-m2d[i];	// - s2d?
		   if(s_id > s_md && s_id > s_dd){
			gdel[jj]=1; del_i[jj]=s_id;
			if(s_id > s_max){ s_max=s_id; t=1; }
		   }
		}
#endif

		// 3. Choose optimum among m->m, i->m or d->m transitions.
		score = ValSMatrix(k,seq2[jj],smx);  // match score at i,j cell.
	        // if(k==1)
	        if(!Internal[j])
		{	// no s2m PENALTY for match at start of block!!!!
		  s_mm=mat_im1[jm1] + score;	// == m2m; no s->m penalty.
		  s_im=ins_im1[jm1] + score;
		  s_dm=del_im1[jm1] + score;
		  // s_dm=del_im1[jm1] + score - d2m[i];
	        } else {
		  s_mm=mat_im1[jm1] + score - m2m[i];
		  s_im=ins_im1[jm1] + score - i2m[i];
		  s_dm=del_im1[jm1] + score - d2m[i];
		} 
		if(s_mm >= s_dm){
		   if(s_mm >= s_im){ mat_i[jj] = s_mm; if(s_mm >= s_max){ s_max=s_mm; t=0; }
		   } else { mat_i[jj]=s_im; if(s_im >= s_max){ s_max=s_im; t=0; } }
		} else if(s_dm >= s_im){ mat_i[jj]=s_dm; if(s_dm >= s_max){ s_max=s_dm; t=0; }
		} else { mat_i[jj]=s_im; if(s_im >= s_max){ s_max=s_im; t=0; } }

		TB[i][jj]=t;	// record optimum route to i,j cell.
		if(!Internal[j] && t <= 0){ Internal[j] = TRUE; }
	   } k++;
	   if(k > LenSMatrix(M[m])) { k = 1; m++; }
	} free(gdel); free(Internal);
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
          TB[i][j] == 0            TB[i][j] == -2          TB[i][j] == 3
             s --> M                s --> I                 s --> D 
              
    model[i]   PT                     W--A                  PhwT
	       :.                     :  :                  :  :
    seq[j]     PS                     WpsA                  P--S

#endif	//********************************************************************
	Int4	s,t,t0,J0,max_i,max_j,end;
	max_i=hmm_len;		// global with respect to profile.
if(global){ max_j=sq_len; score=MAT[max_i][max_j];}
else {
	for(score=INT4_MIN, j=1; j<= sq_len; j++){
	   if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
	}
	// Determine whether optimum alignment ends in a match or a deletion.
	for(j=1; j<= sq_len; j++){
	   if((s=DEL[max_i][j]) > score){ score=s; max_j=j; }
	}
}
#if 0	// make alignment semi-local & symmetric by eliminating initial & ending m2d and d2m penalties.
	for(j=1; j<= sq_len; j++){
	   for(i=max_i; i > 0 && TB[i][j] > 0; i--) ;
	   if(TB[i][j] == 0){
	     if((s=MAT[i][j]) > score){ // then wipe out deletion penalties...
#if 1		// no m2d penalty.
		i++; DEL[i][j]= s-d2d[i];	// this should be m2d[i]
		for(i++; i <= max_i; i++) DEL[i][j]= DEL[i-1][j]-d2d[i];
	        if((s=DEL[max_i][j]) > score){ score=s; max_j=j; }
#else
	   	score=s; max_j=j; 
		for(i++; i <= max_i; i++) DEL[i][j]=s;
#endif
	     }
	   }
	}
#endif
#if debugger_on   //******************** debug using matrix output. *********************
#if 0
	FILE *fptr=stderr;
	Int4 startM=hmm_len-12,endM=hmm_len,endS=sq_len,startS=endS-15,blk=1;
	this->PutDPMatrix(fptr,MAT,TB,DEL,INS,sq_len,seq2,nmod,M,blk,startM,endM,startS,endS);
#endif
#endif  //*****************************************************************************
	// NOTE: insertions at end can't improve the score...
	NEW(operation,hmm_len+sq_len+3,char); operation[0]='E';
	m=nmod; k=LenSMatrix(M[m]);
	for(J0=0,j=sq_len; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
	// for(i=hmm_len; i > 0 || j > 0; ){  // full alignment...
	for(i=max_i; i > 0; ){  // full alignment...Stop at end of profile
	  t0 = TB[i][j];
	  // fprintf(stderr,"TB[%d][%d]=%d\n",i,j,t0);
	  do {
	    if(t0 > 0){ t=1; t0--; } else if(t0 < 0){ t=-1; t0++; } else { t=0; }
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
	} *start=j+1; J0++; operation[J0]='E';
#if 0	// DEBUG.
	if(sq_len==471){
		fprintf(stderr,"hmm: start = %d; max_j=%d; j=%d; operation=%s\n",*start,max_j,j,operation);
	}
#endif
	/*** 5. free allocated memory ***/
	for(i=0; i<=hmm_len; i++) {free(MAT[i]);free(TB[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(TB); free(DEL); free(INS); 
	*alnscore = score; *J = J0;
        return operation;
}

