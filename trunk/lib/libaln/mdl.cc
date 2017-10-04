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

#include "mdl.h"

mdl_typ::mdl_typ(Int4 N,Int4 *length,Int4 maxlen,Int4 *cnts,
	double npseudo,BooLean **null,a_type A)
{ 
    Int4	total,r;
    NumModels=N; AB=A; pseudo=npseudo;
    model = new fm_type [N+1];
    observedNS = new UInt4 [nAlpha(AB)+3];
    freq = new double [nAlpha(AB)+3];
    counts = new UInt4 [nAlpha(AB)+3];
    for(total=0,r=0; r<=nAlpha(AB); r++){
	counts[r]=observedNS[r]=cnts[r]; total+=cnts[r];
    }
    for(r=0; r<=nAlpha(AB); r++){ freq[r] =(double)counts[r]/(double)total; }
    for(Int4 t=1; t<= N; t++){
        model[t]=MkFModel(null[t],length[t],maxlen,npseudo,observedNS,counts,freq,A);
//  Eliminate: M->Ps[] From FModel()???
    } 
}

void    mdl_typ::copy(const mdl_typ& mdl)
{
    Free(); 
    NumModels=mdl.NumModels; AB=mdl.AB; pseudo=mdl.pseudo;
    model = new fm_type [NumModels+1];
    observedNS = new UInt4 [nAlpha(AB)+3];
    counts = new UInt4 [nAlpha(AB)+3];
    freq = new double [nAlpha(AB)+3];
    for(Int4 t=1; t<= NumModels; t++) model[t]=CopyFModel(mdl.model[t]);
    for(Int4 r=0; r<=nAlpha(AB); r++){
	freq[r]=mdl.freq[r]; observedNS[r]=mdl.observedNS[r];
	counts[r]= mdl.counts[r];
    }
}

smx_typ	*mdl_typ::SampleSMX(double pernats, Int4 wt)
{
    smx_typ	*smx;

    NEW(smx,NumModels+2,smx_typ);
    for(Int4 t=1;t<= NumModels;t++) smx[t]=SampleSmatrixFModel(pernats,wt,model[t]);
    return smx;
}

mdl_typ::mdl_typ(const mdl_typ& mdl){ copy(mdl); }
// constructor for 'mdl_typ mdl2=mdl;' or 'mdl_typ mdl2(mdl);

mdl_typ& mdl_typ::operator=(const mdl_typ& mdl){ copy(mdl); return *this; }

void	mdl_typ::Check( )
{
    Int4 r;
// std::cerr << "Checking for consistency of mdl_typ at end of run...\n";
    UInt4 *total = new UInt4 [nAlpha(AB)+1];
    for(r=0; r<=nAlpha(AB); r++) total[r]=0;
    for(Int4 m=1; m<=NumModels; m++){ 
	for(Int4 c=1; c<=LenFModel(model[m]); c++){
	    Int4 *col=SeeColumnFModel(c, model[m]);
	    if(col) for(r=0; r<=nAlpha(AB); r++) total[r]+=col[r];
	}
    }
    for(r=0; r<=nAlpha(AB); r++){
	if(observedNS[r] != (counts[r]-total[r])){
		fprintf(stderr,"observedNS[%d] = %d; counts[%d] = %d; total[%d]=%d\n",
			r,observedNS[r],r,counts[r],r,total[r]);
	} assert(observedNS[r]==(counts[r]-total[r]));
    } delete [] total;
}

void	mdl_typ::Free( )
{
  Check( ); // check at end for inconsistencies in observedNS[].
  if(model!=0){
     for(Int4 t=NumModels; t > 0; t--){
        if(model[t]!= 0) NilFModel(model[t]);
     } delete [] model;
  }
  delete [] observedNS; delete [] freq; delete [] counts;
}

// *********************** Sequence operations: ***********************
void    mdl_typ::Add(Int4 t, unsigned char *seq, Int4 site)
{
	register Int4		c;
	register fm_type	M=model[t];
	unsigned char		*sq=seq;
	for(c=LenFModel(M),sq+=(site+c); c > 0; c--){
		sq--; if(IsColumnFModel(c,M)) observedNS[*sq]--;
	} Add2FModel(seq, site, model[t]); 
}

void    mdl_typ::Remove(Int4 t, unsigned char *seq, Int4 site)
{
	register Int4		c;
	register fm_type	M=model[t];
	unsigned char	*sq=seq;
	for(c=LenFModel(M),sq+=(site+c); c > 0; c--){
		sq--; if(IsColumnFModel(c,M)) observedNS[*sq]++;
	} RmFModel(seq, site, model[t]); 
}

// *********************** Column operations: ***********************
Int4    mdl_typ::MvColumn(Int4 t, Int4 *obs, Int4 lemon, Int4 pos)
{
	Int4 *col=SeeColumnFModel(lemon, model[t]); assert(col);
	for(Int4 r=nAlpha(AB); r>=0; r--) observedNS[r]+=col[r]-obs[r];
	return MvColumnFModel(obs, lemon, pos, model[t]); 
}

Int4    mdl_typ::AddColumn(Int4 t,Int4 *obs,Int4 pos)
{
	for(Int4 r=nAlpha(AB); r>=0; r--) observedNS[r]-=obs[r];
	return AddColumnFModel(obs,pos,model[t]);
}

Int4	mdl_typ::RmColumn(Int4 t, Int4 pos)
{
	Int4 *col=SeeColumnFModel(pos, model[t]); assert(col);
	for(Int4 r=nAlpha(AB); r>=0; r--) observedNS[r]+=col[r];
	return RmColumnFModel2(pos, model[t]);
}

void    mdl_typ::Shift(Int4 t,Int4 *obs, BooLean left)
{
   Int4	*col;
   if(left) col=SeeColumnFModel(1, model[t]);
   else col=SeeColumnFModel(LenFModel(model[t]), model[t]);
   assert(col); 
   for(Int4 r=nAlpha(AB); r>=0; r--) observedNS[r]+=col[r]-obs[r];
   ShiftFModel(obs, left, model[t]); 
}

#if 0 // Potential new functions...
void    mdl_typ::calc_lngamma(cma_typ C)
// static variables are guaranteed to be initialized to zero.
// May want to use this here??? (Instead of cmsa array).
{
        Int4    i,b,N;
        a_type  A;
        double  Ps[30];

	N=NumSeqsCMSA(C); A = AlphabetCMSA(C);
        if(C->lngamma == NULL){
           NEWP(C->lngamma, N+2, double);
           for(i=0; i<=N; i++) NEW(C->lngamma[i], nAlpha(A)+2, double);
        }
        ss_type data = TrueDataCMSA(C); 
        for(b=0; b<=nAlpha(A); b++){ Ps[b]=pseudo*freq[b]; }
        for(i=0; i<=N; i++){
           for(b=0; b<=nAlpha(A); b++){
             if(freq[b] > 0.){
                C->lngamma[i][b]=lngamma((double)i + Ps[b]);
             } else C->lngamma[i][b] = 0.0;
           }
        }
}
#endif

// *********************** General operations: ***********************
void    mdl_typ::SetPseudo(Int4 NumSeq,double new_pseudo)
{
	pseudo=new_pseudo;
	if(NumSeq*pseudo > 20.0) pseudo=20.0/(double) NumSeq;
	for(Int4 t=1; t<=NumModels; t++){
		SetPseudoFModel(NumSeq*pseudo,model[t]);
	}
}

double  mdl_typ::LogProbRatio(Int4 NumSeq, Int4 *SeqLen, double indel_penalty)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (It always uses PSEUDO_CMSA*NumSeq.)
 **********************************************************************/
{
    Int4        n,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,*site,**observed;
    double      MAP,weight,total,Ps[30],v;
    double	npseudo,PSEUDO_MDL=0.1;

    for(len=0, t=1; t <= NumModels; t++) {
	if((n=MaxLenFModel(model[t])) > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    for(weight=0.0, s=1; s<=NumSeq; s++){ // weight for number of blocks.
	T = SeqLen[s] + NumModels;
	for(t=1; t <= NumModels; t++) T -= LenFModel(model[t]);
	weight += lnbico(T, NumModels);
    }
    npseudo = (double) NumSeq*PSEUDO_MDL;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(AB);b++){;Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,MAP=0.0,t=1; t <= NumModels; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	MAP -= lnbico(len - 2, ncol - 2); // column width weight.
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(AB); b++) {
		// MAP += L->lngamma[site[b]][b];  NEW...
	      if(Ps[b] > 0){
		MAP += lngamma((double)site[b] + Ps[b]); 
	      }
	    }
            MAP -= lngamma((double)(totsite-site[0]) + npseudo);
	    // MAP -= lngamma((double)totsite + npseudo);
	  } 
	}
    }
    v = lngamma(npseudo);
    for(total = 0.0, b=1;b<= nAlpha(AB); b++) { // NONSITES MAP.
	   if(Ps[b] > 0){
		MAP += lngamma((double) observedNS[b] + Ps[b]);
		total += (double) observedNS[b];
        	v -= lngamma(Ps[b]);   /** term pulling down MAP **/
	   }
    }
    MAP += v *(double) on_col; 
    MAP -= lngamma((double) total + npseudo);
    for(total=0.0,b=1;b<=nAlpha(AB);b++) { // subtract NULL MAP.
	   if(counts[b] > 0.0){
		MAP -= lngamma(counts[b]+Ps[b]); total+=counts[b];
	   }
    }
    MAP += lngamma((double)total + npseudo);
    if(NumModels > 1){		// WEIGHT FOR COLUMN TRANSFERS.
	weight +=  lnbico(on_col-NumModels-1,NumModels-1);
	if(off_col > 0) weight += lnbico(off_col+NumModels-1,NumModels-1);
    }
    free(observed);
    return (MAP - weight - indel_penalty);
}

///////////////// SW-ALIGNEMNT WITH GAP FUNCTION //////////////////
//
// Perform a Smith-Waterman alignment on mdl. M sequence seq2
// returns optimum subalignment score.  W(k) = a + bk.
// D[0][j]=D[i][0]=0;
//              MAX{ D[i-1][j-1] + S(r1,r2),
// D[i][j] =            D[i-k][j] - W(k) for k=1..i-1,
//                      D[i][j-k]- W(k) for k=1..j-1}.
//
// a = gap open; b = gap extend;
// Find maximum score and trace back to get alignment.
// see T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197.
//
//   Uses O(m*n) Gotoh method for penalties where W(k)=u*k+v where u,v >= 0.
//             .....
//             :DMI:                D = DEL; M = MATCH; I = INS;
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
char    *mdl_typ::gapped_aln_trace(Int4 a,Int4 b,Int4 n2,unsigned char *seq2, 
        Int4 **gapscore,Int4 *J,double *alnscore,Int4 *start)
{
        Int4    m,i,j,k,r1,s,t,v,n1,score,max_i,max_j,mscore;
        Int4    **MAT,**T,**DEL,**INS,**TB;
        Int4    j2,s1,s0,g,gs,g_opt,w=a+b,J0;
        unsigned char   *seq1;
        short   *mtf;
        a_type  A = AB;
        Int4    gINS,*gDEL,t0;
        char    *operation;
	double	pernats=1000.0;

        if(gapscore != NULL) gapscore--; // want gapscore[m-1] for block m.
        /** get total length of profile **/
        for(n1=0, m = 1; m <= NumModels; m++){ n1 += LenFModel(model[m]); }
        /** get concensus sequence for smatrix **/
        NEW(seq1, n1+3, unsigned char); 
        NEW(mtf, n1+3, short); NEWP(TB, NumModels+3, Int4); 
        for(s=0, m=1; m <= NumModels; m++){
            NEW(TB[m], n2+3, Int4); 
            for(i=1; i<= LenFModel(model[m]); i++){ s++; mtf[s] = m; }
        }
        /*** 1. Allocate and initialize memory. ***/
        MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
        MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
        for(i=0; i<= n1; i++) { 
                MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
                MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
        }
        // Make GLOBAL ALIGNMENT with respect to profile,
        MAT[0][0] = 0; DEL[0][0] = 0; INS[0][0] = 0;
        for(i=1; i<= n1; i++) {  
                DEL[i][0] = DEL[i-1][0] - b;
                INS[i][0] = INS[i-1][0] - b;
                MAT[i][0] = INS[i-1][0];  // Local on ends of profile 
                T[i][0] = 1; // for full alignment
        }
        for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
                MAT[0][j] = 0; T[0][j]=-1; // traceback full alignment 
                DEL[0][j] = INS[0][j] = SHRT_MIN;
        }
        // 2. Dynamic programming step. 
        NEW(gDEL,n2+3,Int4);
        for(k=m=i=1; i<= n1; i++) {
           gINS=0;
           if(k==LenFModel(model[m])) {
             for(j=1; j<= n2; j++) {  // Eliminate insert state...
                t=0; s=MAT[i-1][j-1] 
		  + (Int4)floor((pernats*ScoreFModel(seq2[j],k,model[m]) + 0.5));
                if(gapscore == NULL){ // no affine penalty for gapscore.
                   if((s0=MAT[i][j-1]) > (s1=INS[i][j-1])){
                        gINS=-1;  INS[i][j] = s0;
                   } else { gINS--; INS[i][j] = s1; }
                   if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }
                }
                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }
                T[i][j] = t; MAT[i][j] = s;
             }
           } else if(k==1 && m > 1 && gapscore != NULL) {
             for(j=1; j<= n2; j++) {  // Use gap function here.
                // Don't jump over gap function regions! Can't use t=gDEL[j].
                DEL[i][j] = MAXIMUM(Int4,MAT[i-1][j]-w,DEL[i-1][j]-b);
                s = DEL[i][j]; t=1;
                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }
		v=(Int4)floor((pernats*ScoreFModel(seq2[j],k,model[m]) + 0.5));
                for(s1=SHRT_MIN,j2=j-1, g=0; j2 > 0; g++,j2--){  
                    if((gs = gapscore[m][g]) == SHRT_MIN) break; // over max gap?
                    s0 = MAT[i-1][j2] + v + gs;
                    if(s1 < s0){ s1 = s0; g_opt = g; }
                }
                if(s < s1) { s = s1; t=0; TB[m][j] = -g_opt; }
                T[i][j]=t; MAT[i][j] = s;
             }
           } else {                   // Use affine gap penalty
              for(j=1; j<= n2; j++) {
		t=0; s=MAT[i-1][j-1] 
		   + (Int4)floor((pernats*ScoreFModel(seq2[j],k,model[m]) + 0.5));
                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }
                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }
                T[i][j] = t; MAT[i][j] = s;
              }
           }
           k++;
           if(k > LenFModel(model[m])) { k = 1; m++; }
        }
        // 2b. Find optimum global score with respect to profile.
        max_i=n1;               // global with respect to profile.
        for(score=INT4_MIN, j=1; j<= n2; j++)
                if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
        /*** 3. Trace back step. ***/
        // operations: 
        // 'i' = insertion in sequence outside of motifs
        // 'I' = insert in sequence within motif (need to delete this)
        //  'M' = match to start of a motif block
        //  'm' = match to other sites in motif
        //  'D' = deletion of sequence within motif
        //  'd' = deletion of sequence outside of motif (not used)
        NEW(operation,n1+n2+3,char); operation[0]='E';
        m=NumModels; k=LenFModel(model[m]);
        for(J0=0,j=n2; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
        // for(i=n1; i > 0 || j > 0; ){  // full alignment...
        for(i=n1; i > 0; ){  // full alignment...Stop at end of profile
          t0 = T[i][j];
          do {
            if(t0 > 0){ t=1; t0--; }
            else if(t0 < 0){ t=-1; t0++; } else { t=0; }
            switch(t){
                case 0:
                   i--; J0++;
                   if(k==1) operation[J0] = 'M'; else operation[J0] = 'm';
                   if(gapscore != NULL && k==1 && m > 1){
                        for(g=TB[m][j], j--; g < 0; g++){
                            J0++; operation[J0] = 'i'; j--;
                        }
                   } else j--; k--; break;
                case 1:  // Gap ('-') is in sequence; add X's for gaps
                        J0++; 
                        if(k==1) operation[J0] = 'D'; else operation[J0] = 'd';
                        i--;  k--; break;
                case -1:        // Gap ('-') is in profile
                        if(m < 1 || k==LenFModel(model[m]))
                                { J0++; operation[J0] = 'i'; }
                        else { J0++; operation[J0] = 'I'; }
                        j--; break;
                default: print_error("this should not happen"); break;
             }
             if(k==0){ m--; if(m > 0) k=LenFModel(model[m]); }
          } while(t0 != 0);
        }
        *start=j+1; J0++; operation[J0]='E';
        /*** 5. free allocated memory ***/
        for(m=1; m <= NumModels; m++) free(TB[m]); free(TB);
        for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
        free(MAT); free(T); free(DEL); free(INS); 
        free(seq1); free(mtf); free(gDEL);
        *alnscore = score; *J = J0;
        return operation;
}

char    *mdl_typ::GapAlnTrace(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
        Int4 **gapscore, Int4 *start)
// return the trace, the start position of the trace, and the trace length.
{
        Int4    i,j,newlen,lenTrace;
        char    *operation,*tmp;
	double	s;

        tmp=gapped_aln_trace(a,b,len,seq2,gapscore,&lenTrace,&s,start); 
        for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
        NEW(operation,newlen+3,char); operation[0]='E';
        for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
        operation[i]='E'; i++; operation[i]=0; free(tmp);
        return operation;
}

#if 0	// Alex's routine...
double  mdl_typ::DiffLogProbRatio(Int4 NumSeq, double oldindel, double newindel, 
        unsigned char **oldseg, unsigned char **newseg, mdl_typ *mdl)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the mdl->tmodel() pseudocounts.
 (It always uses PSEUDO_CMSA*NumSeq.)
 oldseg[t][i] & newseg[t][i] where t=1...mdl->tNumModels( ) & i = 1..lenModel[t]
 iff oldseg[t][i] != newseg[t][i] then put oldseg[t][i] into nonsites 
 & take newseg[t][i] from nonsite.
 **********************************************************************/
{
    double        DIFF;
    Int4          n,len,i,t,num_del;
    Int4          totsite,b,*site,**observed;
    double        total,Ps[30];
    double        npseudo,PSEUDO_MDL=0.1;
    Int4          *toNonsites;
    unsigned char new_b,old_b;
    fm_type       *model;
    
    model=mdl->tmodel();
    npseudo = (double) NumSeq*PSEUDO_MDL;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(mdl->tAB());b++){;Ps[b]=npseudo*mdl->tfreq( )[b];}

    for(len=0, t=1; t <= mdl->tNumModels( ); t++) {
        if((n=MaxLenFModel(mdl->tmodel()[t])) > len) len = n;
    }
    NEWP(observed,len+2,Int4);
    NEW(toNonsites,nAlpha(mdl->tAB())+1,Int4);
    for(b=0;b<=nAlpha(mdl->tAB());b++) toNonsites[b]=mdl->tobservedNS( )[b];

    for(total = 0.0,DIFF=0.,b=1;b<= nAlpha(mdl->tAB());b++) {
        DIFF -= lngamma((double) toNonsites[b] + Ps[b]);
        total += (double) toNonsites[b];
    } DIFF += lngamma((double) total + npseudo);

    for(t=1; t <= mdl->tNumModels( ); t++) {
        len = ObservedFModel(observed, model[t]);
        totsite = TotSitesFModel(model[t]);
        for(i=1; i <= len; i++) {
          site = observed[i];
          if(site != NULL){
              new_b=newseg[t][i]; old_b=oldseg[t][i];
              if(new_b != old_b){
                num_del=site[0];
                DIFF += lngamma((double)(totsite-num_del) + npseudo);
                if(old_b!=0) DIFF -= lngamma((double)site[old_b] + Ps[old_b]);
                if(new_b!=0) DIFF -= lngamma((double)site[new_b] + Ps[new_b]);
                if(old_b!=0){ toNonsites[old_b]++; } else num_del--;
                if(new_b!=0){ toNonsites[new_b]--; } else num_del++;
                if(new_b!=0) DIFF+=lngamma((double)site[new_b]+1 + Ps[new_b]);
                if(old_b!=0) DIFF+=lngamma((double)site[old_b]-1 + Ps[old_b]);
                DIFF -= lngamma((double)(totsite-num_del) + npseudo);
// std::cerr << DIFF; std::cerr << std::endl;
              } 
          } 
        }
    }
    for(total=0.0, b=1;b<= nAlpha(mdl->tAB()); b++) {
        // if ((double) toNonsites[b] + Ps[b]<=0) print_error("5");
        DIFF += lngamma((double) toNonsites[b] + Ps[b]);
        total += (double) toNonsites[b];
    } DIFF -= lngamma((double) total + npseudo);
    DIFF += oldindel; DIFF -= newindel;
    free(observed);free(toNonsites);
    return (DIFF);
}

#endif
#if 0
//***************************** NEW MAP ROUTINES *************************
double  mdl_typ::DiffLogProbRatio(Int4 NumSeq, Int4 *SeqLen, double indel_penalty,
	Int4 blk, Int4 oldcol)
// test to make sure that column 'col' in block 'blk' corresponds to
// a true column.
{ 
    Int4        n,len,i,t,T,s,ncol,off_col,on_col;
    Int4        totsite,b,*site,**observed;
    double      DIFF,total,Ps[30],v;
    double      npseudo,PSEUDO_MDL=0.1;

    for(len=0, t=1; t <= NumModels; t++) {
        if((n=MaxLenFModel(model[t])) > len) len = n;
    }

    NEWP(observed,len+1,Int4);

    npseudo = (double) NumSeq*PSEUDO_MDL;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(AB);b++){;Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,t=1; t <= NumModels; t++) {
        ncol = nColsFModel(model[t]); on_col += ncol;
        off_col += LenFModel(model[t]) - ncol;
        totsite = TotSitesFModel(model[t]);
    }    

    len = ObservedFModel(observed, model[blk]);

    if (oldcol<1 || oldcol > len) print_error("column not in block");

    DIFF += lnbico(len - 2, ncol - 2);
    if (oldcol != 1 || oldcol != len)
    	DIFF -= lnbico(len - 3, ncol - 3);
    else {if (oldcol = 1) {i=2; while (observed[i]==NULL) i++;}
    	  else {i=1; while (observed[len-i]==NULL) i++;} 

    	  DIFF -= lnbico(len - 2 - i, ncol - 3);

	  for(s=1; s<=NumSeq; s++){
             T = SeqLen[s] + NumModels;
             for(t=1; t <= NumModels; t++)
                 T -= LenFModel(model[t]);
             DIFF += lnbico(T, NumModels);
	     DIFF -= lnbico(T+i, NumModels);
	  }
    }

    if(NumModels > 1){
        DIFF +=  lnbico(on_col-NumModels-1,NumModels-1);
	DIFF -=  lnbico(on_col-NumModels-2,NumModels-1);
        if(off_col > 0) DIFF += lnbico(off_col+NumModels-1,NumModels-1);
	DIFF -= lnbico(off_col+NumModels,NumModels-1);
    }


    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
       DIFF -= lngamma((double) observedNS[b] + Ps[b]);
       total += (double) observedNS[b];
    }
    DIFF += lngamma((double) total + npseudo); 

    DIFF += lngamma((double)totsite + npseudo);

    site = observed[oldcol];
    for(b=1;b<= nAlpha(AB); b++) {
        DIFF -= lngamma((double)site[b] + Ps[b]); 
	totsite -= site[b]; observedNS[b] += site[b];
    }
    DIFF -= lngamma((double)totsite + npseudo);

    v = lngamma(npseudo);
    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
        DIFF += lngamma((double) observedNS[b] + Ps[b]);
        total += (double) observedNS[b];
        v -= lngamma(Ps[b]);
    }
    DIFF -= v;
    DIFF -= lngamma((double) total + npseudo);

    free(observed);
    return (DIFF);
}

double  mdl_typ::DiffLogProbRatio(Int4 NumSeq, Int4 *SeqLen, double indel_penalty,
	Int4 *newcol, Int4 blk, Int4 col)
// make sure that blk = 1..NumModels and col < 1 || col > LenFModel[blk] or
// (col > 0 && col <= LenFModel[blk] && col is turned off (i.e. site == NULL)
// col == 0 means the position just before start of block.
{ 
    Int4        n,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4        totsite,b,*site,**observed;
    double      DIFF,total,Ps[30],v;
    double      npseudo,PSEUDO_MDL=0.1;

    for(len=0, t=1; t <= NumModels; t++) {
        if((n=MaxLenFModel(model[t])) > len) len = n;
    }

    NEWP(observed,len+1,Int4);

    npseudo = (double) NumSeq*PSEUDO_MDL;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(AB);b++){;Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,t=1; t <= NumModels; t++) {
        ncol = nColsFModel(model[t]); on_col += ncol;
        off_col += ObservedFModel(observed, model[t]) - ncol;
        totsite = TotSitesFModel(model[t]);
    }
    len = ObservedFModel(observed, model[blk]);
    DIFF += lnbico(len - 2, ncol - 2);
    if (col  >= 1 || col <= len)
        DIFF -= lnbico(len - 2, ncol - 1);
    else {if (col < 1) {i=abs(col)+1;}
          else {i=col-len;} 

          DIFF -= lnbico(len - 2 + i, ncol - 1);

          for(s=1; s<=NumSeq; s++){
             T = SeqLen[s] + NumModels;
             for(t=1; t <= NumModels; t++)
                 T -= LenFModel(model[t]);
             DIFF += lnbico(T, NumModels);
             DIFF -= lnbico(T-i, NumModels);
          }
    }

    if(NumModels > 1){
        DIFF +=  lnbico(on_col-NumModels-1,NumModels-1);
        DIFF -=  lnbico(on_col-NumModels,NumModels-1);
        if(off_col > 0) DIFF += lnbico(off_col+NumModels-1,NumModels-1);
        DIFF -= lnbico(off_col+NumModels-2,NumModels-1);
    }


    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
       DIFF -= lngamma((double) observedNS[b] + Ps[b]);
       total += (double) observedNS[b];
    }
    DIFF += lngamma((double) total + npseudo); 

    DIFF += lngamma((double)totsite + npseudo);

    site = newcol;
    for(b=1;b<= nAlpha(AB); b++) {
        DIFF += lngamma((double)site[b] + Ps[b]); 
        totsite += site[b]; observedNS[b] -= site[b];
    }
    DIFF -= lngamma((double)totsite + npseudo);

    v = lngamma(npseudo);
    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
        DIFF += lngamma((double) observedNS[b] + Ps[b]);
        total += (double) observedNS[b];
        v -= lngamma(Ps[b]);
    }
    DIFF += v;
    DIFF -= lngamma((double) total + npseudo);

    free(observed);
    return (DIFF);
}

double  mdl_typ::DiffLogProbRatio(Int4 NumSeq, double oldindel, double newindel, 
	unsigned char **oldseg, unsigned char **newseg)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (It always uses PSEUDO_CMSA*NumSeq.)
 oldseg[t][i] &	newseg[t][i] where t=1...NumModels & i = 1..lenModel[t]
 iff oldseg[t][i] != newseg[t][i] then put oldseg[t][i] into nonsites 
 & take newseg[t][i] from nonsite.
 **********************************************************************/
{
    double	  DIFF;
    Int4          n,len,i,t;
    Int4          totsite,b,*site,**observed;
    double        total,Ps[30];
    double        npseudo,PSEUDO_MDL=0.1;
    unsigned char new_b,old_b;
    
    npseudo = (double) NumSeq*PSEUDO_MDL;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(AB);b++){;Ps[b]=npseudo*freq[b];}

    for(len=0, t=1; t <= NumModels; t++) {
        if((n=MaxLenFModel(model[t])) > len) len = n;
    }

    NEWP(observed,len+1,Int4);

    for(total = 0.0,DIFF=0.,b=1;b<= nAlpha(AB);b++) {
                DIFF -= lngamma((double) observedNS[b] + Ps[b]);
                total += (double) observedNS[b];
    }

    DIFF += lngamma((double) total + npseudo);

    for(t=1; t <= NumModels; t++) {
        len = ObservedFModel(observed, model[t]);
        totsite = TotSitesFModel(model[t]);

	DIFF+=lngamma((double)totsite + npseudo);

        for(i=1; i <= len; i++) {
          site = observed[i];
          if(site != NULL){
	      new_b=newseg[t][i]; old_b=oldseg[t][i];
	      if (new_b != old_b){
		if(new_b!=0){
                  DIFF += lngamma((double)site[new_b] + 1 + Ps[new_b]) 
			  - lngamma((double)site[new_b] + Ps[new_b]); 
		  observedNS[new_b]--;totsite++;
		}
		if(old_b!=0){
		  DIFF += lngamma((double)site[old_b] - 1 + Ps[old_b]) 
                          - lngamma((double)site[old_b] + Ps[old_b]);
		  observedNS[old_b]++;totsite--;
		}
              }
          } 
        }
           DIFF -= lngamma((double)totsite + npseudo);

    }

    for(total = 0.0, b=1;b<= nAlpha(AB); b++) {
                DIFF += lngamma((double) observedNS[b] + Ps[b]);
                total += (double) observedNS[b];
    }
    DIFF -= lngamma((double) total + npseudo);
    DIFF += newindel-oldindel;
    free(observed);
    return (DIFF);
}

#endif

