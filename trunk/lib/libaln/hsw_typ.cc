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

#include "hsw_typ.h"

void    FWriteHSW(FILE *fp,hsw_typ hsw)
{
      Int4	s,rtn;
      rtn=fwrite(hsw,sizeof(henikoff_wts_type),1,fp);
      if(rtn != 1) print_error("FWriteHSW() error");
      // rtn=fwrite(&hsw->Length,sizeof(Int4),1,fp);
      // rtn=fwrite(&hsw->NWtSq,sizeof(Int4),1,fp);
      rtn=fwrite(&nAlpha(hsw->AB),sizeof(Int4),1,fp);
      if(rtn != 1) print_error("FWriteHSW() error");
      for(s=1; s <= hsw->Length; s++){
          rtn=fwrite(hsw->Weight[s],sizeof(double),hsw->NWtSq+1,fp);
          if(rtn != hsw->NWtSq+1) print_error("FWriteHSW( ) error");

          rtn=fwrite(hsw->WtCnts[s],sizeof(double),nAlpha(hsw->AB)+1,fp);
          if(rtn != nAlpha(hsw->AB)+1) print_error("FWriteHSW( ) error");

          rtn=fwrite(hsw->WtFreq[s],sizeof(double),nAlpha(hsw->AB)+1,fp);
          if(rtn != nAlpha(hsw->AB)+1) print_error("FWriteHSW( ) error");
      }
      rtn=fwrite(hsw->fract_seq_aln,sizeof(double),hsw->Length+1,fp);
      if(rtn != hsw->Length+1) print_error("FWriteHSW( ) error");
}

hsw_typ	FReadHSW(FILE *fp,a_type AB, cma_typ cma)
{
      Int4	s,x,rtn;
      hsw_typ	hsw;
      MEW(hsw,1,henikoff_wts_type);
      rtn=fread(hsw,sizeof(henikoff_wts_type),1,fp);
      if(rtn != 1) print_error("FReadHSW( ) input error 0");
      hsw->AB=AB; hsw->mcma=cma;
      rtn=fread(&x,sizeof(Int4),1,fp);
      if(rtn != 1) print_error("FReadHSW( ) input error 0");
      if(x != nAlpha(hsw->AB)) print_error("FReadHSW( ) input error 1");
      if(hsw->NWtSq != NumSeqsCMSA(cma)){
         fprintf(stderr,"hsw->Length = %d; LengthCMSA(1,cma) = %d\n",hsw->Length,LengthCMSA(1,cma));
         fprintf(stderr,"hsw->NWtSq = %d; NumSeqsCMSA(cma) = %d\n",hsw->NWtSq,NumSeqsCMSA(cma));
	 print_error("FReadHSW( ) input error 2");
      }
      if(hsw->Length != LengthCMSA(1,cma)) print_error("FReadHSW( ) input error 3");
      NEWP(hsw->Weight,hsw->Length + 3,double);
      NEWP(hsw->WtCnts,hsw->Length + 3,double);
      NEWP(hsw->WtFreq,hsw->Length + 3,double);
      for(s=1; s <= hsw->Length; s++){
          NEW(hsw->Weight[s],hsw->NWtSq + 3,double);
          rtn=fread(hsw->Weight[s],sizeof(double),hsw->NWtSq+1,fp);
          if(rtn != hsw->NWtSq+1) print_error("FReadHSW( ) error 4");

          NEW(hsw->WtCnts[s],x + 3,double);
          rtn=fread(hsw->WtCnts[s],sizeof(double),nAlpha(hsw->AB)+1,fp);
          if(rtn != nAlpha(hsw->AB)+1) print_error("FReadHSW( ) error 5");

          NEW(hsw->WtFreq[s],x+ 3,double);
          rtn=fread(hsw->WtFreq[s],sizeof(double),nAlpha(hsw->AB)+1,fp);
          if(rtn != nAlpha(hsw->AB)+1) print_error("FReadHSW( ) error 6");
      }
      NEW(hsw->fract_seq_aln,hsw->Length+ 3,double);
      rtn=fread(hsw->fract_seq_aln,sizeof(double),hsw->Length+1,fp);
      if(rtn != hsw->Length+1) print_error("FReadHSW( ) error 7");
      return hsw;
}

hsw_typ AddRandomHSW(hsw_typ HSW, cma_typ mcma,cma_typ rcma, cma_typ ccma)
// Add the (presumed) random set of sequences in rcma to end of HSW.
// Calling environment has (presumably) merged the two sets to get the complete set: ccma.
{
        hsw_typ hsw=0;
	Int4	r,s,j,Length,NWtSq;
	a_type	AB;

	assert(HSW->Length == TotalLenCMSA(rcma));
	assert(HSW->Length == TotalLenCMSA(mcma));
	assert(HSW->Length == TotalLenCMSA(ccma));
	assert(nBlksCMSA(rcma)==1);
	assert(nBlksCMSA(mcma)==1);
	assert(nBlksCMSA(ccma)==1);

	Length=HSW->Length; AB=HSW->AB;

	NEW(hsw,1,henikoff_wts_type);
        hsw->Length=HSW->Length;
        hsw->AB=HSW->AB; hsw->mcma=ccma;
        NWtSq=NumSeqsCMSA(mcma) + NumSeqsCMSA(rcma);
	assert(NWtSq == NumSeqsCMSA(ccma));
        hsw->NWtSq=NWtSq;
	NEWP(hsw->Weight,Length+2,double); NEWP(hsw->WtCnts,Length+2,double);
	NEWP(hsw->WtFreq,Length+2,double); NEW(hsw->fract_seq_aln,Length+2,double);
	for(s=1; s <= Length; s++){
		NEW(hsw->WtCnts[s],nAlpha(AB)+2,double);
		NEW(hsw->WtFreq[s],nAlpha(AB)+2,double);
                NEW(hsw->Weight[s],hsw->NWtSq+2,double);
		Int4 sq,rsq=0;
		for(sq=1; sq <= NWtSq; sq++){
		   e_type cE=TrueSeqCMSA(sq,ccma);
		   if(sq <= HSW->NWtSq){
			e_type mE=TrueSeqCMSA(sq,mcma);
			if(!FastIdentSeqs(mE,cE)){
			   print_error("FATAL: non-identical sequences within AddRandomHSW().\n");
			}
			double	wt=HSW->Weight[s][sq];
			hsw->Weight[s][sq]=wt;
			r = ResidueCMSA(1,sq,s,mcma);
			if(r > 0) hsw->WtCnts[s][r] += wt;
		   } else {
			rsq++;
			e_type rE=TrueSeqCMSA(rsq,rcma);
			if(strcmp(PhylumSeq(rE),"simulated") != 0 || 
							strcmp(PhylumSeq(cE),"simulated") != 0){
			   print_error("FATAL: non-random sequence within AddRandomHSW().\n");
			}
			if(!FastIdentSeqs(rE,cE)){
			   print_error("FATAL: non-identical sequences within AddRandomHSW().\n");
			}
			hsw->Weight[s][sq]=1.0;	// random sequences get a weight == 1.0.
			r = ResidueCMSA(1,rsq,s,rcma);
			if(r > 0) hsw->WtCnts[s][r] += 1.0;
		   }
		} assert(rsq == NumSeqsCMSA(rcma));
		double total_wt_cnts=0.0;
          	for(r=0; r<=nAlpha(AB); r++){ total_wt_cnts+=hsw->WtCnts[s][r]+0.01; }
          	for(r=0; r<=nAlpha(AB); r++){
                	hsw->WtFreq[s][r] = (hsw->WtCnts[s][r]+0.01)/total_wt_cnts;
		}
		hsw->fract_seq_aln[s] = HSW->fract_seq_aln[s]; // doesn't matter for mcBPPS...
	} return hsw;
}

hsw_typ GetSubHSW(hsw_typ HSW,set_typ Set, cma_typ cma,cma_typ mcma)
// pass in main cma file and Set == subset of those sequences for rtn_hsw.
{
        hsw_typ hsw=0;
	Int4	r,R,s,sq,i,j,Length,NWtSq;
	a_type	AB;

	assert(HSW->Length == TotalLenCMSA(cma));
	assert(HSW->Length == TotalLenCMSA(mcma));
	assert(nBlksCMSA(cma)==1); assert(nBlksCMSA(mcma)==1);

	Length=HSW->Length; AB=HSW->AB;

	NEW(hsw,1,henikoff_wts_type);
        hsw->Length=HSW->Length;
        hsw->AB=HSW->AB; hsw->mcma=cma;
        NWtSq=NumSeqsCMSA(mcma);
        hsw->NWtSq=CardSet(Set);
	assert(NWtSq >= hsw->NWtSq && CardSet(Set) == NumSeqsCMSA(cma));
#if 1	// DEBUG: make sure that sequence sets match!
	e_type	mE,sE;
	for(sq=1,i=0; sq <= NWtSq; sq++){
	   if(MemberSet(sq,Set)){ 
	        mE = TrueSeqCMSA(sq,mcma);
		i++; sE=TrueSeqCMSA(i,cma); assert(IdentSeqs(mE,sE)); 
	   }
	}
#endif

	NEWP(hsw->Weight,Length+2,double); NEWP(hsw->WtCnts,Length+2,double);
	NEWP(hsw->WtFreq,Length+2,double); NEW(hsw->fract_seq_aln,Length+2,double);
	for(s=1; s <= Length; s++){
		NEW(hsw->WtCnts[s],nAlpha(AB)+2,double);
		NEW(hsw->WtFreq[s],nAlpha(AB)+2,double);
                NEW(hsw->Weight[s],hsw->NWtSq+2,double);
		for(sq=1,i=0; sq <= NWtSq; sq++){
		   if(MemberSet(sq,Set)){
			double	wt=HSW->Weight[s][sq];
			i++; hsw->Weight[s][i]=wt;
			r = ResidueCMSA(1,sq,s,mcma);
			// R = ResidueCMSA(1,i,s,cma); assert(r==R);
			if(r > 0) hsw->WtCnts[s][r] += wt;
		   }
		} assert(i == hsw->NWtSq);
		double total_wt_cnts=0.0;
          	for(r=0; r<=nAlpha(AB); r++){ total_wt_cnts+=hsw->WtCnts[s][r]+0.01; }
          	for(r=0; r<=nAlpha(AB); r++){
                	hsw->WtFreq[s][r] = (hsw->WtCnts[s][r]+0.01)/total_wt_cnts;
		}
		hsw->fract_seq_aln[s] = HSW->fract_seq_aln[s]; // doesn't matter for mcBPPS...
	} return hsw;
}
	
void	NilHSW(hsw_typ HSW)
{
	Int4	s;
	for(s=1; s <= HSW->Length; s++){
		free(HSW->WtCnts[s]); free(HSW->WtFreq[s]); free(HSW->Weight[s]);
	} free(HSW->WtCnts); free(HSW->WtFreq); free(HSW->Weight); free(HSW->fract_seq_aln);
	free(HSW);
}



