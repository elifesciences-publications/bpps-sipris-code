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

#include "lpr_typ.h"

#define LPR_USAGE_START "FATAL: lpr_typ cma_file options\n\
   USAGE: GroupName [options}\n\
     -A<real>:<real> - alpha hyperparameters A0:B0 (default: A0=B0=1.0)\n\
     -global       - Use global sequence weighting\n\
     -P=<str>      - seed pattern string\n\
     -Ri=<real>    - Set prior probability that a row (seq) is in the foreground (default: 0.5)\n\
     -rho=<real>   - set prior probability (rho) that a column is a pattern position (default: 0.5)\n\
                        input is log(rho); e.g., rho=100 --> rho=e^-100 = column 'penalty' of -100 nats\n\
   \n\n"

lpr_typ::lpr_typ(cma_typ in_cma,swt_typ *in_swt)
{ del_as_random=0; cma=in_cma; swt=in_swt; AB=AlphabetCMSA(cma); Init(); }

lpr_typ::lpr_typ(cma_typ in_cma,swt_typ *in_swt,BooLean DelAsRandom)
{ if(DelAsRandom) del_as_random=1; else del_as_random=0; cma=in_cma; swt=in_swt; AB=AlphabetCMSA(cma); Init(); }

void	lpr_typ::usage( ){ print_error(LPR_USAGE_START); }

void	lpr_typ::Free()
{
        for(Int4 s=1; s <= Length; s++){ free(CntFG[s]); free(CntBG[s]); free(SqWt[s]); }
        free(CntFG); free(CntBG); free(SqWt); free(AveSqIWt);
	delete RST;
}

void	lpr_typ::Init()
{
	Length=LengthCMSA(1,cma); assert(nBlksCMSA(cma) ==1);
	UseGlobalSqWts=TRUE;
	RST=new rst_typ('L'); LegalSST=RST->LegalResSets();

        Int4    N=NumSeqsCMSA(cma);
        NEWP(CntFG,Length +3,UInt4); NEWP(CntBG,Length +3,UInt4);
        for(Int4 s=1; s <= Length; s++){
              NEW(CntFG[s],N+3,UInt4); NEW(CntBG[s],N+3,UInt4);
        } AveSqIWt=swt->GetIntegerWts(&SqWt);   // creates SqWt[Length][N] matrix.
        assert(swt->TheLength() == Length); assert(swt->NumWtSeqs() == N);
}

//==================== BPPS routines =======================

sst_typ *lpr_typ::GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,BooLean B,double &L,Int4 x, double pRi)
{ unsigned char *csq; sst_typ *xsst=GetOptPttrnLPR(f,S1,S2,B,L,x,csq,pRi); free(csq); return xsst;  }

sst_typ *lpr_typ::GetOptPttrnLPR(FILE *fp,set_typ SetFG, set_typ SetBG,BooLean Negate,double &lpr,
				Int4 MaxCols,unsigned char *&rtn_csq,double pRi)
{ return GetOptPttrnLPR(fp,SetFG,SetBG,Negate,lpr,MaxCols,rtn_csq,'G',0,pRi); }

sst_typ *lpr_typ::GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,BooLean B,double &L,Int4 x,char Type,
				e_type qE,double pRi)
{ unsigned char *csq; sst_typ *xsst=GetOptPttrnLPR(f,S1,S2,B,L,x,csq,Type,qE,pRi); free(csq); return xsst;  }

sst_typ	*lpr_typ::GetOptPttrnLPR(FILE *fp,set_typ SetFG, set_typ SetBG,BooLean Negate, double &lpr,
		Int4 MaxCols, unsigned char *&rtn_csq, char Type,e_type qE,double pRi)
{
	Int4		s,r,r_max;
	UInt4		cnt,cnt_max;
	double		d_max,d,*p;
	sst_typ		*qsst,best_sst;
	unsigned char	*seq;

	InitCntsFG_BG(SetFG,SetBG,Negate,UseGlobalSqWts); // Get weighted counts for FG & BG seqs to construct qsst.

	// 1. Get the consensus pattern for FG...
	NEW(seq,Length + 3, unsigned char); NEW(qsst,Length + 3, sst_typ);
	for(s=1; s <= Length; s++){
	   for(cnt_max=0,r_max=0,r=1; r <= nAlpha(AB); r++){
		if((cnt=CntFG[s][r]) > cnt_max){ cnt_max=cnt; r_max = r; }
	   } if(r_max > 0) {
		seq[s]=r_max; qsst[s]=SsetLet(r_max); assert(MemSset(seq[s],qsst[s]));
	   } // else assert(r_max > 0); // only occurs if FG has lots of gap residues '-' at site s.
        } // PrintSeq(stderr,seq);
	// 2. Get the initial pps structure;
        bpps_typ *pps=0;
	if(qE){ assert(Length == LenSeq(qE)); pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi,SeqPtr(qE)); }
	else pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi);	
pps->TreatDeletionAsRandom();
	sst_typ *opt_sst=pps->RtnOptPattern(CntBG, CntFG, LegalSST,seq,MaxCols,lpr);
	if(fp) pps->PutSubLPR(fp,0,CntBG,CntFG);
	// 3. free pps memory and return pattern and lpr;
	FreeBPPS(pps); rtn_csq=seq; free(qsst);
	return opt_sst;
}

double	lpr_typ::CalcSetvsPttrnLPR(FILE *fp,set_typ FG, set_typ BG,sst_typ *qsst, BooLean Negate,
			e_type qE,double pRi)
{ return CalcSetvsPttrnLPR(fp,FG, BG,qsst, Negate,'S',qE,pRi); }

double	lpr_typ::CalcSetvsPttrnLPR(FILE *fp,set_typ FG, set_typ BG,sst_typ *qsst, BooLean Negate,char Type,
				e_type qE,double pRi)
{
        bpps_typ *pps=0;
	if(qE){ assert(Length == LenSeq(qE)); pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi,SeqPtr(qE)); }
	else pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi);	
	InitCntsFG_BG(FG,BG,Negate,UseGlobalSqWts);	// CntBG & CntFG initialized );
	double	  lpr=pps->LPR(CntBG,CntFG);
	if(fp && lpr > 0.0) pps->PutSubLPR(fp,CntBG,CntFG);
	FreeBPPS(pps);
	return lpr;
}

void	lpr_typ::PutParameters(FILE *fp,set_typ SetFG, set_typ SetBG,sst_typ *qsst, BooLean Negate,
		char Type, double pRi)
{
        bpps_typ *pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi);
        InitCntsFG_BG(SetFG,SetBG,Negate,UseGlobalSqWts);       // CntBG & CntFG initialized );
        double    lpr=pps->LPR(CntBG,CntFG);
        // if(fp && lpr > 0.0) pps->PutSubLPR(fp,CntBG,CntFG);
	pps->PutParameters(stderr);
	pps->PutRho(stderr);
	FreeBPPS(pps);
}

void	lpr_typ::CompareBPPS(FILE *fp,bpps_typ *PPS,set_typ SetFG, set_typ SetBG,sst_typ *qsst,
								BooLean Negate,char Type,double pRi)
{
        bpps_typ *pps=MakeBPPS(qsst,UseGlobalSqWts,Type,pRi,PPS->RtnQuery());
        InitCntsFG_BG(SetFG,SetBG,Negate,UseGlobalSqWts);       // CntBG & CntFG initialized );
        double    d=pps->LPR(CntBG,CntFG),D=PPS->LPR(CntBG,CntFG);
	if(D != d) fprintf(stderr,"Lpr differ (%.2f vs %.2f)\n",D,d);
	else fprintf(stderr,"Lpr the same (%.2f = %.2f)\n",D,d);
        // if(fp && lpr > 0.0) pps->PutSubLPR(fp,CntBG,CntFG);
	pps->Compare(stderr,PPS);
	FreeBPPS(pps);
}

void	lpr_typ::CompareSqWts(FILE *fp,set_typ SetFG,set_typ SetBG,UInt4 **WtSqFG,UInt4 **WtSqBG,BooLean Neg)
{
	Int4	n=0;
        InitCntsFG_BG(SetFG,SetBG,Neg,UseGlobalSqWts);       // CntBG & CntFG initialized );
	fprintf(fp,"Comparing SqWts...");
	for(Int4 i=1; i <= Length; i++){
	   for(Int4 r=1; r <= nAlpha(AB); r++){
		if(CntFG[i][r] != WtSqFG[i][r]){
		   fprintf(fp,"\n%d: SqWts differ ",i); n++;
	   	   for(Int4 x=1; x <= nAlpha(AB); x++){
			fprintf(fp,"%c%d:%d ",AlphaChar(x,AB),CntFG[i][x],WtSqFG[i][x]);
		   } break;
		}
	   }
	} if(n==0) fprintf(fp,"they're the same.\n"); else fprintf(fp,"\n  %d different!\n",n); 
}

bpps_typ *lpr_typ::MakeBPPS(sst_typ *qsst,BooLean UseGlobalSqWts, char Type, double pRi, unsigned char *Qry)
{ bpps_typ *pps = new bpps_typ(qsst,Qry,Length,Type,AB,UseGlobalSqWts,pRi); return pps; }

void	lpr_typ::FreeBPPS(bpps_typ *pps) { delete pps; }

void	lpr_typ::InitCntsFG_BG(set_typ SetFG, set_typ SetBG,BooLean Negate,BooLean UseGlobalSqWts)
{
#if 1	// Use global weights only.
    register Int4 s,r,sq;
    register unsigned char *seq;
    for(s=1; s <= Length; s++){ for(r=0; r <= nAlpha(AB); r++){ CntFG[s][r]=0; CntBG[s][r]=0; } }
    for(sq=1; sq<=NumSeqsCMSA(cma); sq++){
	seq=GetAlnResInSiteCMSA(1,sq,cma);
	if(MemberSet(sq,SetFG)){		// == FG set...
	   for(s=1; s <= Length; s++) CntFG[s][seq[s]]+=AveSqIWt[sq];
	} else if((!Negate && MemberSet(sq,SetBG)) || (Negate && !MemberSet(sq,SetBG))){
	   for(s=1; s <= Length; s++) CntBG[s][seq[s]]+=AveSqIWt[sq];
	}
    }
#else
    register Int4 s,r,sq;
    register unsigned char *seq;
    for(s=1; s <= Length; s++){ for(r=0; r <= nAlpha(AB); r++){ CntFG[s][r]=0; CntBG[s][r]=0; } }
    if(UseGlobalSqWts){
      for(sq=1; sq<=NumSeqsCMSA(cma); sq++){
	seq=GetAlnResInSiteCMSA(1,sq,cma);
	if(MemberSet(sq,SetFG)){		// == FG set...
	   for(s=1; s <= Length; s++) CntFG[s][seq[s]]+=AveSqIWt[sq];
	} else if((!Negate && MemberSet(sq,SetBG)) || (Negate && !MemberSet(sq,SetBG))){
	   for(s=1; s <= Length; s++) CntBG[s][seq[s]]+=AveSqIWt[sq];
	}
      }
    } else {
      for(sq=1; sq<=NumSeqsCMSA(cma); sq++){
	seq=GetAlnResInSiteCMSA(1,sq,cma);
	if(MemberSet(sq,SetFG)){		// == FG set...
	   for(s=1; s <= Length; s++) CntFG[s][seq[s]]+=(UInt4) SqWt[s][sq];
	} else if((!Negate && MemberSet(sq,SetBG)) || (Negate && !MemberSet(sq,SetBG))){
	   for(s=1; s <= Length; s++) CntBG[s][seq[s]]+=(UInt4) SqWt[s][sq];
	}
      }
    }
#endif
}

double	lpr_typ::WtCardFG_BG_Sets(double &WtCntsFG, double &WtCntsBG)
// compute the number of weighted foreground & background counts.
// WARNING: assumes that CalcSetvsPttrnLPR( ) was called to initialize CntBG & CntFG.
{
	double	CntsFG=0.0,CntsBG=0.0;
	for(Int4 s=1; s <= Length; s++){
	   for(Int4 r=0; r <= nAlpha(AB); r++){
		CntsFG += (double) CntFG[s][r]; CntsBG += (double) CntBG[s][r];
	   }
	} 
	CntsFG = CntsFG / (double)(100*Length); CntsBG = CntsBG / (double)(100*Length);
	WtCntsFG=CntsFG; WtCntsBG=CntsBG;
	return (CntsFG+CntsBG);
}

