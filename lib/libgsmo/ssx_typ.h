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

#if !defined (_SSX_TYP_)
#define _SSX_TYP_
#include "smatrix.h"
#include "cmsa.h"
#include "dms_typ.h"
#include "jlh_typ.h"
#include "blosum62.h"
#include "p2p_typ.h"
#include "ndl_typ.h"
#include "lgm_typ.h"

class ssx_typ {		// sequence sampling extensed type.
public:
	ssx_typ(Int4 aa_per_io,Int4 aa_per_do,Int4 exp_ie, Int4 exp_de, double pn,cma_typ c,
				char dm,ssx_typ *twin_ssx=0);
	~ssx_typ(){ Free(); }
	double	AdjstdBildLLR(double targetWt, BooLean add_gp=FALSE){
		// adjust LLR to compare sequences with different WtNumSeqs.
		   if(ReCalcSMX) CalcSMX();
		   double	sqwt=this->RtnWtNumSeqs();
		   double	d,dd=0.0,Ratio=targetWt/sqwt;
		   UInt4	x,xWtCnts[50];
		   for(Int4 b=1; b <= nBlksCMSA(CMA); b++) {
        		Int4 i,r,len = LengthCMSA(b,CMA);
        		for(i=1; i <= len; i++) {
			   for(r=0; r <= nAlpha(AB); r++){
				x=WtCnts[b][i][r];
				if(x!=0){
				  d=(double)x*Ratio; x=(UInt4) ceil(d-0.49);
				} xWtCnts[r]=x;
				//fprintf(stderr,"%c%d: %d\n",AlphaChar(r,AB),i,WtCnts[b][i][r]);
			   }
			   // dd += DMS->bild(WtCnts[b][i]); 
			   dd += DMS->bild(xWtCnts); 
			}
		   }
		   double gp=0.0; if(add_gp) gp=NDL->IndelPenalty(0,'M',Wt);
		   return dd+gp;
		}
	double	GapMap(double Temp=0.0,char method=' ')
		{ 
		   char mode='S'; if(Temp <= 10.0) mode='M';
		   if(ReCalcSMX) CalcSMX();
		   if(method =='B') return this->BildLLR('M');
		   else if(method =='b') return this->RawBildLLR();
		   else {
#if 0
		     // if(ModeSMX=='W') return VirtualMap(mode,Blk,Col);
		     if(ModeSMX=='W' && Blk > 0) return VirtualMap(mode,Blk,Col);
		     else if(ModeSMX=='W') return DirichletLLR(mode);
#else
		     if(ModeSMX=='W') return DirichletLLR(mode);
#endif
		     else return this->RelMap(mode);
		   } 
		}
	ndl_typ	*RtnNDL(){ return NDL; }
	Int4	RtnFloorI2I(){ return NDL->RtnFloorI2I(); }
	Int4	RtnFloorD2D(){ return NDL->RtnFloorD2D(); }
	void    PutPenalties(FILE *fp, Int4 blk, Int4 I){ NDL->PutPenalties(fp,blk,I); }
	void    PutPenalties(FILE *fp){ this->GapMap(); NDL->PutPenalties(fp); }
	double	RtnIndelPenalty(){ if(ReCalcSMX) CalcSMX(); return NDL->IndelPenalty(0,'S',Wt); }
	double  VirtualMap2(char mode='S');
	double  VirtualMap(char mode='S',Int4 Blk=0,Int4 Col=0);
	double  BildLLR(char mode='S'){ return CoreBildLLR(TRUE,mode); }
	double  RawBildLLR(char mode='S'){ return CoreBildLLR(FALSE,mode); }
	double  CoreBildLLR(BooLean add_indels, char mode='S');
	double	BildScore(Int4 blk, Int4 col){ return DMS->bild(WtCnts[blk][col]); }
	double	BildScore(UInt4 *xWtCnts){ return DMS->bild(xWtCnts); }
	double  ContextBildScore(Int4 blk, Int4 i, double *BldSc=0,double blk_cut=0.20);
	// double	BildScore(Int4 blk, Int4 col){ return DMS->bild(VrtlCnts[blk][col]); }
	double	RelEntropy(Int4 blk, Int4 col);
	double	Entropy(Int4 blk, Int4 col);
	void	SetPriorWt(double D){ NDL->SetPriorWt(D); }
	double	GetPriorWt( ){ return NDL->GetPriorWt( ); }
	double	GetParameters(Int4 &io,Int4 &d_o,Int4 &ie, Int4 &de, double &pw, char &dm)
		  { return NDL->GetParameters(io,d_o,ie,de,pw); dm=dms_mode; }
	void	SetSqWtAdjust(double D)
			{ assert(D < 1.0 && D > 0.0); SqWtAdjst=D;  this->UpdateSeqWts(); }
	double  DirichletLLR(char mode='S',FILE *efp=0);
private:
	double	SqWtAdjst;
	ssx_typ	*twin;
public:
	double  AdjstDirichletLLR(double targetWt,FILE *efp=0);
	//============ New =====================
	double  *EstVirtualColLLR(Int4 Blk);

        //======================== ssx_typ.cc ============================
	Int4	**RemoveFromAlign(Int4 *cluster);
	Int4	*RemoveFromAlign(Int4 sq);
	void	RmFromAlign(Int4 *cluster);
	void	RmFromAlign(Int4 sq);
        void    AddToAlign(Int4 sq,Int4 *site);
        void    ReplaceSeq(Int4 sq, gsq_typ *gsq[]);
        void    UpdateSeqWts(){ SeqWts(); InitCnts(); NDL->InitIndels(Wt); }
	void    SampleAlignSeq(FILE *fp,char mode, e_type qE,double temp=300.0);
	void    AlignSeq(FILE *fp,char mode, e_type qE, double temp=300.0);
	char	*AlignSeq(e_type qE, Int4 &start,double temp=0.0,BooLean global=FALSE);
	BooLean ChangePriors(Int4 pio,Int4 pdo,Int4 eie,Int4 ede)
			{ NDL->ChangePriors(pio,pdo,eie,ede); }
#if 0	// core dumping using the following...
	void	RestoreBest(){ InitMAPCMSA(CMA); NDL->InitIndels(Wt); }	// Don't change weights!
#else	// also update sequence weights, as sequences are realigned!!!
	// void	RestoreBest(){ InitMAPCMSA(CMA); UpdateSeqWts(); }
	void	RestoreBest(){ InitMAPCMSA(CMA); InitCnts(); NDL->InitIndels(Wt); }
	void	ReInitCnts(){ InitCnts(); NDL->InitIndels(Wt); }
#endif
	// void	RestoreBest(){ InitMAPCMSA(CMA); SeqWts(); InitCnts(); }
	// void	UpdateWts(){ SeqWts(); InitCnts(); }
        //======================== ssx_scores.cc ============================
	double  **FractionDeleted( );
	double  FractionDeleted(Int4 blk, Int4 site);
	char	*GapAlnTrace(e_type sbjE, double Temp, Int4 &start,Int4 &score);
	cma_typ	RtnCMA( ){ return CMA; }
	double	TotalWtSeq( ){ return (double)WtN/(double)wt_factor; }
	smx_typ	*RtnSMX( ){ return SMX; }
	smx_typ	*RtnStraightSMX( ){ return StraightSMX(); }
	smx_typ	*RtnWtCntsSMX( ){ return WtCntsSMX(); }
	double	RtnPercentGaps(Int4 blk, Int4 col){
			UInt4 r,T=0,*xWtCnts=WtCnts[blk][col];
			for(r=0; r <= nAlpha(AB); r++){ T+= xWtCnts[r]; }
			return (100*(double)xWtCnts[0]/(double)T);
		}
	void	PutWtCnts(FILE *fp,Int4 blk, Int4 col){ this->PutWtCnts(fp,WtCnts[blk][col]); }
	void	PutWtCnts(FILE *fp,UInt4 *xWtCnts){
			UInt4 r,T=0;
#if 0
			for(r=0; r <= nAlpha(AB); r++){ T+= WtCnts[blk][col][r]; }
			for(r=0; r <= nAlpha(AB); r++){
				fprintf(fp,"%c%.0f ",AlphaChar(r,AB),
					100*(double)WtCnts[blk][col][r]/(double)T);
#elif 0
			for(r=0; r <= nAlpha(AB); r++){
				fprintf(fp,"%c%u ",AlphaChar(r,AB), xWtCnts[r]);
#else
			for(r=0; r <= nAlpha(AB); r++){ T+= xWtCnts[r]; }
			for(r=0; r <= nAlpha(AB); r++){
				fprintf(fp,"%c%.1f ",AlphaChar(r,AB),
					100*(double)xWtCnts[r]/(double)T);
#endif
			} fprintf(fp,"\n");
		}
	void	PutVrtlWtCnts(){
			for(Int4 r=1; r <= nAlpha(AB); r++)
				fprintf(stderr,"%c: %u\n",AlphaChar(r,AB),TtlVrtlCnts[r]);
		}
        //======================== for gmb_p2p.cc ============================
	char	RtnDmsMode(){ return dms_mode; }
	UInt4	*RtnWtCnts(Int4 blk, Int4 i){ 	// may not need this...
			Int4 nblks=nBlksCMSA(CMA);
			assert(blk > 0  && blk <= nblks);
			assert(i > 0  && i <= LengthCMSA(blk,CMA));
			return WtCnts[blk][i];
		}
	char	*AlignP2P(ssx_typ *that,Int4 &start,Int4 &Score);
	double	ColScoreP2P(ssx_typ *that,Int4 i,Int4 j);
private:
	char	dms_mode;
	void	CalcSMX(double temp=0)
		   { if(ModeSMX=='W') CalcWtSMX(temp); else OldSMXs(temp); ReCalcSMX=FALSE; }
	void	CalcWtSMX(double temp=0);
	double	RelMap(char mode='S');	// OLD MAP.
	BooLean	ReCalcSMX;
	static const char ModeSMX='W';	// New, weighted matrix and MAP.
	// static const char ModeSMX='O';	// Old matrix and MAP.
        //======================== ssx_init.cc ============================
	void	Init(double pn,cma_typ);
	void	init(double pn,cma_typ);
	void	Free();
	void	InitCnts();
	void    ComputeTotalWtCnts();
	void	SeqWts();
public:
	void	SetTempSqWt(Int4 sq){ assert(sq > 0 && sq <= N && Wt[sq]==0); Wt[sq]=1; }
	void	UnSetTempSqWt(Int4 sq){ assert(sq > 0 && sq <= N && Wt[sq]==1); Wt[sq]=0; }
	// void	ReInit(){ }
	void	InitNDL(double Temp){ 
		   char mode='S'; if(Temp <= 0.0) mode='M';
		   if(ReCalcSMX) CalcSMX(Temp); 
		   NDL->IndelPenalty(stderr,mode,Wt); }
        //======================== ssx_init.cc ============================
private:
        //======================== ssx_typ.cc ============================
	Int4    TempToWt(double Temp);
	void    CntsFromSMX();
	smx_typ *StraightSMX( );
	smx_typ *WtCntsSMX( );
public:
	smx_typ *SampleSMX(double Temp);
	smx_typ *SampleDirichletSMX(double temp);
	UInt4	RtnSeqWt(Int4 sq) { assert(sq > 0 && sq <= NumSeqsCMSA(CMA)); return Wt[sq]; }
	double	RtnWtNumSeqs() {
		    Int4 total,sq; 
		    for(total=0,sq=1; sq <= NumSeqsCMSA(CMA); sq++) total+=(Int4) Wt[sq];
		    return ((double)total/(double)wt_factor); 
		}
	void	OldSMXs(double temp);
private:

	// temporary routines for testing...
	void    SampleSmatrix(double pernats, Int4 blk, Int4 wt);
        //======================== ssx_typ.cc ============================
	cma_typ	CMA;
	smx_typ	*SMX;
	ndl_typ *NDL;		// Jun Liu's HMM.
	dms_typ	*DMS;		// Dirichlet mixture scores.
	a_type	AB;
	Int4	N,totlen;
	UInt4	WtN,VrtlN;	// total weighted counts...
	double	pernats;
	char	*Wt;		// integer sequence weights; range = 1..100. 
	static const char wt_factor=100;
	static const UInt4 UNDERFLOW_TRIGGER=UINT4_MAX-(UInt4)wt_factor;
	lgm_typ	*LGM;
	UInt4	***WtCnts;	// Weighted residue counts at each position WtCnts[blk][i][r].
	UInt4	*TtlWtCnts;	// Total of weighted residue counts for each residue over all sequence.
	double	*TtlWtFreq;	// Total weighted frequency of residue counts.
	UInt4	***VrtlCnts;	// VrtlCnts[blk][s][r]: virtual weighted residue counts at each position.
	UInt4	**TtlVrtlCnts;	// Total of virtual weighted residue counts at each position.
	static const double pseudo=0.05;

        //======================== ssx_scores.cc ============================
#define testing_fix_for_DMS 1
        // void    mcalc(double *back, int **pam, double **M);
        // double  **mcalc(char **pam);
        double  **mcalc(Int4 **pam);
	void	calc_lngamma( );
	static const UInt4 WtPseudo=100;
	UInt4	TotalWtPseudo;	// Total Pseudocounts.
        // void    psiscore(double *count,double *score);
        void    psiscore(UInt4 *count,double *score);
        double  **Matrix;
        double  *BackGrnd;
        static const double     Pseudo=11.0;
        // static const double  Scale=693.147;  // 2^scale = e^1;   1/693-th bit = 1/1000-th nat.
        // use pernats...don't need scale.
        //======================== ssx_scores.cc ============================

        //======================== ssx_junk.cc ============================
	void	init_freq();
	void	ResFreqs();
	double	***Freq,**TotCnts,*tFreq;
	// ssx_typ(const ssx_typ *ssx);	// copy constructor.
};

#endif

