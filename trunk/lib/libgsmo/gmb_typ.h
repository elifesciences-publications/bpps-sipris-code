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

#if !defined (GAMBIT)
#define	GAMBIT

// #include <random>
#include "cmsa.h"
#include "residues.h"
#include "histogram.h"
#include "gibbs.h"
#include "msaheap.h"
#include "sequence.h"
#include "smatrix.h"
#include "cluster.h"
#include "ssx_typ.h"
#include "jlh_typ.h"
#include "my_ncbi.h"
#include "blosum62.h"
#include "cma_gmb.h"
#include "str_typ.h"
#include "rst_typ.h"

/*************************** ADT GAMBIT ***************************
Gambit:	"A chess 'opening' in which a player risks one or more
        minor pieces to gain an advantage in position."

    here: Sequence 'openings'in which one or more residues are
          inserted or deleted in order to gain a statistical advantage in 
          multiple alignment sequence positioning.

	GAMBIT (Gapped Alignment with MCMC-Based Indel Tempering)

Tempered: "Having the elements mixed in satisfying proportions" (TEMPERATE).

Tempering: "Mixing the elements in satisfying proportions" (TEMPERATE).

/********************************************************************/
class gmb_typ {
public: 
		gmb_typ(Int4 pio,Int4 pdo,Int4 eie,Int4 ede,cma_typ cma,
			   char dms_md=' ',double prr_wt=1.0, ssx_typ *issx=0)
			{ init(pio,pdo,eie,ede,cma,dms_md,prr_wt,issx); }
		~gmb_typ( ){ Free(); }
	BooLean	ChangePriors(Int4 pio,Int4 pdo,Int4 eie,Int4 ede)
		{ SSX->ChangePriors(pio,pdo,eie,ede); }
	void	SetPriorWt(double D){ SSX->SetPriorWt(D); }
	void	SetSqWtAdjust(double D){ SSX->SetSqWtAdjust(D); }
	double	ScoreBILD(Int4 blk, Int4 i){ return SSX->BildScore(blk,i); }
	double	RtnWtNumSeqs(){ return SSX->RtnWtNumSeqs();}
	void	RestoreBest(){ SSX->RestoreBest(); }
	//======================== gmb_put.cc ============================
	void	SampleAlignSeq(FILE *fp,char mode,e_type qE, double temp=300.0)
		{ SSX->SampleAlignSeq(fp,mode,qE,temp); }
	void	AlignSeq(FILE *fp,char mode,e_type qE, double temp=300.0);
	char	*AlignSeq(e_type qE,Int4 &start,double temp=0.0)
			{ return SSX->AlignSeq(qE,start,temp); }
	char	*GlobalAlignSeq(e_type qE,Int4 &start,double temp=0.0)
			{ return SSX->AlignSeq(qE,start,temp,TRUE); }
	cma_typ CreateAlign(ss_type data);
	cma_typ CreateFullCMSA(ss_type Data, set_typ Set);
	cma_typ CreatePurgedCMSA(Int4 percent_ident,e_type *Seqs, set_typ &Set, ss_type &data);
	//======================== gmb_smpl.cc ============================
	cma_typ Sample(char *arg,char c, Int4 sim, double startTemp=300.0,FILE *rfp=0);
	cma_typ Sample(char *name, char mode,Int4 similarity,double startTemp,
                	double endTemp, double incTemp,FILE *rfp=0);
	void	SetMaxIter(Int4 value){ assert(value > 0); MaxIter=value; }
	//======================= hieraln.cc ==============================
	cma_typ Optimize();
	BooLean AlignP2P(FILE *fp,cma_typ child);
	char    *AlignP2P_Operation(cma_typ child, Int4 &start);
	ssx_typ	*RtnSSX(){ return SSX; }
	double	RtnMap(){ return SSX->GapMap( ); }
	void	InitHMM(){ SSX->InitNDL(0.0); }
	//======================== gmb_sticky.cc ============================
	double  SimilarTogether(Int4 MinSticky, double MaxFrctn,double Temp);
	BooLean SampleTogether(set_typ Set,double Temp,double &MAP,dh_type dH);
	set_typ	GambitTogether(Int4 II, Int4 MinSticky, double MaxFrctn, double temp,
			set_typ LstSet,Int4 min_ins=1);
	cma_typ SampleStickyTogether(FILE *fp,Int4 MinSticky, double MaxFrctn,
						char *name,double Temp,char mode);
	cma_typ	foo(double *f){ printf("%.2f\n",f[0]); printf("%.2f\n",f[1]); printf("%.2f\n",f[6]); } 
	cma_typ SampleStickyTogether(FILE *fp,Int4 MinSticky, char m[20], double f[20],
				double T=0.0, char mode='R');
	double  RandomTogether(Int4 MinSticky, double MaxFrctn, double temp);
	double  WorstTogether(Int4 MinSticky, double MaxFrctn, double temp);
	double	ConservedTogether(double MaxFrctn, double Temp);
	BooLean	SamplePurged(Int4 percent_id, double MaxFrctn,double Temp,
			double &sm, double &em, Int4 stage=1);
	BooLean SampleByLayers(Int4 MinSize,double Temp,double &strt_map,double &end_map,Int4 stage);
	void	SetExcludeSet(set_typ set){ SetEx=set; }
	void	UnSetExcludeSet( ){ SetEx=0; }
private:
	set_typ	SetEx;
public:
	//======================== gmb_sticky.cc ============================

	//======================== gmb_cols.cc ============================
	cma_typ AddColumns(double bild_cut=0.0,BooLean Widen=TRUE,set_typ InSet=0,
			double MinMatchFrq=0.80,BooLean	EndsOnly=FALSE);
	cma_typ RmColumns(double bild_cut=0.0, BooLean Shrink=TRUE,double MinDelFrq=0.0);
	cma_typ RmShortBlks(double bild_cut=0.0, BooLean Shrink=TRUE);
	cma_typ	RtnCMSA(){ return SSX->RtnCMA();}
	char	*RtnAddOp(){ char *rtn=AddOp; AddOp=0; return rtn; }
	void	SetAddOp(char *op){ if(AddOp) free(AddOp); AddOp=AllocString(op); }
private:
	char	*AddOp;	// Operations need to revert back to orginal columns after AddColumns();
	char    *FindBlocks(ssx_typ *ssx, char mode, double cutoff =-1.0, double bild_cut=0.0, Int4 Limit=3);
	//======================== gmb_cols.cc ============================

	Int4	aa_per_io,aa_per_do,exp_ie, exp_de;
	char	dms_mode;
	//======================== gmb_operate.cc ============================
	// column and block operations for gambit... move most of this to cma_typ...modify cma routines.
	BooLean SlideColLeft(Int4 t, double *oldMap, double *newMap,cma_typ L, double temperature);
	BooLean SlideColRight(Int4 t, double *oldMap, double *newMap,cma_typ L, double temperature);

	//======================== gmb_operate.cc ============================

	//======================== gmb_typ.cc ============================
	void	init(Int4 eio,Int4 edo,Int4 eie,Int4 ede, cma_typ cma,char dm,
			double prr_wt,ssx_typ *issx=0);
	void	Free( );
	void    DeleteSSX(ssx_typ *ssx);
	//======================== gmb_typ.cc ============================
	ssx_typ	*SSX;
	Int4	MaxIter;
	double	pernats;
	a_type	AB;
	char	str[205];

	//======================== gmb_smpl.cc ============================
	double  GambitSingles(char mode,double Temperature);
public:
	BooLean SampleSingle(Int4 s,double Temperature,double &map);
private:
	double	GambitClusters(char mode,Int4 percent_ident,double Temperature);
	BooLean SampleCluster(Int4 *cluster,Int4 *offset, e_type csq,double Temperature,double &MAP);
	BooLean SampleTransition(double nmap, double omap, double temp);
	//======================== gmb_smpl.cc ============================

	//======================== gmb_clocl.cc ============================
	double  ClustersOfClustersSampling(char mode, double Temp,Int4 percent_ident=60);
	BooLean SampleClustersOfClusters(set_typ ClstSet,Int4 **cluster, Int4 NumClst, Int4 **offset,
							double Temp,double &MAP);
	//======================== gmb_clocl.cc ============================

	//======================== gmb_clstr.cc ============================
	char    *ConvertOperation(char *operation, e_type sE, Int4 mst, Int4 offset);
	e_type  MultiAlnToConsensus(FILE *fp, Int4 *OffSet, Int4 *ListE);
	Int4    **GetClusters(FILE *fp, Int4 **&OffSet, Int4 &Nset);
	e_type	DelimitSeqAln(Int4 i, Int4 &strt, Int4 &end, Int4 &len, Int4 &del);
	Int4    *MultiAlnToLongest(FILE *fp, Int4 *ListE);
	Int4    Extend(register unsigned char  *p1,register unsigned char  *p2,
                	register Int4 len, register char **R, Int4 dropoff);
	class	ops_typ {	// operation type...
public:
		  static const Int4 size=100;
		  ops_typ( ){ this->N=0; }
		  Int4 Add(char *op,Int4 strt,Int4 scr){
			assert(N < size); N++;
			this->operation[N]=op;
			this->trace_len[N]=strlen(op);
			this->start[N]=strt;
			this->score[N]=scr;
		  }
		  void Sort( ){	// using bubble sort algorithm...
			Int4 i,j,t,n;
			char *op;
			for(i=2; i <= N; i++){
			   for(j=i; j > 1 && (score[j-1] < score[j]); j--){
                               t=score[j]; score[j]=score[j-1]; score[j-1]=t;
                               t=start[j]; start[j]=start[j-1]; start[j-1]=t;
                               t=trace_len[j]; trace_len[j]=trace_len[j-1]; trace_len[j-1]=t;
                               op=operation[j]; operation[j]=operation[j-1]; operation[j-1]=op;
                	   }
			}

		       }
		  void Put(FILE *fp)
		     { for(Int4 i=1; i <= N; i++) Put(fp,i); fprintf(fp,"\n"); }
		  void Put(FILE *fp,Int4 i){
			if(i > 0 && i < size){
			  fprintf(fp,"%d: start=%d; score=%d; op=%s\n",
				i,start[i],score[i],operation[i]);
			}
		  }
		  ~ops_typ( )
		     { for(Int4 i=1; i <= N; i++) free(operation[i]); }
		  char *operation[size];
		  Int4 N,start[size], score[size], trace_len[size];
	};
private:
	ops_typ	*OPS;
	Int4    AlignLenOperation(char *Op);
	Int4	GetMultipleTraces(e_type csqE,double Temp);
	Int4    ScoreSeqAlnSMatrix(char *operation, e_type sbjE, Int4 strt);
	//======================== gmb_clstr.cc ============================
	Int4	**CsqCnts;
	static const Int4 MinNumClstrs=25;
	e_type	*CsqQuery;

	//======================== gmb_subaln.cc ============================
public:
	// void    PutSampleOutSqIn(FILE *fp,double Temp,set_typ Set,set_typ tstSet);
	void    PutSampleOutSqIn(FILE *fp,double Temp,set_typ Set,Int4 *tstSet);
	BooLean	SampleOutSqIn(Int4 sq,double Temp,double &MAP,e_type CsqE=0,Int4 min=INT4_MAX);
	// set_typ RtnTestSet(set_typ Set);
	Int4	*RtnTestSet(set_typ Set);
private:
	//======================== gmb_subaln.cc ============================

	//======================== gmb_debug.cc ============================
public:
	void    TestAddRemove();
	void	PutGappedAln(FILE *fp,e_type sbjE,char *operation, Int4 start);
private:
	void	PutDebugClocl(Int4 strt,Int4 os,Int4 sq,char *op, e_type sbjE);
	Int4    PutMultiAlign(FILE *fp, Int4 *ListE,Int4 *offset,e_type cE);
	sap_typ	MultiAlign(Int4 *ListE,Int4 *offset,e_type cE);
	Int4	*PutAlnToLongest(FILE *fp, Int4 **ListE, Int4 s);
	Int4    GetDiagonalEnds(Int4 os,e_type E1,e_type E2,Int4 &S1,Int4 &S2,Int4 &e1,Int4 &e2);
	void    PutListAsAln(FILE *fp, Int4 **ListE, Int4 *Strt, Int4 *Lngth, Int4 s);
	Int4    PutDiagonalSeq(FILE *fptr, Int4 offset, e_type E1, e_type E2);
	void    PutClusterAln(char *name,Int4 *cluster, cma_typ cma);
	void    PutMasterSlaveAln(Int4 s, e_type qE, Int4 start, char *qOp,
                                e_type E, Int4 strt, char *op,Int4 os, cma_typ cma);
	void    PutOperations(FILE *fp,Int4 ii, Int4 s, char *operation, char *op);
	Int4    MaxSegSMatrix(unsigned char *maxseg, smx_typ M);
	Int4	PutGappedSeqAlnSMatrix(FILE *fp, char *operation, Int4 offset, Int4 n2,
                        unsigned char *seq2, Int4 nmod, smx_typ *M);
	void    PutClusterIDs(FILE *fp,Int4 n, Int4 *cluster);
	void	PutCMA_SqAln(FILE *fp, Int4 sq);
	void    PutCMA_SqIndels(FILE *fp, Int4 sq);
	void    PutCMA_SqIndels(FILE *fp, Int4 sq, Int4 *pos, gsq_typ *gsq);
	//======================== gmb_debug.cc ============================

	//======================== gmb_junk.cc ============================
	Int4	**ClusterGPSI_CMSA(char method, double cutoff, Int4 &Nset);
	Int4	**ClusterCMSA(Int4 percent,Int4 *Nset);
	Int4	**GetRepSetCMSA2(Int4 similarity,Int4 *Nset,cma_typ cma);
	// void	ReplaceSSX(ssx_typ *ssx); // Replace SSX and CMA file;
	// BooLean SampleSingle0(Int4 s,double Temperature,double &map);
	// BooLean SampleCluster0(Int4 *cluster,Int4 *offset, e_type csq,double Temperature,double &MAP);
	// cma_typ ProgressiveSampling(char *,Int4 , Int4, Int4);
	//======================== gmb_junk.cc ============================
};

#endif

