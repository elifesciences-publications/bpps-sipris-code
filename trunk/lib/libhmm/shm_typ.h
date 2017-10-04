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

#if !defined(_SHM_TYP_)
#define _SHM_TYP_

#include "idp_typ.h"
#include "stpb_typ.h"
#include "pfam.h"
#include "my_ncbi.h"

class shm_typ {
public:
		shm_typ(){ name=NULL; desc = NULL; nule = NULL; mat_emit=NULL; ins_emit = NULL; stpb = NULL; }
		shm_typ(char *Name, char *Desc, Int4 *Nule, Int4 **Matemit, Int4 **Insemit, 
			Int4 len, idp_typ *idp, double FreqOfMatchAtPositionOne,a_type A);
		shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, 
			stpb_typ *Stpb, a_type A);
		shm_typ(char *Name, char *Desc, double *Nule, double **Matemit, double **Insemit, 
			stpb_typ *Stpb, a_type A);
		shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, 
			stpb_typ *Stpb, a_type A, double Lambda, double U, double SLambda, double SU, 
			Int4 Wthreshold, Int4 Uxparameter, Int4 Uthreshold, Int4 Xparameter, 
			Int4 Gthreshold, Int4 Vthreshold, e_type Cons);
		shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, 
					a_type A, Int4 Rpts);
		shm_typ(char file_type, FILE *fp, a_type A);
		~shm_typ(){ Free(); }
	void    Put(FILE *);
	char    *Align(e_type,Int4 *,Int4 *,Int4 *);
	char    *loc_Align(e_type,Int4 *,Int4 *,Int4 *,Int4 *);
	char    *loc_AlignWOTB(e_type);
	void 	Write(FILE *fp);
	shm_typ *ReverseSHM();
	shm_typ *ConcatenateSHM(Int4 rpts);
	e_type	GetConsensus() { 
			MkConsensus( );
			return consensus; 
		}
	void	SetConsensus(e_type E) { consensus = E; }
	unsigned char 	*ConsensusSeqPtr();
        Int4    *MatToMat() { return stpb->MatToMat(); }
        Int4    *MatToIns() { return stpb->MatToIns(); }
        Int4    *MatToDel() { return stpb->MatToDel(); }
        Int4    *InsToIns() { return stpb->InsToIns(); }
        Int4    *InsToMat() { return stpb->InsToMat(); }
        Int4    *DelToDel() { return stpb->DelToDel(); }
        Int4    *DelToMat() { return stpb->DelToMat(); }
        Int4    *BegToMat() { return stpb->BegToMat(); }
        Int4    *MatToEnd() { return stpb->MatToEnd(); }
        Int4    NB() { return stpb->NB(); }
        Int4    NN() { return stpb->NN(); }
        Int4    EC() { return stpb->EC(); }
        Int4    EJ() { return stpb->EJ(); }
        Int4    CT() { return stpb->CT(); }
        Int4    CC() { return stpb->CC(); }
        Int4    JB() { return stpb->JB(); }
        Int4    JJ() { return stpb->JJ(); }
        Int4    BM() { return stpb->BM(); }
        Int4    BD() { return stpb->BD(); }
	Int4	GG() { return stpb->GG(); }
	Int4    GF() { return stpb->GF(); }
	Int4 	MatToMat(Int4 position) { return stpb->MatToMat(position); }
	Int4 	MatToIns(Int4 position) { return stpb->MatToIns(position); }
	Int4 	MatToDel(Int4 position) { return stpb->MatToDel(position); }
	Int4 	InsToIns(Int4 position) { return stpb->InsToIns(position); }
	Int4 	InsToMat(Int4 position) { return stpb->InsToMat(position); }
	Int4 	DelToDel(Int4 position) { return stpb->DelToDel(position); }
	Int4 	DelToMat(Int4 position) { return stpb->DelToMat(position); }
	Int4 	BegToMat(Int4 position) { return stpb->BegToMat(position); }
	Int4 	MatToEnd(Int4 position) { return stpb->MatToEnd(position); }

	Int4    GetLength() { return shm_len; }
	Int4	Bits() { return bits; }
	Int4	*NULE() { return nule; }
	Int4    *MatEmit(Int4 position) { return mat_emit[position]; }
	Int4    *InsEmit(Int4 position) { return ins_emit[position]; }
	Int4 	**Matrix() { return mat_emit; }
	Int4    **InsMatrix() { return ins_emit; }
	void	SetMatrix(Int4 **mtx) { mat_emit = mtx; }
	void    SetInsMatrix(Int4 **imtx) { ins_emit = imtx; }
	Int4	WordScore(unsigned char *word, Int4 word_len, Int4 position);
	Int4	WordHits(unsigned char *word, Int4 word_len, 
			Int4 trashold, unsigned short *position);
	double  GetLambda() { return lambda; }
	void	SetLambda(double Lambda) { lambda = Lambda; }
	double  GetU() { return u; }
	void	SetU(double U) { u = U; }
	double  GetSLambda() { return Slambda; }
	void	SetSLambda(double Lambda) { Slambda = Lambda; }
	double  GetSU() { return Su; }
	void	SetSU(double U) { Su = U; }
	double  GetChi() { return chi; }
        Int4    GetWthreshold() { return wthreshold; }
	Int4    GetUxparameter() { return uxparameter; }
	Int4    GetUthreshold() { return uthreshold; }
	Int4    GetXparameter() { return xparameter; }
	Int4    GetGthreshold() { return gthreshold; }
	Int4    GetVthreshold() { return vthreshold; }
	void	SetWthreshold(Int4 Wthreshold) { wthreshold = Wthreshold; }
	void	SetUxparameter(Int4 Uxparameter) { uxparameter = Uxparameter; }
	void	SetUthreshold(Int4 Uthreshold) { uthreshold = Uthreshold; }
	void	SetXparameter(Int4 Xparameter) { xparameter = Xparameter; }
	void	SetGthreshold(Int4 Gthreshold) { gthreshold = Gthreshold; }
	void	SetVthreshold(Int4 Vthreshold) { vthreshold = Vthreshold; }
	a_type  GetAlphabet() { return AB; }
	void	SetAlphabet(a_type A) { AB = A; }
	Int4	Getf_open() { return f_open; }
	Int4    Getf_ox() { return f_ox; }
	char 	*GetName() { return name; }
	void	SetName(char *Name) { name = Name; }
	stpb_typ	*GetStpb() { return stpb; }
	char	*GetDesc() { return desc; }
	Int4	Repeats() { return rpts; }
	double 	**Mat_emit_Prob();
	double 	**Ins_emit_Prob();
	double 	*Mat_emit_Prob(Int4 pos);
	double 	*Ins_emit_Prob(Int4 pos);
	e_type 	*EmitSeqSHM(Int4 nseq, double emiss_tmpr, double trans_tmpr);
	void 	EmitSeqSHM(FILE *fp,Int4 nseq, double emiss_tmpr, double trans_tmpr);
	Int4	ScoreLocalSAP(e_type sE, sap_typ sap);
	shm_typ *RenormalizeSHM();	
	shm_typ *TempSHM(double temp, double temp_trans);
	double	ExpectedMainModelScore();
	double	ExpectedScore();
	double	ExpectedMatchStateScore(Int4 i);
	double	ExpectedMatchScore(Int4 i);
	double	ExpectedDeletionScore(Int4 i);
	double	ExpectedInsertionStateScore(Int4 i);
	double	ExpectedInsertionScore(Int4 i);
	double	ExpectedBeginScore();
	double	ExpectedCNJScore();
	double	ExpectedNScore();
	double	ExpectedJScore();
	double	ExpectedCScore();
//	Int4 	*RunningScore(dsp_typ dsp, Int4 *length);
	Int4 	*RunningScore(dsp_typ dsp, Int4 *length, char stype);
private:
	BooLean		reversed;
	double 		lambda,u,Slambda,Su,chi;
	void		readShmB(FILE *fp, a_type A);
	void		MkConsensus( );
	a_type 		AB;
	Int4		shm_len,rpts,bits,f_open,f_ox;
	Int4 		wthreshold,uxparameter,uthreshold;
	Int4 		xparameter,gthreshold,vthreshold;
	e_type		consensus;
	void		init(char *Name, char *Desc, Int4 *NULE, 
				Int4 **Mat_emit, Int4 **Ins_emit, 
				stpb_typ *Stpb, a_type A);
	void		init(char *Name, char *Desc, Int4 *NULE, 
				Int4 **Mat_emit, Int4 **Ins_emit, 
				stpb_typ *Stpb, a_type A, double Lambda, 
				double U, double SLambda, 
				double SU, Int4 Wthreshold, 
				Int4 Uxparameter, Int4 Uthreshold, 
				Int4 Xparameter, Int4 Gthreshold, 
				Int4 Vthreshold, e_type Cons);
	char		*Viterbi(e_type,Int4 *,Int4 *,Int4 *);
	char		*loc_Viterbi(e_type,Int4 *,Int4 *,Int4 *,Int4 *);
	char		*loc_ViterbiWOTB(e_type);
	void		Free();
	stpb_typ	*stpb;
        Int4		**mat_emit,**ins_emit;
	Int4 		*nule;
	Int4		**MAT,**DEL,**INS;
	char		*name,*desc,**traceM,**traceI,**traceD;
	void		inner_loop(Int4 j, Int4 *matj, Int4 *matjm1,
				Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
				Int4 m2m, Int4 d2m, Int4 i2m, Int4  m2d, 
				Int4  d2d, Int4  m2i, 
				Int4 i2i, Int4 seq_len, unsigned char *seq, 
				Int4 *scorej,Int4 *ins_emit);
	char		*get_traceback(Int4, Int4 *, Int4 *);
	char		*loc_get_traceback(Int4, Int4 *, Int4 *, Int4 *, Int4 *);
	void		loc_inner_loop(Int4 j, Int4 *matj, Int4 *matjm1,
				Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
				Int4 m2m, Int4 d2m, Int4 i2m, Int4  m2d, 
				Int4  d2d, Int4  m2i, Int4 i2i, Int4 seq_len, 
				unsigned char *seq, Int4 *scorej,Int4 *ins_emit);
};

shm_typ **ReadShm(FILE *fp, Int4 bits, Int4 *nmbr, a_type A);
shm_typ *ReadProbShm(FILE *fp, Int4 bits, a_type A);
void    WriteSHM(FILE *fp, shm_typ *shm, shm_typ *rshm);
sap_typ AlignToSAP(e_type qE, e_type sE, char *operation, Int4 start_qE, Int4 start_sE, Int4 sc);
sap_typ AlignToSAP2(e_type qE, e_type sE, char *operation, Int4 start_qE, Int4 start_sE, Int4 sc);

#endif

