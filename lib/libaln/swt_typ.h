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

#if !defined(_SWT_TYP_)
#define _SWT_TYP_

#include "cmsa.h"
#include "dheap.h"
#include "alphabet.h"
#include "histogram.h"
#include "random.h"
#include "set_typ.h"
#include "hsw_typ.h"

class swt_typ {         // sequence weighting type
public:
                swt_typ( ){ assert(!"Illegal constructor"); }
		swt_typ(cma_typ,cma_typ,BooLean,BooLean);
		swt_typ(cma_typ,cma_typ,BooLean,BooLean,char *);
		swt_typ(hsw_typ h){ Init(h); }
		swt_typ(cma_typ,BooLean);	// constructor for creating & writing HSW.
		swt_typ(cma_typ,cma_typ,BooLean,e_type,BooLean);
		swt_typ(cma_typ,cma_typ,BooLean,e_type,BooLean,char *);
                ~swt_typ( ){ Free( ); }
	void	Put(FILE *);
	cma_typ	RtnHSW_CMA( ){ if(hsw) return hsw->mcma; else return 0; }
	hsw_typ	RtnHSW( ){
		   // assert(own_hsw);	// don't need to require this...
		   if(own_hsw){	
		      if(compute_wts) HenikoffWeights( );
		      if(rtn_hsw==0){
			NEW(rtn_hsw,1,henikoff_wts_type);
			rtn_hsw->Weight=Weight;
			rtn_hsw->WtCnts=WtCnts;
			rtn_hsw->WtFreq=WtFreq;	
			rtn_hsw->fract_seq_aln=fract_seq_aln;
			rtn_hsw->Length=Length;
			rtn_hsw->NWtSq=NWtSq;
			rtn_hsw->mcma=mcma;
			rtn_hsw->AB=AB;
		      } 
		   } return rtn_hsw; // else simply return the weights.
		 }
	void	Test(FILE *);
	void	FWrite(char *name){
		   assert(own_hsw); if(compute_wts) HenikoffWeights( );
		   FILE *fp = open_file(name,".swt","w");
		   FWriteHenikoffWeights(fp); fclose(fp);
		}
	double	**Weights( ){ if(compute_wts) HenikoffWeights( ); return Weight; }
	double	**WeightedFreq( ){ if(compute_wts) HenikoffWeights( ); return WtFreq; }
	double	*MargProbAln( ){
		   assert(own_hsw); if(compute_mp) ComputeMargProb( );
		   return MargProb;
		}
	double	**MargProbWtFreq( ){
		   assert(own_hsw); if(compute_mp) ComputeMargProb( );
		   return MPWtFreq;
		}
	double	**ObsWtCnts( ){ if(compute_wts) HenikoffWeights( ); return WtCnts; }
	double	*FractSeqAln( ){ if(compute_wts) HenikoffWeights( ); return fract_seq_aln; }
	void	Silent() { Verbose=FALSE; }
	UInt4   *GetIntegerWts(unsigned char ***RtnSqWt);
	UInt4	WtFactor(){ return 100; }
	Int4	TheLength(){ return Length; }
	Int4	NumWtSeqs(){ return NWtSq; }
	BooLean	OwnHSW(){ return own_hsw; }
private:
	void	Init(cma_typ,cma_typ,BooLean,e_type,BooLean,char *);
	void	Init(hsw_typ);
        void    init(cma_typ);
	void	Put(FILE *,Int4);	// use only for debugging...
	BooLean	Verbose;
	//******************** HenikoffWeights information: ***********************
	hsw_typ	rtn_hsw;		// use for passing sequence weights in & out.
	hsw_typ	hsw;		// passed in data on sequence weights (do not own this).
	BooLean	own_hsw;	// should these be freed here or from a calling environment?
	void	HenikoffWeights( );
	void	FWriteHenikoffWeights(FILE *fp);
	void	FReadHenikoffWeights(FILE *fp);
	char	*WtsFileName;	// Filename for precomputed HenikoffWeights( ).
	double	**Weight;	// position specific sequence weights
	double	**WtCnts;	// Weighted counts at each position.
	double	**WtFreq;	
	double	*fract_seq_aln;
	//------------------- temporary arrays NOT NEEDED FOR BINARY READ/WRITE. ----------------
	Int4	*numParticipants;	// temporary array of length LenSeq(csq);
	Int4	*info_order;	// temporary array for ordering s; NOT NEEDED FOR BINARY READ.
	Int4	*order;		// temporary array for ordering s
	Int4	**Observed;	// temporary (s x r) array.
	double	*info;		// temporary array for relative entropy at s
	unsigned char	*nResTyp;
	//******************** end HenikoffWeights information: ***********************
	BooLean	SortColumns(Int4);
	double  relative_entropy(Int4);
	void	ComputeMargProb(Int4 pernats,Int4 gap_open,Int4 gap_extend);
	void	ComputeMargProb(){ ComputeMargProb(100,347,35); }
	double	*MargProb;
        void    Free( );
	double	scale;
        cma_typ ocma,mcma;
	a_type	AB;
	double	minPart,maxPart;
	double	mininfo,maxinfo;
	unsigned char **PtrWtSq;// pointer to FakeSeqs in mcma.
	double	*freq;		// background frequency
	double	**ofreq;	// ortholog frequency
	double	**oObs;		// ortholog Observed residues
        e_type  csq;            // consensus sequence used as query for cma
        set_typ *Participant;   // sets of participating sequences at each site
	set_typ	**RawCnts;	// RawCnts[s][r] = residues r in position s
	set_typ	*IS,is;		// sets for taking intersections.
	double	**MPWtFreq;	
	Int4	NWtSq;
	Int4	Length;
	double	UpperLimitInfo;
	unsigned char StartAlpha;
	BooLean	compute_wts,compute_mp,use_ocma_pseudo;	
	Int4	rm_high_info,keep_many_partic;
	// Speed up attempt...
	void    SumSeqWts(BooLean *skip,Int4 S,double minfract,Int4 to_use);
};

#endif

