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

#if !defined (_LPR_TYP_)
#define _LPR_TYP_

#include "sset.h"
#include "alphabet.h"
#include "dheap.h"
#include "set_typ.h"
#include "sequence.h"
#include "rst_typ.h"
#include "bpps_typ.h"
#include "cmsa.h"	// in libaln
#include "swt_typ.h"	// in libaln

class lpr_typ {	// log-probability ratio type for bpps (later put in for cma as well).
public:
	lpr_typ( ){ usage(); } 
        lpr_typ(cma_typ ,swt_typ *);
        lpr_typ(cma_typ in_cma,swt_typ *in_swt, BooLean);
        ~lpr_typ( ){ Free( ); }
        sst_typ *GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,BooLean B,double &L,Int4 x,double priorRi=0);
        sst_typ *GetOptPttrnLPR(FILE *fp,set_typ SetFG, set_typ SetBG,BooLean Negate,double &lpr,
                                        Int4 MaxCols,unsigned char *&rtn_csq,double priorRi=0);
        sst_typ *GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,BooLean B,double &L,Int4 x,char typ,
			e_type qE=0,double priorRi=0);
        sst_typ *GetOptPttrnLPR(FILE *fp,set_typ SetFG, set_typ SetBG,BooLean Negate, double &lpr,
                        Int4 MaxCols, unsigned char *&rtn_csq, char Type,e_type qE=0,double priorRi=0);
	void    PutParameters(FILE *fp,set_typ SetFG, set_typ SetBG,sst_typ *qsst, BooLean Negate,char Type,
			double priorRi=0);
	void    CompareBPPS(FILE *fp,bpps_typ *PPS,set_typ FG, set_typ BG,sst_typ *qsst,
			BooLean Negate,char Type,double priorRi=0);
	void	CompareSqWts(FILE *fp,set_typ SetFG,set_typ SetBG,UInt4 **WtSqFG,UInt4 **WtSqBG,BooLean Neg);
        double  WtCardFG_BG_Sets(double &WtCntsFG, double &WtCntsBG);
        double  CalcSetvsPttrnLPR(FILE *fp,set_typ FG, set_typ BG,sst_typ *qsst,BooLean neg,
				e_type qE=0,double priorRi=0);
        double  CalcSetvsPttrnLPR(FILE *fp,set_typ FG, set_typ BG,sst_typ *qsst,BooLean neg,char typ,
				e_type qE=0,double priorRi=0);
private: 
	unsigned char del_as_random;
	void	Init();
	void	Free();
	void	FreeBPPS(bpps_typ *pps);
	void	usage();
        void    InitCntsFG_BG(set_typ SetFG, set_typ SetBG,BooLean Negate,BooLean UseGlobalSqWts);
	bpps_typ *MakeBPPS(sst_typ *qsst,BooLean UseGlobalSqWts,char Type, double priorRi=0,
			unsigned char *Qry=0);

	Int4    Length;         // Alignment length.
	cma_typ cma;            // Input alignment.
        rst_typ *RST;
        sst_typ **LegalSST;
	BooLean	UseGlobalSqWts;
        unsigned char   **SqWt;
        UInt4	**CntBG,**CntFG;
        UInt4   *AveSqIWt;
        swt_typ *swt;
	a_type	AB;
};

#endif


