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

#ifndef _PSG_H_
#define _PSG_H_
#include "cmsa.h"

class psg_typ{

public:
			psg_typ();
			~psg_typ();
			psg_typ(cma_typ , Int4 );
			psg_typ(cma_typ , char , char , double *, Int4 );
			psg_typ(cma_typ , char , char , double *, Int4 , char );
	smx_typ		*BlcksMatrices(Int4 nRpts);
	double 		*RelEntropy(Int4 Rpts);
	smx_typ 	*SampleBlcksMatrices(Int4 nRpts, Int4 wt);
	Int4		Nblks() { return nBlks; }
	Int4 		Nseq() {return nSeq; }
	Int4 		Ncol() { return nCol; }
	double          **TargetFreq();	
	Int4		**LogOddMatrix();
	Int4		**LogOddMatrix(Int4);
	double		*HenikoffWeights();
	double		**HenikoffWeightedCounts();
private:
	void		Free();
	void		Mkpsg_typ(cma_typ cma, char w, char p, double *m, Int4 pn);
	double 		**WeightedCounts();
	double 		**PseudoCounts();
	double 		**HenikoffPseudoCounts();
	double          **OnePseudoCount();
	Int4 		*RelEnPosit(double cutoff, Int4 Rpts, Int4 *length);
	char		weightMeth, psdcntMeth;
	double 		*h_mult;
	a_type 		A;
	ss_type		SeqSet;
	Int4		**blPos, pernats;
	Int4 		nSeq, nBlks, nCol;
	Int4		*blLen, *nResTyp, **rawCnts;
	double		**weightedCnts, **psdCnts;
};
#endif
