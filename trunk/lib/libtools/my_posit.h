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

#ifndef __MY_POSIT__
#define __MY_POSIT__
// Taken from posit.h by Alejandro Schaffer
// Contents: header file for position-based GBLAST.
/* $Revision: 6.8 $ * * $Log: posit.h,v $ */
/* Revision 6.8  1999/03/21 19:40:30  madden */

#include "my_ncbi.h"
#include "my_ncbimath.h"
#include "my_blastkar.h"
#include "brm_typ.h"
#include "gpxdrop.h"
#include "bsb_typ.h"
#include "alphabet.h"
#include "histogram.h"
#include "set_typ.h"	// fast convergence checking...

#define UNUSED (-1)

typedef struct {
  signed char	letter;  	// what is the preferred letter here
  Boolean	used;  		// is there any letter here 
  double	e_value; 	// score of highest hsp including this position 
  Int4		leftExtent; 	// How far left do same sequences match?
  Int4		rightExtent; 	// How far right do same sequences match?
} pde_typ; 	// == posDesc == description of position 

typedef struct {
  double	*posWtCount;	// AFN Addition: wieghted counts at each position.
				// Used for computing binomial tails.
  Int4		*posCount; 	// count of how many sequences match at
                  		// each query position, default value is 
                  		// 1 to include query.
  Int4		**posC; 	// position-specific occurrence counts
  double	**posMatchWeights;	// AFN: weight of each residue type 
					// at each position in the query.
  Int4		**posMatrix;
  Int4		**posPrivateMatrix;
  double	**posFreqs;	// AFN: position-specific freqs with pseudocounts.
  Int4		posNumSequences;
  double	*posA;
  double	*posRowSigma;
  Int4		pde_typMatrixLength;	// Length of pde_typMatrix, for deallocation.
  pde_typ	**pde_typMatrix;
  pde_typ	*posExtents;
  double	*posSigma;
  Int4		*posIntervalSizes;  	// interval size used for this column
  Boolean	*posUseSequences;
  double	*posInformation;
  set_typ	HitSet;			// set containing database hits.
  a_type	AB;			// AFN alphabet pointer.
} psim_type,*psim_typ;  // == posSearchItems

typedef struct {
  Uint1Ptr	query;
  Int4		qlength,alphabetSize,pseudoCountConst;
  double	lambda_ideal,K_ideal,ethresh,lambda,*standardProb;
  Int4 		**matrix;
  kbp_typ	*kbp_std, *kbp_psi, *kbp_gap_std, *kbp_gap_psi;
} csi_typ;  // == compactSearchItems
  

psim_typ MakePSIM(Int4 size_dbs, a_type AB);
void    NilPSIM(psim_typ posSearch);
csi_typ *compactSearchNew(bsb_typ search);
void	compactSearchDestruct(csi_typ *compactSearch);

void	outputPosMatrix(FILE *fp,psim_typ posSearch,csi_typ *compactSearch,a_type AB);
Int4	**CposComputation(psim_typ posSearch,bsb_typ search, csi_typ * compactSearch, 
	  sap_typ listOfGSeqAligns, Char *ckptFileName);
Int4Ptr * CposComputation(psim_typ posSearch, bsb_typ search,
        csi_typ *compactSearch, sap_typ listOfGSeqAligns, FILE *chkfp);
void	posPrintInformation(psim_typ posSearch, 
	  bsb_typ search, Int4 passNum);
void	posInitializeInformation(psim_typ posSearch, bsb_typ search);
void	posFreeInformation(psim_typ posSearch);
void	posConvergencedTest(psim_typ posSearch, bsb_typ search,
	  sap_typ listOfGSeqAligns);
	// Cleanup position-specific  data structures after one pass.
void	posCleanup(psim_typ posSearch, csi_typ *compactSearch);

Boolean posTakeCheckpoint(psim_typ posSearch,csi_typ * compactSearch, CharPtr fileName);
Boolean posReadCheckpoint(psim_typ posSearch,csi_typ * compactSearch, CharPtr fileName);
#if 1	// new, multifile read and write routines (afn: 12_8_08)
Boolean posTakeCheckpoint(psim_typ posSearch,csi_typ * compactSearch, FILE *checkFile);
Boolean posReadCheckpoint(psim_typ posSearch,csi_typ * compactSearch, FILE *checkFile);
#endif
void	posCheckpointFreeMemory(psim_typ posSearch, Int4 querySize);

void	posFreqsToMatrix(psim_typ posSearch,csi_typ *compactSearch);

//**************************** MARG PROB ****************************
Int4    **CposComputation(psim_typ posSearch, bsb_typ search, csi_typ *compactSearch,
        sap_typ listOfGSeqAligns, Char *ckptFileName,double ***MatrixMP);
void    posCalcMatchWeights(double *MatrixMP,double **posMatchWeights,
        Int4 qplace, double posA,a_type A);
void    posComputeSequenceWeights(psim_typ posSearch, csi_typ *compactSearch,
        double ***MatrixMP);
//**************************** MARG PROB ****************************

#endif /* __MY_POSIT__ */

