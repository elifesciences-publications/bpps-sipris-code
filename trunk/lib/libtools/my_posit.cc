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

/*****************************************************************************
File name: posit.c   Author: Alejandro Schaffer
Contents: utilities for position-based GBLAST.
*****************************************************************************/
/* $Revision: 6.25 $ */ /* $Log: posit.c,v $ */
/* Revision 6.25  1999/04/05 14:45:40  madden */

#include "my_posit.h"

/* Front-ends to retrieve numbers. */
#define  getCkptDouble(d, ckptFile)  	(getCkptNumber(&(d),sizeof(double),ckptFile))
#define  getCkptInt4(i, ckptFile)       (getCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  getCkptChar(c, ckptFile)       (getCkptNumber(&(c),sizeof(Char),ckptFile))

/*small constants to test against 0*/
#define posEpsilon 0.0001
#define posEpsilon2 0.0000001

#if 0	// OLD with problems....
#define XChar	21	// Alejandro's mapping...  Xchar == 0 changed to AFN's mapping
// Note that Xchar != XChar; Xchar was  21 but AFN changed to 0 (see brm_typ.h)
#define GAP_CHAR 0 	// Representation of a gap in a motif
// #define GAP_CHAR (21) 	// old AFN Representation of a gap in a motif
// WARNING: GAP_CHAR and Xchar values have been swapped!!!
//     MAY NEED TO CHANGE THIS TO BRING INTO CONSISTENCY WITH NCBI PSI-BLAST  -AFN.
#endif

#if 1	// This seems to fix problem...
#define XChar	 0	// AFN: My mapping
#define GAP_CHAR 21 	// Representation of a gap in a motif
#endif
#define GAP_HERE (-1)	// Used inside a seqAlign to reprsent the presence of a gap

/*Used to check that diagnostics printing routine will work*/
#define POSIT_PERCENT 0.05
#define POSIT_NUM_ITERATIONS 10

#define POSIT_SCALE_FACTOR 1000

#define POS_RESTING 0
#define POS_COUNTING 1

#define IDENTITY_RATIO 0.98

static void	FixposMatrix(psim_typ posSearch, csi_typ *compactSearch)
// AFN: 11/7/01 fix for 'x' residues
// WARNING: I still don't understand what is going on here! 
// Need to think this through some more...!!!
{
   Int4 c,length = compactSearch->qlength;
   for(c=0; c < length; c++) {
     if(posSearch->posMatrix[c][XChar] == GBLAST_SCORE_MIN)
		posSearch->posMatrix[c][XChar] = -2; 
   }
}

// Allocate memory for  data structures inside posSearch used in
// position-specific caculations
// posSearch -- to be filled in 
// alphabetSize -- number of distinct characters used in the sequences
// querySize -- number of characters in the query sequence
// numSequences -- number of matching sequences potentially in the model */
static void posAllocateMemory(psim_typ posSearch, 
		       Int4 alphabetSize, Int4 querySize, Int4 numSequences)
{
  Int4	i,j;

// std::cerr << "==================== ENTERING posAllocateMemory( )...====================\n\n"; 
  NEW(posSearch->posCount,querySize+1,Int4);
  NEW(posSearch->posWtCount,querySize+1,double);
  NEWP(posSearch->posC,querySize+1,Int4);
  NEWP(posSearch->posMatrix,querySize+1,Int4);
  NEWP(posSearch->posPrivateMatrix,querySize+1,Int4);
  NEWP(posSearch->posMatchWeights,querySize+1,double);
  NEW(posSearch->posSigma,querySize+1,double);
  NEW(posSearch->posIntervalSizes,querySize+1,Int4);
  for(i = 0; i <= querySize; i++) {
  	NEW(posSearch->posC[i],alphabetSize+2,Int4);
  	NEW(posSearch->posMatrix[i],alphabetSize+2,Int4);
  	NEW(posSearch->posPrivateMatrix[i],alphabetSize+2,Int4);
  	NEW(posSearch->posMatchWeights[i],alphabetSize+2,double);
  }
  posSearch->pde_typMatrixLength = numSequences;
  NEWP(posSearch->pde_typMatrix,numSequences+1,pde_typ);
  for (i = 0; i <= numSequences; i++) {
    NEW(posSearch->pde_typMatrix[i],querySize+1,pde_typ);
    for(j = 0; j < querySize; j++) {
      posSearch->pde_typMatrix[i][j].letter = UNUSED;
      posSearch->pde_typMatrix[i][j].used = FALSE;
      posSearch->pde_typMatrix[i][j].e_value = 1.0;
      posSearch->pde_typMatrix[i][j].leftExtent = -1;
      posSearch->pde_typMatrix[i][j].rightExtent = querySize;
    }
  }  
  NEW(posSearch->posExtents,querySize+1,pde_typ);
  for(j = 0; j < querySize; j++) {
    posSearch->posExtents[j].used = FALSE;
    posSearch->posExtents[j].leftExtent = -1;
    posSearch->posExtents[j].rightExtent = querySize;
  }
  NEW(posSearch->posA,numSequences+ 1,double);
  NEW(posSearch->posRowSigma,numSequences+ 1,double);
}

void posCheckpointFreeMemory(psim_typ posSearch, Int4 querySize)
// Deallocate memory allocated in posReadCheckpoint
// posSearch -- pointer to record used in building the position-specific model
// querySize -- number of characters in the query sequence
{
  for(Int4 i = 0; i <= querySize; i++){
    free(posSearch->posFreqs[i]); free(posSearch->posMatrix[i]); 
    free(posSearch->posPrivateMatrix[i]);
  } free(posSearch->posFreqs); 
  free(posSearch->posMatrix); free(posSearch->posPrivateMatrix);
  posSearch->posMatrix=NULL;
}

void posCleanup(psim_typ posSearch, csi_typ * compactSearch)
// Cleanup position-specific  data structures after one pass
// Deallocate memory allocated in posAllocateMemory
// posSearch -- pointer to record used in building the position-specific model
// querySize -- number of characters in the query sequence
{
  Int4	i,querySize=compactSearch->qlength;
  
  if(posSearch->posCount) free(posSearch->posCount); 
  if(posSearch->posWtCount) free(posSearch->posWtCount); 
  if(posSearch->posExtents) free(posSearch->posExtents);
  if(posSearch->posSigma) free(posSearch->posSigma); 
  for(i = 0; i <= querySize; i++){
    if(posSearch->posC) free(posSearch->posC[i]); 
    free(posSearch->posMatrix[i]);
    free(posSearch->posPrivateMatrix[i]);
    if(posSearch->posMatchWeights) free(posSearch->posMatchWeights[i]);
    free(posSearch->posFreqs[i]);
  } free(posSearch->posFreqs); posSearch->posFreqs=0;
  if(posSearch->posC) free(posSearch->posC);
  if(posSearch->pde_typMatrix){
    for(i = 0; i <= posSearch->pde_typMatrixLength; i++)
    	free(posSearch->pde_typMatrix[i]);
    free(posSearch->pde_typMatrix);
  }
  free(posSearch->posMatrix); posSearch->posMatrix=NULL;
  free(posSearch->posPrivateMatrix);
  if(posSearch->posMatchWeights) free(posSearch->posMatchWeights);
  if(posSearch->posA) free(posSearch->posA);
  if(posSearch->posRowSigma) free(posSearch->posRowSigma);
  if(posSearch->posIntervalSizes) free(posSearch->posIntervalSizes);
  if(posSearch->posUseSequences) free(posSearch->posUseSequences);
}

static double minEvalueForSequence(sap_typ curGSeqAlign, sap_typ listOfGSeqAligns) 
// Find the lowest e-value among all seqAligns for the sequence represented by curGSeqAlign
{
   sap_typ	sap; 		// Index into listOfSeqALigns
   dsp_typ	curSegs, dsp; 	// Used to extract ids from curGSeqAlign, sap 
   double returnValue; 		// stores current best e-value
   Boolean	seen;   	// have we seen a seqAlign matching the sequence yet

   returnValue = curGSeqAlign->evalue;
   curSegs = (dsp_typ) curGSeqAlign->segs;
   UInt4 curId = curSegs->subject_id;  // AFN...
   seen = FALSE;
   for(sap=listOfGSeqAligns; sap; sap=sap->next) {
     dsp = (dsp_typ) sap->segs;
     if(curId == dsp->subject_id){	// AFN
         seen = TRUE;
	 if(sap->evalue < returnValue) returnValue = sap->evalue;
     } else if(seen)  break; // moved on to another sequence; stop looking.
   } return(returnValue);
}

// Count the number of seqAligns in a list (returned) and count the number 
// of distinct target sequences represented (passed back in numSequences);
// if useThreshold is TRUE, only those sequences with e-values below the 
// threshold are counted.
// IMPORTANT ASSUMPTION: All GSeqAligns with the same target sequence
// are consecutive in the list.
static Int4 countGSeqAligns(sap_typ listOfGSeqAligns, Int4 * numSequences, Boolean useThreshold, double threshold)
{
   sap_typ	curGSeqAlign, prevGSeqAlign;
   Int4		seqAlignCounter;
   dsp_typ	curSegs;
   UInt4 curId,prevId=0;  // AFN: Ids of target sequences in current & previous GSeqAlign

   seqAlignCounter = 0;
   *numSequences = 0;
   curGSeqAlign = listOfGSeqAligns; prevGSeqAlign = NULL;
   while (NULL != curGSeqAlign) {
     seqAlignCounter++;
     curSegs = (dsp_typ) curGSeqAlign->segs;
     curId = curSegs->subject_id; 		// AFN
     if((NULL == prevGSeqAlign) || (curId != prevId))	// AFN
       if (!useThreshold || 
		(threshold > minEvalueForSequence(curGSeqAlign,listOfGSeqAligns)))
	 	(*numSequences)++;
     prevGSeqAlign = curGSeqAlign;
     prevId = curId;
     curGSeqAlign = curGSeqAlign->next;
   }
   return(seqAlignCounter);
}

psim_typ MakePSIM(Int4 size_dbs, a_type AB)
{
	psim_typ posSearch;
	NEW(posSearch,1,psim_type);
	posSearch->AB = AB;
	posSearch->HitSet= MakeSet(size_dbs+2);
	ClearSet(posSearch->HitSet);	// initialize to no hits...
	return posSearch;
}

void	NilPSIM(psim_typ posSearch)
{
	if(posSearch != NULL){
	   NilSet(posSearch->HitSet);
           // if(posSearch->posMatrix) posCleanup(posSearch,compactSearch);
           posFreeInformation(posSearch); // AFN: solves purify memory problem.
           free(posSearch); 
        }
}

void	posConvergencedTest(psim_typ posSearch, bsb_typ search, sap_typ sap)
//  Determines if the search has converged from round to round.
//  Checks whether every new sequence found is in posSearch->HitSet.
//  posSearch is the data structure representing the parameters of the position-specific part
//  search represents the overall  GBLAST search
//  listOfGSeqAligns is one representation of the results of the current round.
//  If thissPassNum is 1, then it checks only to see that some sequence
//  distinct from the query was found.
{

  search->posConverged = TRUE; 
  if(sap == NULL) return;
  if(sap && sap->next==NULL){  // only one hit...
	dsp_typ dsp = sap->segs;
	if(IdentSeqs(dsp->sE,dsp->qE)) return; // query only found itself (converged)
  } 
  for( ; sap; sap=sap->next) {
        if(sap->evalue < search->ethresh) {
	  UInt4 thisId = sap->segs->subject_id; 
          if(!MemberSet(thisId,posSearch->HitSet)) { // Speed up convergence test...
              AddSet(thisId,posSearch->HitSet); search->posConverged = FALSE;
          }
        }
  }
}

static void posCancel(psim_typ posSearch, csi_typ *csi, Int4 first, Int4 second,
	Int4 matchStart, Int4 intervalLength)
{
// Eliminate the matches from sequence second starting at position
// matchStart and extending for intervalLength characters.
  Int4 c, i;
  Boolean stillNeeded;

  for(c = matchStart, i = 0; i < intervalLength; i++, c++) {
    posSearch->pde_typMatrix[second][c].used = FALSE;
    // posSearch->pde_typMatrix[second][c].letter = 0;		// Original == gap
    posSearch->pde_typMatrix[second][c].letter = GAP_CHAR; // AFN old mapping == 21 
    // Switched these values after adding char '-' == 21 to alphabet.
  } stillNeeded = FALSE;
  for(c = 0; c < csi->qlength; c++)
    if(posSearch->pde_typMatrix[second][c].used) { stillNeeded = TRUE; break; }
  if(!stillNeeded) posSearch->posUseSequences[second] = FALSE;
}

static void  posPurgeMatches(psim_typ posSearch, csi_typ * compactSearch)
// Eliminate sequences that are identical to the query and partial alignments
// that are identical in two matching sequences
{
  Int4 i, j; // index over sequences
  Boolean matchesQuery; /*Is a matching sequence identical to the query?*/
  Int4 c; /*index over demographics of matching sequence*/
  Int4 state; // state of checking for a match
  Int4 intervalLength, matchStart; /*Length and start of a matching region*/
  Int4 matchNumber; /*number of characters matching*/

  MEW(posSearch->posUseSequences,posSearch->posNumSequences + 1,Boolean);
  for(i=0; i <= posSearch->posNumSequences; i++) posSearch->posUseSequences[i]=TRUE;
  for(i = 1; i <= posSearch->posNumSequences; i++){
    matchesQuery = TRUE;
    for(c=0; c < compactSearch->qlength; c++) {
      if ((!posSearch->pde_typMatrix[i][c].used) ||
          (posSearch->pde_typMatrix[i][c].letter != posSearch->pde_typMatrix[0][c].letter)) {
        matchesQuery = FALSE;
        break;
      }
    }
    if(matchesQuery) posSearch->posUseSequences[i] = FALSE;
  }
  pde_typ       *pde_typMtrx0 = posSearch->pde_typMatrix[0];
  for(j = 1; j <= posSearch->posNumSequences; j++) {
    if (!posSearch->posUseSequences[j]) continue;
    state = POS_COUNTING;
    matchStart = 0; intervalLength = 0; matchNumber = 0;
    pde_typ       *pde_typMtrx = posSearch->pde_typMatrix[j];
    for(c=0; c < compactSearch->qlength; c++) {
      if(pde_typMtrx[c].used) {
	if((pde_typMtrx0[c].letter != XChar) && (pde_typMtrx[c].letter != XChar)) { 
	  if(state == POS_RESTING) {
	    matchStart = c; intervalLength = 1;
	    state = POS_COUNTING; matchNumber = 0;
	  } else intervalLength++;
	  if (pde_typMtrx[c].used &&
		(pde_typMtrx0[c].letter == pde_typMtrx[c].letter)) matchNumber++;
	}
      } else {
	if (state == POS_COUNTING) {
	  if ((intervalLength > 0) && (matchNumber == intervalLength))
	    posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength);
	  state = POS_RESTING;
	}
      }
    }
    if (state == POS_COUNTING) /*at end of sequence i*/
      if ((intervalLength > 0) && (matchNumber == intervalLength))
	posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength);
  }
  for (i = 1; i < posSearch->posNumSequences; i++) {
    if (!posSearch->posUseSequences[i]) continue;
    pde_typ *pde_typMtrxI = posSearch->pde_typMatrix[i];
    for(j = i+1; j <= posSearch->posNumSequences; j++) {
      if (!posSearch->posUseSequences[j]) continue;
      state = POS_COUNTING;
      matchStart = 0; intervalLength = 0; matchNumber = 0;
      pde_typ *pde_typMtrxJ = posSearch->pde_typMatrix[j];
      for(c=0; c < compactSearch->qlength; c++){
	pde_typ *pmIC = &pde_typMtrxI[c];
	pde_typ *pmJC = &pde_typMtrxJ[c];
	if(pmIC->used || pmJC->used) {
	  if((pmIC->letter != XChar) && (pmJC->letter != XChar)) { 
	    if(state == POS_RESTING) {
	      matchStart = c; intervalLength = 1;
	      state = POS_COUNTING; matchNumber = 0;
	    } else intervalLength++;
	    if(pmIC->used && pmJC->used && (pmIC->letter == pmJC->letter)) matchNumber++;
	    }
	} else {
	  if(state==POS_COUNTING) {
	    if((intervalLength > 0) && ((((double ) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	      posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength);
	    state = POS_RESTING;
	  }
	}
      }
      if (state == POS_COUNTING) /*at end of sequence i*/
	if ((intervalLength > 0) && ((((double ) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	  posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength);
    }
  }
}

static void posDemographics(psim_typ posSearch, csi_typ * compactSearch, 
	sap_typ listOfGSeqAligns)
// Compute general information about the sequences that matched on the i-th 
// pass such as how many matched at each query position and what letter matched.
{
   Uint1Ptr q; /*pointers into query */
   Uint1Ptr s; /*pointer into a matching string */
   unsigned char q_use,s_use; 	// AFN fixes
   Int4 length, subjectLength;  /*length of query and subject*/
   Int4 c; /*index into a string*/
   Int4 numseq, numGSeqAligns;  /*number of matching sequences and GSeqAligns*/
   Int4 seqIndex;  /*index for the array of matching sequences*/
   Int4 matchLength; /*length of a match*/
   Int4  queryOffset, subjectOffset, retrievalOffset;  /*offsets needed to make a match align*/
   Int4 qplace, splace; /*index into query string and matching string*/
   sap_typ curGSeqAlign, prevGSeqAlign; /*pointers into listOfGSeqAligns*/
   dsp_typ curSegs, prevSegs;  /*used to extract alignments from curGSeqAlign*/
   Int4 startQ, startS; /*Indices into array of starting positions*/
   Int4 numsegs; /*Number of pieces in the gapped alignment*/
   Int4 segIndex; /*Index for which piece we are at*/
   double thisEvalue;  /*evalue of current partial alignment*/

   q = compactSearch->query;
   length = compactSearch->qlength;
   for(c = 0; c < length; c++) {
     if(q[c] == Xchar) q_use=XChar; else q_use=q[c];	// AFN fix for mapping differences
     posSearch->pde_typMatrix[0][c].letter = (signed char) q_use;
     posSearch->pde_typMatrix[0][c].used = TRUE;
     posSearch->pde_typMatrix[0][c].leftExtent = 0;
     posSearch->pde_typMatrix[0][c].rightExtent = length;
     posSearch->pde_typMatrix[0][c].e_value = compactSearch->ethresh/2;
     posSearch->posC[c][q[c]]++;
     posSearch->posCount[c]++;
   }
   numGSeqAligns = countGSeqAligns(listOfGSeqAligns, &numseq, TRUE, compactSearch->ethresh);
   posSearch->posNumSequences = numseq;
   /*use only those sequences below e-value threshold*/
   seqIndex = 0;
   curGSeqAlign = listOfGSeqAligns;
   prevGSeqAlign = NULL;
   while (curGSeqAlign != NULL) {
     if((thisEvalue = curGSeqAlign->evalue) < compactSearch->ethresh) {
       curSegs = (dsp_typ) curGSeqAlign->segs;
       if (NULL != prevGSeqAlign) {
	 prevSegs = (dsp_typ) prevGSeqAlign->segs;
	 if(curSegs->subject_id != prevSegs->subject_id) seqIndex++;
       }
       s = GetSequenceWithGDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
       startQ = 0; startS = 1;
       numsegs = curSegs->numseg;
       for(segIndex = 0; segIndex < numsegs; segIndex++) {
	 queryOffset = curSegs->starts[startQ];
         if (curSegs->starts[startS] != GAP_HERE)
	   subjectOffset = curSegs->starts[startS] - retrievalOffset;
	 else subjectOffset = GAP_HERE;
	 matchLength = curSegs->lens[segIndex];
	 if ((GAP_HERE ) == queryOffset) { ; /*do nothing, gap in query*/ }
	 else if ((GAP_HERE) == subjectOffset) {
	     for(c = 0, qplace = queryOffset; c < matchLength; c++, qplace++) {
	       posSearch->pde_typMatrix[seqIndex + 1][qplace].used = TRUE;
	       posSearch->pde_typMatrix[seqIndex + 1][qplace].letter = GAP_CHAR;
	       posSearch->pde_typMatrix[seqIndex + 1][qplace].e_value = 1.0;
	     }
	 } else {  /*no gap*/
	     for(c=0, qplace=queryOffset, splace=subjectOffset; c < matchLength;
			c++, qplace++, splace++) {
	       if ((!posSearch->pde_typMatrix[seqIndex+1][qplace].used) ||
		   (thisEvalue < posSearch->pde_typMatrix[seqIndex+1][qplace].e_value)) 
		 if (!posSearch->pde_typMatrix[seqIndex+1][qplace].used) {
		     if(s[splace]==Xchar) s_use= XChar; else s_use=s[splace];	// AFN fix
		     posSearch->pde_typMatrix[seqIndex+1][qplace].letter 
				= (signed char) s_use; 
				// = (signed char) s[splace];  // OLD (AFN)
		     posSearch->pde_typMatrix[seqIndex+1][qplace].used = TRUE; 
		     posSearch->pde_typMatrix[seqIndex+1][qplace].e_value = 
		       thisEvalue;
		 }
	     }
	 } startQ += 2; startS += 2;
       } prevGSeqAlign = curGSeqAlign;
     } curGSeqAlign = curGSeqAlign->next;
   } // closes the while loop over seqAligns
 }

static void posComputeExtents(psim_typ posSearch, csi_typ * compactSearch)
{
   Int4 seqIndex; /*index of sequence*/
   Int4 length; /*length of query*/
   Int4 qplace, qplace2; /*place in query*/
   Int4 numseq; /*number of sequences including query*/
   Uint1Ptr q; /*pointers into query */

   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   q = compactSearch->query;
   for(seqIndex = 0; seqIndex < numseq; seqIndex++) {
     if (!posSearch->posUseSequences[seqIndex+1])
       continue;
     if ((posSearch->pde_typMatrix[seqIndex+1][0].used)
	 && (posSearch->pde_typMatrix[seqIndex+1][length-1].letter != GAP_CHAR))
       posSearch->pde_typMatrix[seqIndex+1][0].leftExtent = 0;
     for(qplace = 1; qplace < length; qplace++)
       if(posSearch->pde_typMatrix[seqIndex+1][qplace].used)
	 if(posSearch->pde_typMatrix[seqIndex+1][qplace-1].used)
	   posSearch->pde_typMatrix[seqIndex+1][qplace].leftExtent =
	     posSearch->pde_typMatrix[seqIndex+1][qplace -1].leftExtent;
	 else
	   posSearch->pde_typMatrix[seqIndex+1][qplace].leftExtent = qplace;
     if ((posSearch->pde_typMatrix[seqIndex+1][length-1].used)
	 && (posSearch->pde_typMatrix[seqIndex+1][length-1].letter != GAP_CHAR))
       posSearch->pde_typMatrix[seqIndex+1][length-1].rightExtent = length -1;
     for(qplace = length -2; qplace >= 0; qplace--)
       if(posSearch->pde_typMatrix[seqIndex+1][qplace].used)
	 if(posSearch->pde_typMatrix[seqIndex+1][qplace+1].used)
	   posSearch->pde_typMatrix[seqIndex+1][qplace].rightExtent =
	     posSearch->pde_typMatrix[seqIndex+1][qplace + 1].rightExtent;
	 else
	   posSearch->pde_typMatrix[seqIndex+1][qplace].rightExtent = qplace;
     for(qplace = 0; qplace < length; qplace++) 
       if (posSearch->pde_typMatrix[seqIndex+1][qplace].used) {
	 posSearch->posExtents[qplace].leftExtent = 
		MAX(posSearch->posExtents[qplace].leftExtent,
		  posSearch->pde_typMatrix[seqIndex+1][qplace].leftExtent);
	 posSearch->posExtents[qplace].rightExtent = 
		MIN(posSearch->posExtents[qplace].rightExtent,
		  posSearch->pde_typMatrix[seqIndex+1][qplace].rightExtent);
	 
       }

     for(qplace = 0; qplace < length; qplace++) // used to check qplace for GAP_CHAR here
       if (posSearch->pde_typMatrix[seqIndex+1][qplace].used) {
	 posSearch->posC[qplace][posSearch->pde_typMatrix[seqIndex+1][qplace].letter]++;
	 posSearch->posCount[qplace]++; /*Add to number of matches in this query position*/
       }
   }
   for(qplace = 0; qplace < length; qplace++)
     posSearch->posIntervalSizes[qplace] = posSearch->posExtents[qplace].rightExtent - 
       posSearch->posExtents[qplace].leftExtent + 1;
   for(qplace =0; qplace < length; qplace++) {
     if(Xchar == q[qplace]) {
       posSearch->posIntervalSizes[qplace] = 0;
       for(qplace2 = 0; qplace2 <qplace; qplace2++) {
	 if(posSearch->posExtents[qplace2].rightExtent >= qplace)
	   posSearch->posIntervalSizes[qplace2]--;
       }
     }
   }
}
 
static void posComputeSequenceWeights(psim_typ posSearch, csi_typ *compactSearch)
// Compute Henikoff weight of each sequence and letter in each position
// WARNING (AFN): THIS IS NOT HANDLING X's IN SEQUENCES PROPERLY.
//   THIS NEEDS TO BE FIXED!!!!!  (DUE TO MY ALPHABET REPRESENTATION.)
{
   Int4 length; /*length of query*/
   Int4 numseq, seqIndex; /*number of matches, index for them*/
   Int4  i; /*index over a multi-alignment block*/
   Int4 qplace; /*index into query*/
   double Sigma; /*Number of different characters occurring in matches within
                   a multi-alignment block, excluding identical columns*/
   double intervalSigma; /*Same as Sigma but includes identical columns*/
   Int4 alphabetSize; /*number of characters in alphabet*/
   Int4 *participatingSequences; /*array of participating sequences at a position*/
   Int4 *oldParticipatingSequences; /*array of participating sequences at a position*/
   Int4 posLocalVariety;  /*number of different characters at a position*/
   Int4 *posLocalC; /*counts of how many of each letter in this column*/
   Int4 c;
   Int4 thisSeq;
   Int4 numParticipating; /*number of sequences in this alignment block*/
   Int4 oldNumParticipating; /*number of sequences in this alignment block*/
   Boolean newSequenceSet; 
   Int4 p; /*index on sequences*/

   alphabetSize = compactSearch->alphabetSize;
   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   NEW(participatingSequences, numseq+1,Int4);
   NEW(oldParticipatingSequences,numseq+1,Int4);
   NEW(posLocalC,alphabetSize +2, Int4);
   for(qplace=0; qplace < length; qplace++) posSearch->posSigma[qplace] = 0.0;
   numParticipating = 0;
   for(qplace = 0; qplace < length; qplace++) {
     if((posSearch->posCount[qplace] > 1) 
	&& (posSearch->posIntervalSizes[qplace] > 0)) {
       oldNumParticipating = numParticipating;
       for(p =0; p < numParticipating; p++)
           oldParticipatingSequences[p] = participatingSequences[p];
       numParticipating = 0;
       for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
         if (!posSearch->posUseSequences[seqIndex]) continue;
	 if ((posSearch->pde_typMatrix[seqIndex][qplace].used) &&
	     (posSearch->pde_typMatrix[seqIndex][qplace].letter != GAP_CHAR)) {
	   participatingSequences[numParticipating] = seqIndex;
	   numParticipating++;
	 }
       }
       newSequenceSet = TRUE;
       if (numParticipating == oldNumParticipating) {
         for(p = 0; p < numParticipating; p++)
           if (oldParticipatingSequences[p] != participatingSequences[p])
             break;
         if (p == numParticipating) newSequenceSet = FALSE;
       }
       if(newSequenceSet) {
	 Sigma = 0; intervalSigma = 0;
	 for(seqIndex = 0; seqIndex <= numseq; seqIndex++) {
	   if (!posSearch->posUseSequences[seqIndex]) continue;
	   posSearch->posRowSigma[seqIndex] = 0.0;
	   posSearch->posA[seqIndex] = 0.0;
	 }
	 for(i = posSearch->posExtents[qplace].leftExtent;
	      i <= posSearch->posExtents[qplace].rightExtent; i++) {
	   posLocalVariety = 0;
	   for(c = 0; c <= alphabetSize; c++) posLocalC[c] = 0;
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     // used to check for GAP here
	     if(0 == posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter])
	       // letter (not a gap) not seen before in this query pos.
	       posLocalVariety++;  
	     posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]++;
	   }
	   intervalSigma += posLocalVariety;
	   if (posLocalVariety > 1) {
	     Sigma += posLocalVariety;
	     for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	       thisSeq = participatingSequences[seqIndex];
	       // used to check for gap here
	       posSearch->posRowSigma[thisSeq] += 
		 ( 1.0 / (double ) posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]);
	     }
	   }
	 }
       }
       // NOTE that 1.0/(# seqs with s[i] @ position i) summed over all seqs ==
       //    		== # of residue types at position i.
       // e.g. AAABB = 1/3 + 1/3 + 1/3 + 1/2 + 1/2 = 2 residue types at i.
       if (Sigma > 0) {
	 for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = posSearch->posRowSigma[thisSeq]/Sigma;
	 }
       } else {
         for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = ((double ) 1 / (double ) numParticipating);
         }
       }
       posSearch->posSigma[qplace] = intervalSigma;
#if 0	// Original code.
       for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 posSearch->posMatchWeights[qplace][posSearch->pde_typMatrix[thisSeq][qplace].letter] 
		+= posSearch->posA[thisSeq];
       }
#endif
#if 1	// AFN: Compute relative weights and effective number of seqs...
       double maxWt=0.0;
       for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 posSearch->posMatchWeights[qplace][posSearch->pde_typMatrix[thisSeq][qplace].letter] 
		+= posSearch->posA[thisSeq];
	 if(maxWt < posSearch->posA[thisSeq]) maxWt = posSearch->posA[thisSeq];
       }
       double	effective_num=0;
       for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 effective_num += posSearch->posA[thisSeq]*(1.0/maxWt); // set max wt seq to 1.0.
       }
       posSearch->posWtCount[qplace]=effective_num;
#if 0	// AFN DEBUG:
	if(qplace == 96){
	   fprintf(stderr,"96: leftExtent = %d; rightExtent = %d; posIntervalSizes = %d; Sigma = %g\n",
		posSearch->posExtents[qplace].leftExtent,
		posSearch->posExtents[qplace].rightExtent,
		posSearch->posIntervalSizes[qplace],Sigma);
		// assert(qplace != 96);
	}
#endif
       // fprintf(stderr,"effective_num_seq=%.2f;",effective_num);
       // fprintf(stderr,"actual number(%d)=%d\n",qplace+1,posSearch->posCount[qplace]);
#endif
#if 0	// DEBUG code AFN.
	double sumX=0;
	for(c = 0; c < alphabetSize; c++){
	  sumX+=posSearch->posMatchWeights[qplace][c];
	}
	if(sumX > 1.01){
	   if(Sigma > 0){
	 	fprintf(stderr,"\nleftExtent =%d; rightExtent = %d\n",
			posSearch->posExtents[qplace].leftExtent,
			posSearch->posExtents[qplace].rightExtent);
		fprintf(stderr,"%d: Sigma = %f; total = %f\n",qplace,Sigma,sumX);
		sumX=0;
		for(seqIndex = 0; seqIndex < numParticipating; seqIndex++){
			thisSeq = participatingSequences[seqIndex];
			fprintf(stderr,"seq %4d: ", thisSeq);
			sumX+=posSearch->posRowSigma[thisSeq];
		       for(i=posSearch->posExtents[qplace].leftExtent; 
				i <= posSearch->posExtents[qplace].rightExtent; i++){
				fprintf(stderr,"%c",
				  AlphaChar(posSearch->pde_typMatrix[thisSeq][i].letter,
				    posSearch->AB));	
			} fprintf(stderr," rowsigma =%f\n",
				posSearch->posRowSigma[thisSeq]);
		}
		fprintf(stderr,"Sum posRowSigma = %f\n",sumX);
#if 1	// AFN: recompute...
	 fprintf(stderr,"pos  ",i);
	 for(c = 0; c <= alphabetSize; c++){
		fprintf(stderr,"  %c", AlphaChar(c,posSearch->AB));
	 }  fprintf(stderr,"\n");
	 Sigma = 0; sumX=0;
	 for(i = posSearch->posExtents[qplace].leftExtent;
	      i <= posSearch->posExtents[qplace].rightExtent; i++) {
	   posLocalVariety = 0;
	   for(c = 0; c <= alphabetSize; c++) posLocalC[c] = 0;
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     if(0 == posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter])
	       posLocalVariety++;  
	     posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]++;
	   }
	   if(posLocalVariety > 1) {
	     Sigma += posLocalVariety;
	     for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	       thisSeq = participatingSequences[seqIndex];
		if(posSearch->pde_typMatrix[thisSeq][i].letter != GAP_CHAR)
		sumX += ( 1.0 / (double ) posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]);
	     }
	   }
	   fprintf(stderr,"%3d: ",i);
	   for(c = 0; c <= alphabetSize; c++){
		fprintf(stderr,"%3d",posLocalC[c]);
	   } fprintf(stderr," (%d)\n",posLocalVariety);
	 }
	 fprintf(stderr,"Sigma = %f; Sum posRowSigma = %f\n",Sigma,sumX);
#endif
	   } assert(sumX <= 1.01);
	}
#endif
     }
   }
   free(participatingSequences);
   free(oldParticipatingSequences);
   free(posLocalC);
}

static double countsFunction(double Sigma, Int4 intervalLength)
{ return(Sigma / intervalLength - 1); }

static double posit_rounddown(double value){ return (double) GNlm_Nint(value); }

static BooLean posCheckWeights(psim_typ posSearch, csi_typ * compactSearch)
// check that weights add to 1 in each column
{
   Int4		r,length,alphabetSize; // length of query & # characters in alphabet
   double	runningSum; 	     // partial total for a column

   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;
   unsigned char *q=compactSearch->query;	// pointer to query
   for(Int4 i = 0; i < length; i++) {
     if ((posSearch->posCount[i] > 1) && (q[i] != Xchar)) {
       runningSum = 0;
       for(Int4 a=0; a < alphabetSize; a++) 
           runningSum += posSearch->posMatchWeights[i][a];
       if((runningSum < 0.99) || (runningSum > 1.01)) {
   	  a_type	AB=posSearch->AB;
	  fprintf(stderr,"alphabetSize = %d; runningSum = %g; column = %d\n",
			alphabetSize,runningSum,i);
	  fprintf(stderr,"posMatchWeights matrix\npos  ");
          for(r=0; r<= NAlpha(AB); r++) fprintf(stderr,"  %c",AlphaChar(r,AB));
	  fprintf(stderr,"\n");
   	  for(Int4 j = 0; j < length; j++) {
		fprintf(stderr,"%3d: ",j);
		runningSum=0.0;
		for(r=0; r <= NAlpha(AB); r++){
		    runningSum+=posSearch->posMatchWeights[j][r];
		    fprintf(stderr,"%3d",(Int4)(100*posSearch->posMatchWeights[j][r]));
	        } fprintf(stderr," %.3f\n",runningSum);
	  } fprintf(stderr,"\n");
	  outputPosMatrix(stderr,posSearch,compactSearch,AB);
	  return FALSE;
       }
     }
   } return TRUE;
}

static void  posFreqsToInformation(psim_typ posSearch, csi_typ *compactSearch)
// Fill in information content per position from pseudo-count frequencies
{
   double	qOverPEstimate;   // intermediate term
   double	infoSum; 	  // information content sum for this position
  
   Int4 length = compactSearch->qlength; // length of the query
   Int4 alphabetSize = compactSearch->alphabetSize;
   for(Int4 c = 0; c < length; c++) {
     infoSum = 0;
     for(Int4 a=0; a < alphabetSize; a++) {
       if(compactSearch->standardProb[a] > posEpsilon) {
         qOverPEstimate=posSearch->posFreqs[c][a]/compactSearch->standardProb[a];
         if (qOverPEstimate > posEpsilon)
	   infoSum+=posSearch->posFreqs[c][a]*log(qOverPEstimate)/NCBIMATH_LN2;
       }
     } posSearch->posInformation[c]=infoSum;
   }
}

void posFreqsToMatrix(psim_typ posSearch, csi_typ *compactSearch)
// Convert pseudo-count frequencies to a score matrix.
{
   Int4		c,a,alphabetSize;
   double	qOverPEstimate, value; 	// intermediate terms
   Boolean	allZeros; 		// are all frequencies in a column 0?

   unsigned char *q = compactSearch->query; 	// pointer to the query
   Int4 length = compactSearch->qlength; 	//length of the query
   alphabetSize = compactSearch->alphabetSize;
   double lambda = compactSearch->lambda_ideal; // Karlin-Altschul parameter
   for(c = 0; c < length; c++) {
     allZeros = TRUE;
     for(a = 0; a < alphabetSize; a++) {
       // Division compensates for multiplication in posComputePsedoFreqs
       if(compactSearch->standardProb[a] > posEpsilon)
	  qOverPEstimate = posSearch->posFreqs[c][a]/compactSearch->standardProb[a];
       else qOverPEstimate = 0.0;
       if(qOverPEstimate != 0.0) allZeros = FALSE;
       if(0.0 == qOverPEstimate || (compactSearch->standardProb[a] < posEpsilon))
	 posSearch->posPrivateMatrix[c][a] = GBLAST_SCORE_MIN;
       else {
	 value = log(qOverPEstimate)/lambda;
	 posSearch->posPrivateMatrix[c][a] = (Int4) posit_rounddown(POSIT_SCALE_FACTOR*value);
       }
     }    
     if(allZeros) {
       for(a = 0; a < alphabetSize; a++) {
         posSearch->posMatrix[c][a] = compactSearch->matrix[q[c]][a];
	 if(compactSearch->matrix[q[c]][a] == GBLAST_SCORE_MIN)
		posSearch->posPrivateMatrix[c][a] = GBLAST_SCORE_MIN;
	 else posSearch->posPrivateMatrix[c][a] = 
		POSIT_SCALE_FACTOR*compactSearch->matrix[q[c]][a];
       }
     }
   }
   for(a = 0; a < alphabetSize; a++) {
     posSearch->posPrivateMatrix[length][a] = GBLAST_SCORE_MIN;
   }
}

#include "blosum62.h"
static double **posComputePseudoFreqs(psim_typ posSearch, csi_typ * compactSearch,
	Boolean Cpos)
{
   Uint1Ptr q;  		// pointer to the query
   Int4 length;  		// length of the query
   Int4 c; 
   Int4 a, aSub, alphabetSize; 	// loop indices and size of alphabet
   double lambda; 		// Karlin-Altschul parameter
   double Sigma;  		// number of characters in an interval
   Int4 intervalLength;  	// length of a block
   double pseudo, numerator, denominator, qOverPEstimate; // intermediate terms
   double infoSum; 		// sum used for information content
   double **posFreqs; 		// store frequencies*/

   q = compactSearch->query;
   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;
   lambda = compactSearch->lambda_ideal;
   NEWP(posFreqs,length+1,double);
   for(c=0; c <= length; c++) NEW(posFreqs[c],alphabetSize+2,double);
   for(c = 0; c < length; c++) {
     if((posSearch->posCount[c] > 1) && (Xchar != q[c])) {
       infoSum = 0;
       for(a = 0; a < alphabetSize; a++) {
         if(compactSearch->standardProb[a] > posEpsilon) {
	   pseudo = 0;
	   for (aSub = 0; aSub < alphabetSize; aSub++){
#if 1
	     if(compactSearch->matrix[a][aSub] != GBLAST_SCORE_MIN){
	       pseudo += (posSearch->posMatchWeights[c][aSub] *
			exp(lambda * compactSearch->matrix[a][aSub]));
	     }
#else	// NEW: this is nearly equivalent to the old method (using freq ratios...)
	     if(compactSearch->matrix[a][aSub] != GBLAST_SCORE_MIN){
	       pseudo += (posSearch->posMatchWeights[c][aSub] * blosum62L[a][aSub]);
	     }
#endif
	   }
	   pseudo *= (compactSearch->pseudoCountConst);
           Sigma = posSearch->posSigma[c];
           intervalLength = posSearch->posIntervalSizes[c];
	   numerator = pseudo + 
             (countsFunction(Sigma, intervalLength) * posSearch->posMatchWeights[c][a]/
                compactSearch->standardProb[a]);
	   denominator = 
		countsFunction(Sigma, intervalLength) + (compactSearch->pseudoCountConst); 
	   qOverPEstimate = numerator / denominator;
	   // Note artificial multiplication by standard probability to normalize
           posFreqs[c][a] = qOverPEstimate * compactSearch->standardProb[a];
	 if (0.0 != qOverPEstimate && (compactSearch->standardProb[a] > posEpsilon))
	   infoSum += qOverPEstimate * 
		compactSearch->standardProb[a] * log(qOverPEstimate)/ NCBIMATH_LN2;
	 } else posFreqs[c][a] = 0.0;
       } if (Cpos) posSearch->posInformation[c] = infoSum;
     } else for(a = 0; a < alphabetSize; a++) { posFreqs[c][a] = 0; }
   }
  return(posFreqs);
}

static void posScaling(psim_typ posSearch, csi_typ * compactSearch)
{
	brm_typ *brm = new brm_typ(compactSearch->alphabetSize, 
		compactSearch->qlength, compactSearch->query,
		compactSearch->standardProb, posSearch->posMatrix,
		posSearch->posPrivateMatrix, compactSearch->kbp_std,
		compactSearch->kbp_psi, compactSearch->kbp_gap_std,
		compactSearch->kbp_gap_psi, compactSearch->lambda_ideal,
			compactSearch->K_ideal);
	brm->Scale( ); delete brm;
}

Int4Ptr * CposComputation(psim_typ posSearch, bsb_typ search,
	csi_typ *compactSearch, sap_typ listOfGSeqAligns, Char *ckptFileName)
{
	Int4Ptr *rtnMtx;
	
	FILE *chkfp=0; /*file in which to take the checkpoint*/
	if(ckptFileName) chkfp = open_file(ckptFileName,"","w");
	rtnMtx = CposComputation(posSearch, search,compactSearch, listOfGSeqAligns, chkfp);
	if(chkfp) fclose(chkfp);
}

Int4Ptr * CposComputation(psim_typ posSearch, bsb_typ search,
	csi_typ *compactSearch, sap_typ listOfGSeqAligns, FILE *chkfp)
{
  search->posConverged = FALSE;
  posAllocateMemory(posSearch, compactSearch->alphabetSize,
		compactSearch->qlength, search->hitlist_count);
  posDemographics(posSearch, compactSearch, listOfGSeqAligns);
  posPurgeMatches(posSearch, compactSearch);
  posComputeExtents(posSearch, compactSearch);
  posComputeSequenceWeights(posSearch, compactSearch);
  if(!posCheckWeights(posSearch, compactSearch)){
	  PutMultiGSeqAlign(stderr,listOfGSeqAligns,60,posSearch->AB);
	  assert(!"ERROR IN WEIGHTS"); // print_error("ERROR IN WEIGHTS");
  }
  posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, TRUE);
  if(NULL != chkfp) posTakeCheckpoint(posSearch,compactSearch,chkfp);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  FixposMatrix(posSearch, compactSearch); // AFN 11/7/01 fix...
  return posSearch->posMatrix;
}

#define POSIT_DEBUG

void outputPosMatrix(FILE *fp,psim_typ posSearch, csi_typ *compactSearch, a_type AB)
// Print out the position-specific matrix
{
   Uint1Ptr q; /*query sequence*/
   Int4 i, index; /*loop indices*/
   Int4	c; /*index over alphabet*/
   Int4 length; /*length of query*/

   q = compactSearch->query; length = compactSearch->qlength;
   index = 0;   
#if 1	// DEBUG
  if(posSearch->posFreqs){
   for(i = 0; i < length; i++) {
    if(posSearch->posFreqs[i]){
	fprintf(fp,"%d:  ",i+1);
  	for(c = 0; c < compactSearch->alphabetSize; c++) {
	   fprintf(fp,"%.4f ",posSearch->posFreqs[i][c]);
  	} fprintf(fp,"\n");
    }
   } fprintf(fp,"\n\n");
  }
#endif
// #ifdef POSIT_DEBUG
   fprintf(fp,"\nCharacter Frequencies by position\n");
   fprintf(fp,"         ");
   for(c = 0; c<= nAlpha(AB); c++) fprintf(fp,"  %c",AlphaChar(c,AB));
   for(i=0; i < length; i++) {
     fprintf(fp,"\n%5d %c   ", i + 1, AlphaChar(q[i],AB));
     for(c = 0; c<= nAlpha(AB); c++) fprintf(fp,"%2d ", posSearch->posC[i][c]);
   }
   fprintf(fp,"\n\n");
#if 1	// ORIGINAL MATRIX
   // fprintf(fp,"\nposition counts used. multiplied by 10 and rounded and");
   fprintf(fp,"\nposition counts used. multiplied by 1 and");
   fprintf(fp,"\nposition character weights used, multiplied by 10 and rounded\n");
   fprintf(fp,"        Counts  ");
   for(c = 0; c<= nAlpha(AB); c++) fprintf(fp,"  %c",AlphaChar(c,AB));
   fprintf(fp," Extent ");
   for(i=0; i < length; i++) {
     fprintf(fp,"\n%5d %c   ", i + 1, AlphaChar(q[i],AB));
     if ((posSearch->posCount[index+i] > 1) && (Xchar != q[index+i])){
       // fprintf(fp,"%4d ", (Int4) posit_rounddown(10 * countsFunction
		// (posSearch->posSigma[i],posSearch->posIntervalSizes[i])));
     // else fprintf(fp,"      ");
       fprintf(fp,"%6.4f ", countsFunction
		(posSearch->posSigma[i],posSearch->posIntervalSizes[i]));
     } else fprintf(fp,"-1.000 ");
     for(c = 0; c<= nAlpha(AB); c++) 
       if((posSearch->posMatrix[i][c] == GBLAST_SCORE_MIN) ||
             (0.0 == posSearch->posMatchWeights[i][c])) fprintf(fp," - ");
         else fprintf(fp,"%2d ",
	   (Int4) posit_rounddown(10 * posSearch->posMatchWeights[i][c]));
     	   fprintf(fp," %4d",
		posSearch->posExtents[i].rightExtent 
		- posSearch->posExtents[i].leftExtent +1);
   } fprintf(fp,"\n\n");
#endif
#if 0	// NEW FLOAT MATRIX
   fprintf(fp,"\nposition counts used. multiplied by 1 and");
   fprintf(fp,"\nposition character weights used, multiplied by 1\n");
   fprintf(fp,"        Counts");
   for(c = 0; c<= nAlpha(AB); c++) fprintf(fp,"  %c",AlphaChar(c,AB));
   fprintf(fp," Extent ");
   for(i=0; i < length; i++) {
     fprintf(fp,"\n%5d %c   ", i + 1, AlphaChar(q[i],AB));
     if ((posSearch->posCount[index+i] > 1) && (Xchar != q[index+i]))
       fprintf(fp,"%.4f ", (countsFunction
		(posSearch->posSigma[i],posSearch->posIntervalSizes[i])));
     else fprintf(fp,"     ");
     for(c = 0; c<= nAlpha(AB); c++) 
       if((posSearch->posMatrix[i][c] == GBLAST_SCORE_MIN) ||
             (0.0 == posSearch->posMatchWeights[i][c])) fprintf(fp,"  -  ");
         else fprintf(fp,"%.3f ",(posSearch->posMatchWeights[i][c]));
     	   fprintf(fp," %4d",
		posSearch->posExtents[i].rightExtent 
		- posSearch->posExtents[i].leftExtent +1);
   } fprintf(fp,"\n\n");
#endif
// #endif
     fprintf(fp,"\nLast position-specific scoring matrix computed\n");
     fprintf(fp,"         ");
     for(c = 0; c<= nAlpha(AB); c++) fprintf(fp,"  %c",AlphaChar(c,AB));
     for(i=0; i < length; i++) {
       fprintf(fp,"\n%5ld %c   ", (Int4) (i + 1), AlphaChar(q[i],AB));
       for(c = 0; c<= nAlpha(AB); c++) 
	 if(posSearch->posMatrix[i][c] == GBLAST_SCORE_MIN)
	   	fprintf(fp,"-I ");
	 else fprintf(fp,"%2ld ", (Int4) posSearch->posMatrix[i][c]);
     }
     fprintf(fp,"\n\n");
     fprintf(fp,"                      K         Lambda\n");
     fprintf(fp,"Standard Ungapped    %6.4f     %6.4f\n",
		compactSearch->kbp_std[0]->K,compactSearch->kbp_std[0]->Lambda);
     fprintf(fp,"Standard Gapped      %6.4f     %6.4f\n",
		compactSearch->kbp_gap_std[0]->K,compactSearch->kbp_gap_std[0]->Lambda);
     fprintf(fp,"PSI Ungapped         %6.4f     %6.4f\n",
		compactSearch->kbp_psi[0]->K,compactSearch->kbp_psi[0]->Lambda);
     fprintf(fp,"PSI Gapped           %6.4f     %6.4f\n",
		compactSearch->kbp_gap_psi[0]->K,compactSearch->kbp_gap_psi[0]->Lambda);
}

void posPrintInformation(psim_typ posSearch, bsb_typ search, Int4 passNum)
{
  Int4 querySize;
  Int4 c,n;

#if 0
  querySize = search->query_length;
  h_type HG=Histogram("Information content by position",0,querySize,2.0);
  // h_type HG=Histogram("Information content by position",0,querySize,1.0);
  for(c = 0; c < querySize; c++) {
	n = (Int4) (posSearch->posInformation[c]*100.0);
	if(n > 0) IncdMHist((double)c+1,n ,HG);
  } PutHist(stderr,60,HG); NilHist(HG);
#ifdef POSIT_DEBUG	// Used ifdef until final decision is made on output.
  printf("\nInformation content by position for pass %d\n", passNum);
  for(c = 0; c < querySize; c++) printf(" %5d", c); 
  printf("\n");
  for(c = 0; c < querySize; c++) printf(" %5.2lf", posSearch->posInformation[c]); 
  printf("\n");
#endif
#endif
}   
 
void posInitializeInformation(psim_typ posSearch, bsb_typ search)
{
  Int4		c, a, alphabetSize;
  sbp_typ	sbp;
  double	*stdrfp; 	// standard frequencies
  double	lambda;
  double	term1, term2, term3, term4;
  double	infoSum;
 
  Int4 querySize = search->query_length;
  unsigned char *query = search->query_sequence;
  NEW(posSearch->posInformation,querySize+1,double);
  alphabetSize = search->sbp->alphabet_size;
  // Compute standard frequencies as in GBlastScoreBlkFill in blastkar.c
  sbp = search->sbp;
  stdrfp = GBlastResFreqNew(sbp);
  GBlastResFreqStdComp(sbp,stdrfp); 
  lambda = search->sbp->kbp[0]->Lambda;
  for(c = 0; c < querySize; c++) {
    for(infoSum = 0,a=0; a < alphabetSize; a++){
      if(stdrfp[a] > posEpsilon) {
        term1 = search->sbp->matrix[query[c]][a];
	term2 = term1 * lambda;
	term3 = exp(term2); term4 = stdrfp[a] * term3;
	infoSum += term4 * log(term4/stdrfp[a])/NCBIMATH_LN2;
      }
    } posSearch->posInformation[c] = infoSum;
  }
  free(stdrfp);
}

/* Is this function used? */
void posFreeInformation(psim_typ posSearch) { free(posSearch->posInformation); }

csi_typ	*compactSearchNew(bsb_typ search)
// Copy a few fields from the large record search into the small record
// compactSearch, so that a small amount of information is passed into posit.c
{
   double	*stdrfp;	// gets standard frequencies in prob field
   csi_typ	*compactSearch;

   NEW(compactSearch,1,csi_typ);
   compactSearch->query = search->query_sequence;
   compactSearch->qlength = search->query_length;
   compactSearch->alphabetSize = search->sbp->alphabet_size;
   compactSearch->pseudoCountConst = search->pseudoCountConst;
   compactSearch->ethresh = search->ethresh;
   compactSearch->lambda =  search->sbp->kbp[0]->Lambda;
   compactSearch->matrix = search->sbp->matrix;
   compactSearch->kbp_psi = search->sbp->kbp_psi;
   compactSearch->kbp_gap_psi = search->sbp->kbp_gap_psi;
   compactSearch->kbp_std = search->sbp->kbp_std;
   compactSearch->kbp_gap_std = search->sbp->kbp_gap_std;
   compactSearch->lambda_ideal = search->sbp->kbp_ideal->Lambda;
   compactSearch->K_ideal = search->sbp->kbp_ideal->K;
   stdrfp = GBlastResFreqNew(search->sbp);
   GBlastResFreqStdComp(search->sbp,stdrfp); 
   NEW(compactSearch->standardProb,compactSearch->alphabetSize+1,double);
   for(Int4 a=0; a < compactSearch->alphabetSize; a++)
     compactSearch->standardProb[a] = stdrfp[a];
   stdrfp = GBlastResFreqDestruct(stdrfp);
   return(compactSearch);
}

void compactSearchDestruct(csi_typ * compactSearch)
// De-allocate memory for a record of type csi_typ
{ free(compactSearch->standardProb); free(compactSearch); }

// Some of the following checkpointing code is taken and adapted from
// code written by K. Shriram for FASTLINK.
// Reference:
//  A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. 
//  Avoiding Recomputation in Linkage Analysis,
//  Human Heredity 44(1994), pp. 225-237. */

#define  putCkptDouble(d,ckptFile)	(putCkptNumber(&(d),sizeof(double),ckptFile))
#define  putCkptInt4(i,ckptFile)	(putCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  putCkptChar(c,ckptFile)	(putCkptNumber(&(c),sizeof(Char),ckptFile))
 
#define MAXALLOC        0x40000000	// Largest permissible memory request

static size_t FileWrite(const void *ptr, size_t size, size_t n, FILE *stream)
//   FileWrite(buf, size, fp)
{
  if(n && (MAXALLOC/n) < size) print_error("FileWrite: size > SIZE_MAX");
  if(!ptr || !stream) return 0;
  size_t cnt; cnt = fwrite(ptr,size,n,stream);
  if (cnt != n) print_error("FileWrite error");
  return cnt;
}

static void  putCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
// General routine for putting the internal representation of a number. 
{ FileWrite(numberPtr,numberSize,1,ckptFile) ;  }

static void    putCkptFreqMatrix (double **theMatrix, Int4 length, FILE * ckptFile,a_type AB)
// Code to put a matrix, vector-by-vector. 
{
  for(Int4 i=0; i < length ; i++ ){
   for(Int4 vectorRef = 1; vectorRef <= nAlpha(AB) ; vectorRef++)
     putCkptDouble(theMatrix[i][vectorRef],ckptFile);
    // Code to put a vector of frequencies; put only the interesting entries
  }
}

static size_t FileRead(void *ptr, size_t size, size_t n, FILE *stream)
//   FileRead(buf, size, fp)
{
  if(n && (MAXALLOC/n) < size) print_error("FileRead: size > SIZE_MAX"); 
  if(!ptr || !stream) return 0;
  return fread(ptr,size,n,stream);
}

static void  getCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
// General routine for getting the internal representation of a number.
{ FileRead(numberPtr,numberSize,1,ckptFile); }

static void    getFreqVector(double * theVector, Int4 length, FILE * ckptFile,a_type AB)
{
  Int4	vectorRef ;
  for(vectorRef = 0; vectorRef < length; vectorRef++) theVector[vectorRef] = 0;
  for(vectorRef = 1; vectorRef <= nAlpha(AB); vectorRef++)
    getCkptDouble(theVector[vectorRef],ckptFile) ;
}

static void getCkptFreqMatrix(double **theMatrix, Int4 length, Int4 width, 
		FILE * ckptFile,a_type AB)
// Code to frequency matrix, vector-by-vector.
{ for(Int4 i=0; i < length; i++) getFreqVector(theMatrix[i],width,ckptFile,AB); }

Boolean posTakeCheckpoint(psim_typ posSearch, csi_typ *compactSearch, CharPtr fileName)
{
  FILE *checkFile; /*file in which to take the checkpoint*/
  checkFile = open_file(fileName,"","w");
  Boolean rtn=posTakeCheckpoint(posSearch, compactSearch, checkFile);
  fclose(checkFile);
  return rtn;
}

/* Take a checkpoint at the end of the current PSI-GBLAST round, stores
 query length, query, and position-specific target frequencies.
 Returns TRUE if checkpoint was sucessful and FALSE otherwise. */
Boolean posTakeCheckpoint(psim_typ posSearch, csi_typ *compactSearch, FILE *checkFile)
{
  Int4 length; /*length of query sequence, and an index for it*/
  Int4 i; /*indices to position and alphabet */
  Char localChar; /*temporary character*/

  length = compactSearch->qlength;
  putCkptInt4(length,checkFile);
  for(i = 0; i < length; i++) {
    localChar = AlphaChar(compactSearch->query[i],posSearch->AB);
    putCkptChar(localChar, checkFile);
  }  
#if 0	// DEBUG
  for(i = 0; i < length; i++) {
 std::cerr << i; std::cerr << ":  ";
  	for(Int4 c = 1; c < compactSearch->alphabetSize; c++) {
	 std::cerr << posSearch->posFreqs[i][c]; std::cerr << " ";
  	} std::cerr << std::endl;
  } std::cerr << std::endl; std::cerr << std::endl;
#endif
  putCkptFreqMatrix(posSearch->posFreqs,length,checkFile, posSearch->AB);
  return(TRUE);
}

Boolean posReadCheckpoint(psim_typ posSearch, csi_typ * compactSearch, CharPtr fileName)
{
  FILE *checkFile; 
  Boolean rtn;
  checkFile = open_file(fileName,"","r");
  rtn = posReadCheckpoint(posSearch, compactSearch, checkFile);
  fclose(checkFile);
  return rtn;
}

// Read a checkpoint from the end of a previous PSI-GBLAST round, get
// query length, query, and position-specific target frequencies.
// Returns TRUE if checkpoint was read sucessfully and FALSE otherwise.
Boolean posReadCheckpoint(psim_typ posSearch, csi_typ * compactSearch, FILE *checkFile)
{
  Int4 length1, length2, c; // length of query sequence, and an index for it
  Char  nextRes; 	    // next residue in stored copy of the query sequence
  Uint1Ptr oldQuery; 	    // array to hold the query sequence*/

  // fprintf(stderr,"Attempting to recover data from previous checkpoint\n");
  length1 = compactSearch->qlength;
  getCkptInt4(length2,checkFile);
  if(length1 != length2) {
    fclose(checkFile); 
    fprintf(stderr,"length query(%d) != length checkpoint(%d)\n",length1,length2);
    print_error("posReadCheckpoint: Failed to recover data\n");
  }
  NEW(oldQuery,length1+1,unsigned char);
  for(c=0; c < length1; c++) {		// check query sequence...
    getCkptChar(nextRes, checkFile);
    oldQuery[c] = AlphaCode(nextRes,posSearch->AB);
    if((oldQuery[c] != compactSearch->query[c]) && (oldQuery[c] != Xchar)) {
      fprintf(stderr,"Checkpoint[%d] = '%c' vs '%c' (alignment)\n",
			c,AlphaChar(oldQuery[c],posSearch->AB),
			AlphaChar(compactSearch->query[c],posSearch->AB));
				
      fprintf(stderr,
	"(oldQuery[%d] != compactSearch->query[%d]) && (oldQuery[%d] != Xchar)\n",c,c,c);
      free(oldQuery); fclose(checkFile);
      print_error("posReadCheckpoint: Failed to recover data\n");
    }
  }
  NEWP(posSearch->posMatrix,length1 + 1,Int4);
  NEWP(posSearch->posPrivateMatrix,length1 + 1,Int4);
  NEWP(posSearch->posFreqs,length1 + 1,double);
  for(Int4 i = 0; i <= length1; i++) {
    NEW(posSearch->posMatrix[i],compactSearch->alphabetSize+2,Int4);
    NEW(posSearch->posPrivateMatrix[i],compactSearch->alphabetSize+2,Int4);
    NEW(posSearch->posFreqs[i],compactSearch->alphabetSize+2,double);
  }
  getCkptFreqMatrix(posSearch->posFreqs,length1,compactSearch->alphabetSize,checkFile,
	posSearch->AB);
#if 0	// DEBUG
  for(i = 0; i <= length1; i++) {
 std::cerr << i; std::cerr << ":  ";
  	for(c = 1; c < compactSearch->alphabetSize; c++) {
	 std::cerr << posSearch->posFreqs[i][c]; std::cerr << " ";
  	} std::cerr << std::endl;
  } std::cerr << std::endl; std::cerr << std::endl;
#endif
  posFreqsToInformation(posSearch,compactSearch);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  FixposMatrix(posSearch, compactSearch); // AFN 11/7/01 fix...
  // fprintf(stderr,"posReadCheckpoint: Data recovered successfully\n");
  free(oldQuery); 
  return(TRUE);
}

//********************* NEW: Marginal Probability PSI-BLAST ************************

Int4	**CposComputation(psim_typ posSearch, bsb_typ search, csi_typ *compactSearch,
	sap_typ listOfGSeqAligns, Char *ckptFileName,double ***MatrixMP)
{
  search->posConverged = FALSE;
  posAllocateMemory(posSearch, compactSearch->alphabetSize,
                compactSearch->qlength, search->hitlist_count);
  posDemographics(posSearch, compactSearch, listOfGSeqAligns);
  posPurgeMatches(posSearch, compactSearch);
  posComputeExtents(posSearch, compactSearch);

  posComputeSequenceWeights(posSearch, compactSearch,MatrixMP);

  posCheckWeights(posSearch, compactSearch);
  posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, TRUE);
  if(NULL != ckptFileName) posTakeCheckpoint(posSearch,compactSearch,ckptFileName);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  FixposMatrix(posSearch, compactSearch); // AFN 11/7/01 fix...
  return posSearch->posMatrix;
}

void	posCalcMatchWeights(double *MatrixMP,double **posMatchWeights,Int4 qplace,
	double posA,a_type A)
// WARNING: assumes that MatrixMP has been renormalized after noise filtering
// That is, assumes that all columns sum to 1.0.
{
	unsigned char	letter;
	// 1. find all residues in subject seq. with prob >= 0.10 
#if 1
	double	total=0.0;
	for(letter=1; letter <= nAlpha(A); letter++){
	  if(MatrixMP[letter] > 0.0){
	        posMatchWeights[qplace][letter] += (MatrixMP[letter])*posA;
		total+= MatrixMP[letter];
	  }
	}
	assert(total < 1.01 && total > 0.99);
	// if(total > 1.01 || total < 0.99) print_error("posCalcMatchWeights() error");
#endif
#if 0
	for(letter=1; letter <= nAlpha(A); letter++){
	  if(MatrixMP[letter] > 0.0){
	        posMatchWeights[qplace][letter] += (MatrixMP[letter])*posA;
	  }
	}
#endif
}

void	posComputeSequenceWeights(psim_typ posSearch, csi_typ *compactSearch,
	double ***MatrixMP)
// Compute Henikoff weight of each sequence and letter in each position
// Call as:
//	posComputeSequenceWeightsMP(posSearch, compactSearch, (void *)posMatrices,
//		posCalcMatchWeights);
{
   Int4 length; /*length of query*/
   Int4 numseq, seqIndex; /*number of matches, index for them*/
   Int4  i; /*index over a multi-alignment block*/
   Int4 qplace; /*index into query*/
   double Sigma; /*Number of different characters occurring in matches within
                   a multi-alignment block, excluding identical columns*/
   double intervalSigma; /*Same as Sigma but includes identical columns*/
   Int4 alphabetSize; /*number of characters in alphabet*/
   Int4 *participatingSequences; /*array of participating sequences at a position*/
   Int4 *oldParticipatingSequences; /*array of participating sequences at a position*/
   Int4 posLocalVariety;  /*number of different characters at a position*/
   Int4 *posLocalC; /*counts of how many of each letter in this column*/
   Int4 c;
   Int4 thisSeq;
   Int4 numParticipating; /*number of sequences in this alignment block*/
   Int4 oldNumParticipating; /*number of sequences in this alignment block*/
   Boolean newSequenceSet; 
   Int4 p; /*index on sequences*/

   alphabetSize = compactSearch->alphabetSize;
   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   NEW(participatingSequences, numseq+1,Int4);
   NEW(oldParticipatingSequences,numseq+1,Int4);
   NEW(posLocalC,alphabetSize +2, Int4);
   for(qplace=0; qplace < length; qplace++) posSearch->posSigma[qplace] = 0.0;
   numParticipating = 0;
   for(qplace = 0; qplace < length; qplace++) {
     if((posSearch->posCount[qplace] > 1) 
	&& (posSearch->posIntervalSizes[qplace] > 0)) {
       oldNumParticipating = numParticipating;
       for(p =0; p < numParticipating; p++)
           oldParticipatingSequences[p] = participatingSequences[p];
       numParticipating = 0;
       for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
         if (!posSearch->posUseSequences[seqIndex]) continue;
	 if ((posSearch->pde_typMatrix[seqIndex][qplace].used) &&
	     (posSearch->pde_typMatrix[seqIndex][qplace].letter != GAP_CHAR)) {
	   participatingSequences[numParticipating] = seqIndex;
	   numParticipating++;
	 }
       }
       newSequenceSet = TRUE;
       if (numParticipating == oldNumParticipating) {
         for(p = 0; p < numParticipating; p++)
           if (oldParticipatingSequences[p] != participatingSequences[p])
             break;
         if (p == numParticipating) newSequenceSet = FALSE;
       }
       if(newSequenceSet) {
	 Sigma = 0; intervalSigma = 0;
	 for(seqIndex = 0; seqIndex <= numseq; seqIndex++) {
	   if (!posSearch->posUseSequences[seqIndex]) continue;
	   posSearch->posRowSigma[seqIndex] = 0.0;
	   posSearch->posA[seqIndex] = 0.0;
	 }
	 for(i = posSearch->posExtents[qplace].leftExtent;
	      i <= posSearch->posExtents[qplace].rightExtent; i++) {
	   posLocalVariety = 0;
	   for(c = 0; c <= alphabetSize; c++) posLocalC[c] = 0;
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     // used to check for GAP here
	     if(0 == posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter])
	       // letter (not a gap) not seen before in this query pos.
	       posLocalVariety++;  
	     posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]++;
	   }
	   intervalSigma += posLocalVariety;
	   if (posLocalVariety > 1) {
	     Sigma += posLocalVariety;
	     for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	       thisSeq = participatingSequences[seqIndex];
	       // used to check for gap here
	       posSearch->posRowSigma[thisSeq] += 
		 ( 1.0 / (double ) posLocalC[posSearch->pde_typMatrix[thisSeq][i].letter]);
	     }
	   }
	 }
       }
       if (Sigma > 0) {
	 for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = posSearch->posRowSigma[thisSeq]/Sigma;
	 }
       } else {
         for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = ((double ) 1 / (double ) numParticipating);
         }
       }
       posSearch->posSigma[qplace] = intervalSigma;
       for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 if(MatrixMP==0){
	   posSearch->posMatchWeights[qplace][posSearch->pde_typMatrix[thisSeq][qplace].letter]
                += posSearch->posA[thisSeq];
	 } else { 
	   posCalcMatchWeights(MatrixMP[thisSeq][qplace],
		posSearch->posMatchWeights,qplace,posSearch->posA[thisSeq],
		posSearch->AB);
	 }
       }
     }
   }
   free(participatingSequences);
   free(oldParticipatingSequences);
   free(posLocalC);
}

