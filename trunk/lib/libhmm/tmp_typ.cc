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

// tmp_typ.cpp: implementation of the tmp_typ class.
//
//////////////////////////////////////////////////////////////////////
/**************************
By Trevor Carlson
   01/13/00
   tcarlson@andrew.cmu.edu

With assistance from 
Dr. Neuwald neuwald@cshl.org and Aleksandar Poleksic poleksic@cshl.org
**************************/
#include "tmp_typ.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

tmp_typ::tmp_typ(a_type A)
{
  AB=A;
  transSet = new Int4 * [2];
  transSet[START] = new Int4[MAX_NUM_HELIX];
  transSet[END] = new Int4[MAX_NUM_HELIX];

  transLen = MAX_NUM_HELIX;
  
  transSet[START][0] = 0; 
  transSet[END][0] = 0;

  theArray=NULL;
}

tmp_typ::~tmp_typ()
{
  if ( transSet != NULL ) { delete [] transSet[0]; delete [] transSet[1]; delete [] transSet; }
  if ( theArray != NULL ) { delete [] theArray; }
}

tmp_typ::tmp_typ( Int4 maxNumHelices ,a_type A) 
{

  AB=A;
  transSet = new Int4 * [2];
  transSet[START] = new Int4[maxNumHelices];
  transSet[END] = new Int4[maxNumHelices];

  transLen = maxNumHelices;
  
  transSet[START][0] = 0; 
  transSet[END][0] = 0;

  theArray = NULL;

}

Int4** tmp_typ::getTrans (Int4 &numRegions, float *scores ,e_type E)
{
  Int4 i;

#if 1	// AFN
  char *seq;
  Int4 seqLen=LenSeq(E);
  char    *buffer = new char[seqLen+4];
  // float   *scores = new float[transLen+3];
  i=SeqToString(buffer, E, AB);
  // buffer[i]=buffer[i+1]=0; // AFN modification to attempt to fix bug.
  seq=buffer;
#endif

  curTransPos = -1;
  if ( theArray != NULL ) { delete [] theArray; }
  aSeq = seq;
  memScores = scores;
  if ( memScores != NULL ) {
    for (i=0;i<transLen;i++){ memScores[i] = 0.0; }
  }

#ifdef DEBUG
  printf("\nSeq Length: %d\n", seqLen);
#endif

  inputTrans( seq, seqLen); 
  computeTrans();

  Int4 ** membraneRegions = new Int4 *[2];

  membraneRegions[START] = new Int4 [curTransPos +1];
  membraneRegions[END] = new Int4 [curTransPos +1];

  for ( i = 0 ; i < curTransPos + 1 ; i++ ) {
      membraneRegions[START][i] = transSet[START][i] + 1;  // Adding one to convert from indexes to residue number.
      membraneRegions[END][i] = transSet[END][i] + 1;
  }

  numRegions = curTransPos + 1;  //Converting index to actual number.
#if 1	// AFN
  delete [] buffer;
#endif
  return membraneRegions;

}

void tmp_typ::inputTrans( char *seq, Int4 seqLen )
{
  Int4 i;
  arrayLen = seqLen;
#if 0
  theArray = new BooLean[seqLen]; 
#else // AFN: fixes Array bounds read error?
  theArray = new BooLean[seqLen+3];
  theArray[seqLen]=0; theArray[seqLen+1]=0; theArray[seqLen+2]=0; 
#endif

  for(i=0 ; i< arrayLen ; i++ ) {
    switch (seq[i]) {
    case 'A': theArray[i] = TRUE; break;
    case 'D': theArray[i] = FALSE; break;
    case 'C': theArray[i] = TRUE; break;
    case 'E': theArray[i] = FALSE; break;
    case 'F': theArray[i] = TRUE; break;
    case 'G': theArray[i] = TRUE; break;
    case 'H': theArray[i] = FALSE; break;
    case 'I': theArray[i] = TRUE; break;
    case 'K': theArray[i] = FALSE; break;
    case 'L': theArray[i] = TRUE; break;
    case 'M': theArray[i] = TRUE; break;
    case 'N': theArray[i] = FALSE; break;
    case 'P': theArray[i] = FALSE; break;
    case 'Q': theArray[i] = FALSE; break;
    case 'R': theArray[i] = FALSE; break;
    case 'S': theArray[i] = TRUE; break;
    case 'T': theArray[i] = TRUE; break;
    case 'V': theArray[i] = TRUE; break;
    case 'W': theArray[i] = TRUE; break;
    case 'Y': theArray[i] = TRUE; break;
    case 'X': theArray[i] = FALSE; break;
    case 'Z': theArray[i] = FALSE; break;
    default:
      printf("Error: invalid input character type, aborting.\n");
      exit(2);
      break;
    }
  }

#ifdef DEBUG
  Int4 j = 1;
  printf("\nBinary Output:\n");
  for ( i = 0 ; i < arrayLen ; i++ ) {
    printf("%d", theArray[i]);
    if ( ( j % 10) == 0 ) { printf(" "); }
    if ( ( j % 60) == 0 ) { printf("\n"); j = 0; }
    j++;   

  }
  printf("\n");
#endif
}

void tmp_typ::computeTrans()
{

  Int4 i;

  makeSetArray();
  
#ifdef DEBUG
  printf("\nOutput from Step #1\n");
  for ( i = 0 ; i <= curTransPos ; i++ ) {
    printf("Pos: %-4d  %d->%d\n", i, transSet[START][i], transSet[END][i] );
  } 
#endif

  trimHelix();

#ifdef DEBUG
  printf("\nOutput from Step #2\n");
  for ( i = 0 ; i <= curTransPos ; i++ ) {
    printf("Pos: %-4d  %d->%d\n", i, transSet[START][i], transSet[END][i] );
  } 
#endif

  if ( memScores != NULL ) {
    setupScores();
  }


  cleanUpTrans();

}

void tmp_typ::makeSetArray () 
{

  Int4 i, j, zeroCount;

  // Count through all of the Residues
  for ( i = 0 ; i <= arrayLen-18 ; i++ ) {

    // If there is a Int4 batch of "0"s, increment past them until a "1" is found
    for ( ; i < arrayLen && !theArray[i] ; i++ ) {}
    zeroCount = 0;
    for(j=i ; j<i+18 ; j++ ) { // Increments through 18 resides to see if it passes the 1st test
      if (!theArray[j]){               // Checks for a "0"
	if (!theArray[j+1]) { 	         // If there are consecutive "0"s, fail.
	  zeroCount = 4;
	  //i = j+1;                       // Sets i because any further checks in this area will result in an additional fail. Also, i will be incremented at the end of the loop.
	  break;
	}
	zeroCount++;	                 // If there aren't consecutive zero's continue.
      }
      if ( zeroCount > 3 ){            // Testing for more than 3 "0"s per 18 unit length
	break;
      }
    }
    if ( zeroCount <= 3 ) {          // Only add to the list if it meets the <= 3 zero criteria.
      if ( !unionSets(i) ) {           // Stop processing if it exceeds the max. number of transmembrane helices.
	return;
      }
    }
  }
}

BooLean tmp_typ::unionSets ( Int4 startPos ) {

  if ( curTransPos < 0 ) {
    curTransPos = 0;
    transSet[START][curTransPos] = startPos;
    transSet[END][curTransPos]   = startPos + 17;
  } else {
    if ( startPos <= transSet[END][curTransPos] ) {
      transSet[END][curTransPos] = startPos + 17;
    } else {
      curTransPos++;
      if ( curTransPos >= transLen ) {
	return FALSE;
      }
      transSet[START][curTransPos] = startPos;
      transSet[END][curTransPos]   = startPos + 17;
    }
  }
  return TRUE;

}

float tmp_typ::getVal ( char c ) {

  switch ( c ) {
    case 'A': return 0.2881819; break;
    case 'D': return -1.766092; break;
    case 'C': return 0.0601539; break;
    case 'E': return -1.783791; break;
    case 'F': return 0.5347374; break;
    case 'G': return -0.002002; break;
    case 'H': return -0.634878; break;
    case 'I': return 0.5894519; break;
    case 'K': return -2.162823; break;
    case 'L': return 0.4842763; break;
    case 'M': return 0.3457151; break;
    case 'N': return -0.705220; break;
    case 'P': return -0.574475; break;
    case 'Q': return -1.070025; break;
    case 'R': return -1.777866; break;
    case 'S': return -0.202116; break;
    case 'T': return -0.176737; break;
    case 'V': return 0.4693784; break;
    case 'W': return 0.2421656; break;
    case 'Y': return 0.1362776; break;
    case 'X': return 0; break;
    case 'Z': return 0; break;
    default:
      printf("Error: invalid input character type(2): %c, aborting.\n", c );
      return 0;
      //exit(3);
      break;
    }

}


float tmp_typ::getWindowVals( Int4 i ) {

  switch ( i ) {

  case 9: return -3.089584;
  case 10: return -3.809399;
  case 11: return -3.636127;
  case 12: return -3.706745;
  case 13: return -4.088112;
  case 14: return -3.469073;
  case 15: return -2.719837;
  case 16: return -2.308329;
  case 17: return -0.233416;
  case 18: return -0.577634;
  case 19: return 0.185573;
  case 20: return 0.648399;
  case 21: return 3.034351;
  case 22: return 0.523821;
  case 23: return 0.211600;
  case 24: return 0.257547;
  case 25: return -0.003938;
  case 26: return -0.740718;
  case 27: return -1.230001;
  case 28: return -1.506814;
  case 29: return -2.249833;
  case 30: return -3.038290;
  case 31: return -5.022422;
  case 32: return -4.162220;
  case 33: return -6.321705;
  case 34: return -5.022422;
  case 35: return -6.034022;
  case 36: return -5.810879;
  case 37: return -5.340875;
  case 38: return -4.375794;
  case 39: return -7.420317;
  case 40: return -6.727170;
  default: return -10.0;

  }
  
}

void tmp_typ::trimHelix() 
{
  Int4 cursorPos; //startPos, endPos, cursorPos, midPoint;

  for ( Int4 i = 0 ; i <= curTransPos ; i++ ) {
    //midPoint = ( transSet[START][i] + transSet[END][i] ) / 2;
    cursorPos = transSet[START][i];
    while ( cursorPos + 3 < transSet[END][i] ) {
      while ( theArray[cursorPos] && (cursorPos + 3 < transSet[END][i]) ) { cursorPos++; }
      if ( cursorPos + 3 > transSet[END][i] ) { continue; }
      if ( !theArray[cursorPos+1] || !theArray[cursorPos+2] || !theArray[cursorPos+3] ) {
	trimHelixHelper( i, cursorPos );
      } cursorPos++;
    }
  }
}


void tmp_typ::trimHelixHelper( Int4 transNum, Int4 cuttingReg ) 
{

  float holder;

#ifdef DEBUG
  printf("Membrane Reg #:%d, CuttingReg:%d\n", transNum, cuttingReg);
#endif

  if ( transSet[END][transNum] - transSet[START][transNum] <= 17 ) {
 

    // Mark region for deletion.
 
    transSet[START][transNum] = transSet[END][transNum];


  } else if ( transSet[END][transNum] - transSet[START][transNum] + 1 >= 40 ) {

    cutLongSegmentsHelper( transNum );

  } else {

    // This section is for normal length segments that need to be trimmed.

    holder = getBestLength( transSet[START][transNum], transSet[END][transNum] );
    
    if ( memScores != NULL ) {
      memScores[transNum] = holder;
    }
  }

}

float tmp_typ::getBestLength( Int4 &startPos, Int4 &endPos ) {

    Int4 j, k;
    float curScore, maxScore=0;

    Int4 bestStart, bestEnd;

    bestStart = startPos;
    bestEnd = endPos; 


    for ( k = 9 ; k <= 40 && k <= endPos; k++ ) { 
      curScore = 0;

      for ( j = startPos ; j < k + startPos && j <= endPos ; j++ ) {
	curScore+= getVal( aSeq[j] );
      }

      curScore+= getWindowVals( k );
      
      for ( j = startPos + 1 ; j <= endPos - k ; j++ ) {
	curScore-= getVal( aSeq[j - 1] );
	curScore+= getVal( aSeq[j + k] );
	if ( curScore > maxScore ) {
	  maxScore = curScore;
	  
	  bestEnd = j+k;
	  bestStart = j;
	}
	/*	
#ifdef DEBUG
	printf("Max Score:%f, StartPos:%d, EndPos:%d, Wind:%d, CurrentScore:%f\n", maxScore, bestStart, bestEnd, k, curScore );
#endif
*/
      }
    }

    startPos = bestStart;
    endPos = bestEnd;

    return maxScore;
}

void tmp_typ::cutLongSegments() {

  Int4 i;
  
  for ( i = 0 ; i <= curTransPos ; i++ ) {
    
    if ( transSet[END][i] - transSet[START][i] + 1 >= 40 ) {
      cutLongSegmentsHelper(i);
    }

  }

}

void tmp_typ::cutLongSegmentsHelper( Int4 transNum ) {

  float bestScoreLeft=0, bestScoreRight=0, tempScoreLeft, tempScoreRight;
  Int4 tempStart, tempEnd, tempMidL, tempMidR, i,
       bestStart, bestEnd, bestMidL, bestMidR;

  tempStart = transSet[START][transNum];
  tempEnd   = transSet[END]  [transNum];

  for (i=tempStart+17;i<tempEnd-17;i++) {

    tempMidL = i;
#ifdef DEBUG
    printf("\nLONG1:: %d, %d\n",tempStart, tempMidL);
#endif
    tempScoreLeft=getBestLength( tempStart, tempMidL );
#ifdef DEBUG
    printf("LONG2:: %d, %d\n",tempStart, tempMidL);
#endif
 
    tempMidR = i+1;
#ifdef DEBUG
    printf("\nLONG1:: %d, %d\n", tempMidR, tempEnd);
#endif
    tempScoreRight=getBestLength( tempMidR, tempEnd );
#ifdef DEBUG
    printf("LONG2:: %d, %d\n",tempMidR, tempEnd);
#endif

    if((tempScoreLeft+tempScoreRight)>(bestScoreLeft+bestScoreRight)) {

      bestStart = tempStart;
      bestEnd   = tempEnd;
      bestMidL  = tempMidL;
      bestMidR  = tempMidR;

      bestScoreLeft =tempScoreLeft;
      bestScoreRight=tempScoreRight;

    }

#ifdef DEBUG
    printf("bs:%d bml:%d bmr:%d be:%d\n", bestStart, bestMidL, bestMidR, bestEnd);
#endif

    tempStart = transSet[START][transNum]; //Reset the starting and ending positions
    tempEnd   = transSet[END]  [transNum];
  }

  
  for(i=(transLen-1<curTransPos?transLen-1:curTransPos);i>transNum;i--){
    transSet[START][i+1] = transSet[START][i];
    transSet[END][i+1] = transSet[END][i];
    
    if ( memScores != NULL ) {
      memScores[i+1] = memScores[i];
    }
  }


  transSet[START][transNum] = bestStart;
  transSet[END]  [transNum] = bestMidL;
  transSet[START][transNum+1]=bestMidR;
  transSet[END]  [transNum+1]=bestEnd;

  if ( memScores != NULL ) {
    memScores[transNum] = bestScoreLeft;
    memScores[transNum+1]=bestScoreRight;
  }

  curTransPos++;
  
  
}

/* // An older attempt to use a running score and large drops to determine the 2 TM_helix's
void cutLongSegmentsHelper2( Int4 transNum )
{

  Int4 i;
  float localMax=0, localMin=0, Score=0, oldMin=0, oldMax=0;

  i = transSet[START][transNum];


  for ( i = transSet[START][transNum] ; i <= transSet[END][transNum] ; i++ ) {

    Score += getVal( aSeq[i] );

    if ( Score > localMax ) {

      if ( (localMax-localMin) > (oldMax-oldMin) // && i > transSet[START][i]+9 && i < transSet[END][i] - 9  ) {
	if (i > transSet[START][i]+9 && i < transSet[END][i] - 9) {
	  oldMin = localMin;
	  localMax = localMin = Score;
	  oldMax = localMax;
	}
      } else {
	localMax = localMin = Score;
      }

    } else if ( localMin > Score ) {
      localMin = Score;
    }

    printf("Loc:%d Score:%5f LMin:%5f LMax:%5f OMin:%5f, OMax:%5f\n", i, Score, localMin, localMax, oldMin, oldMax );
    
  }

}
*/

void tmp_typ::setupScores() {

  Int4 i, j;
  float Score=0;

  for ( i=0;i<=curTransPos;i++ ) {
    if ( memScores[i] == 0.0 ) {
      for (j=transSet[START][i];j<=transSet[END][i];j++) {
	Score+= getVal( aSeq[j] );
      }
      memScores[i]=Score; Score=0;
    }
  }
}

void tmp_typ::cleanUpTrans() {

  Int4 i, j;
  
  for ( i = 0 ; i <= curTransPos ; i++ ) {
    if ( transSet[START][i] == transSet[END][i]  ) {
      for ( j = i + 1 ; j <= curTransPos && j < transLen ; j++ ) {
	transSet[START][j-1] = transSet[START][j];
	transSet[END]  [j-1] = transSet[END]  [j];
      } curTransPos--;
    }
  }
}
