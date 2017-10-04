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

// tmp_typ.h: interface for the tmp_typ class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __PREDTRANSMEMBRANE_H__
#define __PREDTRANSMEMBRANE_H__

#include "stdinc.h"
#include "sequence.h"
#include "alphabet.h"

//#define DEBUG 1

#define MAX_NUM_HELIX 50
#define START 0
#define END 1

class tmp_typ  {
public:
	tmp_typ(a_type);
	tmp_typ( Int4 ,a_type);
	Int4 ** getTrans(Int4 &, float * ,e_type);
	virtual ~tmp_typ();

private:
	void computeTrans();

	void inputTrans( char *, Int4);
	void makeSetArray();
	BooLean	unionSets( Int4 );
	float getVal( char );
	float getWindowVals( Int4 );
	float getBestLength( Int4&, Int4& ); 
	void trimHelix();
	void trimHelixHelper( Int4, Int4 );
	void cutLongSegments();
	void cutLongSegmentsHelper( Int4 );
	void setupScores();
	void cleanUpTrans();
	
	// Pointer to the TM Helices and the number of TM Helices.
	Int4 ** transSet;
	Int4    transLen;
	a_type	AB;

	// Reference to the character array of the sequence.
	char*   aSeq;
	float * memScores;
	// Number of TM Helix regions that exist.
	Int4 curTransPos;
	// Array for determining if the A.A. is likely for the region.
	BooLean	*theArray;
	Int4 arrayLen;
};

#endif 
