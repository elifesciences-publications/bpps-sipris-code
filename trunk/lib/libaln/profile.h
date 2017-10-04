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

#if !defined (PROFILE)
#define	PROFILE
#include "stdinc.h"
#include "alphabet.h"
#include "blosum62.h"
#include "smatrix.h"
/*************************** ADT PROFILE ***************************
	
	profile model data type

**********************************************************************/

/*************************** Profile Type **************************/
typedef struct {
	smx_typ	smx;
	a_type	A;
	Int4	*max;			/* maximum */
	Int4	*min;			/* minimum */
	Int4	K;			/* length */
	char	*maxseg;		/* maximum scoring segment */ 
	Int4	**cnts;			/* residue counts */
	Int4	nlet;
	Int4	nseqs;
	BooLean	calc_prof;		/* calculate profile */
} profile_type;
typedef profile_type *pfl_typ;

/******************************* private *****************************/
void	profile_error(char *s);
void	calc_profile(pfl_typ M);

/******************************* Public ******************************/
/***************************** operations ****************************/
pfl_typ	MkProfile(Int4 N, Int4 K, double *freq, a_type A);
void	NilProfile(pfl_typ M);
void	Add2Profile(unsigned char *seq, Int4 start, pfl_typ M);
void	RmProfile(unsigned char *seq, Int4 start, pfl_typ M);
void    PutProfile(FILE *fptr, pfl_typ M);
Int4     ScoreProfile(register unsigned char *seq,register Int4 start,register pfl_typ M);
double  ProfileProb(Int4 score, pfl_typ M);
/**************************** macro operations **********************/
#define kProfile(P)	((P)->K)
#define meanProfile(P)	(meanSMatrix((P)->smx))
#define sdProfile(P)	(sdSMatrix((P)->smx))
#define nSeqProfile(P)	((P)->nseqs)

#endif

