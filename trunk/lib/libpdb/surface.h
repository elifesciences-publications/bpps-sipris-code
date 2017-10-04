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

/****************** surface.h ***************/
#if !defined(SURFACE)
#define SURFACE
#include "afnio.h"
#include "stdinc.h"
#include "alphabet.h"
#include "residues.h"
#include "dheap.h"
#include "sequence.h"

/************************** surface datatype *****************************/
typedef struct {
	a_type	A;		/** residue alphabet **/
	Int4	len;		/** lenght of protein **/
	unsigned char	*seq;		/** sequence of exposed residues **/
	char	*exp;		/** exp[1..3]=buried,inter,surface **/
	char	filename[51];	/** filename of input **/
	double	*exposed;	/** solvent accessible area of residues **/
	double	*percent;	/** percent accessible area **/
} surface_type;

typedef	surface_type	*sfc_typ;

/******************************** private **********************************/
unsigned char    residue_surface(a_type A, char aa[]);
Int4    bubble_sort_surface(Int4 *L);
Int4    *LeeRichExposed(BooLean surface, sfc_typ S);
void    surface_error(char *s);

/******************************** PUBLIC ***********************************/
sfc_typ MkLeeRichards(Int4 maxres, char *filename, a_type A);
void    PutLeeRich(FILE *fptr, sfc_typ S);
void    NilLeeRich(sfc_typ S);
Int4    *SurfaceLeeRich(sfc_typ S);
Int4    *BuriedLeeRich(sfc_typ S);
Int4    CntExposedLeeRich(Int4 *surface, Int4 *buried, sfc_typ S);
double  ScoreDiffLeeRich(Int4 site, Int4 New, Int4 old, e_type E, void *X);
double  TotalScoreLeeRich(e_type E, sfc_typ S);
double  NegScoreDiffLeeRich(Int4 site, Int4 New, Int4 old, e_type E, void *X);
double  NegTotalScoreLeeRich(e_type E, sfc_typ S);
/******************************** MACROS ***********************************/
#define	PresentSurface(s,X)	((X)->seq[(s)] > 0)

#endif

