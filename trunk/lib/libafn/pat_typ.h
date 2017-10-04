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

#if !defined (_PAT_TYP_)
#define _PAT_TYP_

#include "stdinc.h"
#include "afnio.h"
#include "sset.h"
#include "alphabet.h"

#define MAX_NUM_PATTERNS 500
#define MAX_NUM_PATTERN_POS 200

class pat_typ {         // pattern type 
public:
		pat_typ(){ print_error("pat_typ( ) constructor is disallowed"); }
		pat_typ(FILE *fp){ OpenPttrnFile(fp); }
		~pat_typ( ){ Free(); }
	void	Put(FILE *fp);
	pat_typ *SubPatterns(char *sets);
	void	SimilarPatterns(a_type AB);
	void	SimilarPatterns(FILE *fp, pat_typ *pat, a_type AB);
	double	SimilarPatterns(FILE *fp,double cutoff,Int4 I,Int4 J, pat_typ *pat, a_type AB);
	Int4	NumPttrns;
	Int4	*NumPttrnRes;	// NumPttrnRes[pttrn] = 15;
	char	***PttrnRes;	// PttrnRes[pttrn][pos] = "ILV";
	Int4	**PosPttrn;	// PosPttrn[pttrn][i] = 36;
	Int4	MaxPttrnPos;	// maximum value in PosPttrn matrix.
	Int4	*PttrnCategory;	// PttrnCategory[pttrn] = 13; // == column in hpt.
private:
	void	Free();
	void    OpenPttrnFile(FILE *fp);
};

#endif

