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

#if !defined (_PTRN_TYP_)
#define _PTRN_TYP_

#include "alphabet.h"
#include "afnio.h"
#include "sset.h"

#if 0
class ptrn_typ {         // hpt pattern type.
public:
        	ptrn_typ( ){ assert(!"Illegal constructor"); }
        	ptrn_typ(Int4 n,Int4 len,sst_typ **sst, double **lpr){ 
			
		}
	BooLean	Empty( ){ return emptyHeap(dH); }
        ~ptrn_typ( ){ Free( ); }
private:
	Int4    N;      // number Columns.
        Int4    Len;    // number positions.
        double  **LPR;
        sst_typ **SST;  // SST[column][position];
};

#endif

typedef struct {
        Int4    N;      // number Columns.
        Int4    Len;    // number positions.
        double  **LPR;
        sst_typ **SST;  // SST[column][position];
} pattern_type;

typedef pattern_type        *ptrn_typ;

ptrn_typ MakePtrn(Int4 N, Int4 Len, sst_typ **sst, double **lpr);
void	 WritePtrn(FILE *fp, ptrn_typ ptrn);
void	 PutPtrn(FILE *fp,ptrn_typ ptrn, a_type AB);
void    PutPtrn(FILE *fp,ptrn_typ ptrn, Int4 i, a_type AB);
ptrn_typ ReadPtrn(FILE *fp);
void	 NilPtrn(ptrn_typ ptrn);


#endif
