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

#if !defined (_PCH_TYP_)
#define _PCH_TYP_

#include "alphabet.h"
#include "afnio.h"
#include "dheap.h"

class pch_typ {         // pattern column heap type.
public:
        	pch_typ( ){ assert(!"Illegal constructor"); }
        	pch_typ(Int4 N,Int4 C){ Init(N,C); }
	Int4	Insert(Int4 i, keytyp key, Int4 Cat, Int4 Col) {
			assert(Cat > 0 && Cat <= NumCat);
			assert(Col > 0 && Col <= NumCol);
			assert(i > 0 && i <= hpsz);
			assert(!memHeap(i,dH));		// make sure that i not already in heap.
			insrtHeap(i,-key,dH); 
			Category[i]=Cat; Column[i]=Col;
		}
	Int4    DeleteMax(keytyp &key, Int4 &Cat, Int4 &Col) {
			if(emptyHeap(dH)){ key=0.0; Cat=Col=0; return 0; }
			key=-minkeyHeap(dH);
			Int4 i = delminHeap(dH); 
			assert(i > 0 && i <= hpsz);
			Cat = Category[i]; Col = Column[i];
			Category[i]=Column[i]=0;
			return i;
		}
	Int4	NumItems( ){ return ItemsInHeap(dH); }
	BooLean	Empty( ){ return emptyHeap(dH); }
        ~pch_typ( ){ Free( ); }
private:
        void    Init(Int4 N,Int4 C){
			NumCat=N; NumCol=C;
			hpsz = NumCat * NumCol;
			dH = dheap(hpsz+2,4);
			NEW(Category, hpsz + 2, Int4); NEW(Column, hpsz + 2, Int4);
		}
        void    Free( ){ free(Category); free(Column); Nildheap(dH); }
	dh_type	dH;
	Int4	hpsz;
	Int4	*Category;
	Int4	*Column;
	Int4	NumCat,NumCol;
	
};

#endif


