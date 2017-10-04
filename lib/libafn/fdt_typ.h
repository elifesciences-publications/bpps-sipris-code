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

#if !defined (_FDT_TYP_)
#define _FDT_TYP_

#include "dheap.h"
#include "stdinc.h"
#include "set_typ.h"
#include "afnio.h"

class fdt_typ {		// Functional divergence table type.
public:
          fdt_typ( ){ assert(!"Illegal constructor"); }
	  fdt_typ(Int4 Ni,char **iR,char **iName,char *tosi,Int4 No, char **oR,char **oName,char *toso){
		NumORows=No; NumIRows=Ni; irow=iR; orow=oR; NameORow=oName; NameIRow=iName;
		NEWP(Hits,NumORows+3,Int4); 
		for(Int4 i=1; i<=NumORows; i++){ NEW(Hits[i],NumIRows+3,Int4); }
		TypeOfSetO=toso;
		TypeOfSetI=tosi;
		FindMisc();
	  }
	void	AssignHits(Int4 o, Int4 i,Int4 hits){
		assert(i >= 0 && i <= NumIRows);
		assert(o > 0 && o <= NumORows);
		// if(i==0) fprintf(stderr,"DEBUG: hits = %d\n",hits);
		Hits[o][i]=hits;
	}
	void	Put(FILE *fp);
	~fdt_typ( ){ Free( ); }
private:
	void	FindMisc();
	set_typ	*ParentNodeI;	
	set_typ	*ParentNodeO;	
	Int4	NumORows,NumIRows;
	char	*TypeOfSetO,*TypeOfSetI;
	char	**orow;	// output hpt rows.
	char	**irow;	// input hpt rows.
	char	**NameORow,**NameIRow;
	Int4	**Hits;
	void	Free();
};

#endif

