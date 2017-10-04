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

#if !defined(_TXN_TYP_)
#define _TXN_TYP_

#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "seqset.h"
#include "residues.h"

#if 0	//*************************************************************

	[1|2|3|4|...|i|...|max]

	max == 500,000? (12 + 4 + 1 = 17 bytes) (~200,000 currently).
		17 * 500,000 = 8.5 Mbytes.
#endif  //*************************************************************

class txn_typ {		// Taxonomy node type.
public:
                txn_typ( ){ assert(!"Illegal constructor"); }
                txn_typ(UInt4 Tax_id,char *Rank,UInt4 Parent);
                ~txn_typ( ){ ; }
private:
	UInt4	tax_id;		// node id in GenBank taxonomy database
	char		*rank;		// superkingdom, kingdom, ...
	UInt4	parent;		// parent node id in GenBank taxonomy database.
	char		**name;		// names for node.
	unsigned char	gc_code;	// genetic code id.

        // NEW: P = Protista; B = Bacteria; M = Metazoan; G = Green Plants;
        // F = Fungi; O = Other; U = Unknown
};

#endif

