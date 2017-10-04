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

#if !defined (_ESC_TYP_)
#define _ESC_TYP_

#include "alphabet.h"
#include "histogram.h"
// #include "dheap.h"
// #include "random.h"
// #include "residues.h"
// #include "blosum62.h"
#include "atom.h"
#include "dsets.h"
#include "pdb.h"
#include "sequence.h"
#include "mps_typ.h"
#include "swaln.h"

class esc_typ {		// equivalent (protein) structural chains type
// ME
public:
			esc_typ(){ print_error("esc_typ( ) constructor is disallowed"); }
			esc_typ(mps_typ *mpdb,Int4 min_seq_overlap){
				pdb=mpdb->pdb; NumFiles=mpdb->NumPDB; pdbSeq=mpdb->pdbSeq;
				AB=AminoAcidAlphabetPDB(pdb[1]); MinSeqOverlap=min_seq_overlap;
				CreateDsets(); 
			}
			~esc_typ(){ Free(); }
	Int4	**PDB_SetI;	// PDB_SetI[S][i] = I; 
	Int4	**PDB_SetC;	// PDB_SetC[S][i] = C; 
	Int4	*NumPDB_Set;	// NumPDB_Set[S] = number in each set.
	Int4	NumPDB_Sets;	// NumPDB_Sets = number of sets.
        e_type  *FullSeq;	// FullSeq[S].

	Int4	*BestInSet;	// number of best chain in set s

        e_type  **pdbSeq;	// created here, but passed back to calling environment.

	Int4	MinSeqOverlap;	// passed in 
private:
	void	Free();
	void	CreateDsets();	// Disjoint sets.  Clustering of pdb chains into equivalence sets.
	void	GetPDB_Seq();	// get the sequences corresponding to structures.
	
	// passed in as arguments but not owned by this object.
	pdb_typ *pdb;		// protein structures; Dont need to store this here?
	a_type	AB;
	Int4	NumFiles;	// number of structural files.
};

#endif

