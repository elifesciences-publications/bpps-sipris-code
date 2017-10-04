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

#if !defined (_MPS_TYP_)
#define _MPS_TYP_

#include "alphabet.h"
#include "histogram.h"
// #include "residues.h"
#include "pdb.h"
#include "atom.h"
#include "sequence.h"

#define MAX_NUM_PDB_INPUT_FILES 5000
#define	MAX_NUMBER_INTERNAL_RPTS 100

class mps_typ {         // multiple protein structure type
public:
		mps_typ(){ print_error("mps_typ( ) constructor is disallowed"); }
		mps_typ(char *infile){
			UseBeta=FALSE;
			pdb_infile=AllocString(infile);
			GetFileNamesPDB(); Read();
			GetResidues( );
		}
		~mps_typ( ){ Free(); }
	char	**pdb_file;	// names of file (with path).
	Int4	NumPDB;
	pdb_typ *pdb;
        e_type  **pdbSeq;	// direct pointer to esc->pdbSeq; ???
	char	*pdb_infile;	// name of pdb input file.
	atm_typ ***Calpha;
	res_typ ***ResALL;
	Int4    **num_resALL;
private:
	BooLean	UseBeta;
	void	GetResidues( );
	void	Read();
	void	GetFileNamesPDB();
	void	Free();
};

#endif

