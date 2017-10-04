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

#if !defined(_TXS_TYP_)
#define _TXS_TYP_

#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "sequence.h"
#include "seqset.h"
#include "residues.h"

#define MAX_TAX_GROUPS  500

class txc_typ {		// taxonomic category type.
public:
	txc_typ( );
	char      	**Group[MAX_TAX_GROUPS];
	char            Kingdom[MAX_TAX_GROUPS];
	Int4    	GroupNum;
	UInt4   GroupHits[MAX_TAX_GROUPS];
	~txc_typ( ){ Free(); }
private:
	void	Free();
	char	**CopyStr(const char **str);
};

class txs_typ {
public:
		txs_typ( ){ assert(!"Illegal constructor"); }
		txs_typ(char *,a_type);
		txs_typ(char *,Int4 *,Int4,const char **,char *,UInt4,e_type*,a_type);
		~txs_typ( ){ Free( ); }
	void	RmSubseq( );
	Int4	TotalSeqs( ){ return total_seqs; }
	Int4	NumTaxGroups( ){ return num_groups; }
	char	*GroupName(Int4 i){
		  if(i < 1 || i > num_groups) return 0;
		  else return group_name[i];
		}
	ss_type	Group(Int4 i){
		if(i < 1 || i > num_groups) return 0;
		else return group[i];
	}
	char	Kingdom(Int4 g){
		  if(g < 1 || g > num_groups) return 0;
		  else return kingdom[g];
		}
	e_type	QuerySeq(){ return SeqSetE(1,group[1]); }
	void	WriteFile(char *);
	void	PutSeqs(FILE *);
private:
	void	Free( );
	Int4	ReadFile( );
	void	init();
	void	Load(e_type *Seq);
	
	a_type	AB;
	Int4	num_groups;
	char	*filename;
	ss_type	*group;
	char	**group_name;
	// NEW: P = Protista; B = Bacteria; M = Metazoan; G = Green Plants; 
	// F = Fungi; O = Other; U = Unknown
	char	*kingdom;
	Int4	*group_size;
	Int4	total_seqs;
};

#endif

