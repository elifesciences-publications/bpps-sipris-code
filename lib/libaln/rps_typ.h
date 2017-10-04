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

#if !defined(_RPS_TYP_)
#define _RPS_TYP_
#include "dom_typ.h"

/************************ Rich Protein Sequence Type ************************
Data type for storing information about protein sequences, such as:
    1. Domain locations.
    2. Subseq locations within full sequences (for repeats).
    3. Taxonomic information.
    4. Protein family tree(?) (phylip).
    5. PDB structural information (including: *.pdb, *.ss, *.exp, *.rrc *.rpc)
    6. ???

// FullSeq information for SubSeq repeats.
void    AddFullCountsCMSA(ss_type FullSeq,unsigned short *FullRpts,cma_typ cma);
void    CopyFullCountsCMSA(cma_typ cma2, cma_typ cma);
void    RmFullCountsCMSA(cma_typ cma);

// Domain information for FullSeqs
void    AddDomainsCMSA(dom_typ *domains, cma_typ cma);
void    CopyDomainsCMSA(cma_typ cma2, cma_typ cma);
void    RmDomainsCMSA(cma_typ cma);

 *******************************************************************************/

class rps_typ {
public: 
	rps_typ( ){ init(); }
	rps_typ(FILE *);
	rps_typ& operator=(const rps_typ&);     // assignment operator
	~rps_typ( ){ Free(); }
	void	Put(FILE *,BooLean *,char *);
	void	Put(FILE *fp,BooLean *skip){ Put(fp,skip,0); }
	void	Put(FILE *fp){ Put(fp,0,0); }
	Int4	SqStart(UInt4 s, unsigned short d);
	Int4	SqEnd(UInt4 s, unsigned short d);
	float	Evalue(UInt4 s, unsigned short d);
	char	*DomainName(UInt4 s, unsigned short d);
	Int4	NumDomains(UInt4 s);
	UInt4	NumSeq( ){ return NumSeqs; }
private:
	void		copy(const rps_typ& dom);
	void		init( );
	void		Free();		// free memory...
//      Information about full sequences for subseqs.
        ss_type         FullSeq;        // full size sequences;
        unsigned short  *FullRpts;      // number of repeats in each full seq.
        Int4            *SubToFull;     // map subseq to fullseq.
//      Information about other domains withing sequences.
        dom_typ         *Domains;
};

#endif

