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

#if !defined(_N2A_TYP_)
#define _N2A_TYP_

#include "sequence.h"
#include "alphabet.h"
#include "residues.h"
#include "gpsi_typ.h"

class n2a_typ {
public:
		n2a_typ(){ assert(!"Illegal n2a_typ constructor"); }
		n2a_typ(Int4 , Int4 , a_type , a_type );
		~n2a_typ(){ Free(); }
	Int4	Translate(e_type);
	e_type	BestTranslate(e_type , double *, e_type );
	e_type  *AllTranslate(Int4 *NumRFs,e_type DnaE);
	e_type	Sequence(Int4 ,Int4 );
	Int4	NumReadFrames( ){ return 6; }
	void	ResetCode(Int4 c){ reset( ); code = c; set_table( ); }
	Int4	SeqPerReadFrame(Int4 rf){
			return (((rf) > 0) && ((rf) <= 6) ? (num[(rf)]):0);
		}
private:
	void	Init(Int4 , Int4 , a_type , a_type );
	void	Free( );
	Int4	***set_table( );
	void	reset( );
	void	reverse_compl();
	Int4	translate(BooLean,e_type);
	unsigned char A,T,C,G,N;
	Int4	***table;
	a_type	AB, dnaAB;
	Int4	min_len;
	char	code;
	Int4	num[7];
	e_type	*E[7],dnaE;
};

#endif
