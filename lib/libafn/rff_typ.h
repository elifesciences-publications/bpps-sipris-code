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

#if !defined(_RFF_TYP_)
#define _RFF_TYP_
#include "stdinc.h"

/************************* Rich Fasta Format Type ***************************

 *******************************************************************************/
class rff_typ {
public: 
			rff_typ( );
			rff_typ(FILE *);
			rff_typ(char *);
			rff_typ(const rff_typ *rff){ init(); copy2(rff); }
			rff_typ(const rff_typ *,UInt4,UInt4);
			rff_typ(const rff_typ& rff){ init(); copy(rff); }
		   	// constructor for 'rff_typ rff2=rff;' or 'rff_typ rff2(rff);
			~rff_typ( ){ Free(); }
			rff_typ& operator=(const rff_typ&);
	void		Put(FILE *fp){ Put(fp,1,len_ss); }
	void		Put(FILE *,UInt4,UInt4);
	UInt4	LengthS(){ return len_ss; }
	BooLean		IsStruct( ){ return (mode == 'S'); }
	char		Struct(UInt4 r){
			   if(r <= len_ss && ss) return ss[r]; else return 0; 
			}
private:
	void		Read(FILE *);
	void		Read(char *);
	void		init( );
	void    	copy(const rff_typ& rff);
	void    	copy2(const rff_typ *rff);
	void    	copy2(const rff_typ *rff,UInt4,UInt4);
	void		Free();		// free memory...
	char		*ss;    // secondary structure (h,s,c) state.
	UInt4   len_ss;  // length of secondary structure string.
	char		mode;	// S = secondary structure;
};

char	ScanOverRFF(FILE *fptr);

#endif

