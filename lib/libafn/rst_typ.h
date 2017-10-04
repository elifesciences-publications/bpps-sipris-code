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

#if !defined (_RST_TYP_)
#define _RST_TYP_
#include "residues.h"
#include "alphabet.h"
#include "set_typ.h"
#include "afnio.h"
#include "sset.h"
#include "blosum62.h"

class rst_typ { // Sequence dependence type...
public:
                rst_typ( ){ assert(!"Illegal constructor"); }
                rst_typ(char);
                rst_typ(char,a_type);
                ~rst_typ( ){ Free(); }
	sst_typ	**LegalResSets( ){ return ResidueSSet; }
	BooLean	IsLegalSet(sst_typ sst){
		// Is residue set sst legitimate?
		   sst_typ	isst;
		   Int4	rset;
		   for(unsigned char r =1; r <= nAlpha(AB); r++){
		      if(MemSset(r,sst)){
		   	for(rset=1; ResidueSSet[r][rset] != 0; rset++){
			    isst = IntersectSset(ResidueSSet[r][rset],sst);
			    if(isst == sst && isst == ResidueSSet[r][rset]){
				return TRUE;
			    }
			}
		      }
		   } return FALSE;
		}
	void	Put(FILE *,a_type ab=0);
	unsigned char	*NumResSets;
	sst_typ	**ResidueSSet;
	Int4	MaxResSet(){ return max_residue_set; }
	BooLean	**EqualsPrevious;	// Is Set[r][s] == Set[r'][s'] where r > r'?
private:
        void    Free(){
			for(char r = 0; r <= nAlpha(AB); r++){
				if(ResSetString[r]) free(ResSetString[r]);
				free(ResidueSSet[r]);
				if(EqualsPrevious[r]) free(EqualsPrevious[r]);
			} free(ResidueSSet); free(NumResSets); free(ResSetString);
			free(EqualsPrevious); 
			if(own_AB) NilAlpha(AB);	
		}
	a_type	AB;
	BooLean	own_AB;
	void	FindEquivalent();
	void	init_res_sets(char mode );
	void	init_res_sets(){ init_res_sets('G'); }
	void	ValidityCheck(FILE *fp);
	static const int	res_array_size=26;
	static const int	max_residue_set=24;
	void	Init_Res_Sets(const char *ResidueSets[21][res_array_size]);
	char	**ResSetString;
};

#endif

