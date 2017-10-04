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

#if !defined (_LGM_TYP_)
#define _LGM_TYP_
#include "stdinc.h"

class lgm_typ {		// precomputed lgamma values for digitized real values.
public:
	lgm_typ(UInt4 MWCs, char WtFactor){
		MaxWtCnts=MWCs;
		assert(WtFactor==wt_factor);
		NEW(this->LnGamma, MaxWtCnts+3, double);
#if 1
		for(UInt4 i=1; i<=MaxWtCnts; i++) this->LnGamma[i]=CalcLnGamma(i);
#else	// fix issue with deletions...
		UInt4	i,ZeroOut=1*wt_factor;
		for(i=1; i<=ZeroOut; i++) this->LnGamma[i]=0.0;
		for( ; i<=MaxWtCnts; i++) this->LnGamma[i]=CalcLnGamma(i);
#endif
	}
	~lgm_typ(){  free(this->LnGamma); }
	double	*RtnArray( ){ return LnGamma; }
	double	LGamma(UInt4 x){
			if(x <= MaxWtCnts) return LnGamma[x]; else return this->CalcLnGamma(x);
		}
	double  CalcLnGamma(UInt4 cnts) { return lgamma((double)cnts/(double)wt_factor); }
private:
	double  *LnGamma;       // array for computing lngamma quickly.
	static const char wt_factor=100;
	UInt4	MaxWtCnts;
};

#endif
