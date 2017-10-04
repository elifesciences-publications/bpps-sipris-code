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

#if !defined (_P2P_TYP_)
#define _P2P_TYP_
#include "stdinc.h"

class p2p_typ {		// profile-to-profile score type.
public:
	p2p_typ(Int4 lenC, Int4 n, Int4 *lenP, Int4 totWtC){
		LenC=lenC; nBlks=n;
		TotWtC=totWtC;
		NEW(LenP,nBlks+3,Int4); 
		for(Int4 b=1; b <= nBlks; b++) LenP[b]=lenP[b];
		NEWPP(MTX,LenC +3,Int4);
		for(Int4 c=1; c <= LenC; c++){
		    NEWP(MTX[c],nBlks+3,Int4);
		    for(Int4 b=1; b <= nBlks; b++) NEW(MTX[c][b],LenP[b] +3,Int4);
		}
	}
	~p2p_typ(){ 
		for(Int4 c=1; c <= LenC; c++){
		    for(Int4 b=1; b <= nBlks; b++) free(MTX[c][b]);
		    free(MTX[c]);
		} free(MTX); free(LenP);
	}
	void	SetScore(Int4 i, Int4 b, Int4 j, Int4 score){ 
			assert(b > 0 && b <= nBlks);
			assert(i > 0 && i <= LenC && j > 0 && j <= LenP[b]);
			MTX[i][b][j]=score;
		}
	Int4	RtnNumBlks(){ return nBlks; }
	Int4	RtnTotLenP(){
		   Int4 x=0;
		   for(Int4 b=1; b <= nBlks; b++) x+=LenP[b];
		   return x; 
		}
	Int4	RtnWtC(){ return TotWtC; }
	Int4	RtnLenC(){ return LenC; }
	Int4	RtnLenP(Int4 b){ return LenP[b]; }
	Int4	RtnScore(Int4 i, Int4 b, Int4 j){
			assert(b > 0 && b <= nBlks);
			assert(i > 0 && i <= LenC && j > 0 && j <= LenP[b]);
			return MTX[i][b][j];
		}
private:
	Int4	***MTX,TotWtC;
	Int4	LenC;
	Int4	*LenP,nBlks;
};

#endif
