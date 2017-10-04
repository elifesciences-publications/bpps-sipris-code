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

#if !defined(_STPB_TYP_)
#define _STPB_TYP_

#include "stdinc.h"
#include "afnio.h"


#define BigNegative 			-99999999
#define exp2(x) 			(exp((x) * 0.69314718))
#define log2(x) 			(1.44269504*log(x))
#define BitsToProb(sc,nll,bts)		((nll)*exp2(((double) (sc))/(bts)))
#define ProbToBits(prob,nll,bits)	(prob != 0) ? ((Int4)(floor(0.5+(bits)*(log2((prob)/(nll)))))) : BigNegative



class stpb_typ {
public: 
			stpb_typ( ){ assert(!"Illegal constructor"); }
			stpb_typ(Int4, double *, double *, double *, double *,
			        double *, double *, double *, double *, double *, 
				double, double, double, double, double, double, double, 
				double, double, double, double, double,Int4);
			stpb_typ(Int4,Int4 *,Int4 *,Int4 *,Int4 *,Int4 *,Int4 *,
				Int4 *,Int4 *,Int4 *,Int4,Int4,Int4,Int4,Int4,Int4,Int4,
					Int4,Int4,Int4,Int4,Int4,Int4);
			stpb_typ(Int4, Int4 *M2M, Int4 *M2I, Int4 *M2D, Int4 *I2I,
				Int4 *I2M,Int4 *D2D,Int4 *D2M,double FreqOfMatchAtPositionOne);
			stpb_typ(stpb_typ *,Int4);
			stpb_typ& operator=(const stpb_typ&);	// assignment operator.
			~stpb_typ(){ Free(); }
	void		Put(FILE *);
	Int4		Length( ){ return length; }
	Int4		*MatToMat( ){ return m2m; }
	Int4		*MatToIns( ){ return m2i; }
	Int4		*MatToDel( ){ return m2d; }
	Int4		*InsToIns( ){ return i2i; }
	Int4		*InsToMat( ){ return i2m; }
	Int4		*DelToDel( ){ return d2d; }
	Int4		*DelToMat( ){ return d2m; }
	Int4		*BegToMat( ){ return b2m; }
	Int4		*MatToEnd( ){ return m2e; }
	Int4		MatToMat(Int4 i){ return m2m[i]; }
	Int4		MatToIns(Int4 i){ return m2i[i]; }
	Int4		MatToDel(Int4 i){ return m2d[i]; }
	Int4		InsToIns(Int4 i){ return i2i[i]; }
	Int4		InsToMat(Int4 i){ return i2m[i]; }
	Int4		DelToDel(Int4 i){ return d2d[i]; }
	Int4		DelToMat(Int4 i){ return d2m[i]; }
	Int4		BegToMat(Int4 i){ return b2m[i]; }
	Int4		MatToEnd(Int4 i){ return m2e[i]; }
	Int4		NB() { return nb; }	// N -> B = XT = 8 special transitions.
	Int4		NN() { return nn; }	// N -> N  (Rean Plan 7 documentation)
	Int4		EC() { return ec; }	// E -> C
	Int4		EJ() { return ej; }	// E -> J
	Int4		CT() { return ct; }	// C -> T
	Int4		CC() { return cc; }	// C -> C
	Int4		JB() { return jb; }	// J -> B
	Int4		JJ() { return jj; }	// J -> J
	Int4		BM() { return bm; }	// B -> M = begin to match; B -> I always *
	Int4 		BD() { return bd; }	// B -> D = begin to deletion.
	Int4 		GG() { return gg; }	// G -> G = NULT = null model trans. prob.
	Int4 		GF() { return gf; }	// G -> F
	Int4		Bits(){ return bits; }
	stpb_typ 	*ReverseStpb();
	stpb_typ 	*RenormalizeStpb();
	stpb_typ 	*ConcatenateStpb(Int4 rpts);
	stpb_typ 	*TempStpb(double temp);
private:
	void		Init(Int4 ,Int4 *,Int4 *,Int4 *,Int4 *,Int4 *,Int4 *,
				Int4 *,Int4 *,Int4 *,Int4,Int4,Int4,Int4,Int4,Int4,Int4,
					Int4,Int4,Int4,Int4,Int4,Int4);
	void		Free();
	unsigned short 	length;
	Int4		*m2m,*m2i,*m2d;
	Int4		*i2i,*i2m;
	Int4		*d2d,*d2m;
	Int4		*b2m,*m2e;
	Int4		nb,nn,ec,ej,ct,cc,jb,jj;
	Int4		bm,bd;
	Int4 		gg,gf;
	Int4		bits;
};

#endif

