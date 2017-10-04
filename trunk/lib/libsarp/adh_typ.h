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

/* adh_typ.h - . */
#if !defined (ADH_TYP)
#define ADH_TYP
#include "cmsa.h"
#include "atom.h"
#include "sequence.h"
#include "histogram.h"
#include "residues.h"
#include "mheap.h"

class adv_typ {		// atomic distance variance type
public:
		adv_typ(){ print_error("not allowed"); }
		adv_typ(Int4 i, Int4 j, double v, double m, Int4 N){
			Init(i,j,v,m,N);
		}
		~adv_typ(){ Free(); }
	double	Variance(){ return variance; }
	double	Mean(){ return mean; }
	Int4	Number(){ return Num; }
	Int4	ColI(){ return I; }
	Int4	ColJ(){ return J; }
	void	Put(FILE *fp);
private:
	void	Init(Int4 i, Int4 j, double v, double m, Int4 N);
	void	Free();
	double	variance,mean;
	Int4	Num,I,J;
};

class adh_typ {		// adh_typ = atomic distance heap type
public:
		adh_typ(){ print_error("not allowed"); }
		adh_typ(Int4 h){ Init(h); }
		~adh_typ(){ Free(); }
	Int4    Insert(adv_typ *adv,double key);
	adv_typ	*DelMin(double *key, Int4 *Item);
	adv_typ	*DelMax(double *key, Int4 *Item);
	Int4	NumItems(){ return ItemsInMheap(mheap); }
private:
	mh_type mheap;
	Int4	hpsz;
	adv_typ	**ADV;
	void	Init(Int4 hs);
	void	Free();
};

#endif

