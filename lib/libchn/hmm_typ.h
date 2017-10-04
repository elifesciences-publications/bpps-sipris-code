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

#if !defined(_HHM_TYP_)
#define _HHM_TYP_

#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"

#define BigNegative                     -99999999
#define log2(x)                         (1.44269504*log(x))
// #define exp2(x)                         (exp((x) * 0.69314718))
// #define BitsToProb(sc,nll,bts)          ((nll)*exp2(((double) (sc))/(bts)))
// #define ProbToBits(prob,nll,bits)       (prob != 0) ? ((Int4)(floor(0.5+(bits)*(log2((prob)/(nll)))))) : BigNeg

class hmm_typ {
public:
	hmm_typ(FILE *fp,a_type A);
	hmm_typ(char *Name, Int4 len, Int4 **matemit, Int4 **insemit, Int4 *mm, Int4 *mi, Int4 *md, 
			Int4 *ii, Int4 *im, Int4 *dd, Int4 *dm, Int4 *bm, Int4 *me, a_type A, char mode='N'){
		  name=AllocString(Name); desc=0; AB=A;
		  Init(len, matemit, insemit,mm,mi,md,ii,im,dd,dm,bm,me,mode); 
		}
	~hmm_typ(){ Free(); }
	hmm_typ *Copy(){
		   hmm_typ *hmm= new hmm_typ(name,length,mat_emit,ins_emit,
						m2m,m2i,m2d,i2i,i2m,d2d,d2m,b2m,m2e,AB);
		   return hmm;
		}
	void	Put(FILE *fp);
private:
	void	Init(Int4 len, Int4 **matemit, Int4 **insemit, Int4 *mm, Int4 *mi, Int4 *md, 
			Int4 *ii, Int4 *im, Int4 *dd, Int4 *dm, Int4 *bm, Int4 *me,char mode);
	void	InitAsNull();
	void	Free();
	a_type	AB;
	char	*name,*desc;
	Int4	*DefaultInsEmit;
	Int4	**mat_emit,**ins_emit,*nule;
	Int4	length;
        Int4    *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        Int4    nb,nn,ec,ej,ct,cc,jb,jj;
        Int4    BM,BD;
        Int4    gg,gf;
        Int4    bits;
};

#endif

