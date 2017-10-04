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

/* rih_typ.h - . */
#if !defined (RIH_TYP)
#define RIH_TYP
#include "cmsa.h"
#include "atom.h"
#include "sequence.h"
#include "histogram.h"
#include "residues.h"
#include "mheap.h"

class rpi_typ {		// residue pairwise interation (rpi) type
public:
		rpi_typ(){ print_error("not allowed"); }
#if 0
		rpi_typ(Int4 i, Int4 j, double scr, Int4 n){ Init(i,j,scr,n,0,0,0); }
		rpi_typ(Int4 i, Int4 j, double scr, Int4 n,char *pI,char *pJ){ Init(i,j,scr,n,pI,pJ,0); }
		rpi_typ(Int4 i, Int4 j, double scr, Int4 n,char *pI,char *pJ,Int4 ca){ Init(i,j,scr,n,pI,pJ,ca); }
#else
		rpi_typ(Int4 i, Int4 j, double scr, Int4 n,char *pI=0,char *pJ=0,Int4 ca=0)
			{ Init(i,j,scr,n,pI,pJ,ca); }
		rpi_typ(char *Iname,char *Janame,Int4 i,Int4 j,double scr,Int4 n,char *pI=0,char *pJ=0,Int4 ca=0,
			char line=' ')
			{ Init(i,j,scr,n,pI,pJ,ca,Iname,Janame); Lineage=line; }
#endif
		~rpi_typ(){ Free(); }
	char	Lineage;
	Int4	ColumnI(){ return I; }
	Int4	ColumnJ(){ return J; }
	Int4	ContrastAlignment(){ return CA; }
	Int4	Number(){ return Num; }
	Int4	SetNumber(Int4 n){ assert(n >= 0); Num=n; return Num; }
	BooLean	TheSame(rpi_typ *rpi);
	void	Put(FILE *fp);
private:
	void	Init(Int4 i,Int4 j,double scr,Int4 n,char *pI,char *pJ,Int4 ca,char *Iname=0,char *Jname=0);
	void	Free();
	double	Score;
	Int4	Num,I,J,CA;
	char	*pttrnI,*pttrnJ;
	char	*nameI,*nameJ;
};

class rih_typ {		// rih_typ = atomic distance heap type
public:
		rih_typ(){ print_error("not allowed"); }
		rih_typ(Int4 h){ Init(h); }
		~rih_typ(){ Free(); }
	Int4    Insert(rpi_typ *rpi,double key);
	rpi_typ	*DelMin(double *key, Int4 *Item);
	rpi_typ	*DelMax(double *key, Int4 *Item);
private:
	mh_type mheap;
	Int4	hpsz;
	rpi_typ	**RPI;
	void	Init(Int4 hs);
	void	Free();
};

#endif

