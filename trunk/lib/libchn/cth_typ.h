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

#if !defined (_CTH_TYP_)
#define	_CTH_TYP_

#include "alphabet.h"
#include "table.h"
#include "dsets.h"
#include "set_typ.h"
#include "probability.h"
#include "mheap.h"

Int4    SsetToStr(sst_typ sst, char *str,a_type AB);

class cti_typ { 	// contingency table item
public:
		cti_typ(Int4 i, Int4 j,sst_typ sst[2],set_typ Set[2], 
			unsigned char id_ij[2],
			unsigned char ri, unsigned char rj, t_type T,
			Int4 node, Int4 **Cell, Int4 Factor);
		cti_typ( ){ assert(!"Illegal constructor"); }
		~cti_typ( ){ Free(); }
	t_type	Table( ) { return Tab; }
	double	Obs00( ) { return valTab(0,0,Tab); }
	double	Exp00( ) { if(!expect) ExpectTable(E,Tab); return E[0][0]; }
	double	ExactTest() { return ExactTab(Tab); }
	set_typ	SetI() { return setI; }
	set_typ	SetJ() { return setJ; }
	unsigned char	ResI,ResJ;
	unsigned char rsidI,rsidJ;	// residue set identifiers
	Int4	SiteI,SiteJ;
	void	PutTable(FILE *fp){ PutTableMod(fp,Tab); }
	void	PutCells(FILE *fp){
                    fprintf(fp," %5d | %5d\n",cell[1][1],cell[1][2]);
                    fprintf(fp," %5d | %5d",cell[2][1],cell[2][2]);
                    fprintf(fp," (total: %d)\n",cell[1][1]+cell[1][2]+cell[2][1]+cell[2][2]);
                    fprintf(fp," Fisher's Exact = %g\n",
                        FisherExactTest(cell[1][1],cell[1][2],cell[2][1],cell[2][2]));
		}
	void	StrResSetI(char *str,a_type AB){ 
#if 1
		SsetToStr(sstI, str,AB);
#else 
		  Int4 i,j;
		  for(j=0,i=1;i<=nAlpha(AB);i++) {
           		if(MemSset(i,sstI)){ str[j]=AlphaChar(i,AB); j++; } 
		  } str[j]=0;
#endif
                }
	void	StrResSetJ(char *str,a_type AB){ 
#if 1
		SsetToStr(sstJ, str,AB);
#else 
		  Int4 j,i;
		  for(j=0,i=1;i<=nAlpha(AB);i++) {
           		if(MemSset(i,sstJ)){ str[j]=AlphaChar(i,AB); j++; } 
		  } str[j]=0;
#endif
                }
	double  FisherExactTest(register Int4 ro, register Int4 bo,
        			register Int4 ri, register Int4 bi);
	sst_typ	sstI,sstJ;
private:
	void	Free();
	BooLean	expect;
	double	E[2][2];
	t_type	Tab;
	double	FisherPval;
	double	ChiSq;
	Int4	*cell[3];
	Int4	Node;
	Int4	factor;	// multiplication factor for table.
	set_typ	setI,setJ;	// intersection set.
};

#if 0
	Node:  Node = (set_num * N) + Site); // range 1...(set_num +1) *N
   Note that:
	(Node - 1) % N + 1 == Site.
	Node / N == subset.
#endif

class cth_typ { 	// contingency table heap
public:
		cth_typ( ){ assert(!"Illegal constructor"); }
		cth_typ(Int4 hpsz);
		~cth_typ( ){ Free(); }
	BooLean	Insert(double, cti_typ *);
	cti_typ	*DeleteMin( );
	BooLean	Empty( ){ return EmptyMheap(mH); }
	double	MinKey( ){ return MinKeyMheap(mH); }
	Int4	NumInHeap( ){ return ItemsInMheap(mH); }
private:
	void	Free();
	mh_type	mH;
	Int4	hpsz;
	cti_typ	**table;
};


#endif

