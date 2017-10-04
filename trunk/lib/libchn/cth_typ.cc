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

#include "cth_typ.h"

Int4	SsetToStr(sst_typ sst, char *str,a_type AB)
{
	Int4 i,j;
        for(j=0,i=1;i<=nAlpha(AB);i++) {
           if(MemSset(i,sst)){ str[j]=AlphaChar(i,AB); j++; }
        } str[j]=0;
	return j;
}

//************************ cti_typ.h ***************************

cti_typ::cti_typ(Int4 i, Int4 j,sst_typ sst[2],set_typ Set[2], unsigned char id_ij[2],
		unsigned char ri, unsigned char rj, t_type tab,
		Int4 node, Int4 **Cell, Int4 Factor)
{
	cell[0]=0;
	NEW(cell[1],4,Int4); NEW(cell[2],4,Int4);
        ResI = ri; ResJ=rj;
   if(Cell){
	cell[1][1]=Cell[1][1]; cell[1][2]=Cell[1][2];
	cell[2][1]=Cell[2][1]; cell[2][2]=Cell[2][2];
   }
	expect=FALSE;
        Node=node;
        SiteI=i;SiteJ=j;
        factor=4;	// multiplication factor for table.
        setI=Set[0]; 
        setJ=Set[1];
	sstI=sst[0];
	sstJ=sst[1];
	rsidI=id_ij[0];
	rsidJ=id_ij[1];
	Tab=tab;
	factor=Factor;
	// PutTableMod(stderr,Tab);
}

double  cti_typ::FisherExactTest(register Int4 ro, register Int4 bo,
        register Int4 ri, register Int4 bi)
/*************** Fisher's Exact Test for a 2x2 contingency table *****
 *  observed:
 *
 *              out     in                        red    black
 *      -----+-------+-------+-------      -----+-------+-------+-------
 *      red  |   ro  |   ri  | ro+ri        out |   ro  |   bo  | ro+bo
 *      -----+-------+-------+-------  or  -----+-------+-------+-------
 *      black|   bo  |   bi  | bo+bi        in  |   ri  |   bi  | ri+bi
 *      -----+-------+-------+-------      -----+-------+-------+-------
 *           | ro+bo | ri+bi |                  | ro+ri | bo+bi |
 *
 *   Randomly draw out ro+bo balls from an urn containing ro+ri red balls
 *   and bo+bi black balls.  The FisherExactTest() returns the chance of
 *   obtaining >= ro red balls among the drawn balls.
 **********************************************************************/
// WARNING: about 50x's slower than cntab1() routine below!
{
#if 0
        return CumHypGeomProb(ro+ri,bo+bi,ro+bo,ro);
#else
        register double p;
        register Int4   end,C;

        if(ro == 0) return 1.0;
        end = MINIMUM(Int4,ro+ri,ro+bo); // minimum of red or out
        C = (lnfact(ro+ri)+lnfact(bo+bi)+lnfact(ro+bo)+lnfact(ri+bi));
        C -= lnfact(ro+ri+bo+bi);
        for(p=0.0; ro <= end; ro++,bi++, ri--,bo--){
           p += exp(C-(lnfact(ro)+lnfact(ri)+lnfact(bo)+lnfact(bi)));
        } return p;
#endif
}

void	cti_typ::Free()
{
	if(setI) NilSet(setI);
	if(setJ) NilSet(setJ);
	NilTab(Tab);
	free(cell[1]); free(cell[2]);
}


//************************ cth_typ.h ***************************

cth_typ::cth_typ(Int4 HeapSize)
{
	hpsz=HeapSize;
	assert(hpsz > 0);
	NEWP(table,hpsz+3,cti_typ);
	mH = Mheap(hpsz,3);
}

void	cth_typ::Free()
{
	cti_typ *cti;
	while(cti=this->DeleteMin( )){ delete cti; }
	// for(Int4 i=1; i<= hpsz; i++) if(table[i]) delete table[i];
	free(table);
	NilMheap(mH);
}

BooLean cth_typ::Insert(double key,cti_typ *cti)
{
	Int4	item = InsertMheap(key,mH);
	if(item > 0){
		if(table[item]) delete table[item];
		table[item] = cti;
		return	TRUE;
	} else return FALSE;
}

cti_typ *cth_typ::DeleteMin( )
{
	cti_typ *cti;
	if(EmptyMheap(mH)) return 0;
	Int4 item = DelMinMheap(mH);
	if(item){ cti = table[item]; table[item]=0; return cti; }
	else return 0;
}


