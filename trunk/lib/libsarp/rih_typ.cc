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

#include "rih_typ.h"

//*********************************** rpi_typ ***********************************

void    rpi_typ::Init(Int4 i, Int4 j, double scr, Int4 n,char *pI, char *pJ,Int4 ca, char *Iname, char *Jname)
{
	Score=scr; Num=n; I=i; J=j; pttrnI=pI; pttrnJ=pJ; CA=ca; nameI=Iname; nameJ= Jname; Lineage=' ';
}

BooLean	rpi_typ::TheSame(rpi_typ *rpi )
{
	if(this->I != rpi->I) return FALSE;
	if(this->J != rpi->J) return FALSE;
	if(strcmp(this->pttrnI,rpi->pttrnI) != 0) return FALSE;
	if(strcmp(this->pttrnJ,rpi->pttrnJ) != 0) return FALSE;
	if(this->Num != rpi->Num) return FALSE;
	return TRUE;
}

void    rpi_typ::Free( )
{
	// do nothing...
}

void    rpi_typ::Put(FILE *fp)
{
	if(pttrnI && pttrnJ) fprintf(fp," %s%d vs %s%d: %.2f contact score (%d)", pttrnI,I,pttrnJ,J,Score,Num); 
	else fprintf(fp," %d vs %d: %.2f contacts (%d)", I,J,Score,Num); 
	if(nameI && nameJ){
		if(nameI == nameJ) fprintf(fp,"\t--> %s %c\n",nameI,Lineage); 
		else fprintf(fp,"\t--> %s + %s %c\n",nameI,nameJ,Lineage); 
	} else fprintf(fp," %c\n",Lineage);
}

//*********************************** rih_typ ***********************************

void    rih_typ::Init(Int4 hs)
{
	hpsz=hs;
	mheap=Mheap(hpsz, 3);
	NEWP(RPI,hpsz +5, rpi_typ);
}

void    rih_typ::Free()
{
	double	key;
	Int4	item;
	rpi_typ *rpi;

	while((rpi=DelMin(&key,&item)) != NULL){ delete rpi; }
	free(RPI); NilMheap(mheap); 
}

Int4	rih_typ::Insert(rpi_typ *rpi, double key)
{
	Int4	item,i;

#if 1	// If already in heap then don't add it!
	for(i=1; i <= hpsz; i++){
	   if(RPI[i] != 0){
		if(rpi->TheSame(RPI[i])) return 0; 
	   } 
	}
#endif
	item=InsertMheap(key, mheap);
	if(item==0) return 0;
	if(RPI[item]!= NULL) delete RPI[item];
	RPI[item]=rpi; 
	return item;
}

rpi_typ	*rih_typ::DelMin(double *key, Int4 *Item)
{
	Int4	item;
	rpi_typ *rpi;

	*key = MinKeyMheap(mheap);
	item = DelMinMheap(mheap); *Item = item;
	if(item==0) return NULL;
	rpi = RPI[item]; RPI[item]=0;
	return rpi;
}

rpi_typ	*rih_typ::DelMax(double *key, Int4 *Item)
{
	Int4	item;
	rpi_typ *rpi;

	*key = MaxKeyMheap(mheap);
	item = DelMaxMheap(mheap); *Item = item;
	if(item==0) return NULL;
	rpi = RPI[item]; RPI[item]=NULL;
	return rpi;
}

