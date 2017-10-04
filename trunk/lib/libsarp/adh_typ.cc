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

#include "adh_typ.h"

//*********************************** adv_typ ***********************************

void    adv_typ::Init(Int4 i, Int4 j, double v, double m, Int4 N)
{ variance = v; mean=m; Num=N; I=i; J=j; }

void    adv_typ::Free( )
{
	// do nothing...
}

void    adv_typ::Put(FILE *fp)
{
	fprintf(fp,"columns %d vs %d: var = %.2f; mean = %.2f; Num = %d\n",
		I,J,variance,mean,Num);
}

//*********************************** adh_typ ***********************************

void    adh_typ::Init(Int4 hs)
{
	hpsz=hs;
	mheap=Mheap(hpsz, 3);
	NEWP(ADV,hpsz +5, adv_typ);
}

void    adh_typ::Free()
{
	double	key;
	Int4	item;
	adv_typ *adv;

	while((adv=DelMin(&key,&item)) != NULL){
		delete adv;
	}
	free(ADV);
	NilMheap(mheap); 
}

Int4	adh_typ::Insert(adv_typ *adv, double key)
{
	Int4	item;

	item=InsertMheap(key, mheap);
	if(item==0) return 0;
	if(ADV[item]!= NULL) delete ADV[item];
	ADV[item]=adv; 
	return item;
}

adv_typ	*adh_typ::DelMin(double *key, Int4 *Item)
{
	Int4	item;
	adv_typ *adv;

	*key = MinKeyMheap(mheap);
	item = DelMinMheap(mheap); *Item = item;
	if(item==0) return NULL;
	adv = ADV[item]; ADV[item]=0;
	return adv;
}

adv_typ	*adh_typ::DelMax(double *key, Int4 *Item)
{
	Int4	item;
	adv_typ *adv;

	*key = MaxKeyMheap(mheap);
	item = DelMaxMheap(mheap); *Item = item;
	if(item==0) return NULL;
	adv = ADV[item]; ADV[item]=NULL;
	return adv;
}

