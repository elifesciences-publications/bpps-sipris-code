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

#include "sph_typ.h"

sph_typ  MakeSapHeap(Int4 hpsz)
{
	sph_typ	H;

	NEW(H,1,sap_heap_type);
	NEW(H->SAP,hpsz+2,sap_typ);
	H->mH=Mheap(hpsz, 3);
	return H;
}

void    NilSapHeap(sph_typ H)
{
	NilMheap(H->mH); free(H->SAP);
	free(H);
}

Int4	InsertSapHeap(sap_typ sap, double key, sph_typ H)
{
	Int4	item;

	item=InsertMheap(key, H->mH);
	if(item==NULL) return NULL;	// failed to insert...?
	H->SAP[item]=sap; 
	return item;
}

sap_typ DelMinSapHeap(double *key, Int4 *Item, sph_typ H)
{
	Int4	item;
	sap_typ sap;

	*key = MinKeyMheap(H->mH);
	item = DelMinMheap(H->mH); *Item = item;
	if(item==NULL) return NULL;
	sap = H->SAP[item]; H->SAP[item]=NULL;
	return sap;
}

sap_typ DelMaxSapHeap(double *key, Int4 *Item, sph_typ H)
{
	Int4	item;
	sap_typ sap;

	*key = MaxKeyMheap(H->mH);
	item = DelMaxMheap(H->mH); *Item = item;
	if(item==NULL) return NULL;
	sap = H->SAP[item]; H->SAP[item]=NULL;
	return sap;
}

