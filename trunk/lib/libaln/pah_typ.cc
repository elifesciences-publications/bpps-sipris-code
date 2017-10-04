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

#include "pah_typ.h"

pah_typ::pah_typ(Int4 Size)
{
	NEW(info, Size+3, psaln_info_type);
	hpsz = Size;
	mH=Mheap(hpsz, 3);
	HG = Histogram("Scores",-100000,100000,200.0);
}

void	pah_typ::Free()
{
	for(Int4 i=1; i <= hpsz; i++) NilInfo(i);
	free(info); NilMheap(mH); NilHist(HG);
}

void	pah_typ::MkInfo(Int4 item,e_type E, char *operation,
	Int4 Score, unsigned char rpts, Int4 start)
{
	assert(item > 0 && item <= hpsz);
	NilInfo(item);
	info[item].E = E;
	info[item].operation = operation;
	info[item].Score = Score;
	info[item].start= start;
	info[item].rpts = rpts;
}

void	pah_typ::NilInfo(Int4 item)
{
	if(info[item].E) NilSeq(info[item].E);
	if(info[item].operation) free(info[item].operation);
}

Int4    pah_typ::Insert(e_type E, char *operation, Int4 Score, Int4 start)
{ return Insert(E, operation,Score,1,start); }

Int4    pah_typ::Insert(e_type E, char *operation, Int4 Score, unsigned char rpts,
	Int4 start)
{
        Int4    item,i;

        IncdHist(Score,HG);
        item=InsertMheap(-Score,mH);
        if(item==NULL) return NULL;
        NilInfo(item);
	MkInfo(item,E,operation,Score,rpts,start);
        return item;
}

Int4	pah_typ::DelMin(e_type *E, char **operation, unsigned char *rpts, Int4 *start,
	Int4 *score)
{
	Int4    item;

        item = DelMaxMheap(mH);
        if(item==NULL) return NULL;
        *E = info[item].E; info[item].E=0;
        *operation = info[item].operation; info[item].operation=0;
	*rpts = info[item].rpts; 
	*score = info[item].Score;
	*start=info[item].start;
        return item;
}

Int4	pah_typ::DelMax(e_type *E, char **operation, unsigned char *rpts,Int4 *start,
	Int4 *score)
{
	Int4    item;

        item = DelMinMheap(mH);
        if(item==NULL) return NULL;
        *E = info[item].E; info[item].E=0;
        *operation = info[item].operation; info[item].operation=0;
	*rpts = info[item].rpts; 
	*score = info[item].Score;
	*start=info[item].start;
        return item;
}

void    pah_typ::Put(FILE *fp)
{
	PutHist(fp,60,HG);
}


