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
#include "str_typ.h"

// Initialize all items to be in lists containing only themselves.
str_typ::str_typ(char c) { label = c; next = 0; prev = 0; }

str_typ::~str_typ() {
#if 0
	if(prev != 0) prev->next=next;
	if(next != 0) next->prev=prev;
#else
	if(next != 0) delete next;
#endif
}

str_typ	*str_typ::Append(char label) {
	str_typ *nxt=this->next,*str = new str_typ(label);
	str->next=nxt; this->next = str; str->prev=this; 
	if(nxt) nxt->prev=str;
	return str;
}

char	*str_typ::Return()
{
	str_typ *str;
	Int4	n;
	for(n=0,str = this; str != 0; str=str->next) n++;
	char	*rtn; NEW(rtn,n+9,char);
	for(n=0,str = this; str != 0; str=str->next){ rtn[n]=str->label; n++; }
	return rtn;
}

// Print the contents of list starting from 'this'.
void str_typ::Print(FILE *fp)
{
	for(str_typ *str = this; str != 0; str=str->next) fprintf(fp,"%c",str->label);
	fprintf(fp,"\n");
}

