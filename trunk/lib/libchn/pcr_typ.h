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

#if !defined(_PCR_TYP_H_)
#define _PCR_TYP_H_

#include "stdinc.h"
#include "alphabet.h"

#define	MAX_FUNCTIONAL_CATEGORIES 26

class pcr_typ {	// printable colored (categorized) residue type
public:
	pcr_typ( ){
		for(char c='A'; c <= 'Z'; c++){ head[c]=tail[c]=NULL; }
	}
	~pcr_typ( ){ Free( ); }
	void print(FILE *fp){ print(fp,0); }
	void print(FILE *fp,Int4 offset){
	   for(char c='A'; c <= 'Z'; c++){
		pcr_cell *tmp=head[c];
		while(tmp != NULL){
		     fprintf(fp,"%c%d.%c",tmp->res,tmp->site+offset,c);
		     if(tmp != tail[c]){ fprintf(fp,","); }
		     else fprintf(fp,"\n\n");
		     tmp=tmp->next;
		}
	   }
	}
	void append(char r, Int4 s, char c, Int4 scr){
		assert(c >= 'A' && c <= 'Z');
		assert(r >= 'A' && r <= 'Z');
		pcr_cell *tmp = new pcr_cell;
		tmp->next=NULL;
		tmp->res = r; tmp->site = s; tmp->score=scr;
		if(head[c] == NULL){
		    head[c] = tail[c] = tmp;
		} else {
		   assert(tail[c] != NULL);
		   tail[c]->next = tmp; tail[c] = tmp;
		}
	}
private:
	struct pcr_cell // list cell
	{
		char		res;
		Int4		site;
		Int4		score;
		pcr_cell	*next;
	};
	pcr_cell	*head[128];
	pcr_cell	*tail[128];
	void	Free( ){
	   for(char c='A'; c <= 'Z'; c++){
	        pcr_cell *tmp=head[c];
		while(tmp != NULL){
	           pcr_cell *nxt=tmp->next;
		   delete tmp;
		   tmp=nxt;
		}
	   }
	}
};

#endif

