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

#include "guide.h"

gd_type MkGuide(char *infile, a_type A)
{
	gd_type	G;

	NEW(G,1,guide_type);
	G->E = ReadGuideSeqFA(infile, 1, &G->ignore, A);
	G->len = LenSeq(G->E);
	return G;
}

gd_type CopyGuide(gd_type G)
{
	gd_type	G2;
	Int4	s;

	NEW(G2,1,guide_type);
	G2->E = CopySeq(G->E);
	G2->len = G->len;
	NEW(G2->ignore, G2->len +2, BooLean);
	for(s=1; s <= G2->len; s++) G2->ignore[s] = G->ignore[s];
	return G2;
}

void    PutGuide(FILE *fptr,gd_type G, a_type A)
{
	Int4	i,s,len;
	unsigned char	*seq;
	char	c;
	BooLean	low;

	PutSeqInfo2(fptr,G->E);
	seq = SeqPtr(G->E);
	len = LenSeq(G->E);
	for(low=FALSE,i=s=1; s <= len; s++){
	   c = AlphaChar(seq[s],A);
	   if(G->ignore[s]){
		 if(!low) { 
			low = TRUE; i++;
	         	if(i == 60){ fprintf(fptr,"\n"); i=2; }
			fprintf(fptr,"/"); 
		 }
		 fprintf(fptr,"%c", tolower(c));
	   } else {
		 if(low) { 
			i++; 
			fprintf(fptr,"\\"); low = FALSE;
	         	if(i == 60){ fprintf(fptr,"\n"); i=1; } 
		 }
		 fprintf(fptr,"%c",c);
	   }
	   i++; if(i == 60) { fprintf(fptr,"\n"); i=1; }
	}
	fprintf(fptr,"\n\n");
}

void    NilGuide(gd_type G)
{
	NilSeq(G->E);
	free(G->ignore);
	free(G);
}

void    guide_error(char *s){fprintf(stderr,"guide error: %s\n",s);exit(1);}


