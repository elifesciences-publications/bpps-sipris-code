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

#include "dft_typ.h"

void    dft_typ::Init(Int4 N, Int4 *Array)
{
	num_nodes=N;
	NEW(Level,N+3,Int4);
	for(Int4 i=0; i < N; i++) Level[i]=Array[i];
	ValidityCheck();
}

BooLean	dft_typ::Read(FILE *fp)
// Read a 
// e.g., "0,1,2,2,2,3,3,3,1."  (spaces allowed anywhere)
{
	char	c,str[200];
	Int4	i,j,N=0;
	do {
		c = fgetc(fp);
		if(isspace(c)) continue;
		else if(isdigit(c)){
		   while(isdigit(c)){ c = fgetc(fp);}
		   ungetc(c,fp); N++; 
		} else if(c != ',' && c != '.') return FALSE;
        } while(c != '.'); 
	// syntax is okay if reached this point.
	rewind(fp);
	NEW(Level,N+3,Int4);
	i=0;
	do {
		c = fgetc(fp);
		if(isspace(c)) continue;
		else if(isdigit(c)){
		   j=0; str[j]=c;
		   while(isdigit(c=fgetc(fp))){ j++; str[j]=c; }
		   j++; str[j]=0;
		   if(sscanf(str,"%d",&Level[i])!=1) print_error("dft_typ::Read(): this should not happen");
		   else i++;
		} else if(c != ',' && c != '.') return FALSE;
        } while(c != '.'); 
	if(i != N) print_error("dft_typ::Read(): this should not happen");; 	
	num_nodes=N;
	return ValidityCheck();
}

void	dft_typ::Write(FILE *fp)
{
	for(Int4 i=0; i < num_nodes; i++){
		fprintf(fp,"%d",Level[i]);
		if(i < num_nodes-1) fprintf(fp,",");
	}
	fprintf(fp,".\n");
}

void    dft_typ::Free()
{
	free(Level);
}

BooLean	dft_typ::ValidityCheck()
// L[0]=0. (generates the root node)
// 1 <= L[i] <= L[i-1] +1 for i >= 1.
// Adds an ith node as a child along the lineage leading to the i -1 node.

{
	if(Level[0] != 0) return FALSE;
	for(Int4 i=1; i < num_nodes; i++){
		if(Level[i] < 1 || Level[i] > (Level[i-1] + 1)){
		// if(Level[i] > (Level[i-1] + 1)){
			Valid=FALSE;
			return FALSE;
		} 
	}
	Valid=TRUE;
	return TRUE;
}


