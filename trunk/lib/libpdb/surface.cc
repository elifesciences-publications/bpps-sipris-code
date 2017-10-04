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

#include "surface.h"
#include "exposed.h"

sfc_typ	MkLeeRichards(Int4 maxres, char *filename, a_type A)
{
	sfc_typ	S;
	double	exp,percent;
	Int4	s,i,*res,n;
	char	str[200],aa[10],residue=0;
	keytyp	k;
	FILE	*fptr;

	fptr = open_file(filename,".exp","r");
	NEW(S,1,surface_type);
	NEW(S->exposed, maxres+3,double);
	NEW(S->seq, maxres+3,unsigned char);
	NEW(S->exp, maxres+3,char);
	NEW(S->percent, maxres+3,double);
	strncpy(S->filename,filename,45);
	S->A = A;
	S->len = maxres;
	while(fgets(str,195,fptr) != NULL){
	  if(str[0] == '#') residue = residue_surface(A, str);
	  else if(residue > 0){
	      if(sscanf(str," %d %lf", &s,&exp)!=2) {
		surface_error("Lee & Richards file (*.exp) input error");
 	      }
	      if(s > maxres){
		fprintf(stderr,"site: %d; exp: %f; length: %d\n",
			s, exp, maxres);
	      } else if(s > 0) {
		S->exposed[s] = exp;
		S->percent[s] = 100*exp/maximum_exposed[residue];
		if(S->percent[s] < percent_exposed[1]) S->exp[s] = 1;
		else if(S->percent[s] < percent_exposed[2]) S->exp[s] = 2;
		else S->exp[s] = 3;
		S->seq[s] = residue;
	      }
	  }
	}
	fclose(fptr);
	/** PutLeeRich(stdout, S); /*****/
	return S;
}

double  NegScoreDiffLeeRich(Int4 site, Int4 New, Int4 old, e_type E, void *X)
{ return -ScoreDiffLeeRich(site, New, old, E, X); }

double  NegTotalScoreLeeRich(e_type E, sfc_typ S)
{ return -TotalScoreLeeRich(E, S); }

double	TotalScoreLeeRich(e_type E, sfc_typ S)
{
	Int4	s,v,r;
	double	score;
	unsigned char 	*seq = SeqPtr(E);

	if(LenSeq(E) != S->len) surface_error("TotalScoreLeeRich( ) error");
	for(score = 0, s=1; s <= (Int4) S->len; s++) {
	   if(S->seq[s] > 0){
		v = S->exp[s];
		r = seq[s];
		score += exposed_info[r][v];
	   }
	}
	return score;
}

double	ScoreDiffLeeRich(Int4 site, Int4 New, Int4 old, e_type E, void *X)
{
	register Int4   v;
	register sfc_typ S= (sfc_typ) X;

	if(site < 1 || site > S->len)
                surface_error("ScoreDiffLeeRich( ) input error");
	v = S->exp[site];
	return (exposed_info[New][v] - exposed_info[old][v]);
}

void	PutLeeRich(FILE *fptr, sfc_typ S)
{
	Int4	s;
	
	fprintf(fptr," File %s\n",S->filename);
	fprintf(fptr," Lee & Richards solvent accessible area\n");
	fprintf(fptr,"pos  res loc   area    percent\n");
	for(s=1; s<= S->len; s++) {
		fprintf(fptr,"%4d  %c  (%c) %7.3f  %7.3f\n", s, 
			AlphaChar(S->seq[s],S->A),
			exposed_character[S->exp[s]],
			S->exposed[s],
			S->percent[s]);
	}
	fprintf(fptr,"\n");
}

void	NilLeeRich(sfc_typ S)
{
	free(S->seq); free(S->exp);
	free(S->exposed); free(S->percent);
	free(S);
}

unsigned char	residue_surface(a_type A, char str[])
{
	char	aa[10];

	if(sscanf(str,"# %s", aa)!=1) return 0;
	if(strncmp(aa,"ALA",3)==0) return AlphaCode('A',A);
	if(strncmp(aa,"CYS",3)==0) return AlphaCode('C',A);
	if(strncmp(aa,"ASP",3)==0) return AlphaCode('D',A);
	if(strncmp(aa,"GLU",3)==0) return AlphaCode('E',A);
	if(strncmp(aa,"PHE",3)==0) return AlphaCode('F',A);
	if(strncmp(aa,"GLY",3)==0) return AlphaCode('G',A);
	if(strncmp(aa,"HIS",3)==0) return AlphaCode('H',A);
	if(strncmp(aa,"ILE",3)==0) return AlphaCode('I',A);
	if(strncmp(aa,"LYS",3)==0) return AlphaCode('K',A);
	if(strncmp(aa,"LEU",3)==0) return AlphaCode('L',A);
	if(strncmp(aa,"MET",3)==0) return AlphaCode('M',A);
	if(strncmp(aa,"ASN",3)==0) return AlphaCode('N',A);
	if(strncmp(aa,"PRO",3)==0) return AlphaCode('P',A);
	if(strncmp(aa,"GLN",3)==0) return AlphaCode('Q',A);
	if(strncmp(aa,"ARG",3)==0) return AlphaCode('R',A);
	if(strncmp(aa,"SER",3)==0) return AlphaCode('S',A);
	// if(strncmp(aa,"SEP",3)==0) return AlphaCode('S',A);
	if(strncmp(aa,"THR",3)==0) return AlphaCode('T',A);
	if(strncmp(aa,"VAL",3)==0) return AlphaCode('V',A);
	if(strncmp(aa,"TRP",3)==0) return AlphaCode('W',A);
	if(strncmp(aa,"TYR",3)==0) return AlphaCode('Y',A);
	return 0;
}

Int4	*SurfaceLeeRich(sfc_typ S)
{ return LeeRichExposed(TRUE, S);}

Int4	*BuriedLeeRich(sfc_typ S)
{ return LeeRichExposed(FALSE,S);}

Int4    bubble_sort_surface(Int4 *L)
{
        Int4 i,j,t,n;

	L++;
        for(n=i=0; L[i]!= -1; i++) n++;
        for(i=1; i < n; i++){
                for(j=i;j > 0 && (L[j-1] > L[j]);j--){
                        t=L[j]; L[j]=L[j-1]; L[j-1]=t;
                }
        }
	return n;
}

Int4	*LeeRichExposed(BooLean surface, sfc_typ S)
/** sample input format: ' PHE    18      76.029'  **/
/** use heap to save the half most exposed  residues **/
{
	Int4	s,i,*res,n,maxres;
	keytyp	k;
	dh_type H;

	H = dheap(S->len+1,4);
	for(s=1; s <= (Int4) S->len; s++) {
	   if(S->seq[s] > 0){
 		k = (keytyp)S->exposed[s];
		if(surface) k = -k;
		insrtHeap(s,k,H);
	   }
	}
	maxres = ItemsInHeap(H);
	n = 1 + (maxres/2); /*** half ***/
	/** n = 1 + (2*maxres/3); /*** 2/3rds ***/
	/** n = 1 + (4*maxres/5); /*** 4/5ths ***/
	NEW(res, maxres+3,Int4);
	for(i=1; (s = delminHeap(H)) != NULL; i++){ 
		res[i] = s; 
		if(i >= n) break; 
	}
	res[i] = -1;
	Nildheap(H); 
	bubble_sort_surface(res);
	return res;
}

Int4	CntExposedLeeRich(Int4 *surface, Int4 *buried, sfc_typ S)
/*******************************************************************
 Count the number of each residue type that are found at the surface 
 versus in the protein interior.
 (uses heap to save the half most exposed  residues)
 ********************************************************************/

{
	Int4	s,i,n,maxres = S->len;
	keytyp	key;
	char	r;
	dh_type H;

	H = dheap(maxres+1,4);
	for(s=1; s <= (Int4) S->len; s++) {
	   if(S->seq[s] > 0){
		/**** Option A: percent surface exposed ****
		key = (keytyp)S->percent[s]
		/**** Option B: actual surface exposed *****/
		key = (keytyp)S->exposed[s];
		/*******************************************/
		insrtHeap(s,key,H); 
	   }
	}
	n = ItemsInHeap(H);
	if(n%2==0){	n = n/2; }
	else {	n = 1 + n/2; }
	for(i=1; (s = delminHeap(H)) != NULL; i++){ 
		buried[i] = s; 
		if(i >= n) break; 
	} buried[i] = -1;
	for(i=1; (s = delminHeap(H)) != NULL; i++){ surface[i] = s; }
	surface[i] = -1;
	Nildheap(H); 
	return n;
}

void    surface_error(char *s){fprintf(stderr,"surface: %s\n",s);exit(1);}

