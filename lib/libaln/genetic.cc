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

#include "genetic.h"

ga_typ  MakeGenetic(Int4 size, Int4 aveblk, Int4 avecol, Int4 minlenseq)
{
	Int4	Nc,b;
	ga_typ	GA;

	NEW(GA,1,genetic_type);
	GA->size=size;
        GA->maH = MkMSAHeap(size);

	Nc = minlenseq;
	/**** NEW ****/
	if(Nc < 10) print_error("MakeGenetic( ): seq.too short");
	if(aveblk*avecol > Nc){ /** adjust for short sequences **/
	   for(b = aveblk; TRUE; b--){
		avecol = (Int4) floor((double)Nc/((double)b*2.0));
		if(avecol >= 4) { aveblk = b; break; }
	   } 
	   fprintf(stderr,"aveblk = %d; avecol =%d\n",aveblk,avecol);
	   if(aveblk < 2) print_error("MakeGenetic( ): seq.too short");
	}
	/**** OLD ****
	if(Nc<35 || avecol>Nc) print_error("MakeGenetic( ): seq.too short");
	/**** OLD ****/
	if(aveblk > 30) print_error("MakeGenetic( ) input error");
	NEW(GA->ncols,Nc+3,Int4);
	GA->totcols = aveblk*avecol;
	GA->Bc = MakeBinomial(Nc, avecol,"number of columns");
	GA->Bb = MakeBinomial(30, aveblk,"number of blocks");
	PutBinomial(stderr, 10000, GA->Bb);
	PutBinomial(stderr, 10000, GA->Bc);
	return GA;
}

cma_typ	NilGenetic(ga_typ GA)
{
	cma_typ ma;
	double	map;

	ma = DelMinMSAHeap(&map, GA->maH);
	free(GA->ncols);
	NilMSAHeap(GA->maH);
	NilBinomial(GA->Bc); NilBinomial(GA->Bb);
	free(GA);
	return ma;
}

#if 0
double	BreedGenetic(Int4 narg, char *argument[], double bestmap, 
	char *name, ga_typ GA, gd_type Gd)
{
   gs_type	G;
   double	map,map1,map2,*mapn,*mapo,temperature;
   cma_typ	ma,ma1,ma2,*New,*old;
   Int4		N,n,nn,no,i,j,k,o,c=1,x;
   mah_typ	H = HeapGenetic(GA);
   BooLean	*recomb,flag;

   N = 2*SizeMSAHeap(H); 
   NEW(New,N+2,cma_typ); NEW(old,N+2,cma_typ);
   NEW(mapn,N+2,double); NEW(mapo,N+2,double); 
   NEW(recomb,N+2,BooLean);
   PurgeMSAHeap(H); // remove all but one of those items with identical keys
   for(n=0,c=1; !EmptyMSAHeap(H); c++){
	for(no=0,o=1; o <= n; o++){  // n==0 --> skipped the first time around 
	   if(recomb[o]){ NilCMSA(New[o]); } // discard parent alignments 
	   else { ++no; old[no] = New[o]; mapo[no]=mapn[o]; }
	}
        for(nn=0; (ma=DelMinMSAHeap(&map, H)) != NULL; )
		{ ++nn; New[nn]=ma; mapn[nn] = map; recomb[nn]=FALSE; }
	for(n=nn,o=1; o <= no; o++) // also add previous recombinants to heap
		{ ++n; New[n]=old[o]; mapn[n]=mapo[o]; recomb[n]=FALSE; }
	fprintf(stderr,
	   "******** %d: %d new, %d old, %d total alignments ********\n",
		c,nn,no,n);
	if(no==0) nn--;
	for(i=1; i <= nn; i++){  /** try every combination of recombinants **/
	   ma1 = New[i]; map1=mapn[i];
	   for(j=i+1; j <= n; j++){
	     ma2 = New[j]; map2=mapn[j];
	     fprintf(stderr,"x");
	     if((ma=RecombineCMSA(ma1,ma2)) != NULL){
		   // if get a recombinant then try to refine it. 
		   recomb[i]=recomb[j]=TRUE;
		   fprintf(stderr,"RecombineCMSA( ) map = %g\n",RelMapCMSA(ma));
		   map=SimulatedAnnealingGibbs(narg,argument,&ma,'S',150);
		   if(map > bestmap){
			bestmap=map; PutAlnCMSA(name,ma,Gd); SetGenetic(ma, GA);
			fprintf(stderr,"bestmap = %g\n",bestmap);
		   } 
		   if(!KeyInMSAHeap(map, H)){
			if(InsertMSAHeap(ma, map, H) == NULL) NilCMSA(ma);
			else { fprintf(stderr,"\thybrid map=%f\n",map); }
		   } else NilCMSA(ma);
	     } 
#if 0 // NEED TO ADD "SAME SPECIES" TEST FOR CMSAs.
             else if(nBlksCMSA(ma1) == nBlksCMSA(ma2)){
                fprintf(stderr,"TEST INTERSECTION PROCEDURE (WARNING THIS IS FULL OF BUGS!!!!)\n");
                ma = IntersectionCMSA(ma1,ma2);
                if(ma != NULL){
                  map=SimulatedAnnealingGibbs(narg,argument,&ma,'S',250);
                  if(map > bestmap){
                        bestmap=map; PutAlnCMSA(name,ma,Gd); SetGenetic(ma, GA);
                        fprintf(stderr,"bestmap = %g\n",bestmap);
                        fprintf(stderr,"IntersectionCMSA( ) worked!\n");
                        fprintf(stderr,"...but this is still not really implemented!!!!\n");
                  }
                  if(!KeyInMSAHeap(map, H)){
                        if(InsertMSAHeap(ma, map, H) == NULL) NilCMSA(ma);
                        else { fprintf(stderr,"\thybrid map=%f\n",map); }
                  } else NilCMSA(ma);
                }
             }
#endif
   	   }
	}
	fprintf(stderr,"\n");
	if(ConvergedMSAHeap(H) < 2.0) break; // I haven't tested this at all!!!!
	// ConvergedMSAHeap(H);
   }
   for(o=1; o <= n; o++) {
	/** after breeding move the parent population to the heap **/
	if(InsertMSAHeap(New[o],mapn[o],H)==NULL) NilCMSA(New[o]);
   }
   PurgeMSAHeap(H);
   // get the map for the best alignment.
   ma=DelMinMSAHeap(&map, H); 
   if(InsertMSAHeap(ma,map,H)==NULL) NilCMSA(ma);
   // get the map for the best alignment.
   free(New); free(old); free(mapo); free(mapn); free(recomb);
   return map;
}
#endif

Int4	BestBlksGenetic(ga_typ GA) { return AveBinomial(GA->Bb); }

Int4	BestColsGenetic(Int4 t, ga_typ GA)
{
	Int4	b=AveBinomial(GA->Bb);

	if(t < 1 || t > b) print_error("BestColsGenetic( ) input error");
	return GA->ncols[t];
}

void	SetGenetic(cma_typ msa, ga_typ GA)
{
	Int4	i,b,c,*cols=GA->ncols;

	b=ConfigCMSA2(cols, msa); /**** NEW ****/
	/** b=ConfigCMSA(cols, msa); /**** OLD ****/
	for(c=0,i=1; i<=b; i++) { c+=cols[i]; }
	GA->totcols=c;
	c = (Int4) floor( ((double)c/(double)b) + 0.5);
	SetBinomial(c,GA->Bc);
	SetBinomial(b,GA->Bb);
	fprintf(stderr,"\t%d blocks; %d columns\n",b,GA->totcols);
	PutBinomial(stderr, 10000, GA->Bb);
	PutBinomial(stderr, 10000, GA->Bc);
}

Int4	GeneticSampC2(Int4 t, Int4 m,Int4 max,ga_typ GA)
{ 
	Int4	c;
	
	if(t < 1 || t > AveBinomial(GA->Bb)) 
		print_error("GeneticSampC2( ) input error");
	c = GA->ncols[t];
	SetBinomial(c,GA->Bc);
	if(c < m || c > max)
		print_error("GeneticSampC2( ) min/max input error");
	return SampleBinomial(m,max,GA->Bc);
}

Int4	GeneticSampC(Int4 b, Int4 m,Int4 max,ga_typ GA)
{ 
	Int4	c;
	
	c = GA->totcols;
	c = (Int4) ceil( ((double)c/(double)b) );
	SetBinomial(c,GA->Bc);
	return SampleBinomial(m,max,GA->Bc);
}

Int4	GeneticSampB(Int4 m,Int4 max,ga_typ GA)
{ return SampleBinomial(m,max,GA->Bb); }

