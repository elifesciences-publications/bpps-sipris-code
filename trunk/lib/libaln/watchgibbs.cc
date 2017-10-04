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

#include "watchgibbs.h"

wg_type MakeWatchGibbs(Int4 maxblk, Int4 maxcol)
{
	wg_type	WG;
	Int4	i,j;

	NEW(WG,1,watchgibbs_type);
	NEWP(WG->map,maxblk +2,float);
	for(i=1; i<=maxblk; i++){
	  NEW(WG->map[i],maxcol+2,float);
	  for(j=1; j<=maxcol; j++){
		WG->map[i][j] = -FLT_MAX;
	  }
	}
	NEW(WG->n,maxcol+2,Int4);
	WG->N = 0;
	WG->maxcol = maxcol;
	WG->lowcol = maxcol;
	WG->hicol = 0;
	WG->maxblk = maxblk;
	return WG;
}

void	PutContourWatchGibbs(FILE *fp, wg_type WG)
// output a contour map of gibbs sampler...
// scale:  0  1  2  3  4  5  6  7  8  9    (9 >= 0.9 * max_map)
// symbol:' ' .  -  +  4  5  6  7  8  9 
{
	

}

void    PutWatchGibbs(FILE *fp, wg_type WG)
{
	Int4	i,j,maxblk,k,k0,s;
	float	inc,max;

   if(WG->N > 0){	
	inc = WG->bestmap/50.0;
	fprintf(fp,"LPR as a function of the number of colums:\n");
	for(j=WG->lowcol; j<=WG->hicol; j++){
	    max = 0.0; maxblk = 0;
	    for(i=1; i<=WG->maxblk; i++){
		if(WG->map[i][j] > max){
			max = WG->map[i][j]; maxblk = i;
		}
	    }
	    fprintf(fp,"%4d | ",j);
	    if(max > 0){
	      k = (Int4)(max/inc);
	      for(s = 0; s < k; s++) fprintf(fp," ");
	      fprintf(fp,"o %d (%.1f) #%d\n",maxblk,max,WG->n[j]);
	    } else fprintf(fp,"\n");
	}
	fprintf(fp,"\n\n");
   }
}

void    NilWatchGibbs(wg_type WG)
{
	Int4	i;

	for(i=1; i<=WG->maxblk; i++) free(WG->map[i]);
	free(WG->map);
	free(WG->n);
	free(WG);
}

BooLean	AddWatchGibbs(float map, Int4 blk, Int4 col, wg_type WG)
{
	if(blk < 1 || blk > WG->maxblk || col < 1 || col > WG->maxcol){
		fprintf(stderr,"blk=%d (%d max); col=%d (%d max)\n",
			blk,WG->maxblk,col, WG->maxcol);
		fprintf(stderr,"\n!!AddWatchGibbs( ) input error\n\n");
		// print_error("AddWatchGibbs() input error");
		return FALSE;
	}
	if(map > WG->bestmap) { WG->bestmap = map; WG->bestcol= col; }
	if(map > WG->map[blk][col]){
	   if(WG->lowcol > col) WG->lowcol = col;
	   if(WG->hicol < col) WG->hicol = col;
	   WG->N++;
	   WG->n[col] = WG->N;
	   WG->map[blk][col] = map; 
	   return TRUE;
	} else return FALSE;
}

