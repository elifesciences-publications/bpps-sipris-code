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

#include "gmb_typ.h"

BooLean	SlideColLeftCMSA(Int4 t, double *oldMap, double *newMap, cma_typ L)
{
}

BooLean	SlideColRightCMSA(Int4 t, double *oldMap, double *newMap, cma_typ L)
{
}

BooLean	gmb_typ::SlideColLeft(Int4 t, double *oldMap, double *newMap,
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

   --[**.*.*]--[*.:**.*...***]--
               |		This operation allows escape from some 
               V                 nasty traps in alignment space.
   --[**.*.*.*]---[**.*...***]--

  returns TRUE if succeeds in sampling ends, otherwise returns FALSE.
 *************************************************************************/
{ 
    st_type		sites = SitesCMSA(L);
    fm_type		M,*model=ModelsCMSA(L);
    Int4		t1,t2,n,s,end,N,ntyp;
    Int4 		d0,d,i,j,e,start,length,offset,diff;
    double		map0,map,rand_no,*prob,total,ratio;
    double		*tmp_map;
    unsigned char	*seq;
    BooLean		right;
    
    ntyp=nTypeSites(sites);
    if(ntyp < 2) return FALSE; else if(t < 1 || t >= ntyp) return FALSE;
    t1=t; t2 = t1+1; // slide left: source = t2, denstination = t1.

    // 2. Make sure there are enough columns.
    if(nColsFModel(model[t2]) <= 3) return FALSE;
    map0 = RelMapCMSA(L); /** 1.b. Determine the alignment map. **/
    *oldMap=map0;

        // 3. Remove the last column from the source model. 
	M = model[t2]; length = LenFModel(M);
	d0=RmColumnMSA(t2,1,L);
	diff = length - LenFModel(M);
	length = LenFModel(model[t1]);

	// 4. Add the column to destination model.
	end = length + diff;
	tmp_map = new double [end+3];
	NEW(prob,end+3,double);
	for(total=1.0, i=length+1; i<=end; i++){
		AddColumnMSA(t1,i,L);
		map = RelMapCMSA(L);
		tmp_map[i]=map;
		ratio = exp(map-map0);	/** ratio of new to old **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		prob[i] = ratio; total += ratio; 
		s = LenFModel(model[t1]); RmColumnMSA(t1,s,L);
	}
	// 5. Sample an alignment.
	rand_no = ((double)Random()/(double)RANDOM_MAX)*total;
	for(i=length+1; i<=end; i++){
		rand_no -= prob[i];
		if(rand_no <= 0.0){
			AddColumnMSA(t1,i,L);
// fprintf(stderr,"slide columns right to left (diff = %d)\n",diff);
			*newMap=tmp_map[i]; delete [] tmp_map;
			free(prob); return TRUE;
		} 
	}
	delete [] tmp_map;
	AddColumnMSA(t2,1-diff,L); free(prob); return FALSE;
}

BooLean	gmb_typ::SlideColRight(Int4 t, double *oldMap, double *newMap, 
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

   --[**.*.*.*]---[**.*...***]--
               |		This operation allows escape from some 
               V                 nasty traps in alignment space.
   --[**.*.*]--[*.:**.*...***]--

  returns TRUE if succeeds in sampling ends, otherwise returns FALSE.
 *************************************************************************/
{ 
    st_type		sites = SitesCMSA(L);
    fm_type		M,*model=ModelsCMSA(L);
    Int4		t1,t2,n,s,end,N,ntyp;
    Int4 		d0,d,i,e,start,length,offset,diff;
    double		map0,map,rand_no,*prob,total,ratio;
    unsigned char	*seq;
    BooLean		right;
    double		*tmp_map;
    
    ntyp=nTypeSites(sites);
    if(ntyp < 2) return FALSE; else if(t < 1 || t >= ntyp) return FALSE;
    t1=t; t2 = t1+1;
    if(nColsFModel(model[t1]) <= 3) return FALSE;
    map0 = RelMapCMSA(L); /** 1.b. Determine the alignment map. **/
    *oldMap=map0;
    // source = t1, denstination = t2.
    	M = model[t1]; length = LenFModel(M);
        /** 4. Remove the last column from the source model. **/
	RmColumnMSA(t1,length,L);
	diff = length - LenFModel(M);
	NEW(prob,diff+3,double);
	tmp_map = new double [diff+2];

	/** 5. Add the columns to destination model **/
	for(total=1.0, i=1; i<=diff; i++){
		d=AddColumnMSA(t2,-(i-1),L);
		map = RelMapCMSA(L);
		tmp_map[i]=map;
		ratio = exp(map-map0);	/** ratio of new to old **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		prob[i] = ratio; total += ratio; 
		RmColumnMSA(t2,1,L);
	}
	/** sample an alignment **/
	rand_no = ((double)Random()/(double)RANDOM_MAX)*total;
	for(i=1; i <= diff; i++){
		rand_no -= prob[i];
		if(rand_no <= 0.0){
			AddColumnMSA(t2,-(i-1),L);
#if 1
fprintf(stderr,"slide columns from left to right adjacent models (diff = %d)\n",diff);
#endif
			*newMap=tmp_map[i]; delete [] tmp_map;
			free(prob); return TRUE;
		} 
	} delete [] tmp_map;
	AddColumnMSA(t1,length,L); free(prob);
	return FALSE;
}

