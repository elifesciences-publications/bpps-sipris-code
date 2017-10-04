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

#include "cmsa.h"

static Int4 colmultfactor = 3;
// static Int4 MAX_FLANK_CMSA = 5;

char	*GapAlignSeqCMSA(FILE *ftpr,Int4 a, Int4 b, Int4 *Score, e_type  E,
	Int4 **gapscore,cma_typ cmsa)
// aligns a sequence against an fmodel of the input cmsa object.
{
        Int4    n,t,ntyps,score=0,*pos;
        smx_typ *smx;

        ntyps=nBlksCMSA(cmsa);
        NEW(smx,ntyps+1,smx_typ);
        NEW(pos,ntyps+1,Int4);
#if 1	// weighted model
	a_type	A=AlphabetCMSA(cmsa);
        Int4 maxrpts=1;
	char method='H';
	char mode='O';
        Int4 maxLength=LenSeq(E)+10;
        BooLean weight=TRUE;
        float minmap=-99999990.0;	// Don't eliminate poorly conserved blocks.
        double pseudo=0.5;
	double minprob=1.0;

	FILE *fptr=tmpfile(); assert(fptr);
	PrintAlnCMSA(fptr,"tempfile",cmsa,NULL);	// output scan file...
	rewind(fptr);
        ptm_typ PM=MakePrtnModel(fptr, A, maxrpts, method, mode, maxLength,
                        weight, minmap,pseudo,minprob,blosum62freq);
	fclose(fptr);
        for(t=1; t<=ntyps; t++){
	   smx[t]=SMatrixPrtnModel(t,PM);
	   if(smx[t]==0){
		fprintf(stderr,"block %d\n",t);
		 print_error("Fatal: negative map block in cma file");
	   }
	}
#else
        double  pernats=1000.0;
        for(t=1; t<=ntyps; t++){
                smx[t]=GetSmatrixFModel(pernats,ModelCMSA(t,cmsa));
        }
#endif
        // score=AlnSeqSMatrix(LenSeq(E),SeqPtr(E),ntyps,smx,pos);
	char    *operation=GapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),ntyps,smx,gapscore);
        if(ftpr) PutSeqInfo2(ftpr,E);
	fprintf(stderr,"length operation = %d; length seq = %d\n",strlen(operation),
		LenSeq(E));
#if 1
	if(ftpr) PutSeqAlnSMatrixSW(ftpr,a,b,LenSeq(E),SeqPtr(E),ntyps,smx,0);
#else
        for(t=1; t<=ntyps; t++){
	   // NilSMatrix(smx[t]);
        } if(ftpr) fprintf(ftpr,"score = %d\n",score);
#endif
        for(t=1; t<=LenSeq(E); t++){
		fprintf(stderr,"%c%d: %c\n",AlphaChar(ResSeq(t,E),AlphabetCMSA(cmsa)),
			t+OffSetSeq(E),operation[t-1]);
	}
        *Score=score;
        free(smx);
        return operation;
}

double	ReAlignBestCMSA(cma_typ cmsa)
// relalign sequences against an fmodel of the rest of the sequences
{
	Int4	n,t,ntyps,score,*pos;
	smx_typ	*smx;
#if 0
	double	pernats=1000.0;
#else
	gss_typ	*gss=gssCMSA(cmsa);
	double	pernats=gss->PerNats();;
#endif
	e_type	E;

   for(n=NumSeqsCMSA(cmsa); n > 0;  n--){
	ntyps=nBlksCMSA(cmsa);
	E = SeqSetE(n,DataCMSA(cmsa));
	NEW(smx,ntyps+1,smx_typ);
	NEW(pos,ntyps+1,Int4);
	VacateSitesCMSA(n,cmsa);
	for(t=1; t<=ntyps; t++){
		smx[t]=GetSmatrixFModel(pernats,ModelCMSA(t,cmsa));
	}
	score=AlnSeqSMatrix(LenSeq(E),SeqPtr(E),ntyps,smx,pos);
	// score=AlnSeqSMatrix(LenSeq(E),XSeqPtr(E),ntyps,smx,pos);
	// PutSeqInfo2(stderr,E);
	for(t=1; t<=ntyps; t++){
		// PutSeqRegion(stderr,pos[t],LengthCMSA(t,cmsa),E,AlphabetCMSA(cmsa));
		// fprintf(stderr,"\n");
		AddSiteCMSA(t,n,pos[t], cmsa);
		NilSMatrix(smx[t]);
	}
	// fprintf(stderr,"\nmap = %g\n",RelMapCMSA(cmsa));
	free(pos); free(smx);
   }
	return RelMapCMSA(cmsa);
}

BooLean	comp_map_cmsa(double nmap, double omap, double temperature)
/** short routine to be used for sampling **/
{
	double	ratio;

	// if(fabs(nmap-omap) > 500) print_error("error in comp_map_cmsa( )");
	if((nmap-omap) > 40)  return TRUE;
	else if((omap-nmap) > 40)  return FALSE;
	if(omap < nmap){
		ratio = exp(nmap-omap);	/** ratio of new to old **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		ratio = ((double)Random()/(double)RANDOM_MAX)*(ratio+1);
		if(ratio >= 1.) return TRUE;
		else return FALSE;
	} else {	/** omap >= nmap **/
		ratio = exp(omap-nmap);	/** ratio of old to new **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		ratio = ((double)Random()/(double)RANDOM_MAX)*(ratio+1);
		if(ratio > 1.) return FALSE;
		else return TRUE;
	}
}

/******************************** columns.c **********************************/

BooLean	SampAddColCMSA(Int4 t, cma_typ L)
{ return SampAddColTempCMSA(t, 1.0, L); }

BooLean	SampAddColTempCMSA(Int4 t, double temperature, cma_typ L)
/*************************************************************************
  Sample an additional column for block t.  If one is added then return 
  TRUE; otherwise return FALSE.
/*************************************************************************/
{
	Int4	New,d;
	double	omap,nmap;
	FILE 	*efp=0;

	omap=RelMapCMSA(L); New=SampNewColMSA(t,L);
	if(New==1) return FALSE;	/** no room to add a column **/
	d=AddColumnMSA(t, New, L);	/** d = displacement of model **/
	nmap=RelMapCMSA(L);
	if(comp_map_cmsa(nmap, omap,temperature)){
		if(efp) fprintf(stderr,"Sampled new column %d (%d) in block %d: delta LLR = %.2f\n",
			New,New-d,t,nmap-omap);
		return TRUE;
	}
	RmColumnMSA(t, New-d, L);
	return FALSE;
}

BooLean	SampRmColCMSA(Int4 t, cma_typ L)
{ return SampRmColTempCMSA(t,1.0,L); }

BooLean	SampRmColTempCMSA(Int4 t, double temperature,cma_typ L)
/*************************************************************************
  Sample a column to remove from block t.  If one is removed then return 
  TRUE; otherwise return FALSE.

   * . * * . *            . . * * . *          
  |1|2|3|4|5|6|...  -->  |x|x|3|4|5|6|...  removed column 1

	lemon = 1; d = +2.  then revert to original:

  lemon-d = 1 - 2 = -1.	   (Note: 0 = column just to left of start of block)

/*************************************************************************/
{
	Int4	lemon,d;
	double	omap,nmap;
	fm_type	M = ModelCMSA(t,L);

	if(nColsFModel(M) < 4) return FALSE;
	omap = RelMapCMSA(L);
	if(temperature < 1.0) lemon = OrangeFModel(M);	/** 1..len = block **/
	else lemon = LemonFModel(M);	/** 1..len = block **/
	if(lemon == INT4_MIN) return FALSE;  // temporary fix for overflow.
	d = RmColumnMSA(t, lemon, L);	/** d = displacement of model **/
	nmap=RelMapCMSA(L);
	if(comp_map_cmsa(nmap, omap,temperature)) return TRUE;
	AddColumnMSA(t, lemon-d, L); /** revert to original **/ 
	return FALSE;
}

BooLean	SampSlideColLtCMSA(Int4 t, double *oldMap, double *newMap,
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

   --[**.*.*]--[*.:**.*...***]--
               |
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

BooLean	SampSlideColRtCMSA(Int4 t, double *oldMap, double *newMap, 
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

   --[**.*.*.*]---[**.*...***]--
               |
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


Int4	SampNewColMSA(Int4 t, cma_typ L)
/*******************************************************************
  Sample a new column (to add to model); return the new column sampled.
  If the sampling fails 1 is returned.
********************************************************************/
{
	BooLean	debug=0;
	Int4	lemon;
	st_type S=L->sites;
	Int4	**site_freq,ncol,start,fend,end,i,j,d,k,flank,leng,maxlen;
	double	*ratio,total,rand_no,weight;
	fm_type M=ModelCMSA(t,L);
	a_type	AB=AlphabetCMSA(L);
	FILE	*efp=0;

	lemon = LemonFModel(ModelCMSA(t,L)); /** proportional to prob **/
#if 0
	lemon = OrangeFModel(ModelCMSA(t,L)); /** at random **/
#endif
	if(lemon == INT4_MIN) return 1;  // temporary fix for overflow.
	
	ncol = nColsFModel(M);
	// OLD: flank = (L->maxlen[t] - LenFModel(M))/2; 
	// NEW: don't allow model to get larger than colmultfactor x numcolumns.
	maxlen = MINIMUM(Int4, L->maxlen[t], colmultfactor*ncol);
	flank = (maxlen - LenFModel(M))/2;
	flank = MINIMUM(Int4, flank, MAX_FLANK_CMSA);
	if(flank < 0) flank = 0;
	/**** end NEW. ****/
	fend = end =  LenFModel(M) + flank;
	leng = LenFModel(M); 
        NEWP(site_freq, end+MAX_FLANK_CMSA+3,Int4);
        NEW(ratio, end+MAX_FLANK_CMSA+3,double);
	site_freq += flank; ratio += flank; 
	for(total=0.0,i=0; i>= -flank; i--){
		site_freq[i] = GetSiteFreq(S,t,i-1);
		if(site_freq[i] == NULL) break;
		ratio[i] = LnRatioFModel(site_freq[i],lemon, M);
		weight = add_col_weight(ncol, leng, i);
		ratio[i] += log(weight); 
// the following fixes overflow problems...
		if(ratio[i] > 27.6){ 	// e^27.6 == one trillion..
		  if(debug){
		   Int4 r,sum=0;
		   for(r=1; r <= nAlpha(AB); r++){
			fprintf(stderr,"%c=%d ",AlphaChar(r,AB),site_freq[i][r]);
			if(r%10 == 0) fprintf(stderr,"\n");
			sum+= site_freq[i][r];
		   } fprintf(stderr," sum = %d\n",sum);
		   fprintf(stderr,"blk=%d(%d/%d); i=%d; lemon=%d; weight=%g\n",t,ncol,leng,i,lemon,weight);
		   Int4 *lemon_freq = GetSiteFreq(S,t,lemon-1);
		   sum=0;
		   for(r=1; r <= nAlpha(AB); r++){
			fprintf(stderr,"%c=%d ",AlphaChar(r,AB),lemon_freq[r]);
			if(r%10 == 0) fprintf(stderr,"\n");
			sum+= site_freq[i][r];
		   } fprintf(stderr," sum = %d\n",sum);
		   free(lemon_freq);
		  }
		  if(efp) fprintf(stderr,
		    "SampNewColMSA( ): %g = ratio for blk %d column %d reset to 1 trillion to 1.\n",
		     ratio[i],t,i);
		  ratio[i] = exp(27.6) + ratio[i];	// make slightly more probable.
		} else ratio[i] = exp(ratio[i]);
		// i.e. new column is MUCH better than old...
		// No need to go beyond 100 million million to 1 odds!
// NEED BETTER SOLUTION TO THIS!!!
		total += ratio[i];
	} start = i+1;
	site_freq[1] = NULL;
	for(i=2; i<= end; i++){
	   if(NullSiteFModel(i,M)){
		site_freq[i] = GetSiteFreq(S,t,i-1);
		if(site_freq[i] == NULL) { end = i-1; break; }
		ratio[i] = LnRatioFModel(site_freq[i],lemon, M);
		weight = add_col_weight(ncol, leng, i);
		ratio[i] += log(weight); 
// the following fixes overflow problems...
		if(ratio[i] > 27.6){ 	// e^27.6 == one trillion..
		  if(debug){
		   Int4 r,sum=0;
		   for(r=1; r <= nAlpha(AB); r++){
			fprintf(stderr,"%c=%d ",AlphaChar(r,AB),site_freq[i][r]);
			if(r%10 == 0) fprintf(stderr,"\n");
			sum+= site_freq[i][r];
		   } fprintf(stderr," sum = %d\n",sum);
		   fprintf(stderr,"blk=%d(%d/%d); i=%d; lemon=%d; weight=%g\n",t,ncol,leng,i,lemon,weight);
		   Int4 *lemon_freq = GetSiteFreq(S,t,lemon-1);
		   sum=0;
		   for(r=1; r <= nAlpha(AB); r++){
			fprintf(stderr,"%c=%d ",AlphaChar(r,AB),lemon_freq[r]);
			if(r%10 == 0) fprintf(stderr,"\n");
			sum+= site_freq[i][r];
		   } fprintf(stderr," sum = %d\n",sum);
		   double *Ps=PseudoFModel(M);
		   for(r=1; r <= nAlpha(AB); r++){
			fprintf(stderr,"%c=%.1f ",AlphaChar(r,AB),Ps[r]);
			if(r%10 == 0) fprintf(stderr,"\n");
		   } 
		   free(lemon_freq);
		  }
		  if(efp) fprintf(stderr,
		    "SampNewColMSA( ): %g = ratio for blk %d column %d reset to 1 trillion to 1.\n",
		     ratio[i],t,i);
		  ratio[i] = exp(27.6) + ratio[i];	// make slightly more probable.
		} else ratio[i] = exp(ratio[i]);
		// i.e. new column is MUCH better than old...
// NEED BETTER SOLUTION TO THIS!!!
		total += ratio[i];
	   } else { site_freq[i] = NULL; }
	}
	rand_no = ((double)Random()/(double)RANDOM_MAX)*total;
	for(i=start; i <= end; i++){
	   if(NullSiteFModel(i,M)){
	      if((rand_no -= ratio[i]) <= 0.0) {
		for(j=-flank;j<=fend;j++){
		    if(site_freq[j]!=NULL){
			free(site_freq[j]);site_freq[j]=NULL;
		    }
		}
		site_freq-=flank; ratio-=flank;
		free(site_freq); free(ratio);
		return i;
	      }
	   }
	}
	if(total==0.0){
	   site_freq-=flank; ratio-=flank;
	   free(site_freq); free(ratio);
	   return 1;	/** no room to add a column **/
	}
	/** else an error has occurred **/
	fprintf(stderr,"total = %g; flank = %d; start=%d; end=%d\n",
		total,flank,start,end);
	for(i=start; i <= end; i++){
	   if(NullSiteFModel(i,M)){
	      fprintf(stderr,"ratio[%d]=%g\n",i,ratio[i]);
	   }
	} NullSitesFModel(L->null, M);
	PutSites(stderr,t,S,NULL,L->null);
	print_error("error in SampNewColMSA(): this should not happen");
	return 0;  // to stop compiler from complaining.
}

cma_typ DeleteBlkCMSA(Int4 b,cma_typ L){ return delete_blk_cmsa(L,b); }

cma_typ	DeleteBlkCMSA(cma_typ L)
/*********************************************************************
  If it improves the map to delete a block then delete it and return
  the resultant msa.
/*********************************************************************/
{
	Int4	tdel,t;
	double	map,bmap;

	bmap = RelMapCMSA(L); tdel=0;
	for(t=1; t <= nBlksCMSA(L); t++){
		map=RelMapMinusCMSA(t, L);
		if(map > bmap){ tdel=t; bmap=map; }
	}
	if(tdel==0) return NULL; else return delete_blk_cmsa(L, tdel);
}

cma_typ	RefineCMSA(cma_typ msa)
/*********************************************************************
 Refine model by deleting bogus blocks.
 MAP < 0.0 -> model is no more likely than background.
 i.e., Eliminate insignificant columns?
/*********************************************************************/
{
	Int4	tdel,t;
	double	map,wmap;

	for(wmap=999,tdel=0,t=1; t <= nBlksCMSA(msa); t++){
		map=FieldRelMapCMSA(msa, t);
		if(map < 0.0 && map < wmap){ tdel=t; wmap=map; }
	}
	if(tdel==0) return NULL;
	else return delete_blk_cmsa(msa, tdel);
}

#if 0	/** tested refinement based on a weight matrix; replace smodel with fwmodel. **/
void	RefineCMSA2(char method, cma_typ msa)
/*********************************************************************
 Refine model using a blosum62 based profile.
/*********************************************************************/
{
   e_type	E;
   ss_type	data=DataCMSA(msa);
   st_type	S=SitesCMSA(msa);
   fm_type	*model;
   sm_type	*M;
   smx_typ	*sM;
   a_type	A=AlphabetCMSA(msa);
   unsigned char	*seq;
   char		c;
   Int4		score,tot_len,n,N,*len,nseq,nblks,cycle;
   Int4		i,j,m,s,*oldsite,*site,maxcycle=20;
   BooLean	changed,done,use_null=TRUE,**null;

   nseq = NSeqsSeqSet(data);
   nblks = nBlksCMSA(msa);
   if(isupper(method)){ use_null=FALSE; method=tolower(method); }
   else {
	model = ModelsCMSA(msa);
	NEWP(null,nblks+2,BooLean);
        for(m=1; m <= nblks; m++){
	    NEW(null[m],MaxSeqSeqSet(data)+2,BooLean);
	    NullSitesFModel(null[m], model[m]);
	}
   }
   fprintf(stderr,"creating profile\n");
   /*** Create profile of alignment. ***/
   NEW(M,nblks+2,sm_type); NEW(len,nblks+2,Int4); 
   for(tot_len=0,m=1; m <= nblks; m++){
	len[m] = SiteLen(m,S); tot_len += len[m];
	if(use_null) M[m]=MkSModel(null[m],len[m],PSEUDO_CMSA,CntsSeqSet(data),A);
	else M[m]=MkSModel(NULL,len[m],PSEUDO_CMSA,CntsSeqSet(data),A);
	SetMethodSModel(method,M[m]);
	for(s=1; s<=nseq; s++){ 
		E = SeqSetE(s,data); seq = SeqPtr(E);
		i = SitePos(m,s,1,S);
		Add2SModel(seq,i,M[m]); 
	}
   }
   fprintf(stderr,"checking for errors\n");
   /*** check for errors ***/
   for(s=1; s<=nseq;s++) {
	  E = SeqSetE(s,data);
	  if(tot_len > (Int4) LenSeq(E)) print_error("RefineCMSA( ) error");
   }
   fprintf(stderr,"refining model\n");
   /*** Iteratively refine model and new msa **/
   NEW(site,nblks+2,Int4); NEW(oldsite,nblks+2,Int4); 
   NEW(sM,nblks+2,smx_typ);
   for(done=FALSE,n=0,cycle=0; cycle <=maxcycle && !done; cycle++){
        for(s=1; s<=nseq;s++) {
	  E = SeqSetE(s,data); seq = SeqPtr(E);
	  for(m=1; m<=nblks; m++) { 
	        oldsite[m]=SitePos(m,s,1,S);
		RmSiteCMSA(m,s,oldsite[m],msa);
		RmSModel(seq, oldsite[m], M[m]);
		sM[m]=GetSMatrixSModel(M[m]); 
	  }
	  score=AlnSeqSMatrix(LenSeq(E),SeqPtr(E),nblks,sM,site);
	  for(changed=FALSE,m=1; m<=nblks; m++){ 
		AddSiteCMSA(m,s,site[m],msa); 
		Add2SModel(seq,site[m],M[m]); 
		if(oldsite[m] != site[m]){ changed=TRUE; }
	  }
	  if(changed) n=0;
	  else n++;
	  if(n > nseq) { done=TRUE; break; }
	}
	fprintf(stderr,"cycle %d; n = %d\n",cycle,n);
    }
    for(m=1; m <= nblks; m++) NilSModel(M[m]); free(M);
    free(oldsite); free(site); free(len); free(sM); 
    if(use_null){
        for(m=1; m <= nblks; m++) free(null[m]); free(null);
    }
    fprintf(stderr,"done with refining model\n");
}
#endif

BooLean	MoveMultiColsMSA(cma_typ L, Int4 t, Int4 num)
{
	Int4	i,lemon;
	BooLean	result=FALSE;

	for(i=1; i<=num; i++){
	   lemon = ChoiceLemonFModel(ModelCMSA(t,L));
	   if(lemon == INT4_MIN) continue;  // temporary fix for overflow.
	   // lemon = ChoiceOrangeFModel(ModelCMSA(t,L)); 
	   if(move_column_msa(L, lemon, t)) { result=TRUE; } else break;
	}
	return result;
}

BooLean MoveColumnCMSA(cma_typ L, Int4 t)
{
	Int4	lemon;
	
	lemon = LemonFModel(ModelCMSA(t,L)); /** proportional to prob **/
	if(lemon == INT4_MIN) return FALSE;  // temporary fix for overflow.
	// lemon = OrangeFModel(ModelCMSA(t,L)); /*** randomly ***/
	return move_column_msa(L, lemon, t);
}

BooLean	move_column_msa(cma_typ L, Int4 lemon, Int4 t)
/*******************************************************************
  Sample a column to remove and then sample a column to replace it.
  if lemon = new site then don't bother to move ??? 
********************************************************************/
{
	st_type S=L->sites;
	Int4	**site_freq,ncol,start,fend,end,i,j,d,k,flank,leng,maxlen;
	double	*ratio,total,rand_no,weight;
	fm_type M=ModelCMSA(t,L);
	mdl_typ	*mdl=mdlCMSA(L);

	ncol = nColsFModel(M);
	NullSitesFModel(L->null, M);
	// OLD: flank = (L->maxlen[t] - LenFModel(M))/2;
	// NEW: don't allow model to get larger than colmultfactor x numcolumns.
	maxlen = MINIMUM(Int4, L->maxlen[t], colmultfactor*ncol);
	flank = (maxlen - LenFModel(M))/2;
	flank = MINIMUM(Int4, flank, MAX_FLANK_CMSA);
	if(flank < 0) flank = 0;
	/**** end NEW. ****/
	fend = end =  LenFModel(M) + flank;
	leng = LenFModel(M); 
        NEWP(site_freq, end+MAX_FLANK_CMSA+3,Int4);
        NEW(ratio, end+MAX_FLANK_CMSA+3,double);
	site_freq += flank; ratio += flank; 
	for(total=1.0,i=0; i>= -flank; i--){ 
		site_freq[i] = GetSiteFreq(S,t,i-1); 
		if(site_freq[i] == NULL) break;
		ratio[i] = RatioFModel(site_freq[i],lemon, M);
		weight = mv_col_weight(ncol, leng, i, lemon, L->null);
		ratio[i] *= weight; total += ratio[i];
	} start = i+1;
	site_freq[1] = NULL;
	for(i=2; i<= end; i++){
	   if(NullSiteFModel(i,M)){
		site_freq[i] = GetSiteFreq(S,t,i-1);
		if(site_freq[i] == NULL) { end = i-1; break; }
		ratio[i] = RatioFModel(site_freq[i],lemon, M);
		weight = mv_col_weight(ncol, leng, i, lemon, L->null);
		ratio[i] *= weight; total += ratio[i];
	   } else { site_freq[i] = NULL; }
	}
	rand_no = ((double)Random()/(double)RANDOM_MAX)*total;
	for(i=start; i <= end; i++){
	   if(NullSiteFModel(i,M)){
	      if((rand_no -= ratio[i]) <= 0.0) {
		d = mdl->MvColumn(t,site_freq[i], lemon,i); 
		site_freq[i]= NULL;
		/**** fprintf(stderr," (new=%d; old=%d; w=%d) %g",
			i,lemon,LenFModel(M), ratio[i]); /****/
		if(LenFModel(M) < leng){	/* model shrunk */
			k = leng - LenFModel(M);
			for(j = 1; j <= k; j++) { ShrinkSites(t,S); }
			if(d != 0) ShiftSitesM(S, t, d);
		}else if(LenFModel(M) > leng){	/* model has grown */
			k = LenFModel(M) - leng;
			if(d != 0) ShiftSitesM(S, t, d);
			for(j = 1; j <= k; j++) { GrowSites(t,S); }
		} else if(d!=0) ShiftSitesM(S,t,d);
		for(j=-flank;j<=fend;j++){
			if(site_freq[j]!=NULL) {
				free(site_freq[j]);site_freq[j]=NULL;
			}
		}
	        site_freq-=flank; ratio-=flank;
		free(site_freq); free(ratio);
		return TRUE; 
	      }
	   }
	}
	for(j=-flank; j<=fend; j++){
	   if(site_freq[j]!=NULL){free(site_freq[j]);site_freq[j]=NULL;}
	}
	site_freq-=flank; ratio-=flank;
	free(site_freq); free(ratio);
	return FALSE;
}

double	transfer_weight(double wd, double cd, double nwr, double owr, 
	double cr)
// WARNING: nwr is never referenced!!!  Check this out.
{
	double weight;

	/******************* weight adjustment ********************
	weight donor(M1) < 1; weight recipient(M2) > 1.
	donor: (w-2 choose c-2-1)/(w-2 choose c-2) = (c-2)/(w-2-c+3).
	recipient: (w-2 choose c-2+1)/(w-2 choose c-2) = (w-2-c+2)/(c-1).
	net = donor * recipient
	 ******************* weight adjustment ********************/
	weight = (cd-2)/(wd-cd+1);
	weight *= (owr-cr)/(cr-1);
	return weight;
}

BooLean TransferColumnCMSA(cma_typ L)
{
	st_type S=L->sites;
	Int4	**site_freq,t1,t2,start,leng1,leng2;
	Int4	fend,end,i,j,d,k,flank,lemon,ntyp;
	double	*ratio,total,rand_no,weight,c1,c2;
	fm_type	M1,M2,*model;
	mdl_typ *mdl=mdlCMSA(L);

	ntyp=nTypeSites(S);
	model = ModelsCMSA(L);
	do{	/* select a model for lemon */
		rand_no = (double)Random()/(double)RANDOM_MAX;
		rand_no *= (double) ntyp;
		t1 = 1 + (Int4) rand_no;
	} while(t1 < 1 || t1 > ntyp);
	M1 = model[t1];

	do{	/* select a model for new column */
		rand_no = (double)Random()/(double)RANDOM_MAX;
		rand_no *= (double) ntyp;
		t2 = 1 + (Int4) rand_no;
	} while(t2 < 1 || t2 > ntyp || t2 == t1);
	M2 = model[t2];

	leng1 = LenFModel(M1);
	if(nColsFModel(M1) == 3) return FALSE;
#if 1
	lemon = LemonFModel(M1); /** sample proportional to prob. **/
#endif
#if 0
	lemon = OrangeFModel(M1); /** sample at random **/
	if(lemon == 1 || lemon == leng1) return FALSE;
#endif
	if(lemon == INT4_MIN) return FALSE;  // temporary fix for overflow.

	leng2 = LenFModel(M2);
	// OLD: flank = ((L->maxlen[t2] - LenFModel(M2))/2); 
	flank = 0; /*** need to weight columns if changes model size ***/
	/******************* weight adjustment ********************
	weight donor(M1) < 1; weight recipient(M2) > 1.
	donor: (w-2 choose c-2-1)/(w-2 choose c-2) = (c-2)/(w-2-c+3).
	recipient: (w-2 choose c-2+1)/(w-2 choose c-2) = (w-2-c+2)/(c-1).
	net = donor * recipient
	 ******************* weight adjustment ********************/
	c1 = nColsFModel(M1); c2 = nColsFModel(M2);
	if(c2 == leng2) return FALSE;
	if((leng1 < (Int4)c1) || (leng2 < (Int4)c2)) {
		fprintf(stderr,"w1=%d; w2=%d; c1 =%g; c2 = %g\n",
			leng1,leng2,c1,c2);
		print_error("what the...?");
	}
	weight = (c1-2)/(leng1-c1+1);
	weight *= (leng2-c2)/(c2-1);
	/*** weight=1; /******
	fprintf(stderr,"weight = %g\n",weight);
	/*********************************************************/
	fend = end = leng2 + flank;
        NEWP(site_freq, end+MAX_FLANK_CMSA+3,Int4);
        NEW(ratio, end+MAX_FLANK_CMSA+3,double);
	site_freq += flank; ratio += flank; site_freq[1] = NULL;
	for(start=1,total=1.0,i=0; flank && i>= -flank; i--){
		site_freq[i] = GetSiteFreq(S,t2,i-1);
		if(site_freq[i] == NULL) { start = i+1; break; }
		ratio[i] = weight*RatioFModel(site_freq[i],lemon, M1);
		total += ratio[i];
	}
	site_freq[1] = NULL;
	for(i=2; i<= end; i++){
	   if(NullSiteFModel(i,M2)){
		site_freq[i] = GetSiteFreq(S,t2,i-1);
		if(site_freq[i] == NULL) { end = i-1; break; }
		ratio[i] = weight*RatioFModel(site_freq[i],lemon, M1);
		total += ratio[i];
	   } else { site_freq[i] = NULL; }
	}
	rand_no = (double)Random()/(double)RANDOM_MAX;
	rand_no *= total;
	for(total=0.0, i=start; i <= end; i++){
	   if(NullSiteFModel(i,M2)){
	      if((rand_no -= ratio[i]) <= 0.0) {
	if(t1!=t2) {  /* t1 can't be == t2  - see above. */
	   /*** FIRST REMOVE OLD COLUMN ***/
		d = mdl->RmColumn(t1,lemon);
#if 0
		fprintf(stderr," (old=%d; length=%d) %g model %c\n",
			lemon,LenFModel(M1), ratio[i], t1+'A'-1); 
#endif
		if(LenFModel(M1) < leng1){	/* model shrunk */
			k = leng1 - LenFModel(M1);
			for(j = 1; j <= k; j++) { ShrinkSites(t1,S); }
			if(d != 0) ShiftSitesM(S, t1, d);
		} else if(LenFModel(M1) > leng1){	/* model has grown */
			print_error("this should not happen");
		} else if(d!=0) ShiftSitesM(S, t1, d);

	  /*** THEN ADD NEW COLUMN ***/
		d = mdl->AddColumn(t2,site_freq[i], i);
		site_freq[i]= NULL;
#if 0
		fprintf(stderr," (new=%d; length=%d) %g model %c\n",
			i,LenFModel(M2), ratio[i],'A' + t2 - 1);
#endif
		if(LenFModel(M2) < leng2){	/* model shrunk */
			print_error("this should not happen");
		} else if(LenFModel(M2) > leng2){	/* model has grown */
			k = LenFModel(M2) - leng2;
			if(d != 0) ShiftSitesM(S, t2, d);
			for(j = 1; j <= k; j++) { GrowSites(t2,S); }
		} else if(d!=0) ShiftSitesM(S, t2, d);
	/*** FREE POTENTIAL COLUMNS ***/
		for(j=-flank; j<=fend; j++){
			if(site_freq[j]!=NULL) {
				free(site_freq[j]); site_freq[j] = NULL;
			}
		}
		site_freq-=flank; ratio-=flank;
		free(site_freq); free(ratio);
		return TRUE; 
	} else print_error("this should not happen!?");
	      }
	   }
	}
	for(j=-flank; j<=fend; j++){
		if(site_freq[j]!=NULL) { 
			free(site_freq[j]); site_freq[j] = NULL;
		}
	}
	site_freq-=flank; ratio-=flank;
	free(site_freq); free(ratio);
	return FALSE;
}

BooLean	ShiftCMSA(cma_typ L, Int4 t)
/* shift element to the left or right (allows shifts into other elements of 
   same type since they will also shift over in the same way).  */
{
	ss_type	P=DataCMSA(L);
	st_type S=L->sites;
	Int4	*site_freq[2][15],*stfreq,b,p,len;
	Int4	end[2],i,j,o,d,d2,n,k,shift,t2;
	double	ratio[2][15], total,rand_no;
	BooLean	left[2],*overlap;
	fm_type	M;
	a_type	A=AlphabetCMSA(L);
	e_type	E;
	mdl_typ	*mdl=mdlCMSA(L);
	/** char	c; /****/
	
	M = ModelCMSA(t,L);
	shift = MINIMUM(Int4,LenFModel(M)/2,10);
	shift = MAXIMUM(Int4,shift,1);
	NEW(overlap,shift+3,BooLean);
	end[0] = shift; left[0]=FALSE;	/* note: right == 0 == FALSE */
	end[1] = shift; left[1]=TRUE;	/* note: left == 1 ==TRUE */
	len = SiteLen(t,S);
	for(total=1.0,i=0; i<2; i++){
	   ratio[i][0] = 1.0;
	   for(d=1; d<= end[i]; d++){
		/*** = GetSiteFreq ***/
		MEW(stfreq, nAlpha(A)+2,Int4);
		for(b=0;b<= nAlpha(A); b++) stfreq[b] = 0;
		if(left[i]){	/*** [__site__] ===> ***/
		  d2 = d + len - 1;
		  for(n=1; n <=NSeqsSeqSet(P) && stfreq!= NULL; n++){
		    E = SeqSetE(n,P);
		    for(k=1; k<=nSites(t,n,S); k++){
			p = SitePos(t,n,k,S) + d2;
			t2 = TypeSite(n,p,S);
			overlap[d]=OpenPos(n,p,S);
			if(p<=(Int4)SqLenSeqSet(n,P) && 
			    (overlap[d] || t2 == t || t2==BlockedSite(t,S))){
			   // b=XSeq(p,E); stfreq[b]++;  /** ERROR **/
			   b=ResSeq(p,E); stfreq[b]++; 
			} else { free(stfreq); stfreq = NULL; break; }
		    }
		  }
		} else {	/*** <=== [__site__] ***/
		  d2 = -d;
		  for(n=1; n <=NSeqsSeqSet(P) && stfreq!= NULL; n++){
		    E = SeqSetE(n,P);
		    for(k=1; k<=nSites(t,n,S); k++){
			p = SitePos(t,n,k,S) + d2;
			t2 = TypeSite(n,p,S);
			overlap[d]=OpenPos(n,p,S);
			if(p>=1 && (overlap[d] || t2==BlockedSite(t,S))){
			   // b=XSeq(p,E); stfreq[b]++;  /** ERROR **/
			   b=ResSeq(p,E); stfreq[b]++; 
			} else { free(stfreq); stfreq = NULL; break;}
		    }
		  }
		}
		/*** end GetSiteFreq ***/
		if(stfreq == NULL) { end[i] = d-1; break; }
		else {
		  site_freq[i][d] = stfreq;
		  if(left[i]) d2 = d;
		  else d2 = SiteLen(t,S) - d + 1;
		  ratio[i][d]=RatioFModel(site_freq[i][d],d2, M);
		  ratio[i][d] *= ratio[i][d-1];
		  total += ratio[i][d];
		}
	   }
	}
	rand_no = ((double)Random()/(double)RANDOM_MAX)*total;
	for(total=0.0, i=0; i < 2; i++){
	   // if(left[i]) c= '+'; else c = '-'; 
	   for(d=1; d<= end[i]; d++){
	      if((total += ratio[i][d]) >= rand_no) {
		// fprintf(stderr,"\r["); 
		for(j=1; j<=d; j++){
		    // fprintf(stderr,"%c", c);
		    mdl->Shift(t,site_freq[i][j], left[i]); 
		    if(overlap[d]){
			ShrinkSites(t,S);
			ShiftSites(S,t,left[i]);;
			GrowSites(t, S);
		    } else { ShiftSites(S,t,left[i]);; }
		}
		// fprintf(stderr,"] %g ",ratio[i][d]); 
		for(j=d+1; j<=end[i]; j++){
			if(site_freq[i][j]!=NULL) free(site_freq[i][j]);
		}
		for(o=(i+1)%2,j=1; j<=end[o]; j++){
			if(site_freq[o][j]!=NULL) free(site_freq[o][j]);
		}
		free(overlap);
		return TRUE; 
	      }
	   }
	}
	for(i=0; i<2; i++){
	   for(d=1; d<= end[i]; d++){
		if(site_freq[i][d]!=NULL) free(site_freq[i][d]);
	   }
	}
	free(overlap);
	return FALSE;
}

/************************* SAMPLING INDELs **************************/

Int4	**gap_score_cmsa(cma_typ cma)
// obtain gaps except for sequence s
{
	Int4	**observedGap,N=NumSeqsCMSA(cma);
	Int4	g,maxLength=MaxTrueSeqCMSA(cma),M=nBlksCMSA(cma);
	a_type	A=AlphabetCMSA(cma);

        NEWP(observedGap, M+2, Int4);
        NEW(observedGap[0], maxLength+3, Int4);;
        for(Int4 m=1; m < M; m++) {
            NEW(observedGap[m], maxLength+3, Int4);;
	    for(Int4 s=1; s <= N; s++){
	       g=GapBetweenSites(s,m,SitesCMSA(cma)); observedGap[m][g]++;
	    }
        }
        NEW(observedGap[M], maxLength+3, Int4);;
	return observedGap;
}

BooLean	SampleGapsGibbsCMSA(Int4 s,cma_typ *oldcma,double Temperature)
// sample insertions and deletions in sequence s.
{
	Int4	m,n,i,t,start,trace_length,indels,gapopen,gapextend;
	Int4    score,**gapfnct,**observedGap;
	Int4	*newpos,*oldpos;
	smx_typ	*smx;
	a_type	A=AlphabetCMSA(*oldcma);
	e_type	oldE,E;
	gss_typ	*gss=gssCMSA(*oldcma);
	double	pernats;
	char	*operation;
	Int4	wt;
	cma_typ	cma;
	double	oldmap,newmap;

	cma=CopyCMSA(*oldcma);
	oldmap=RelMapCMSA(*oldcma);
	observedGap= gap_score_cmsa(cma); // 0. Obtain Gap Function.
	gapopen=gss->GapOpen(); gapextend=gss->GapExtend(); pernats=gss->PerNats();

        if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
        else if(Temperature > 275.0) wt=2;
        else if(Temperature > 250.0) wt=3;
        else if(Temperature > 225.0) wt=4;
        else if(Temperature > 200.0) wt=5;
        else if(Temperature > 175.0) wt=6;
        else if(Temperature > 150.0) wt=7;
        else if(Temperature > 125.0) wt=8;
        else if(Temperature > 100.0) wt=9;
        else if(Temperature > 75.0) wt=10;
        else wt=0;  // Use maximum likelihood.

	// 1. Remove sequence s from alignment.
	NEW(oldpos,nBlksCMSA(cma)+2,Int4);  // save old sites.
	NEW(newpos,nBlksCMSA(cma)+2,Int4);  // for new sites.
	for(t=1; t<=nBlksCMSA(cma); t++){
		PosSiteCMSA(t, s, cma->pos, cma); oldpos[t]=cma->pos[1];
	} VacateSitesCMSA(s,cma);

	// ===== 2. Obtain scoring matrix from rest of alignment.
	NEW(smx,nBlksCMSA(cma) + 2,smx_typ);
	for(n=0,m=1; m<= nBlksCMSA(cma); m++){
           if(wt == 0) smx[m]=GetSmatrixFModel(pernats,ModelCMSA(m,cma));
	   else smx[m]=SampleSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
	   // else smx[m]=SampleDirichletSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
	   n+=LenSMatrix(smx[m]);
// PutSMatrix(stderr,smx[m]);
	}

#if 0	// ========== 3. Smooth gap function. ===========
	js_type	Spg= MkSpouge(nAlpha(A)+1,tFreqSeqSet(TrueDataCMSA(cma)),
			nBlksCMSA(cma),smx,MaxTrueSeqCMSA(cma)+1,
			observedGap,'g');
        gapfnct = GapScoreSpouge(Spg);
// fprintf(stderr,"max_gap_score = %d\n",MaxGapScoreSpouge(Spg));
// for(m=1; m<= nBlksCMSA(cma); m++) PutHistSpouge(stderr,m,5.0,Spg);
	NilSpouge(Spg); 
#endif

	// ========== 4. Sample a gapped alignment for sequence s. ===========
        E = gss->TrueSeq(s); 	
        oldE = gss->FakeSeq(s); 	
	// PutSeq(stderr,E,AlphabetCMSA(cma)); fflush(stderr);

#if 0	// ================== SBA PROGRAM. =======================
        Int4    *io,*ie,*od,*ed;
        NEW(io,n+2,Int4);NEW(ie,n+2,Int4);NEW(od,n+2,Int4);NEW(ed,n+2,Int4);
        for(i=1; i<=n; i++){
             io[i]=od[i]=-gapopen-gapextend; ie[i]=ed[i]=-gapextend;
        }
        sba_typ sba;
        sba.initialize(LenSeq(E),SeqPtr(E),io,ie,od,ed,nBlksCMSA(cma),
                        smx,gapfnct,pernats);
        // SAMPLE AT LOWER TEMPERATURES BY DEFAULT.
        operation=sba.TraceBack(&score,&trace_length,&start,
          pernats,1000.0/Temperature); // normal temperature -> pow(x,1.0);
	free(io);free(ie);free(od);free(ed);
#endif // ================== SBA PROGRAM. =======================
#if 1
	operation=GapAlnTraceSMatrix(gapopen,gapextend,LenSeq(E),
			SeqPtr(E),nBlksCMSA(cma),smx, NULL,&start);
#endif
#if 0
	operation=cma->mdl->GapAlnTrace(gapopen,gapextend,LenSeq(E),
		SeqPtr(E),gapfnct,&start);
	operation=GapAlnTraceSMatrix(gapopen,gapextend,LenSeq(E),
			SeqPtr(E),nBlksCMSA(cma),smx, gapfnct,&start);
#endif
	trace_length=strlen(operation);
#if 0
	PutSeqInfo2(stderr,E);
	PutGappedSeqAlnSMatrix(stderr, operation, start-1,
		LenSeq(E)-start+1, SeqPtr(E)+start-1, nBlksCMSA(cma), smx);
#endif
	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,E,newpos);
#if 0
if(s==968){
	fprintf(stderr,"sq=%d; start=%d; pos[1]=%d; operation=%s\n",s,start,newpos[1],operation);
	gsq->Put(stderr,A);
}
#endif
#if 0
if(s==22){
std::cerr << "newpos("; std::cerr << s; std::cerr << "):\n";
for(t=nBlksCMSA(cma); t > 0; t--){
	fprintf(stderr,"new: blk %d; site = %d (%d)\n",
		t,newpos[t],LenSeq(gsq->FakeSeq()));
PutSeq(stderr,gsq->FakeSeq(),A);
	fprintf(stderr,"old: blk %d; site = %d (%d)\n",
		t,oldpos[t],LenSeq(oldE));
PutSeq(stderr,oldE,A);
	// fprintf(stderr,"new: blk %d; site = %d\n",t,newpos[t]+OffSetSeq(E));
	// fprintf(stderr,"old: blk %d; site = %d\n",t,oldpos[t]+OffSetSeq(E));
}
std::cerr << "\n\n";
}
#endif

	// ========== 6. Deallocate memory. ===========
	for(m=1; m<= nBlksCMSA(cma); m++) NilSMatrix(smx[m]); free(smx);
        for(m=0; m<=nBlksCMSA(cma); m++) free(observedGap[m]); free(observedGap);
        free(operation); 

	// ========== 7. If sequence changed then replace s with gsq. =========
	if(gss->Identical(s,*gsq)){ 
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ(); NilCMSA(cma); 
		// for(t=1 ; t <= nBlksCMSA(cma); t++) AddSiteCMSA(t,s,oldpos[t], cma);
// std::cerr << "SampleGapsGibbsCMSA( ) failed\n\n";
		free(oldpos); free(newpos); delete []gsq; return FALSE; 
	}
	// gsq->Put(stderr, 60, A);
// std::cerr << "SampleGapsGibbsCMSA( ) succeeded!!!!\n\n";
	ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & fmodel.
	for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	free(oldpos); free(newpos);
//PutTypeSites(stderr, SitesCMSA(cma));
	newmap=RelMapCMSA(cma);
// fprintf(stderr,"oldmap = %.2f; newmap = %.2f\n",oldmap,newmap);
	// Now Sample one of these alignments.
	if(comp_map_cmsa(newmap, oldmap, 300./Temperature)){
		gss_typ     *gss0=gssCMSA(*oldcma);
		gss0->~gss_typ();
		NilCMSA(*oldcma); 
		*oldcma=cma; return TRUE;
		
	} else {
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ();
		NilCMSA(cma); return FALSE; 
	}
}

Int4	sites_from_operations_cmsa(char *operation, Int4 start, Int4 *pos,
	Int4 *indels,gsq_typ *gsq)
// get sites in fake sequence (i) from alignment trace for real sequence (j).
// 'start' is position in Real sequence for start of trace.
// gsq->FakeToReal(1) is position in 
{
   Int4   t,o,i,gaps=0,offset;

   if(gsq->FakeToReal(1) == 0){	offset=0; } // deletion at start of fakeSeq
   else { offset = start - gsq->FakeToReal(1); }
   for(t=i=0,o=1 ;operation[o] != 'E'; o++){
       switch(operation[o]){
            case 'M': i++; t++; pos[t]=i+offset; break;
	    case 'm': i++; break;  // matches.
            case 'D': i++; t++; pos[t]=i+offset; gaps++; break;
	    case 'd': i++; gaps++; break;       // deletions relative to profile.
            case 'i': i++; break;  // insertion between blocks.
            case 'I': while(operation[o]=='I'){ o++; gaps++; }
		    o--; break;		// insertion within a block.
            default: fprintf(stderr,"operations = %s\n",operation);
               print_error("sites_from_operations_cmsa( ): input error"); break;
       }  
   } *indels=gaps;
   return t;
}

double	GappedSimAnnealCMSA(cma_typ *M,double StartTemp,double EndTemp,double inc)
// sample with gaps using simulated annealing...
{
	cma_typ	cmsa=*M;
	Int4	n,iter,i,j,r,s,t,N=NumSeqsCMSA(cmsa);
	double	temp,TotLike,best,map,oldMap,newMap;
	dh_type dH;

	assert(cmsa != NULL); assert(StartTemp >= EndTemp);
	dH = dheap(N,3); fprintf(stderr,"."); 
	SetPseudoToMapCMSA(cmsa); best=map=RelMapCMSA(cmsa); 
	SaveBestCMSA(cmsa); // don't discard previous map for this one
	for(temp=StartTemp; temp >= EndTemp; temp-=5.0) {
	   for(n = 1; n <= NumColumnsCMSA(cmsa); n++){
		t=random_integer(nBlksCMSA(cmsa))+1; 
		if(nBlksCMSA(cmsa) <= 1) s = random_integer(3);
		else s = random_integer(6);
		switch(s){
		  case 0: if(ContigBlkCMSA(t,cmsa)) ShiftCMSA(cmsa,t);
		     MoveColumnCMSA(cmsa,t); break;
		  case 1: SampAddColTempCMSA(t,temp,cmsa); break;
		  case 2: SampRmColTempCMSA(t,temp,cmsa); break;
		  // more than one block 
		  case 3: TransferColumnCMSA(cmsa); break;
		  case 4: SampSlideColRtCMSA(t,&oldMap,&newMap,cmsa,temp); break;
		  default: SampSlideColLtCMSA(t,&oldMap,&newMap,cmsa,temp); break;
		}
	  }
       	  for(n = 1; n <= N; n++) insrtHeap(n,((keytyp)Random()),dH);
       	  while((n=delminHeap(dH)) != NULL){ 
	    SampleGapsGibbsCMSA(n,M,temp); 
            if((TotLike=RelMapCMSA(cmsa)) > best){ 
std::cerr << "map improves from"; std::cerr << best; std::cerr << " to "; std::cerr << TotLike; std::cerr << std::endl;
         	   best=TotLike; 
         	   if(map < best){ map=best; SaveBestCMSA(cmsa); }
            }
       	  }
	}
	fprintf(stderr,"final temperature = %2f K\n %.1f\n", temp,map);
	Nildheap(dH); InitMAPCMSA(cmsa);
	return map;
}

//======================== NEW ROUTINES ==========================

Int4	slide_left_dfs_cmsa(Int4 t, double *new_map, cma_typ cma, Int4 depth)
// depth-first-search for the lpr of slided alignments.
// WARNING: Assumes that columns are not fragmented.
{ 
	Int4		t1,t2,s,ntyp,i;

	t1=t; t2=t+1; // for slide left: source = t2, denstination = t1.
	if(LengthCMSA(t2,cma) <= 3) return depth; else depth++;
// if(NullSiteCMSA(t,LengthCMSA(t2,cma)-1,cma)) return depth;

	RmColumnMSA(t2,1,cma); // Remove first column from source block and...
	i=LengthCMSA(t1,cma); AddColumnMSA(t1,i+1,cma);  // add to target block.
	new_map[depth]= RelMapCMSA(cma);
	return slide_left_dfs_cmsa(t,new_map,cma,depth);
}

BooLean	dfsSlideColLtCMSA(Int4 t, double *oldMap, double *newMap,
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

      t1         t2
   --[******]---[**:**********]--
               |		This operation allows escape from some 
               V                 nasty traps in alignment space.
   --[******:**]---[**********]--

  returns TRUE if succeeds in sampling ends, otherwise returns FALSE.
 *************************************************************************/
{ 
    Int4		t1,t2,n,s,ntyp,depth,d,i,j;
    double		map0,rand_no,*prob,total,ratio;
    double		*new_map;
    
    ntyp=nBlksCMSA(L);
    if(ntyp < 2) return FALSE; else if(t < 1 || t >= ntyp) return FALSE;
    t1=t; t2 = t1+1; // slide left: source = t2, denstination = t1.
    if(!(ContigBlkCMSA(t1,L) && ContigBlkCMSA(t2,L))) return FALSE;
    if(LengthCMSA(t2,L) <= 3) return FALSE;

    map0 = RelMapCMSA(L); // 1.b. Determine the old alignment map.
    *oldMap=map0;
    new_map = new double [LengthCMSA(t2,L) +2];
    depth=slide_left_dfs_cmsa(t,new_map,L,0);
    assert(depth > 0);

	prob = new double [depth+3]; 
	for(total=1.0,d=1; d <= depth; d++){
		ratio = exp(new_map[d]-map0);	// ratio of new to old.
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		prob[d] = ratio; total += ratio; 
	} // 5. Sample an alignment.
	rand_no = total*SampleUniformProb();
	for(d=depth; d > 0; d--,i++){
		if((rand_no -= prob[d]) <= 0.0){
			*newMap=new_map[d]; delete [] new_map;
// fprintf(stderr,"%d columns moved from block %d to %d\n",d,t2,t1);
			delete [] prob; return TRUE;
		} else {  // move columns back where they were...
			i=LengthCMSA(t1,L); RmColumnMSA(t1,i,L);
			AddColumnMSA(t2,0,L);
	        }
	}
	delete [] new_map; delete [] prob; return FALSE;
}

Int4	slide_right_dfs_cmsa(Int4 t, double *new_map, cma_typ cma, Int4 depth)
// depth-first-search for the lpr of slided alignments.
// WARNING: Assumes that columns are not fragmented and nBlksCMSA(cma) > 1.
{ 
	Int4		t1,t2,s,ntyp,i;

	t1=t; t2=t+1; // for slide left: source = t1, denstination = t2.
	if(LengthCMSA(t1,cma) <= 3) return depth; else depth++;
	i=LengthCMSA(t1,cma);
	RmColumnMSA(t1,i,cma);  // Remove last column from source block and...
	AddColumnMSA(t2,0,cma); // add to target block.
	new_map[depth]= RelMapCMSA(cma);
	return slide_right_dfs_cmsa(t,new_map,cma,depth);
}


BooLean	dfsSlideColRtCMSA(Int4 t, double *oldMap, double *newMap, 
	cma_typ L, double temperature)
/*************************************************************************
  slide the end of one model onto the next model.

   --[******:**]--[**********]--
               |		This operation allows escape from some 
               V                 nasty traps in alignment space.
   --[******]--[**:**********]--

  returns TRUE if succeeds in sampling ends, otherwise returns FALSE.
 *************************************************************************/
{ 
    Int4		t1,t2,n,s,ntyp,depth,d,i,j;
    double		map0,map,rand_no,*prob,total,ratio;
    double		*new_map;
    
    ntyp=nBlksCMSA(L);
    if(ntyp < 2) return FALSE; else if(t < 1 || t >= ntyp) return FALSE;
    t1=t; t2 = t1+1; // slide left: source = t1, denstination = t2.
    if(!(ContigBlkCMSA(t1,L) && ContigBlkCMSA(t2,L))) return FALSE;
    if(LengthCMSA(t1,L) <= 3) return FALSE;

    map0 = RelMapCMSA(L); // 1.b. Determine the old alignment map.
    *oldMap=map0;
    new_map = new double [LengthCMSA(t1,L)+2];
    depth=slide_right_dfs_cmsa(t,new_map,L,0);
    assert(depth > 0);

	NEW(prob,depth+3,double);
	for(total=1.0,d=1; d <= depth; d++){
		ratio = exp(new_map[d]-map0);	// ratio of new to old.
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		prob[d] = ratio; total += ratio; 
	} // 5. Sample an alignment.
	rand_no = total*SampleUniformProb();
	for(d=depth; d > 0; d--,i++){
		if((rand_no -= prob[d]) <= 0.0){
// fprintf(stderr,"%d columns moved from block %d to %d\n",d,t1,t2);
			*newMap=new_map[d]; delete [] new_map;
			free(prob); return TRUE;
		} else {  // move column back...
			RmColumnMSA(t2,1,L);
			i=LengthCMSA(t1,L); AddColumnMSA(t1,i+1,L);
	        }
	}
	delete [] new_map; free(prob); return FALSE;
}

/****************** Squeeze in Tight colums ************************/

#if 0
BooLean	SampInsertColCMSA(Int4 t, double temperature, cma_typ L)
/*************************************************************************
  Samples either a right inserted column or a left inserted column 
  (returns TRUE) or no inserted colum (returns FALSE) for block t.
/*************************************************************************/
{
	Int4	New,d;
	double	omap,nmap;

	omap=RelMapCMSA(L); New=SampNewColMSA(t,L);
	if(New==1) return FALSE;	/** no room to add a column **/
	d=AddColumnMSA(t, New, L);	/** d = displacement of model **/
	nmap=RelMapCMSA(L);
	if(comp_map_cmsa(nmap, omap,temperature)) return TRUE;
	RmColumnMSA(t, New-d, L);
	return FALSE;
}

BooLean	SampleEndColCMSA(Int4 t, double temperature, cma_typ L)
/*******************************************************************
   probability ratio for RatioFModel(site_freq[i],1, M);
  Sample a new column (to add to model); return the new column sampled.

  1. Get frequencies of columns (use 'X' residues for blocked sites).
   1a. At the same time get positions in seqs. of new column (pos[1..NumSeq]).
   1b. Find out how many and which of these sites are blocked by 
	calling the sites function: OpenPos(n,p,S) for each sequence.

Int4    *AlwaysGetSiteFreq(st_type S,Int4 t,Int4 d,Int4 *pos,char *blocked);

  4. Find out gap penalties for an insertion at each blocked site.
	Use new function call: gss->GapCost(n,p,gap) for this
	which calls gsq.(n,p,gap).

Int4    GapCost(Int4,UInt4,unsigned short);


  5. Get total likelihoods for each of these columns and for
	no new column at all and sample one of these three options.

p = RatioFModel(site_freq[i],1, M);
  
  4. For those sites that are not open, call another new sites
	function InsertGapSites(Int4 n, Int4 pos, cma_typ L);
	which calls gss->InsertGap(n,p,gap);

void    InsertGap(Int4, UInt4,unsigned short);

  
  NEED to know positions of columns in 

BooLean OccupiedSite(register Int4 t, register Int4 n, register Int4 site,
        register st_type S);
void    ReplaceSeqSites(Int4 n, gsq_typ *gsq, st_type S);
********************************************************************/
#endif

BooLean	SampleEndColCMSA(Int4 t, double temperature, cma_typ cma)
{
	st_type S=SitesCMSA(cma);
	Int4	lemon,*site_freq[3],ncol,leng,*location;
	double	ratio[3],total,rand_no,weight;
	fm_type M=ModelCMSA(t,cma);
	char	*blocked;

	lemon = LemonFModel(M); /** proportional to prob **/

if(lemon == INT4_MIN) return FALSE;  // temporary fix for overflow.

	location = new Int4[NumSeqsCMSA(cma) + 2];
	blocked = new char[NumSeqsCMSA(cma) + 2];
	ncol = nColsFModel(M); leng = LenFModel(M); 
	total=1.0;	// add one for no sampling at all.
	site_freq[1]=AlwaysGetSiteFreq(S,t,FALSE, location,blocked);
	ratio[1] = RatioFModel(site_freq[1],lemon, M);
	weight = add_col_weight(ncol, leng, 0);
	ratio[1] *= weight; total += ratio[1];

	site_freq[2]=AlwaysGetSiteFreq(S,t,TRUE,location,blocked);
	ratio[2] = RatioFModel(site_freq[2],lemon, M);
	weight = add_col_weight(ncol, leng, leng+1);
	ratio[2] *= weight; total += ratio[2];

// WARNING: NEED TO ADD gap penalties...
// CALL: Int4    gss.GapCost(Int4,UInt4,unsigned short);

	delete [] location; delete [] blocked;
	free(site_freq[1]); free(site_freq[2]);
	rand_no = total*SampleUniformProb();;
	if((rand_no -= ratio[1]) <= 0.0) {	// add left column.
// std::cerr << "Sampled insertion of left column...\n";
		InsertColCMSA(t, FALSE, cma); return TRUE;
	} 
	if((rand_no -= ratio[2]) <= 0.0) {	// add right column.
// std::cerr << "Sampled insertion of right column...\n";
		InsertColCMSA(t, TRUE, cma); return TRUE;
	}
	return FALSE;
}

BooLean	ReAlignGSqCMSA(Int4 s,char *operation, Int4 start, cma_typ *oldcma)
// realign sequence s using operation string.
{
	Int4	m,n,i,t,trace_length,gapopen,gapextend;
	Int4    score,*newpos,*oldpos;
	a_type	A=AlphabetCMSA(*oldcma);
	e_type	oldE,E;
	gss_typ	*gss=gssCMSA(*oldcma);
	double	pernats;
	cma_typ	cma;
	double	oldmap,newmap;

	cma=CopyCMSA(*oldcma);
	gapopen=gss->GapOpen(); gapextend=gss->GapExtend(); pernats=gss->PerNats();

	// 1. Remove sequence s from alignment.
	NEW(oldpos,nBlksCMSA(cma)+2,Int4);  // save old sites.
	NEW(newpos,nBlksCMSA(cma)+2,Int4);  // for new sites.
	for(t=1; t<=nBlksCMSA(cma); t++){
		PosSiteCMSA(t, s, cma->pos, cma); oldpos[t]=cma->pos[1];
	} VacateSitesCMSA(s,cma);

	// ========== 4. Sample a gapped alignment for sequence s. ===========
        E = gss->TrueSeq(s); 	
        oldE = gss->FakeSeq(s); 	
	// PutSeq(stderr,E,AlphabetCMSA(cma)); fflush(stderr);
	trace_length=strlen(operation);
	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,E,newpos);
	// ========== 7. If sequence changed then replace s with gsq. =========
	if(gss->Identical(s,*gsq)){  // This should be FALSE with new calling routine...
		// gsq->Put(stderr,A);
		gss_typ     *gss0=gssCMSA(cma); gss0->~gss_typ(); NilCMSA(cma); 
		free(oldpos); free(newpos); delete []gsq; return FALSE; 
	}
	// gsq->Put(stderr, 60, A);
	ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & in fmodel.
	for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	free(oldpos); free(newpos);
#if 0
//PutTypeSites(stderr, SitesCMSA(cma));
	oldmap=RelMapCMSA(*oldcma);
	newmap=RelMapCMSA(cma);
// fprintf(stderr,"oldmap = %.2f; newmap = %.2f\n",oldmap,newmap);
#endif
	gss_typ     *gss0=gssCMSA(*oldcma);
	gss0->~gss_typ();
	NilCMSA(*oldcma); 
	*oldcma=cma; return TRUE;
}

BooLean	ClusterSampleGapsGibbsCMSA(Int4 *sqset,cma_typ *oldcma,double Temperature)
// sample insertions and deletions in sequence s.
// Int4 *sqset is an array of sequence ids and ends with 0.
// e.g. [23,45,98,0] --> Nsqset = 3.
{
	Int4	m,n,i,t,start,trace_length,indels,gapopen,gapextend;
	Int4	**newpos,**oldpos,s,Nsqset;
	smx_typ	*smx;
	a_type	A=AlphabetCMSA(*oldcma);
	e_type	oldE,E;
	gss_typ	*gss=gssCMSA(*oldcma);
	double	pernats;
	char	*operation;
	Int4	wt;
	cma_typ	cma;
	double	oldmap,newmap;

	cma=CopyCMSA(*oldcma);
	oldmap=RelMapCMSA(*oldcma);
	gapopen=gss->GapOpen(); gapextend=gss->GapExtend(); pernats=gss->PerNats();

        if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
        else if(Temperature > 275.0) wt=2;
        else if(Temperature > 250.0) wt=3;
        else if(Temperature > 225.0) wt=4;
        else if(Temperature > 200.0) wt=5;
        else if(Temperature > 175.0) wt=6;
        else if(Temperature > 150.0) wt=7;
        else if(Temperature > 125.0) wt=8;
        else if(Temperature > 100.0) wt=9;
        else if(Temperature > 75.0) wt=10;
        else wt=0;  // Use maximum likelihood.

//**************************************************************
	// 1. Remove sequences in sqset from alignment.
  	for(Nsqset=i=0; (s=sqset[i]); i++) Nsqset++;
	NEWP(oldpos,Nsqset+2,Int4);  // save old sites.
	NEWP(newpos,Nsqset+2,Int4);  // for new sites.
	for(i=0; (s=sqset[i]); i++){
	  NEW(oldpos[i],nBlksCMSA(cma)+2,Int4);  // save old sites.
	  NEW(newpos[i],nBlksCMSA(cma)+2,Int4);  // for new sites.
	  for(t=1; t<=nBlksCMSA(cma); t++){
		PosSiteCMSA(t, s, cma->pos, cma); oldpos[i][t]=cma->pos[1];
	  } VacateSitesCMSA(s,cma);
	}
//**************************************************************

	// ===== 2. Obtain scoring matrix from rest of alignment.
	NEW(smx,nBlksCMSA(cma) + 2,smx_typ);
	for(n=0,m=1; m<= nBlksCMSA(cma); m++){
           if(wt == 0) smx[m]=GetSmatrixFModel(pernats,ModelCMSA(m,cma));
	   else smx[m]=SampleSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
	   // else smx[m]=SampleDirichletSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
	   n+=LenSMatrix(smx[m]);
// PutSMatrix(stderr,smx[m]);
	}

//**************************************************************
      BooLean	AllIdentical=TRUE;
      for(i=0; (s=sqset[i]); i++){
	// ========== 4. Sample a gapped alignment for sequence s. ===========
        E = gss->TrueSeq(s); 	
        oldE = gss->FakeSeq(s); 	
	// PutSeq(stderr,E,AlphabetCMSA(cma)); fflush(stderr);

	operation=GapAlnTraceSMatrix(gapopen,gapextend,LenSeq(E),
			SeqPtr(E),nBlksCMSA(cma),smx, NULL,&start);
	trace_length=strlen(operation);

	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,E,newpos[i]);

	// ========== 7. If sequence changed then replace s with gsq. =========
	if(gss->Identical(s,*gsq)){ delete []gsq; }
	else {
		AllIdentical=FALSE;
		// gsq->Put(stderr, 60, A);
		// std::cerr << "ClusterSampleGapsGibbsCMSA( ) succeeded!!!!\n\n";
		ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & fmodel.
		for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[i][t], cma);
	} free(oldpos[i]); free(newpos[i]); free(operation); 
      } free(oldpos); free(newpos); 
//**************************************************************

	// ========== 6. Deallocate memory. ===========
	for(m=1; m<= nBlksCMSA(cma); m++) NilSMatrix(smx[m]); free(smx);

	if(AllIdentical){	// then return without doing anything.
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ();
		NilCMSA(cma); return FALSE; 
	}

	//PutTypeSites(stderr, SitesCMSA(cma));
	newmap=RelMapCMSA(cma);
	// fprintf(stderr,"oldmap = %.2f; newmap = %.2f\n",oldmap,newmap);
	// Sample either old or new alignment.
	if(comp_map_cmsa(newmap, oldmap, 300./Temperature)){
		gss_typ     *gss0=gssCMSA(*oldcma);
		gss0->~gss_typ();
		NilCMSA(*oldcma); 
		*oldcma=cma; return TRUE;
	} else {
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ();
		NilCMSA(cma); return FALSE; 
	}
}



