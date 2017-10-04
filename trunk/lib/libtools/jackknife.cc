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

#include "jackknife.h"

Int4	Jackknife(ss_type data, ss_type msadata, sma_typ MA, sma_typ *MA2, 
	ss_type *data2)
{	
	double	cutoff=0.01;
	double	Cluster_cutoff=0.01;
	Int4	dbs_size=300000,MaxBadBlks=0;
	char	method='e'; 
	return Jackknife(method, Cluster_cutoff,cutoff,dbs_size,data, msadata, MA, MA2,
		data2,MaxBadBlks);
}

Int4	Jackknife(char method,double Cluster_cutoff, double cutoff, Int4 dbs_size, 
	ss_type data, ss_type msadata, sma_typ MA, sma_typ *MA2, ss_type *data2,
	Int4 MaxBadBlks)
{
	sma_typ	ma1,ma2;
	e_type	*ListE,**Esets,dummy[3];
	Int4	Nsets,*counts,*cnts1,*cnts2,i,j,n,s,g,N,maxseq;
	a_type	A=SeqSetA(data);
	Int4	*setsize,maxset,max;
	unsigned short  *nsize;
	char	str[220];
	snh_typ sH;
	BooLean *select,*rm_grp,*fullselect;
	FILE	*fp;
	double	adj_cutoff;

	maxseq = MaxSeqSeqSet(data);
	i = MaxSeqSeqSet(msadata);
	maxseq = MAXIMUM(Int4,i,maxseq);
	/************************ MAIN ALGORITHM ********************************/
	//  1. Merge msafa and srch files removing identities from srch_file.
	ListE = MergeSeqSet(&N, msadata, data);
	fprintf(stdout,"  cutoff = %g; database size = %d\n",Cluster_cutoff,dbs_size);
	fprintf(stdout,"Adjusted cutoff = %g\n",Cluster_cutoff/(double)dbs_size);
	adj_cutoff= -log10(Cluster_cutoff * (double)N/(double)dbs_size);
	// multiply by N because this is divided out again within ClusterSeqs( ) routine.

	//  2. If sequences need to ignore motifs (e.g., Walker A or B) then  
     	//     mask out these regions. (wait on this...) 

	//  3. Cluster these sequences into distinct transitively-related groups. 
	cnts1 = CntsSeqSet(data);
	cnts2 = CntsSeqSet(msadata);
	NEW(counts, nAlpha(A) +2, Int4);
	for(i = 0 ; i<= nAlpha(A); i++) counts[i] = cnts1[i] + cnts2[i];
#if 1
	Esets = ClusterSeqs(ListE,method, adj_cutoff, counts, A, &Nsets);
#else
        if(isupper(method)){
          ListE1 = ClusterGPSI(method,adj_cutoff,P,A,&Nset,T);
        } else Esets =ClusterSeqs(ListE,method,adj_cutoff,CntsSeqSet(P),A,&Nset);
#endif
	NEW(setsize, Nsets +2, Int4);
	for(max=0,g = 1; g <= Nsets; g++){
		fprintf(stdout,"group %d:\n",g);
		for(s=1; Esets[g][s] != NULL; s++){
			PutSeqInfo(stdout,Esets[g][s]);
		}
		setsize[g] = s-1;
		if(max < setsize[g]){ max = setsize[g]; maxset = g; }
		fprintf(stdout,"\n");
	}

/**  4. Remove sequences not in main group from the alignment. **/
	select = IntersectSeqs(Esets[maxset],msadata,&n); 
	ma1 = RmSMA(select, n, MA);
	assert(strlen(NameSeqSet(msadata)) < 200);
	sprintf(str,"%s.tmp.msa",NameSeqSet(msadata));
	fp = open_file(str,"","w"); PutSMA(fp,ma1); fclose(fp);
	free(select); 

/**  5. For each group: **/
	NEW(rm_grp, Nsets + 2, BooLean);
#if 0
	osn_typ F;
	F=MakeOScan(0.1,str,counts,A,1,'M','G',cutoff,1.0,maxseq,0,FALSE);
	NoMaskOScan(F);
	for(g = 1; g <= Nsets; g++){
	 if(g != maxset){
	  fprintf(stdout,"group %d out of %d\n",g,Nsets);

	  /** b. Scan other groups against the main alignment model. **/
	  fullselect = IntersectSeqs(Esets[g],data,&n);
	  NEW(nsize, NSeqsSeqSet(data)+2, unsigned short);
	  fp = tmpfile();
	  for(i=0,s=1; s<=NSeqsSeqSet(data); s++){
	    if(fullselect[s]){
		PutSeqSetE(fp,s, data);
		i++; nsize[i] = LenSeq(SeqSetE(s,data));
	    }
	  } rewind(fp);
	  adj_cutoff = cutoff*(double)setsize[g]/dbs_size;
	  fprintf(stdout,"Group %d oscan adjusted cutoff = %g\n",g,adj_cutoff);
	  SetEvalOScan(adj_cutoff,1.0,F);
	  sH=OScanScan(fp,i, nsize, F);
	  fclose(fp);

	  /** c. Remove groups that fail the scan cutoff score.  **/
	  if(nScanHeap(sH) <= 0) rm_grp[g] = TRUE;
	  else {
		fprintf(stdout,"    %d hits\n",nScanHeap(sH));
		PutInfoScanHeap(stdout, sH, A);
	  }
	  NilScanHeap(sH); 
	  free(nsize); 
	  free(fullselect);
	 }
	}
	NilOScan(F);
#else
	double	*freq;
	// BooLean	weights=FALSE;
	BooLean	weights=TRUE;
	Int4	total;
        NEW(freq,nAlpha(A)+2,double);
	for(total=0, i=0; i<=nAlpha(A); i++) total += counts[i];
        for(i=0; i<=nAlpha(A); i++) freq[i]= (double) counts[i]/(double)total;

	gsn_typ F;
	// F=MakeGOScan(0.1,str,counts,A,1,'M','G',cutoff,1.0,maxseq,0,FALSE);
	// F=MakeGOScan(argv[2],A,maxrpts,method,mode,expect,singleEval,m,weights,
        //         minmap,pseudo,freq);
	F=MakeGOScan(str,A,1,'H','O',cutoff,1.0,maxseq,weights,0.0,0.1,freq);
	NoMaskGOScan(F);
	for(g = 1; g <= Nsets; g++){
	 if(g != maxset){
	  fprintf(stdout,"group %d out of %d\n",g,Nsets);

	  /** b. Scan other groups against the main alignment model. **/
	  fullselect = IntersectSeqs(Esets[g],data,&n);
	  NEW(nsize, NSeqsSeqSet(data)+2, unsigned short);
	  fp = tmpfile();
	  for(total=i=0,s=1; s<=NSeqsSeqSet(data); s++){
	    if(fullselect[s]){
		PutSeqSetE(fp,s, data);
		i++; nsize[i] = LenSeq(SeqSetE(s,data));
		total += nsize[i];
	    }
	  } rewind(fp);
	  adj_cutoff = cutoff*(double)setsize[g]/dbs_size;
	  fprintf(stdout,"Group %d oscan adjusted cutoff = %g\n",g,adj_cutoff);
	  SetEvalGOScan(adj_cutoff,1.0,F);
#if 1
	  // sH=GOScanScan(fp,i, nsize, F);
	  sH=GOScanScan(fp,i,total,nsize,1,MaxBadBlks,F);
#else
	  sH=GOScanScan(fp,i,total,nsize,F);
snh_typ GOScanScan(fp,i,Int4 total,unsigned short *nsize,
                Int4 min_rpt, gsn_typ F);
snh_typ GOScanScan(fp,i,Int4 total,unsigned short *nsize,
                Int4 min_rpt, Int4 MaxBadBlks, gsn_typ F);
#endif
	  fclose(fp);

	  /** c. Remove groups that fail the scan cutoff score.  **/
	  if(nScanHeap(sH) <= 0) rm_grp[g] = TRUE;
	  else {
		fprintf(stdout,"    %d hits\n",nScanHeap(sH));
		PutInfoScanHeap(stdout, sH, A);
	  }
	  NilScanHeap(sH); 
	  free(nsize); 
	  free(fullselect);
	 }
	}
	NilGOScan(F);
#endif
	NilSMA(ma1);
fprintf(stderr,"done with jackknife test.\n");

/**  6. return the retained seqset aInt4 with a corresponding msa model for that set. **/
	dummy[1]=NULL;
	select = IntersectSeqs(ListE,msadata,&n); 
	ma1 = RmSMA(select, n, MA);
#if 0
	PutSMA(stderr,ma1);
#endif
	free(select);
	for(j=0,g = 1; g <= Nsets; g++){
	   if(rm_grp[g]){
fprintf(stdout,"set %d failed.\n",g);
		select = DifferSeqs(Esets[g],msadata,&n); 
		ma2 = RmSMA(select, n, ma1);
		NilSMA(ma1); ma1 = ma2;
		for(i = NSeqsSeqSet(msadata); i > 0; i--){
		  if(!select[i]) {
			RemoveSeqSet(i, msadata);
		  }
	        }
		free(select); 
#if 1
		sprintf(str,"%s.%dsq",NameSeqSet(msadata),g);
		FILE *tfp= open_file(str,"","w"); 
		for(i=1; Esets[g][i] != NULL; i++) {
			PutSeq(tfp,Esets[g][i],A);
			NilSeq(Esets[g][i]);
		} fclose(tfp);
#else
		for(i=1; Esets[g][i] != NULL; i++) NilSeq(Esets[g][i]);
#endif
	   } else {
fprintf(stdout,"set %d retained.\n",g);
		for(i=1; Esets[g][i] != NULL; i++){
		   j++; ListE[j]=Esets[g][i]; 
		}
	   }
	   free(Esets[g]);
	}
	N = j;
fprintf(stderr,"creating new alignment model...\n");
	*MA2 = ma1;
	*data2 = Array2SeqSet(ListE, N, NameSeqSet(data),A);
	free(counts);
	free(rm_grp);
	free(Esets);
	free(setsize);
#if 0
/** later loop over this function recursively if merged some but not all groups. **/
/*************************************************************************
	
 *************************************************************************/
     if(some_merged && Nset > 1){
	N = Jackknife(data, msadata, MA, &MA2, &data2);
     } else return Nset;
#endif
	return N;
}

BooLean	*DifferSeqs(e_type *ListE, ss_type data, Int4 *Nnew)
{
	Int4	N,n,i,j;
	BooLean	*select;
	e_type	E;

	N = NSeqsSeqSet(data);
	NEW(select, N+2, BooLean);
	for(n=0,i=1; i <= N; i++){
	   E = SeqSetE(i,data);
	   select[i] = TRUE;
	   for(j=1; ListE[j] != NULL; j++){ 
		if(IdentSeqs(E,ListE[j])) { select[i] = FALSE; break; }
	   }
	   if(select[i]) n++;
	}
	*Nnew = n;
	return select;
}

BooLean	*IntersectSeqs(e_type *ListE, ss_type data, Int4 *Nnew)
{
	Int4	N,n,i,j;
	BooLean	*select;
	e_type	E;

	N = NSeqsSeqSet(data);
	NEW(select, N+2, BooLean);
	for(n=0,i=1; i <= N; i++){
	   E = SeqSetE(i,data);
	   for(j=1; ListE[j] != NULL; j++){ 
		if(IdentSeqs(E,ListE[j])) { select[i] = TRUE; break; }
	   }
	   if(select[i]) n++;
	}
	*Nnew = n;
	return select;
}

e_type	*MergeSeqSet(Int4 *N, ss_type data1, ss_type data2)
/** merge sequences in both sets, preferentially retaining sequences in data1 **/
{
	Int4	N1,N2,i,j;
	e_type  *EList,E;

	N1 = NSeqsSeqSet(data1);
	N2 = NSeqsSeqSet(data2);
	NEW(EList, N1 + N2 +2, e_type);
	for(i=0; i < N1; ){
		i++; E = SeqSetE(i,data1); EList[i] = CopySeq(E);
	}
	for(j=1; j <= N2; j++) {
		i++; E = SeqSetE(j,data2); EList[i] = CopySeq(E);
	}
	*N = RmIdent(EList);
	return EList;
}


