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

#include "scanheap.h"

snh_typ  MakeScanHeap(Int4 hpsz,Int4 nblk, Int4 *len, double Ecut)
{ return MakeScanHeap(hpsz,nblk, len,Ecut,0); }

snh_typ  MakeScanHeap(Int4 hpsz,Int4 nblk, Int4 *len, double Ecut,
	Int4 MaxNotOK)
{
	snh_typ	H;
	Int4	i;
	char	str[40];

	NEW(H,1,scanheap_type);
	NEW(H->info,hpsz+2,sni_typ);
	H->mH=Mheap(hpsz, 3);
	H->hpsz=hpsz;
	H->nblk=nblk;
	H->HG = Histogram("-log10(E-values)", -100, 200,2.0);
	H->Ecut=Ecut;
	H->use_cutoff=-100.0;
	// H->use_cutoff=1.3; // old cutoff
	H->neglog10Ecut=-log10(Ecut);
	NEW(H->hg,nblk+2,h_type);
	NEW(H->len,nblk+2,Int4);
	for(i=1; i<=nblk; i++){
	   H->len[i]=len[i];
	   sprintf(str,"block %d -log10(E-values)",i);
	   H->hg[i]=Histogram(str, -100, 200,2.0);
	}
	H->not_okay=NULL; H->sort=NULL;
	H->max_not_ok=MaxNotOK;
	H->number=-1;
	return H;
}

void    NilScanHeap(snh_typ H)
{
	Int4	i;

	process_scanheap(H);
	for(i=1; i<=H->number; i++) NilScanInfo(H->sort[i]);
	NilMheap(H->mH); free(H->info); free(H->len);
	NilHist(H->HG); free(H->not_okay); free(H->sort); 
	for(i=1; i<=H->nblk; i++) {
		NilHist(H->hg[i]); 
	}
	free(H->hg);
	free(H->use);
	free(H);
}

sni_typ	*RtnNilScanHeap(snh_typ H)
/** destroy scanheap & return sorted list of hits: from 1..number (end=NULL) **/
{
	Int4	i,j,n;
	sni_typ	*sort;

	process_scanheap(H);
	for(n=0,i=1; i<=H->number; i++) if(H->not_okay[i] <= H->max_not_ok) n++; 
	NEW(sort, n+2, sni_typ);
	for(j=0,i=1; i<=H->number; i++){
		if(H->not_okay[i] <= H->max_not_ok){ j++; sort[j] = H->sort[i]; } 
		else NilScanInfo(H->sort[i]);
	}
	j++; sort[j] = NULL;
	free(H->sort); free(H->info); free(H->len);
	NilMheap(H->mH); NilHist(H->HG); free(H->not_okay); 
	for(i=1; i<=H->nblk; i++) NilHist(H->hg[i]); free(H->hg);
	free(H->use); free(H);
	return sort;
}

Int4	InsertScanHeap(sni_typ I, snh_typ H)
{
	Int4	item,i;

	double x=ScoreScanInfo(I);
	if(!isfinite(x)){	// this is not a valid floating point number.
	  if(isnan(x)) print_error("InsertScanHeap() error: NaN"); 
	  else if(isinf(x)){
		  if(isinf(x) == 1) print_error("InsertScanHeap() error: pos. infinite"); 
		  else if(isinf(x) == -1){ 	// negative infinity
			x = -DBL_MAX;
	  		IncdHist(-x,H->HG); item=InsertMheap(x, H->mH);
			// print_error("InsertScanHeap() error: neg. infinite"); 
		  } else print_error("InsertScanHeap() error: score is infinite");
	  } else print_error("InsertScanHeap() error: score input error");
	} else {		// insert score as is.
	  IncdHist(-ScoreScanInfo(I),H->HG);
	  item=InsertMheap(ScoreScanInfo(I), H->mH);
	}
	if(item==NULL) return NULL;
	if(H->info[item]!= NULL) NilScanInfo(H->info[item]);
	H->info[item]=I; 
	for(i=1; i<=H->nblk; i++){ 
		IncdHist(scoreScanInfo(i,I),H->hg[i]); 
	}
	return item;
}

sni_typ	DelMinScanHeap(snh_typ H)
{
	Int4	item;
	sni_typ	I;

	item = DelMinMheap(H->mH);
	if(item==NULL) return NULL;
	I = H->info[item]; H->info[item]=NULL;
	return I;
}

sni_typ	DelMaxScanHeap(snh_typ H)
{
	Int4	item;
	sni_typ	I;

	item = DelMaxMheap(H->mH);
	if(item==NULL) return NULL;
	I = H->info[item]; H->info[item]=NULL;
	return I;
}

sni_typ *InfoScanHeap(Int4 *n, snh_typ H)
{
	Int4	i,N;
	sni_typ *info,I;

	process_scanheap(H);
	NEW(info,H->number+2,sni_typ);
	for(N=0,i=1; i<=H->number; i++){
	  I=H->sort[i];
	  if(H->not_okay[i] <= H->max_not_ok){ N++; info[N]=I; }
	}
	*n=N;
	return info;
}

sni_typ	*DomainsScanHeap(Int4 *n, Int4 Nt_flank, Int4 Ct_flank, snh_typ H)
/*** return a subsequence corresponding to the domain delineated by I ***/
{
	Int4	s,i,j,j0,m,m2,m0,N,N0,nblks,gap,s2,rpts;
	Int4	start,end,r,len,total;
	float	p;
	sni_typ	I,*listI;
	BooLean	okay=TRUE;

	process_scanheap(H);
	N0 = nblks=nblksScanHeap(H);
	for(total=0,i=1; i<=H->number; i++){
	   I=H->sort[i];
	   total += RptsScanInfo(I);
	}
	NEW(listI,total+2,sni_typ);
	for(N=0,i=1; i<=H->number; i++){
	  I=H->sort[i];
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
	  if((rpts*N0) != nblks) print_error("DomainsScanHeap( ) input error");

         if(rpts > 1 || (H->not_okay[i] <= H->max_not_ok)){
	  for(r=1, m=1; r<=rpts; r++, m+=N0){
		start=SiteScanInfo(m,I);
		m2 = m + N0 - 1;
		s=SiteScanInfo(m2,I);
		m0 = ((m2-1) % N0) + 1;
		end = m+N0;
		for(okay=TRUE, j=m; j < end; j++){
		   p = scoreScanInfo(j,I);
		   j0 = ((j-1) % N0) + 1;
		   if(H->use[j0] && p < H->neglog10Ecut){ okay=FALSE; break; }
		}
		if(okay){
		  N++;
		  listI[N] = MakeRptScanInfo(r,Nt_flank,H->len[m0]+Ct_flank-1,I);
		}
	  }
	 }
	}
	*n=N;
	return listI;
}

e_type  *DomainsScanInfo2(Int4 *n, Int4 Nt_flank, Int4 Ct_flank, snh_typ H)
/*** return a subsequence corresponding to the domain delineated by I ***/
{
        e_type  *E,E1;
        Int4    i,s,e,N,nblk,lastblklen;
	sni_typ	I;

	process_scanheap(H);
	nblk = nblksScanHeap(H);
	lastblklen = BlkLenScanHeap(nblk,H);
	NEW(E,H->number+2,e_type);
	for(N=0,i=1; i<=H->number; i++){
	  if(H->not_okay[i] <= H->max_not_ok){ 
	  	I=H->sort[i];
        	E1 = SeqScanInfo(I);
        	s = SiteScanInfo(1,I);
        	s = MAXIMUM(Int4, s - Nt_flank, 1);
        	e = SiteScanInfo(nblk,I) + lastblklen - 1;
		e = MINIMUM(Int4, e + Ct_flank, LenSeq(E1));
		N++; 
		E[N] = MkSubSeq(s,e,E1); 
	  }
	}
	*n=N;
	return E;
}

double	*WeightsScanHeap(snh_typ H, a_type A)
{
	sni_typ *info,I;
        Int4    i,j,k,r,s,m,n,nblks,totlen;
        double  w,max,min,*wt,N,d;
        e_type  E;
        short   **nres,*ntyp;
        unsigned char    *seq;
	h_type	HG;

	process_scanheap(H);
	NEW(info,H->number+1,sni_typ);
	for(n=0,i=1; i<=H->number; i++){
	  I=H->sort[i];
	  if(H->not_okay[i] <= H->max_not_ok){ n++; info[n]=I; }
	}
        nblks=nblksScanHeap(H);
        for(totlen=0,m=1; m<=nblks; m++){ totlen+=BlkLenScanHeap(m,H); }

        NEWP(nres,totlen+2,short); NEW(ntyp,totlen+2,short);
        for(i=0; i<=totlen; i++) { NEW(nres[i],nAlpha(A)+2,short); }
        NEW(wt,n+2,double);

        /*** 1. determine the number of residues at each position ***/
        for(s=1; s<=n; s++){
           I = info[s]; E = SeqScanInfo(I);
           for(j=1,m=1; m<=nblks; m++){
                seq = XSeqPtr(E) + SiteScanInfo(m,I) - 1;
                for(i=1; i<=BlkLenScanHeap(m,H); i++){
                    r=seq[i]; 
		    if(nres[j][r] == 0) { ntyp[j]++; } 
		    nres[j][r]++; j++;
                }
           }
        }

        /*** 2. determine the sequence weights. ***/
        for(max=0.,min=DBL_MAX,s=1; s<=n; s++){
           I = info[s]; E = SeqScanInfo(I);
	   w = 0.0; 
           for(j=1,m=1; m<=nblks; m++){
                seq = XSeqPtr(E) + SiteScanInfo(m,I) - 1;
                for(i=1; i<=BlkLenScanHeap(m,H); i++){
                    w+= 1.0/ (double) (ntyp[j]*nres[j][seq[i]]); j++; 
                }
           }
	   /*** w /= (double) (j-1); /****/
	   /***/ w /= (double) totlen; /****/
	   wt[s]=w;
	   max = MAXIMUM(double,w,max);
	   min = MINIMUM(double,w,min);
        }

	HG = Histogram("number residue types", 0, 30,1.0); /****/
        for(j=1; j<=totlen; j++){ IncdHist(ntyp[j],HG); }
	/*** PutHist(stdout,60,HG);  /****/
	NilHist(HG);

        /*** 3. normalize the weights. ***/
	HG = Histogram("sequence weights", 0, 1,0.02); /****/
	/** HG = Histogram("sequence weights", 0, 100,1.0); /*****/
        for(N=0.0,s=1; s<=n; s++){
	   wt[s] /= max; 
	   IncdHist(wt[s],HG);
	   N += wt[s];
           // printf("seq %d: weight = %g\n",s,wt[s]);
        }
	/*** PutHist(stdout,60,HG);  /****/
	NilHist(HG);

fprintf(stdout,"effective number of sequences N = %f\n",N);

	/** 4. print out Dirichlet parameters **/
#if 0
        for(i=1; i<=totlen; i++) {
	   printf("(");
	   for(r=0; r<=nAlpha(A); r++){
		if(nres[i][r] > 0) {
			d = N*(double)nres[i][r]/(double) n;
			printf("%.2f,",d);
		} else printf("0,");
		if(r==10) printf("\n ");
	   }
	   printf("%.4f)\n",N);
	}
#endif

        for(i=0; i<=totlen; i++) { free(nres[i]); }
        free(nres); free(ntyp); free(info);
	return wt;
}

void	process_scanheap(snh_typ H)
{
	e_type	E;
	Int4	end,i,s,m,n,nblks;
	float	p,score;
	double	mean,var;
	sni_typ	I;
	static double cut=-1;

	if(H->sort != NULL) return ;
	n=nScanHeap(H);
	nblks=nblksScanHeap(H);
	NEW(H->not_okay,nScanHeap(H)+2,Int4);
	NEW(H->sort,nScanHeap(H)+2,sni_typ);
	NEW(H->use,nblks+2,BooLean);
	for(i=1; (I=DelMinScanHeap(H))!=NULL; i++){ H->sort[i]=I; }
	end=i;
	for(m=1; m<=nblks; m++){
		mean=MeanHist(H->hg[m]);
		var=VarianceHist(H->hg[m]);
		if(n < 10 || mean >= H->use_cutoff) H->use[m]=TRUE;
	}
	for(i=1; i<end; i++){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  H->not_okay[i]=0;
	  for(m=1; m<=nblks; m++){
	     if(H->use[m]){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		if(p < H->neglog10Ecut) H->not_okay[i]++; 
	     }
	  }
	}
	H->number=end-1;
}

Int4	AdjustBlksScanHeap(snh_typ H)
{
	Int4	m,n,nblks;

	process_scanheap(H);
	nblks=nblksScanHeap(H);
	for(n=0,m=1; m<=nblks; m++) if(H->use[m]) n++;
	return n;
}


void	PutRptsScanHeap(FILE *fptr, Int4 left, Int4 right, snh_typ H,
	a_type A)
{ PutRptsScanHeapFull(fptr, left, right, H->number, 1,H, A); }

void	PutTopRptsScanHeap(FILE *fptr, Int4 left, Int4 right, Int4 number, 
	snh_typ H, a_type A)
{ PutRptsScanHeapFull(fptr, left, right, number, 1,H, A); }

Int4	PutMinRptsScanHeap(FILE *fptr, Int4 left, Int4 right, Int4 min_rpt, 
	snh_typ H, a_type A)
{ PutRptsScanHeapFull(fptr, left, right, H->number, min_rpt,H, A); return 0; }

void	PutRptsScanHeapFull(FILE *fptr, Int4 left, Int4 right, Int4 number, 
	UInt4 min_rpt, snh_typ H, a_type A)
{
	e_type	E;
	Int4	s,i,j,j0,m,m2,m0,N0,n,nblks,gap,s2,rpts;
	Int4	start,end,r,len;
	float	score,p;
	sni_typ	I;
	Int4	bad=0;

	assert(min_rpt > 0);
	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);
	if(number > H->number) number=H->number;
	for(i=1; i<=number; i++){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
	  if((rpts*N0) != nblks) print_error("PutRptsScanHeap( ) input error");
          // if(rpts > 1 || H->not_okay[i] <= H->max_not_ok)  // OLD
          if(rpts >= min_rpt) // NEW
	  {
	    for(r=1, m=1; r<=rpts; r++, m+=N0){
		start=SiteScanInfo(m,I);
		m2 = m + N0 - 1;
		s=SiteScanInfo(m2,I);
		m0 = ((m2-1) % N0) + 1;
		end = m+N0;
		for(bad=0, j=m; j < end; j++){
		   p = scoreScanInfo(j,I);
		   j0 = ((j-1) % N0) + 1;
		   if(H->use[j0] && p < H->neglog10Ecut){ bad++; }
		}
		if(bad <= H->max_not_ok){
		  end = s + H->len[m0] - 1;
		  PutSubSeq(fptr,start-left,end+right,E,A);
		}
	    }
	  }
	}
}

void	PutSelexAlnScanHeap(FILE *fptr, snh_typ H, a_type A)
{
	e_type	E;
	Int4	s,i,j,m,m0,N0,n,nblks,gap,s2,rpts;
	Int4	start,end,r;
	Int4	*maxgap;
	float	score,p;
	sni_typ	I;
	char	c;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);
	NEW(maxgap,N0+3,Int4);
	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  for(m=1; m<=nblks; m++){
	    for(m=1; m<=nblks; m++){
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		if(m0 < nblks){
		   s2=SiteScanInfo(m+1,I);
		   gap = s2 - (s + H->len[m0]);
		   if(gap > maxgap[m0]) maxgap[m0]=gap;
		}
	    }
	  }
	}
	for(i=1; i<=H->number; i++){
	 if(H->not_okay[i] <= H->max_not_ok){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  PutSeqID2(fptr,E,20); 

/**** TEST: REPEATS ****/
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
/**** TEST: REPEATS ****/

	  for(m=1; m<=nblks; m++){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		end = s + H->len[m0] - 1;
		for(j = s; j <= end; j++){
		   r = ResSeq(j,E); 
		   fprintf(fptr,"%c",AlphaChar(r,A));
		}
		/***** put insert lengths here *****/
		if((m%N0) != 0){
		   start = s + H->len[m0];
		   end=SiteScanInfo(m+1,I);
		   for(gap=0,j = start; j < end; j++){
			r = ResSeq(j,E); 
			c = AlphaChar(r,A); c = tolower(c);
			fprintf(fptr,"%c",c); gap++;
		   }
		   for(j=gap; j <= maxgap[m0]; j++) {
			fprintf(fptr,"-"); 
		   }
		} 
		if((m%N0) == 0 && m < nblks) {
			fprintf(fptr,"\n"); PutSeqID2(fptr,E,20); 
		} 
	  }
	  fprintf(fptr,"\n");
	 }
	}
	fprintf(fptr,"\n\n");
	free(maxgap);
}

void	PutFAAlnScanHeap(FILE *fptr, snh_typ H, a_type A)
{
	e_type	E;
	Int4	s,i,j,m,m0,N0,n,nblks,gap,s2,rpts;
	Int4	start,end,r,len;
	Int4	*maxgap;
	float	score,p;
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);
	NEW(maxgap,N0+3,Int4);
	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  for(len=0,m=1; m<=nblks; m++){
	    for(m=1; m<=nblks; m++){
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		if(m0 < nblks){
		   s2=SiteScanInfo(m+1,I);
		   gap = s2 - (s + H->len[m0]);
		   if(gap > maxgap[m0]) maxgap[m0]=gap;
		}
	    }
	  }
	}
	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  fprintf(fptr,">"); PutSeqInfo(fptr,E); 

/**** TEST: REPEATS ****/
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
/**** TEST: REPEATS ****/

	  for(len=0,m=1; m<=nblks; m++){
/**** NEW ******/
		if((m%N0) != 1){
		   start = s + H->len[m0];
		   end=SiteScanInfo(m,I);
		   start = start + ((end - start)/2);
		   for(gap=0,j = start; j <= end; j++){
			r = ResSeq(j,E); 
			fprintf(fptr,"%c",AlphaChar(r,A)); len++; gap++;
		   	if(len % 70 == 0) fprintf(fptr,"\n");
		   }
		} 
/**************/
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		end = s + H->len[m0] - 1;
		for(j = s; j <= end; j++){
		   r = ResSeq(j,E); 
		   fprintf(fptr,"%c",AlphaChar(r,A)); len++;
		   if(len % 70 == 0) fprintf(fptr,"\n");
		}
		/***** put insert lengths here *****/
		if((m%N0) != 0){
		   start = s + H->len[m0];
		   end=SiteScanInfo(m+1,I);
/*** modify to print out have of gap ***/
		   end= start + ((end - start)/2);
/*** modify to print out have of gap ***/
		   for(gap=0,j = start; j <= end; j++){
			r = ResSeq(j,E); 
			fprintf(fptr,"%c",AlphaChar(r,A)); len++; gap++;
		   	if(len % 70 == 0) fprintf(fptr,"\n");
		   }
/******/
		   for(j=gap; j <= (maxgap[m0]/2); j++) {
			fprintf(fptr,"-"); len++;
		   	if(len % 70 == 0) fprintf(fptr,"\n");
		   }
/******/
		} 
		if((m%N0) == 0 && m < nblks) {
			fprintf(fptr,"\n\n>"); PutSeqInfo(fptr,E); len=0;
		} 
	  }
	  fprintf(fptr,"\n\n");
	}
	fprintf(fptr,"\n\n");
	free(maxgap);
}

void	PutGapsScanHeap(FILE *fptr, snh_typ H)
{
	e_type	E;
	Int4	s,i,m,m0,N0,n,nblks,gap,s2,rpts;
	float	score,p;
	sni_typ	I;
	float	p2;
	h_type	*HG;
	char	str[100];

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);

	NEW(HG,nblks+2,h_type);
	for(m=0; m<=nblks; m++){
	   sprintf(str,"gaps: motif %d",m);
	   HG[m] = Histogram(str, 0, 200,5.0);
	}

	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);

/**** TEST: REPEATS ****/
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
/**** TEST: REPEATS ****/

	  for(m=1; m<=nblks; m++){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		/***** put insert lengths here *****/
		if(m0 < N0){
		   s2=SiteScanInfo(m+1,I);
		   gap = s2 - (s + H->len[m0]);
		} else {
		   s2=LenSeq(E)+1;
		   gap = s2 - (s + H->len[m0]);
		} 
	
		IncdHist(gap,HG[m0]);
		if(m==1) IncdHist(s-1,HG[0]);
	  }
	}
	for(m=0; m<=N0; m++){ PutHist(fptr,60,HG[m]); NilHist(HG[m]); }
	free(HG);
}

void	PutMSAScanHeap(FILE *fptr, snh_typ H, a_type A)
{
	e_type	E;
	Int4	s,i,j,m,r,m0,N0,n,nblks,gap,s2,rpts,cols,Nseq,b;
	float	p;
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);
	for(cols=0,m=1; m<=N0;m++) cols+=H->len[m];
	for(Nseq=0,i=1; i<=H->number; i++) {
	        I=H->sort[i];
		rpts=RptsScanInfo(I); Nseq+=rpts;
	        nblks = nBlksScanInfo(I);
	        if(rpts*N0 != nblks) print_error("PutMSAScanHeap( ) error");
	}
	fprintf(fptr,"//\nID   XXX\nAC   A00000\nDE   xxx\nLPR  100.00.\n");
	fprintf(fptr,"NU   %d blocks; %d columns.\n",N0,cols);
	fprintf(fptr,"SQ   %d sequences.\n",Nseq);
	for(i=1; i<=H->number; i++){
	  I=H->sort[i]; E = SeqScanInfo(I);
	  rpts = RptsScanInfo(I);
	  for(r=1; r <= rpts; r++){
	  	fprintf(fptr,"     "); PutSeqInfo2(fptr,E);
	  }
	}
	for(b=1; b<=N0; b++){
	  fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
			b,H->len[b],H->len[b]);
	  fprintf(fptr,"MAP  100.0 (50.0 only; 50.0 w/o).\n");
	  fprintf(fptr,"AL   ");
	  for(j=1; j <= H->len[b]; j++) fprintf(fptr,"*"); fprintf(fptr,"\n");
	  for(i=1; i<=H->number; i++){
	     I=H->sort[i];
	     E = SeqScanInfo(I);
	     rpts = RptsScanInfo(I);
	     for(r=1; r <= rpts; r++){
		m = b + ((r - 1) * N0);
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		PutSeqRegionFormatMSA(fptr,s, H->len[b], E, p, A);
	     }
	  }
	}
	fprintf(fptr,"\n");
}

void	PutInfoScanHeap(FILE *fptr, snh_typ H, a_type A)
{ put_info_scan_heap(fptr, H, FALSE, A); }

void	PutFormatMSAScanHeap(FILE *fptr, snh_typ H, a_type A)
{ put_info_scan_heap(fptr, H, TRUE, A); }

void	put_info_scan_heap(FILE *fptr, snh_typ H, BooLean msaformat, a_type A)
{
	e_type	E;
	Int4	s,i,m,m0,N0,n,nblks,gap,s2,rpts,gapped;
	float	score,p;
	char	method;
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);

	gapped=0;
	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
#if 0	/** NEW **/
	  PutScanInfo(fptr, I, H->len, msaformat, A);
#endif
#if 1	/** OLD **/
	  method = MethodScanInfo(I);
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  if(!msaformat) { 
		fprintf(fptr,"[%3.2f] ", -score);
	  } else fprintf(fptr,"     ");
	  PutSeqInfo2(fptr,E);

/**** TEST: REPEATS ****/
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
	  if(rpts > 1) fprintf(fptr,"  (%d repeats)\n",rpts);
/**** TEST: REPEATS ****/

	  for(m=1; m<=nblks; m++){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
	        if(!msaformat) fprintf(fptr," %4.1f: ",p);
		m0 = ((m-1) % N0) + 1;
		if(msaformat) {
		   PutSeqRegionFormatMSA(fptr,s, H->len[m0], E, p, A);
		} else {
		   PutSeqRegion2(fptr,s,H->len[m0],E,0,A);
		   /***** put insert lengths here *****/
		  if(m < nblks){
		   s2=SiteScanInfo(m+1,I);
		   gap = s2 - (s + H->len[m0]);
		   if(m0 < N0 && gap >= 200) fprintf(fptr," (%d)<--LARGE_DOMAIN_INSERT",gap);	
		   else if(m0 < N0 && gap >= 150) fprintf(fptr," (%d)<--DOMAIN_INSERT",gap);	
		   else if(m0 < N0 && gap >= 100) fprintf(fptr," (%d)<--VERY_LARGE_INSERT",gap);	
		   else if(m0 < N0 && gap >= 50) fprintf(fptr," (%d)<--LARGE_INSERT",gap);	
		   else fprintf(fptr," (%d)",gap);	
		  } else {
		   s2=LenSeq(E)+CtermExtendSeq(E)+1;
		   gap = s2 - (s + H->len[m0]);
		   if(gap < 0) gap = 0;
		   fprintf(fptr," (%d)",gap);	
		  } 
		  if(H->use[m0]) fprintf(fptr,"\n");
		  else fprintf(fptr,"*\n");
		} 
		if((m%N0) == 0 && m < nblks) fprintf(fptr,"\n");
	  }
	  if(islower(method)) {
		gapped++;
		fprintf(fptr,"  (gapped E-value used)");
	  }
#endif
	  if(H->not_okay[i] <= H->max_not_ok) fprintf(fptr,"\n");
	  else fprintf(fptr,"  (rejected)\n");
	}
	for(m=1; m<=nblks; m++){
	    m0 = ((m-1) % N0) + 1;
	   if(!H->use[m0]){ fprintf(fptr," (note: '*' = block rejected)"); break; }
	}
	if(gapped > 0) fprintf(fptr,
		"\n (%d additional sequences detected by gap function)\n\n",gapped);
	else fprintf(fptr,"\n\n");
}

void	PutMaskSeqScanHeap(FILE *fptr, snh_typ H, BooLean neg, a_type A)
{
	e_type	E;
	Int4	s,i,j,k,m,m0,N0,n,nblks,rpts;
	unsigned char *seq;
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);

	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  E = SeqScanInfo(I);
	  seq = SeqPtr(E);
	  fprintf(fptr,">");	
	  PutSeqInfo(fptr,E);

	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);

	  for(j=1,m=1; m<=nblks; m++){
		s=SiteScanInfo(m,I);
		while(j < s) {
		   if(neg) fprintf(fptr,"%c",tolower(AlphaChar(seq[j],A)));
		   // if(neg) fprintf(fptr,"x");
		   else fprintf(fptr,"%c",AlphaChar(seq[j],A));
		   if(j %60 == 0) fprintf(fptr,"\n");
		   j++;
		}
		m0 = ((m-1) % N0) + 1;
		for(k=1; k <= H->len[m0]; k++){
		   if(neg) fprintf(fptr,"%c",AlphaChar(seq[j],A));
		   else fprintf(fptr,"x");
		   if(j %60 == 0) fprintf(fptr,"\n");
		   j++;
		}
	  }
	  s = LenSeq(E);
	  while(j <= s) {
		if(neg) fprintf(fptr,"%c",tolower(AlphaChar(seq[j],A)));
		// if(neg) fprintf(fptr,"x");
		else fprintf(fptr,"%c",AlphaChar(seq[j],A));
		if(j %60 == 0) fprintf(fptr,"\n");
		j++;
	  }
	  fprintf(fptr,"\n\n");
	}
}

void	PutSeqBlkScanHeap(char *name, snh_typ H, Int4 Nflank, Int4 Cflank, a_type A)
/*********************************************************************
 output fasta files for each of the conserved blocks.
 *********************************************************************/
{
	e_type	E;
	Int4	s,i,j,k,m,m0,N0,n,nblks,rpts;
	unsigned char *seq;
	FILE	**fp;
	char	str[20];
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);
	NEWP(fp, N0 + 2, FILE);
	for(m=1; m <= N0; m++){ 
		sprintf(str,".%dsq",m); fp[m] = open_file(name,str,"w"); 
	}

	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  E = SeqScanInfo(I);
	  seq = SeqPtr(E);
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);

	  for(m=1; m<=nblks; m++){
		m0 = ((m-1) % N0) + 1;
	  	fprintf(fp[m0],">");	
	  	PutSeqInfo(fp[m0],E);
		s=SiteScanInfo(m,I);
		j = MAXIMUM(Int4,1,s-Nflank);
		for( ; j < s; j++){
		   fprintf(fp[m0],"%c",AlphaChar(seq[j],A));
		}
		for(k=1; k <= H->len[m0]; k++){
		   fprintf(fp[m0],"%c",AlphaChar(seq[j],A)); j++;
		}
	        s = MINIMUM(Int4,j+Cflank,LenSeq(E));
		for( ; j < s; j++){
		   fprintf(fp[m0],"%c",AlphaChar(seq[j],A));
		}
	  	fprintf(fp[m0],"\n\n");
	  }
	}
	for(m=1; m <= N0; m++) fclose(fp[m]);
	free(fp);
}

void 	PutRasmolScanHeap(FILE *fptr, snh_typ H)
/** write rasmol scripts (in one file) for sequences in heap **/
{
	e_type	E;
	Int4	s,i,j,m,m0,N0,n,nblks,gap,s2,rpts;
	float	score,p;
	const char	*colors[7]={"magenta","red","orange","yellow","green","cyan","blue"};
	char	*id,str[30],su;
	Int4	color;
	sni_typ	I;

	process_scanheap(H);
	n=nScanHeap(H);
	N0 = nblks=nblksScanHeap(H);

	for(i=1; i<=H->number; i++){
	  I=H->sort[i];
	  score = ScoreScanInfo(I);
	  E = SeqScanInfo(I);
	  id = SeqKey(E);
	  id = strstr(id,"pdb|");
	  if(id == NULL) continue;
	  sscanf(id,"pdb|%[^|]|%c",str,&su);
	  for(j=0; str[j] != 0; j++) if(isupper(str[j])) str[j]=tolower(str[j]);
	  fprintf(fptr,"\nload ~/fasta/pdb/full/pdb%s.ent\n",str);
	  /*****/
	  fprintf(fptr,"set background white\nrenumber 1\n");
	  fprintf(fptr,"select all\nbackbone off\nwireframe off\n");
	  if(isalpha(su)) fprintf(fptr,"select *%c\ncenter *%c\n",su,su);
	  fprintf(fptr,"set strands 1\nstrands on\ncolor gray\n");

/**** TEST: REPEATS ****/
	  nblks = nBlksScanInfo(I);
	  rpts = RptsScanInfo(I);
/**** TEST: REPEATS ****/

	  for(m=1; m<=nblks; m++){
		p=scoreScanInfo(m,I);
		s=SiteScanInfo(m,I);
		m0 = ((m-1) % N0) + 1;
		
		color = m0%7;
		if(isalpha(su)) {
	  	  fprintf(fptr,"select *%c and %d-%d\n",su,s,s+H->len[m0]-1);
		} else {
	  	  fprintf(fptr,"select %d-%d\n",s,s+H->len[m0]-1);
		}
	  	fprintf(fptr,"strands off\ncartoons on\n");
		if(p >= H->neglog10Ecut){
			fprintf(fptr,"color %s\n",colors[color]);
		} else {
			fprintf(fptr,"color gray\n");
		}
	  }
	  fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
}

Int4	PutSeqScanHeap(FILE *fptr, snh_typ H, a_type A)
{ return PutTopSeqScanHeap(fptr, H->number, H, A); }

Int4	PutTopSeqScanHeap(FILE *fptr, Int4 number, snh_typ H, a_type A)
{
	e_type	E;
	Int4	i,n;
	sni_typ	I;

	process_scanheap(H);
	if(number > H->number) number=H->number;
	for(n=0,i=1; i<=number; i++){
	  I=H->sort[i];
	  E=SeqScanInfo(I);
	  if(H->not_okay[i] <= H->max_not_ok) { n++; PutSeq(fptr,E,A); }
	}
	return n;
}

Int4	AddSetScanHeap(snh_typ H, set_typ S, BooLean AddAll)
/** set bits in block B corresponding to the sequences in hit list **/
{
	Int4	i,n;
	sni_typ	I;

	process_scanheap(H);
	for(n=0,i=1; i<=H->number; i++){
	  I=H->sort[i];
	  if(AddAll || (H->not_okay[i] <= H->max_not_ok)){
		n++; AddSet(SeqIDScanInfo(I),S); 
	  }
	}
	return n;
}

Int4	PutFailSeqScanHeap(FILE *fptr, snh_typ H, a_type A)
{
	e_type	E;
	Int4	i,n;
	sni_typ	I;

	process_scanheap(H);
	for(n=0,i=1; i<=H->number; i++){
	  I=H->sort[i];
	  E=SeqScanInfo(I);
	  if(H->not_okay[i] > H->max_not_ok) { n++; PutSeq(fptr,E,A); }
	}
	return n;
}

Int4	NumFailSeqScanHeap(snh_typ H)
{
	e_type	E;
	Int4	i,n;
	sni_typ	I;

	process_scanheap(H);
	for(n=0,i=1; i<=H->number; i++){ if(H->not_okay[i] > H->max_not_ok) n++; }
	return n;
}

