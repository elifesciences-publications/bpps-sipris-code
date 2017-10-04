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

#include "cluster.h"

e_type	**ClusterSeqs(e_type *ListE, char method, double cutoff, Int4 *counts, 
	a_type A, Int4 *Nset)
{ return ClusterSeqs(ListE, method, cutoff, counts, A, Nset, 11); }

e_type	**ClusterSeqs(e_type *ListE, char method, double cutoff, Int4 *counts, 
	a_type A, Int4 *Nset, Int4 T)
/** cluster the sequences in ListE into distinct sets **/
{ 
	Int4	i,j,n,m,N,score,minimum=1,v;
	Int4	total,cut,len,low,high,s,s1,s2;
	Int4	sizedbs=0,max_set,num_set,max_i,num_i;
	double	lambda,K,H,*pr,*freq,avelen,p,Eval;
	e_type	E1,E2,**ListE2;
	gb_typ  B;
	BooLean	flag,related;
	ds_type	sets;

	for(N=0; ListE[N+1] != NULL; ) N++;
	NEWP(ListE2, N +2, e_type);
	if(sizedbs <1) sizedbs = N-1;
	sets = DSets(N);
	if(method == 'e'){
          for(total=0, i=0; i<= nAlpha(A); i++) total += counts[i];
          avelen = (double) total/(double)N;
          NEW(freq,nAlpha(A)+2,double);
          for(i=0;i<=nAlpha(A);i++) freq[i]=(double)counts[i]/(double)total;
          low=lowAlphaR(A); high=highAlphaR(A); len=high - low + 1;
          NEW(pr,len+1,double);
          for(i=0; i<=nAlpha(A); i++){
            for(j=0; j<=nAlpha(A); j++){
                v=(Int4)valAlphaR(i,j,A)-low; pr[v]+=freq[i]*freq[j];
            }
          }
          if(!karlin(low,high,pr,&lambda,&K,&H)){print_error("fatal");}
          free(pr); free(freq);
          // T=(Int4) floor((4.676 + 1.3*log(H))/lambda - 0.5);
          // T=11;
          /** end compute karlin-altschul parameters **/
	}
	// PutDSets(stdout,sets); 
	for(n = 1; n < N; n++){
	    E1 = ListE[n];
	    s1 = findDSets(n,sets);	
	    // PutXSeq(stdout,E1,A);
	    if(method=='e'){
	      len=LenSeq(E1);
	      cut=(Int4)((log((double)(len*avelen))+log(K/0.5))/lambda);
	      B = MkGBlast(50, cut, T, E1, A);
	      for(m = n+1; m <= N; m++){
	        E2=ListE[m];
		score=MatcherGBlastStr(LenSeq(E2),SeqPtr(E2),B,AlphaR(A));
		p=SumStatGBlastStr(len,LenSeq(E2),lambda,H,K,B); 
		if(p > 0.0){
		   Eval = -log10(p*(double)(sizedbs));
		} else Eval = DBL_MAX;
		// fprintf(stderr,"n=%d; m=%d: %g vs %g\n",n,m,Eval,cutoff); 
		if(Eval >= (double) cutoff){
	           s2 = findDSets(m,sets);	
		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
		}
	      } NilGBlast(B);
	    } else {
#if 0	// old slower 'b' method
	        if(method=='b') B = MakeGBlast(T, E1, A);
	        for(m = n+1; m <= N; m++){
	          E2 = ListE[m];
		  if(method=='b'){
			score = MatcherGBlast(NULL,E2,B);
		  } else { score = FastAlnSeqSW(12, 4, E1, E2, A); }
		  if(score >= (Int4) cutoff){ /** then merge sets **/
	           s2 = findDSets(m,sets);	
		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
		  }
	        } if(method=='b') NilGBlast(B);
#else
	        if(method=='b') B = MakeGBlast(T, E1, A);
	        for(m = n+1; m <= N; m++){
	          E2 = ListE[m];
		  if(method=='b'){
			related= FastMatcherGBlast(E2,B,cutoff);
		  } else { related = RelatedSeqs(cutoff, E1, E2, A); }
		  if(related){ /** then merge sets **/
#if 0		  // in purge I am using Sets datatype for edges:
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
#endif
	           s2 = findDSets(m,sets);	
		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
		  }
	        } if(method=='b') NilGBlast(B);
#endif
	   }
	} // PutDSets(stdout,sets); 
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    if(i == findDSets(n,sets)){	
		if(!flag){ num_set++; flag=TRUE; }
		num_i++;
	    }	
	  } if(flag) NEW(ListE2[num_set],num_i+2,e_type);
	}
	*Nset = num_set;
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    E2 = ListE[n]; s1 = findDSets(n,sets);	
	    if(i==s1){
		if(!flag){ num_set++; flag=TRUE; s++; }
		num_i++; ListE2[num_set][num_i] = E2;
	    }	
	  }
	} NilDSets(sets);
	return ListE2;
}

int	PutRepSeq(FILE *fptr, char c, BooLean pval, e_type *EList,
		BooLean UseLabeled,a_type A)
// c = method
{ 
	Int4	i,j,n,m,sum,N,score,minimum=1,max,v;
	Int4	total,cut,len,low,high,ave,maxseq,maxlen;
	double	lambda,T,K,H,*pr,avelen,p,maxscore,*avescore;
	ss_type	P = NULL;
	e_type	E1,E2,*NewEList;
	gb_typ  B;
        static double freq[21] = {  0.0, /*  X
            C     G     A     S     T     N     D     E     Q     K */
        0.025,0.074,0.074,0.057,0.051,0.045,0.054,0.054,0.034,0.058,
        /*  R     H     W     Y     F     V     I     L     M     P */
        0.052,0.026,0.013,0.032,0.047,0.073,0.068,0.099,0.025,0.039 };


	for(n=1; EList[n] != NULL; n++) ; 
	NEW(NewEList,n+3,e_type);
	for(n=1; EList[n] != NULL; n++) NewEList[n]=EList[n]; 
	n--; 
	NEW(avescore,n+3,double);
	P = Array2SeqSet(NewEList, n, "null",A);
	// PutSeqSetPIDs(stderr,P); 
	// fprintf(stderr,"     ");
	for(total=0, n = 1; n <= NSeqsSeqSet(P); n++) {
	  	E1 = SeqSetE(n,P); total += LenSeq(E1);
		// fprintf(stderr," %4d",n);
	}
	avelen = (double) total/(double) NSeqsSeqSet(P);
	// fprintf(stderr," %4s\n","ave");
	if(pval){
          low=lowAlphaR(A); high=highAlphaR(A); len=high - low + 1;
          NEW(pr,len+1,double);
          for(i=0; i<=nAlpha(A); i++){
            for(j=0; j<=nAlpha(A); j++){
                v=(Int4)valAlphaR(i,j,A)-low; pr[v]+=freq[i]*freq[j];
            }
          }
          if(!karlin(low,high,pr,&lambda,&K,&H)){print_error("fatal");}
          free(pr); 
	  T=(Int4) floor((4.676 + 1.3*log(H))/lambda - 0.5);
          /** end compute karlin-altschul parameters **/
	}
	for(max=0,sum = N = 0, n = 1; n <= NSeqsSeqSet(P); n++){
	    E1 = SeqSetE(n,P);
	    // fprintf(stderr," %4d",n);
	    if(pval){
	      len=LenSeq(E1);
	      cut=(Int4)((log((double)(len*avelen))+log(K/0.5))/lambda);
	      B = MkGBlast(50, cut, T, E1, A);
	      for(m = 1; m < n; m++){
	        E2=SeqSetE(m,P);
		score=MatcherGBlastStr(LenSeq(E2),SeqPtr(E2),B,AlphaR(A));
		p=SumStatGBlastStr(len,LenSeq(E2),lambda,H,K,B); p = -log(p);
		// fprintf(stderr," %4d",(Int4)p);
	      	avescore[n] += score; avescore[m] += score;
	      }
	      NilGBlast(B);
	    } else {
	      if(c=='b') B = MakeGBlast(11, E1, A);
	      for(m = 1; m < n; m++){
	        E2 = SeqSetE(m,P);
		if(c=='b') score = MatcherGBlast(NULL,E2,B);
		else { score = FastAlnSeqSW(12, 4, E1, E2, A); }
		// fprintf(stderr," %4d",score);
	      	avescore[n] += score; avescore[m] += score;
	      }
	      if(c=='b') NilGBlast(B);
	   }
	   // fprintf(stderr,"\n");
	}

	maxscore=-99999999.;
	for(n = 1; n <= NSeqsSeqSet(P); n++){
	      E1 = SeqSetE(n,P); ave = avescore[n];
#if 0	// DEBUG: UseLabeled option...
	      fprintf(stderr,"%c ",kingdomSeq(E1));
	      // PutSeqID(stderr,E1);
#endif
#if 1	// UseLabeled option...
	      if(UseLabeled){
		if(PdbSeq(E1)) ave += 100000.0;
		else if(LabeledSeq(E1)) ave += 50000.0;
	      }
#else
	      if(UseLabeled && LabeledSeq(E1)) ave += 50000.0;
#endif
	      if(maxscore < ave) { maxscore=ave; maxseq = n; maxlen = LenSeq(E1); }
	      else if(maxscore == ave) {
		if(maxlen < LenSeq(E1)){
		     maxscore=ave; maxseq = n; maxlen = LenSeq(E1);
		}
	      }
		// fprintf(stderr,"len = %d\n",LenSeq(E1));
	}
	// fprintf(stderr,"\n\tmaximum sum score = %g\n", (double)maxscore);
	E1 = SeqSetE(maxseq,P);
#if 0	// DEBUG: UseLabeled option...
	fprintf(stderr,"\n***** maximum sum score = %g: %c ****\n",
		(double)maxscore,kingdomSeq(E1));
	// PutSeqID(stderr,E1);
	// fprintf(stderr,"\n*****************************\n", (double)maxscore);
#endif
	PutSeq(fptr,E1,A);
	free(NilSeqSetRtnSeqs(P));
	free(avescore);
	return 1;
}

Int4	RepSetCluster(FILE *fp, char *DBS_NAME, Int4 sizedbs, BooLean xnu, 
	char mode, Int4 cutoff, a_type  A)
{ return RepSetCluster(fp,DBS_NAME,sizedbs,xnu,mode, cutoff,FALSE,A,11); }

Int4	RepSetCluster(FILE *fp, char *DBS_NAME, Int4 sizedbs, BooLean xnu, 
	char mode, Int4 cutoff,BooLean UseLabeled, a_type  A)
{ return RepSetCluster(fp,DBS_NAME,sizedbs,xnu,mode,cutoff,UseLabeled,A,11); }

Int4	RepSetCluster(FILE *fp, char *DBS_NAME, Int4 sizedbs, BooLean xnu, 
	char mode, Int4 cutoff,BooLean UseLabeled, a_type  A,Int4 T)
{
	Int4	n,N,minimum=1,Nset,s;
	ss_type	P = NULL;
	e_type	*ListE,**ListE2;
	BooLean	verbose = TRUE,shuffle=FALSE,pval=FALSE;

	if(xnu) P = MkXSeqSet(DBS_NAME,A);
	else P = SeqSet(DBS_NAME,A);
	if(isupper(mode)){
	  ListE2 = ClusterGPSI(mode,cutoff,P,A,&Nset);
	} else {
	  N = NSeqsSeqSet(P);
	  NEW(ListE, N+2,e_type);
	  for(n = 1; n <= N; n++) ListE[n] = SeqSetE(n,P); 
	  if(sizedbs <1) sizedbs = N-1;
	  ListE2 = ClusterSeqs(ListE,mode,cutoff,CntsSeqSet(P),A,&Nset,T);
	  free(ListE);
	}
	for(s=1; s<=Nset; s++)
	   { PutRepSeq(fp,'b',pval,ListE2[s],UseLabeled,A); free(ListE2[s]); }
	free(ListE2); 
	NilSeqSet(P);
	return Nset;
}

e_type	**ClusterGPSI(char method, double cutoff, ss_type data, 
	a_type A, Int4 *Nset)
{ return ClusterGPSI(method, cutoff, data, A, Nset,11); }

e_type	**ClusterGPSI(char method, double cutoff, ss_type data, 
	a_type A, Int4 *Nset,Int4 Threshold)
// THIS FUNCTION IS REDUNDANT WITH ONE BELOW AND CAN BE ELIMINATED!!!
/** cluster the sequences in ListE into distinct sets **/
{ 
	Int4	i,n,m,N;
	Int4	s,s1,s2;
	Int4	max_set,num_set,max_i,num_i;
	double	*scores;
	e_type	E1,E2,**ListE2;
	BooLean	flag;
	ds_type	sets;

        Int4    T=11,gap_open=11,gap_extend=1;
        double  x_parameter=15.0;
        Int4    hpsz=10;
        // double  Ethresh=0.01,Ecutoff=0.01;
        double  Ethresh=1.0,Ecutoff=1.0;
        BooLean top_only=TRUE;
	T = Threshold;

	N = NSeqsSeqSet(data);
	sets = DSets(N); // PutDSets(stdout,sets); 
	for(n = 2; n <= N; n++){
	    E1 = SeqSetE(n,data);
	    s1 = findDSets(n,sets);	
	    // PutXSeq(stdout,E1,A); 
	    gpsi_type gpsi(E1,data,Ethresh,Ecutoff,x_parameter,hpsz,1);
// fprintf(stderr,"%4d: ",n);
	    scores=gpsi.FastSearch(top_only,n-1,T,method);
	    for(m = 1; m < n; m++){
// fprintf(stderr,"%4.0f ",scores[m]);
		  if(scores[m] >= cutoff){ /** then merge sets **/
	           s2 = findDSets(m,sets);	
		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
	          }
	    } free(scores);
// std::cerr << std::endl;
	} // PutDSets(stdout,sets); 
	NEWP(ListE2, N +2, e_type);
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    if(i==findDSets(n,sets)){	
		if(!flag){ num_set++; flag=TRUE; }
		num_i++;
	    }	
	  }
	  if(flag) NEW(ListE2[num_set], num_i + 2, e_type);
	} *Nset = num_set;
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    E2=SeqSetE(n,data); s1=findDSets(n,sets);	
	    if(i==s1){
		if(!flag){ num_set++; flag=TRUE; s++; }
		num_i++; ListE2[num_set][num_i] = E2;
	    }	
	  }
	} NilDSets(sets);
	return ListE2;
}

ds_type	ClusterGPSI(char method, double cutoff, ss_type data)
/** cluster the sequences in ListE into distinct sets **/
{ 
	Int4	n,m,N,s1,s2;
	double	*scores;
	e_type	E1;
	ds_type	sets;
	a_type	A=SeqSetA(data);

	N = NSeqsSeqSet(data);
        Int4    T=11,gap_open=11,gap_extend=1;
        double  x_parameter=15.0;
        Int4    hpsz=N;
        // double  Ethresh=0.01,Ecutoff=0.01;
        double  Ethresh=1.0,Ecutoff=1.0;
        BooLean top_only=TRUE;

	sets = DSets(N); // PutDSets(stdout,sets); 
	for(n = 2; n <= N; n++){
	    E1 = SeqSetE(n,data);
	    s1 = findDSets(n,sets);	
	    // PutXSeq(stdout,E1,A); 
	    gpsi_type gpsi(E1,data,Ethresh,Ecutoff,x_parameter,hpsz,1);
// fprintf(stderr,"%4d: ",n);
	    scores=gpsi.FastSearch(top_only,n-1,T,method);
	    for(m = 1; m < n; m++){
// fprintf(stderr,"%4.0f ",scores[m]);
		  if(scores[m] >= cutoff){ /** then merge sets **/
	            s2 = findDSets(m,sets);	
		    if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
	          }
	    } free(scores);
// std::cerr << std::endl;
	} 
	// PutDSets(stdout,sets); 
	return sets;
}

e_type  **ClusterGPSI(char method, double cutoff, ss_type data, Int4 *Nset)
{
	a_type	A=SeqSetA(data);
	Int4	i,n,N,s,s1,max_set,num_set,max_i,num_i;
	e_type	E2,**ListE2;
	BooLean	flag;

	N = NSeqsSeqSet(data);
	ds_type	sets = ClusterGPSI(method, cutoff, data);
	NEWP(ListE2, N +2, e_type);
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    if(i==findDSets(n,sets)){	
		if(!flag){ num_set++; flag=TRUE; }
		num_i++;
	    }	
	  }
	  if(flag) NEW(ListE2[num_set], num_i + 2, e_type);
	} *Nset = num_set;
	for(num_set=0,max_i=max_set=s=0,i=1; i <= N; i++){
	  for(num_i=0,flag=FALSE,n=1; n <= N; n++){
	    E2=SeqSetE(n,data); s1=findDSets(n,sets);	
	    if(i==s1){
		if(!flag){ num_set++; flag=TRUE; s++; }
		num_i++; ListE2[num_set][num_i] = E2;
	    }	
	  }
	} NilDSets(sets);
	return ListE2;
}

Int4	*ClusterIntoSets(double cutoff, ss_type data, Int4 *Nset)
{
	Int4	i,j,s,set,N;
	Int4	*Set;

	N = NSeqsSeqSet(data);
	ds_type	sets = ClusterGPSI('E', cutoff, data);
	NEW(Set, N +2, Int4);
        for(set=0,i=1;i<=N;i++){
           s = findDSets(i,sets);
           if(s == i){
	      set++;
              for(j=1;j<=N;j++){
                 if(s==findDSets(j,sets)) Set[j]=set;
              }
           }
        } NilDSets(sets);
	*Nset = set;
	return Set;
}

