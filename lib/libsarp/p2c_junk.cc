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

#include "psc_typ.h"

void	p2c_typ::GetStructAlnScores(FILE *scrfp, Int4 Row,BooLean PerformControl)
// check quality of MAPGAPS alignment based on Calpha variances
{
#if 1
	this->RoundTwo(scrfp, Row,PerformControl); return;
#endif
   Int4		i,j,n,I,C,R,k1,k2,cmaSeqID;
   double	d,m;
   adv_typ	*adv=0;
   e_type	pdbE,cmaE;
   char		c1,c2,str[505];
   double	**Distance;     // distances for each |i-j| separation in seq.
   Int4		*NumDist,NumCol;       // number of distances for |i-j|
   Int4		Start,Stop;
   cma_typ	cma=0;
   FILE		*fp=0;

   if(IN_CMA[Row] == 0) return; // No sequences were in this set.
   NumCol=LengthCMSA(1,IN_CMA[1]);
   if(Row == 0){
      Start=1; Stop=Number - 1;
      fp=open_file(OutFile,"_All.ca_scores","w");
   } else {	// run on one cma file only.
      assert(Row <= Number && Row > 0); Start=Stop=Row;
      if(scrfp==0){
        sprintf(str,"%s_%d",OutFile,Row); fp=open_file(str,".ca_scores","w"); 
      } else {
	fp=scrfp; 
	fprintf(fp,"================== %s_%d ca_scores ===================\n",OutFile,Row);
      }
   }
   if(begin==0 && end ==0){ begin=1; end=NumCol; }

   //********************** 1. Get differences in Calpha distances *********************
   NEW(NumDist, MaxSqDist +9, Int4);	// Save Number of distances for each separation...
   // create dheap to store most deviant hits.
   // adh_typ = atomic distance heap typ
   adh_typ *adh = new adh_typ(5000);
   h_type HG0=Histogram("Variances in distance between Ca atoms.",0,MaxMeanDist,bin_size);
   h_type HG2=Histogram("Numbers of data points for computing variances.",0,mpdb->NumPDB,10);

   double Fraction=0.50;
   Int4 Num,MinNumData=(Int4) floor(Fraction*(double)NumSeqsCMSA(IN_CMA[1]));
   fprintf(fp,"============> %d pdb files; Number=%d; MinNumData = %d\n",
		mpdb->NumPDB,NumSeqsCMSA(IN_CMA[1]),MinNumData);
MaxMeanDist=5000;

   for(k1=begin; k1 < end; k1++){
     for(k2 = k1+ MinDist; k2 <= end; k2++){
       if(KeyCol && !(KeyCol == k1 || KeyCol == k2)) continue;
       if(Begin && !((Begin <= k1 && k1 <= End) || (Begin <= k2 && k2 <= End))) continue;
       h_type HG=Histogram("Distances between Ca atoms (Angstroms)",0,50,1.0);
       h_type HG_bst=Histogram("Distances between Ca atoms (Angstroms)",0,50,1.0);

       for(I=1; I <=mpdb->NumPDB; I++){
	BooLean	found=FALSE;
	double d_bst=99999999.0;
        for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
	  for(R=1; R <= NumRpts[I][C]; R++){
	    if(RelevantSeq[I][C][R] < Start || RelevantSeq[I][C][R] > Stop) continue;
	    if(Col2pdbSeq[I][C][R]==0) continue;
	    i = Col2pdbSeq[I][C][R][k1];
	    j = Col2pdbSeq[I][C][R][k2];
	    if(i == 0 || j == 0) continue;
#if 1		// what about indels?  This will throw things off...
	    Int4 DistInSeq=abs(i-j);
	    if(DistInSeq < MinDistInSeq) continue;
#endif
	    if(mpdb->Calpha[I][C][i] && mpdb->Calpha[I][C][j]){ 
	        d=DistanceAtoms(mpdb->Calpha[I][C][i],mpdb->Calpha[I][C][j]);
		D=abs(i-j); assert(D > 0);
   	        // D=abs(k1-k2); assert(D > 0);
		if(D > MaxSqDist) D = MaxSqDist+1;
		NumDist[D]++;	// Store this for use below...
		if(d < d_bst) d_bst=d;
	        IncdHist(d,HG); found=TRUE;
	        if((KeyCol && I == 1) || (K1 && k1==K1 && k2 == K2)){
	          pdbE=mpdb->pdbSeq[I][C];
	          c1=AlphaChar(ResSeq(i,pdbE),AB); c2=AlphaChar(ResSeq(j,pdbE),AB);
		  StrSeqID(str, 100, pdbE);
	          printf("(%d,%d) %c%d(%d)...%c%d(%d) = %.2f Angstroms.  %s\n",
			I,C,c1,FakeToRealCMA(n,k1,IN_CMA[Row]),k1,
			c2,FakeToRealCMA(n,k2,IN_CMA[Row]),k2,d,str);
	        }
	    }
	  }
        } if(found) IncdHist(d_bst,HG_bst);
       }
#if 1
       h_type HG_use=HG_bst; 
#else
       h_type HG_use=HG;
#endif
// MinNumData=80;
// if(NumDataHist(HG_use) >= 40)
       if(MeanHist(HG_use) <= MaxMeanDist){
	   // PutHist(fp,60,HG_use);
           d=sqrt(VarianceHist(HG_use));
	   if(d >= MinVar){
	     m = MeanHist(HG_use);
	     d = 10*d/m;
       	     // VAR[Begin]=d;	// for scatter diagram
	     if(K1 && k1==K1 && k2 == K2) PutHist(fp,60,HG_use); 
	     Num = NumDataHist(HG_use);
	     if(Num < MinNumData) continue;
             IncdHist(d,HG0); IncdHist(Num,HG2);
	     adv = new adv_typ(k1,k2,d,m,Num);
	     if(FindBest) adh->Insert(adv,d);	// best superimposed...
	     else adh->Insert(adv,-d);		// worst superimposed...
           }
       }
       // PutHist(stdout,60,HG);
       NilHist(HG); NilHist(HG_bst);
     }
   }
   Int4	Item;
   double key;
   h_type HG1=Histogram("Weighted variances in distance between Ca atoms.",0,10,0.25);
   h_type HG3=Histogram("Top variances in distance between Ca atoms.",0,MaxMeanDist,bin_size);
   h_type HG4=Histogram("Numbers of data points for computing top variances.",0,mpdb->NumPDB,10);
   for(i=1; adv=adh->DelMin(&key,&Item); i++){
	d=adv->Variance(); Num=adv->Number(); IncdMHist(d,Num,HG1); IncdHist(d,HG3); IncdHist(Num,HG4);
	if(i <= 200){
	  adv->Put(fp);
	  if(i % 50 == 0){ fprintf(fp,"============= %d ===============\n",i); }
	  else if(i % 10 == 0){ fprintf(fp,"------------- %d ---------------\n",i); }
	} delete adv;
   } delete adh;
   PutHist(fp,60,HG0); PutHistArray(fp,HG0);
   PutHist(fp,60,HG2); NilHist(HG2);
   PutHist(fp,60,HG3); PutHistArray(fp,HG3); NilHist(HG3);
   PutHist(fp,60,HG4); NilHist(HG4);
   PutHist(fp,60,HG1); NilHist(HG1);

//******************************************* NEW **********************************************

   double **ResFreq=ColResFreqsCMSA(1,IN_CMA[Start]);
   for(Int4 AA=Start+1; AA <= Stop; AA++){
	double **tmp_freq=ColResFreqsCMSA(1,IN_CMA[AA]);
	for(i=1; i<=NumCol; i++){
	   for(Int4 r=0; r<=nAlpha(AB); r++) ResFreq[i][r] += tmp_freq[i][r];
	} free(tmp_freq);
   }
   if(Start < Stop){
     double d=Stop-Start+1;
     for(i=1; i<=NumCol; i++){
	for(Int4 r=0; r<=nAlpha(AB); r++) ResFreq[i][r] = ResFreq[i][r]/d;
     }
   }
   double *re; NEW(re,NumCol +4, double);
   double AveRe=0.0;
   for(i=1; i<=NumCol; i++){
	double Re=0.0;
        if(Begin && !((Begin <= i && i <= End))) continue;
	for(Int4 r=1; r<=nAlpha(AB); r++){
		if(ResFreq[i][r] > 0.0){
			Re += ResFreq[i][r]*log(ResFreq[i][r]/blosum62freq[r]);
	        }
		// fprintf(stdout,"freq[%d][%c]=%.3f\n",i,AlphaChar(r,AB),ResFreq[i][r]);
	}
	// fprintf(stdout,"%3d: %.2f\n",i,Re);
	re[i]=Re;
	AveRe += Re; 
   } 
   // RE[Begin]=AveRe;
   fprintf(fp,"AveRe=%.4f\n",AveRe/(double)(End-Begin+1));
   fprintf(fp,"\n");
   
   if(!PerformControl){ NilHist(HG0); return; }
//********************************** Test versus Control relative entropy ***************
   double *xvalT, *yvalT,*ProbT=0, P=0.0;
   double *xvalC, *yvalC,*ProbC;
   Int4		TotalT;
   Int4		TotalC;
   Int4 NumBinsT = RtnHistArray(&xvalT,&yvalT,&TotalT,HG0);
 if(TotalT > 0 && NumBinsT > 0){
   double tiny=0.0000000001;
   NEW(ProbT,NumBinsT+5,double); 
   for(i=0; i<=NumBinsT+1; i++){
	ProbT[i] = ((double) yvalT[i] + tiny)/((double)TotalT + (NumBinsT+2)*tiny);
	P+=ProbT[i];
	// if(ProbT[i] > 0.0) fprintf(stdout,"%.2f: prob[%d]=%.3f\n",xvalT[i],i,ProbT[i]);
   }
   // fprintf(stdout,"total prob = %.3f\n",P);
 }
   
   NilHist(HG0);

//******************************************* NEW **********************************************
   //**************** 5. Store Calpha distances as function of seq leng ***************
   //************************* & Randomly shuffle distances ***************************

   dh_type *dH=NULL;		// one heap for every seq. distance.
   NEWP(Distance, MaxSqDist +9, double);
   NEW(dH, MaxSqDist +9, dh_type);
   for(D=1; D <= MaxSqDist +1; D++){	// arrays for storing info.
	NEW(Distance[D], NumDist[D] +9, double);
	dH[D]=dheap(NumDist[D]+9,4);
	NumDist[D]=0;			// zero out for use as a pointer.
   }
   Int4 DistInSeq=abs(i-j);
   h_type HG_D=0;
   // HG_D=Histogram("Distances between residues in sequences",0,500,5.0);
   for(k1=begin; k1 < end; k1++){
     for(k2 = k1+ MinDist; k2 <= end; k2++){
       if(KeyCol && !(KeyCol == k1 || KeyCol == k2)) continue;
       if(Begin && !((Begin <= k1 && k1 <= End) || (Begin <= k2 && k2 <= End))) continue;
       for(I=1; I <=mpdb->NumPDB; I++){
        for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
	 for(R=1; R <= NumRpts[I][C]; R++){
	  if(RelevantSeq[I][C][R] < Start || RelevantSeq[I][C][R] > Stop) continue;
	  if(Col2pdbSeq[I][C][R]){
	    i = Col2pdbSeq[I][C][R][k1];
	    j = Col2pdbSeq[I][C][R][k2];
	    if(i == 0 || j == 0) continue;
#if 1		// what about indels?  This will throw things off...
	    Int4 DistInSeq=abs(i-j);
	    if(DistInSeq < MinDistInSeq) continue;
#endif
	    if(mpdb->Calpha[I][C][i] && mpdb->Calpha[I][C][j]){ 
	       d=DistanceAtoms(mpdb->Calpha[I][C][i],mpdb->Calpha[I][C][j]);
#if 1
   	       D=abs(i-j); assert(D > 0);
#else
   	       D=abs(k1-k2); assert(D > 0);
#endif
	       if(D > MaxSqDist) D = MaxSqDist+1;
	       if(HG_D) IncdHist(D,HG_D);
	       NumDist[D]++; Distance[D][NumDist[D]] = d; 
	       insrtHeap(NumDist[D],(keytyp)Random(),dH[D]);
	    }
         }
	}
       }
      }
     }
   }
   if(HG_D) PutHist(fp,60,HG_D);
   if(HG_D) NilHist(HG_D);
   
   //****************** 7. Perform control ******************
   HG0=Histogram("Control for variance in Ca distances.",0,MaxMeanDist,bin_size);
   for(k1=begin; k1 < end; k1++){
     for(k2 = k1+ MinDist; k2 <= end; k2++){
       if(KeyCol && !(KeyCol == k1 || KeyCol == k2)) continue;
       if(Begin && !((Begin <= k1 && k1 <= End) || (Begin <= k2 && k2 <= End))) continue;
       h_type HG=Histogram("Distances between Ca atoms (Angstroms)",0,50,1.0);
       for(I=1; I <=mpdb->NumPDB; I++){
	BooLean	found=FALSE;
	double d_bst=99999999.0;
        for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
	 for(R=1; R <= NumRpts[I][C]; R++){
	  if(RelevantSeq[I][C][R] < Start || RelevantSeq[I][C][R] > Stop) continue;
	  if(Col2pdbSeq[I][C][R]){
	    i = Col2pdbSeq[I][C][R][k1];
	    j = Col2pdbSeq[I][C][R][k2];
	    if(i == 0 || j == 0) continue;
#if 1		// what about indels?  This will throw things off...
	    Int4 DistInSeq=abs(i-j);
	    if(DistInSeq < MinDistInSeq) continue;
#endif
	    if(mpdb->Calpha[I][C][i] && mpdb->Calpha[I][C][j]){ 
#if 1		// use sequence distance or column distance.
   	       D=abs(i-j); assert(D > 0);
#else
   	       D=abs(k1-k2); assert(D > 0);
#endif
	       if(D > MaxSqDist) D = MaxSqDist+1;
	       n=delminHeap(dH[D]);
	       assert(n > 0);
	       d=Distance[D][n];
	       if(d < d_bst){ d_bst=d; found=TRUE; }
	    }
	  }
	 }
        } if(found) IncdHist(d,HG);
       }
       if(MeanHist(HG) <= MaxMeanDist){
           d=sqrt(VarianceHist(HG));
	   if(d >= MinVar){
	     m = MeanHist(HG);
             IncdHist(d,HG0);
           }
       } NilHist(HG);
     }
   }
   PutHist(fp,60,HG0);
   PutHistArray(fp,HG0);
   fflush(fp);

//********************* Compute relative entropy ***********************
   Int4 NumBinsC = RtnHistArray(&xvalC,&yvalC,&TotalC,HG0);
 if(ProbT && TotalC > 0){
   NEW(ProbC,NumBinsC+5,double); 
   double RelEnt=0.0,tiny=0.0000000001;
   for(P=0.0,i=0; i<=NumBinsC+1; i++){
	// ProbC[i] = (double) yvalC[i]/(double)TotalC;
	ProbC[i] = ((double) yvalC[i] + tiny)/((double)TotalC + (NumBinsC+2)*tiny);
	P+=ProbC[i];
	// if(ProbC[i] > 0.0) fprintf(stdout,"%.2f: prob[%d]=%.3f\n",xvalC[i],i,ProbC[i]);

	if(ProbC[i] > 0.0 || ProbT[i] > 0.0){
		// RelEnt += (ProbT[i] + tiny)*log((ProbT[i]+tiny)/(ProbC[i]+tiny));
		RelEnt += ProbT[i]*log(ProbT[i]/ProbC[i]);
	}

   }
   // fprintf(stdout,"total prob = %.3f\n",P);
   fprintf(fp,"Rel. Entropy=%.3f\n",RelEnt);
 }
//********************* Compute relative entropy ***********************
   
   NilHist(HG0);
   for(D=1; D <= MaxSqDist +1; D++){	// arrays for storing info.
	Nildheap(dH[D]); free(Distance[D]);
   } if(scrfp==0) fclose(fp); else fprintf(fp,"\n"); 
}


