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

// double	**p2c_typ::RtnCaMeanDist(FILE *scrfp, Int4 Row)
adh_typ	*p2c_typ::RoundOne(FILE *scrfp, Int4 Row, Int4 OmitCol, Int4 OmitSeq)
// check quality of an alignment based on Calpha variances
{
   Int4		i,j,n,I,C,R,k1,k2,cmaSeqID;
   Int4		Start,Stop,*NumDist,NumCol;       // number of distances for |i-j|
   double	d,m;
   e_type	pdbE,cmaE;
   char		c1,c2,str[505];
   cma_typ	cma=0;
   FILE		*fp=0;

   if(IN_CMA[Row] == 0) return 0; // No sequences were in this set.
   NumCol=LengthCMSA(1,IN_CMA[1]);
   if(Row == 0){ 	// shouldn't need this....??  Used for SARP;
	assert(Row > 0);
	Start=1; Stop=Number - 1; fp=open_file(OutFile,"_All.ca_scores","w"); 
   } else {	// run on one cma file only.
      assert(Row <= Number && Row > 0); Start=Stop=Row;
      fp=scrfp; 
   } if(begin==0 && end ==0){ begin=1; end=NumCol; }
#if 0
   double **MeanDistMtrx; NEWP(MeanDistMtrx, NumCol + 20, double);
   for(i=1; i <= NumCol; i++) NEW(MeanDistMtrx[i], NumCol + 20, double);
#endif

   //********************** 1. Get differences in Calpha distances *********************
   NEW(NumDist, MaxSqDist +9, Int4);	// Save Number of distances for each separation...
   adv_typ	*adv=0;
   // create dheap to store most deviant hits.
   // adh_typ = atomic distance heap typ
   adh_typ *adh = new adh_typ(MaxDataPoints);
   // adh_typ *adh = new adh_typ(2*MaxDataPoints);
   // adh_typ *adh = new adh_typ(NumCol*NumCol);
   h_type HG0=Histogram("variances in distance between Ca atoms.",0,100,0.1);
   // h_type HG0=Histogram("variances in distance between Ca atoms.",0,MaxMeanDist,bin_size);
   h_type HG2=Histogram("numbers of data points for computing variances.",0,mpdb->NumPDB,10);
   h_type HG_D=Histogram("number of residues between Ca atoms within the chain.",MinDist,500,20);

   double Fraction=0.50;
   Int4 Num,MinNumData=(Int4) floor(Fraction*(double)NumSeqsCMSA(IN_CMA[1]));
   if(Diagnostic && fp) fprintf(fp,"Read Faults: %s\n",Diagnostic);
   if(fp) fprintf(fp,"============> %d pdb files; Number=%d; MinNumData = %d; NumCompare = %d\n",
		mpdb->NumPDB,NumSeqsCMSA(IN_CMA[1]),MinNumData,(NumCol*(NumCol-1))/2);
MaxMeanDist=5000;
// MinCaDist=0.0; MaxCaDist=8.0;
// MinCaDist=8.0; MaxCaDist=14.0;
// MinCaDist=14.0; MaxCaDist=20.0;
// MinCaDist=20.0; MaxCaDist=26.0;
// MinCaDist=26.0; MaxCaDist=32.0;
// MinCaDist=32.0; MaxCaDist=500.0;

   enum fault { off=0, nil=1, zero=2, close=3, noatom=4, lost=5, miss=6, okay=7 };
   Int4	Fault[9]; for(i=0; i < 8; i++) Fault[i]=0; 
this->LastTotalData=0;

   for(k1=begin; k1 < end; k1++){
     if(k1==OmitCol && OmitSeq < 0) continue;	// new: find contribution to ARSD score for each column!!
     for(k2=k1+1; k2 <= end; k2++){
       if(k2==OmitCol && OmitSeq < 0) continue;
       if(Begin && !((Begin <= k1 && k1 <= End) || (Begin <= k2 && k2 <= End))) continue;
       h_type HG=Histogram("Distances between Ca atoms (Angstroms)",0,50,1.0);
       h_type HG_bst=Histogram("Distances between Ca atoms (Angstroms)",0,50,1.0);

       for(I=1; I <=mpdb->NumPDB; I++){
	if(I == OmitSeq && OmitCol < 0) continue;
	if(I == OmitSeq && (OmitCol == k1 || OmitCol == k2)) continue;
	BooLean	found=FALSE;
	double d_bst=99999999.0;
        for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
#if 1	// Column Pair Sets...Skip sequences that don't correspond...
          if(ColPairSet && ColPairSet[k1][k2]){
		// pdbIC=mpdb->pdbSeq[I][C]; assert(pdbIC);
		if(MapPdb2cmaSq[I] == 0) continue;
		if(!MemberSet(MapPdb2cmaSq[I][C],ColPairSet[k1][k2])) continue;
          }
#endif
	  for(R=1; R <= NumRpts[I][C]; R++){
	    if(RelevantSeq[I][C][R] < Start || RelevantSeq[I][C][R] > Stop)
		{ Fault[off]++; continue; }

	    if(Col2pdbSeq[I][C][R]==0){ Fault[nil]++; continue; }
	    i = Col2pdbSeq[I][C][R][k1]; j = Col2pdbSeq[I][C][R][k2];
	    if(i == 0 || j == 0){ Fault[zero]++; continue; }
	    // what about indels?  This will throw things off...
	    if(abs(i-j) < MinDistInSeq){ Fault[close]++;  continue; }
	    if(mpdb->Calpha[I][C][i] && mpdb->Calpha[I][C][j]){ 
	        d=DistanceAtoms(mpdb->Calpha[I][C][i],mpdb->Calpha[I][C][j]);
#if 0
if(d < MinCaDist){ continue; }
if(d > MaxCaDist){ continue; }
#endif
		D=abs(i-j); assert(D > 0);
   	        // D=abs(k1-k2); assert(D > 0);
		IncdHist(D,HG_D);
		if(D > MaxSqDist) D = MaxSqDist+1;
		NumDist[D]++;	// Store this for use below...
		if(d < d_bst) d_bst=d;
	        IncdHist(d,HG); found=TRUE;
	    }
	  }
        } if(found) IncdHist(d_bst,HG_bst); else Fault[lost]++; 
       }
       h_type HG_use=HG_bst; // HG_use=HG;
       // if(MeanHist(HG_use) <= MaxMeanDist)
       // if(MeanHist(HG_use) <= MaxMeanDist && NumDataHist(HG_use) >= Target)
this->LastTotalData += NumDataHist(HG_use);
       if(NumDataHist(HG_use) >= Target)
       {
       	   Fault[okay]++;
	   // MeanDistMtrx[k1,k2]=MeanDistMtrxp[k2][k1]=MeanHist(HG_use);
	   // if(fp) PutHist(fp,60,HG_use);
           d=sqrt(VarianceHist(HG_use));
	   m = MeanHist(HG_use); 
#if 1
	   d = d/m; 
#else
	   d = pow(d,1.2)/m; 
#endif
	   Num = NumDataHist(HG_use);
           IncdHist(d,HG0); IncdHist(Num,HG2);
	   adv = new adv_typ(k1,k2,d,m,Num);
           adh->Insert(adv,-fabs(1.0-d));   	// key is minimum d.
	   // ^ save those with the least variability from 1.0.
       } else Fault[miss]++;
	NilHist(HG); NilHist(HG_bst);
     }
   } free(NumDist);
   if(fp){
	PutHist(fp,60,HG0); PutHistArray(fp,HG0); 
   	fprintf(fp," mean stdev = %.5f (%d data points).\n",MeanHist(HG0),NumDataHist(HG0));
   }
   this->LastARSD=MeanHist(HG0); this->LastNumData=NumDataHist(HG0);
   if(fp){
	fprintf(fp,"Faults: off=%d; nil=%d; zero=%d; close=%d; noatom=%d; lost=%d; miss=%d; okay=%d.\n",
		Fault[0],Fault[1], Fault[2],Fault[3],Fault[4],Fault[5],Fault[6],Fault[7]);
	PutHist(fp,60,HG2); PutHist(fp,60,HG_D); 
   } NilHist(HG0); NilHist(HG2); NilHist(HG_D);
   return adh;
}

adv_typ **p2c_typ::GetAdvList(Int4 &Num,adh_typ *adh)
{
	double	key;
	Int4	item,i;
	adv_typ	*adv,**ADV;
	NEW(ADV,adh->NumItems() +5, adv_typ*);
	for(i=0; (adv=adh->DelMin(&key,&item)) != 0; ){
	    i++; ADV[i]=adv; 
	} delete adh; Num=i;
	return ADV;
}

adv_typ	**p2c_typ::RunAndRtnADV(FILE *scrfp,Int4 &NumADV)
// check quality of MAPGAPS alignment based on Calpha variances
{

   Int4		ii,jj,i,j,n,I,C,R,k1,k2,cmaSeqID;
   double	dd,d,m;
   e_type	pdbE,cmaE;
   char		c1,c2,str[505];
   Int4		Start,Stop,Row=1;
   cma_typ	cma=0;
   FILE		*fp=0;

   assert(IN_CMA[Row] != 0);
   adv_typ **ADV=this->GetAdvList(NumADV,this->RoundOne(scrfp, Row));
   if(Row == 0){
      Start=1; Stop=Number - 1;
      fp=open_file(OutFile,"_All.ca_scores","w");
   } else {	// run on one cma file only.
      assert(Row <= Number && Row > 0); Start=Stop=Row;
      if(scrfp==0){
        // sprintf(str,"%s_%d",OutFile,Row); fp=open_file(str,".ca_scores","w"); 
      } else {
	fp=scrfp; 
	fprintf(fp,"================== %s_%d ca_scores ===================\n",OutFile,Row);
      }
   }

   //********************** 1. Get differences in Calpha distances *********************
   Int4 NumCol=LengthCMSA(1,IN_CMA[1]);
   adv_typ	*adv=0;
   // create dheap to store most deviant hits.
   // adh_typ = atomic distance heap typ
   // adh_typ *adh = new adh_typ(MaxDataPoints);
   // adh_typ *adh = new adh_typ(MaxDataPoints);
   adh_typ *adh = new adh_typ(NumCol*NumCol);

   double Fraction=0.50,*Distance;
   Int4 Num,MinNumData=(Int4) floor(Fraction*(double)NumSeqsCMSA(IN_CMA[1]));
   if(fp) fprintf(fp,"============> %d pdb files; Number=%d; MinNumData = %d\n",
		mpdb->NumPDB,NumSeqsCMSA(IN_CMA[1]),MinNumData);
// MaxMeanDist=5000; MaxMeanDist=2.0;

   //=============== Store the 'target' # of the least variant sequences ===============

   // find most misaligned...
   Int4 hpsz=mpdb->NumPDB;
   double *SumDev;   NEW(SumDev,mpdb->NumPDB +3,double);
   Int4	  *NumData;  NEW(NumData,mpdb->NumPDB +3,Int4);
   dh_type misaln_dH=dheap(hpsz*NumADV + 10,4);
   Int4	   *NumDt,xx=0,*JJ,*II;
   NEW(JJ,hpsz*NumADV +3, Int4); NEW(II,hpsz*NumADV +3, Int4);
   NEW(NumDt,NumADV +3, Int4);

   enum fault { off=0, nil=1, zero=2, close=3, noatom=4, lost=5, miss=6, okay=7 };
   Int4	Fault[9]; for(i=0; i <= 8; i++) Fault[i]=0; 

   for(xx=0,jj=1; jj <= NumADV; jj++){
       k1=ADV[jj]->ColI(); k2=ADV[jj]->ColJ(); 
       if(Begin && !((Begin <= k1 && k1 <= End) || (Begin <= k2 && k2 <= End))) continue;
       h_type HG=Histogram("distances between Ca atoms (Angstroms)",0,50,1.0);
       h_type HG_bst=Histogram("distances between Ca atoms (Angstroms)",0,50,1.0);

       dh_type dH=dheap(hpsz + 10,4); NEW(Distance,hpsz+10,double);;
       for(I=1; I <=mpdb->NumPDB; I++){
	BooLean	found=FALSE;
	double d_bst=99999999.0;
        for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
#if 1	// Column Pair Sets...Skip sequences that don't correspond...
          if(ColPairSet && ColPairSet[k1][k2]){
		// pdbIC=mpdb->pdbSeq[I][C]; assert(pdbIC);
		if(MapPdb2cmaSq[I] == 0) continue;
		Int4 sq=MapPdb2cmaSq[I][C];
		if(!MemberSet(sq,ColPairSet[k1][k2])) continue;
          }
#endif
	  for(R=1; R <= NumRpts[I][C]; R++){
	    if(RelevantSeq[I][C][R] < Start || RelevantSeq[I][C][R] > Stop){ Fault[off]++; continue; }
	    if(Col2pdbSeq[I][C][R]==0){ Fault[nil]++; continue; }
	    i = Col2pdbSeq[I][C][R][k1]; j = Col2pdbSeq[I][C][R][k2];
	    if(i == 0 || j == 0){ Fault[zero]++; continue; }
	    if(abs(i-j) < MinDistInSeq){ Fault[close]++;  continue; }
	    if(mpdb->Calpha[I][C][i] && mpdb->Calpha[I][C][j]){ 
	        d=DistanceAtoms(mpdb->Calpha[I][C][i],mpdb->Calpha[I][C][j]);
		// m=ADV[jj]->Mean(); d=10*d/m;
		m=ADV[jj]->Mean(); 
#if 1	// sqrt
		d=d/m;
#else
		d=pow(d,1.2)/m;
#endif

		if(d < d_bst) d_bst=d;
	        IncdHist(d,HG); found=TRUE;
	    } else Fault[noatom]++;
	  }
        }
	if(found){
	    insrtHeap(I,(keytyp)ADV[jj]->Variance(),dH); Distance[I]=d_bst;

	    SumDev[I] += d_bst; NumData[I]++;
	    xx++; insrtHeap(xx,-(keytyp)d_bst,misaln_dH); II[xx]=I; JJ[xx]=jj; NumDt[jj]++;
	} else Fault[lost]++;
       } // PutHeap(fp,dH); 
       // for(ii=0; ii < Target && !emptyHeap(dH); ii++)
       for(ii=0; !emptyHeap(dH); ii++)
       {
	   dd=(double)minkeyHeap(dH); I=delminHeap(dH); d=Distance[I]; IncdHist(d,HG_bst); 
       } Nildheap(dH); free(Distance);  dH=0; Distance=0;
       h_type HG_use=HG_bst; // HG_use=HG; // PutHist(fp,60,HG_bst);
       // if(MeanHist(HG_use) <= MaxMeanDist)
       {
	   // MeanDistMtrx[k1,k2]=MeanDistMtrxp[k2][k1]=MeanHist(HG_use);
	   // PutHist(fp,60,HG_use);
           d=sqrt(VarianceHist(HG_use));
	   m = MeanHist(HG_use); 
#if 1	// don't change this...
	   d = d/m;
#else
	   d = pow(d,1.2)/m;
#endif
	   Num = NumDataHist(HG_use);
	   assert(Num >= Target);
	   Fault[okay]++;
	   adv = new adv_typ(k1,k2,d,m,Num);
           adh->Insert(adv,-fabs(1.0-d));   	// key is minimum d.
       } NilHist(HG); NilHist(HG_bst);
   }
#if 1	// find most misaligned...
   h_type all_HG=Histogram("all normalized Ca distances (Angstroms)",0,10,0.10);
   for(ii=0; !emptyHeap(misaln_dH); ii++){
	dd=-(double)minkeyHeap(misaln_dH); xx=delminHeap(misaln_dH);
	jj=JJ[xx]; I=II[xx]; k1=ADV[jj]->ColI(); k2=ADV[jj]->ColJ(); 
#if 0
	if(ii < Target) fprintf(stderr,"%d('%s'): d_bst[%d][%d]=%.2f.\n",xx,mpdb->pdb_file[I],k1,k2,dd);
#endif
	if(NumDt[jj] > Target) IncdHist(dd,all_HG);
	
   } Nildheap(misaln_dH); free(JJ); free(II); free(NumDt);
   if(fp) fprintf(fp," full stdev = %.5f (%d data points; %d cols).\n",
		sqrt(VarianceHist(all_HG)),NumDataHist(all_HG),LengthCMSA(1,IN_CMA[1]));
#if 1
   fprintf(stdout,"%s: full stdev = %.5f (%d data points; %d cols).\n",
		OutFile,sqrt(VarianceHist(all_HG)),NumDataHist(all_HG),LengthCMSA(1,IN_CMA[1]));
   this->FinalARSD=sqrt(VarianceHist(all_HG));
   this->FinalPnts=NumDataHist(all_HG);
   this->FinalCols=LengthCMSA(1,IN_CMA[1]);
#endif
   // fprintf(fp,"mean = %.5f; stdev = %.5f.\n",MeanHist(all_HG),sqrt(VarianceHist(all_HG)));
   if(fp) PutHist(fp,60,all_HG); NilHist(all_HG);
   for(jj=1; jj <= NumADV; jj++){ delete ADV[jj]; ADV[jj]=0; } free(ADV);
   for(dd=0,jj=0,I=1; I <=mpdb->NumPDB; I++){
	if(NumData[I] > 1){
	   if(0) fprintf(stderr,"%d ('%s'): %.2f (%d)\n",I,mpdb->pdb_file[I],
					SumDev[I]/(double)NumData[I],NumData[I]);
	    dd+=SumDev[I]/(double)NumData[I]; jj++;
	}
   } // fprintf(stderr,"   averge = %.2f (%d)\n",dd/(double)jj,jj);
   // enum fault { off=0, nil=1, zero=2, close=3, noatom=4 };
   if(fp) fprintf(fp,"Faults: off=%d; nil=%d; zero=%d; close=%d; noatom=%d; lost=%d; okay=%d.\n",
		Fault[0],Fault[1], Fault[2],Fault[3],Fault[4],Fault[5],Fault[7]);
   free(SumDev); free(NumData); // exit(1);
#endif
   ADV=this->GetAdvList(NumADV,adh);	// deletes adh
   return ADV;
}

void    p2c_typ::PutAdvList(FILE *fp,Int4 NumADV, adv_typ **ADV,Int4 NumShown, char *message)
{
   Int4		ii,jj,i,j,Num;
   double	dd,d,m;
   adv_typ	*adv;

   // h_type HG0=Histogram("Variances in distance between Ca atoms.",0,MaxMeanDist,bin_size);
   h_type HG0=Histogram("Variances in distance between Ca atoms.",0,100,0.1);
   h_type HG2=Histogram("Numbers of data points for computing variances",0,mpdb->NumPDB,10);
   for(i=1,jj=1; jj <= NumADV; jj++){
	adv=ADV[jj]; d=adv->Variance(); Num=adv->Number();
        IncdHist(d,HG0); IncdHist(Num,HG2);
   } 
   if(message){
      fprintf(fp,"%s: mean variance = %.5f (%d data points).\n",
		message,MeanHist(HG0),NumDataHist(HG0));
   }
   PutHist(fp,60,HG0); PutHistArray(fp,HG0); NilHist(HG0); 
   PutHist(fp,60,HG2); NilHist(HG2); 
   for(i=1,jj=1; jj <= NumADV; i++,jj++){
        if(i <= NumShown){
	  ADV[jj]->Put(fp);
          if(i % 50 == 0){ fprintf(fp,"============= %d ===============\n",i); }
          else if(i % 10 == 0){ fprintf(fp,"------------- %d ---------------\n",i); }
        } // delete ADV[jj]; ADV[jj]=0;
   } // free(ADV);
}

void	p2c_typ::RoundTwo(FILE *scrfp, Int4 Row,BooLean PerformControl)
// check quality of MAPGAPS alignment based on Calpha variances
{
   Int4		ii,jj,i,j,NumADV,Num;
   FILE		*fp=scrfp; 
   double	dd,d,m;

   adv_typ **ADV=this->RunAndRtnADV(scrfp,NumADV);
   if(scrfp) this->PutAdvList(scrfp, NumADV,ADV);
   for(jj=1; jj <= NumADV; jj++){ delete ADV[jj]; ADV[jj]=0; } free(ADV);
   return;
}

