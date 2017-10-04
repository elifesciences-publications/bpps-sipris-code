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

#include "sqd_typ.h"
#include "blosum62.h"

#define	TEST_ResSetAll	1
#define	MaxResSets	16

sqd_typ::sqd_typ(cma_typ CMA,cma_typ KeyCMA, double MinFreq,double max_gap_frq,char mode)
{
	init(CMA,KeyCMA,MinFreq,max_gap_frq,mode); 
	// rst->Put(stderr);	
}

// I may want to define this based on the actual alignment at some point...

void    sqd_typ::Free()
{
	Int4	r,s;
	for(s=0; s <= Length; s++){
	   if(Set[s]){
              for(r=0; r<=nAlpha(AB); r++){
#if 0	// debug....
		if(Set[s][r]){
fprintf(stderr,"sqd->Free(): %c%d (%d)\n",AlphaChar(r,AB),s,CardSet(Set[s][r]));
		   NilSet(Set[s][r]);
		}
#else
		if(Set[s][r]) NilSet(Set[s][r]);
#endif
	      }
	      free(Set[s]);
	   } 
#if 1	// new data...
	   if(key_obs[s]) free(key_obs[s]); 
	   if(keyFreq[s]) free(keyFreq[s]); 
#endif
	} free(Set);
#if 1	// new...
	free(key_obs); free(keyFreq); 
#endif
	NilSet(SetI); NilSet(SetJ);
	NilSet(NotSetI); NilSet(NotSetJ); NilSet(USet);
	Int4	r1,r2,rs1,rs2;
	for(r1=1; r1 <= nAlpha(AB); r1++){
	     for(r2=1; r2 <= nAlpha(AB); r2++){
	       for(rs1=1; rs1 <= NumResSets[r1]; rs1++){
	         for(rs2=1; rs2 <= NumResSets[r2]; rs2++) {
			t_type	tab=table[r1][r2][rs1][rs2];
			if(tab) NilTab(tab);
		 } free(table[r1][r2][rs1]);
	       } free(table[r1][r2]);
	     }
	}
	free(cell[1]); free(cell[2]);
	if(CnTabList){
	  for(r=1; r<= NumCnTab; r++){
		if(CnTabList[r]) delete CnTabList[r];
	  } free(CnTabList);
	}
	if(NumSsetStr > 0 && sst_str){
        	for(s=0; s <= NumSsetStr; s++){
			if(sst_str[s]) free(sst_str[s]);
		} free(sst_str);
	}
	delete ctn;
	delete rst;
#if TEST_ResSetAll
	free(ResSetAll);
	free(BigSetArray);
	free(ResSetProb);
	free(BigProbArray);
#endif
	if(keyE) NilSeq(keyE);
	free(ResSetLegal);
	free(BigSetLegal);
	free(NumGaps);
}

void	sqd_typ::init(cma_typ CMA,cma_typ KeyCMA, double MinFreq,double max_gap_frq,
		char mode)
{
	Int4	j,s,r,r1,r2,rs;
	sst_typ	sst;

	table_fp=0;
	assert(CMA); assert(KeyCMA);
	keycma=KeyCMA;
	cma=CMA;
	ctn =0;
	NumClust=0; sst_str=0; NumSsetStr=0;
	CnTabList=0; NumCnTab=0;
	cth=0;
	HG=0;
#if 1	// Newest: Fix both indels and extensions within query seq. afn: 10-30-08
        Int4 start = TruePosCMSA(1,1,KeyCMA);      // ignores offset...
        Int4 end = start+LengthCMSA(1,KeyCMA)-1;   // will jump over gaps...
        keyE=MkSubSeq(start,end,SeqSetE(1,DataCMSA(KeyCMA)));  // Fake sequence...
	assert(LenSeq(keyE) == LengthCMSA(1,KeyCMA));
	for(s=1; s <= LenSeq(keyE); s++){	// set deletions to 'X'
		if(IsDeletedCMSA(1,s,KeyCMA)){ EqSeq(s,0,keyE); EqXSeq(s,0,keyE); }
	}	// NOT yet sure if this really fixes the problem...but works for now.
#elif 0	// older method: WARNING: Don't free keyE in Free() if using this option!!!
	keyE=SeqSetE(1,DataCMSA(KeyCMA)); // Fake seq!!
     	// Try fixing bug relating to conversion of X to A residues...
	// This routine was using A's as if they were real instead of deletions....
	keyE=CopySeq(keyE);	// not freed up!!!
	for(s=1; s <= LenSeq(keyE); s++){
		if(IsDeletedCMSA(1,s,KeyCMA)){ EqSeq(s,0,keyE); EqXSeq(s,0,keyE); }
	}	// NOT yet sure if this really fixes the problem...but works for now.
#else	// next == concensus...
	keyE=SeqSetE(1,DataCMSA(CMA));
#endif
// double  **keyFreq=ColResFreqsCMSA(1, Int4 ***observed, keycma);
	keyFreq=ColResFreqsCMSA(1,&key_obs,keycma);
	assert(MinFreq >= 0.0 && MinFreq <= 1.0);
	MinKeyFrq=MinFreq;
	MaxGapFrq=max_gap_frq;
	assert(MaxGapFrq >= 0.0 && MaxGapFrq < 1.0);
	Tab=0;
	AB = AlphabetCMSA(cma);
	NEW(cell[1],3,Int4);
	NEW(cell[2],3,Int4);
	NumSeqs=NumSeqsCMSA(cma);
	Length=LenSeq(keyE);
	if(Length != LengthCMSA(1,cma)){
		fprintf(stderr,"LenSeq(keyE) = %d; LengthCMSA = %d\n",
			Length,LengthCMSA(1,cma));
		assert(Length == LengthCMSA(1,cma));
		// print_error("Fatal: length inconsistency");
	}
	rst=new rst_typ(mode);
	ResidueSSet = rst->ResidueSSet;
	NumResSets= rst->NumResSets;
	max_residue_set = rst->MaxResSet();
	// rst->Put(stderr);

#if 0	// OLD:
	if(KeyE) keyE=KeyE;
	else keyE=SeqSetE(1,DataCMSA(cma)); // == 'Fake' sequence...
#endif
	// fprintf(stderr,"Creating contingency tables: \n");
	// *********************** Allocate sets... ***************************
	// fprintf(stderr,"Allocating fundamental sets...\n");
	NEWP(Set,Length +3, set_typ);
	for(s=1; s <= Length; s++){
	   NEW(Set[s],nAlpha(AB)+2,set_typ);	// 
	   r = ResSeq(s,keyE);
	   //********* first create non-related set...
	   Set[s][0] = MakeSet(NumSeqs+1);	// create 'X' set...
	   if(r) {
		sst = ResidueSSet[r][0];
	   	for(r2=1; r2 <= nAlpha(AB); r2++){
		   if(MemSset(r2,sst)) Set[s][r2]=MakeSet(NumSeqs+1);
		}
	   }
	}
	// WARNING: Sets contain elements 0..NumSeqs-1
	SetI=MakeSet(NumSeqs+1);
	SetJ=MakeSet(NumSeqs+1);
	NotSetI=MakeSet(NumSeqs+1);
	NotSetJ=MakeSet(NumSeqs+1);
	USet=MakeSet(NumSeqs+1);
	FillSet(USet); DeleteSet(0,USet);
	// *********************** Fill sets... ***************************
#if 1	// Speed up...
	{
	  Int4 sq,s0;
	  unsigned char *qsq=SeqPtr(keyE);
	  for(sq=2; sq <=NumSeqsCMSA(cma); sq++){
	     s0=SitePos(1,sq,1,SitesCMSA(cma)); // start of fake sequence.
	     unsigned char *fsq = SeqSeqSet(sq,DataCMSA(cma));
	     for(s=1; s<=LengthCMSA(1,cma); s++,s0++){
		rs = fsq[s0];
	     	sst = ResidueSSet[qsq[s]][0];
		if(MemSset(rs,sst)) AddSet(sq,Set[s][rs]);
		else if(rs==0) AddSet(sq,Set[s][0]);
	     }
	  }
	}
#else
	for(s=1; s<=LengthCMSA(1,cma); s++){
	   	Int4 rq = ResSeq(s,keyE);
		// r = ResidueCMSA(1,1,s0,cma);
		for(Int4 sq=2; sq <=NumSeqsCMSA(cma); sq++){
               		rs=ResidueCMSA(1,sq,s,cma);
	     		sst = ResidueSSet[rq][0];
			if(MemSset(rs,sst)) AddSet(sq,Set[s][rs]);
			else if(rs==0) AddSet(sq,Set[s][0]);
		}
	}
#endif
	NEW(NumGaps,LengthCMSA(1,cma)+3,Int4);
	for(s=1; s<=LengthCMSA(1,cma); s++) NumGaps[s] = CardSet(Set[s][0]);
	// *********************** Get Tables... ***************************
	sst_typ	sstL,sstR;
	Int4	rs1,rs2;
	for(r1=1; r1 <= nAlpha(AB); r1++){
	     for(r2=1; r2 <= nAlpha(AB); r2++){
	       NEWP(table[r1][r2],NumResSets[r1]+3,t_type);	// 
	       for(rs1=1; rs1 <= NumResSets[r1]; rs1++){
	          NEW(table[r1][r2][rs1],NumResSets[r2]+3,t_type);	// 
	       }
	     }
	}
	Int4 total_tab=0;
	for(r1=1; r1 <= nAlpha(AB); r1++){
	   for(rs1=1; rs1 <= NumResSets[r1]; rs1++){
	     sstL=ResidueSSet[r1][rs1];;
	     sst_typ sst1=SsetLet(r1);
	     for(r2=1; r2 <= nAlpha(AB); r2++){
	        sst_typ sst2=SsetLet(r2);
		for(rs2=1; rs2 <= NumResSets[r2]; rs2++){
	          sstR=ResidueSSet[r2][rs2];
		  t_type T=Table(r1,sst1,sst2,AB); // create a 2x2 table.
		  AddColumnsTab(r2,sstL,sstR,T);
		  // ClearTab(T);
#if 0
Int4 unit=20;
AddTab(0,0,unit,T); AddTab(1,0,unit,T);
AddTab(0,1,unit,T); AddTab(1,1,unit,T);
#endif
		  table[r1][r2][rs1][rs2]=T;
		  // PutTableMod(stderr,table[r1][r2][rs1][rs2]);
		  total_tab++;
		}
	     }
	   }
	} // fprintf(stderr,"\t%d tables made.\n",total_tab);
	sst_typ *sst_ptr;
#if TEST_ResSetAll
	NEW(BigSetArray,Length*(MaxResSets+2),sst_typ);
	NEW(BigProbArray,Length*(MaxResSets+2),double);
	NEWP(ResSetAll,Length+3,sst_typ);
	NEWP(ResSetProb,Length+3,double);
	sst_ptr=BigSetArray;
	double	*prb_ptr=BigProbArray;
	for(s=1; s<=LengthCMSA(1,cma); s++){
	   ResSetAll[s] = sst_ptr;
	   ResSetProb[s] = prb_ptr;
	   for(rs1=1; rs1 <= MaxResSets; rs1++){
		ResSetProb[s][rs1]=DBL_MAX;
	   }
	   sst_ptr += MaxResSets+1;
	   prb_ptr += MaxResSets+1;
	}
#endif
#if 1	// TEST_ResSetLegal
	Int4 ri,ris;
	sst_typ sstI;
	NEW(BigSetLegal,Length*(MaxResSets+3),sst_typ);
	NEWP(ResSetLegal,Length+3,sst_typ);
	sst_ptr=BigSetLegal;
	for(j=1; j <= Length; j++){
	   ResSetLegal[j]=sst_ptr;
	   ri = ResSeq(j,keyE); 
	   for(ris=0; ris <= NumResSets[ri]; ris++){
	        ResSetLegal[j][ris]=ResidueSSet[ri][ris];;
	   }
	   ResSetLegal[j][ris]=0;
	   sst_ptr += NumResSets[ri]+2;
	}
#endif
	// rst->Put(stderr);	
}

#if 1	// return all legal residue sets...
sst_typ **sqd_typ::LegalResSets( )
{
	return ResSetLegal;
}
#endif

sst_typ **sqd_typ::SortedResSets( )
// Sort 
{
	Int4	s,rs,rs0;
	double	key;
	sst_typ	*tmp;
	dh_type	dH=dheap(MaxResSets+1,3);
	NEW(tmp,MaxResSets+2,sst_typ);
	
	for(s=1; s<=Length; s++){
	  for(rs=1; rs <= MaxResSets; rs++){
		tmp[rs]=0;
		if(ResSetProb[s][rs] < 1.0){
			if(ResSetProb[s][rs] == 0.0) key = (double) INT4_MIN;
			else key=log(ResSetProb[s][rs]);
			insrtHeap(rs,key,dH);
		}
	  }
	  for(rs0=0; (rs=delminHeap(dH)) != NULL; ) {
		rs0++; tmp[rs0] = ResSetAll[s][rs];
	  }
	  for(rs=1; rs <= rs0; rs++) ResSetAll[s][rs]=tmp[rs];
	  for( ; rs <= MaxResSets; rs++) ResSetAll[s][rs]=0;
	}
	free(tmp); Nildheap(dH);
	return ResSetAll;
}

#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-30

void cntab1(Int4 **nn, Int4 ni, Int4 nj, double *chisq, double *df, double *prob,
	double *cramrv, double *ccc)
{
	double gammq(double a, double x);
	Int4	nnj,nni,j,i,minij;
	double sum=0.0,expctd,*sumi,*sumj,temp;

        MEW(sumi,(ni)+2,double);            /* = Ni. values */
        MEW(sumj,(nj)+2,double);            /* = N.j values */
	nni=ni;
	nnj=nj;
	for (i=1;i<=ni;i++) {
		sumi[i]=0.0;
		for (j=1;j<=nj;j++) {
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
		if (sumi[i] == 0.0) --nni;
	}
	for (j=1;j<=nj;j++) {
		sumj[j]=0.0;
		for (i=1;i<=ni;i++) sumj[j] += nn[i][j];
		if (sumj[j] == 0.0) --nnj;
	}
	*df=nni*nnj-nni-nnj+1;
	*chisq=0.0;
	for (i=1;i<=ni;i++) {
		for (j=1;j<=nj;j++) {
			expctd=sumj[j]*sumi[i]/sum+1;
			temp=nn[i][j]-expctd;
#if 0
if(nn[1][1] > 3000){
fprintf(stderr,"Obs[%d][%d]=%d; Exp=%g; temp=%g\n",i,j,nn[i][j],expctd,temp);
}
#endif
			*chisq += temp*temp/(expctd+TINY);
		}
	}
// if(nn[1][1] > 3000) fprintf(stderr,"chisq=%g\n",*chisq);
	if(*chisq < 0.0 || *df <= 0.0) *prob=1.0; 
	else *prob=gammq(0.5*(*df),0.5*(*chisq));
	minij = nni < nnj ? nni-1 : nnj-1;
	*cramrv=sqrt(*chisq/(sum*minij));
	*ccc=sqrt(*chisq/(*chisq+sum));
	free(sumi);
	free(sumj);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

#if 0
Int4	CardSset(sst_typ sst)
{
}

#endif

Int4    sqd_typ::expected( )
/*      E[i][j] = M[0][i]*M[1][j]/N*N  */
{
        double  M[3][3],N;	// M = margins...
        Int4    i,j;

        for(N=0.0,i=1;i<=2;i++) {                /* Get the row totals */
                for(M[1][i]=0.0,j=1;j<=2;j++) M[1][i] += cell[i][j];
                N += M[1][i];                   /* Get total number */
        }
        for(j=1;j<=2;j++) {
                for(M[2][j]=0.0,i=1;i<=2;i++) M[2][j]+= cell[i][j];
        }       /*** Calculate frequencies for each cell */
        if(N)for(i=1;i<=2;i++)for(j=1;j<=2;j++) E[i][j]=((M[1][i]*M[2][j])/N);
        return  (Int4) N;
}

BooLean	sqd_typ::TableItem(unsigned short i,unsigned short j,double Fisher_pcut)
#if 0	//***********************************************************************
	// Take intersect and determine the number of each type of
	// residue and the fraction of ri (or rj) for each of these
	// types.  Use (blosum62L[ri][r]/(blosum62L[ri][r]+1.0)
	// Set[j][r] == # r residues at j.
		   // IntersectSet1(set_typ B1, set_typ B2,set_typ IB);
		   IntersectSet1(SetI,SetJ,IB);
		   for(r=1; r <= nAlpha(AB); r++){
		      if(r == ri){
		         V += CardSet(IB)*(blosum62L[ri][r]/(blosum62L[ri][r]+1.0));
		      } else if(MemSset(r,sstI)){
		         IntersectSet1(IB,Set[j][r],IB);
		         V += CardSet(IB)*(blosum62L[ri][r]/(blosum62L[ri][r]+1.0));
			
		      } 
		   } 
	// Alternatively compute each intersect for all ris, rjs combinations.
	// Weigh these by
#endif	//***********************************************************************
{
	t_type	BestTab=0;
	double	BestInfo,BestPrb;
	Int4	r,ri,rj;
	sst_typ	sstI,sstJ;
	sst_typ	bestsstI,bestsstJ;
	double 	chisq,df,prob,cramrv,ccc;
	double	keyFrqI,keyFrqJ;	// minFrq
	set_typ	BestSetI=0,BestSetJ=0;
	unsigned char ris,rjs,set_ij[2],best_ij[2];
	cti_typ	*cti=0;

	assert(i > 0 && j > 0 && i <= Length && j <= Length);
	// Do either i or j contain too many gaps.
	ccc= (double) NumGaps[i]/(double)NumSeqs;
	if(ccc > MaxGapFrq) return FALSE;	
	ccc= (double) NumGaps[j]/(double)NumSeqs;
	if(ccc > MaxGapFrq) return FALSE;	

	ri = ResSeq(i,keyE); rj = ResSeq(j,keyE);
#if 1	// Temporary cth_typ
	Int4 hpsz0=NumResSets[ri]*NumResSets[rj]+5;
	cth_typ *cth0=new cth_typ(hpsz0);
#endif
	BestInfo=-99999;
#if 0	// First test for marginally significant correlation between 
	// single related residues.
	
#endif
	for(ris=1; ris <= NumResSets[ri]; ris++){
		Int4	V;
		set_ij[0]=ris;
	        sstI=ResidueSSet[ri][ris];;
		// 1. obtain set of all sequences with residues in sstI at i.
		ClearSet(SetI);
		BooLean	UselessSet = FALSE;
		for(keyFrqI=0.0,r=1; r <= nAlpha(AB); r++){
		   if(MemSset(r,sstI)){
			if(!Set[i][r]){
				fprintf(stderr,"i=%d; r=%d (%c)\n", i,r,AlphaChar(r,AB));
				print_error("sqd_typ error1");
			}
			if(r != ri && CardSet(Set[i][r]) < 3){ UselessSet=TRUE; break; }
			// Want an average of > 1 count per cell.
			UnionSet(SetI,Set[i][r]);
			keyFrqI += keyFreq[i][r];
		   }
		}
		if(UselessSet) continue;	// residue not found at this position.
		if(keyFrqI < MinKeyFrq) continue; // too few of these in subfamily.
		IntersectNotSet(USet,SetI,NotSetI);
		IntersectNotSet(NotSetI,Set[i][0]); // Remove 'X's from NotSetI
	        for(rjs=1; rjs <= NumResSets[rj]; rjs++){
		    set_ij[1]=rjs;
		    Tab=table[ri][rj][ris][rjs];
	            sstJ=ResidueSSet[rj][rjs];;
		    // 1. Find useless sets...
		    UselessSet = FALSE;
		    {
		      Int4 rJ,n1,n2;
		      for(rJ=1; rJ <= nAlpha(AB); rJ++){
		        if(rJ==rj || !MemSset(rJ,sstJ)) continue;

			// Make sure rJ contributes:
			if(CardSet(Set[j][rJ]) < 3){ UselessSet=TRUE; break; }
			// Want an average of > 1 count per cell.
		    	n1=CardInterSet(SetI,Set[j][rJ]); 
			if(n1 < 1){ UselessSet=TRUE; break; }
		    	n2=CardInterSet(NotSetI,Set[j][rJ]); 
		    	// n2=CardInterSetNotIJ(SetI,Set[j][rJ]); 
			// NOTE: something wrong with CardInterSetNotIJ()!!!
#if 0
			fprintf(stdout,"J(%c%d): n1 = %d; n2 = %d; debug\n",
				AlphaChar(rJ,AB),j,n1,n2);
#endif
			if(n1 < n2){ UselessSet=TRUE; break; }
		      }
		    }
		    if(UselessSet) continue;	// few related residues at this position.
		    // 2. Obtain set of all sequences with residues in sstJ at j.
		    ClearSet(SetJ);
		    for(keyFrqJ=0.0,r=1; r <= nAlpha(AB); r++){
		      if(MemSset(r,sstJ)){
			if(!Set[j][r]){
				fprintf(stderr,"j=%d; r = %d (%c)\n", j,r,AlphaChar(r,AB));
				print_error("sqd_typ error1");
			}
			UnionSet(SetJ,Set[j][r]);
			keyFrqJ += keyFreq[j][r];
		      }
		    }
		    if(keyFrqJ < MinKeyFrq) continue; // too few of these in subfamily.
		    IntersectNotSet(USet,SetJ,NotSetJ);
		    IntersectNotSet(NotSetJ,Set[j][0]); // Remove 'X's from NotSetJ
		    // Check sets...  // F82 vs S89 if(i==37 && j==56)
		    {
		      Int4 rI,rJ,n1,n2;
		      for(rI=1; rI <= nAlpha(AB); rI++){
			if(rI==ri || !MemSset(rI,sstI)) continue;
			// 1. make sure rI contributes:
		    	n1=CardInterSet(Set[i][rI],SetJ); 
		    	n2=CardInterSet(Set[i][rI],NotSetJ); 
		    	// n2=CardInterSetINotJ(Set[i][rI],SetJ); 
			// something wrong with CardInterSetINotJ()!!!  // now fixed...afn: 12_09_08
#if 0
			fprintf(stdout,"I(%c%d): n1 = %d; n2 = %d; debug\n",
				AlphaChar(rI,AB),i,n1,n2);
#endif
			if(n1 < n2){ UselessSet=TRUE; break; }
		      } if(UselessSet) continue;
		    }

		    // 3. Obtain sets for contingency table.
		    cell[1][1]=V=CardInterSet(SetI,SetJ); 
		    EqTab(0,0,ceil(V*0.25),Tab);

		    cell[2][1]=V=CardInterSet(SetI,NotSetJ);
		    EqTab(1,0,ceil(V*0.25),Tab);

		    cell[1][2]=V=CardInterSet(NotSetI,SetJ); 
		    EqTab(0,1,ceil(V*0.25),Tab);

		    cell[2][2]=V=CardInterSet(NotSetI,NotSetJ);
		    EqTab(1,1,ceil(V*0.25),Tab);
		    if(V == 0.0) continue;
		    this->expected( );
		    if(cell[1][1] <= E[1][1]) continue;

		    double I;
		    cntab1(cell,2,2,&chisq, &df, &prob, &cramrv, &ccc);
		    // this is ~50 x's faster than current Exact test routine!
		    I=chisq;
		    // fprintf(stderr,"I = %g; BestInfo = %g\n",I,BestInfo);
#if 1	// insert all significant hits....
		    if(prob > 0.0) IncdHist(-log10(prob),HG);
		    else IncdHist(9999999,HG);
		    if(prob <= Fisher_pcut){
#if 0	// FIX with x's...
		    if(i >= 210 || j >= 210){
			 fprintf(stderr,"%d(%d x's) vs %d(%d x's): intersect x's = %d\n",i,
				CardSet(Set[i][0]),j,CardSet(Set[j][0]),
				CardInterSet(Set[i][0],Set[j][0]));
			 PutTableMod(stderr,Tab);
		    }
#endif
			BestInfo=I; 
			t_type T=Table(0,SsetLet(ri),SsetLet(rj),AB); // create a 2x2 table.
			AddColumnsTab(0,sstI,sstJ,T);
			sst_typ	tmp_sst[2]; tmp_sst[0]=sstI; tmp_sst[1]=sstJ;
#if TEST_ResSetAll
			if(ResSetProb[i][ris] > prob){
				ResSetAll[i][ris] = sstI;
				ResSetProb[i][ris] = prob;
			}
			if(ResSetProb[j][rjs] > prob){
				ResSetAll[j][rjs] = sstJ;
				ResSetProb[j][rjs] = prob;
			}
#endif
			V = valTab(0,0,Tab); EqTab(0,0,V,T);
			V = valTab(1,0,Tab); EqTab(1,0,V,T);
			V = valTab(0,1,Tab); EqTab(0,1,V,T);
			V = valTab(1,1,Tab); EqTab(1,1,V,T);
			set_typ	tmp_set[2];
			tmp_set[0]=MakeSet(SetN(SetI)); CopySet(tmp_set[0],SetI);
		        tmp_set[1]=MakeSet(SetN(SetJ)); CopySet(tmp_set[1],SetJ);
			cti=new cti_typ(i,j,tmp_sst,tmp_set,set_ij,ri,rj,T,0,cell,4);
			assert(cti);
#if 0
			cell[1][1]-=CardInterSet(SetI,Set[i][AlphaCode('M',AB)]);
			cell[1][1]+=CardInterSet(NotSetI,Set[i][AlphaCode('M',AB)]);
			CardInterSet(SetI,Set[i][AlphaCode('M',AB)]);
#endif
			double prb=cti->ExactTest();
#if 1
			if(!cth0->Insert(prb,cti)){ delete cti; cti=0; } 
#else
			if(!cth->Insert(prb,cti)){ delete cti; cti=0; } 
#endif
		    }
#else
		    if(I > BestInfo){
			bestsstI = sstI; bestsstJ = sstJ;
			BestTab=Tab; BestInfo=I; 
			if(!BestSetI) BestSetI=MakeSet(SetN(SetI)); 
			CopySet(BestSetI,SetI);
			if(!BestSetJ) BestSetJ=MakeSet(SetN(SetJ)); 
			CopySet(BestSetJ,SetJ);
			best_ij[0]=set_ij[0]; best_ij[1]=set_ij[1];
			BestPrb=prob;
		    }
#endif
		}
	}
#if 1	// temporary heap
	// 1. Remove hits from the heap saving only those hits 
	BooLean	*UsedI,*UsedJ;
	NEW(UsedI,NumResSets[ri]+3,BooLean);	// have these nodes been assigned?
	NEW(UsedJ,NumResSets[rj]+3,BooLean);	// have these nodes been assigned?
	// Binary optimum matching of edges...
	while(!cth0->Empty( )){
	   cti=cth0->DeleteMin( );
	   // ris=cti->rsidI; rjs=UsedJ[cti->rsidJ];
	   ris=cti->rsidI; rjs=cti->rsidJ;
	   if(UsedI[ris] || UsedJ[rjs]){ // if already used...
		delete cti; cti=0;
	   } else {
		double prb=cti->ExactTest();
	   	if(!cth->Insert(prb,cti)){ delete cti; cti=0; } 
		else { UsedI[ris] = TRUE; UsedJ[rjs] = TRUE; }
	   }
	}
	free(UsedI); free(UsedJ);
	if(cth0) delete cth0;
#endif
	if(!BestTab) return FALSE;
	else return TRUE;
}

Int4	sqd_typ::TotalNodes()
{
	Int4	nodes=0,i;
	for(i=1; i <= Length; i++){
	   nodes+=NumResSets[ResSeq(i,keyE)];
	} return nodes;
}

double	sqd_typ::SearchSpace()
// Compute the size of the search space.
// only allow pairs
{
	double	size=0.0;
	Int4	i,j,ri,rj,min;
	for(i=1; i < Length; i++){
	   ri = ResSeq(i,keyE); 
	   for(j=i+1; j <= Length; j++){
		rj = ResSeq(j,keyE);
#if 0
		min = MINIMUM(Int4,NumResSets[ri],NumResSets[rj]);
		size += (double) min;	// only allow one pair between nodes.
#else
		size += (double) NumResSets[ri]*(double)NumResSets[rj];
#endif
	   }
	} return size;
}

double	ChiSqTable(double *prob,t_type T)
{
        Int4	*M[3];
	double 	chisq,df,cramrv,ccc;

	MEW(M[1],3,Int4);
	MEW(M[2],3,Int4);
	M[1][1] = T->cell[0][0];
	M[1][2] = T->cell[0][1];
	M[2][1] = T->cell[1][0];
	M[2][2] = T->cell[1][1];
	cntab1(M,2,2,&chisq, &df, prob, &cramrv, &ccc);
#if 0
	fprintf(stdout,"chisq=%g; df= %.2f; prob = %g; cramrv=%g; ccc=%g\n",
		chisq,df,prob,cramrv,ccc);
#endif
	free(M[2]); free(M[1]);
        return  chisq;
}

#define	TEST_CTH_TYPE	1

Int4	sqd_typ::MaxNumEdges()
// assuming binary matching of residue sets.
{
	Int4	size=0;
	Int4	i,j,ri,rj,min;
	for(i=1; i < Length; i++){
	   ri = ResSeq(i,keyE); 
	   for(j=i+1; j <= Length; j++){
		rj = ResSeq(j,keyE);
		min = MINIMUM(Int4,NumResSets[ri],NumResSets[rj]);
		size += min;	// only allow one pair between nodes.
	   }
	} return size;
}

char    **sqd_typ::CliqueStrings(Int4 MinClique,Int4 *NumClusters,double pcut,
		double Fisher_pcut)
{
	this->CliqueClusters(MinClique,NumClusters,pcut,Fisher_pcut);
	*NumClusters=NumSsetStr;
	return sst_str;
}
                
ctn_typ	*sqd_typ::CliqueClusters(Int4 MinClique,Int4 *NumClusters,double pcut,
	double Fisher_pcut)
{
	Int4	i,j,k,JunN=0,hpsz;
	char    str[300];

	HG=Histogram("-log10 of exact test tail probabilities",-10,300,5.0);
	Int4 r0q,r1q,s1,s0;

	// keyE=KeySeq( );
	// JunN=LenSeq(keyE);
	JunN=this->TotalNodes();
	Int4 minHeapCut;
	if(JunN > 0){	//**************************************************
#if 0	//******************************************************************
 Jun says:
 I found some interesting results regarding your question, but I hope to
 find some more delicate ones. The theorem states that for a random graph
 with n vertices and  each pair of vertices has probability p to become an
 edge (your  problem can be converted into this form to let p=M/(N choose
 2) and the results hold for both  cases), the size of the maximal clique
 of almost  all the graphs is d-1, d or d+1 where
 d=2*log(n)/log(1/p) + O(loglog n).

 In fact, d is the greatest natural number for which
 (n choose d) * p ^ (d choose 2)  >= 1.

#endif 	//******************************************************************
	   Int4	N=JunN,d,D,M,end;	// M == number of edges.
	   double	R,bestR;
	   Int4 best_d=0,minM,maxM;

// fprintf(stderr,"DEBUG 0\n");
#if 1
	   end=MaxNumEdges();
	   maxM=0; minM=end;
#else 
	   end = N*(N-1)/2;
	   maxM=0; minM= N*(N-1)/2;
#endif
// fprintf(stderr,"DEBUG 1\n");
	   assert(end > MinClique);
	   for(M=MinClique-1; M <= end; M++){
	     BooLean	Success=FALSE;
	     for(d=2; d <= N; d++){
#if 0
		double p = M/(double) MaxNumEdges();
#else
		double p = M/bico(N,2);
#endif
		R=bico(N,d)*pow(p,bico(d,2));
		if(R < 1.0) break;
		else { Success=TRUE; D = d; bestR=R; }
	     }
	     if(Success && D >= best_d){
		// fprintf(stderr,"N = %d; d = %d; M = %d; R = %g\n",N,D,M,bestR);
		best_d=D;
	     }
	     if(Success && D == MinClique-1){
		if(M < minM) minM=M;
		if(M > maxM) maxM=M;
	     } else if(D == MinClique) break;
	   }
// fprintf(stderr,"DEBUG 2\n");
#if 0
	   M = minM + (maxM-minM)/2;
	   hpsz=M;
	   minHeapCut=minM;
#else
	   M = minM;
	   hpsz=minM;		// make it harder to get an edge...
	   minHeapCut=minM;
#endif
#if 0	// Jun's new formula:
	   // P(cl(G_n) >= r)  <=  n^{d(n)-r}
	   // (2*pi)^{-1/2}* n^{n+.5}*(n-d)^{-n+d-.5}*d^{-d-.5}*p^{d(d-1)/2} =1
	   // where p = Edges/(N choose 2)
	   double d_n = pow((2.0*3.141592653589793238),-0.5);
	   d_n *= pow(N,N+ 0.5);
	   d_n *= pow(N-,N+ 0.5);
	   double P_clique = pow(N,d_n - MinClique);
#endif
	   // fprintf(stderr,"N = %d; d = %d; M = %d; hpsz = %d\n",N,MinClique-1,M,hpsz);
	   // exit(0);
	} //**************************************************

	if(cth) delete cth;
	cth = new cth_typ(hpsz);
	cti_typ *cti=0;
	for(s0=1; s0 <= LengthCMSA(1,cma); s0++){
	  // r0q = ResidueCMSA(1,1,s0,cma);
	  r0q = ResSeq(s0,keyE);
	  // fprintf(stderr,"%c%d \n",AlphaChar(r0q,AB),s0);
	  for(s1=s0+1; s1 <= LengthCMSA(1,cma); s1++){
		TableItem(s0,s1,Fisher_pcut);
	  }
	}

	// Clique analysis routines...
	Int4	setlim=2; // minimum size of clique to find
	Int4	N=LenSeq(keyE)+2; // vertex numbering is from 0.
	NumCnTab=cth->NumInHeap();
	// cti_typ	**CnTabList;
	assert(CnTabList==0);
	NEWP(CnTabList,NumCnTab +3,cti_typ);

	assert(ctn==0);
	ctn= new ctn_typ(LenSeq(keyE),max_residue_set,this->TotalNodes());
	Int4	*Rank,*Edges;
	double	*MaxCorrelate;
	sst_typ	*sstList;
	NEW(MaxCorrelate,ctn->NumNodes( ) +3,double);
	NEW(Rank,ctn->NumNodes( ) +3,Int4);
	NEW(Edges,ctn->NumNodes( ) +3,Int4);
	NEW(sstList,ctn->NumNodes( ) +3,sst_typ);
	Int4 res,rank;
	Int4	n1,n2;
	for(n1=1; n1 <= ctn->NumNodes( ); n1++){ MaxCorrelate[n1]=1.0; }

// cti->SetI()
	cti_typ	***edge;
#if 0
	NEWPP(edge,LenSeq(keyE)+3,cti_typ);
	for(i=1; i <= LenSeq(keyE); i++) NEWP(edge[i],LenSeq(keyE)+3,cti_typ);
#else
	NEWPP(edge,ctn->NumNodes( )+3,cti_typ);
	for(i=1; i <= ctn->NumNodes( ); i++) NEWP(edge[i],ctn->NumNodes( )+3,cti_typ);
#endif

	for(rank=1; !cth->Empty( ); rank++){
	   double prb = cth->MinKey( );
	   cti=cth->DeleteMin( );

	   r0q=cti->ResI; r1q=cti->ResJ;
	   s0=cti->SiteI; s1=cti->SiteJ;
	   n1=ctn->Node(cti->SiteI,cti->rsidI);
	   n2=ctn->Node(cti->SiteJ,cti->rsidJ);
	   edge[n1][n2]=edge[n2][n1]=cti;
	   if(sstList[n1]==0) sstList[n1]=cti->sstI;
	   if(sstList[n2]==0) sstList[n2]=cti->sstJ;
	   
#if 1
	   ctn->AddEdge(cti);
#else 
	   ctn->AddEdge(n1,n2);
#endif
	   CnTabList[rank]=cti;

	   if(Rank[n1] == 0){ Rank[n1]=rank; } 
	   if(Rank[n2] == 0){ Rank[n2]=rank; }
	   Edges[n1]++; Edges[n2]++;
	   MaxCorrelate[n1] = MINIMUM(double,MaxCorrelate[n1],prb);
	   MaxCorrelate[n2] = MINIMUM(double,MaxCorrelate[n2],prb);
#if 1	// only output of results: create verbose or Verbose option and Silent()
	if(table_fp){
	   sprintf(str,"(%d) %c%d vs %c%d: %.3f\n",rank,
		AlphaChar(r0q,AB),s0+OffSetSeq(keyE),
		AlphaChar(r1q,AB),s1+OffSetSeq(keyE),-log10(prb));
	   PutTableMod(table_fp,str,cti->Table()); 
	   cti->PutCells(table_fp);
	   double prob;
	   fprintf(table_fp,"ChiSq=%g\n",ChiSqTable(&prob,cti->Table()));
	   fprintf(table_fp,"----------------------------\n");
	}
#endif
	   // delete cti;
	}
	
// Clique analysis: multiple clique clusters...
  printf("Creating seed patterns.\n");

  NumClust;
  vst_typ **clique=ctn->CreateCliques(MinClique,hpsz,&NumClust,pcut);
  *NumClusters=NumClust;

  char  color[]="MGBCYR";
  Int4	numcolors=6;
#if 0	// Test Scores...
unsigned char *qseq;
NEW(qseq,LenSeq(keyE)+3,unsigned char);
#endif

  if(clique && NumClust > 0){
#if 1
    Int4	node,*nodeList,ptr_node;
    BooLean	*resList;
    NEW(nodeList,ctn->NumNodes()+3,Int4);
#endif
#if 0
    ctn->ConciseCliques(MinClique,1+NumCnTab/5,Rank); // use top 3rd to avoid overspecifying.
#endif
    NEWP(sst_str,NumClust+2,char);
    NumSsetStr=0;
    char tmp_str[100];
#if 0	// NEW: create query subfamily string...
    NEW(sst_str[0],20*4*LenSeq(keyE)+2,char);
    for(j=1; j <= Length; j++){
	sprintf(tmp_str,"%c%d,",AlphaChar(ResSeq(j,keyE),AB),j+OffSetSeq(keyE));
	strcat(sst_str[0],tmp_str);
    } j = strlen(sst_str[0]); sst_str[0][j-1]=0; // eliminate last ','.
    fprintf(stdout,"sst_str[0] = %s\n",sst_str[0]);
#endif
    for(Int4 n=1; n <= NumClust; n++){
     // clique[n]->Put(stdout,keyE,ctn->MaxNumResSets( ),AB); 
     if(clique[n]->Size() >= MinClique + 2){
	NumSsetStr++;
        NEW(sst_str[NumSsetStr],20*4*LenSeq(keyE)+2,char);
	Int4   clr=n%numcolors;
	ptr_node=0;
        NEW(resList,LenSeq(keyE)+3,BooLean);
	for(j=0; j < clique[n]->Size(); j++){
	  node=clique[n]->Vertex(j);
	  res=ctn->Site(node);

	  nodeList[ptr_node] = node; ptr_node++;

#if 0	// Test Scores...
	  if(n==1) qseq[res] = ResSeq(res,keyE);
#endif

// if(Rank[node] <= minHeapCut)
	    SsetToStr(sstList[node],str,AB);
	    // if(2*Rank[node] > NumCnTab) fprintf(stdout,"#");

	    //sprintf(tmp_str,"%s%d,",str,res);	// bug fix
#if 1
	    Int4 realres = FakeToRealCMA(1,res,keycma);
	    sprintf(tmp_str,"%s%d,",str,realres);
#else
	    sprintf(tmp_str,"%s%d,",str,res+OffSetSeq(keyE));
#endif
	    strcat(sst_str[NumSsetStr],tmp_str);

#if 0
	    fprintf(stdout,"%c%d.%c // {%s} %d rank -> %.3f (%d edges)\n",
		AlphaChar(ResSeq(res,keyE),AB),
		res+OffSetSeq(keyE),color[clr],str,Rank[node],
		-log10(MaxCorrelate[node]),Edges[node]);
#endif
	    resList[res]=TRUE;
	} nodeList[ptr_node] = 0;
	free(resList);
	j = strlen(sst_str[NumSsetStr]); sst_str[NumSsetStr][j-1]=0; // eliminate last ','.

#if 1	// Intersections and unions...
        set_typ Inclusive,Exclusive;
        Inclusive=MakeSet(SetN(SetI));
        Exclusive=MakeSet(SetN(SetI));
	// Determine 
	FillSet(Exclusive);
        for(j=0; nodeList[j]!=0; j++){
         for(k=j+1; nodeList[k]!=0; k++){
	   cti=edge[nodeList[j]][nodeList[k]];
	   // assert(cti); // cti==0 means these were clustered from cliques
	   if(cti){
	     UnionSet(Inclusive,cti->SetI());
	     IntersectSet3(Exclusive,cti->SetI());
	   }
	 }
        }
#if 0
	fprintf(stdout,"Cluster %d: Intersect = %d; Union = %d\n",
		n,CardSet(Exclusive),CardSet(Inclusive));
#endif
	// 3 way cliques...
	Int4 ei,ej,ek,res_i,res_j;
        for(i=0; nodeList[i]!=0; i++){
	 ei=nodeList[i];
	 res_i=ctn->Site(ei);
         for(j=i+1; nodeList[j]!=0; j++){
	  ej=nodeList[j];
	  res_j=ctn->Site(ej);
	  if(edge[ei][ej]){
	    if(ei < ej) edge[ei][ej]->StrResSetI(str,AB);
	    else edge[ei][ej]->StrResSetJ(str,AB);
#if 0
	    fprintf(stdout,"%c%d(%d){%s} && ",
		AlphaChar(ResSeq(res_i,keyE),AB),res_i,CardSet(edge[ei][ej]->SetI()),str);
#endif

	    if(ei < ej) edge[ei][ej]->StrResSetJ(str,AB);
	    else edge[ei][ej]->StrResSetI(str,AB);
#if 0
	    fprintf(stdout,"{%s}%c%d(%d)= %d\n", str,
		AlphaChar(ResSeq(res_j,keyE),AB),res_j,CardSet(edge[ei][ej]->SetJ()),
		CardInterSet(edge[ei][ej]->SetI(),edge[ei][ej]->SetJ()));
#endif
	    // edge[ei][ej]->PutTable(stdout);
	  }
          for(k=j+1; nodeList[k]!=0; k++){
	   ek=nodeList[k];
	   if(edge[ei][ej] && edge[ei][ek] && edge[ej][ek]){
	     
#if 0
	     FillSet(Exclusive);
	     ClearSet(Inclusive);
	     UnionSet(Inclusive,edge[ei][ej]->SetI());
	     UnionSet(Inclusive,edge[ei][ek]->SetI());
	     UnionSet(Inclusive,edge[ej][ek]->SetI());
	     IntersectSet3(Exclusive,edge[ei][ej]->SetI());
	     IntersectSet3(Exclusive,edge[ei][ek]->SetI());
	     IntersectSet3(Exclusive,edge[ej][ek]->SetI());
	     fprintf(stdout,"%c%d(%d,%d:%d) U %c%d(%d,%d:%d) U %c%d(%d,%d:%d)=%d\n\
\tIntersect = %d\n",
		AlphaChar(ResSeq(ei,keyE),AB),ei+OffSetSeq(keyE),
		CardSet(edge[ei][ej]->SetI()),CardSet(edge[ei][ek]->SetI()),
		CardInterSet(edge[ei][ej]->SetI(),edge[ei][ek]->SetI()),
		AlphaChar(ResSeq(ej,keyE),AB),ej+OffSetSeq(keyE),
		CardSet(edge[ei][ej]->SetI()),CardSet(edge[ej][ek]->SetI()),
		CardInterSet(edge[ei][ej]->SetI(),edge[ej][ek]->SetI()),
		AlphaChar(ResSeq(ek,keyE),AB),ek+OffSetSeq(keyE),
		CardSet(edge[ei][ek]->SetI()),CardSet(edge[ej][ek]->SetI()),
		CardInterSet(edge[ei][ek]->SetI(),edge[ej][ek]->SetI()),
		CardSet(Inclusive), CardSet(Exclusive));
#endif
#if 0
	     edge[ei][ek]->PutTable(stdout);
	     edge[ej][ek]->PutTable(stdout);
#endif
	   }
	  }
	 }
        }
	NilSet(Exclusive); NilSet(Inclusive);
#endif

     } 
    }
    free(nodeList);
  } 
  for(i=1; i <= ctn->NumNodes( ); i++) free(edge[i]); free(edge);

#if 0	// Test scores...
   {
	h_type HGM=Histogram("Scores against main clique",-10,100,1.0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    // WARNING: Rank needs to be fixed!!!
	    double score=PseudoAlnScoreCMSA(qseq,sq,Rank,cma);
	    IncdHist(score,HGM);
	}
	PutHist(stdout,60,HGM); NilHist(HGM); free(qseq);
   }
#endif

	if(cth) delete cth; cth=0;
    	free(sstList);
	free(Edges); free(Rank);
	// PutHist(stdout,60,HG);
	NilHist(HG); HG=0;
	free(MaxCorrelate);
#if 1
	return ctn;
#else
	delete ctn;
	return clique;
#endif
}

#if 0 // moved to ctn_typ.h...
void	ctn_typ::ConciseCliques(Int4 MinClique,Int4 MaxRank,Int4 *Rank)
// Return clique of sites only...
{
  if(!clique) return;
  else {
    Int4	res,j,n,m,node;
    BooLean	found;
    NEWP(concise,NumCluster+3,vst_typ);
    for(m=0,n=1; n <= NumCluster; n++){
     found=FALSE;
     if(clique[n]->Size() >= MinClique){
	for(j=0; j < clique[n]->Size(); j++){
	  node=clique[n]->Vertex(j);
	  // res=this->Site(node);
	  if(Rank[node] > MaxRank) continue;
	  else {
		if(!found){
			m++; found=TRUE;
			concise[m]=new vst_typ(N+1);
		} concise[m]->Add(node);
	  }
	}
     }
    } NumConcise=m;
  }
}
#endif

#if 0	// OLD resurrected CODE
/**********************************************************************
CnTable() - Given a 2-dimensional contingency table in the form 
of a double array bin[1..ni][1..nj], this routine returns the 
chi-square 'chsq', the number of degrees of freedom 'df' and the 
significance level 'prob' (small values indicating a significant association).
    Assuming no associations between i and j the expected value 
can be derived by the formula for a two dimensional table...

        N(..) = N(i.)*N(.j)/N....                 for 2 dimensions      (1)

    where 
                ---                       ---
                \                         \
        N(.j) = /   N(ij),        N(i.) = /    N(ij), 
                ---                       ---
                 i                         j

    The number of degrees of freedom is equal to the number of entries
in the table - 1 (e.g., dim=2: ni*nj=4-1) minus the number of probabilities
estimated from the data for the particular hypothesis being tested.
Since we are here considering the hypothesis of mutual independence
there are (ni-1)+(nj-1) probabilites estimated from the
data.  Therefore 
        df = (ni*nj) - (ni+nj) + 1.  
**********************************************************************/
double  CnTable(t_type T,double *chsq, double *df)
{
        int     ni,nj,i,j,p,sg,sh;
        double  expctd,denom=1.0,N=0.0, *Ni,*Nj,temp;
        double  small=0.0,total=0.0,**bin;
        boolean okay=TRUE;

        bin = T->bin; ni=2; nj=2;
        NEW(Ni,(ni),double);            /* = Ni. values */
        NEW(Nj,(nj),double);            /* = N.j values */
        for(i=0;i<ni && okay;i++) {     /* Get the row totals */
                Ni[i]=0.0;
                for(j=0;j<nj;j++) Ni[i] += bin[i][j];
                N += Ni[i];                     /* Get total number */
                if(Ni[i] == 0.0) okay=FALSE;
        }
        denom *= N;
        for(j=0;j<nj && okay;j++) {
                Nj[j]=0.0;
                for(i=0;i<ni;i++)Nj[j]+=bin[i][j];
                if(Nj[j] == 0.0) okay=FALSE;
        }
        *chsq=0.0;
        for(i=0;i<ni;i++) {
                for(j=0;j<nj;j++) {
                        expctd=(Ni[i]*Nj[j])/(denom);
                        /****** TEST VALIDITY OF CHI SQUARED VALUE ******
                        total += 1.0;
                        if(expctd < 2.0) { okay = FALSE; break; }
                        if (expctd < 5.0) small += 1.0;
                        /****************** END OF TEST *****************/
                        temp=bin[i][j]-expctd;
                        *chsq += temp*temp/(expctd+TABSTINY); /* Here TINY */
                }       /* guarantees that any eliminated row or column */
        }               /* will not contribute to the N. */
        /****** TEST VALIDITY OF CHI SQUARED VALUE ******
        if(okay && small/total > 0.20) okay = FALSE;
        /****************** END OF TEST *****************/
        if(okay) {
           *df=(double)(ni*nj-(ni+nj)+1);
           free(Ni); free(Nj); 
           return gammq((double)(0.5*(*df)),(double)(0.5*(*chsq)));
        } else { free(Ni);free(Nj);return ILLEGAL; }
}
#endif

