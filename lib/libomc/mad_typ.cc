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

#include "mad_typ.h"
#include "blosum62.h"

#define	MaxResSets	16

void    mad_typ::Free()
{
	Int4	r,s;
	for(s=0; s <= Length; s++){
	   if(Set[s]){
              for(r=0; r<=nAlpha(AB); r++){
		if(Set[s][r]) NilSet(Set[s][r]);
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
	if(sst_str){
           for(s=0; s <= NumSsetStr; s++){ if(sst_str[s]) free(sst_str[s]); } free(sst_str);
	}
	delete ctn;
	delete rst;
	if(keyE) NilSeq(keyE);
	if(bstE) NilSeq(bstE);
	free(skip); free(NumGaps);
}

void	mad_typ::SubGrpCsqBsq(set_typ set)
{
	Int4	i,j,x,s,sq,Sq,**cnts,d,lenM,max,pos[3];
	double	likelihood,best,score;
	st_type	S=SitesCMSA(cma);
	a_type	AB=AlphabetCMSA(cma);

    	assert(nBlksCMSA(cma) == 1);
	NEW(skip,NumSeqsCMSA(cma) + 3, BooLean);
	for(i=1; i <= NumSeqsCMSA(cma); i++) if(!MemberSet(i,set)) skip[i]=TRUE;

	keyFreq=ColResFreqsCMSA(1,skip,&key_obs,cma);

        Int4 N,N0; N=N0=NumSeqsCMSA(cma);
        for(N=0,sq=1; sq <= N0; sq++) if(!skip[sq]) N++; 
	NumSeqs=N;	// # seqs in subset.
	assert(N == CardSet(set));

        lenM=LengthCMSA(1,cma);

	// Get csqE and WtFreq[i][r]
	char 	r,R,*seq; NEW(seq,lenM + 3,char);
	double  total,**WtFrq,**SqWt=swt->Weights();
	NEWP(WtFrq,lenM +4,double);
        for(i=1; i <= lenM; i++){
	   NEW(WtFrq[i],nAlpha(AB) +3, double);
	   for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	   	if(skip[sq]) continue;
#if 0		// This is equivalent to below iff x == 0 --> r == 0;
	   	e_type E=TrueSeqCMSA(sq,cma);
		x=TruePosCMSA(sq,i,cma); r=ResSeq(x,E);
		WtFrq[i][r] += SqWt[i][sq];
#else
	   	s=SitePos(1,sq,1,SitesCMSA(cma)); // start of fake sequence; first blk first repeat.
	   	unsigned char *fsq = SeqSeqSet(sq,DataCMSA(cma));
		r = fsq[s+i-1];
		WtFrq[i][r] += SqWt[i][sq];
#endif
	   } total=0.0;
#if 1
	   for(best=0.0,R=0,r=0; r <= nAlpha(AB); r++){
		if(best < WtFrq[i][r]){ best=WtFrq[i][r]; R=r; }
		total+=WtFrq[i][r];
	   } if(R == 0) seq[i]='A'; else seq[i]=AlphaChar(R,AB);
#else	// Pick Csq based on # observed residues, not on weighted frequencies...
	   Int4 bst;
	   for(bst=0,R=0,r=0; r <= nAlpha(AB); r++){
		if(bst < key_obs[i][r]){ bst=key_obs[i][r]; R=r; }
		total+=WtFrq[i][r];
	   } if(R == 0) seq[i]='A'; else seq[i]=AlphaChar(R,AB);
#endif
	   if(total > 0.0){
	     for(r=0; r <= nAlpha(AB); r++){ WtFrq[i][r] = WtFrq[i][r]/total; }
	     // fprintf(stderr,"best[%d] = %g (%c) (total = %g)\n",best,i,AlphaChar(R,AB),total);
	   } 
	} seq[i]=seq[0]=0; csqE=StringToSeq(seq+1,"Consensus seq", 1, AB); 
	// PutSeq(stderr,csqE,AB);

	// Get bstE 
	// h_type xHG=Histogram("Sequence Scores",0,10000,2);
        for(best=0.0,Sq=0,sq=1; sq <= NumSeqsCMSA(cma) ; sq++){
	   if(skip[sq]) continue;
#if 0
	   e_type E=TrueSeqCMSA(sq,cma);
           for(score=0.0,i=1; i <= lenM; i++){
		if((x=TruePosCMSA(sq,i,cma)) == 0) continue; 
		if((r=ResSeq(x,E)) > 0) score+=WtFrq[i][r];
	   } if(best < score){ best=score; Sq=sq; }
#else
	   s=SitePos(1,sq,1,SitesCMSA(cma)); // start of fake sequence.
	   unsigned char *fsq = SeqSeqSet(sq,DataCMSA(cma));
	   for(score=0.0,i=1; i<=lenM; i++,s++){
		r = fsq[s];
		if(r > 0) score+=WtFrq[i][r];
	   } if(best < score){ best=score; Sq=sq; }
#endif
	   // IncdHist(score,xHG);
	} // PutHist(stderr,60,xHG); NilHist(xHG);
	bstE=GetSeqAsCsqCMSA(Sq,cma);
	char Ala= AlphaCode('A',AB);
        for(i=1; i <= LenSeq(bstE); i++){
	   if(ResSeq(i,bstE) == 0){ EqSeq(i,Ala,bstE); EqXSeq(i,Ala,bstE); }
	} // PutSeq(stderr,bstE,AB);
        for(i=1; i <= lenM; i++){ free(WtFrq[i]); } free(WtFrq); free(seq);
}


// void	mad_typ::init(cma_typ CMA,double MinFreq,double max_gap_frq, char mode)
void	mad_typ::init(swt_typ *inswt, set_typ inset, cma_typ CMA,double MinFreq,double max_gap_frq, char mode)
{
	Int4	j,s,r,r1,r2,rs;
	sst_typ	sst;

	verbose=FALSE; // verbose=1;
	assert(CMA); cma=CMA; AB=AlphabetCMSA(cma); swt=inswt;
	ctn =0; NumClust=0; sst_str=0; NumSsetStr=0; 
	CnTabList=0; NumCnTab=0; cth=0; HG=0;

	SubGrpCsqBsq(inset); keyE=csqE;  // skip & NumSeqs initialized here!

	assert(MinFreq >= 0.0 && MinFreq <= 1.0);
	MinKeyFrq=MinFreq; MaxGapFrq=max_gap_frq;
	assert(MaxGapFrq >= 0.0 && MaxGapFrq < 1.0);
	Tab=0; NEW(cell[1],3,Int4); NEW(cell[2],3,Int4);
	Length=LenSeq(keyE);
	if(Length != LengthCMSA(1,cma)){
		fprintf(stderr,"LenSeq(keyE) = %d; LengthCMSA = %d\n", Length,LengthCMSA(1,cma));
		assert(Length == LengthCMSA(1,cma));
		// print_error("Fatal: length inconsistency");
	}
	rst=new rst_typ(mode);
	ResidueSSet = rst->ResidueSSet; NumResSets= rst->NumResSets;
	max_residue_set = rst->MaxResSet();
	// rst->Put(stderr);

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
	{
	  Int4 sq,s0,Sq;
	  unsigned char *qsq=SeqPtr(keyE);
	  // for(Sq=0,sq=2; sq <=NumSeqsCMSA(cma); sq++)
	  for(Sq=0,sq=1; sq <=NumSeqsCMSA(cma); sq++) {
	     if(skip[sq]) continue; else Sq++;
#if 0	// problems with this !?!? Due to gap '-' residues positions in TrueSeq.
	     for(s=1; s<=Length; s++){
		e_type	E=TrueSeqCMSA(sq,cma);
	     	s0=TruePosCMSA(sq,s,cma); rs=ResSeq(s0,E);
if(s0 == 0){	// for gaps!!!
	e_type fE=FakeSeqCMSA(sq,cma);
	AlnSeqSW(12, 4,E,fE, AB);
	assert(s0 > 0);
}
	     	sst = ResidueSSet[qsq[s]][0];
		if(MemSset(rs,sst)) AddSet(Sq,Set[s][rs]);
		else if(rs==0) AddSet(Sq,Set[s][0]);
	     }
#else	// This seems to work!
	     s0=SitePos(1,sq,1,SitesCMSA(cma)); // start of fake sequence blk==1; rpt==1.
	     unsigned char *fsq = SeqSeqSet(sq,DataCMSA(cma));
	     for(s=1; s<=Length; s++,s0++){
		rs = fsq[s0];
	     	sst = ResidueSSet[qsq[s]][0];
		if(rs==0) AddSet(Sq,Set[s][0]);
		else if(MemSset(rs,sst)) AddSet(Sq,Set[s][rs]);
	     }
#endif
	  } assert(Sq == NumSeqs);
	} 
	NEW(NumGaps,Length+3,Int4);
	for(s=1; s<=Length; s++) NumGaps[s] = CardSet(Set[s][0]);
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
		  table[r1][r2][rs1][rs2]=T;
		  // PutTableMod(stderr,table[r1][r2][rs1][rs2]);
		  total_tab++;
		}
	     }
	   }
	} // fprintf(stderr,"\t%d tables made.\n",total_tab);
	sst_typ *sst_ptr;
	// rst->Put(stderr);	
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
		for (j=1;j<=nj;j++) { sumi[i] += nn[i][j]; sum += nn[i][j]; }
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

Int4    mad_typ::expected( )
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

BooLean	mad_typ::TableItem(unsigned short i,unsigned short j,double Fisher_pcut)
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
				print_error("mad_typ error1");
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
#if 0	// DEBUG
			if(Set[j][rJ] == 0){
				fprintf(stderr," r = %c: ",AlphaChar(rj,AB));
				PutSST(stderr,sstJ,AB);
				fprintf(stderr," --> Set[%d][%c] non-existent\n",j,AlphaChar(rJ,AB));
			}
#endif
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
				print_error("mad_typ error1");
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
		    EqTab(0,0,ceil(V*0.25),Tab);	// divide by 0.25 to avoid overflow; use as key only.

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
			BestInfo=I; 
			t_type T=Table(0,SsetLet(ri),SsetLet(rj),AB); // create a 2x2 table.
			AddColumnsTab(0,sstI,sstJ,T);
			sst_typ	tmp_sst[2]; tmp_sst[0]=sstI; tmp_sst[1]=sstJ;
			V = valTab(0,0,Tab); EqTab(0,0,V,T);
			V = valTab(1,0,Tab); EqTab(1,0,V,T);
			V = valTab(0,1,Tab); EqTab(0,1,V,T);
			V = valTab(1,1,Tab); EqTab(1,1,V,T);
			set_typ	tmp_set[2];
			tmp_set[0]=MakeSet(SetN(SetI)); CopySet(tmp_set[0],SetI);
		        tmp_set[1]=MakeSet(SetN(SetJ)); CopySet(tmp_set[1],SetJ);
			cti=new cti_typ(i,j,tmp_sst,tmp_set,set_ij,ri,rj,T,0,cell,4);
			assert(cti);
			double prb=cti->ExactTest();	// prb is used as a key only...
			if(!cth0->Insert(prb,cti)){ delete cti; cti=0; } 
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
	// 1. Remove hits from the heap saving only the most significant patterns for each i,j.
	BooLean	*UsedI,*UsedJ;
	NEW(UsedI,NumResSets[ri]+3,BooLean);	// have these nodes been assigned?
	NEW(UsedJ,NumResSets[rj]+3,BooLean);	// have these nodes been assigned?
	// Binary optimum matching of edges...
	while(!cth0->Empty( )){
	   cti=cth0->DeleteMin( );
	   ris=cti->rsidI; rjs=cti->rsidJ;
	   if(UsedI[ris] || UsedJ[rjs]){ // if already used...
		delete cti; cti=0;
	   } else {
		double prb=cti->ExactTest();
	   	if(!cth->Insert(prb,cti)){ delete cti; cti=0; } 
		else { UsedI[ris] = TRUE; UsedJ[rjs] = TRUE; }
	   }
	} free(UsedI); free(UsedJ);
	if(cth0) delete cth0;
#endif
	if(!BestTab) return FALSE;
	else return TRUE;
}

Int4	mad_typ::TotalNodes()
{
	Int4	nodes=0,i;
	for(i=1; i <= Length; i++){
	   nodes+=NumResSets[ResSeq(i,keyE)];
	} return nodes;
}

double	mad_typ::SearchSpace()
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

Int4	mad_typ::MaxNumEdges()
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

char    **mad_typ::CliqueStrings(Int4 MinClique,Int4 *NumClusters,double pcut,
		double Fisher_pcut,sst_typ **&xsst)
{
	this->CliqueClusters(MinClique,NumClusters,pcut,Fisher_pcut/SearchSpace(),xsst);
	*NumClusters=NumSsetStr;
	return sst_str;
}
                
ctn_typ	*mad_typ::CliqueClusters(Int4 MinClique,Int4 *NumClusters,double pcut,
	double Fisher_pcut,sst_typ **&sst_rtn)
{
	Int4	i,j,k,JunN=0,hpsz;
	char    str[300];

	HG=Histogram("-log10 of exact test tail probabilities",-10,300,5.0);
	Int4 r0q,r1q,s1,s0;

  	sst_rtn=0; JunN=this->TotalNodes();
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
	for(s0=1; s0 <= Length; s0++){
	  r0q = ResSeq(s0,keyE);
	  // fprintf(stderr,"%c%d \n",AlphaChar(r0q,AB),s0);
	  for(s1=s0+1; s1 <= Length; s1++){
		this->TableItem(s0,s1,Fisher_pcut);
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
	NEWPP(edge,ctn->NumNodes( )+3,cti_typ);
	for(i=1; i <= ctn->NumNodes( ); i++) NEWP(edge[i],ctn->NumNodes( )+3,cti_typ);
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
	   ctn->AddEdge(cti);
	   CnTabList[rank]=cti;

	   if(Rank[n1] == 0){ Rank[n1]=rank; } 
	   if(Rank[n2] == 0){ Rank[n2]=rank; }
	   Edges[n1]++; Edges[n2]++;
	   MaxCorrelate[n1] = MINIMUM(double,MaxCorrelate[n1],prb);
	   MaxCorrelate[n2] = MINIMUM(double,MaxCorrelate[n2],prb);
	if(verbose){ // output results: create verbose or Verbose option and Silent()
#if 0
	   fprintf(stderr,"(%d) %c%d vs %c%d: %.3f\n",rank,
		AlphaChar(r0q,AB),s0+OffSetSeq(keyE),
		AlphaChar(r1q,AB),s1+OffSetSeq(keyE),-log10(prb));
#else
	   sprintf(str,"(%d) %c%d vs %c%d: %.3f\n",rank,
		AlphaChar(r0q,AB),s0+OffSetSeq(keyE),
		AlphaChar(r1q,AB),s1+OffSetSeq(keyE),-log10(prb));
	   PutTableMod(stderr,str,cti->Table()); 
	   cti->PutCells(stderr);
	   double prob;
	   fprintf(stderr,"ChiSq=%g\n",ChiSqTable(&prob,cti->Table()));
	   fprintf(stderr,"----------------------------\n");
#endif
	}
	   // delete cti;
	}
	
// Clique analysis: multiple clique clusters...
  fprintf(stderr,"Creating seed patterns.\n");

  NumClust=0;
  FILE *efp=0; if(verbose) efp=stderr;
  vst_typ **clique=ctn->CreateCliques(efp,MinClique,hpsz,&NumClust,pcut);
  *NumClusters=NumClust;

  char  color[]="MGBCYR";
  Int4	numcolors=6;
    NumSsetStr=0;
#if 0	// Test Scores...
unsigned char *qseq;
NEW(qseq,LenSeq(keyE)+3,unsigned char);
#endif

  if(clique && NumClust > 0){
    Int4	node,*nodeList,ptr_node;
    NEW(nodeList,ctn->NumNodes()+3,Int4);
#if 0
    ctn->ConciseCliques(MinClique,1+NumCnTab/5,Rank); // use top 3rd to avoid overspecifying.
#endif
    NEWP(sst_str,NumClust+2,char);
    NEWP(sst_rtn,NumClust+2,sst_typ);
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
        NEW(sst_rtn[NumSsetStr],LenSeq(keyE)+2,sst_typ);
	Int4   clr=n%numcolors;
	ptr_node=0;
#if 0
	for(j=0; j < clique[n]->Size(); j++){
	  node=clique[n]->Vertex(j);
	  res=ctn->Site(node);
	  nodeList[ptr_node] = node; ptr_node++;
	  // if(Rank[node] <= minHeapCut)
	  SsetToStr(sstList[node],str,AB);
	  // if(verbose) assert(sst_rtn[NumSsetStr][res] == 0);	// duplicate sites are showing up; fix this!!
	  sst_rtn[NumSsetStr][res] = sstList[node]; // residue set at position 'res'.
	  sprintf(tmp_str,"%s%d,",str,res);
	  strcat(sst_str[NumSsetStr],tmp_str);
#if 0
	  fprintf(stdout,"%c%d.%c // {%s} %d rank -> %.3f (%d edges)\n",
		AlphaChar(ResSeq(res,keyE),AB),
		res+OffSetSeq(keyE),color[clr],str,Rank[node],
		-log10(MaxCorrelate[node]),Edges[node]);
#endif
	}
#else	//
	for(j=0; j < clique[n]->Size(); j++){
	  node=clique[n]->Vertex(j);
	  res=ctn->Site(node);
	  nodeList[ptr_node] = node; ptr_node++;
	  sst_rtn[NumSsetStr][res] = sstList[node]; // residue set at position 'res'.
	}
	for(j=1; j <= LenSeq(keyE); j++){
	  if(sst_rtn[NumSsetStr][j]){
		SsetToStr(sst_rtn[NumSsetStr][j],str,AB);
	  	sprintf(tmp_str,"%s%d,",str,j);
	  	strcat(sst_str[NumSsetStr],tmp_str);
	  }
	}
#endif
	nodeList[ptr_node] = 0;
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
	    if(ei < ej) edge[ei][ej]->StrResSetJ(str,AB);
	    else edge[ei][ej]->StrResSetI(str,AB);
	  }
          for(k=j+1; nodeList[k]!=0; k++){
	   ek=nodeList[k];
	   if(edge[ei][ej] && edge[ei][ek] && edge[ej][ek]){
	     
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

	if(cth) delete cth; cth=0;
    	free(sstList); free(Edges); free(Rank);
	if(verbose) PutHist(stderr,60,HG);
	NilHist(HG); HG=0;
	free(MaxCorrelate);
	return ctn;
}


