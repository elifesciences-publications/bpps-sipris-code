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

#include "omc_typ.h"
#include "hmm_typ.h"

//================ for developing & testing 'B' option ====================

double  LogCumHypGeomProb(double N1,double N2, double n,double x)
#if 0   //****************************************************
N total balls with N1 red balls and N2 = N-N1 black balls.
Choose n balls at random.  The probability that the
group so chosen will contain x or more red balls is
given by: p=CumHypGeomProb(N1,N2,n,x).
#endif  //****************************************************
{
        double  p,K,end;

        if(x == 0) return 1.0;
        end = MINIMUM(double,N1,n);
        K = (lngamma(N1+1)+lngamma(N2+1)-lngamma(N2+N1+1)+lngamma(n+1)+lngamma(N2+N1-n+1));
        for(p=0.0; x <= end; x += 1.0){
           p += K-lngamma(x+1)-lngamma(N1-x+1)-lngamma(n-x+1)-lngamma(N2-n+x+1);
        } return p;
}

double	***BoltzScore;	// BoltzScore[n][j][r];
double	***BoltzStd;	// BoltzStd[n][j][r];
double	**Entropy;	// Entropy[n][j];

double	bpps_typ::PutBoltzmannLike(FILE *fp,double *StdFrq, Int4 j, UInt4 **CntBG,UInt4 **CntFG, Int4 n)
// for saving information on pattern residues and 
// read this in using csp->ReadBinomialBPPS(-cutoff,res_evals);
{
	Int4 i,r;
	FILE	*efp=0; // efp=stderr;
	double N,M,p,q,n1,n2,m1,m2,P,StdP;
//	fprintf(fp,"#             __Foreground___    __Background___     ____Information____    WtNumSeq\n");
//	fprintf(fp,"#   Pattern:  Match  Diverge     Match  Diverge \n");
	 
	double	map=0.0,map0=0.0,RefSubJ,StdRefSubJ,StdSubJ,d,D,DD,dd,RefMap;
	sst_typ QstJ=qst[j];
	unsigned char QueryJ=query[j];
	{ qst[j]=0; query[j]=0; }
	SubLPR(CntBG,CntFG); map0=subLPR[0];

#if 1	// Get standard and reference FG and BG counts.
	UInt4 **RefCntFG,TotalCntFG,TotalCntBG,**StdCntFG,**StdCntBG;
	NEWP(RefCntFG,k+3,UInt4); NEWP(StdCntFG,k+3,UInt4); NEWP(StdCntBG,k+3,UInt4);
	// compute total FG and BG counts at position j.
	for(TotalCntFG=0,TotalCntBG=0,r=0; r <= nAlpha(AB); r++)
		{ TotalCntFG += CntFG[j][r]; TotalCntBG += CntBG[j][r]; }
	for(i=1; i <= k; i++){
	    if(i != j){	// for positions != j use the regular counts (i.e., we want the difference in LLR at j).
		StdCntBG[i]=CntBG[i]; RefCntFG[i]=CntFG[i]; StdCntFG[i]=CntFG[i]; continue; 
	    } // for jth position compute...
	    NEW(RefCntFG[j],nAlpha(AB)+3,UInt4);
	    NEW(StdCntFG[j],nAlpha(AB)+3,UInt4); NEW(StdCntBG[j],nAlpha(AB)+3,UInt4);
	    for(r=0; r <= nAlpha(AB); r++){
		// m1=(double) CntBG[j][r]*wt_factor; M=(double) TotalCntBG*wt_factor;
		m1=(double) CntBG[j][r]; M=(double) TotalCntBG;
		d = m1/M;	// = the mean of beta distribution.
		D = floor( (d*(double)TotalCntFG) + 0.5);  // make FG[r] proportional to BG[r] count...
		if(D < 0.0) D = 0.0; RefCntFG[j][r]= (UInt4) D;

		d = StdFrq[r];		// = the overall frequency in the input alignment.
		D = floor( (d*(double)TotalCntFG) + 0.5); 
		if(D < 0.0) D = 0.0; StdCntFG[j][r]= (UInt4) D;

		D = floor( (d*(double)TotalCntBG) + 0.5); 
		if(D < 0.0) D = 0.0; StdCntBG[j][r]= (UInt4) D;

		// fprintf(fp,"%c=%.1f; ",AlphaChar(r,AB),wt_factor*(double)RefCntFG[j][r]);
	    } // fprintf(fp,"\n  Ref=%d; TotalFG=%d\n",C,TotalCntFG);
	}
#endif
	double	max=0,min=999999,entropy=0;
	double	maxFrq=0,minFrq=0;
	unsigned char maxR=0,minR=0;
	char	PttrnPos[3]; PttrnPos[0]=0;
	// for(unsigned char R=0; R <= nAlpha(AB); R++)
	for(unsigned char R=1; R <= nAlpha(AB); R++)
	{
	   if(R==0){ qst[j] = QstJ; query[j]=QueryJ; }
	   else { qst[j] = SsetLet(R); query[j]=R; } 

#if 1
	   SubLPR(StdCntBG,CntFG); RefMap=subLPR[0]; StdSubJ=subLPR[j];
	   SubLPR(StdCntBG,StdCntFG); RefMap=subLPR[0]; StdRefSubJ=subLPR[j];

	   SubLPR(CntBG,RefCntFG); RefMap=subLPR[0]; RefSubJ=subLPR[j];
#endif
	   SubLPR(CntBG,CntFG); map=subLPR[0];

	   m1 = MatchBG[j]; m2 = MisMatchBG[j]; M=m1+m2; p=(m1+1)/(M+2);
	   n1=MatchFG[j];  n2=MisMatchFG[j]; N=n1+n2; q=(n1+1)/(N+2);
	   // if(n1 <= 0.0) continue;
	   p = n1/N;	// p = probability of match in foreground...
	   if(n1 > 0.0) entropy -= p*log(p); else n1=0.001;

#if 1
	   StdP = StdSubJ - StdRefSubJ; StdP=10*StdP/N;
	   P = subLPR[j] - RefSubJ; P=10*P/N;
#elif 1
	   P = subLPR[j]/(wt_factor*N);
#else
	   P = (map-map0)/(wt_factor*N);
#endif
#if 0
	   // P = LnCBP(n1,N,q) - LnCBP(n1,N,p); P = 100*P/N;
	   double Red=m1+n1,Black=m2+n2,mode,mean;
	   mean = (N+M)*Red/(Red+Black);
	   mode = n1*(N+M+1)*(Red+1)/(Red+Black+1);
	   dd=LogCumHypGeomProb(Red,Black,N,n1) - LogCumHypGeomProb(Red,Black,M,mode);
	   // dd=map-RefMap;
	   DD=-log((n1/N)/(m1/M));
#endif
	   if(efp && qst[j]){	// if a pattern position?  This is only needed for printout...
	     char tmp[30],*ptr; ptr=tmp;
	     for(r=StartAlpha; r <= nAlpha(AB); r++){
		if(MemSset(r,qst[j])){ sprintf(ptr,"%c",AlphaChar(r,AB)); ptr++; }
	     }
	     if(R == 0){
		PttrnPos[0]='!'; PttrnPos[1]=0;
	        sprintf(ptr,"%d",j);
	        fprintf(fp,"%12s: %5.0f %6.0f (%2.0f%c) %5.0f %6.0f (%2.0f%c)    %5.1f  %5.1f   %5.1f",
		  tmp,n1,n2, 100.0*n1/N,'%',m1,m2, 100.0*m1/M,'%', subLPR[j],P,StdP); 
	        if(R > 0) fprintf(fp,"      %5.1f\n", M+N);
	        else fprintf(fp,"      %5.1f\n---------\n", M+N);
	     } else {
	       fprintf(fp,"%s %.0f %.3f %.1f %.2f %.3f\n",tmp,n1,100*n1/N,-P,-StdP,100*m1/M); 
	     }
	   } // D=m1/M; d=n1/N;
#if 1
	   if(!isfinite(P)) continue;
	   if(R > 0){
		BoltzScore[n][j][R]=P; BoltzStd[n][j][R]=StdP; 
	        fprintf(stderr,"BoltzScore[%d][%d][%c] = %.3f\n",n,j,AlphaChar(R,AB),P); 
	        fprintf(stderr,"BoltzStd[%d][%d][%c] = %.3f\n",n,j,AlphaChar(R,AB),StdP); 
		if(-P > max){ max = -P;  maxR=R; maxFrq=100*m1/M; }
		if(-P < min){ min = -P;  minR=R; minFrq=100*n1/N; }
	   }
#else
	   if(!isfinite(P)) continue;
	   else if(R > 0){ 
		BoltzScore[n][j][R]=P; BoltzStd[n][j][R]=StdP; 
		if(P > max){ max=P;  maxR=R; minR=R; maxFrq=n1/N; minFrq=m1/M; }
		// if(P > max){ max = P;  maxR=R; maxFrq=100*n1/N; }
		// if(P < min){ min = P;  minR=R; minFrq=100*n1/N; }
	   }
#endif
	}
	if(efp) fprintf(fp,"Range(%d)%s=%.1f: %c%d(%.2f) --> %c%d(%.2f); entropy=%.3f\n\n",
		n,PttrnPos,max,AlphaChar(maxR,AB),j,maxFrq,AlphaChar(minR,AB),j,minFrq,entropy);
	Entropy[n][j]=entropy;
	qst[j] = QstJ; query[j]=QueryJ;	// return this to original state.
#if 1
	free(RefCntFG[j]); free(RefCntFG); 
	free(StdCntFG[j]); free(StdCntFG); free(StdCntBG[j]); free(StdCntBG); 
	return max-min;
#else
	// return max;
	return max*maxFrq;
#endif
}

double	che_typ::PutBoltzmannLike(FILE *fptr,double *frq,Int4 c, Int4 n)
{
        return pps->PutBoltzmannLike(fptr,frq,c,ResWt[2],ResWt[1],n);
}

double	mcs_typ::PutBoltzmannLike(FILE *fp,Int4 n,Int4 c)
{
	fprintf(stdout,"====== %d: %s (column %d) ========\n",n,Hpt->ElmntSetName(n),c);
	double  *freq = tFreqSeqSet(TrueDataCMSA(TrueMainCMA));
	// double  *freq = tFreqSeqSet(DataCMSA(TrueMainCMA));
        return che[n]->PutBoltzmannLike(fp,freq,c,n); 
}

double  GetSqProbCMSA(cma_typ cma, Int4 n, cma_typ CMA)
{
    assert(TotalLenCMSA(cma) == TotalLenCMSA(CMA));
    assert(nBlksCMSA(cma) == 1); assert(nBlksCMSA(CMA) == 1);
    ss_type             data=DataCMSA(CMA);
    st_type             sites = SitesCMSA(CMA);
    fm_type             *model=ModelsCMSA(cma);
    Int4                s,end;
    double              prob;
    unsigned char       *seq;

    if(n < 1 || n > NSeqsSeqSet(data)) print_error("GetSqProbCMSA( ) input error");
    end = SqLenSeqSet(n,data) + 1;
    end += 1 - LenFModel(model[1]);
    seq = SeqSeqSet(n,data);
    s = SitePos(1,n,1,sites);
    prob = (double) LikelihoodFModel(seq, s, model[1]);
    if(prob > 0.0) prob=log10(prob); else assert(prob > 0.0);
    return prob;
}

smx_typ	*omc_typ::ProfileHiHMM(Int4 nn)
// Create the optimal HMM for node nn.
{
	//=============== 1. Check input for validity. =================
	hpt_typ *Hpt=mcs->GetHpt();
	assert(nn > 1 && nn < Hpt->NumSets());
	char	dms_mode='F',wt_factor=100;	// 
	Int4	pernats=1443;
	//=============== 2. Obtain main datatypes. =================
	dms_typ	*dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	che_typ **che=mcs->RtnChe( );	
        // mcs->PutHyperPartition(stdout); 
	double	dd,d,Dd,FG,BG,DD,dD,RootFG,aveFG,**rtnDD=0;
	double	max=-INT4_MAX,min=INT4_MAX; 
	double	Max=-INT4_MAX,Min=INT4_MAX; 
	Int4	i,m,n,x,*Parent; Hpt->IsTree(Parent);
	Int4	X,parent=Parent[nn];
	dh_type dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	dh_type ave_dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	double *fg_freq,*bg_freq,*std_freq=tFreqCMSA(TrueMainCMA);
	UInt4  **WtCntsPrntFG=che[parent]->GetResWtsFG();
	UInt4  **WtCntsPrntBG=0;
	if(parent != 1) WtCntsPrntBG=che[parent]->GetResWtsBG();
	UInt4  **WtCntsFG=che[nn]->GetResWtsFG();
	UInt4  **WtCntsBG=che[nn]->GetResWtsBG();
	smx_typ	*smx; NEW(smx,5,smx_typ);
	smx[0]=MkSMatrix(2.0, LengthCMSA(1,TrueMainCMA),tFreqCMSA(TrueMainCMA),AB);
	smx[1]=MkSMatrix(2.0, LengthCMSA(1,TrueMainCMA),tFreqCMSA(TrueMainCMA),AB);
	smx[2]=MkSMatrix(2.0, LengthCMSA(1,TrueMainCMA),tFreqCMSA(TrueMainCMA),AB);
	
	double	diff,scoreFG[50],scoreBG[50];
	for(i=1; i <= LengthCMSA(1,TrueMainCMA); i++){
	    DD=dms->bild(WtCntsPrntFG[i]); dd=dms->bild(WtCntsFG[i]); Dd=dms->bild(WtCntsBG[i]); 
	    diff = DD-(dd + Dd);
	    if(diff < 0.0){		// subgroup is better...
	    	fprintf(stderr,"bild %d: %.2f + %.2f = %.2f > %.2f (%.2f)\n",i,Dd,dd,(Dd+dd),DD,diff); 
#if 0
		dms->score(WtCntsFG[i],scoreFG); dms->score(WtCntsBG[i],scoreBG);
		for(x=1; x <= nAlpha(AB); x++){
			X=(Int4) floor((scoreFG[x]-scoreBG[x]) + 0.5); SetSMatrix(x,i,X,smx[1]);
		}
#else		// use standard BG score.
		dms->score(WtCntsFG[i],scoreFG);
		for(x=1; x <= nAlpha(AB); x++) SetSMatrix(x,i,(Int4)floor(scoreFG[x] + 0.5),smx[1]);
#endif
	    } else {			// parent group is no worse.
	    	fprintf(stderr,"bild %d: %.2f > %.2f (%.2f)\n",i,DD,Dd+dd,diff); 
	        if(1 || parent == 1){
			dms->score(WtCntsPrntFG[i],scoreFG); 
			for(x=1; x <= nAlpha(AB); x++) SetSMatrix(x,i,(Int4)floor(scoreFG[x] + 0.5),smx[1]);
		} else {
			dms->score(WtCntsPrntFG[i],scoreFG); dms->score(WtCntsPrntBG[i],scoreBG); 
			for(x=1; x <= nAlpha(AB); x++){
				X = (Int4) floor((scoreFG[x]-scoreBG[x]) + 0.5);
				SetSMatrix(x,i,X,smx[1]);
			}
		}
	    }
	    dms->score(WtCntsFG[i],scoreFG); 	// FG only = smx[2].
	    for(x=1; x <= nAlpha(AB); x++) SetSMatrix(x,i,(Int4)floor(scoreFG[x] + 0.5),smx[2]);
	    dms->score(WtCntsPrntFG[i],scoreFG); 	// parent FG only = smx[0].
	    for(x=1; x <= nAlpha(AB); x++) SetSMatrix(x,i,(Int4)floor(scoreFG[x] + 0.5),smx[0]);
	    
	    // if(DD > max) max=DD; if(DD < min) min=DD; insrtHeap(i,-(keytyp)DD,dH); fprintf(stdout," ==> %.1f\n",DD); 

	    // aveFG = aveFG/(double)m; insrtHeap(i,-(keytyp)aveFG,ave_dH); if(aveFG > Max) Max=aveFG; if(aveFG < Min) Min=aveFG;
	} delete dms;
#if 0
	h_type	HG1=Histogram("BILD scores",-50,100,5);
	while(!emptyHeap(dH)){
		dd=-minkeyHeap(dH); assert((i=delminHeap(dH)) != 0);
		Dd = 100*(dd - min)/(max - min); IncdHist(Dd, HG1);
		fprintf(stdout,"column %d: %.1f average dF (%.1f)\n",i,dd,Dd);
	} fprintf(stdout,"\n"); Nildheap(dH);
	fprintf(stdout,"max=%.3f; min=%.3f\n",max,min); 
	PutHist(stdout,60,HG1); NilHist(HG1);

	HG1=Histogram("aveBILD scores",-50,100,5);
	while(!emptyHeap(ave_dH)){
		dd=-minkeyHeap(ave_dH); assert((i=delminHeap(ave_dH)) != 0);
		Dd = 100*(dd - Min)/(Max - Min); IncdHist(Dd, HG1);
		fprintf(stdout,"column %d: %.1f average dF (%.1f)\n",i,dd,Dd);
	} fprintf(stdout,"\n"); Nildheap(ave_dH);
	fprintf(stdout,"Max=%.3f; Min=%.3f\n",Max,Min); 
	PutHist(stdout,60,HG1); NilHist(HG1);
#endif
	return smx;
}

Int4	omc_typ::TestHiHMM( )
{
	Int4	id,i,j,n,N,s,iter=0,result=0,No,mode=0;
	double	best_lpr,last_lpr,d,D,dd,DD,EE,e;
	mcs_typ *xmcs=0,*rmcs=0;

	SetStringency(stringency);
	assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
	mcs->NoFailureMode=TRUE;  // if node configuration subLPR <= 0.0, then reject it.
	mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
	if(Evolve) mcs->DoEvolve(); else mcs->DoNotEvolve();
	Int4 c,nn,*Parent=0;	
	// for files in /home/aneuwald/working/PIPELINE/GISMO/AAAplus/FULL_SET/JUNK
	   // nn=14; // == moxR = Set41. // misaligned
	nn=57; // == dynein = Set12.
	nn=56; // == DnaA = Set10.
	nn=2; // == exeR = Set14.
	nn=50; // == AFG1 = Set40.
	nn=51; // == gamma + delta' clamp loaders = Set6.
	nn=47; // == SARP = Set33.
	nn=27; // == clpX = Set23.
	nn=34; // == lon = Set44.

	smx_typ	*smx=this->ProfileHiHMM(nn);
	PutSMatrix(stderr,smx[1]); PutSMatrix(stderr,smx[2]);
	h_type bgHG=Histogram("Background scores",-2000,2000,10);
	h_type stdHG=Histogram("Standard scores for node",-2000,2000,10);
	sprintf(str,"HiHMM scores for node %d (%s)",nn,Hpt->ElmntSetName(nn));
	h_type	HG=Histogram(str,-2000,+2000,10);
	for(Int4 sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
		for(DD=dd=EE=0.0,c=1; c <= LengthCMSA(1,TrueMainCMA); c++){
	   	    Int4 r=ResidueCMSA(1,sq,c,TrueMainCMA);
		    D=ValSMatrix(c,r,smx[2]); 	// standard HMM for subgroup 
		    d=ValSMatrix(c,r,smx[1]);	// hiHMM
		    e=ValSMatrix(c,r,smx[0]);
		    D=(Int4) floor((D/1000) + 0.5); 
		    d=(Int4) floor((d/1000) + 0.5);
		    e=(Int4) floor((e/1000) + 0.5);
		    DD += D; dd += d; EE+=e; 
	   	} IncdHist(dd,HG); IncdHist(DD,stdHG); IncdHist(EE,bgHG);
	}
#if 1
	Int4 pn=1443;  // This is what Sean Eddy is using: ln(x)*1.443 * 1000; 1/1443 nats == 1/1000 bits.
	cma_typ	cma=TrueMainCMA;
        // dms_mode='O';
        char dms_mode='T';       // nearly as good as 'F'.
        dms_mode='f';       // nearly as good as 'F'.
        dms_mode='F';       // best setting...
        // Int4   aa_per_io=200,aa_per_do=200,exp_ie=1,exp_de=1;
        Int4    aa_per_io=20,aa_per_do=150, exp_ie=1,exp_de=1;      // best setting...
#if 0
        ssx_typ *ssx = new ssx_typ(aa_per_io,aa_per_do,exp_ie, exp_de, pn,cma,dms_mode);
        // if(PriorWt > 0) ssx->SetPriorWt(PriorWt);
        // if(SqWtAdj > 0) ssx->SetSqWtAdjust(SqWtAdj);
#endif

	Int4 **mat_emit= ValuesSMatrix(smx[1]),**ins_emit=0;
	char *name= NameCMSA(cma);
        Int4 len=LengthCMSA(1,cma);
{
        Int4 *mm=0,*mi=0,*md=0,*ii=0,*im=0,*dd=0,*dm=0,*bm=0,*me=0;
#if 0
        ndl_typ *ndl=ssx->RtnNDL();
        ndl->GetTransProb(&mm,&mi,&md,&ii,&im,&dd,&dm,&bm);
	// m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
	//  -2 -10254 -11296   -894  -1115   -701  -1378    -42     0
#endif
        hmm_typ *hmm = new hmm_typ("hiHMM",len,mat_emit,ins_emit,mm,mi,md,ii,im,dd,dm,bm,me,AB,'R');
        hmm->Put(stdout);
        delete hmm;
	mat_emit= ValuesSMatrix(smx[2]);
        hmm = new hmm_typ("stdHMM",len,mat_emit,ins_emit,mm,mi,md,ii,im,dd,dm,bm,me,AB,'R');
        hmm->Put(stdout);
        delete hmm;
        fprintf(stderr,"Use hmmcalibrate to add EVD parameters. Use hmmsearch to search for hits.\n");
// for(Int4 x =1; smx[x] != 0; x++) NilSMatrix(smx[x]); free(smx);
}
#endif
	PutHist(stderr,60,bgHG); NilHist(bgHG);
	PutHist(stderr,60,HG); NilHist(HG);
	PutHist(stderr,60,stdHG); NilHist(stdHG);
	NilSMatrix(smx[0]); NilSMatrix(smx[1]); NilSMatrix(smx[2]); free(smx);
}

Int4	omc_typ::TestHiHMM2( )
// older code for delta-F analysis.
{
	Int4	id,i,j,n,N,s,iter=0,result=0,No,mode=0;
	double	best_lpr,last_lpr,d,D,dd,DD;
	mcs_typ *xmcs=0,*rmcs=0;

	SetStringency(stringency);
	assert(Hpt->NumBPPS() == Hpt->NumSets()-1);
	mcs->NoFailureMode=TRUE;  // if node configuration subLPR <= 0.0, then reject it.
	mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
	if(Evolve) mcs->DoEvolve(); else mcs->DoNotEvolve();
        // Boltzmann-like distributions.
	{
	mcs->PutHyperPartition(stdout);
	dh_type dH=dheap(LengthCMSA(1,TrueMainCMA)+2,4);
	Int4 c;
	NEWPP(BoltzScore, Hpt->NumBPPS() + 3, double);
	NEWPP(BoltzStd, Hpt->NumBPPS() + 3, double);
	for(n=1; n<= Hpt->NumBPPS(); n++){
	       NEWP(BoltzScore[n], LengthCMSA(1,TrueMainCMA) + 3, double);
	       NEWP(BoltzStd[n], LengthCMSA(1,TrueMainCMA) + 3, double);
               for(c=1; c <= LengthCMSA(1,TrueMainCMA); c++){
	          NEW(BoltzScore[n][c], nAlpha(AB)+3, double);
	          NEW(BoltzStd[n][c], nAlpha(AB)+3, double);
	       }
	}
	double *aveFC,**FC;  
	NEW(aveFC,Hpt->NumBPPS() + 3,double);
	NEWP(FC,Hpt->NumBPPS() + 3,double);
	NEWP(Entropy,Hpt->NumBPPS() + 3,double);
	for(n=1; n<= Hpt->NumBPPS(); n++){
	      NEW(FC[n],LengthCMSA(1,TrueMainCMA) + 3,double);
	      NEW(Entropy[n],LengthCMSA(1,TrueMainCMA) + 3,double);
	}
	Int4 m,*P,*NumData; Hpt->IsTree(P);
	set_typ Leaves=Hpt->MkLeafSet();
	NEW(NumData,LengthCMSA(1,TrueMainCMA) +4,Int4);
	h_type HG,aveHG=Histogram("Ave range of functional capacity",0,200,2);
        for(c=1 ; c <= LengthCMSA(1,TrueMainCMA); c++){
		// if(c==47) c=73; else if(c > 84) break; // c=35;
	   	HG=Histogram("Ranges of functional capacity",0,200,1);
	        for(m=0,D=0.0,n=2; n <= Hpt->NumBPPS(); n++){  // Skip root node!!
		   // if(P[n] != 1) continue;
		   // if(!MemberSet(n,Leaves)) continue;
	   	   d=mcs->PutBoltzmannLike(stdout,n,c);
		   if(!isfinite(d)) continue;
		   if(d <= 0) d=0.0;
		   else { IncdHist(d,HG); D += d; FC[n][c]=d; m++; }
		}
		// if(m > 0){ D = D/(double)(m); IncdHist(D,aveHG); }
		if(m > 0){ D = D/(double)(n-1); IncdHist(D,aveHG); }
		NumData[c]=NumDataHist(aveHG);
		insrtHeap(c,-(keytyp)D,dH);
		fprintf(stdout,"======= Column %d: %.1f ave =====\n",c,-D);
		PutHist(stdout,60,HG); NilHist(HG); 
	} free(P); NilSet(Leaves);
	PutHist(stdout,60,aveHG); NilHist(aveHG);
	while(!emptyHeap(dH)){
		D=minkeyHeap(dH);
	        assert((j=delminHeap(dH)) != 0);
		fprintf(stdout,"column %d: %.1f average dF (%d pts)\n",j,-D,NumData[j]);
	} fprintf(stdout,"\n"); Nildheap(dH);

	for(n=1; n<= Hpt->NumBPPS(); n++){
              for(D=0.0,c=1; c <= LengthCMSA(1,TrueMainCMA); c++) D+= FC[n][c];
	      aveFC[n]=D/(double)LengthCMSA(1,TrueMainCMA);
	}

	set_typ *set=mcs->CopyOfSeqSets();
	Int4 mm,nn,*Parent=0;	
	// for files in /home/aneuwald/working/PIPELINE/GISMO/AAAplus/FULL_SET/JUNK
	   // nn=14; // == moxR = Set41. // misaligned
	nn=57; // == dynein = Set12.
	nn=56; // == DnaA = Set10.
	nn=2; // == exeR = Set14.
	nn=50; // == AFG1 = Set40.
	nn=51; // == gamma + delta' clamp loaders = Set6.
	nn=47; // == SARP = Set33.
	nn=27; // == clpX = Set23.
	nn=34; // == lon = Set44.
	double FC_min=40; FC_min=0;

	// smx_typ	*smx=this->ProfileHiHMM(nn);
	// print out histograms for each set...
	h_type stdHG=0,prbHG=0;
	h_type scrHG=Histogram("BoltzScores",-200,200,2);
#if 0
	stdHG=Histogram("Standard FG functional capacity",-2000,5000,25);
	sprintf(str,"BoltzScores for node %d (%s)",nn,Hpt->ElmntSetName(nn));
	HG=Histogram(str,-2000,+5000,25);
	// if(Hpt->TypeOfSet(n) == '*') continue; // skip internal nodes...
	for(c=1 ; c <= LengthCMSA(1,TrueMainCMA); c++){ mcs->PutBoltzmannLike(stderr,nn,c); }
	for(Int4 sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
		for(DD=dd=0.0,c=1; c <= LengthCMSA(1,TrueMainCMA); c++){
	   	    Int4 r=ResidueCMSA(1,sq,c,TrueMainCMA);
		    d = BoltzScore[nn][c][r]; D = BoltzStd[nn][c][r];

#if 0
	            if(d != 0) fprintf(stderr,"BoltzScore[%d][%d][%c] = %.3f\n",nn,c,AlphaChar(r,AB),d); 
	            if(D != 0) fprintf(stderr,"BoltzStd[%d][%d][%c] = %.3f\n",nn,c,AlphaChar(r,AB),D); 
#endif
		    DD += D; dd += d; IncdHist(d,scrHG);
	   	} IncdHist(dd,HG); IncdHist(DD,stdHG);
	} 
	PutHist(stdout,60,HG); NilHist(HG);
	PutHist(stdout,60,stdHG); NilHist(stdHG);
	PutHist(stdout,60,scrHG); NilHist(scrHG);
	exit(1);
#endif
	for(mm=1; mm < Hpt->NumSets(); mm++){
	    mm=nn;
	    set_typ subtreeFG=MakeSet(MaxNumNodesPlus+2);
	    ReSetSubTree(subtreeFG,mm); fprintf(stderr,"%d: ",nn); PutSet(stderr,subtreeFG);
	    stdHG=Histogram("Standard FG functional capacity",-2000,5000,25);
	    sprintf(str,"BoltzScores for node %d (%s)",mm,Hpt->ElmntSetName(mm));
	    HG=Histogram(str,-2000,+5000,25);
	    for(n=1; n < Hpt->NumSets(); n++){
	     if(!MemberSet(n,subtreeFG)) continue;
	     // if(Hpt->TypeOfSet(n) == '*') continue; // skip internal nodes...
	     for(Int4 sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
		// if(!MemberSet(sq,set[n])) continue;
		for(DD=dd=0.0,c=1 ; c <= LengthCMSA(1,TrueMainCMA); c++){
		    // if(FC[n][c] < aveFC[n]) continue;
		    // if(FC[n][c] < FC_min) continue;
	   	    Int4 r=ResidueCMSA(1,sq,c,TrueMainCMA);
		    // d = BoltzScore[n][c][r]; D = BoltzStd[n][c][r];
		    d = BoltzScore[mm][c][r]; D = BoltzStd[mm][c][r];
#if 0
	        fprintf(stderr,"BoltzScore[%d][%d][%c] = %.3f\n",n,c,AlphaChar(r,AB),d); 
	        fprintf(stderr,"BoltzStd[%d][%d][%c] = %.3f\n",n,c,AlphaChar(r,AB),D); 
#endif
		    DD += D; dd += d; 
		    IncdHist(d,scrHG);
	   	} IncdHist(dd,HG); IncdHist(DD,stdHG);
	     }
	    } PutHist(stdout,60,HG); 
	    NilHist(HG); NilSet(subtreeFG);  subtreeFG=0;
	    PutHist(stdout,60,stdHG); NilHist(stdHG);
	break;
	} PutHist(stdout,60,scrHG); NilHist(scrHG);
exit(1);

	assert(Hpt->IsScrambledTree(Parent));

	cma_typ KeyCMA=GetInSetCMSA(set[nn],TrueMainCMA);
#if 0
	sprintf(str,"LLR FG scores for node %d (%s)",nn,Hpt->ElmntSetName(nn));
	prbHG=Histogram(str,-2000,5000,25);
	for(Int4 sq=1; sq <= NumSeqsCMSA(KeyCMA); sq++){
		d=20*GetTotalProbCMSA(sq,KeyCMA);
		IncdHist(d,prbHG); IncdHist(d,prbHG);
	} PutHist(stdout,60,prbHG); NilHist(prbHG);
#endif

	FILE *xfp=0; // xfp=open_file("junk",".cma","w"); 
	double BS_cut=10;
	set_typ	TopSet=0;
	   if(xfp){ TopSet=MakeSet(NumSeqsCMSA(TrueMainCMA)+9); ClearSet(TopSet); }
	   set_typ subtree=MakeSet(MaxNumNodesPlus+2);
	   ReSetSubTree(subtree,Parent[nn]); fprintf(stderr,"%d: ",Parent[nn]); PutSet(stderr,subtree);
	   sprintf(str,"BoltzScores for node %d (%s) BG seqs",nn,Hpt->ElmntSetName(nn));
	   HG=Histogram(str,-2000,+5000,25);
	   stdHG=Histogram("Standard functional capacity",-2000,5000,25);
	   sprintf(str,"LLR BG scores for node %d (%s)",nn,Hpt->ElmntSetName(nn));
	   prbHG=Histogram(str,-3500,5000,25);
#if 0
	   cma_typ RandCMA=ShuffleSeqCMSA(1.0, TrueMainCMA);
	   TotalNilCMSA(TrueMainCMA); TrueMainCMA=RandCMA;
#endif
	   for(n=1; n < Hpt->NumSets(); n++){
#if 0
	     if(MemberSet(n,subtree)) continue;
#else
	     if(!MemberSet(n,subtree)) continue;
	     if(n == nn){ continue; }
	     if(n == Parent[nn]) continue;	// likely to be some in here!
	     // if(Parent[n] != Parent[nn]) continue;	// Don't get too far removed.
#endif
	     for(Int4 sq=1; sq <= NumSeqsCMSA(TrueMainCMA); sq++){
		if(!MemberSet(sq,set[n])) continue;
		for(DD=dd=0.0,c=1; c <= LengthCMSA(1,TrueMainCMA); c++){
		        if(FC[n][c] < FC_min) continue;
	   		Int4 r=ResidueCMSA(1,sq,c,TrueMainCMA);
			// r = random_integer(nAlpha(AB)) + 1;
			d = BoltzScore[nn][c][r]; D = BoltzStd[nn][c][r];
			dd += d; DD += D;
	   	} IncdHist(dd,HG); IncdHist(DD,stdHG); 
		if(xfp && dd > BS_cut) AddSet(sq,TopSet);
		d=20*GetSqProbCMSA(KeyCMA,sq,TrueMainCMA);
		// d=GetTotalProbCMSA(sq,KeyCMA);
		IncdHist(d,prbHG); IncdHist(d,prbHG);
	     }
	   } PutHist(stdout,60,HG); PutHist(stdout,60,stdHG); 
	   NilHist(stdHG); NilHist(HG); 
	   PutHist(stdout,60,prbHG); NilHist(prbHG);
	   NilSet(subtree); TotalNilCMSA(KeyCMA);
	   if(xfp){ PutInSetCMSA(xfp,TopSet,TrueMainCMA);  NilSet(TopSet); fclose(xfp); }

	   free(Parent);
	   for(n=1; n<= Hpt->NumSets(); n++) if(set[n]) NilSet(set[n]); free(set);
	   for(n=1; n<= Hpt->NumBPPS(); n++){
               for(c=1; c <= LengthCMSA(1,TrueMainCMA); c++) free(BoltzScore[n][c]);
	       free(BoltzScore[n]);
	       free(FC[n]);
	   } free(BoltzScore); free(FC); free(aveFC);
#if 0
	   while(!emptyHeap(dH)){
		D=minkeyHeap(dH);
	        assert((j=delminHeap(dH)) != 0);
		fprintf(stdout,"column %d: %.1f average range\n",j,D);
	   } fprintf(stdout,"\n"); Nildheap(dH);
#endif
   } return 0;
}


