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

#include "scl_typ.h"

long double  LnFact(Int4 n)
/* static variables are guaranteed to be initialized to zero */
{
        static long double lnft[501];

        if (n <= 1) return 0.0;
        if (n <= 500) return lnft[n] ? lnft[n] : (lnft[n]=lgammal(n+1.0));
        else return lgammal(n+1.0);
}

long double  CumHyperGeomProb(Int4 N1,Int4 N2, Int4 n,Int4 x)
#if 0   //****************************************************
N total balls with N1 red balls and N2 = N-N1 black balls.
Choose n balls at random.  The probability that the
group so chosen will contain x or more red balls is
given by: p=CumHypGeomProb(N1,N2,n,x).
#endif  //****************************************************
{
        Int4    end;
        long double  p,K;

        if(x == 0) return 1.0;
        end = MINIMUM(Int4,N1,n);
        K = (LnFact(N1)+LnFact(N2)-LnFact(N2+N1)+LnFact(n)+LnFact(N2+N1-n));
        for(p=0.0,x; x<=end; x++){
           p += expl(K-LnFact(x)-LnFact(N1-x)-LnFact(n-x)-LnFact(N2-n+x));
        } return p;
}

long double	scm_typ::DCAvsPDB(FILE *fp,Int4 NumDistinguished,float **MtrxPDB, float **MtrxDCA,
		BooLean **Used)
// find the pvalue for DCA versus PDB.
{
	Int4	i,j,n,m,N,M=num_resC*num_resC,II,JJ;
	Int4	ri=0,ro=0,bi=0,bo=0,Nb,Nr,nO,nI,x;
	dh_type dH_pdb,dH_dca;

	// Insert pairs into dH_pdb heap using PDB distances as keys.
	// Insert pairs into dH_dca heap using DCA scores as keys.
	dH_pdb = dheap(M+5,4); dH_dca = dheap(M+5,4);
	for(n=0,i=1; i < num_resC; i++){
	   II=ResidueID(ResALL[i]);
	   for(j=i+1; j <= num_resC; j++){	
	      JJ=ResidueID(ResALL[j]);
	      if(abs(II-JJ) <= 5) continue;	// EV-fold requires |JJ-II| > 5.
	      // fprintf(stderr,"%d: %d vs %d\n",i,II,JJ);
	      if(MtrxPDB[i][j] > 0 && Used[i][j]){	// is this a legitimate pair?
	        n++;	// pair identifier == n;
		insrtHeap(n,(keytyp)MtrxPDB[i][j],dH_pdb); 
	        if(Simulate > 0) insrtHeap(n,(keytyp)Random(),dH_dca); 
	        else insrtHeap(n,(keytyp)MtrxDCA[i][j],dH_dca); 
	      }
	   }
	} N=ItemsInHeap(dH_pdb); assert(n == N); 
	if(NumDistinguished >= N/2){	// then stop here & return -1.0
	    fprintf(stderr,"NumDistinguished=%d; N=%d\n",NumDistinguished,N);
	    Nildheap(dH_pdb); Nildheap(dH_dca); return -1.0;	
	}

	// Pick D column pairs with the best DCA scores as the distinguished set.
	set_typ SetDCA=MakeSet(M);
	for(i=1; !emptyHeap(dH_dca); i++) {	
		double dD = (double) minkeyHeap(dH_dca);
		assert((n=delminHeap(dH_dca)) != 0); AddSet(n,SetDCA);
		if(i >= NumDistinguished) break;
	} Nr=CardSet(SetDCA); Nildheap(dH_dca); 
	
	// Sort all of the corresponding residue pairs based on their structural distances.
	Int4	*SortedByDistPDB; NEW(SortedByDistPDB,N+5,Int4);
	for(i=1; !emptyHeap(dH_pdb); i++) {
		double dP = (double) minkeyHeap(dH_pdb);
		assert((n=delminHeap(dH_pdb)) != 0); SortedByDistPDB[i]=n;
	} Nildheap(dH_pdb); 

	// Create an array of residue/column pairs 
	Int4 L,D,*pos,X;
	long double pval,prob=1.0;
	BooLean jf=FALSE,cf=TRUE,pf=TRUE;
	L = N;  D = Nr; NEW(pos,N+5,Int4);
	for(x=n=0,i=1; i <= N; i++){
	     n=SortedByDistPDB[i];
	     if(MemberSet(n,SetDCA)){
		x++; pos[x]=i; assert(x <= D); 
	        // fprintf(stderr,"pos[%d]=%d (%d)\n",x,pos[x],n);
	     }
	} X=pvcalcD(L,D,pos,pval,jf,cf,pf);   // Call ICA routine to compute the pval.

	if(fp){	// Compute probability using a Ball-In-Urn model.
	  set_typ	SetOut=MakeSet(M);
	  // Select drawn out pos[X] total balls & X red balls..
	  for(i=1; i <= pos[X]; i++){ n=SortedByDistPDB[i]; AddSet(n,SetOut); }
	  ro=CardInterSet(SetOut,SetDCA); ri=Nr-ro;	// compute # red balls in & out of urn.
	  Nb = N -Nr;		// compute # of black balls in & out of urn.
	  nO=CardSet(SetOut); 	// Total # of balls drawn out.
	  bo=nO-ro; bi=Nb-bo;	// # black balls in and out.
	  if(fp) fprintf(fp,"Nr=%d; Nb=%d; ro=%d; ri=%d; bo=%d; bi=%d; Out=%d; In=%d\n",
			Nr,Nb,ro,ri,bo,bi,ro+bo,ri+bi);
	  prob = CumHyperGeomProb(ro+ri,bo+bi,ro+bo,ro);
	  if(fp) fprintf(fp,"CumHypGeomProb(Nr=%d,Nb=%d,No=%d,ro=%d)=%.3Lg\n",ro+ri,bo+bi,ro+bo,ro,prob);
	  NilSet(SetOut); 
	}
	// Print the result if fp != NULL.
	if(fp) fprintf(fp,"pvcalc:  L=%d, D=%d; X=%d; #D before X=%d (%.1f%c); d/X=%.3f.\n",
			L,D,pos[X],X,100*(double)X/(double)D,'%',(double)X/(double)pos[X]);
	if(fp) fprintf(fp,"pvcalc p-value = %.3Lg (%.2Lf).\n",pval,-log10l(pval));
	free(SortedByDistPDB); free(pos); NilSet(SetDCA); 
	return pval; 
}

Int4	scl_typ::RunDCA( )
{
	Int4	i,j,k,n,x;
	Int4	resStart=0,resEnd=0,**PttrnSites,nAln,II=0; 
	char	c=' ',*chain=0,**pdb_file=0,*res_class;
	FILE	*infp=open_file(infile,"","r");
	pdb_typ	PDB;
	set_typ RtnSet=0;
	do {
	  n=ReadPttrnRes(infp,resStart,resEnd,RtnSet,nAln,pdb_file,chain,res_class,PttrnSites);
#if 0	// Unused data structures by scm_typ.
	  nAln=0;
	  free(res_class); res_class=0;
#endif
#if 0
	  for(i=1; i <= n; i++){ 
		for(j=1; PttrnSites[i][j]; j++) fprintf(stderr,"%d.%d: %d\n",i,j,PttrnSites[i][j]); 
	  } fprintf(stderr,"n=%d\n",nAln); 
#endif
	  if(n > 0){
	    for(j=1; pdb_file[j] != 0; j++){
	      assert(chain[j] != 0);
	      PDB = MakePDB(pdb_file[j]);	// PutPDB(stdout,PDB);
	      II++;
	      if(ChainI){
		for(k=0; ChainI[k]; k++){
	          if(!ChainExistsPDB(ChainI[k],PDB)){
		    free(pdb_file[j]); NilPDB(PDB); PDB=0; break;
		  }
		} if(PDB==0) continue;
		// Don't include query chain!
	      	if(strchr(ChainI,chain[j]) != NULL){ free(pdb_file[j]); NilPDB(PDB); continue; }
	      }
	      fprintf(stdout,"\n=========== %d. %s:%c ============\n",II,pdb_file[j],chain[j]);
	      if(dcaE){	// check for consistency...
		e_type E=GetPDBSeq(GetChainNumberPDB(PDB,chain[j]),PDB);
		a_type ab=AminoAcidAlphabetPDB(PDB); AlnSeqSW(stderr,11, 1,dcaE, E, ab);
        	if(!this->OverlappingSeqs(dcaOS,30,0,dcaE,E)){
           	  fprintf(stderr,"Sequences mismatch (%s)\n",pdb_file[j]);
           	  // PutSeq(stderr,dcaE,ab); PutSeq(stderr,E,ab);
		  free(pdb_file[j]); NilPDB(PDB); continue; 
        	} 
		Int4 pdbOS=OffSetSeq(E),os=OffSetSeq(dcaE); dcaOS=os-pdbOS+dcaOS;
		fprintf(stderr,"pdbOS= %d; os=%d; dcaOS=%d\n",pdbOS,os,dcaOS);
		NilSeq(E);
	      }

	      // Compute the significance of DCA-structural overlap.
	      scm_typ scm(resStart,resEnd,RtnSet,nAln,chain[j],res_class,PttrnSites,PDB,this);
	      assert(dcaE);
	      long double	pMin=1.0,pmin=1.0,score,pp;
	      Int4	xx,Nr,bstNr=0,BstNr=0,start=LenSeq(dcaE),end=15*start;
	      Int4	inc=10,Inc=25,Sim=0;
	      float    **MtrxPDB=0,**MtrxDCA=0;
	      BooLean	**Used=0;
	      MtrxPDB=scm.CalcDistMtrx(); MtrxDCA=scm.GetMtrxDCA(dca_file,Used); 
	      Int4 num_resC=scm.RtnNumResC();
#if 0	// compute max # of possible pairwise contacts within 10 angstroms...
	      start=MINIMUM(Int4,LenSeq(dcaE),50);
	      for(end=0,i=1; i < num_resC; i++){
		for(j=i+1; j <= num_resC; j++){	
		   if(MtrxPDB[i][j] <= 10.0 && Used[i][j]) end++;
	   	}
	      } end = MAXIMUM(Int4,end,(start*2));
#endif
	      h_type	HG=0;
	      if(Simulate > 0){ HG=Histogram("-log(p)",0,20,0.25); Sim=Simulate; }
	      else HG=Histogram("-log(p)",0,500,5.0);
	      do {
		 for(xx=0,Nr=start; Nr <= end; Nr+= inc){ 
		   score=1.0; pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp < 0) break;
		   if(pp > 0){ score=-log10l(pp); IncdHist(score,HG); }
		   if(pp < pmin){
			fprintf(stderr,"%d. %.1Lf (%d)\n",Nr,score,xx);
			pmin=pp; bstNr=Nr; xx=0; 
		   } else xx++;
		 } BstNr=bstNr;
		 for(Nr=BstNr+Inc; Nr > BstNr; Nr--){ 
		   pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp > 0){ score=-log10l(pp); IncdHist(score,HG); }
		   if(pp < pmin){ pmin=pp; bstNr=Nr; }
		 }
		 for(Nr=BstNr-Inc; Nr < BstNr; Nr++){ 
		   pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp > 0){ score=-log10l(pp); IncdHist(score,HG); }
		   if(pp < pmin){ pmin=pp; bstNr=Nr; }
		 }
		 if(Sim > 0){ 
		   fprintf(stderr,"======== %d ==========\n",Sim); Sim--; 
		   fprintf(stderr,"pmin=%.3Lg; bstNr=%d\n",pmin,bstNr);
		 }
	      } while(Sim > 0);
	      PutHist(stderr,60,HG); NilHist(HG);
	      if(1 || Simulate == 0){
	           fprintf(stdout,"\n=========== %d. %s:%c ============\n",II,pdb_file[j],chain[j]);
	           fprintf(stdout,"start=%d; end=%d\n",start,end);
		   scm.DCAvsPDB(stdout,bstNr,MtrxPDB,MtrxDCA,Used); 
	      } else fprintf(stderr,"pmin=%.3Lg; bstNr=%d\n",pmin,bstNr);
	      for(i=1; i <= num_resC; i++){ free(MtrxPDB[i]); } free(MtrxPDB);
	      for(i=1; i <= num_resC; i++){ free(MtrxDCA[i]); } free(MtrxDCA);
	      for(i=1; i <= num_resC; i++){ free(Used[i]); } free(Used);
	      free(pdb_file[j]); NilPDB(PDB);
	    }
	    for(i=1; i <= n; i++){ if(PttrnSites[i]) free(PttrnSites[i]); }
	    free(chain); free(pdb_file); free(PttrnSites); if(res_class) free(res_class); NilSet(RtnSet);
	  }
	} while(n); fclose(infp);
	return 0;
}

#define DUSAGE "Usage: evalDCA klst_file dca_score_file [options]\n\
       Computes p-value for overlap between DCA scores and structural contacts.\n\
       -J          Use (i.e., don't correct for) Jeffreys priors.\n\
       -S          Randomize input to perform <int> simulations.\n\
       -seed=<int> Provide a random seed.\n\
\n\n"

void	scl_typ::InitDCA(int argc,char *argv[])
{
	Int4	arg;
	char	str[150],*map_file=0;
	BooLean	openVSI=FALSE,OpenTable=FALSE;

	InitDefaults( );
	if(argc < 3) print_error(DUSAGE);
	program=AllocString(argv[0]);
	infile=AllocString(argv[1]);
	dca_file=AllocString(argv[2]);
	AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
        for(arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(DUSAGE);
           switch(argv[arg][1]) {
                case 'J': { if(argv[arg][2] != 0) print_error(DUSAGE); else UseJeffreys=TRUE; } break;
                case 'S': if(argv[arg][2] != 0) print_error(DUSAGE); else Simulate=1; break; 
		case 's': if(sscanf(argv[arg],"-seed=%u",&seed) != 1) print_error(DUSAGE);; break;
                default: print_error(DUSAGE);;
            }
        }
        FILE *fp=0,*ifp=0;
	char	Str[105],cI,cJ;
	double Dd;
	// if using GaussDCA output as input...then need cma file to compute true positions.
	if((fp=fopen(dca_file,"r")) == NULL){	// then look for dca_file.cma & dca_file.dca as input.
	     fp=open_file(dca_file,".cma","r");
	     cma_typ cma=ReadCMSA(fp,AB); fclose(fp);
	     assert(nBlksCMSA(cma) == 1);
	     Int4 ii,jj,i,j,len=LengthCMSA(1,cma);
	     char *Seq; NEW(Seq, len+5, char); 
	     e_type qE=TrueSeqCMSA(1,cma);
	     ifp=open_file(dca_file,".dca","r");
	     fp=open_file(dca_file,"","w");	// dca_file as output...adjusts pos and adds residues.
	     Int4 pc[12],os=OffSetSeq(qE);
	     double *X; NEW(X,LenSeq(qE)+os+5,double);
	     set_typ SetAln=0;
#if 1	// ignore inserts
	     SetAln=MakeSet(LenSeq(qE)+len+9); ClearSet(SetAln);
	     for(i=1; i <= len; i++){
		j=TruePosCMSA(1,i,cma); //  + OffSetSeq(qE);
	        AddSet(j,SetAln); 
	     } PutSet(stderr,SetAln);
#endif
             while(fgets(Str,100,ifp) != NULL){
	      // PSICOV input format.
              if(sscanf(Str,"%d %d %d %d %lf",&i,&j,&pc[1],&pc[2],&Dd) == 5){	// PSICOV input
		assert(i > 0 && i <= LenSeq(qE)); assert(j > 0 && j <= LenSeq(qE));
		ii=i; jj=j; // ii=ii+os; jj=jj+os; 
if(SetAln && (!MemberSet(ii,SetAln) || !MemberSet(jj,SetAln))) continue;
		// ii=TruePosCMSA(1,i,cma)+os; jj=TruePosCMSA(1,j,cma)+os;
		cI=AlphaChar(ResSeq(ii,qE),AB); cJ=AlphaChar(ResSeq(jj,qE),AB);
                fprintf(fp,"%d,%d,%.5lf,0,0,0,0,0,0,0,%c,%c\n",ii,jj,Dd,cI,cJ);
		if(Dd > 0) { X[ii] += Dd; X[jj] += Dd; }
	      // GaussDCA input format 
              } else if(sscanf(Str,"%d %d %lf",&i,&j,&Dd) == 3){		// gDCA input
		assert(i > 0 && i <= len); assert(j > 0 && j <= len);
		// ii=TruePosCMSA(1,i,cma)+os; jj=TruePosCMSA(1,j,cma)+os;
		ii=TruePosCMSA(1,i,cma); jj=TruePosCMSA(1,j,cma);
		assert(ii > 0 && ii <= LenSeq(qE)); assert(jj > 0 && jj <= LenSeq(qE));
if(SetAln && (!MemberSet(ii,SetAln) || !MemberSet(jj,SetAln))) continue;
		cI=AlphaChar(ResidueCMSA(1,1,i,cma),AB); cJ=AlphaChar(ResidueCMSA(1,1,j,cma),AB);
                fprintf(fp,"%d,%d,%.5lf,0,0,0,0,0,0,0,%c,%c\n",ii,jj,Dd,cI,cJ);
		if(Dd > 0) { X[ii] += Dd; X[jj] += Dd; }
	      } else print_error("-EV option input error 0.");
	     } fclose(ifp); fclose(fp); 
if(SetAln) NilSet(SetAln);
	     fp=open_file(dca_file,".clst","w");	// Order DCA residues...
	     dh_type dH=dheap(LenSeq(qE)+5,4); 
	     for(i=1; i <= LenSeq(qE); i++) if(X[i] > 0.0) insrtHeap(i,-(keytyp)X[i],dH);;
	     for(j=0; !emptyHeap(dH); j++){
		i=delminHeap(dH); 
		// fprintf(stdout, "%d = %.5lf\n",i,X[i]); 
		if(j > 0) fprintf(fp, ",%d",i); else fprintf(fp, "C=%d",i); 
		if(j > 50) break;
	     } free(X); Nildheap(dH); fprintf(fp,"\n"); fclose(fp);
	     TotalNilCMSA(cma); fp=open_file(dca_file,"","r");
	}
	// Read in EVcouplings-formatted input file...
        Int4  max_ij=0,min_ij=9999;
        unsigned char *seq;     NEW(seq,1005,unsigned char); // 30-500 is evfold max...
        while(fgets(Str,100,fp) != NULL){
                Int4	i,j,dm[12];    // dummy variables.
                if(sscanf(Str,"%d,%d,%lf,%d,%d,%d,%d,%d,%d,%d,%c,%c,",
                   &i,&j,&Dd,&dm[1],&dm[2],&dm[3],&dm[4],&dm[5],&dm[6],&dm[7],&cI,&cJ) != 12)
                        print_error("GetMtrxDCA() input error 1.");
                assert(i < 1000 && j < 1000 && i > 0 && j > 0);
                if(i > max_ij) max_ij=i; if(j > max_ij) max_ij=j;
                if(i < min_ij) min_ij=i; if(j < min_ij) min_ij=j;
                seq[i]=AlphaCode(cI,AB); seq[j]=AlphaCode(cJ,AB);
                // EqSeq(i,AlphaCode(cI,AB),dcaE); EqSeq(j,AlphaCode(cJ,AB),dcaE);
        } fclose(fp);
        dcaE=MkSeq("DCA seq",max_ij-min_ij +1,seq + min_ij -1); free(seq);
        SetOffSetSeq(min_ij -1,dcaE);
        PutSeq(stderr,dcaE,AB); 
	time1=time(NULL);
        if(seed == 18364592) seed = (unsigned int) (time(NULL)/2);
	sRandom(seed);
        for(arg=0; arg < argc; arg++) fprintf(stdout,"%s ",argv[arg]); 
	if(Simulate) fprintf(stderr,"seed = %u\n",seed); else fprintf(stdout,"\n");
}

