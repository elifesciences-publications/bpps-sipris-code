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

#include "ssx_typ.h"
#include "blosum62.h"

ssx_typ::ssx_typ(Int4 aa_per_io,Int4 aa_per_do,Int4 exp_ie, Int4 exp_de, double pn,cma_typ c,
				char dm, ssx_typ *twin_ssx)
{
	twin=twin_ssx;
	dms_mode=dm;
	LGM=0;
	// JLH= new jlh_typ(aa_per_io,aa_per_do,exp_ie,exp_de,c,pn,wt_factor);
	// Init(pn,c); JLH->InitIndels(Wt);
	Init(pn,c); 
	assert(LGM != 0);
	NDL= new ndl_typ(aa_per_io,aa_per_do,exp_ie,exp_de,c,LGM,pn,wt_factor);
	NDL->InitIndels(Wt);
#if 1	// Fix for uninitialized gap penalties problem... afn: 11-16-2014.
	this->InitNDL(0.0);
#endif
}

void	ssx_typ::Init(double pn,cma_typ in_cma){ init(pn,in_cma); SeqWts(); InitCnts(); }

#if 1
static double Blsm62frq[21] = {  0.0, /*  X
    C      G        A       S       T       N       D       E       Q      K */
0.02469,0.07415,0.07422,0.05723,0.05089,0.04465,0.05363,0.05431,0.03426,0.05816,
/*  R       H       W       Y       F       V       I       L       M      P */
0.05161,0.02621,0.01303,0.03228,0.04742,0.07292,0.06792,0.09891,0.02499,0.03854};
#endif

void	ssx_typ::init(double pn,cma_typ in_cma)
{
	Int4	r,len,i,blk;
	CMA=in_cma;
	gss_typ *gss=gssCMSA(CMA); gss->SetLeftFlank(500); gss->SetRightFlank(500);
	AB=AlphabetCMSA(CMA); pernats=pn;
	N=NumSeqsCMSA(CMA); 
	SqWtAdjst=1.0;
	TotalWtPseudo=WtPseudo*(UInt4) nAlpha(AB);
	ReCalcSMX=TRUE;
assert(pernats==1000 || pernats == 693 || pernats == 1443);
	WtN=0; VrtlN=0; Matrix=0;
	NEW(Wt,N+3,char); 
	NEW(TtlWtCnts,nAlpha(AB)+3,UInt4); NEW(TtlWtFreq,nAlpha(AB)+3,double);
	NEW(BackGrnd,nAlpha(AB)+5,double);
	if(dms_mode != 'P'){
	  DMS = new dms_typ(wt_factor,pernats,dms_mode,AB);	// sets up dmscore();
	  double *BG=DMS->BackGrnd();
          for(r=0; r <= nAlpha(AB); r++){ BackGrnd[r]=BG[r]; }
	} else {
	  DMS=0;
#if 1
	  // DMS = new dms_typ(wt_factor,pernats,'T',AB);	// sets up dmscore();
	  Int4	s,**R; NEWP(R,nAlpha(AB)+3,Int4);
          for(r=0; r <= nAlpha(AB); r++){
	   NEW(R[r],nAlpha(AB)+3,Int4);
	   BackGrnd[r]=Blsm62frq[r];
           for(s=0; s <= nAlpha(AB); s++){
		double d=blosum62[r][s]; 
		R[r][s]=(Int4) (10*floor(d+0.5));
	   }
	  }
	  Matrix=mcalc(R);       // Stephen Altschul's routine...
          for(r=0; r <= nAlpha(AB); r++) free(R[r]); free(R);
	
#elif 1
          // for(r=0; r <= nAlpha(AB); r++) BackGrnd[r]=blosum62freq[r];
          for(r=0; r <= nAlpha(AB); r++) BackGrnd[r]=Blsm62frq[r];
	  Matrix=mcalc(AlphaR(AB));       // Stephen Altschul's routine...
#else
	  double	*frq=tFreqSeqSet(TrueDataCMSA(CMA));
          for(r=0; r <= nAlpha(AB); r++) BackGrnd[r]=frq[r];
	  Matrix=mcalc(AlphaR(AB));       // Stephen Altschul's routine...
#endif
	}
	NEWPP(WtCnts,nBlksCMSA(CMA)+3,UInt4);
	NEWPP(VrtlCnts,nBlksCMSA(CMA)+3,UInt4);
	NEWP(TtlVrtlCnts,nBlksCMSA(CMA)+3,UInt4);
	NEW(SMX, nBlksCMSA(CMA) +3, smx_typ);
        for(totlen=0,blk=1; blk<=nBlksCMSA(CMA); blk++){
		len=LengthCMSA(blk,CMA); totlen+=len;
	    	SMX[blk]=MkSMatrix(2.0,len,BackGrnd,AB);
	        NEWP(VrtlCnts[blk],len+3,UInt4);
		NEWP(WtCnts[blk],len +3,UInt4);
		NEW(TtlVrtlCnts[blk],len+3,UInt4);
                for(i=1; i<=len; i++){
		   NEW(WtCnts[blk][i],nAlpha(AB) +3,UInt4); 
		   NEW(VrtlCnts[blk][i],nAlpha(AB) +3,UInt4); 
		}
	}
}

void	ssx_typ::Free()
{
	Int4	i,blk;
	if(Matrix){ for(i=0; i<= nAlpha(AB); i++) free(Matrix[i]); free(Matrix); }
	free(BackGrnd); 
        for(blk=1; blk<=nBlksCMSA(CMA); blk++){
                for(i=1; i<=LengthCMSA(blk,CMA); i++){
			 free(WtCnts[blk][i]); free(VrtlCnts[blk][i]);
		} free(WtCnts[blk]); free(VrtlCnts[blk]); free(TtlVrtlCnts[blk]);
		NilSMatrix(SMX[blk]);
	} free(Wt); free(WtCnts); free(VrtlCnts); free(TtlVrtlCnts); free(SMX);
	// delete JLH;
	delete NDL;
	if(LGM) delete LGM;
	free(TtlWtFreq); free(TtlWtCnts);
	if(DMS) delete DMS;
}

void    ssx_typ::InitCnts()
// compute total counts for smatrix.
{
	Int4	i,j,len,st,pos[5],blk,sq,N=NumSeqsCMSA(CMA);
	unsigned char   r,*isq;

	for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	    UInt4	**wtCnts=WtCnts[blk];
	    len=LengthCMSA(blk,CMA);
	    for(i=1; i <= len; i++){ for(r=0; r <= nAlpha(AB); r++) wtCnts[i][r] = 0; }
	    for(sq=1; sq <= N; sq++){
	      isq=SeqPtrCMSA(sq,CMA);
	      if(PosSiteCMSA(blk,sq,pos,CMA)){
		 st=pos[1];
	         for(i=1; i <= len; st++,i++){ r=isq[st]; wtCnts[i][r] += Wt[sq]; }
	      }
	    }
	}
}

void	ssx_typ::ComputeTotalWtCnts()
// update only when updating SeqWts.
{
	UInt4	Total,*num; NEW(num,nAlpha(AB)+3,UInt4);
	Int4	sq,r;
	for(r=0; r <= nAlpha(AB); r++) TtlWtCnts[r] = 0;
        for(sq=1; sq<=N; sq++){
	   NumResSeq(TrueSeqCMSA(sq,CMA),num,AB);
	   for(r=0; r <= nAlpha(AB); r++) TtlWtCnts[r] += Wt[sq]*num[r];
	}
	for(Total=0,r=0; r <= nAlpha(AB); r++) Total+=TtlWtCnts[r];
	for(r=0; r <= nAlpha(AB); r++) TtlWtFreq[r] = (double)TtlWtCnts[r]/(double) Total;
	free(num);
}

void	ssx_typ::SeqWts()
// return an integer weight for each sequence in the cmsa.
// NOTE!!! Eventually call ComputSeqWtsCMSA() in cma_gmb.cc to replace below...
{
        Int4    i,j,k,r,sq,b,n;
        double  w,max,N2,d,*wt;
        unsigned char    *seq;
	UInt4	**nres,*ntyp;

#if 1	// for calculating comparable LLRs with variable # columns...
	if(twin){
		assert(twin->N == this->N);
		this->WtN=twin->WtN;
		char *tWt=twin->Wt;
		for(sq=1; sq <= N; sq++){ this->Wt[sq] = tWt[sq]; }
		ComputeTotalWtCnts();
		this->calc_lngamma();
		return;
	}
#endif
	NEW(wt,N+2,double); 
        NEWP(nres,totlen+2,UInt4); NEW(ntyp,totlen+2,UInt4);
        for(i=0; i<=totlen; i++) { NEW(nres[i],nAlpha(AB)+2,UInt4); }
        /*** 1. determine the number of residues at each position ***/
        for(sq=1; sq<=N; sq++){
           for(j=1,b=1; b <= nBlksCMSA(CMA); b++){
		seq=GetAlnResInSiteCMSA(b,sq,CMA);
		if(seq == 0) continue;
                for(i=1; i<=LengthCMSA(b,CMA); i++){
                    r=seq[i]; 
		    if(r==0) continue;
		    if(nres[j][r] == 0) { ntyp[j]++; } 
		    nres[j][r]++; j++;
                }
           }
        }

        /*** 2. determine the sequence weights. ***/
        for(max=0.,sq=1; sq<=N; sq++){
           for(w=0.0,j=1,b=1; b <= nBlksCMSA(CMA); b++){
		seq=GetAlnResInSiteCMSA(b,sq,CMA);
		if(seq == 0) continue;
                for(i=1; i<=LengthCMSA(b,CMA); i++){
                    r=seq[i]; 
		    if(r==0) continue;
                    w+= 1.0/ (double) (ntyp[j]*nres[j][r]); j++; 
                }
           } w /= (double) totlen; 
	   wt[sq]=w; max = MAXIMUM(double,w,max); 
        }
	h_type	HG=0;
#if 0	// For debugging code...
	HG = Histogram("number residue types", 0, 30,1.0);
        for(j=1; j<=totlen; j++){ IncdHist(ntyp[j],HG); }
	PutHist(stderr,60,HG);  NilHist(HG);

        /*** 3. normalize the weights. ***/
	HG = Histogram("sequence weights", 0, 100,2); 
	// HG = Histogram("sequence weights", 0, 100,1.0); 
#endif
        for(N2=0.0,sq=1; sq<=N; sq++){
	   wt[sq] /= max; 
	   if(SqWtAdjst < 1.0) wt[sq] = pow(wt[sq],SqWtAdjst);
	   N2 += wt[sq];
           // if(wt[s] > 0.9) fprintf(stderr,">%s (%g)\n",infoSMA(s,MA),wt[s]);
	   // fprintf(stderr,"seq %d: weight = %g\n",s,wt[s]); 
	   char tmp=(char) ceil((double)wt_factor*wt[sq]);
#if 0	// allow tmp == 0;	--> implies sequences excluded from alignment.
	   if(!(tmp <= wt_factor && tmp > 0)){
		// fprintf(stderr,"tmp=%d; wt_factor=%d\n",tmp,wt_factor);
		if(tmp==0) tmp=1;
	   	assert(tmp <= wt_factor && tmp > 0);
	   } else if(HG) IncdHist(tmp,HG);
#else
	   assert(tmp <= wt_factor && tmp >= 0);
	   if(HG) IncdHist(tmp,HG);
#endif
#if 1	// Fix for AAA+ gismo core dump with XXXXXXXXXXXXX pdb seqs.
	   if(tmp==0) tmp=1;
#endif
	   Wt[sq]= tmp;
        }
	for(WtN=0,sq=1; sq <= N; sq++){ WtN += (UInt4) Wt[sq]; }
	if(HG){ PutHist(stderr,60,HG);  NilHist(HG); }
	// fprintf(stderr,"%d seqs: effective number of sequences N = %f\n",N,N2);
        for(i=0; i<=totlen; i++) { free(nres[i]); }
        free(nres); free(ntyp); free(wt); 
	ComputeTotalWtCnts();
	this->calc_lngamma();
}

