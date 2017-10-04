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

double  ssx_typ::ContextBildScore(Int4 blk, Int4 i, double *BldSc,double blk_cut)
{
	Int4	j,k,n,x,ins,ins_res,nins,Ncols=LengthCMSA(blk,CMA);
	// double  SumWt=0,wt[9]={ 1.0, 0.750, 0.50, 0.25, 0.125 };
	double  SumWt=0,wt[9]={ 1.0, 0.67, 0.33, 0.0, 0.0 };
        double  d,BS,bs,bs_rt=0,bs_lt=0,wt_lt=0,wt_rt=0;
	Int4	flank=2;
        for(SumWt=0,n=0,bs=0.0,k=-flank; k <= flank; k++){
           j=i+k; x=abs(k);
           if(j > 0 && j <= Ncols){
               if(k < 0) nins=NumInsertsCMSA(1,j,CMA,ins_res);
               else if(k > 0) nins=NumInsertsCMSA(1,j-1,CMA,ins_res);
               else nins=0;
               double dd=(double)nins/(double)NumSeqsCMSA(CMA);
               if(dd > blk_cut){
                    if(k > 0) break; // stop here...end of block
                    else if(k < 0){ bs_lt=0; wt_lt=0; } // start over...
                    continue;
               }
#if 1
	       if(BldSc){	// 3 x's faster to compute...
                 if(k < 0){ bs_lt += BldSc[j]*wt[x]; wt_lt +=wt[x]; }
                 else if(k > 0){ bs_rt += BldSc[j]*wt[x]; wt_rt +=wt[x]; }
                 else { BS = BldSc[j]*wt[x]; bs+=BS; SumWt+=wt[x]; }
	       } else {
                 if(k < 0){ bs_lt += DMS->bild(WtCnts[blk][j])*wt[x]; wt_lt +=wt[x]; }
                 else if(k > 0){ bs_rt += DMS->bild(WtCnts[blk][j])*wt[x]; wt_rt +=wt[x]; }
                 else { BS = DMS->bild(WtCnts[blk][j])*wt[x]; bs+=BS; SumWt+=wt[x]; }
	       }
#else	// relative entropy...
               if(k < 0){ bs_lt += this->RelEntropy(1,j)*wt[x]; wt_lt +=wt[x]; }
               else if(k > 0){ bs_rt += this->RelEntropy(1,j)*wt[x]; wt_rt +=wt[x]; }
               else { BS = this->RelEntropy(1,j)*wt[x]; bs+=BS; SumWt+=wt[x]; }
#endif
           }
        }
#if 0	// pick one side or the other.
        if(wt_lt > 0 && bs_lt < bs_rt){ bs += 2*bs_lt*wt_lt; SumWt+=2*wt_lt; }
        else if(wt_rt > 0) { bs += 2*bs_rt*wt_rt; SumWt+=2*wt_rt; }
#elif 0	// Take the maximum of the three bild scores.
        if(wt_rt > 0 && bs_rt > bs){ bs = bs_rt; SumWt=wt_rt; }
        if(wt_lt > 0 && bs_lt > bs){ bs = bs_lt; SumWt=wt_lt; }
#else	// take the weighted sum of three scores.
        bs += bs_rt*wt_rt + bs_lt*wt_lt; SumWt += wt_rt + wt_lt;
#endif
	nins=NumInsertsCMSA(1,i,CMA,ins_res);
	d=bs/SumWt;
#if 0
        fprintf(stderr,"%d: BILD = %.3f; bs_lt = %.3f; bs_rt = %.3f; cBILD = %.3f; SumWt=%.3f\n",
             i,BS,bs_lt/wt_lt,bs_rt/wt_rt,d,SumWt);
        if(nins > 0) fprintf(stderr,"    %d insertions\n",nins);
#endif
	return d; 
}

void    ssx_typ::AddToAlign(Int4 sq,Int4 *site)
{
	Int4    r,st,i,len,blk,pos[5];

	assert(!SeqIsInCMSA(sq,CMA));
        for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	      AddSiteCMSA(blk,sq,site[blk], CMA);
	      // add sq to WtCnts...
              UInt4     **Cnts=WtCnts[blk];
              len=LengthCMSA(blk,CMA);
              unsigned char   *isq=SeqPtrCMSA(sq,CMA);
              assert(PosSiteCMSA(blk,sq,pos,CMA) > 0); st=pos[1];
	      assert(site[blk]==st);
              for(i=1; i <= len; st++,i++){
		r=isq[st];
#if 0
		if(r) Cnts[i][r] += (UInt4) Wt[sq]; 
#else
		Cnts[i][r] += (UInt4) Wt[sq]; 
#endif
	      }
        } ReCalcSMX=TRUE;
	// JLH->AddIndels(sq,Wt[sq]);
	NDL->AddIndels(sq,Wt[sq]);
}

Int4	*ssx_typ::RemoveFromAlign(Int4 sq)
// Remove a sequence from the alignment and return the old positions.
// Return null if sequence not in alignment!
{
	Int4	st,r,j,len,blk,pos[5],*site;
	
#if 0	// check for vacant site...afn: 11-26-2014.
	if(!SeqIsInCMSA(sq,CMA)) return 0;
#else 
	assert(SeqIsInCMSA(sq,CMA));
#endif
	NEW(site,nBlksCMSA(CMA) + 5, Int4);
	// JLH->RmIndels(sq,Wt[sq]);
	NDL->RmIndels(sq,Wt[sq]);
	for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	      UInt4	**Cnts=WtCnts[blk];
	      len=LengthCMSA(blk,CMA);
	      unsigned char   *isq=SeqPtrCMSA(sq,CMA);
	      assert(PosSiteCMSA(blk,sq,pos,CMA) > 0); st=pos[1]; site[blk]=st;
	      for(j=1; j <= len; st++,j++){
		r=isq[st]; 
#if 0
		if(r) Cnts[j][r] -= (UInt4) Wt[sq]; 
#else
		Cnts[j][r] -= (UInt4) Wt[sq]; 
#endif
		if(Cnts[j][r] >= UNDERFLOW_TRIGGER){
		    fprintf(stderr,"Wt[%d]=%d; Cnts[%d][%c]=%u\n",
				sq,Wt[sq],j,AlphaChar(r,AB),Cnts[j][r]);
		    FILE *efp=open_file("junk_RmFromAln_err",".cma","w");
		    PutAlnSeqsCMSA(efp,CMA); fclose(efp);
		    assert(Cnts[j][r] < UNDERFLOW_TRIGGER);  // indicative of underflow!!!
		}
	      }
	} VacateSitesCMSA(sq,CMA); ReCalcSMX=TRUE;
	return site;
}

Int4	**ssx_typ::RemoveFromAlign(Int4 *cluster)
// Remove a sequence from the alignment and return the old positions.
{
	Int4	sq,i,**site;
	NEWP(site,cluster[0] + 5, Int4);
	for(i=1; (sq=cluster[i]) != 0; i++){ site[i]=this->RemoveFromAlign(sq); }
        return site;
}

void	ssx_typ::RmFromAlign(Int4 sq)
// Remove a sequence from the alignment.
{
	Int4	st,r,j,len,blk,pos[5];
	
#if 1
	assert(SeqIsInCMSA(sq,CMA));
#endif
	NDL->RmIndels(sq,Wt[sq]);
	// JLH->RmIndels(sq,Wt[sq]);
	for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	      UInt4	**Cnts=WtCnts[blk];
	      len=LengthCMSA(blk,CMA);
	      unsigned char   *isq=SeqPtrCMSA(sq,CMA);
	      PosSiteCMSA(blk,sq,pos,CMA); st=pos[1];
	      for(j=1; j <= len; st++,j++){
		  r=isq[st]; 
#if 0
		  if(r) Cnts[j][r] -= Wt[sq]; 
#else
		  Cnts[j][r] -= Wt[sq]; 
#endif
	      }
	} VacateSitesCMSA(sq,CMA);
        ReCalcSMX=TRUE; 
}

void	ssx_typ::RmFromAlign(Int4 *cluster)
// Remove a sequence from the alignment.
{ Int4	sq,i; for(i=1; (sq=cluster[i]) != 0; i++) this->RmFromAlign(sq); }

Int4    ssx_typ::TempToWt(double Temperature)
{
        Int4    wt;
        if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
	else if(Temperature < 30) wt = 1000;
	else { 
		double	d = 300.0/Temperature; d = pow(d,3.0);
		wt=(Int4)floor(0.5+(double)d);
		assert(wt <= 1000);
	} return wt;
}

void	ssx_typ::SampleSmatrix(double pernats, Int4 blk, Int4 wt)
// return an smatrix with parameters sampled from the beta distribution ...
// Temporary for now.
{
        Int4	i,j,c,o,*observed,total,a,b;
	double	s,*likelihood,*Ps,*targetNS,totalS,obs,totalNS;
	static Int4 Seed=0;
	fm_type M=ModelCMSA(blk,CMA);
	smx_typ	smx = SMX[blk];

	if(Seed==0) { 
	   Seed=-Random();  // need to initialize with a negative number.
	   fprintf(stderr,"INITIALIZING SEED (%d)\n",Seed);
	}
	assert(wt > 0 && wt <= 20);
	Ps = M->Ps; 
	for(totalNS=0.0,c=nAlpha(M->A);c > 0; c--){totalNS+=M->observedNS[c];}
	totalNS+=M->npseudo;
	for(c=nAlpha(M->A); c > 0; c--){
	   if(Ps[c] > 0.0) M->targetNS[c]=(M->observedNS[c]+Ps[c])/totalNS;  
	   else M->targetNS[c]=1.0;  // this should not be used...
	} targetNS = M->targetNS;
	total=M->totsites;
	totalS = (double) (M->totsites + M->npseudo);
        for(i=1,j=M->start; j<=M->end; j++,i++){
	    // if(IsColumnFModel(i,M))
	    if(M->observed[j]!=NULL){
	      observed=M->observed[j];
               for(c=nAlpha(M->A); c>0; c--){
		if(Ps[c] > 0.0){
		  o=observed[c];
		  if(o > 0 && o < total){
		    a = MAXIMUM(Int4,1,o*wt),b = MAXIMUM(Int4,1,(total-o)*wt);
// std::cerr << "Debug" << c << "." << j << std::endl;
		    assert(a > 0 && b > 0);
		    // obs=(double)total*betadev(o*wt,(total-o)*wt,&Seed);
		    obs=(double)total*betadev(a,b,&Seed);
		    if(obs > 0.0){
		      s=((obs+Ps[c])/(totalS+(obs-o)))/targetNS[c];
		    } else s=((o+Ps[c])/totalS)/targetNS[c];
		  } else if(o) s=1.0/targetNS[c]; else s=(Ps[c]/totalS)/targetNS[c];
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
                } else SetSMatrix(c,i,0,smx);
	       }
	       if(Ps[0] == 0.0) SetSMatrix(0,i,0,smx);
	       else SetSMatrix(0,i,(Int4)floor((pernats*log(M->likelihood[j][0]))+0.5),smx);
             } else for(c=0; c<=nAlpha(M->A); c++) SetSMatrix(c,i,0,smx);
        } // PutFModel(stderr, M); PutSMatrix(stderr, smx); 
}

void	ssx_typ::OldSMXs(double temp)
{
	Int4    i,j,r,n,blk,wt=0; 
	if(temp==0) wt = 0; else wt=TempToWt(temp);;
        for(blk=1; blk<= nBlksCMSA(CMA); blk++){
#if 0
	   NilSMatrix(SMX[blk]);
           if(wt == 0) SMX[blk]=GetSmatrixFModel(PerNatsCMSA(CMA),ModelCMSA(blk,CMA));
           else SMX[blk]=SampleSmatrixFModel(PerNatsCMSA(CMA),wt,ModelCMSA(blk,CMA));
#elif 0
           if(wt == 0){
	   	NilSMatrix(SMX[blk]);
		SMX[blk]=GetSmatrixFModel(PerNatsCMSA(CMA),ModelCMSA(blk,CMA));
	   } else this->SampleSmatrix(pernats, blk, wt);
#else
	   if(wt == 0){
	     smx_typ smx = SMX[blk];
	     fm_type M=ModelCMSA(blk,CMA);
	     double  *Ps=PseudoFModel(M);
	     for(j=1; j<=LenFModel(M); j++){
	       if(IsColumnFModel(j,M)){
                 for(r=0; r <= nAlpha(AB); r++){
                    if(Ps[r] > 0.0){
			double score=ScoreFModel(r,j,M);
			Int4 s = (Int4) floor((pernats*score + 0.5));
			SetSMatrix(r,j,s,smx);
                    } else SetSMatrix(r,j,0,smx);
                  }
               } else for(r=0; r <= nAlpha(AB); r++) SetSMatrix(r,j,0,smx);
	     } // if(m==1) PutSMatrix(stderr,smx[m]); exit(1);
	   } else this->SampleSmatrix(pernats, blk, wt);
#endif
        }
}

#define debug_vrtl_cnts 0

void	ssx_typ::CntsFromSMX()
// score[i][r]/pernats = ln(freq_FG[i][r]/freq_BG[r])
// freq_FG[i][r]/freq_BG[r] = exp(score/pernatrs)
// freq_FG[i][r] = freq_BG * exp(score/pernatrs)
// Obs_FG[i][r] = (Total[i][r]/wt_factgor) * freq_BG * exp(score/pernatrs)
// eObs[i][r] = total_Observed[i]*BG[r]*exp(score[i][r])
// Calculate the virtual number of observed residue counts at each position.
{
	Int4	len,i,r,blk,score;
	double	dd,d,D;
	UInt4	Total,total;
	for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	    UInt4   **wtCnts=WtCnts[blk];
	    smx_typ smx=SMX[blk];
	    len=LengthCMSA(blk,CMA);
	    for(i=1; i <= len; i++){
		for(Total=0,r=1; r <= nAlpha(AB); r++) Total+=wtCnts[i][r]; 
// for(Total=0,r=1; r <= nAlpha(AB); r++) Total+=wtCnts[i][r]+WtPseudo; 
// Total+=TotalWtPseudo;
		D=(double) Total;		// total observed counts for residue r.
		for(total=0,r=1; r <= nAlpha(AB); r++){
	   		score=ValSMatrix(i,r,smx);
			dd= (double)score/(double)pernats;	// convert score into nats.
			d = BackGrnd[r]*exp(dd);	// implied observed frequency of r.
	   		VrtlCnts[blk][i][r] = (UInt4) ceil(D*d);  // round up...
// VrtlCnts[blk][i][r] *= 2;
// VrtlCnts[blk][i][r] += WtPseudo;
// if(VrtlCnts[blk][i][r] < WtPseudo) VrtlCnts[blk][i][r] = WtPseudo;
			total += VrtlCnts[blk][i][r];
			// assert(VrtlCnts[blk][i][r] <= WtN); // or a bit higher due to rounding up?
#if debug_vrtl_cnts	// this will be affected by testing_fix_for_DMS == 1
	fprintf(stderr,"score = %d; FG=%.2f; BG=%.2f; Total=%.1f; VrtlCnts[%d][%d][%c] = %.2f; Obs=%.2f\n",
			score,d,BackGrnd[r],D/(double)wt_factor,blk,i, AlphaChar(r,AB),
			(double)VrtlCnts[blk][i][r]/(double)wt_factor,
			(double)wtCnts[i][r]/(double)wt_factor);
#endif
		} TtlVrtlCnts[blk][i]=total;
#if debug_vrtl_cnts
	fprintf(stderr," *** Total Cnts = %d; Total VrtlCnts=%d; diff=%d; ratio=%.3f ***\n\n",
		Total,total,total-Total,(double)total/(double)Total);
#endif
	    }
	}
}

void	ssx_typ::CalcWtSMX(double temp)
{
	smx_typ smx; 
	Int4	score,c,i,j,r,len,st,pos[5],blk,sq,N=NumSeqsCMSA(CMA);
	UInt4	**wtCnts;
	double	*frq,*Score;

	NEW(frq,nAlpha(AB)+5,double); NEW(Score,nAlpha(AB)+5,double);
	for(blk=1; blk <= nBlksCMSA(CMA); blk++){
	    len=LengthCMSA(blk,CMA);
	    wtCnts=WtCnts[blk];
	    // NEWP(Cnts,len+5,double);
	    // for(i=1; i <= len; i++) NEW(Cnts[i],nAlpha(AB)+5,double);
	    // for(r=1; r <= nAlpha(AB); r++) 
	    smx=SMX[blk]; 
	    for(i=1; i <= len; i++){
		// if(DMS) DMS->score(wtCnts[i],Score); else psiscore(wtCnts[i],Score);
		if(dms_mode != 'P') DMS->score(wtCnts[i],Score); else psiscore(wtCnts[i],Score);
                for(c=0; c<=nAlpha(AB); c++){
#if 0
fprintf(stderr,"%d: Score[%c]=%.2f; wtCnts[%d]=%d\n", i,AlphaChar(c,AB),Score[c],i,wtCnts[i]);
#endif
		   score = (Int4) floor(Score[c] + 0.5);
		   SetSMatrix(c,i,score,smx);
		} // free(Cnts[i]);
	    } // free(Cnts);
	} free(Score); free(frq);
	CntsFromSMX();	
	// fprintf(stderr,"done calculating WtSMX\n");
}

smx_typ	*ssx_typ::SampleDirichletSMX(double temp)
/***********************************************************************
 Returns an smatrix with parameters sampled from the Dirichlet distribution
 corresponding to the observed counts with uninformed priors (i.e., with
 one pseudocount per residue type).
 ***********************************************************************/
{ return 0; }  // see fmodel.cc to adapt this to ssx_typ.

#if 0	// Stephen's code for Dirichlet deviates.
/***********************************************************************
I think the code below is all you need.  In short, if you want to sample a 
multinomial from a Dirichlet distribution with "alpha"
parameters a[i], you just sample 20 gamma variates, and normalize.
The code below samples gamma variates, but requires normal variates 
for any alphas>100.  The normal variate code requires one to set
nflag=0 before calling it the first time.  It then generates pairs 
of normal variates.  I do not remember where I found these algorithms, 
but I could probably track it down if necessary.  
Let me know if this deems to work.  Best wishes, Stephen
 ***********************************************************************/

void	DirichletDeviate(double *a, double *p, Int4 dim)
//  Sample from a Dirichlet distribution.
{
	double	sum;
	Int4	i;
        for (sum=0.0,i=0;i < dim; ++i) sum+=p[i]=gd(a[i]);
        for (i=0; i < dim; ++i) p[i]/=sum;
}


double  gd(double x)
// Gamma variate code.
{
         double  a,b,e,v,sum,bound;

         bound=100;
         if (x>bound) return(x+sqrt(x)*normal());
         for (sum=0;x>=1;x-=1) sum-= log(drand48());
         e=exp(1.0);
         v=e/(e+x);
         do {
                 if (drand48()<=v) { a=exp(log(drand48()/x)); b=exp(log(a)*(x-1)); }
                 else { a=1-log(drand48()); b=exp(-a); }
         } while (b*drand48()>exp(log(a)*(x-1)-a));
         return(sum+a);
}

int     nflag;                  /* Normal variate flag    */
double  extra;                  /* Extra normal variate   */

double  normal()
// Normal variate code 
{
         double  v,w,r,f;

         nflag=1-nflag;
         if (nflag) {
                 do { v=2*drand48()-1; w=2*drand48()-1; r=v*v+w*w; } while (r>=1);
                 f=sqrt(-2*log(r)/r); extra=v*f; return(w*f);
         } else return(extra);
}

#endif

smx_typ	*ssx_typ::WtCntsSMX( )
// return an SMX based on weighted counts plus pseudocounts...
{
    Int4 blk;
    smx_typ *rtn_smx; NEW(rtn_smx,nBlksCMSA(CMA)+3, smx_typ);
    for(blk =1; blk <= nBlksCMSA(CMA); blk++){
        Int4	i,j,c,o;
	UInt4	*observed,pseudo=1;	// 1/100 pseudocounts only
	double	s,O,d,obs,*aa,total; NEW(aa,nAlpha(AB) +3, double);

	smx_typ	smx = MkSMatrix(2.0,LenSMatrix(SMX[blk]),this->BackGrnd,AB);
        for(i=1; i <= LenSMatrix(SMX[blk]); i++){
              assert(this->WtCnts[blk][i] != NULL);
	      observed=this->WtCnts[blk][i];
              for(total=0,c=nAlpha(AB); c > 0; c--){
		// O=(double)(observed[c] + pseudo)/(double)wt_factor; o = (Int4) ceil(O); 
		aa[c]=(double)(observed[c] + pseudo); 
		total += aa[c];
	      } 
              for(c=nAlpha(AB); c>0; c--){
		  d=(double) aa[c]/(double) total;
// fprintf(stderr,"%d%c = %.1f/%.1f = %.3f\n",i,AlphaChar(c,AB),aa[c],total,d);
		  assert(d > 0.0 && d <= 1.0);
		  s=d/this->BackGrnd[c];
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
	      } SetSMatrix(0,i,0,smx); 
        } // if(blk == 1) PutSMatrix(stderr, smx); // if(TRUE) exit(1);
	rtn_smx[blk]=smx; free(aa);
    } return rtn_smx;
}

smx_typ	*ssx_typ::StraightSMX( )
// return the SMX without sampling parameters...
{
    Int4 sx;
    smx_typ *rtn_smx; NEW(rtn_smx,nBlksCMSA(CMA)+3, smx_typ);
    for(sx =1; sx <= nBlksCMSA(CMA); sx++){
        Int4	i,j,c,o,total;
	UInt4	*observed;
	double	s,O,d,totalD,obs;

	smx_typ	smx = MkSMatrix(2.0,LenSMatrix(SMX[sx]),this->BackGrnd,AB);
        for(i=1; i <= LenSMatrix(SMX[sx]); i++){
	      // totalD=(double)TtlVrtlCnts[sx][i]/(double)wt_factor;
              assert(this->VrtlCnts[sx][i] != NULL);
	      observed=this->VrtlCnts[sx][i];
	      // total = (Int4) ceil(totalD);
	      Int4 *aa; NEW(aa,nAlpha(AB) +3, Int4);
              for(total=0,c=nAlpha(AB); c > 0; c--){
		O=(double)observed[c]/(double)wt_factor; o = (Int4) ceil(O); 
		aa[c] = MAXIMUM(Int4,1,o); 
		total += aa[c];
	      } 
              for(c=nAlpha(AB); c>0; c--){
		  d=(double) aa[c]/(double) total;
		  assert(d > 0.0 && d <= 1.0);
		  s=d/this->BackGrnd[c];
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
	      } SetSMatrix(0,i,0,smx); free(aa);
        } // if(sx == 1) PutSMatrix(stderr, smx); // if(TRUE) exit(1);
	rtn_smx[sx]=smx;
    } return rtn_smx;
}

smx_typ	*ssx_typ::SampleSMX(double Temp)
// return an smatrix with parameters sampled from the beta distribution ...
{
    Int4 sx,wt=this->TempToWt(Temp);	// sample by temperature...
    smx_typ *rtn_smx; NEW(rtn_smx,nBlksCMSA(CMA)+3, smx_typ);
    for(sx =1; sx <= nBlksCMSA(CMA); sx++){
        Int4	i,j,c,o,total;
	UInt4	*observed;
	double	s,O,d,totalD,obs;
	static Int4 Seed=0;

	if(Seed==0) { 
	   Seed=-Random();  // need to initialize with a negative number.
	   fprintf(stderr,"ssx_typ::SampleSMX() INITIALIZING SEED (%d)\n",Seed);
	}
	smx_typ	smx = MkSMatrix(2.0,LenSMatrix(SMX[sx]),this->BackGrnd,AB);
        for(i=1; i <= LenSMatrix(SMX[sx]); i++){
	    totalD=(double)TtlVrtlCnts[sx][i]/(double)wt_factor;
            if(this->VrtlCnts[sx][i] != NULL){
	      observed=this->VrtlCnts[sx][i];
	      total = (Int4) ceil(totalD);
#if 1
	      Int4 *aa; NEW(aa,nAlpha(AB) +3, Int4);
	      double	tt,*pp; NEW(pp,nAlpha(AB)+3,double);
// if(i==37) fprintf(stderr,"37:    ");
              for(c=nAlpha(AB); c > 0; c--){
		O=((double)observed[c]*(double)wt)/(double)wt_factor;
		o = (Int4) ceil(O); 
		aa[c] = MAXIMUM(Int4,1,o); 
// if(i==37) fprintf(stderr,"%c%d ",AlphaChar(c,AB),o);
	      } tt=DirichletDev(aa, pp, nAlpha(AB),&Seed);
// if(i==37) fprintf(stderr,"\n");
	      //   a[c] = # of residue c; p[c] = returned probability!
              for(c=nAlpha(AB); c>0; c--){
		  d=pp[c];
		  assert(d > 0.0 && d <= 1.0);
		  // obs=(double)total*d;
		  s=d/this->BackGrnd[c];
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
		  // SetSMatrix(c,i,(Int4)floor(total*d+0.5),smx);
	      } SetSMatrix(0,i,0,smx); free(pp); free(aa);
#else 
              for(c=nAlpha(AB); c > 0; c--){
		  O=(double)observed[c]/(double)wt_factor;
		  o = (Int4) ceil(O); 
		  if(o > 0 && o <= total){
		    Int4 a = MAXIMUM(Int4,1,o*wt),b = MAXIMUM(Int4,1,(total-o)*wt);
		    d=betadev(a,b,&Seed);
// if(i==37) fprintf(stderr,"Seed = %d; beta = %.5f.; a=%d; b=%d\n",Seed,d,a,b);
		    if(d == 0.0) d = 0.000001;
		    assert(d > 0.0 && d <= 1.0);
		    // obs=(double)total*d;
		    s=d/this->BackGrnd[c];
		  } else if(o) s=1.0/this->BackGrnd[c]; 
		  else {	// all gaps in this position!!!
	       		SetSMatrix(c,i,0,smx); continue;
			fprintf(stderr,"O=%g; o=%d; total = %d\n",O,o,total);
			PutSMatrix(stderr,SMX[sx]);
			assert(o > 0);
		  }
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
	       } SetSMatrix(0,i,0,smx);
	       // SetSMatrix(0,i,-1000,smx);
#endif
             } else for(c=0; c<=nAlpha(AB); c++) SetSMatrix(c,i,0,smx);
             // else for(c=0; c<=nAlpha(AB); c++) SetSMatrix(c,i,-1000,smx);
        }
	// if(sx == 1) PutSMatrix(stderr, smx); // if(TRUE) exit(1);
	rtn_smx[sx]=smx;
    }
	return rtn_smx;
}

char    *ssx_typ::GapAlnTrace(e_type seq, double Temp, Int4 &start,Int4 &score)
{
	if(ReCalcSMX) this->CalcSMX(Temp);
	smx_typ	*smx=0;
	if(Temp >= 30.0){ smx=this->SampleSMX(Temp); } else smx=SMX;
	char	*operation;
// if(Temp < 30.0) PutSMatrix(stderr,smx[1]);
#if 0	// need to fix this; random number in sampling has issues...
	if(nBlksCMSA(CMA) == 1 && Temp >= 250.0){	// Careful: temperature not yet incorporated!!!
		operation=NDL->SampleGapAlnTrace(seq,nBlksCMSA(CMA), smx, &start,&score);
	} else 
#endif
	{
        	// operation=JLH->GapAlnTrace(seq,nBlksCMSA(CMA),smx,&start,&score);
        	operation=NDL->GapAlnTrace(seq,nBlksCMSA(CMA),smx,&start,&score);
	}
	if(smx != SMX){ for(Int4 x =1; x <= nBlksCMSA(CMA); x++) NilSMatrix(smx[x]); free(smx); }
	return operation;
}

void    ssx_typ::SampleAlignSeq(FILE *fp,char mode, e_type qE,double temp)
{
	Int4	m,wt=1,score,start;
        double	map;
	char	*operation=0;
	smx_typ	*smx=0;
        map=this->GapMap(temp);	// if any site is removed will core dump!!
	fprintf(fp,"\n=================================================\n"); 
PutSeqInfo2(fp,qE); // fprintf(fp," ==========\n");
	if(ReCalcSMX){ this->CalcSMX(); }
	if(temp > 0.0){ smx=this->SampleSMX(temp); } else smx=SMX;
        // operation=JLH->GapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score);
        operation=NDL->GapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score);
        // operation=NDL->SampleGapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score);
        // fprintf(fp,"start=%d; score=%d\n",start,score); fprintf(fp,"operation=%s\n",operation);
        put_seqaln_smatrixSW(fp,operation,LenSeq(qE)-start+1,
                SeqPtr(qE)+start-1,OffSetSeq(qE)+start-1,strlen(operation),nBlksCMSA(CMA),smx);
	if(smx != SMX){ for(Int4 x =1; x <= nBlksCMSA(CMA); x++) NilSMatrix(smx[x]); free(smx); }
        free(operation);
}

void    ssx_typ::AlignSeq(FILE *fp,char mode, e_type qE,double temp)
{
	Int4	m,wt=1,score,start;
        double	map, penalty;
	char	*operation=0;
	smx_typ	*smx=0;
#if 0	// This appears necessary to initialize JLH paramters to model.
        map=this->RelMap( );
        if(mode == 's' || mode == 'S'){ penalty=JLH->IndelPenalty(stderr,'S',CMA); }
	else if(mode=='m'){ wt=0; penalty=JLH->IndelPenalty(stderr,mode,CMA); }
	else { wt=0; penalty=JLH->IndelPenalty(stderr,'m',CMA); }
	printf("map = %g; penalty = %g; total map = %g\n",map,penalty,map+penalty);
	map=map+penalty;
#elif 1
        map=this->GapMap(temp);	// if any site is removed will core dump!!
#endif
	fprintf(fp,"=========== "); PutSeqInfo2(fp,qE); // fprintf(fp," ==========\n");
#if 1
	if(ReCalcSMX){ this->CalcSMX(); }
	if(temp > 0.0){ smx=this->SampleSMX(temp); } else smx=SMX;
        // operation=JLH->GapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score);
        operation=NDL->GapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score);
        fprintf(fp,"start=%d; score=%d\n",start,score); fprintf(fp,"operation=%s\n",operation);
        put_seqaln_smatrixSW(fp,operation,LenSeq(qE)-start+1,
                SeqPtr(qE)+start-1,OffSetSeq(qE)+start-1,strlen(operation),nBlksCMSA(CMA),smx);
#if 0	// DEBUG...
	// for(Int4 x =1; x <= nBlksCMSA(CMA); x++) PutSMatrix(stderr, smx[x]);
	PutSMatrix(stderr, smx[1]);
#endif
	if(smx != SMX){ for(Int4 x =1; x <= nBlksCMSA(CMA); x++) NilSMatrix(smx[x]); free(smx); }
#else
	if(ReCalcSMX) this->CalcSMX();
        // operation=JLH->GapAlnTrace(qE,nBlksCMSA(CMA),SMX,&start,&score);
        operation=NDL->GapAlnTrace(qE,nBlksCMSA(CMA),SMX,&start,&score);
        fprintf(fp,"start=%d; score=%d\n",start,score); fprintf(fp,"operation=%s\n",operation);
        put_seqaln_smatrixSW(fp,operation,LenSeq(qE)-start+1,
                SeqPtr(qE)+start-1,OffSetSeq(qE)+start-1,strlen(operation),nBlksCMSA(CMA),SMX);
#endif
        free(operation);
}

char	*ssx_typ::AlignSeq(e_type qE,Int4 &start,double temp,BooLean global)
{
        Int4    score;
	a_type  AB = AlphabetCMSA(CMA);
        char    *operation=0;
        double  map=this->GapMap(temp); // if any site is removed will core dump!!
        if(ReCalcSMX){ this->CalcSMX(); }
        smx_typ *smx=SMX; // if(temp > 0.0){ smx=this->SampleSMX(temp); } 
        operation=NDL->GapAlnTrace(qE,nBlksCMSA(CMA),smx,&start,&score,global);
	assert(nBlksCMSA(CMA) == 1);
	Int4	sites[3]; sites[1]=start; sites[2]=0;
#if 0
	FILE	*fp=stdout;
        fprintf(fp,"start=%d; score=%d\n",start,score); fprintf(fp,"operation=%s\n",operation);
        fprintf(fp,"=========== "); PutSeqInfo2(fp,qE); // fprintf(fp," ==========\n");
	fflush(fp);
#endif
#if 0
	gsq_typ *gsq=new gsq_typ[1];
	// gsq->initialize(nBlksCMSA(CMA),LengthsCMSA(CMA),sites,qE);
	gsq->initialize(operation,qE,sites);
	gsq->Put_cma_format(fp,0,nBlksCMSA(CMA),sites,LengthsCMSA(CMA),AB);
	delete []gsq;
        put_seqaln_smatrixSW(fp,operation,LenSeq(qE)-start+1, SeqPtr(qE)+start-1,
			OffSetSeq(qE)+start-1,strlen(operation),nBlksCMSA(CMA),smx);
	fflush(fp);
#endif
#if 0   // DEBUG...
        // for(Int4 x =1; x <= nBlksCMSA(CMA); x++) PutSMatrix(stderr, smx[x]);
        PutSMatrix(stderr, smx[1]);
#endif
        // if(smx != SMX){ for(Int4 x =1; x <= nBlksCMSA(CMA); x++) NilSMatrix(smx[x]); free(smx); }
        return operation;
}

