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

#include "hma_typ.h"

hma_typ::hma_typ(Int4 maxrpts,cma_typ cma,double *exp_rpt_gap)
// -G200..1300,20..100/0..50:500,35..250
{
	init(maxrpts,"-P200..1300,20..100:500,35..250",200, cma,exp_rpt_gap); 
	char *best_arg = GetBestTransProbPairs( );
	Free( );
	init(maxrpts,best_arg,200,cma,exp_rpt_gap);
}

char    *hma_typ::GetBestTransProbPairs( )
// Use Newton-Raphson Method to find optimum parameters.
{
	char *best_arg=0; 
	best_arg=AllocString(hmm->ArgTP( ));
	return best_arg;
}

hma_typ::hma_typ(Int4 maxrpts,char *pssm_arg,cma_typ cma,Int4 pernats,double *exp_rpt_gap,
	BooLean use_wt)
{ init(maxrpts,pssm_arg,pernats,cma,exp_rpt_gap); weightMAP=use_wt; }

hma_typ::hma_typ(Int4 maxrpts,char *pssm_arg,cma_typ cma,Int4 pernats,double *exp_rpt_gap)
{ init(maxrpts,pssm_arg,pernats,cma,exp_rpt_gap); }

void	hma_typ::init(Int4 maxrpts, char *hma_arg, Int4 pernats, cma_typ cma, 
	double *exp_rpt_gap)
{
	weightMAP=FALSE;
	NEWP(operation,NumSeqsCMSA(cma)+4,char);
	NEW(start,NumSeqsCMSA(cma)+4,Int4);
	NEW(skipseq,NumSeqsCMSA(cma)+4,BooLean);
	NEW(oper_len,NumSeqsCMSA(cma)+4,Int4);
	hmm = new HMM_typ(maxrpts,hma_arg,cma,pernats,0);
        hmm->Put(stderr);
	cmsa = cma;
	PerNats=pernats;
}

void	hma_typ::Free( )
{
	delete hmm; 
        for(Int4 i=1;i<=NumSeqsCMSA(cmsa);i++) if(operation) free(operation[i]);
	free(operation); free(start); free(oper_len); free(skipseq);
}

void	hma_typ::Put(FILE *fp)
{
}

double	hma_typ::SampleRelMap(double Temperature)
{
	Int4	wt;
	double	inc=10.0;

	if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
	else if(Temperature > 75.0){
	   wt = 1 + (Int4) floor((300.0+inc-Temperature)/inc);
	   // for inc = 10.0 --> wt for T=290-300 is 2;  280-290 is 3; 270-280 is 4; etc...
	   // down to 75-80 = 24.
	}  else wt=0;  // Use maximum likelihood.

	Int4	Score,Start,OperLen;
	double	map=0.0;
	double *weight;
	if(weightMAP) weight=HenikoffWeightsCMSA(cmsa);
	for(Int4 sq=1;sq<=NumSeqsCMSA(cmsa);sq++){
		gss_typ *gss=gssCMSA(cmsa);
		e_type E=gss->TrueSeq(sq);
		if(operation[sq] != 0){
		  Score=hmm->GetScore(E,operation[sq],start[sq],oper_len[sq]);
		} else {
		  if(wt > 0) SampleAlign(sq,1, &Score, &Start, &OperLen,wt);
		  else Align(sq,1, &Score, &Start, &OperLen);
		  // fprintf(stderr,"Score(%d) = %d\n",sq,Score);
		}
		if(weightMAP) map+=(double) Score*weight[sq]/(double) PerNats;
		else map+=(double) Score/(double) PerNats;
		// if(Score <= 0) hmm->PutAlign(stderr,E,Operation,oper_len,Start,1);
	} 
	if(weightMAP) free(weight);
	return map; 
}

double	hma_typ::RelMap(Int4 rm_seq)
{
	Int4	Score,Start,OperLen;
	double	map=0.0;
	double *weight;
	if(weightMAP) weight=HenikoffWeightsCMSA(cmsa);
	for(Int4 sq=1;sq<=NumSeqsCMSA(cmsa);sq++){
		if(sq==rm_seq) continue;
		gss_typ *gss=gssCMSA(cmsa);
		e_type E=gss->TrueSeq(sq);
		if(operation[sq] != 0){
		  Score=hmm->GetScore(E,operation[sq],start[sq],oper_len[sq]);
		} else {
		  operation[sq] = hmm->Align(E,1, &Score, &Start, &OperLen);
		  start[sq]=Start; oper_len[sq]=OperLen;
		  fprintf(stderr,"Score(%d) = %d\n",sq,Score);
		}
		if(weightMAP) map+=(double) Score*weight[sq]/(double) PerNats;
		else map+=(double) Score/(double) PerNats;
		// if(Score <= 0) hmm->PutAlign(stderr,E,Operation,oper_len,Start,1);
	}
	if(weightMAP) free(weight);
	return map; 
}

char	*hma_typ::FindBestRpts(Int4 sq, Int4 start_rpts, Int4 *Rpts, 
	Int4 *Score, Int4 *Start, Int4 *OperLen)
// Stub for now...
{
	char	*best_operation=0;
	return best_operation;
}

char	*hma_typ::Align(Int4 sq, Int4 rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len) 
{ 
	gss_typ *gss=gssCMSA(cmsa);
	e_type E=gss->TrueSeq(sq);
	if(operation[sq]) free(operation[sq]);
	operation[sq]=hmm->Align(E,rpts,Score,Start,Oper_len); 
	start[sq]=*Start; oper_len[sq]=*Oper_len;
	return operation[sq];
}

char	*hma_typ::SampleAlign(Int4 sq, Int4 Rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len,
	Int4 wt) 
// Stub for now...
{
	assert(wt != 0);
	gss_typ *gss=gssCMSA(cmsa);
	e_type E=gss->TrueSeq(sq);
	if(operation[sq]) free(operation[sq]);
	operation[sq]=hmm->SampleAlign(E,Rpts,Score,Start,Oper_len,wt);
	// operation[sq]=hmm->Align(E,Rpts,Score,Start,Oper_len);
	start[sq]=*Start; oper_len[sq]=*Oper_len;

	Int4 score2; 
	score2=hmm->GetScore(E,operation[sq],start[sq],oper_len[sq]);
#if 0
	fprintf(stderr,"score = %d; score2 = %d\n",*Score,score2);
	if(*Score != score2){
		hmm->PutAlign(stderr,E,operation[sq],oper_len[sq],start[sq],1);
		// exit(1);
	} *Score=score2;	// need right parameters for map...
#endif
	return operation[sq];
}

char	hma_comp_map_cmsa3(double rmap, double nmap, double omap, double temperature)
/** short routine to be used for sampling three different maps **/
// MOVE THIS TO BETTER LOCATION AND ELIMINATE REDUNDANCIES!!! (w/ comp_map_cmsa3)
{
	double	ratio,maxmap,minmap,midmap;
	char	max,min,mid;

	if(rmap > nmap){
	   if(rmap > omap){ 	// i.e., rmap > nmap & rmap > omap.
	     maxmap=rmap; max = 'r'; 
	     if(omap > nmap){			// r > o > n.
		midmap = omap; minmap=nmap; 
		mid = 'o'; min = 'n'; 
	     } else { 				// r > n > o.
		midmap = nmap; minmap=omap; 
		mid = 'n'; min = 'o'; 
	     }
	   } else {				// o > r > n.
		maxmap=omap; midmap=rmap; minmap=nmap; 
		max = 'o'; mid = 'r'; min = 'n';
	   }
	} else {		
	   if(nmap > omap){	// i.e., nmap > rmap & nmap > omap.
		maxmap=nmap; max = 'n'; 
		if(omap > rmap){		// n > o > r.
		   midmap = omap; minmap=rmap; 
		   mid = 'o'; min = 'r'; 
		} else {			// n > r > o.
		   midmap = rmap; minmap=omap;
		   mid='r'; min = 'o'; 
		}
	   } else {				// o > n > r.
		maxmap=omap; midmap= nmap; minmap=rmap; 
		max = 'o'; mid = 'n'; min = 'r';
	   }
	}
	if(temperature <= 50.0) return max;
	else temperature = 300./temperature;
#if 0
 fprintf(stderr,"rmap = %g; nmap = %g; omap = %g\n",rmap,nmap,omap);
 fprintf(stderr,"maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
 fprintf(stderr,"max = %c; mid = %c; min = %c\n",max,mid,min);
#endif
	if((maxmap - midmap) > 100) return max;
	else { maxmap -= minmap; midmap -= minmap; }
	if(maxmap > 100) { 		// i.e., discard minmap.
	   maxmap -= midmap;
	   if(maxmap > 100) return max;
	   else {
		ratio = exp(maxmap);	/** ratio of new to old **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		ratio = ((double)Random()/(double)RANDOM_MAX)*(ratio+1);
		if(ratio >= 1.) return max; else return mid;
	   }
	} else {
	   double total,rand;
#if 0
fprintf(stderr,"rel: maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
	   maxmap = exp(maxmap); midmap = exp(midmap); minmap = 1.0;
fprintf(stderr,"exp: maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
#endif
	   if(temperature != 1.0){
		 maxmap = pow(maxmap,temperature);
		 midmap = pow(midmap,temperature);
	   }
	   total = maxmap + midmap + 1.0;
	   rand = total*((double)Random()/(double)RANDOM_MAX);
	   // total = |------max-+----|---mid--|-|
	   // rand  = |----------^
#if 0
fprintf(stderr,"total = %g; rand = %g; maxmap = %g; midmap = %g; minmap = %g\n",
		total,rand,maxmap,midmap,minmap);
 fprintf(stderr,"max = %c; mid = %c; min = %c\n",max,mid,min);
#endif
	   if((rand -= maxmap) <= 0.0) return max; 
	   else if((rand - midmap) <= 0.0) return mid;
	   else return min;
	}
}

char	hma_typ::SampleGSeqGibbs(Int4 s, BooLean *Skipseq, double *OldMap, 
	double Temperature, double prob_accept, UInt4 min_rpt)
// For sampling gapped sequences into and out of cmsa.
{
	Int4	t,Start,trace_length,gapopen,gapextend;
	Int4    score;
	Int4	*newpos;
	a_type	A=AlphabetCMSA(cmsa);
	e_type	oldE,E;
	gss_typ	*gss=gssCMSA(cmsa);
	double	pernats;
	char	*Operation;
	Int4	wt;
	cma_typ	cma;
	double	oldmap,newmap,rm_map;
	double	inc=10.0;
	BooLean	CanDelete=FALSE;

	assert(prob_accept > 0.0);
	cma=CopyCMSA(cmsa);
	oldmap=RelMap( ); 	// calls hma->RelMap
	*OldMap = oldmap;
	gapopen=gss->GapOpen(); gapextend=gss->GapExtend(); pernats=gss->PerNats();

	if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
	else if(Temperature > 75.0){
	   wt = 1 + (Int4) floor((300.0+inc-Temperature)/inc);
	   // for inc = 10.0 --> wt for T=290-300 is 2;  280-290 is 3; 270-280 is 4; etc...
	   // down to 75-80 = 24.
	}  else wt=0;  // Use maximum likelihood.

	if((cma)->FullSeq){
	    if((cma)->FullRpts[(cma)->SubToFull[s]] > min_rpt) CanDelete=TRUE;
	    else CanDelete=FALSE;
	} else CanDelete=TRUE;
	// 1. Remove sequence s from alignment.
	NEW(newpos,nBlksCMSA(cma)+2,Int4);  // for new sites.
	if(!skipseq[s] && CanDelete){
	  VacateSitesCMSA(s,cma); 
	  if((cma)->FullSeq){ (cma)->FullRpts[(cma)->SubToFull[s]]--; }
	  rm_map=RelMap(s); 
#if 1   // NEW prob_accept
          rm_map = rm_map - prob_accept; // scale up rm_map to favor removal by prior prob.
#endif
	  if((cma)->FullSeq){ (cma)->FullRpts[(cma)->SubToFull[s]]++; }
// if(oldmap < rm_map){ fprintf(stderr,"#### oldmap = %g; rm_map = %g\n",oldmap,rm_map); }
	} else { rm_map = -200.0; }

        E = gss->TrueSeq(s); 	
        oldE = gss->FakeSeq(s); 	
	// PutSeq(stderr,E,AlphabetCMSA(cma)); fflush(stderr);

	// ===== 2. Obtain scoring matrix from rest of alignment.

	// ========== 4. Sample a gapped alignment for sequence s. ===========
	// assert(nBlksCMSA(cma) == 1);
	if(wt == 0) Operation=Align(s,1,&score,&Start,&trace_length);
	else Operation=SampleAlign(s,1,&score,&Start,&trace_length,wt);
        /// Operation is save as hma_typ operation[sq]; 
	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				Operation,trace_length,Start,E,newpos);

	// ========== 7. If sequence changed then replace s with gsq. =========
	if(!skipseq[s] && gss->Identical(s,*gsq)){ 
	  free(newpos); delete []gsq; newmap = -200.0; 
	} else {
	  ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & fmodel.
	  for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	  free(newpos);
	  hmm->ReComputeSMX(cma);
	  newmap=RelMap( ); 
	}
	// Now Sample one of these alignments.
	char option;
// fprintf(stderr,"#### oldmap = %g; newmap = %g; rm_map = %g\n",oldmap,newmap,rm_map); 
// oldmap=0; *OldMap = oldmap;
	if(Temperature<=10.0) option = hma_comp_map_cmsa3(rm_map,newmap,oldmap,0.0);
	else option = hma_comp_map_cmsa3(rm_map,newmap,oldmap,Temperature);
	switch(option){
	  case 'r': {   // remove sequence from alignment.
		Skipseq[s]=skipseq[s]=TRUE;   
		ReplaceCMSA(s,NULL,cmsa);
		if((cmsa)->FullSeq){
		  (cmsa)->FullRpts[(cmsa)->SubToFull[s]]--;
		}
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ(); NilCMSA(cma); 
	    } break;
	  case 'o': {  // retain old sequence.
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ(); NilCMSA(cma); 
	    } break;
	  case 'n': {  // convert to new sequence;
		if(skipseq[s]) {
		    Skipseq[s]=skipseq[s]=FALSE;   
		    option = 'a';
		    if(cma->FullSeq){
		      cma->FullRpts[cma->SubToFull[s]]++;
		    }
		}
		gss_typ     *gss0=gssCMSA(cmsa);
		gss0->~gss_typ(); NilCMSA(cmsa); 
		cmsa=cma; 
	    } break;
	  default:  print_error("This should not happen"); break;
	}
	return option;
}

