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

#include "gsm_typ.h"

double  gsm_typ::log_like_ratio(cma_typ cma)
/**********************************************************************
 return the relative map for aligment.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 This ignores block combinatorics.
 **********************************************************************/
{
    Int4        N,n,len,i,j,t,T,s,ncol,on_col;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(cma);
    fm_type	*model = ModelsCMSA(cma);    
    Int4	nblks=nBlksCMSA(cma);
    a_type	A;
    ss_type	truedata = TrueDataCMSA(cma);

    A = SeqSetA(data); N = NSeqsSeqSet(data);
    if(cma->FullSeq){	// then use full_counts and freq;
	counts = CntsSeqSet(cma->FullSeq); freq=tFreqSeqSet(cma->FullSeq);
    } else { // OLD: prior to full_counts.
       counts = CntsSeqSet(truedata); freq=tFreqSeqSet(truedata); 
    }
    for(len=0, t=1; t <= nblks; t++) { n=MaxLenFModel(model[t]); if(n > len) len=n; }
    NEWP(observed,len+1,Int4); npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	totsite = TotSitesFModel(model[t]);
	/***** MAP *******/
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(A); b++) {
		MAP += cma->lngamma[site[b]][b]; 
		nonsite[b] -= site[b];
	    }
	    MAP -= lngamma((double)(totsite-site[0]) + npseudo);
	    // MAP -= lngamma((double)totsite + npseudo);
	  } 
	}
    } v = lngamma(npseudo);
    /***** NONSITES MAP *******/
    if(cma->FullSeq){	// check for problems due to flanking regions...
       for(b=1;b<= nAlpha(A); b++) assert(nonsite[b] >= 0);
    }
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(Ps[b] > 0){
		MAP += lngamma((double) nonsite[b] + Ps[b]);
		total += (double) nonsite[b];
        	v -= lngamma(Ps[b]);   /** term pulling down MAP **/
	   }
    }
    MAP += v *(double) on_col; 
    MAP -= lngamma((double) total + npseudo);
    /***** subtract NULL MAP ******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){ MAP -= lngamma(counts[b] + Ps[b]); total += counts[b]; }
    } MAP += lngamma((double)total + npseudo); free(observed);
    return (MAP);
}


void	gsm_typ::purge() 
{
    Int4	n;

    if(input_msa) return; 
    sprintf(str,"%s",name); sprintf(name,"%s0",Argv[1]);
    if(UseRepset){  // then first generate representative set...
        // FILE *fp2 = open_file(name,".rep","w");
        FILE *fp2 = open_file(name,"","w");
	n=RepSetCluster(fp2,str,0,!noseg,'S',cutoff,UseLabeled,AB);
	// this should call gapped blast with score option...
	fclose(fp2);
        // sprintf(str,"%s.rep",name);
    } else n=PurgeFiles(str,name,cutoff,inc,minseq,maxseq,AB);
    if(n < minseq) print_error("purged set has too few sequences");
}

BooLean	gsm_typ::Refine()
// Not currently used....
{
	BooLean improved=FALSE;
	cma_typ	maD;

	fprintf(stderr,"delete blocks?\n");
        while(nBlksCMSA(cmsa) > 1 && (maD=RefineCMSA(cmsa))!=NULL){
	     fprintf(stderr,"Block deleted\n");
  	     NilCMSA(cmsa); cmsa=maD; 
#if 0
	     gibbs(NULL); improved = TRUE;
#endif
#if 1
	  if(use_gseq) strcpy(options,"-t1 -g "); else strcpy(options,"-t1 ");
  	  map = core_gibbs(&cmsa);
	  // fprintf(stderr," refined map %.1f\n",map);
  	  Record( ); improved=TRUE; 
#endif
	}
        return improved;
}

BooLean	gsm_typ::Record( )
// Record current alignment data.
{
	FILE	*fp;

	if(0) fprintf(stderr,"@map = %.3f; %d blocks; %d columns\n",
		   map, nBlksCMSA(cmsa), NumColumnsCMSA(cmsa)); 
	if(WG && cmsa){
		fp=open_file(name,".see","w");
		AddWatchGibbs(map,nBlksCMSA(cmsa),NumColumnsCMSA(cmsa),WG);
		PutWatchGibbs(fp, WG); fclose(fp);
	}
	if(map > bestmap){
	        // fprintf(stderr,"!!!bestmap improved from %.1f to %.1f\n",bestmap,map);
		bestmap=map; if(write) PutAlnCMSA(name,cmsa,Guide); 
		if(bestmsa) NilCMSA(bestmsa); bestmsa=CopyCMSA(cmsa);
		return TRUE;
	} return FALSE;
}

double  gsm_typ::core_gibbs(cma_typ *M, char mode, double temp, double minmap)
{
        double  map;
        gs_type G;
        cma_typ ma=*M;
        char    *Opt[100];
        Int4    i,nopt=string2argv(Opt,options); assert(nopt < 100);

        assert(ma != NULL);
        G = blank_gibbs( ); G->cmsa=ma;
        G->stage1_minmap=minmap;
        OptionsGibbs(nopt,Opt,G);
        switch (mode){
           case 'N': if(AlignedCMSA(ma)) G->nruns=0; // if ma->best != NULL: save previous best alignment.
                G->mod_temp=FALSE; G->temp0=0; G->temp=1.0; break;
           case 'D': G->nruns=-1; break; // discard previous bestmsa.
           case 'd': G->nruns=-1; break; // discard previous bestmsa.
           case 'S': G->nruns=0; break; // don't discard previous bestmsa.
           case 'H': G->mod_temp=FALSE; G->nruns=-1; G->temp0 = 0;
                G->mod_temp=FALSE; // don't change temperature during sampling
            break;      // discard previous bestmsa.
           default: print_error("SimulatedAnnealingGibbs( ) input error");
        }
        if(temp >= 0.0){ assert(SetTemperatureGibbs(temp,G)); }
	if(mode == 'N' || mode == 'd') G->mod_temp=FALSE; 
	map = this->the_gibbs_sampler(G);
        if(mode == 'H'){ SaveBestCMSA(G->cmsa); map = RelMapCMSA(G->cmsa); }
        else if(map > 0.0) InitMAPGibbs(G);
        *M = NilGibbs(G);
	for(i=0; i<nopt; i++) free(Opt[i]); 
        return map;
}

BooLean	gsm_typ::gibbs(cma_typ ma2,char tweakmode)
{
	cma_typ	recombinant;
	double	map2;
	Int4	b,c,j;
	BooLean	improved = FALSE;
	FILE 	*efp=0;

	if(ma2 != NULL){	//============ Not used with gismo++ !! =============
	    // print_error("gibbs() cma input not used");
	    if(use_gseq) strcpy(options,"-t1 -g ");
	    else strcpy(options,"-t1 ");
  	    if(useSA) map2=core_gibbs(&ma2,'S',150);
  	    else map2 = core_gibbs(&ma2);
	    if((recombinant = RecombineCMSA(ma2, cmsa)) != NULL){
  		if(useSA) map2=core_gibbs(&recombinant,'S',150);
  		else map2 = RelMapCMSA(recombinant);
//	    	else map2 = core_gibbs(&recombinant);
		fprintf(stderr,"!!!recombinant map = %.2f\n", map2);
		NilCMSA(ma2); ma2 = recombinant; 
	    }
	    if(map2 > map){
  		if(!useSA) map2=core_gibbs(&ma2,'S',150);
	        fprintf(stderr,"!!!map improved from %.1f to %.1f\n", map,map2);
		NilCMSA(cmsa); cmsa=ma2; map=map2; 
		Record( ); return TRUE;
	    } else { 
		fprintf(stderr,"map drops ap from %.1f to %.1f\n",map,map2);
		NilCMSA(ma2); 
		return FALSE;
	    }
        } else if(cmsa!=NULL){	//=========== simply tweak the cmsa ==============
// print_error("gibbs() cmsa input not used"); // used by gismo++ !!!
	  switch(tweakmode){
	   case 'N': // use gapped...
	  	if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
	        else sprintf(options,"-t1 -l%d ",limit);
		map = core_gibbs(&cmsa); break; 
	   case 'n': // Normal (Fast) mode
	        sprintf(options,"-t1 -l%d ",limit);
		map = core_gibbs(&cmsa); break;
	   case 'S': // Simulated Annealing mode
  	    	map = core_gibbs(&cmsa,'S',temperature); break;
	   default: 
  	    	if(useSA) map = core_gibbs(&cmsa,'S',temperature);
  	    	else map = core_gibbs(&cmsa);
	    break;
	  }
	  // fprintf(stderr," refined map %.1f (%.1f K)\n",map,temperature);
  	  if(Record( )) improved=TRUE; 
#if 0	// experimental gap function...
if(use_gseq){
  fprintf(stderr,"performing gapped SA version of gibbs...\n");
  map2=GappedSimAnnealCMSA(&cmsa, 200.0, 75.0, 5.0);
  if(map2 > map){ fprintf(stderr,"!!!map improved from %.1f to %.1f\n", map,map2); map=map2; }
  if(Record( )) improved=TRUE; 
}
#endif
	} else {	//================= return a new cmsa... ===============
	 // char	other_options[]="-T1";
	 // char	other_options[]="-T1 -m10 ";	// -m == nconverge.
	 char	other_options[]=" ";
	 if(Run==1 && cmsa_in != NULL){
print_error("gibbs() cmsa_in input not used");
		// sprintf(options,"-t1 -T1 "); // no gaps (-g) for initial alignment...
		if(use_gseq) sprintf(options,"-t1 -g -l%d %s ",limit,other_options); 
		else sprintf(options,"-t1 -l%d %s ",limit,other_options);
		over=1; under=0;
		// b=ConfigCMSA2(num,cmsa_in); c=0;
		b=SampleConfig( ); c=0;
		fprintf(stderr,"%d: %d blocks\n",cycle,b); 
		// for(j=1; j<= b; j++) fprintf(stderr," block %d: %d columns\n",j,num[j]); 
	 } else {
	
	  if(cycle == 1){
	   if(use_gseq) sprintf(options,"-t1 -l%d -g %s ",limit,other_options);
	   else sprintf(options,"-t1 -l%d %s ",limit,other_options);
	   BLK=aveblk; COL=avecol;
	   if(fix){ BLK = MAXIMUM(Int4,BLK,minblk); BLK = MINIMUM(Int4,BLK,maxblk); }
	   assert(BLK > 0); assert(aveblk==BLK);
	   assert(COL <= maxlenModels);
	   assert(COL <= MAX_LENG_SITES);
	   b = SampleConfig( );
	  } else {  // cycle != 1..

	    if(use_gseq) sprintf(options,"-t1 -g -l%d %s ",limit,other_options);
	    else sprintf(options,"-t1 -l%d %s ",limit,other_options);
// std::cerr << options;  std::cerr << std::endl;
	    if(cycle > 15){ over=1; under=0; }
	    else if(cycle > 5){ over=2; under=1; } else { over=3; under=2; }
	    b = SampleConfig();
	  }	// end cycle != 1.
	 }	// end of 'else // i.e., not Run ==1 ...'
	 for(j=1; j<= b; j++){
		// fprintf(stderr," block %d: %d columns\n",j,num[j]); 
		if(num[j] > MAX_LENG_SITES){
		    if(cmsa_in != 0){
	 		PutConfigCMSA(stderr,cmsa_in);
			FILE *Xfp=open_file("core_dump",".cma","w");
			PutCMSA(Xfp,cmsa_in); fclose(Xfp);
		    } assert(num[j] <= MAX_LENG_SITES);
		}
	 }
	 if(StrtLens) cmsa=RandomCMSA(StrtBlks,StrtLens,*gss); 
	 else cmsa=RandomCMSA(b,num, *gss);	// WARNING: model == NULL.
#if 0
	if(TestMode){
	 // SaveBestCMSA(cmsa);	// save this new configuration as the best...
	 return TRUE;
	}
#endif
	 // PutConfigCMSA(stderr,cmsa);
// std::cerr << options; std::cerr << std::endl;
	 // map=RunMinMapGibbs(options, bestmap, &cmsa);  // Quit if not promising...
#if 1	// new routine...
	BooLean	found;
	Int4 iter=0; // MaxIter=4; by default == five rounds maximum!
	SmplBlksCols=TRUE;
	char mode='S';
	double Temp=100;
	do {
	  found=FALSE;
	  if(efp) PutConfigCMSA(stderr,cmsa);
	  if(iter==0){
	     // InitCMSA(cmsa);
	     map = core_gibbs(&cmsa,'N',150); 
	     if(map < 5.0) break; 
	     if(efp) fprintf(stderr,"iter %d of core_gibbs; map=%.3f; %d blks; %d cols.\n",
			iter,map,nBlksCMSA(cmsa),NumColumnsCMSA(cmsa));
	     if(efp) PutConfigCMSA(stderr,cmsa); iter++;
	  } if(AlignedCMSA(cmsa)){ InitMAPCMSA(cmsa); } // save the best...
	  for(b=1; b <= nBlksCMSA(cmsa); b++){
	      if(nBlksCMSA(cmsa) >= MAX_NO_TYPESITES) break;
	      Int4 lft_len,len=LengthCMSA(b,cmsa);
	      if(len >= 12){
	        found=TRUE; lft_len=len/2;
	        cma_typ rcma=SplitBlkCMSA(b,lft_len,5,cmsa);
	        assert(rcma != 0); NilCMSA(cmsa); cmsa=rcma;
	  	if(efp) PutConfigCMSA(stderr,cmsa); b--;
	      }
  	  }
	  if(found) SaveBestCMSA(cmsa);	// save this new configuration as the best...
	  if(iter < 2) found=TRUE;  // go at least 2 iterations.
	  if(iter==MaxIter || !found) this->SmplBlksCols=FALSE;
	  // map = core_gibbs(&cmsa,'N',Temp); this->SmplBlksCols=TRUE;
	  map = core_gibbs(&cmsa); this->SmplBlksCols=TRUE;
	  if(map < 5.0) break;
	  if(efp) fprintf(stderr,"iter %d of core_gibbs; map=%.3f; %d blks; %d cols.\n",
			iter,map,nBlksCMSA(cmsa),NumColumnsCMSA(cmsa));
	  if(efp) PutConfigCMSA(stderr,cmsa);
	  Temp = MAXIMUM(double,Temp-50,0.0);
  	  if(iter >= MaxIter) break;	
	  iter++;
	} while(found);
#else
	 map = core_gibbs(&cmsa);
#endif
	 // PutConfigCMSA(stderr,cmsa);
	 if(map <= 5.0) { NilCMSA(cmsa); cmsa = NULL; return FALSE; }
  	 Record( ); return TRUE; 
	}	// end of if(ma2==cmsa==null).
	return improved;
}

Int4	gsm_typ::SampleConfig()
{
	Int4	b,c;
	if(num) free(num); NEW(num,aveblk*2+2,Int4);
        for(c=0,b=1; b <= aveblk; b++) c+=num[b]=avecol;
	fprintf(stderr,"%d: %d blocks %d columns (ave=%d; COL=%d)\n",cycle,aveblk,c,avecol,COL); 
	return aveblk;
}

cma_typ	*gsm_typ::RtnBestCMA() 
// FINAL OPTIMIZATION OF BEST ALIGNMENTs.
{
     Int4	i,j,min;
     cma_typ	ma2,ma3,*rcma; NEW(rcma,nMSAHeap(maH)+3,cma_typ);
     double	map2;

     if(maH==NULL) return 0;
     if(nMSAHeap(maH) == 0) return 0;
     if(use_gseq) sprintf(options,"-t1 -g "); else sprintf(options,"-t1 ");

   for(i=1; !EmptyMSAHeap(maH); i++){
     // assert(bestmsa != NULL); min = nBlksCMSA(bestmsa);
     // if(!fix){ minblk = MAXIMUM(Int4,min-under,2); maxblk = min+over;}
     cmsa = DelBestMSAHeap(&map2, maH); // NilMSAHeap(maH); maH=NULL;
     // Record( );
     map = core_gibbs(&cmsa,'S',150);
     Record( ); // if(map > bestmap){ bestmap=map; PutAlnCMSA(name,cmsa,Guide); }
     rcma[i]=cmsa; cmsa=0;
   } return rcma;
}

#define gismo_no_speedup 0

cma_typ	*gsm_typ::Align(Int4 ab,Int4 ac,Int4 mnb,Int4 mxb,Int4 len_mode)
// for gismo++: ab=aveblk; ac=avecol.
{
   Int4		i,j,n,sum=0,maxfutile=1,rtn; 
   double	spacing;
   FILE		*efp=0;
   
  BLK=0; read_input_data(); 
#if 1	//**************** randomize parameter settings... ***************
     n = 15*5*100;	// ac; sp; ab.
     const Int4 HpSz=8000,HpSzPlus=HpSz+9;
     dh_type	dH=dheap(HpSzPlus,4);
     Int4	tries,J,xcol[HpSzPlus],xblk[HpSzPlus];
     double	xspc[HpSzPlus];
  for(i=0; i < HpSz; ){
     for(ac=5; ac <= 15; ac++) {
	for(spacing=2.0; spacing < 4.1; spacing +=0.5){	// 2 2.5 3 3.5 4 = 5 settings.
           ab=(Int4) ceil((double)len_mode/(spacing*(double)ac));
	   ab = MINIMUM(Int4,ab,100);	// no more than 100 settings.
           if(ab > mxb) ab=mxb;
	   for(Int4 b=ab; b >= 3; b--){
		i++; xblk[i]=b; xcol[i]=ac; xspc[i]=spacing; 
		// fprintf(stderr,"### %d: spc=%.1f; ab=%d; ac=%d\n",i,spacing,b,ac);
		insrtHeap(i,(keytyp)Random(),dH); 
	        if(i >= HpSz) break;
	   } if(i >= HpSz) break;
	} if(i >= HpSz) break;
     }
  }
  // maxcycle = 2;
  maxcycle = 1;
  for(sum=0,tries=0,rtn=1; !emptyHeap(dH) && sum < mhpsz; ){
     assert((J=delminHeap(dH)) != 0); ac=xcol[J]; spacing=xspc[J]; ab=xblk[J];
     do {
      this->ResetParameters(ab,ac,mnb,mxb); Run++; 
      fprintf(stderr,"************ %d: create_cma() *******************\n",Run);
      rtn=create_cma(n,maxfutile); tries++; 	// tries to find maxcycle # msa
	// sum+=n; 
      if(n > 0){
         sum+=1; 
	 fprintf(stderr,"### %d: n=%d; ab=%d; ac=%d; sum=%d; maxcycle=%d; mhpsz=%d; futile=%d; rtn=%d\n",
		i,n,ab,ac,sum,maxcycle,mhpsz,maxfutile,rtn);
     	 { if(efp) PutConfigCMSA(stderr,SeeMSAHeap(BestItemMSAheap(maH),maH)); }
      }
      if(rtn == 1){		// this configuration failed...try another...
	if((J=delminHeap(dH)) == 0){
#if 1	// let calling environment decide...
	     assert(sum < mhpsz); free_input_data(); return 0;
#else
	     print_error("No conserved regions found in the input sequences");
#endif
	} assert(J != 0); ac=xcol[J]; spacing=xspc[J]; ab=xblk[J];
      } else if(J != 0){	// this configuration succeeded...try it again later.
	insrtHeap(J,(keytyp)Random(),dH); xblk[J]=ab; xcol[J]=ac; xspc[J]=spacing; J=0; 
      }
      // if(tries > 500 && sum == 0) break; // print_error("No conserved regions found in the input sequences");
     } while(rtn == 1);
     // if(tries > 500 && sum == 0) break; // print_error("No conserved regions found in the input sequences");
     // if(tries > 8000) break;
  } Nildheap(dH);
#if 1
  if(sum < mhpsz){ free_input_data(); return 0; }
#else
  if(sum == 0){ free_input_data(); return 0; }
#endif
#else	//**************** original method... **********************8
  for(rtn=1,ac=5; ac <= 15 && rtn==1; ac++)
  // for(rtn=1,ac=15; ac >= 5 && rtn==1; ac--)
  // for(rtn=1,j=1,ac=10; ac >= 5 && ac <= 15 && rtn==1; j++)
  {
   for(spacing=2.0; spacing < 4.1 && rtn==1; spacing +=0.5){	// spacing 2...3...4
     i=0; ab=(Int4) ceil((double)len_mode/(spacing*(double)ac));
     if(ab > mxb) ab=mxb;
     do {
      this->ResetParameters(ab,ac,mnb,mxb); Run++;
      rtn=create_cma(n,maxfutile); sum+=n;
      if(n > 0){
	 fprintf(stderr,"### %d: n=%d; ab=%d; ac=%d; sum=%d; maxcycle=%d; mhpsz=%d; futile=%d; rtn=%d\n",
		i,n,ab,ac,sum,maxcycle,mhpsz,maxfutile,rtn);
     	 if(sum > 0){
		if(efp) PutConfigCMSA(stderr,SeeMSAHeap(BestItemMSAheap(maH),maH)); 
#if 0	// Debug
	    for(Int4 x=1; x <= SizeMSAHeap(maH); x++){
		if(SeeMSAHeap(x,maH) == 0) continue;
		cma_typ xcma=OneBlockCMSA(SeeMSAHeap(x,maH)); 
		WriteCMSA("junk.cma",xcma); NilCMSA(xcma);
	    }
#endif
	 }
      }
      if(rtn == 1){ i++; 
	if(n < 1 && sum < mhpsz){
		ab--; if(ab < 3) ab=3; else i=0; if(ac > 15) ac=15; }
	else if(n >= mhpsz) i=0; 
      } 
      if(sum > maxcycle){ rtn=0;  break; }	// found enough alignments.
      // if(i > 10) break;
      if(i > 5) break;
     } while(rtn == 1);
   }
// break;
   // if(j%2 == 0){ ac += j; }	//  9+2=11; 8+4=12; 7+6=13; 6+8=14; 5+10=15;
   // else { ac -= j; }	      // 10-1=9; 11-3=8; 12-5=7; 13-7=6; 14-9=5; stop.
  }
  if(rtn==1 && sum < mhpsz){ free_input_data(); return 0; }
#endif	//******************* end original method ***************
#if gismo_no_speedup
  this->RefineMSAHeap();
#endif
  // if(Delete('n')) gibbs(NULL);
  cma_typ *rcma=this->RtnBestCMA();
#if 0
  for(i=1; rcma[i]; i++){ core_gibbs(&rcma[i],'S',150); }
#endif
  return rcma;
}

Int4	gsm_typ::create_cma(Int4 &num_success,Int4 max_futile)
// BUILD INITIAL POPULATION (MSAHEAP).
{
   Int4		failed_attempts=0;
   UInt4	time0,time0b;
   cma_typ	xcma;
   double	avelen=0.0,dd,d;
   BooLean	improved;

   num_success=0; item=NULL; time0b=time(NULL); 
   if(cmsa != NULL){ NilCMSA(cmsa); cmsa= NULL; }
   for(cycle=1; cycle <= maxcycle; cycle++){
	time0=time(NULL); 
	cmsa=NULL; gibbs(NULL);
	if(cmsa){ InitMAPCMSA(cmsa); map=UnGappedRelMapCMSA(cmsa); Record( ); }
	if(cmsa && map > 0) {  // THEN ADD ALIGNMENT TO HEAP.
#if gismo_no_speedup
	    if(!dont_delete) improved=Delete('n');  // Delete and tweak alignment.
#endif
	    // if(improved) gibbs(NULL);
   	    while(this->DeletePoor()){
#if gismo_no_speedup
		map=core_gibbs(&cmsa,'S',150); // GISMO speedup: turned off afn2_4_2015.
#endif
		InitMAPCMSA(cmsa); map=UnGappedRelMapCMSA(cmsa);
	     	Record( ); // if(map > bestmap){ bestmap=map; PutAlnCMSA(name,cmsa,Guide); }
	    }
#if 1
	    if(bestmsa != 0 && bestmsa != cmsa){
		NilCMSA(cmsa); cmsa=bestmsa; bestmsa=0; bestmap=-9999999999.;
	    }
#endif
	    // PutConfigCMSA(stderr,cmsa); 
	    // if(1) dd=LogLikeRatio(cmsa); else dd=map;
	    dd=UnGappedRelMapCMSA(cmsa);
	    if((item=InsertMSAHeap(cmsa,dd,maH))==0){ NilCMSA(cmsa); } cmsa=0; 
   	    num_success++; ConvergedMSAHeap(maH);
	    if(BestItemMSAheap(maH) != item){ // i.e., new optimum alignment.
		item = BestItemMSAheap(maH);
		xcma=SeeMSAHeap(item,maH);
	    	aveblk = 3 + nBlksCMSA(xcma); 
#if 0
		avecol=(Int4) floor((double)TotalLenCMSA(xcma)/(double) nBlksCMSA(xcma));
		if(avecol < 5) avecol=5; else if(avecol > 10) avecol=10;
#endif
	    } failed_attempts=0;
	} else {	// found a lousy alignment, so delete it.
		if(cmsa!=NULL) NilCMSA(cmsa); cmsa=NULL; 
		cycle--; failed_attempts++; 
		d = (double)num_success/(double)(num_success + failed_attempts);
		if(failed_attempts > max_futile || d < 0.25){
			return 1; // print_error("failed to find a motif");
		} else if(failed_attempts % 10 == 0) { 
		   if(!fix){
			minblk--; maxblk--; minblk = MAXIMUM(Int4,minblk,1); 
			if(maxblk < minblk) maxblk = minblk;
		   }
		}
	} if(0) fprintf(stderr,"\ttime: %d seconds (%0.2f minutes)\n",
                time(NULL)-time0,(float)(time(NULL)-time0)/60.0);
     }
     fprintf(stderr,"\tTotal time: %d seconds (%0.2f minutes)\n",
                time(NULL)-time0b,(float)(time(NULL)-time0b)/60.0);
     return 0;
}

void	gsm_typ::RefineMSAHeap()
// DO SOME POST PROCESSING HERE...
// GISMO+ speed up: turned this off afn 2_4_2015.
{
     cma_typ	*list,xcma;
     double	*lpr;
     Int4	t,aln,NumAln=nMSAHeap(maH);
     NEW(list,NumAln+3,cma_typ);
     NEW(lpr,NumAln+3,double);
     for(aln=1; nMSAHeap(maH) > 0; aln++){ 
	assert((xcma=DelMinMSAHeap(&map, maH)) != NULL);
fprintf(stderr,"******************** refining alignment %d *********************\n",aln);
	// for(t=1; t<=nBlksCMSA(xcma); t++) avelen+=LengthCMSA(t,xcma);
	// avelen/=(double)nBlksCMSA(xcma);
	if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
	else sprintf(options,"-t1 -l%d ",limit);
        // map=core_gibbs(&xcma,'S',300);
        map=core_gibbs(&xcma); // cmsa=xcma; Record( ); cmsa=0;
	InitMAPCMSA(xcma); map=UnGappedRelMapCMSA(xcma);
	list[aln]=xcma; lpr[aln]=map;
     }
     for(aln=1; aln <= NumAln; aln++){ assert(InsertMSAHeap(list[aln],lpr[aln],maH)!=NULL); }
     free(list); free(lpr);
}

void gsm_typ::Optimize() 
// FINAL OPTIMIZATION OF BEST ALIGNMENT.
{
   Int4		j,min;
   cma_typ	ma2,ma3;
   double	map2;

   if(maH==NULL) print_error("gsm_typ::Optimize() input error");
   if(use_gseq) sprintf(options,"-t1 -g "); else sprintf(options,"-t1 ");

   // fprintf(stderr,"******** breed ********\n"); map=Breed( );

   while(nMSAHeap(maH) > 1){
        if((ma2 = DelMaxMSAHeap(&map, maH)) != NULL) NilCMSA(ma2);
   } assert(bestmsa != NULL); min = nBlksCMSA(bestmsa);
   if(!fix){ minblk = MAXIMUM(Int4,min-under,2); maxblk = min+over;}
   fprintf(stderr,"--> End breed ");
   double runtime=difftime(time(NULL),time1);
   fprintf(stderr," time thus far: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
   cmsa = DelMinMSAHeap(&map2, maH); // NilMSAHeap(maH); maH=NULL;
   // Note: map inserted as -map!!!
   Record( );
#if 0
   if(use_gseq) strcpy(options,"-t1 -g "); else strcpy(options,"-t1 ");
   ma2 = CopyCMSA(cmsa);
   fprintf(stderr,"starting temperature = %.2f K\n", 300.);
   for(j=1; j<=10; j++){
     map2=core_gibbs(&ma2,'S',300);
     fprintf(stderr,"map = %.2f\n", map2);
     if((ma3 = RecombineCMSA(cmsa, ma2)) != NULL){
   	map2=core_gibbs(&ma3,'S',150);
        NilCMSA(ma2); ma2=ma3; ma3=NULL;
     } 
     if(map2 > bestmap){ 
	NilCMSA(cmsa); cmsa=ma2; bestmap=map2; PutAlnCMSA(name,cmsa,Guide); 
	ma2=NULL; break;
     } 
   }
   if(ma2 != NULL) NilCMSA(ma2);
#endif
     map = core_gibbs(&cmsa,'S',150);
     Record( ); // if(map > bestmap){ bestmap=map; PutAlnCMSA(name,cmsa,Guide); }
}

Int4	gsm_typ::align()
{
   Int4 num_success;
   BLK=0; read_input_data(); Run++; // maH=MkMSAHeap(mhpsz);  // created within gsm_init.cc
#if 1
   if(create_cma(num_success) == 1){ free_input_data(); return 1; }
   // if(Delete('n')) gibbs(NULL);
   Optimize();	// removes cma's from maH.
#else
   create_population(); Recombine(); optimize();
#endif
   if(bestmsa != NULL && write) PutAlnCMSA(name,bestmsa,Guide);
   free_input_data(); return 0;
}

void gsm_typ::gapped_align()
// this is not used...
{
   double	d;
   BLK=0;
   read_input_data();
   fprintf(stderr,"blocks columns .. MAP\n"); Run++;
   create_population();
   cmsa = DelMinMSAHeap(&d, maH); NilMSAHeap(maH); maH=NULL;
   free_input_data();
}

cma_typ	gsm_typ::Improve(cma_typ cma0)
{ 
  char		Str[300];
  double	Temperature=200;

  if(cma0 == 0){ sprintf(Str,"%s.cma",name); cma0=ReadCMSA2(Str,AB); }
  SetPenaltyCMSA(gapopen,gapextend,cma0);
  sprintf(name,"%s%d",Argv[1],Run);  
  FILE *fptr=open_file(name,"","w");
  PutSeqSetEs(fptr,TrueDataCMSA(cma0));
  fclose(fptr);
  assert(cmsa==0); cmsa = cma0; map = RelMapCMSA(cmsa);
  // gss_typ *gssX=gssCMSA(cmsa); make_binomials(gssX->FakeSqSet());
  Record( ); cmsa=0; // sets bestmsa to cmsa;
  if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
  else sprintf(options,"-t1 -l%d ",limit);

	h_type H = Histogram("contributions to map",-1000,1000,1.0);
	fprintf(stderr,"full MAP = %.2f\n",map);
	Int4 N = NSeqsSeqSet(DataCMSA(bestmsa));
	double minus_map,change_map;
	for(Int4 s=1; s<=N; s++){
		minus_map = RelMapMinusSeqCMSA(s,bestmsa);
		change_map=minus_map-map;
#if 0
		fprintf(stderr,"minusMAP(%d) = %.2f (diff = %.2f)\n",
			s,minus_map,change_map);
#endif
		IncdHist(change_map,H);
		if(change_map <= improve_cut) PutSeqSetE(stdout,s,TrueDataCMSA(bestmsa));
		// if(change_map > 20) PutSeqSetE(stdout,s,TrueDataCMSA(bestmsa));
	}
	PutHist(stderr,60,H); NilHist(H);
  return cmsa;	// need to free this up in calling environment.
}

//============= New CDD -R=<in_cma> routines below... ====================
BooLean	gsm_typ::SplitUp(Int4 minlength)
{
	BooLean improved=FALSE;
	double	netlen,p;
	cma_typ	maS;
#if 1
	dh_type	dH=dheap(nBlksCMSA(cmsa)+4,3);
	Int4	i,len,end=nBlksCMSA(cmsa),t,blk;
	for(t=1; t <= end; t++){ insrtHeap(t,-(keytyp)LengthCMSA(t,cmsa),dH); }
#endif
map=RelMapCMSA(cmsa);

gss_typ *gssX=gssCMSA(cmsa);
fprintf(stderr,"cmsa: go=%d; gx=%d; pernats=%d\n",gssX->GapOpen(),gssX->GapExtend(),gssX->PerNats());
WriteMtfCMSA("junkC", cmsa, NULL);
	fprintf(stderr,"Try splitting blocks...(%d total blks)\n",nBlksCMSA(cmsa));
	// for(Int4 end=nBlksCMSA(cmsa), t=1; t <= end; t++)
	for(i=1; (blk=delminHeap(dH)) != NULL; i++)
        {
	  if(LengthCMSA(blk,cmsa) < minlength) continue;
	  netlen= (double) (LengthCMSA(blk,cmsa) - 3);
	  p = (netlen*netlen + 1.0)/400.0;
	  // len=8: p = 0.065; len=12: p = 0.20; len=18: p = 0.56; len=22: p = 1.0.
	  if(1 || netlen >= minlength && p >= SampleUniformProb( )){
	   if((maS=SplitBlkCMSA(blk, minlength, cmsa))!= NULL){
	     fprintf(stderr,"******** block %d split in two.\n",blk);
	     double	map2;
             char    *arguments[100];
	     // strcpy(options,"-t1 -g "); 
	     // strcpy(options,"-t3 -g "); 
	     strcpy(options,"-t3 "); 
             Int4    nopt=string2argv(arguments,options); assert(nopt < 100);
// temperature=300.0;
             map2=CoreGibbs(nopt,arguments,&maS,'N',temperature,0.0);
map2=RelMapCMSA(maS); gssX=gssCMSA(maS);
fprintf(stderr,"maS: go=%d; gx=%d; pernats=%d\n",gssX->GapOpen(),gssX->GapExtend(),gssX->PerNats());
	fprintf(stderr,"map changes from %.1f to %.1f\n",map,map2);
WriteMtfCMSA("junkS", maS, NULL);
exit(1);
             for(Int4 j=0; j<nopt; j++) free(arguments[j]); 
	     if(map2 > map){
		fprintf(stderr,"map improves from %.1f to %.1f\n",map,map2);
		if(cmsa) NilCMSA(cmsa); cmsa=maS; map=map2; 
  	        Record( );
		Nildheap(dH); end=nBlksCMSA(cmsa); dH=dheap(end+4,3);
		for(t=1; t <= end; t++){ insrtHeap(t,-(keytyp)LengthCMSA(t,cmsa),dH); }
	    	improved=TRUE; 
	     } else { 
		fprintf(stderr,"map drops ap from %.1f to %.1f\n",map,map2); NilCMSA(maS); 
	     }
           }
	  }
	} Nildheap(dH);
	return improved;
}

void	gsm_typ::RefineInputAln(cma_typ cma)
// assume that this was created using MAPGAPS.
{
	BooLean	improved;
	assert(cmsa==NULL);
	cmsa=cma;
	map=RelMapCMSA(cma);
	data=TrueDataCMSA(cmsa);
	gss=gssCMSA(cmsa);
	SetIndelPenaltySeqSet(indel_penalty,data);
	// make_binomials(data);
	do {
		improved=SplitUp(8);
		if(improved) fprintf(stderr,"******** improved (now %d total blks)\n",nBlksCMSA(cmsa));
		else fprintf(stderr,"failed to improve\n");
		if(Delete()){
		   fprintf(stderr,"******** Deletion improves (now %d total blks)\n",nBlksCMSA(cmsa));
		}
		// gibbs(0);
	} while(improved);
	do { } while(Delete());
	fprintf(stderr,"Running full sampler...\n");
	Gibbs(0);
}

Int4	gsm_typ::Align2(Int4 ab,Int4 ac,Int4 mnb,Int4 mxb,Int4 len_mode)
// for gismo++: ab=aveblk; ac=avecol.
{
   Int4		i,j,n,sum=0,maxfutile=1,rtn; 
   double	spacing;
   
  BLK=0; read_input_data(); 
  for(rtn=1,ac=5; ac <= 15 && rtn==1; ac++){
   for(spacing=2.0; spacing < 4.1 && rtn==1; spacing +=0.5){	// spacing 2...3...4
     i=0; ab=(Int4) ceil((double)len_mode/(spacing*(double)ac));
     if(ab > mxb) ab=mxb;
     do {
      this->ResetParameters(ab,ac,mnb,mxb); Run++;
      rtn=create_cma(n,maxfutile); sum+=n;
      if(n > 0){
	 fprintf(stderr,"### %d: n=%d; ab=%d; ac=%d; sum=%d; maxcycle=%d; mhpsz=%d; futile=%d; rtn=%d\n",
		i,n,ab,ac,sum,maxcycle,mhpsz,maxfutile,rtn);
     	 if(sum > 0){
		PutConfigCMSA(stderr,SeeMSAHeap(BestItemMSAheap(maH),maH)); 
#if 0	// Debug
	    for(Int4 x=1; x <= SizeMSAHeap(maH); x++){
		if(SeeMSAHeap(x,maH) == 0) continue;
		cma_typ xcma=OneBlockCMSA(SeeMSAHeap(x,maH)); 
		WriteCMSA("junk.cma",xcma); NilCMSA(xcma);
	    }
#endif
	 }
      }
      if(rtn == 1){ i++; 
	if(n < 1 && sum < mhpsz){
		ab--; if(ab < 3) ab=3; else i=0; if(ac > 15) ac=15; }
	else if(n >= mhpsz) i=0; 
      } 
      if(sum > maxcycle){ rtn=0;  break; }	// found enough alignments.
      // if(i > 10) break;
      if(i > 5) break;
     } while(rtn == 1);
   }
  } if(rtn==1 && sum < mhpsz){ free_input_data(); return 1; }
  // if(Delete('n')) gibbs(NULL);
  Optimize();	// empties maH and sets cmsa = best on heap...
  // if(bestmsa != NULL) PutAlnCMSA(name,bestmsa,Guide);
  free_input_data(); 
  return 0;
}

