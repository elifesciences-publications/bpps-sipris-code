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

// Gibbs Sampling algorithms for finding multiple sites in multiple sequences 
#include "gibbs.h"

gs_type	blank_gibbs( )
// *************************** initialize *******************************
{
	gs_type G;

	NEW(G,1,gibbs_sampler_type); 
	G->verbose = FALSE; 
        G->nruns=10; G->nconverge = 500; G->limit = 12;
        for(Int4 n=0; n<10;  n++) G->test[n]=FALSE;
        G->matrix=NULL; G->end = NULL; G->strt=NULL; 
	G->oldsite=NULL; G->temp=1.0; G->mode='n';
	// modes = N,n=Normal; S,s = Simulated Annealing; Z,z = Zero degree sampling.
	G->hot_gibbs=FALSE; G->gapped_gibbs=FALSE;
	G->stage1_minmap=0.0;
	return G;
}

double	RunGibbs(char *options,cma_typ *M){ return SimAnnealGibbs(options,M,'N',-1.0);}

double	HotGibbs(char *options,cma_typ *M,double temp)
{ return SimAnnealGibbs(options,M,'H',temp);}

double	SimAnnealGibbs(char *options, cma_typ *M, char mode, double temp)
{
	double	map;
	char	*arguments[100];
	Int4	nopt=string2argv(arguments,options); assert(nopt < 100);
	map=SimulatedAnnealingGibbs(nopt, arguments, M, mode, temp); 
	for(Int4 i=0; i<nopt; i++) free(arguments[i]); return map;
}

double	RunGibbs2(Int4 nopt, char *options[],cma_typ *M)
{ return SimulatedAnnealingGibbs(nopt, options, M, 'N', -1.0); }

double	SimulatedAnnealingGibbs(Int4 nopt, char *options[], cma_typ *M, char mode,
	double temp)
{ return CoreGibbs(nopt, options, M, mode, temp, 0.0); }

double	RunMinMapGibbs(char *options,double minmap, cma_typ *M)
{ 
	double	map;
	char	*arguments[100];
	Int4	nopt=string2argv(arguments,options); assert(nopt < 100);
	map=CoreGibbs(nopt, arguments, M, 'N', -1.0, minmap);
	for(Int4 i=0; i<nopt; i++) free(arguments[i]); return map;
}

double	CoreGibbs(Int4 nopt, char *options[], cma_typ *M, char mode,
	double temp, double minmap)
{
	double	map;
	gs_type G;
	cma_typ	ma=*M;

	assert(ma != NULL);
	G = blank_gibbs( ); G->cmsa=ma; 
	G->stage1_minmap=minmap;
	OptionsGibbs(nopt,options,G);
	switch (mode){  
	   case 'N': if(AlignedCMSA(ma)) G->nruns=0; // save previous best alignment.
		G->mod_temp=FALSE; G->temp0=0; G->temp=1.0; break;
	   case 'D': G->nruns=-1; break; // discard previous bestmsa. 
	   case 'S': G->nruns=0; break; // don't discard previous bestmsa.
	   case 'H': G->mod_temp=FALSE; G->nruns=-1; G->temp0 = 0;
		G->mod_temp=FALSE; // don't change temperature during sampling
	    break;	// discard previous bestmsa.
	   default: print_error("SimulatedAnnealingGibbs( ) input error");
	}
	if(temp >= 0.0) assert(SetTemperatureGibbs(temp,G));
        map = RunPropagate(stderr, G); 
	if(mode == 'H'){ SaveBestCMSA(G->cmsa); map = RelMapCMSA(G->cmsa); }
	else if(map > 0.0) InitMAPGibbs(G); 
        *M = NilGibbs(G); 
	return map;
}

BooLean	SetTemperatureGibbs(double T, gs_type G)
{
	if(T > 5000.) return FALSE;
	if(T <= 75.0) { // set right to zero degrees
		G->mod_temp = 2; G->temp0 = T; G->temp = 300./G->temp0; 
	} else { G->mod_temp=TRUE; G->temp0 = T; G->temp = 300./T; }
	return TRUE;
}

#define GIBBS_USAGE "\nUsage: gibbs file lengths [options] \n\
  lengths = <long>[,<long>]: lengths of elements for each type\n\
  <long>[,<long>] = numbers for each element (e.g. \"10,12,15\" for 3 types)\n\
  options:\n\
\t-D          - print sequence information\n\
\t-d          - DON'T use fragmented model (shut off??)\n\
\t-g          - gapped mode\n\
\t-l<int>     - set rapid convergence limit (higher = longer to converge)\n\
\t-m<int>     - set maximum number of cycles in each run \n\
\t-T<int>     - test modes\n\
\t-t<int>     - maximum number of sampling runs\n\
\t-v          - verbose mode\n\
\n"

BooLean	OptionsGibbs(Int4 argc, char *argv[], gs_type G)
{
   Int4		a,i,t,n,N;
   BooLean	flag=FALSE;

   for(a=0; a < argc; a++){
       if(argv[a][0] != '-') print_error(GIBBS_USAGE);
	switch(argv[a][1]) {
	   case 'D': flag=TRUE; break;
	   case 'd': std::cerr << "gibbs -d option is shut off!\n"; break;
	   case 'g': G->gapped_gibbs = TRUE; break;
	   case 'l': G->limit=IntOption(argv[a],'l',1,1000,GIBBS_USAGE); break;
	   case 'm': G->nconverge=IntOption(argv[a],'m',1,10000,GIBBS_USAGE); break;
	   case 'T': i=IntOption(argv[a],'T',0,10,GIBBS_USAGE); G->test[i]=TRUE; break;
	   case 't': G->nruns=IntOption(argv[a],'t',1,10000,GIBBS_USAGE);
		break;
	   case 'v': G->verbose=TRUE; break;
	   default: print_error(GIBBS_USAGE);
	}
    }
    if(flag){
	for(a = 0; a < argc; a++) printf("%s ",argv[a]); printf("\n");
	PutIDsCMSA(stdout,G->cmsa);
    }
    N = NumSeqsCMSA(G->cmsa);
    // G->maxlen=MaxSeqCMSA(G->cmsa);
    G->maxlen=MaxTrueSeqCMSA(G->cmsa);
    NEW(G->pos,(2*G->maxlen)+2,Int4);
    NEWP(G->matrix,nBlksCMSA(G->cmsa)+3,double);
    NEWP(G->end,nBlksCMSA(G->cmsa)+3,Int4);
    NEWP(G->strt,nBlksCMSA(G->cmsa)+3,Int4);
    NEW(G->oldsite,nBlksCMSA(G->cmsa)+3,Int4);
    for(t=0; t <= nBlksCMSA(G->cmsa); t++) {
	   NEW(G->matrix[t],(2*G->maxlen)+3,double);
	   NEW(G->end[t],N+3,Int4); NEW(G->strt[t],N+3,Int4);
    }
    for(i=0; i<=G->maxlen; i++) G->matrix[0][i] = 1.0;
    return TRUE;
}

BooLean EnlargeGibbs(gs_type G)
// if sequence set lengths change then may need to enlarge arrays.
{
    if(G->maxlen >= MaxSeqCMSA(G->cmsa)) return FALSE;
#if 0
std::cerr << "Enlarging gibbs\n";
fprintf(stderr,"oldmaxlen=%d; newmaxlen=%d\n",G->maxlen,MaxSeqCMSA(G->cmsa));
#endif
    // G->maxlen=MaxSeqCMSA(G->cmsa);
    G->maxlen=MaxTrueSeqCMSA(G->cmsa);
    free(G->pos); NEW(G->pos,(2*G->maxlen)+2,Int4);
    for(Int4 t=0; t <= nBlksCMSA(G->cmsa); t++) {
	   free(G->matrix[t]); NEW(G->matrix[t],(2*G->maxlen)+3,double);
    }
    for(Int4 i=0; i<=G->maxlen; i++) G->matrix[0][i] = 1.0;
    return TRUE;
}

cma_typ	NilGibbs(gs_type G)
/** Destroy Gibbs structure but return cmsa  **/
{
	cma_typ	M;
	Int4	t,n,N,ntyps;

    	N = NumSeqsCMSA(G->cmsa);
	ntyps = nBlksCMSA(G->cmsa);
	free(G->pos);
	if(G->matrix != NULL){ 
    		for(t=0;t<=ntyps;t++) free(G->matrix[t]); free(G->matrix);
	}
	if(G->end != NULL){
    		for(t=0;t<=ntyps; t++) free(G->end[t]); free(G->end);
	}
	if(G->oldsite != NULL) free(G->oldsite);
	if(G->strt != NULL){
    		for(t=0;t<=ntyps; t++) free(G->strt[t]); free(G->strt);
	} 
	M = G->cmsa; free(G);
	return M;
}

void	InitMAPGibbs(gs_type G)
{ InitMAPCMSA(G->cmsa); update_propagate_gibbs(G); }

cma_typ	GibbsCMSA(gs_type G) { return G->cmsa;  }

double	TotLikeGibbs(gs_type G)
{
	double	d,e;

	d = RelMapCMSA(G->cmsa); 
#if 0
	e = DirichletRelMap(G->cmsa);
	fprintf(stderr,"Map = %g; Blsm62MAP = %g\n", d,e);
#endif
	return d;
}

double	GibbsSampler(BooLean (*Update)(char, gs_type), 
	Int4 (*Sample)(Int4, double *, BooLean *,gs_type), gs_type G)
{
	Int4	n,iter,i,j,r,s,t,N=NumSeqsCMSA(G->cmsa);
	Int4	nruns=G->nruns,nconverg=G->nconverge;
	Int4	run,number,num_fail;
	double	TotLike,best,lastbest,map,map2,T;
	BooLean	flag;
	dh_type dH;
	char	stage;
	/**** DEBUG *****/
	static FILE *fp=NULL;
	BooLean	debug=FALSE;
	UInt4	time1;
	Int4	*ncols,c;
	/**** DEBUG *****/
	BooLean	*moved;
	Int4	*sq_moves,sqmv,tot_moves,cmvN,tot_try;
	double	*prob,prb,max_prb,min_prb,ave_prb,max_max_prb;
	double	oldMap,newMap;
	keytyp	key;

#if 0
   if(debug) {
	NEW(ncols, nBlksCMSA(G->cmsa)+2, Int4);
	time1=time(NULL);
	fp = open_file(NameSeqSetCMSA(G->cmsa),".prn","a");
   }
#endif

   //========================== Set up memory ============================
   // fprintf(stderr,"nconv = %d\n", nconverg); 
   assert(N < INT4_MAX && N > 0);
   map=lastbest=-DBL_MAX; 
   Update('c',G);
   dH = dheap(N,3);
   NEW(sq_moves,N+2,Int4); NEW(prob,N+2,double);
   mem_typ memS(N,'S'); // Sites: number of actions == number of sequences.
   mem_typ memG(N,'G'); // Gaps: number of actions == number of gapped sequences
#if 0
fprintf(stderr,"N=%d; N*5 = %d; MAX_USHRT=%d\n",N,N*5,USHRT_MAX);
   mem_typ memB(N,N*5,nBlksCMSA(G->cmsa),0.10,'B'); // BLOCKS: # actions == # blocks.
#elif 0	// N should equal the max number of columnts? Need to understand this better. afn: 3-31-2014.
	assert(MaxTrueSeqCMSA(G->cmsa) < (Int4) USHRT_MAX);
	unsigned short x=(unsigned short)MaxTrueSeqCMSA(G->cmsa);
   mem_typ memB(x*2,x*10,nBlksCMSA(G->cmsa),0.10,'B'); // BLOCKS: # actions == # blocks.
#else
   unsigned short x=0,X=0;
   if(N < (Int4) (USHRT_MAX-20)) X= (unsigned short) N; else X = (USHRT_MAX-20);
   X = MINIMUM(unsigned short,X,nBlksCMSA(G->cmsa)*1000);	// no need to make this so huge!
   x = MAXIMUM(unsigned short,X/5,5);
   mem_typ memB(nBlksCMSA(G->cmsa),'B',x,X,0.10); // BLOCKS: # actions == # blocks.
   // mem_typ memB(nBlksCMSA(G->cmsa),'B',5,35,0.10); // BLOCKS: # actions == # blocks.
#endif
   mem_typ memC(5,'C',3,15,0.33); // Columns: 5 types of operations.
   memC.Clear(4); // set cardinality of these to zero.
   moved = new BooLean [nBlksCMSA(G->cmsa)+2]; 

   //========================== Perform sampling ============================
   for(run = 1; (nruns <= 0 || run <= nruns); run++){
#if 0
	if(debug) {
		if(run<=0) fprintf(fp,"************* tweaking **************\n");  
		else fprintf(fp,"************* new run **************\n");  
	}
#endif
	fprintf(stderr,"."); best=-DBL_MAX; 
	if(nruns <= 0) { 
	    iter=0; G->stage=stage=1; 
            cmvN=NumColumnsCMSA(G->cmsa); 
	    SetPseudoToMapCMSA(G->cmsa);
	    if(nruns == 0){
		best=map=lastbest=TotLikeGibbs(G); TotLike=0.0; 
		SaveBestCMSA(G->cmsa); // don't discard previous map for this one
	    } else { map=0.01; TotLike=0.0; iter=-2; SaveBestCMSA(G->cmsa); }
	} else {
		map=TotLike=0.0; 
		Update('i',G); G->stage=stage=0; iter=-3; cmvN=nBlksCMSA(G->cmsa); 
	}
	Update('p',G); flag=FALSE;
	max_max_prb=-999999999.0; 
/********************* MAIN Gibbs sampling loop ************************/
	for(number=num_fail=0; iter <= nconverg; iter++) {
/********************* SLIDE COLUMNS SAMPLING ROUTINES ************************/
	  if(map > 0.0){  // SLIDE COLUMNS SAMPLING ROUTINES.
	  // if(map > 0.0 || !memC.Empty(4)){  // at start only.
	    for(t=1; t < nBlksCMSA(G->cmsa); t++){
		    newMap=oldMap=0.0;
		    T = 1.0;
		    // T = 0.5;
		    // if(iter > 5) T = G->temp; else T = 0.5;
		    if(dfsSlideColLtCMSA(t,&oldMap,&newMap,G->cmsa,T)){
			memC.Success(4); flag=TRUE;
         		if(newMap > best){ number = 0; best = newMap; }
		    } else {
			memC.Failure(4); 
         		if(oldMap > best){ number = 0; best = oldMap; }
         	        if((nruns >= 0 || stage > 1) && map < best){
			    map=best; SaveBestCMSA(G->cmsa); 
			}
		        newMap=oldMap=0.0;
		       if(dfsSlideColRtCMSA(t,&oldMap,&newMap,G->cmsa,T)){
			  memC.Success(4); flag=TRUE;
         		  if(newMap > best){ number = 0; best=newMap; }
		       } else {
			  memC.Failure(4); 
         		  if(oldMap > best){ number = 0; best=oldMap; }
		       }
		    }
         	    if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
	    }
	  }
/********************* COLUMNS SAMPLING ROUTINES ************************/
	  if(map > 0.0){
	     for(n = 1; n <= cmvN; n++){
		// BooLean success=FALSE; 
   if(!G->test[1]){
	        do { t=random_integer(nBlksCMSA(G->cmsa))+1;
	        } while(memB.Chance(t) < SampleUniformProb());
		switch(random_integer(4)){
		  case 0: 
		    if(ContigBlkCMSA(t,G->cmsa)){
			if(ShiftCMSA(G->cmsa,t)){ memC.Success(1); flag=TRUE; }
			else { memC.Failure(1); }
		    } break;
		  case 1:
		    if(SampAddColTempCMSA(t,G->temp,G->cmsa)) {
				memC.Success(2); flag=TRUE; 
		    } else { memC.Failure(2); }
		    break;
		  case 2: // do more often.
		    if(SampRmColTempCMSA(t,G->temp,G->cmsa)){
				memC.Success(3); flag=TRUE; 
		    } else { memC.Failure(3); }
		    break;
		  default: 
		    if(SampleEndColCMSA(t,G->temp,G->cmsa)){
				memC.Success(5); EnlargeGibbs(G); flag=TRUE;
		    } else { memC.Failure(5); }
		    break;
		}
    }
#if 0
if(G->test[2]) {
   switch(random_integer(4)){
     case 0:
	if(MoveColumnCMSA(G->cmsa,t)) { memC.Success(2); flag=TRUE; }
	else memC.Failure(2); 
      break;
     default: 
	if(nBlksCMSA(G->cmsa) > 1){
	  if(TransferColumnCMSA(G->cmsa)){ memC.Success(2); flag=TRUE; }
	  else memC.Failure(2);
	}
   }
}
#endif
	    }
	  } else for(t = 1; t <= nBlksCMSA(G->cmsa); t++){
	      if(ContigBlkCMSA(t,G->cmsa) && memC.Chance(1) >= SampleUniformProb()){
		if(ShiftCMSA(G->cmsa,t)){ memC.Success(1); flag=TRUE; }
	      }
	  }
	  if(flag) { Update('p',G); flag=FALSE; }
// for(t=1; t<=nBlksCMSA(G->cmsa); t++){ std::cerr << "BLOCKS:";  memB.Put(stderr,t); }
// for(i=1; i<=4; i++) if(memC.Chance(i) > 0.5){ std::cerr << "COLUMNS:";  memC.Put(stderr,i); }
// if(G->test[1]) for(i=1; i<=4; i++)if(!memC.Empty(i)){ std::cerr << "COLUMNS:";  memC.Put(stderr,i); }
if(debug) if(stage > 0 && !memC.Empty(5)){ std::cerr << "END COLUMNS:";  memC.Put(stderr,5); }

/********************* GAPPED-BASED Gibbs SAMPLING ************************/
	  if(G->gapped_gibbs && ( (stage > 0 && map > 0.0 && !G->mod_temp)
		|| (G->mod_temp && stage > 1) )){
if(debug) std::cerr << "starting gapped sampling loop\n";
       	    for(n=1; n <= N; n++) insrtHeap(n,(keytyp) Random(),dH); 
       	    while((n=delminHeap(dH)) != NULL){ 
	     if(memG.Chance(n) >= SampleUniformProb()){
		double previous=TotLike;
		BooLean gapworked;
		if(G->mod_temp && G->temp0 > 0)
		     gapworked=SampleGapsGibbsCMSA(n,&G->cmsa,G->temp0);
		else gapworked=SampleGapsGibbsCMSA(n,&G->cmsa,300);
		if(gapworked){
			// WARNING: always uses zero temp right now...
         		EnlargeGibbs(G); Update('p',G);
         	     	if((TotLike=TotLikeGibbs(G)) > previous){ 
			    memG.Success(n);
if(debug) std::cerr << "gapped sampling worked...\n"; 
			} else memG.Failure(n);
         	     	if(TotLike > best){ 
			  num_fail=number=0; best=TotLike; 
         		  if((nruns >= 0 || stage > 1) && map < best)
				{ map=best; SaveBestCMSA(G->cmsa); }
	                } else {
			  num_fail++;
			  if(map > 0 && TotLike < 0.8*map){  // then start over...
			   if(AlignedCMSA(G->cmsa)){
if(debug) std::cerr << "starting over with previous best alignment...\n\n"; 
				InitMAPCMSA(G->cmsa); Update('p',G);
			   }
			  }
	                }
if(debug){ std::cerr << TotLike; std::cerr << "; map = "; std::cerr << map; std::cerr << std::endl; }
// std::cerr << TotLike; std::cerr << "; map = "; std::cerr << RelMapCMSA(G->cmsa); 
// std::cerr << "; newmap = "; std::cerr << RelMapCMSA2(G->cmsa); std::cerr << std::endl; 
	        } else { num_fail++; memG.Failure(n); }
	     }
// if(TRUE || memG.Chance(n) <= 0.50) { std::cerr << "GAPS:";  memG.Put(stderr,n); }
	    }
	  }
if(G->mod_temp && G->temp0 > 0){
		// if(G->temp0 > 75.) SetTemperatureGibbs(G->temp0-1.0, G);
		if(G->temp0 > 75.) SetTemperatureGibbs(G->temp0-2.0, G);
fprintf(stderr,"*************** temperature = %2f K **************\n", G->temp0);
} else {
/********************* block-based Gibbs sampling ************************/
	  max_prb=-999999999.0; min_prb=+999999999.0;
   	  tot_try=tot_moves=0;
       	  for(n = 1; n <= N; n++) insrtHeap(n,(keytyp) Random(),dH); 
       	  while((n=delminHeap(dH)) != NULL){ 
	    if(memS.Chance(n) >= SampleUniformProb()){
         	if((sqmv=Sample(n,&prb,moved,G)) > 0){
		   for(t=1; t <= nBlksCMSA(G->cmsa); t++){
			if(moved[t]) memB.Success(t);
			else memB.Failure(t);
// memB.Put(stderr,t); 
		   } memS.Success(n);
		} else memS.Failure(n);
// if(memS.Chance(n) <= 0.50) memS.Put(stderr,n); 
// memS.Put(stderr,n); 
		if(max_prb < prb){
         	    if((max_prb=prb) > max_max_prb){
			max_max_prb = max_prb;
         	     	   if((TotLike=TotLikeGibbs(G)) > best){ 
         			num_fail=number=0; best=TotLike; 
         			if((nruns >= 0 || stage > 1) && map < best){ 
					map=best; SaveBestCMSA(G->cmsa); 
if(G->mod_temp && G->temp0 > 0){
   fprintf(stderr,"map = %.2f @ %.2f K (%.2f degrees F)\n",
                	map,G->temp0,((9.*G->temp0)/5.)-459);
}
				}
	                   } else {
			     num_fail++;
			     if(map > 0 && TotLike < 0.8*map){  // then start over...
			      if(AlignedCMSA(G->cmsa)){
if(debug) std::cerr << "starting over with previous best alignment...\n\n"; 
				InitMAPCMSA(G->cmsa); Update('p',G);
			      }
			     }
	                   }
         	    }
         	}
         	min_prb=MINIMUM(double,min_prb,prb);
         	sq_moves[n] = sqmv; prob[n] = prb;
         	tot_moves+=sqmv; tot_try+=nBlksCMSA(G->cmsa);
	    }
          }
       	  ave_prb=((max_prb-min_prb)/4.0) + min_prb;
}
/********************* checkpoint for map ************************/
if(G->mod_temp && G->temp0 > 0){
		if(G->temp0 > 75.) SetTemperatureGibbs(G->temp0-2.0, G);
		// if(G->temp0 > 75.) SetTemperatureGibbs(G->temp0-1.0, G);
fprintf(stderr,"*************** temperature = %2f K **************\n", G->temp0);
}
	  if(iter > 0){ // simulated annealing option 
	     TotLike = TotLikeGibbs(G);
if(debug) { std::cerr << TotLike; std::cerr << "; map = "; std::cerr << map; std::cerr << std::endl; }
#if 0
fprintf(stderr,
"** stage %d: blocks = %d; number = %d; TotLike = %g \n",stage,nBlksCMSA(G->cmsa), number,TotLike);
fprintf(stderr,
"  max = %.2f; min %.2f; ave_prob = %.2f\n", max_prb,min_prb,ave_prb);
fprintf(stderr,
" %.2f seq_moves (%d/%d)\n",
			(double)tot_moves/(double)tot_try,tot_moves,tot_try);
#endif
	     if(TotLike > best){ 
#if 0
	       if(debug) {
		ConfigCMSA(ncols, G->cmsa);
		fprintf(fp,"%4d  %5.2f ",iter,TotLike);
		for(c=0, t=1; t<=nBlksCMSA(G->cmsa); t++){	 
			fprintf(fp,"%6.1f ",FieldRelMapCMSA(G->cmsa,t));
			c+=ncols[t];
		}
		fprintf(fp,"%3d ",c);
		for(t=1; t<=nBlksCMSA(G->cmsa); t++){ fprintf(fp,"%3d ",ncols[t]); }
		fprintf(fp,"\n");
	       }
#endif
		num_fail=number=0; best=TotLike; 
		// if(stage > 0 && map < best){ map=best; SaveBestCMSA(G->cmsa); }
		if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
	     } else {
		num_fail++;
		if(TotLike < 0.8*map){  // then start over...
		 if(AlignedCMSA(G->cmsa)){
if(debug) std::cerr << "starting over with previous best alignment...\n\n"; 
		  InitMAPCMSA(G->cmsa); Update('p',G);
		 }
		}
	     }
	  } 
	  if(debug) fprintf(stderr,"********** num_fail = %d *********\n",num_fail);
/********************* check for convergence ************************/
	  if(iter <= 0) number=0;
	  // if(++number > G->limit || num_fail > NumSeqsCMSA(G->cmsa)) {
	  if(++number > G->limit) {
		if(stage == 3) {
if(debug) { std::cerr << map; std::cerr << " = stage3 map\n\n";  }
			G->stage++; stage++; break; 
		} else if(stage == 2) { 
			// G->stage++; stage++; break; 
if(debug){ std::cerr << map; std::cerr << " = stage2 map\n\n";  }
			G->stage++; stage++; number=G->limit-3; G->mod_temp=2; G->temp0=0;
		} else if(stage==1) { G->stage=stage=2; 
			number=G->limit-3; 
if(debug){ std::cerr << map; std::cerr << " = stage1 map\n\n";  }
			// max_max_prb=-9999999999999.; map = 0.0;
			// number=0; // best=map=0.0;
		} else {  /** stage == 0 **/
		   G->stage=stage=1; number=0; 
		   cmvN=NumColumnsCMSA(G->cmsa);
		   SetPseudoToMapCMSA(G->cmsa);
		   best=TotLikeGibbs(G); 
		   if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
if(debug){ std::cerr << map; std::cerr << " = stage0 map\n\n"; }
		   if(map < G->stage1_minmap){ map=0.0; break; }
if(nruns > 0 && G->test[1]) break;
                   if(best < 0){ /** THEN ABORT **/ break; }
		}
	  }
#if 0
	  if(nruns > 0 && G->test[1]){
		if(memC.Empty(1) && memC.Empty(2) && memC.Empty(3)){
		   SetPseudoToMapCMSA(G->cmsa);
		   best=TotLikeGibbs(G); 
		   if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
		   break; 
		}
	  }
#endif
	}
	if(best==lastbest || nruns <= 0) break; 
	if(lastbest < best) lastbest = best;
    }
if(G->mod_temp && G->temp0 > 0) fprintf(stderr,"final temperature = %2f K\n", G->temp0);
    Nildheap(dH); free(sq_moves); free(prob);
    fprintf(stderr," %.1f\n", map);
#if 0
    if(debug) {
	fclose(fp); free(ncols);
	if(nruns > 0){
	    fp = open_file(NameSeqSetCMSA(G->cmsa),".log","a");
	    fprintf(fp,"%5.1f %3d %3ld %3d %3d\n",
		map,iter,time(NULL)-time1, c,nBlksCMSA(G->cmsa));
	    fclose(fp);
  	}
    }
#endif
    delete [] moved;
    return map;
}

/*********************** propagation routines *****************************/
double	RunPropagate(FILE *fptr, gs_type G)
{ return GibbsSampler(UpdatePropagate, propagate_gibbs, G); }

BooLean UpdatePropagate(char c, gs_type G)
{
   Int4         k,n,t,N;

   switch(c){
     case 'p': /** update for propagation **/
        update_propagate_gibbs(G); return TRUE;
     case 'c': /** check data structure **/
	if(!SamePerSeqCMSA(G->cmsa)) print_error("ordering error");
        return TRUE;
     case 'i': /** initialize gibbs **/
      InitCMSA(G->cmsa); update_propagate_gibbs(G); return TRUE;
     case 's': /** shift occurred **/
      return FALSE;
     default: return FALSE;  /** command ignored **/
  }
}

void	update_propagate_gibbs(gs_type G)
/*** find ends of matrix ***/
{
	Int4	n,t,blocked;

	for(blocked=1, t=1; t<=nBlksCMSA(G->cmsa); t++) {
	    for(n = 1; n <= NumSeqsCMSA(G->cmsa); n++) {
		G->strt[t][n] = blocked;
	    }
	    blocked += LengthCMSA(t,G->cmsa);
	}
	for(blocked=-1, t=nBlksCMSA(G->cmsa); t > 0; t--) {
	    blocked += LengthCMSA(t,G->cmsa);
	    for(n = 1; n <= NumSeqsCMSA(G->cmsa); n++) {
		G->end[t][n] = LenSeqCMSA(n,G->cmsa) - blocked;
	    }
	}
	for(t=nBlksCMSA(G->cmsa); t > 0; t--) {
	    for(n=1; n <= NumSeqsCMSA(G->cmsa); n++) {
		if(G->strt[t][n] > G->end[t][n]){
		  fprintf(stderr,
			"t=%d;n=%d;strt=%d;end=%d;seqlen=%d;modelen=%d\n",
			t,n,G->strt[t][n],G->end[t][n],
			LenSeqCMSA(n,G->cmsa), LengthCMSA(t,G->cmsa));
		  print_error("update_propagate_gibbs( ): seq. too short");
		}
	    }
	}
}

Int4	propagate_gibbs(Int4 n, double *L, BooLean *moved, gs_type G)
// *L = likelihood; * 
{
	Int4	end,k,t,s,ntyp,*oldsite=G->oldsite;
        double  p,q,sum;
	Int4	moves = 0;

	ntyp = nBlksCMSA(G->cmsa);
	// 1: remove sites from sequence n.
	for(t = 1; t <= ntyp; t++) {
	   PosSiteCMSA(t,n,G->pos,G->cmsa);
	   s = G->pos[1]; RmSiteCMSA(t,n,s, G->cmsa);
	   oldsite[t]=s;
	}

	// 2: compute conditional probabilities.
        for(p=0.0,t=1; t<=ntyp; t++) {
             sum=CondProbPropagate(t,n,G); 
             q=G->matrix[t][G->end[t][n]]; 
	     p+=log(q)+log(sum);
        } *L = p;

	// 3: sample backwards through the matrix adding sampled 
	//	sites back into sequence.
	for(s=G->end[ntyp][n], t = ntyp; t > 0;  t--) {
	     if(G->mod_temp == 2) s=BestSitePropagate(s,t,n,G);
	     else s=ChooseSitePropagate(s,t,n,G); 
	     if(s != oldsite[t]) { moved[t]=TRUE; moves++; } else moved[t]=FALSE;
	     AddSiteCMSA(t,n,s, G->cmsa);
	     if(t > 1) s -= LengthCMSA(t-1,G->cmsa);
	}
	return moves;
}

Int4     ChooseSitePropagate(register Int4 end, Int4 t, Int4 n, gs_type G)
// sample a t site in sequence n of S. return site location
// WARNING: Assumes that BLOCKED sites have zero probability. 
{
        register double	rand_no, cum_prob;
        register Int4	site;
        register double	*P=G->matrix[t];

        if(end <= G->end[t][n]){
          do {
		rand_no = (double) Random()/(double) RANDOM_MAX;
          } while(rand_no == 0.0);
          rand_no *= P[end];	/** cumulative probability **/
          for(site=G->strt[t][n],cum_prob = 0.0; site <= end; site++){
           if((cum_prob += (double) (P[site]-P[site-1])) >= rand_no)
			return site;
          }
	} else print_error("ChooseSitePropagate( ): end error");
        fprintf(stderr,
		"site = %d; n = %d; start = %d; end = %d (=%d?)\n", 
		site,n,G->strt[t][n],G->end[t][n],end);
        fprintf(stderr,
		"block %d; total prob = %g; rand_no = %g; cum_prob = %g\n",
                t,P[end],rand_no,cum_prob);
	for(site=G->strt[t][n]; site <=  G->end[t][n]; site++){
		fprintf(stderr," s = %d; P[s] = %g\n",site,P[site]);
	}
	// PutMSA(stderr, G->cmsa);
	PutSeqAlnCMSA(stderr, G->cmsa);
	assert(cum_prob < rand_no);
        print_error("ChooseSiteGibbs( ) - this should not happen!?");
}

Int4	BestSitePropagate(register Int4 end, Int4 t, Int4 n, gs_type G)
// Pick the best type t site in sequence n of S. return site location
//  WARNING: Assumes that BLOCKED sites have zero probability. 
{
        register double	max_prob;
        register Int4	max_site,site;
        register double	*P = G->matrix[t];

        if(end <= G->end[t][n]){
           for(max_site=site=G->strt[t][n],max_prob=0.; site<=end; site++){
		if(max_prob < (double) (P[site]-P[site-1])){
			max_prob=P[site]-P[site-1]; max_site=site;
		}
           }
           return max_site;
	} print_error("BestSiteGibbs( ) - this should not happen!?");
}

double	LogLikePropagate(gs_type G)
{
	Int4	end,k,t,s,ntyp,n;
	double	p,q,sum;

	ntyp = nBlksCMSA(G->cmsa);
	for(p=0,n = 1; n <= NumSeqsCMSA(G->cmsa); n++){
	   for(t=1; t<=ntyp; t++) {
	     sum=CondProbPropagate(t,n,G); end=G->end[t][n];
	     q=G->matrix[t][end]; p+=log(q)+log(sum);
	   }
	}
	return (1.4427*p);
}

double	CondProbPropagate(Int4 t, Int4 n, gs_type G)
/*********************************************************************

       A          B             C

    [P(A1)][       0     ][      0      ]
    [P(A2)][       0     ][      0      ]
       :       :      :      :
       :       :      :      :
    [P(Ai)][P(Bi,A<i-end)][      0      ]
       :       :      :      :
       :       :      :      :
    [P(Ai)][P(Bi,A<i-end)][P(Ci,Bj,Ak)  ]

	P(A,B,C) = P(A) * P(B|A) * P(C|B).

   Note: P(M0,i) = 1.0 & THERE ARE ROUNDING ERRORS but this is dealt with.

   end[t][n] = end for sequence n and model t.

 *********************************************************************/
{
        register unsigned char    *seq;
        register Int4    i,end,len;
        register double  sum,*P = G->matrix[t];
        register double  *Pm1 = G->matrix[t-1];
	double		normL;
	fm_type		M=ModelCMSA(t,G->cmsa);

	/*** before convergence use seg'ed sequences ***/
        if(G->stage < 1) seq=XSeqPtrCMSA(n,G->cmsa);
	/*** but after convergence use 'raw' sequences ***/
	else seq = SeqPtrCMSA(n,G->cmsa);
	/*** set prob = 0 for i = 1..start-1 ****/
        for(i=0; i < G->strt[t][n]; i++) P[i] = 0.0;

	/*** compute P(Mi) ****/
        end = G->end[t][n];
  if(G->mod_temp && G->temp0 > 0){  // WARNING: NOT SURE THAT THIS IS CORRECT
  // if(G->mod_temp){
    normL=NormLikelihoodFModel(M);
    do {
        for(sum=0.0, i=G->strt[t][n]; i<= end; i++){
           P[i] = pow((LikelihoodFModel(seq,i,M)/normL),G->temp);
	   sum += P[i];		
        }
	if(sum >= DBL_MAX){ 
	    fprintf(stderr,
		"!!WARNING: overflow in CondProbPropagate( ) (norm=%g)!!\n", normL);
	    fprintf(stderr,"Sampling Temperature = %.2f K\n",G->temp0);
	    fprintf(stderr,"Adjusting likelihoods to compensate...\n");
		normL *= 100000.0;	/** increase normL **/
	} else if(sum == 0.0){
	    fprintf(stderr,
		"!!WARNING: underflow in CondProbPropagate( ) (norm=%g)!!\n", normL);
	    fprintf(stderr,"Sampling Temperature = %.2f K\n",G->temp0);
	    fprintf(stderr,"Adjusting likelihoods to compensate...\n");
		normL /= 100000.0;	/** decrease normL **/
	} else break;
    } while(TRUE);
  } else {
        for(sum=0.0, i=G->strt[t][n]; i<= end; i++){
           P[i] = (double)LikelihoodFModel(seq, i, M);
	   sum += P[i];		
        }
  }
	if(sum >= DBL_MAX) print_error("sum > DBL_MAX");
	// normalize probabilities.
        for(i=G->strt[t][n]; i<= end; i++) { P[i] /= sum; }

	// compute Markovian conditional probability.
	if(t > 1) len = LengthCMSA(t-1,G->cmsa); else len=0;

        for(i=G->strt[t][n]; i<= end; i++){ P[i] *= Pm1[i-len]; }
#if 0
/**** add gap penalty right here based on distance --v ****/
???????
        for(i=G->strt[t][n]; i<= end; i++){ P[i] *= Pm1[i-len]*gapfunct[i-len]; }
???????
#endif

/*** treat gaps as another residue (insertions or deletions) so
     that in addition to summing probabilities over all residue positions
     also need to sum over gap positions as well.  ************/

	// convert to cumulative probabilities.
        for(P[0] = 0.0,i=G->strt[t][n]; i<= end; i++){ 
		P[i] += P[i-1]; 
#if 0
	if(P[i] >= 1.0) fprintf(stderr,"?P[%d][%d][%d] = %g (%g + %g)\n",
			t,n,i,P[i],P[i]-P[i-1],P[i-1]);
#endif
	}
	return sum;
}

