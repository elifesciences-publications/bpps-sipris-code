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

#include "asset.h"

ast_typ	CreateAsset(char *DBS_NAME, a_type A)
{
	ast_typ	D;

	NEW(D,1,asset_type); D->infile = AllocString(DBS_NAME); D->A=A; 
	D->P = NULL; D->tfreq = NULL; D->H=NULL; D->pheap=NULL; 
	D->Q = NULL; D->eblocks = NULL; D->E = NULL;
	D->query = FALSE; D->shuffle = FALSE; D->xnu = TRUE;
	D->d_min=MINDEPTH; D->d_max=6; D->maxpat=20; D->verbose = FALSE; 
	D->k_max = 15; D->sd=2.0; D->ncalc = 0; D->percent = 0;
	D->hpsz=0; D->c_min=4; D->n_min = 3; 
	D->pout = 0.5; D->min_prob=(-5.0);
	D->minscore= 0; D->s_min = 4; D->mode = 26;
	return D;
}

BooLean	SetAsset(char *command, ast_typ D)
/* set parameters for depth-first-search */
{
	Int4	i;
	float	f;

	if(D->P != NULL) 
		asset_error("set parameters prior to executing operations");;
        if(command[0] != '-') return FALSE;
        switch(command[1]) {
              case 'c':
         	if(sscanf(command,"-c%d", &i) == 1 )
			{ D->c_min=MAXIMUM(Int4,1,i); return TRUE; }
              break;
              case 'D':
		if(sscanf(command,"-D%d",&i)==1){
        	   if(i < MINDEPTH) {
			fprintf(stderr,"d_max must be > %d\n",MINDEPTH-1);
		   } else { D->d_max=i; return TRUE; }
		}
              break;
              case 'd':
		if(sscanf(command,"-d%d",&i)==1){
        	   if(i < MINDEPTH) {
			fprintf(stderr,"d_min must be > %d\n",MINDEPTH-1);
		   } else { D->d_min=i; return TRUE; }
		}
              break;
              case 'E':
                if(sscanf(command,"-E%f", &f)== 1)
			{D->sd=MINIMUM(double,f,20.0); return TRUE; }
              break;
              case 'e':
                if(sscanf(command,"-e%f", &f)== 1)
			{D->pout=MINIMUM(double,f,5000000.0); return TRUE; }
              break;
              case 'f': 
                if(sscanf(command,"-f%d", &i)== 1){
		  if(i < 0 || i > 100){
			fprintf(stderr,"scan file percent must be 0-100");
		  } else { D->percent=i; return TRUE; }
		}
	      break;
              case 'H':
		   D->H=Histogram("log10(E-values)",-10,20,1.0); return TRUE;
              case 'h':    
                if(sscanf(command,"-h%d", &i) == 1 )
		   { D->hpsz=MINIMUM(Int4,100000,MAXIMUM(Int4,1,i)); return TRUE;}
              break;
              case 'k':
                 if(sscanf(command,"-k%d", &i)== 1) {
        		if(i <= MAXSEGMENT && i >= MINSEGMENT) { 
				D->k_max=i; return TRUE; 
			} else {
            		    fprintf(stderr,"k_max must be > %d and < %d\n",
		 		  MINSEGMENT - 1, MAXSEGMENT +1);
			}
		  }
              break;
              case 'l': D->xnu=FALSE; return TRUE; 
              case 'n':
         	if(sscanf(command,"-n%d", &i)== 1)
			{D->n_min=MAXIMUM(Int4,1,i);return TRUE;}
              break;
              case 'P':
                if(sscanf(command,"-P%d", &i)== 1)
        			{ D->mode=MAXIMUM(Int4,0,i); return TRUE; }
              break;
              case 'O':    
                if(sscanf(command,"-O%d", &i) == 1 )
				{ D->minscore=MAXIMUM(Int4,0,i); return TRUE; }
              break;
              case 'o':    
                if(sscanf(command,"-o%d", &i) == 1 )
				{ D->maxpat=MAXIMUM(Int4,1,i); return TRUE; }
              break;
              case 'q': D->query=TRUE; return TRUE; 
              case 'r': D->shuffle=TRUE; return TRUE;
              case 'S': 
		if(sscanf(command,"-S%d",&i)==1)
				{ D->seed=i; sRandom(i); return TRUE; }
		break;
              case 's': 
		if(sscanf(command,"-s%d",&i)==1)
				{ D->s_min = MAXIMUM(Int4,0,i); return TRUE; }
		break;
              case 'v': D->verbose=TRUE; return TRUE;
              case 'x':
               	if(sscanf(command,"-x%f", &f) == 1)
		   { D->min_prob=MAXIMUM(double,-50.0,f); return TRUE; } 
              break;
              default: return FALSE;
        }
	return FALSE;
}

ast_typ	InitAsset(ast_typ D)
{
	Int4	i;
	a_type	A=D->A;

	if(D->P != NULL) asset_error("initialization already done");
        if(D->d_max > D->k_max) asset_error("d_max must be <= k_max");
        if(D->s_min > D->d_max) D->s_min=D->d_max;
        if(D->d_min > D->d_max) D->d_min=D->d_max;
        if(D->s_min < 2) asset_error("s_min must be > 2");
	D->mode = MINIMUM(Int4,D->mode,NPairsAlpha(D->A));
	if(D->shuffle){
		D->P = SeqSet(D->infile,D->A);
		ShuffleSeqSet(D->P);
        } else if(D->xnu) {
           if((D->P=MkXSeqSet(D->infile,D->A)) == NULL){
                asset_error("can not remove low entropy regions");
           }
        } else D->P = SeqSet(D->infile,D->A);
	D->tfreq = tFreqSeqSet(D->P);
	D->Q = Pattern(D->k_max,nAlpha(D->A)); MultiVar2Pattern(D->d_max,D->Q);
	if(D->hpsz<=0) {
		switch(D->d_max) {
                case 4: D->hpsz = 200; break;
                case 5: D->hpsz = 500; break;
                case 6: D->hpsz = 1000; break;
                default: D->hpsz = 2000;
                }
	}
	/*** alloc_asset_blocks ***/
	if(D->n_min > NSeqsSeqSet(D->P)) D->n_min = NSeqsSeqSet(D->P);
	if(D->n_min > D->c_min) D->n_min = D->c_min;
	D->d_max = MINIMUM(Int4, D->k_max,iPattern(D->Q)); 
	D->ncalc=0;
	D->eblocks = MkEBlocks(D->P, D->k_max);
	D->IB=BlockL(BlockN(EBlocksBu(D->eblocks)));
	D->UB=Block(BlockN(EBlocksBu(D->eblocks)));
	ClearBlockL(D->IB); ClearBlock(D->UB);
	NEWP(D->motif,nAlpha(A)+1,char);
	for(i=0;i<=nAlpha(A);i++) NEW(D->motif[i],D->k_max+1,char);
	if(D->s_min < D->d_max && D->mode != 0){
		AddEBlocks2(D->eblocks, D->mode);
        	NEWP(D->Card,D->d_max+1,Int4);
        	for(i=1; i<=D->d_max; i++) { 
			NEW(D->Card[i],nAlpha(D->A)+2,Int4); 
		}
	}
	if(D->s_min > D->d_max) D->s_min = D->d_max;
	MEW(D->B,D->d_max+1,b_type);
	MEW(D->EOB,D->d_max+1,b_type);
	MEW(D->cnts,D->d_max+1,Int4);
	for(i=0;i<=D->d_max;i++) {
		D->cnts[i] = 0;
		D->B[i] =BlockL(BlockN(EBlocksBu(D->eblocks)));
		D->EOB[i] =BlockL(BlockN(EBlocksBu(D->eblocks)));
	}
	D->E=MkEvalue(D->d_max, D->d_min,D->k_max, D->min_prob, log10(D->pout),
		D->s_min , D->mode,EBlocksN(D->eblocks),
			tFreqSeqSet(D->P),D->c_min,D->A);
	/*** end alloc_asset_blocks ***/
	return D;
}

void	PutAssetIDs(FILE *fptr,ast_typ D)
{ if(D->P == NULL) InitAsset(D); PutSeqSetPIDs(fptr,D->P); }

void	PutInfoAsset(FILE *fptr, ast_typ D)
{
	Int4	i;
	double	sum,n;

	if(D->P == NULL) InitAsset(D);
	fprintf(fptr,"\n  input file:\n");
	fprintf(fptr,"\tname: \"%s\"\n\ttotal sequences: %d",
			NameSeqSet(D->P),NSeqsSeqSet(D->P));
	fprintf(fptr,"\n\tsequence lengths: %d-%d residues\n",
		       MinSeqSeqSet(D->P),MaxSeqSeqSet(D->P));
	for(sum=0.0,i=1;i<=NSeqsSeqSet(D->P); i++){
		sum += (double) SqLenSeqSet(i,D->P);
	}
	n=(double) NSeqsSeqSet(D->P);
	fprintf(fptr,"\taverage length = %d\n",
			(Int4)((sum/n)+0.5));
	fprintf(fptr,"\t%d segments of length %d\n\n",
			EBlocksN_k(D->k_max,D->eblocks),D->k_max);
	fprintf(fptr,"  search parameters: \n");
	fprintf(fptr,"\td_max: %d\n",D->d_max);
	fprintf(fptr,"\td_min: %d\n",D->d_min);
	fprintf(fptr,"\ts_min: %d\n", D->s_min);
	fprintf(fptr,"\tc_min: %d\n",D->c_min);
	fprintf(fptr,"\theap size: %d\n",D->hpsz);
	if(D->s_min < D->d_max){
		fprintf(fptr,"\tnumber of pairs: %d\n",D->mode);
	}
	fprintf(fptr,"\tminimum matching sequences: %d\n",D->n_min);
	fprintf(fptr,"\tsearching %g patterns\n", KEvalue(D->E));
	fprintf(fptr,"\tstore in heap if pattern probabitity < 1e%.2f\n",
		pmax0Evalue(D->E));
	fprintf(fptr,"\tprint if E-value < %g (1e%.1f)\n", 
		D->pout,log10(D->pout));
}

/*********************** assets output **************************/
void	PutAssetHeap(FILE *fptr,ast_typ D)
{
  Int4		numG,count,os,npat,g,i,j,k;
  Int4		N,n,m;
  double	p,K,pout,hg_cut,rcut;
  s_type	**list;
  ptn_typ	*motif,*pat,Q;
  char		c;
  aln_typ	*align;
  FILE		*snfptr=NULL;
  char		str[300];

  rcut = pow(10.0,(-(double)D->minscore/100.0)); 
  K=KEvalue(D->E);
  // PutPheap(fptr, D->pheap); // DEBUG...
#if 0
N=EBlocksN_k(D->k_max,D->eblocks); /** OLD **/
#endif
  N=NSegsEBlocks(D->eblocks); 
  pout = pmax0Evalue(D->E);
  while(pout<MaxKeyPheap(D->pheap)&&(Q=DelMaxPheap(NULL,D->pheap))!=NULL)
	   		{ NilPattern(Q); }
  i = ItemsInPheap(D->pheap);
  fprintf(fptr,"Saved %d out of %g patterns\n",i,K);
  hg_cut = i*i*D->k_max; hg_cut = 0.05/hg_cut;
  hg_cut = MAXIMUM(double,hg_cut,0.000001);
  if(!EmptyPheap(D->pheap)) {
	Q = D->Q;
#if 0
npurged = PurgePheap(D->pheap);
#endif
	pat = CombinePheap(&motif, &list, &numG, &npat, N, D->pheap);
	fprintf(stderr,"Finding segments with motifs");
	align = MkAlignsCombine(&numG, motif, list, hg_cut, D->P, rcut,N);
	free(list);
	fprintf(stderr,"\n");
	for(m=g=0; g<numG; g++){
	  if(!EmptyAlign(align[g])){
	   	D->Q = motif[g]; k = kPattern(D->Q);
		/*** PRINT PATTERNS FOR GROUP ***/
	        fprintf(fptr,"\n                ******* MOTIF ");
	        if(numG<26)fprintf(fptr,"%c *******\n\n",m+'A');
		else fprintf(fptr,"%d *******\n\n", m+1);
		m++;
		fprintf(fptr,"%*s %14s\n%-*s %7s %4s%5s  %7s  %7s\n",
                	k+19,"", "__log(E-value)__",k,
                	"PATTERN","EXP","OBS","SEQ","  Mult.","Sampled");
		for(count=j=0; j<npat; j++) { 
		   if(CombinablePatterns(motif[g],pat[j],&os)){ 
		      if(count < D->maxpat){
			D->Q=pat[j]; PutAsset(fptr,D,os,k);
		      } count++;
		   }
		}
		D->Q = motif[g];
        	for(j=0;j<k;j++){
			if(j%5==0)c='+';else c='-';fprintf(fptr,"%c",c);
		} fprintf(fptr,"\n");
        	for(j=0;j<k;j++){if(j%5==0)fprintf(fptr,"%-5d",j);}
		/*** PRINT SEGMENTS FOR GROUP ***/
		if(count > D->maxpat){
		    fprintf(fptr,
			"(%d related patterns not shown)\n",
				count - D->maxpat); 
		} else fprintf(fptr,"\n");
		n = PutAlign(fptr, align[g], D->verbose, D->Q);
		if(D->percent > 0){
		  p = ((double)n*100.0)/(double) NSeqsSeqSet(D->P);
		  if(p >= (double)D->percent){
		    if(!D->query || SeqInAlign(1, align[g])){
	   	      if(snfptr== NULL){
	   		strcpy(str,NameSeqSet(D->P));
           		strcat(str,".sn");
		  	if((snfptr = fopen(str,"w")) == NULL){
      		  	   asset_error("fatal - could not open scan file");
			}
          	      }
		      ScanfileAlign(align[g],D->verbose,snfptr);
		    }
		  }
		}
	   }
	   D->Q = NULL; 
	}
	if(snfptr != NULL) fclose(snfptr);
	fprintf(fptr,"\n\t%d groups.\n", m);
	for(g=0; g<numG; g++){
		NilAlign(align[g]); NilPattern(motif[g]); motif[g] = NULL;
	}
        for(i=0;i<npat;i++) { NilPattern(pat[i]); } free(pat);
	free(motif); free(align);
	D->Q = Q;
	fprintf(fptr,"\t%d hits\n\n", npat); /** new **/
  } else fprintf(fptr,"\tno hits\n\n");
  if(D->H != NULL) PutHist(fptr,60,D->H);
  NilPheap(D->pheap); D->pheap=NULL;
  fprintf(fptr,"\t%d probability calculations (>%-1.2f std.dev.).\n",
                D->ncalc,D->sd );
  for(i=1;i<=D->d_max;i++) {
              fprintf(fptr,"\tdepth %d: %d patterns\n",i,D->cnts[i]);
  }
}

BooLean	PutAsset(FILE *fptr,ast_typ D, Int4 offset, Int4 length)
{
	Int4	C,c,j,i,tmax,N,add;
	double	p,prob,E;

    if(D->Q == NULL || offset < 0) {
	fprintf(fptr,"\tillegal\n"); return FALSE; 
    }
    add = length - LengthPattern(D->Q) - offset;
    C=EvaluateAsset(&p,&i,D);
    N=EBlocksN_k(i+1,D->eblocks);
    if((prob=Log10CBP(C,N,p))==ILLEGAL) return FALSE;
    E=(double)N*p;
    tmax = GetPattern(D->motif,D->Q,D->A);
    for(c=0;c<tmax;c++) {
	if(c==0) for(j=0; j< offset; j++) fprintf(fptr,".");
	else for(j=0; j< offset; j++) fprintf(fptr," ");
	fprintf(fptr,"%s",D->motif[c]);
	if(c==0 && prob <= 1.0) {
	    for(j=0; j < add; j++) fprintf(fptr,".");
	    if(E>2.0) fprintf(fptr," %7.1f %4d %4d ",
		E,C,NumEBlocksP(D->IB,D->eblocks));
	    else fprintf(fptr," %7.1g %4d %4d ",
		E,C, NumEBlocksP(D->IB,D->eblocks));
	   if(D->E != NULL){
		fprintf(fptr,"%+7.1f",MultEvalue(D->E,prob)); 
		if(MultEvalue(D->E,prob) > D->min_prob) {
			fprintf(fptr,"  %+7.1f\n",
				FastEvalue(D->E, prob));
		} else fprintf(fptr,"      nd\n");
	   } else fprintf(fptr,"%7s\n","  -  "); 
	} else fprintf(fptr,"\n");
    }
    return TRUE;
}

void	NilAsset(ast_typ D)
{
	ss_type	data;
	ph_type PH = NilAssetRtnPHeap(&data,D);
	if(PH != NULL) NilPheap(PH);
	if(data !=NULL) NilSeqSet(data); 
}

ph_type	NilAssetRtnPHeap(ss_type *data, ast_typ D)
{
	Int4	c,i;

	if(D!=NULL) {
		ph_type PH = D->pheap;
		*data = D->P;
		/*** nilassetblocks ***/
		for(i=0;i<=D->d_max;i++) {
			NilBlock(D->B[i]); NilBlock(D->EOB[i]);
		}
		if(D->IB != NULL) NilBlock(D->IB);
		if(D->UB != NULL) NilBlock(D->UB);
		free(D->B); free(D->EOB); free(D->cnts); 
		NilEBlocks(D->eblocks);
		if(D->s_min < D->d_max && D->mode != 0){
        		for(i=1; i<=D->d_max; i++) free(D->Card[i]);
        		free(D->Card); D->Card = NULL;
		}
		NilEvalue(D->E); D->E = NULL;
		/*** end nilassetblocks ***/
		free(D->infile);
		if(D->motif != NULL) {
		    for(c=0;c<=nAlpha(D->A);c++) free(D->motif[c]);
		    free(D->motif);
		}
		if(D->H != NULL) NilHist(D->H);
		if(D->Q != NULL) NilPattern(D->Q);
		free(D);
		return PH;
	} return 0;
}

void	asset_error(const char *s) { fprintf(stderr,"asset: %s\n",s); exit(1); }

