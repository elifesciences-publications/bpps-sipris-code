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

#include "dsets.h"
#include "clique.h"
#include "probability.h"
#include "cmsa.h"
#include "chn_typ.h"
#include "sqd_typ.h"
#include "che_typ.h"

// #include "cha_typ.h"

#define	USAGE_START	"USAGE: chn_pps chn_fafile [options]\n\
   options:\n\
     -A<real>:<real> - alpha hyperparameters A0:B0 (default: A0=B0=1.0)\n\
     -B=<char>     - Background options:\n\
                      'A' = Use MainSet as background for all\n\
                      'B' = Use BPPS background partition\n\
                      'M' = merge foreground into background (default)\n\
     -best         - Output the optimum partition only\n\
     -c=<real>,<real> - fraction of hardest to classify seqs to ignore in each partition (default: 0.0)\n\
     -concise	   - Don't create extra files\n\
     -col=<int>:<int>  - Specify the min and max number fo columns allowed\n\
     -E<real>      - Fisher's exact test adjusted p-value cutoff (default: 0.001)\n\
     -F            - Forbid other residue positions with input pattern\n\
     -global       - Use global sequence weighting\n\
     -fixed        - Sample over pattern positions only (keep sequence partition fixed)\n\
     -I=<str>      - User inputed residue positions to be ignored (input string: \"3,5-7,9,11-17\").\n\
     -M<real>      - minimum frequence of conserved residues in foreground (default: 0.5)\n\
     			WARNING: Using this option can yield confusing results\n\
     				unless the corresponding '-c' option is used in chn_see\n\
     -m<int>       - minimum amount of information required at each pattern position (default: 5)\n\
     -mode=<char>  - set mode for initializing sequence partition (default: 'R')\n\
     			'A' = All sequences into foreground (for use with -fixed option)\n\
     			'S' = by score to query\n\
     			'R' = randomly\n\
     			'O' = by order in input file\n\
     -N=<int>      - Maximum number of significant pattern positions to highlight\n\
     		   - (This sets the contrast for the alignment.)\n\
     -P<set><int>,[<set><int>]  - User inputed residue sets (e.g., YF90,ST97,YF98).\n\
     -q<int>       - Minimum clique size (3-10; default 4)\n\
     -R<int>       - Random seed\n\
     -rho=<real>   - set prior probability (rho) that a column is a pattern position (default: 0.5)\n\
     			WARNING: no longer accepts ln(rho) as input in nats! (range: 0.0 < rho <= 0.5)\n\
     -Ri=<real>    - set prior probability that a sequence belongs to the foreground (default: 0.5)\n\
     			input is (Ri); e.g., Ri=0.001 --> Ri = 'penalty' of ln(0.001/0.999) nats per wt_seq\n\
     -sigma=<int>  - set percent stringency of BnmlKL divergence (default: 70% FG/30% BG)\n\
     			(Note: 50 <= sigma <= 100)\n\
     -sets=<char>  - Residue sets options (default 'G')\n\
     			'O' = allow only one residue at pattern positions\n\
     			'G' = \"gold standard\" residue sets\n\
     			'L' = large number of residue sets\n\
     			'R' = Root node residue sets\n\
     			'M' = middle node residue sets\n\
     -U<real>      - Clustering p-value cutoff (default: 0.0001)\n\
     -Z            - run zero run only (can't use with -P option).\n\
     -x            - Examine legal residue sets (corresponds to '-sets=' option).\n\n"

//  -c<real>      - minimum contribution of a seq to the MAP to be included in either partition (default: 0.0)\n

/**************************** Global Variables ******************************/
Int4	ChainBPPS(Int4 argc,char *argv[])
{ 
	Int4	arg;
	Int4	hpsz=200,MinClique=5;
	UInt4   seed=7061950;
	BooLean	UseGlobalSqWts=FALSE;
UseGlobalSqWts=TRUE;	// set by default...
	char	SFBG='B';
	BooLean	Forbid=FALSE,ZeroOnly=FALSE,concise=FALSE;
	double	MinKeyFrq=0.5,MaxGapFrq=0.5,pcut=0.0001,Fisher_pcut=0.001;
	Int4	min_nats=5;
	double	min_contrib=0.0;
	double	fract_ignored[20];
	char	*input_str=0,*forbid_str=0;
	Int4	start_pattern=0;
	double	A0=1.0,B0=1.0;
	char	Mode='S';
	Int4	RtnInt=0,Contrast=0;
	BooLean	fixed_partition=FALSE,best_only=FALSE,output_sst=FALSE;
	char    str[300],sets_mode='G';
	Int4	sigma=80;
	double	Sigma=0.80;
	// double	LnRho= 0.6931471805599452862;	// -log(0.5);
	double	**Rho=0;
	double	rho=0.5;
	double	PriorRi= 0.5;	// -log(0.5);
	Int4	min_num_col=3;
	Int4	max_num_col=30;

	// pcut=0.00001;

	time_t	time1=time(NULL); 
	if(argc < 2) print_error(USAGE_START);
	if(0){
         FILE *fp = open_file(argv[1],".pcmd","w");
         for(arg = 0; arg < argc; arg++) fprintf(fp,"%s ",argv[arg]); 
	 fprintf(fp,"\n"); fclose(fp);
	}
	fract_ignored[1]=0.0; fract_ignored[2]=0.0;
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
             case 'A': if(sscanf(argv[arg],"-A%lf:%lf",&A0,&B0)==2){
			if(A0 <= 0.0 || B0 <= 0.0) print_error(USAGE_START); 
			argv[arg][1] = ' '; 
		       } else print_error(USAGE_START); 
		break;
	     case 'B':
		if(sscanf(argv[arg],"-B=%c",&SFBG)==1){
		  if(!(SFBG=='M' || SFBG=='B' || SFBG=='A')) print_error(USAGE_START);
		} else print_error(USAGE_START); 
		  argv[arg][1] = ' '; 
		break;
	     case 'b': 
		if(strcmp("-best",argv[arg]) == 0){ 
			best_only=TRUE; argv[arg][1] = ' '; 
		} else if(argv[arg][2] == 0){ argv[arg][1] = ' '; }
		break;
	     case 'c': 
              if(sscanf(argv[arg],"-col=%d:%d",&min_num_col,&max_num_col)==2){
			if(min_num_col < 2 || min_num_col > max_num_col){
				fprintf(stderr,"Min(%d)/Max(%d) number columns out of range\n",
					min_num_col,max_num_col);
				print_error(USAGE_START);
			} argv[arg][1] = ' '; 
	      } else if(strcmp("-concise",argv[arg]) == 0){
                        concise=TRUE;
	      } else {
		// OLD method based on contribution in nats.
		// min_contrib=RealOption(argv[arg],'c',-1000.0,20.0,USAGE_START);
		// NEW method based on fraction of poorest seqs to ignore.
		// fract_ignored=RealOption(argv[arg],'c',0.0,0.5000001,USAGE_START);
		// Int4    ParseReals(char *str, double *values, const char *msg);
		// for getting separate values for each partitions
		if(argv[arg][2] != '=') print_error(USAGE_START);
		Int4 n = ParseReals(argv[arg] + 3,fract_ignored,USAGE_START);
		if(n != 2) print_error(USAGE_START);
		fract_ignored[2]=fract_ignored[1];
		fract_ignored[1]=fract_ignored[0];
		if(fract_ignored[2] < 0.0 || fract_ignored[2] > 0.9) print_error(USAGE_START);
		if(fract_ignored[1] < 0.0 || fract_ignored[1] > 0.9) print_error(USAGE_START);
                argv[arg][1] = ' '; 
	      } break;
	     case 'g': 
	      if(strcmp("-global",argv[arg]) == 0) UseGlobalSqWts=TRUE;
	      else print_error(USAGE_START);
		break;
	     case 'q': 
		MinClique=IntOption(argv[arg],'q',3,10,USAGE_START); 
                argv[arg][1] = ' '; break;
	     case 'E': Fisher_pcut=RealOption(argv[arg],'E',0.0,1.0,USAGE_START);
                argv[arg][1] = ' '; break;
	     case 'U': pcut=RealOption(argv[arg],'U',0.0,1.0,USAGE_START);
                argv[arg][1] = ' '; break;
	     case 'P':
		input_str=argv[arg]+2; start_pattern=1; argv[arg][1] = ' '; 
		break;
	     case 'I': if(argv[arg][2] == '='){ 
			forbid_str=argv[arg]+3; argv[arg][1] = ' '; 
		} break;
	     case 'm': 
		if(sscanf(argv[arg],"-mode=%c",&Mode)==1){
			if(Mode != 'A' && Mode != 'S' && Mode != 'R' && Mode != 'O'){
				print_error(USAGE_START);
			}
                  argv[arg][1] = ' '; 
		} else {
		  min_nats=IntOption(argv[arg],'m',0,1000,USAGE_START); 
                  argv[arg][1] = ' '; 
		}
		break;
	     case 'M': MinKeyFrq=RealOption(argv[arg],'M',0.0,1.0,USAGE_START);
                argv[arg][1] = ' '; break;
	     case 'N': 
             	if(sscanf(argv[arg],"-N=%d",&Contrast)==1){
			if(Contrast < 1) print_error(USAGE_START);
			argv[arg][1] = ' ';
                } else print_error(USAGE_START); 
		break;
             case 'R':
		if(sscanf(argv[arg],"-Ri=%lf",&PriorRi)==1){
		   if(PriorRi <= 0.0 || PriorRi >= 1.0) print_error(USAGE_START);
		    argv[arg][1] = ' ';
		} else if(sscanf(argv[arg],"-R%d",&seed)==1){
			 argv[arg][1] = ' ';
                } else print_error(USAGE_START); 
		break;
	     case 'r': 
	        if(sscanf(argv[arg],"-rho=%lf",&rho)==1){
		   // print_error("-rho option no longer valid");
		   if(rho <= 0.0 || rho > 0.5) print_error(USAGE_START);
		   argv[arg][1] = ' '; 
                } else print_error(USAGE_START); 
		break;
	     case 'F': Forbid=TRUE; argv[arg][1] = ' '; break;
	     case 'f': 
	      	if(strcmp("-fixed",argv[arg]) == 0){
			fixed_partition=TRUE; argv[arg][1] = ' '; 
		} else print_error(USAGE_START); 
		break;
	     case 's': 
	        if(sscanf(argv[arg],"-sigma=%d",&sigma)==1){
		   if(sigma < 50 || sigma > 100) print_error(USAGE_START);
		   Sigma = (double) sigma/100.0;
	        } else if(sscanf(argv[arg],"-sets=%c",&sets_mode)==1){
		   if(sets_mode != 'M' && sets_mode != 'G' && sets_mode != 'O' 
			&& sets_mode != 'R' && sets_mode != 'L'){
				print_error(USAGE_START);
		   }
		   // don't need to pass sets along to chn_typ, as it is not used within chn_pps
		} else print_error(USAGE_START);
		argv[arg][1] = ' '; break;
	     case 'Z': ZeroOnly=TRUE; 
		argv[arg][1] = ' '; break;
	     case 0: print_error(USAGE_START); break;
	     case 'x': output_sst=TRUE; argv[arg][1] = ' '; break;
	     default: ; // do nothing.
	   }
	}
	if(seed == 7061950) seed = (UInt4) time(NULL);
	sRandom(seed);
	if(input_str && concise) best_only=TRUE;
        chn_typ *chn = new chn_typ(argc,argv,200);
	cma_typ *IN_CMA=chn->GetIN_CMSA(); 	
	a_type	AB=AlphabetCMSA(IN_CMA[2]);
	sqd_typ	*sqd=0;

	char	*buffer=NULL;
#if 1	// Fix for '-Z' option 
     {
	char	tmp[20],r;
	Int4	res,s;
	Int4	site,os;
	if(ZeroOnly){	// Need to get input_string == query.
	  // WARNING: assumes that there are no flanking sequences in first seq!!!
	  e_type fakeE1=FakeSeqCMSA(1,IN_CMA[1]);
	  e_type trueE1=TrueSeqCMSA(1,IN_CMA[1]);
	  os=OffSetSeq(trueE1);
	  assert((LenSeq(fakeE1)+OffSetSeq(trueE1) + 5) < 100000);
	  NEW(buffer,LenSeq(fakeE1)*8 + 3,char);	// Handles seqs < 100,000 aa
	  for(s=1; s <= LengthCMSA(1,IN_CMA[1]); s++){
		if(IsDeletedCMSA(1, s,IN_CMA[1])) continue;
		res=ResSeq(s,fakeE1); r = AlphaChar(res,AB);
		if(r == 0 || r == 'X') continue;
	  	site = TruePosCMSA(1,s,IN_CMA[1]);
		site += os;
		sprintf(tmp,"%c%d,",r,site);
		strcat(buffer,tmp);
	  }
	  s = strlen(buffer);
	  assert(s > 2);
	  buffer[s-1]=0;
	  input_str=buffer;
	  start_pattern=1;
	  // exit(1);
	}
     }
	if(input_str)  fprintf(stderr,"seed pattern = %s\n",input_str);
#endif
#if 1	// Temporary fix for gaps ('-') in query sequence...need better fix later.
	{
	 e_type fakeE1=FakeSeqCMSA(1,IN_CMA[1]);
	 for(Int4 s=1; s<=LenSeq(fakeE1); s++){
	   Int4 r=ResSeq(s,fakeE1);
	   if(r==AlphaCode('X',AB)){ // set 'X' to 'A'in fake query seq.
	     fprintf(stderr,"WARNING: setting 'X' at position %d in query to 'A'.\n",s);
	     r=AlphaCode('A',AB); EqSeq(s,r,fakeE1);
	   }
	 }
	}
#endif
	sqd=new sqd_typ(IN_CMA[2],IN_CMA[1],MinKeyFrq,MaxGapFrq,sets_mode);
	Int4	NumClust;
	char **sst_str;
	if(input_str){
	  NEWP(sst_str,4,char);
	  sst_str[1]=input_str;
	  NumClust=1;
	} else {
	  if(!concise) fprintf(stderr,"Fisher_pcut: raw = %g",Fisher_pcut);
	  double adjust=sqd->SearchSpace();
	  Fisher_pcut = Fisher_pcut/adjust;
	  FILE *tabfp = 0;
#if 0
	  if(!concise){
	     fprintf(stderr,"; adjusted (%g) = %g\n",adjust,Fisher_pcut);
             tabfp = open_file(argv[1],".tables","w");
	     sqd->OutPutTables(tabfp); // will output nothing if tabfp == 0.
	  }
#endif
	  sst_str = sqd->CliqueStrings(MinClique,&NumClust,pcut,Fisher_pcut);
	  if(tabfp){ fclose(tabfp); }
	}

	// fprintf(stderr,"Returned significant seed patterns.\n");
	// sst_typ	**sstS=sqd->SortedResSets( );
	sst_typ	**sst=sqd->LegalResSets( );
	e_type Query=FakeSeqCMSA(1,IN_CMA[1]);
	if(output_sst){	// use to examine legal residue sets...
          for(Int4 j=1; j <= LenSeq(Query); j++){
                fprintf(stdout,"%d: ",j);
                for(Int4 s=1; sst[j][s]; s++){
                        fprintf(stdout,"{");
                        PutSST(stdout,sst[j][s],AB);
                        fprintf(stdout,"} ");
                } fprintf(stdout,"\n");
          } exit(1);
	}
#if 0	// categorical distribution for rho given # sets...use harmonic weighting...1/2 , 1/4, 1/8,...
        NEWP(Rho,LenSeq(Query) +3, double);
        for(Int4 j=1; j <= LenSeq(Query); j++){
	   // Count the number of sets with 1, 2, 3, ... etc. residues.
	   Int4 *NumSets; NEW(NumSets,nAlpha(AB) + 3, Int4);
           for(Int4 s=1; sst[j][s]; s++) NumSets[CardSST(sst[j][s],AB)]++;
	   NEW(Rho[j],nAlpha(AB) + 3, double);
	   double d=rho,total=0.0;
	   for(Int4 x=1; x <= nAlpha(AB); x++){
		double sum=0.0;
	 	for(Int4 y=1; y <= NumSets[x]; y++){
			d = d*rho; sum += d; 
		} total += sum;
		if(NumSets[x] > 0) Rho[j][x] = sum/(double)NumSets[x];	// take the average
	   } Rho[j][0] = 1.0 - total;
           fprintf(stdout,"%d: ",j);
	   total=0.0;
           for(Int4 s=1; sst[j][s]; s++){
		fprintf(stdout,"{");
  		PutSST(stdout,sst[j][s],AB);
		Int4 card=CardSST(sst[j][s],AB);
		fprintf(stdout,"}(%d(%d): %.3g) ",card,NumSets[card],Rho[j][card]);
		total+=Rho[j][card];
           } fprintf(stdout," total:%.3f\n",total);
	   free(NumSets);
	} // exit(1);
#else	// moved the above to alphabet.cc
	// uses a geometric distribution with p = rho input.  
	Rho=GetRhoCategoricalPriors(stdout, LenSeq(Query), rho, sst, AB);
#endif

	Int4	NumberFound=0;
	if(sst_str){
	  Int4		n,best_n=-1;
	  che_typ	**che=NULL,**che_sort=NULL;
	  double	map,*Map,bestmap=0.0,*map_sort=NULL;
	  NEWP(che,NumClust+4,che_typ);
	  NEW(Map,NumClust+4,double);
	  fprintf(stderr,"Running Bayesian partitioning with pattern selection:\n");
	  dh_type dH = dheap(NumClust+3,4);
#if 1	  // open seeds file
	  if(!concise){
            FILE *sfp = open_file(argv[1],".seeds","w");
            for(n=start_pattern; n <= NumClust; n++){
	  	if(sst_str[n] && !concise) fprintf(sfp,"pattern %d: '%s'\n",n,sst_str[n]);
	    }
	    fclose(sfp);
	  }
#endif
          for(n=start_pattern; n <= NumClust; n++){
	    // if(sst_str[n]) fprintf(stderr," seed pattern %d:\n\t%s\n",n+1,sst_str[n]);
	    fprintf(stderr," seed pattern %d.\n",n);
	    // if(n != 4) continue;	// Temporary for testing...
	    if(sst_str[n]){
	     if(sst_str[n] && !concise) fprintf(stderr,"pattern %d: '%s'\n",n,sst_str[n]);
	     che[n] = new che_typ(sst_str[n],chn,A0,B0,!concise,Mode,Rho,PriorRi,UseGlobalSqWts);
				// NOTE: !concise==verbose
	     if(fixed_partition) che[n]->FixedPartition( );
	     che[n]->SetMinNumColumns(min_num_col);
	     // che[n]->SetMaxNumColumns(max_num_col);
	     if(Contrast > 0) che[n]->SetContrast(Contrast);
	     if(input_str && Forbid){	// if pattern input...then don't use other columns
	       che[n]->Forbid(sst);
	     }
	     if(forbid_str && !Forbid){	// specify positions to ignore...
	       che[n]->Forbid(forbid_str,sst);
	     }
	     // partition into heirarchical sets using Gibbs sampling.
	     // Map[n]=map=che[n]->Sample(n,sst);
	     Map[n]=map=che[n]->Sample(TRUE,sst);
	     Map[n]=map=che[n]->Sample(FALSE,sst);

#if 0		// use -N= option not max # columns...
	     che[n]->SetMaxNumColumns(max_num_col);
	     che[n]->RmExcessColumns( );

	     Map[n]=map=che[n]->Sample(TRUE,sst);
	     Map[n]=map=che[n]->Sample(FALSE,sst);
#endif

	     fprintf(stderr,"%d columns: %d +columns; map = %.2f\n",che[n]->NumColumns( ),
			che[n]->NumSignificantColumns( ),map);
	     che[n]->PutSubLPRs(stderr);
	     if(NumClust == 1 && che[n]->NumSignificantColumns( ) < 3){
	  		print_error("Failed to find a significant divergent subgroup.");
	     } else if(che[n]->NumColumns( ) >= 3){ 	// >= 3 columns required.
	        if(map > bestmap){ bestmap=map; best_n=n; }
		insrtHeap(n+1,(keytyp)-map,dH);
	     }
	    }
	  }
	  // if(best_n < 0) print_error("Failed to find significant divergent subgroups.");
	  if(best_only){
	   if(best_n >= 0){
#if 1
            	bpps_typ *bpps=che[best_n]->BPPS();
		if(bpps){
            	      double alpha,A0,B0;
		      cma_typ mcma = che[best_n]->RtnMainSet();
            	      alpha = bpps->Parameters(A0,B0);
            	      SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,sets_mode,mcma);
		}
#endif
		if(concise) che[best_n]->PutChn("_pps.cha",SFBG,min_nats,MinKeyFrq); 
		else che[best_n]->PutChn("_pps.chn",SFBG,min_nats,MinKeyFrq); 
        	sprintf(str,"%s_pps.ptn",argv[1]);
		che[best_n]->PutSubLPRs(str);
		NumberFound++;
	   	RtnInt=-1;
	   } else RtnInt=0;
	  } else {
	   Int4 Ident=0,m;
	   BooLean	thesame;
           if(FALSE && ZeroOnly) {
#if 0		// Don't need this as treating ZeroOnly as an input pattern...
	     if(che[0] && best_n >= 0){
		if(Map[0] > 0){
		   NumberFound++;
		   if(concise){
        	     sprintf(str,"%s_pps.ptn",argv[1]);
		     che[0]->PutSubLPRs(str);
		     che[0]->PutChn("_pps.cha",SFBG,min_nats,MinKeyFrq); 
		   } else {
		     che[0]->ContribSeqMAP(0,SFBG,fract_ignored,min_nats,MinKeyFrq);
		     che[0]->PutAll(0,SFBG,min_nats);
		   }
		}
	     }
#endif
	   } else {	// sort by map...
	     if(!input_str){
	       NEW(map_sort,NumClust+4,double);
	       NEWP(che_sort,NumClust+4,che_typ);
	       for(m=start_pattern; !emptyHeap(dH); m++){
		  map = -minkeyHeap(dH); n=delminHeap(dH); n--;
	 	  che_sort[m]=che[n]; map_sort[m]=Map[n];
	       }
	       for(n=start_pattern; n <= NumClust; n++){
		  che[n]=che_sort[n]; Map[n]=map_sort[n];
	       } // if deleted set to null
	       free(che_sort); free(map_sort);
	     }
	     for(n=start_pattern; n <= NumClust; n++){
		if(che[n]){
		  thesame=FALSE;
		  for(m=start_pattern; m < n; m++){
		     if(che[m]){
			if(che[n]->TheSameSets(che[m])){ thesame=TRUE; break; }
			else if(che[n]->TheSamePattern(che[m])){ thesame=TRUE; break; }
			// else if(che[n]->TheSame(che[m])){ thesame=TRUE; break; }
		     }
		  } 
		  if(thesame) continue;	// nearly identical partition-pattern pairs found.
#if 0
		  if(Map[n] > 0 && n==0){
			// che[0]->ContribSeqMAP(0,SFBG,min_contrib,min_nats,MinKeyFrq);
			NumberFound++;
			che[0]->ContribSeqMAP(0,SFBG,fract_ignored, min_nats,MinKeyFrq);
			che[0]->PutAll(0,SFBG,min_nats);
		  } else 
#endif
		  // fprintf(stderr,"Map[%d] = %g\n",n,Map[n]);
		  if(Map[n] > 0){
			NumberFound++;
			// che[n]->ContribSeqMAP(Ident,SFBG,min_contrib,min_nats,MinKeyFrq);
			if(concise){
			  if(n == start_pattern){
        		   sprintf(str,"%s_pps.ptn",argv[1]);
			   che[n]->PutSubLPRs(str);
        		   sprintf(str,"_pps.cha");
		 	   che[n]->PutChn(str,SFBG,min_nats,MinKeyFrq); 
			  } else {
        		   sprintf(str,"%s%d_pps.ptn",argv[1],Ident);
			   che[n]->PutSubLPRs(str);
        		   sprintf(str,"%d_pps.cha",Ident);
		 	   che[n]->PutChn(str,SFBG,min_nats,MinKeyFrq); 
			  }
			} else {
			   che[n]->ContribSeqMAP(Ident,SFBG,fract_ignored,min_nats,MinKeyFrq);
#if 1
            		   bpps_typ *bpps=che[n]->BPPS();
			   if(bpps){
            		      double alpha,A0,B0;
			      cma_typ mcma = che[n]->RtnMainSet();
            		      alpha = bpps->Parameters(A0,B0);
            		      SetBPPS_CMA(alpha,(Int4)A0,(Int4)B0,sets_mode,mcma);
			   }
#endif
			   che[n]->PutAll(Ident,SFBG,min_nats,MinKeyFrq);
			   // che[0]->PutAll(0,SFBG,min_nats); // Did use no MinKeyFrq for n=0.
			} Ident++;
		  }
		}
	     }
	   } RtnInt=NumberFound;
	 }
	 Nildheap(dH);
         for(n=0; n <= NumClust; n++){ if(che[n]) delete che[n]; } free(che);
         for(Int4 j=1; j <= LenSeq(Query); j++) free(Rho[j]);
	 free(Rho); free(Map); 
	}
	if(sqd) delete sqd;
	if(buffer) free(buffer);
	delete chn;
	if(!concise){
	   if(NumberFound == 0) fprintf(stderr,"Failed to find a significant pattern-partition\n");
	   double runtime=difftime(time(NULL),time1);
	   fprintf(stderr,
		"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	}
	return RtnInt;
}


