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

#include "mcs_typ.h"
#include "blosum62.h"
#include "rst_typ.h"
#include "swt_typ.h"

#define	MAIN_MCS_USAGE  "USAGE: mcsBPPS file_prefix [options]\n\
      (put input main alignment cma formated file in <file_prefix>.mma)\n\
      (put input seed alignment cma formated files in <file_prefix>.sma)\n\
      (put hyperpartition specification in <file_prefix>.hpt)\n\
   options:\n\
     -A=<real>:<real> - Set default alpha hyperparameters A0:B0\n\
     -Am=<real>:<real> - Set miscellaneous default alpha hyperparameters A0:B0\n\
     -col=<int>:<int>  - Default min and max number of columns\n\
     -col_start=<int> - Iteration on which column sampling starts (default: 2)\n\
     -convergence  - monitor convergence by outputting LPRs in <file_refix>.cnv\n\
     -del          - Don't treat deletions as random frequencies (default: treat as random)\n\
     -noseeds	   - ignore seed patterns in hyperpartition (*.hpt) file\n\
     -nocsq	   - Don't add a consensus sequence to the seed alignments\n\
     -fixed=<str>  - Keep all labeled input sequences within the set <str> of hpt.\n\
     -global       - Use global sequence weighting(not yet implemented)\n\
     -iter=<int>:<int>:<int>  - set IterStart,IterEvolve,IterRounds (default: 1:3:4).\n\
     -maxcol=<int>  - Reset default max number of columns\n\
     -minnats=<real>  - Reset default minimum nats per weighted sequence (default: 2.0)\n\
     -N=<int>      - Set default contrast setting (N = number of columns to highlight)\n\
     -numcol=<int> - Initial number of columns for seed patterns(default: 25)\n\
     -put          - output intermediate files for hyperpartition categories and sets\n\
     -ppb=<int>    - parts per billion increase in LPR required to keep going (default: 100)\n\
                             Larger values lead to shorter (less optimized) runs (range: 1-1000000)\n\
     -Pttrn=<int>  - Specify length of automatically generated seed patterns\n\
     -R=<int>      - Random seed\n\
     -RandomInit   - Randomly partition the sequences to initialize\n\
     -Ri=<real>    - Set global prior probability that a sequence belongs to the foreground\n\
     -Rim=<real>   - Set misc global prior probability that a sequence belongs to the foreground\n\
     -rho=<real>   - Set global rho parameter\n\
     -rtf          - output individual rtf files, one for each hyperpartition set\n\
     -sets	   - Save the sets for comparing hierarchies (evalBPPS program)\n\
     -strict       - Require strict pattern independence between categories\n\
     -syntax       - show hyperpartition syntax\n\
     -tree	   - The input FD-table is in tree format\n\
   \n"

void	mcs_typ::PrintUsage(FILE *fp)
{
	fprintf(stderr,"%s\n",MAIN_MCS_USAGE);
	print_error(HPT_PUBLIC_USAGE);
}

void	mcs_typ::ReadMainArg(Int4 argc, char *argv[])
{
	Int4 x,arg;
	UInt4   seed=7061950;
	double	d;
	// StrictIndepend=FALSE; SaveSets=FALSE; PrintEachRTF=FALSE; NoCSQ=FALSE; cfp=ifp=0;
	if(argc < 2){
		print_error(MAIN_MCS_USAGE); 
	}
	program_name=AllocString(argv[0]);
	infile=AllocString(argv[1]);
	if(argc == 2 && argv[1][0] == '-'){
	      // if(strcmp("-syntax",argv[1]) == 0) print_error(HPT_PUBLIC_USAGE);
	      fprintf(stderr,"%s\n",MAIN_MCS_USAGE);
	      print_error(HPT_PUBLIC_USAGE);
	}
	for(arg = 2; arg < argc; arg++){
// fprintf(stderr,"DEBUG 0: argc=%d\n",argc); assert(argc != 3);
	   if(argv[arg][0] != '-') print_error(MAIN_MCS_USAGE);
	   switch(argv[arg][1]) {
             case 'A':
		if(sscanf(argv[arg],"-Am=%lf:%lf",&MiscGlobalA0,&MiscGlobalB0)==2){
                        if(MiscGlobalA0 <= 0.0 || MiscGlobalB0 <= 0.0) print_error(MAIN_MCS_USAGE);
		} else if(sscanf(argv[arg],"-A=%lf:%lf",&GlobalA0,&GlobalB0)==2){
                        if(GlobalA0 <= 0.0 || GlobalB0 <= 0.0) print_error(MAIN_MCS_USAGE);
                } else print_error(MAIN_MCS_USAGE);
                break;
	     case 'c': 
              if(sscanf(argv[arg],"-col_start=%d",&x)==1){
		if(x < 1) print_error(MCS_USAGE_START);
		else ColSampleStart_DF=x;
              } else if(sscanf(argv[arg],"-col=%d:%d",&DefaultMinCol,&DefaultMaxCol)==2){
			if(DefaultMinCol < 2 || DefaultMinCol > DefaultMaxCol){
				fprintf(stderr,"Default Min(%d)/Max(%d) # columns out of range\n",
					DefaultMinCol,DefaultMaxCol);
				print_error(MCS_USAGE_START);
			} argv[arg][1] = ' '; 
	      } else if(strcmp("-convergence",argv[arg]) == 0){
		cfp=open_file(infile,".cnv","w");
		ifp=open_file(infile,".itr","w");	// simulated annealing temperature
              } else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'd': 
	      if(strcmp("-del",argv[arg]) == 0){ del_as_random=0; } else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'i': 
		{
		  Int4	s,e,r;
                  if(sscanf(argv[arg],"-iter=%d:%d:%d",&s,&e,&r)==3){
			if(s < 1) print_error(MAIN_MCS_USAGE);
			if(e < 1) print_error(MAIN_MCS_USAGE);
			if(r < 2) print_error(MAIN_MCS_USAGE);
			IterStart_DF=s; IterEvolve_DF=e; NumRounds_DF=r;
		  } else print_error(MAIN_MCS_USAGE);
		}
	      break;
	     case 'm': 
               if(sscanf(argv[arg],"-minnats=%lf",&d)==1){
			// fprintf(stderr,"%s: d=%.3f\n",argv[arg],d);
			if(d <= 0.0 || d > 100) print_error(MCS_USAGE_START);
			else { MinimumNatsPerWtSeq=d; argv[arg][1] = ' ';  }
               } else if(sscanf(argv[arg],"-maxcol=%d",&DefaultMaxCol)==1){
		  if(DefaultMinCol > DefaultMaxCol){
			fprintf(stderr,"Default Min(%d)/Max(%d) # columns out of range\n",
					DefaultMinCol,DefaultMaxCol);
			print_error(MCS_USAGE_START);
		  } argv[arg][1] = ' '; 
               } else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'N': 
                if(sscanf(argv[arg],"-N=%d",&GlobalN)==1){
			if(GlobalN < 2) print_error(MAIN_MCS_USAGE);
		} else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'n': 
              if(sscanf(argv[arg],"-numcol=%d",&x)==1){
		if(x < 1) print_error(MCS_USAGE_START);
		else SeedPttrnLen=x;
	      } else if(strcmp("-noseeds",argv[arg]) == 0) NoSeeds=TRUE;
	      else if(strcmp("-nocsq",argv[arg]) == 0) NoCSQ=TRUE;
	      else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'P': 
              if(sscanf(argv[arg],"-Pttrn=%d",&SeedPttrnLen)==1){
		if(SeedPttrnLen < 5) print_error("SeedPttrnLen must be >= 5");
              } else print_error(MAIN_MCS_USAGE);
	      break;
	     case 'p': 
	      if(sscanf(argv[arg],"-ppb=%u",&ppb_increase) == 1){
		if(ppb_increase < 1 || ppb_increase > 1000000) print_error(MAIN_MCS_USAGE);
	      } else if(strcmp("-put",argv[arg]) == 0) PutIntermediateFiles=TRUE;
	      else print_error(MAIN_MCS_USAGE);
	      break;
             case 'R':
	        if(strcmp("-RandomInit",argv[arg]) == 0) PartitionRandomly=TRUE;
		else if(sscanf(argv[arg],"-Rim=%lf",&MiscGlobalRi)==1){
                   if(MiscGlobalRi <= 0.0 || MiscGlobalRi >= 1.0) print_error(MAIN_MCS_USAGE);
		} else if(sscanf(argv[arg],"-Ri=%lf",&GlobalRi)==1){
                   if(GlobalRi <= 0.0 || GlobalRi >= 1.0) print_error(MAIN_MCS_USAGE);
                } else if(sscanf(argv[arg],"-R=%d",&seed)!=1) print_error(MAIN_MCS_USAGE);
                break;
             case 'r':
		if(sscanf(argv[arg],"-rho=%lf",&Global_rho)==1){
                   if(Global_rho <= 0.0 || Global_rho >= 0.5) print_error(MAIN_MCS_USAGE);
		} else if(strcmp("-rtf",argv[arg]) == 0) PrintEachRTF=TRUE;
	        else print_error(MAIN_MCS_USAGE);
                break;
             case 's':
	      if(strcmp("-sets",argv[arg]) == 0) SaveSets=TRUE;
	      else if(strcmp("-strict",argv[arg]) == 0) StrictIndepend=TRUE;
	      else if(strcmp("-syntax",argv[arg]) == 0) print_error(HPT_PUBLIC_USAGE);
	      else print_error(MAIN_MCS_USAGE);
                break;
	     case 't':
              if(strcmp("-tree",argv[arg]) == 0){  
		  // then allow elimination of inappropriate subcategories...as for pmcBPPS.
                  IsTreeHpt=TRUE;
                } else print_error(MAIN_MCS_USAGE);
                break;

             case ' ': break;	// ignore these...
	     default: print_error(MAIN_MCS_USAGE);
	      break;
	  }
	}
#if 0	// seeded from calling routine...
	if(seed == 7061950) seed = (UInt4) time(NULL);
	if(seed != 0) sRandom(seed);
	// if(strcmp(program_name,"cdhBPPS")==0) IsTreeHpt=TRUE;
	fprintf(stderr,"random seed = %u\n",seed);
#endif
	// fprintf(stderr,"Program Name = \"%s\"\n",program_name); exit(1);
}


