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
#include "lha_typ.h"
#include "hierview.h"
#include "c2h_typ.h"

#define USAGE_START "  *** Bayesian Partitioning with Pattern Selection ***\n\
  Usage: bpps <step> <prefix> ... [options]\n\
  Options:\n\
     step = 1: Initial hierarchical partitioning of MSA into subgroups.\n\
        Usage: bpps 1 <prefix> [options]\n\
           type \"bpps 1\" for options\n\
     step = 2: Creation of hiMSA.\n\
        Usage: bpps 2 <prefix> [options]\n\
           type \"bpps 2\" for options\n\
           creates <prefix>_himsa.* files as input to BPPS 3\n\
     step = 3: Create a contrast alignment for <int>th node (& optional pymol script).\n\
        Usage: bpps 3 <prefix> <int> [options]\n\
           type \"bpps 3\" for options\n\
     step = A: Run steps 1-3 using default options with optional structural mapping.\n\
        Usage: bpps A <prefix> [pdb_paths]\n\
          [pdb_paths] is an optional file listing paths to pdb coordinate files\n\
            list one pdb path per line (e.g., /tmp/pdb_nat/4ag9_H.pdb)\n\
     step = E: Evaluate the consistency between hierarchies.\n\
        Usage: bpps E <prefix1> <prefix2> [options]\n\
          type \"bpps E\" for options\n\
          requires MSA (.mma) & checkpoint (.chk) files from BPPS 1 as input\n\
        Note: <prefix1>.mma file must be identical to <prefix2>.mma file.\n\
   Reference:\n\
     Neuwald, A.F., L. Aravind & S.F. Altschul. 2017. Inferring Joint Sequence-Structural\n\
        Determinants of Protein Functional Specificity. eLife In revision.\n\
\n"

static void	PrintError()
{ PrintLicenseStatement("bpps v1.0.8"); print_error(USAGE_START); }

int     main(int argc, char *argv[])
{
	int	i,a,A,Argc;
	char	**Argv,step; 
	if(argc < 2) PrintError();
	else step=argv[1][0];
	if(strchr("A123E",step)==NULL) PrintError();
	TurnOffLicenseStatement();
	NEWP(Argv,argc+5,char);
	Argv[0]=argv[0]; Argv[0]=argv[0]; Argc=argc-1;
	for(a=2,A=1; a < argc; a++,A++){ Argv[A]=argv[a]; }
	switch (step){
	   case 'A':
	     {	
		char str[20],pdb[100];
		if(argc < 3) print_error(USAGE_START);
		if(Argv[2] && Argv[2][0] != '-'){
		   sprintf(pdb,"-pdb=%s",Argv[2]); Argv[3]=pdb;
		} else Argv[3]=0;
		Int4 NumSets=0; Argc=2;
	        { 
		  omc_typ omc(1,Argc,Argv); Int4 rtn=omc.Run(); omc.PrintTime(stderr); 
		  hpt_typ *Hpt=omc.RtnHpt(); NumSets=Hpt->NumSets()-1; 
		} run_hieraln(2,Argc,Argv); 
		if(Argv[3]) Argc=4; else Argc=3;
		sprintf(str,"%s_himsa",Argv[1]); Argv[1]=AllocString(str);
		for(i=2; i <= NumSets; i++){
		    sprintf(str,"%d",i); Argv[2]=str; run_hierview(3,Argc,Argv); 
		} free(Argv[1]);
	     } break;
	   case '1':
	     { omc_typ	omc(1,Argc,Argv); Int4 rtn=omc.Run(); omc.PrintTime(stderr); } break;
	   case '2': { run_hieraln(2,Argc,Argv); } break;
	   case '3': { run_hierview(3,Argc,Argv); } break;
	   case 'E': { c2h_typ c2h(step,Argc,Argv); c2h.Run(); } break;
	   default: print_error(USAGE_START); break;
	} free(Argv); return 0;
}

