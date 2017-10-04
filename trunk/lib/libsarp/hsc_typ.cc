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

#include "hsc_typ.h"

//*************************** hsc_typ ****************************************
//*************************** hsc_typ ****************************************
//*************************** hsc_typ ****************************************

#define USAGE_START \
"USAGE: sarp <pdb_name_file> <mcBPPS_prefix> [-options]\n\
    (files needed: <mcBPPS_prefix>_new.hpt <mcBPPS_prefix>.pttrns <mcBPPS_prefix>_new.mma\n\
    -color=<str>    - Designate repeat colors to use (e.g., <str>=\"MROYGC\")\n\
    -Color=<str>    - Designate sidechain colors to use (e.g., <str>=\"MROYGC\")\n\
    -dummy \n\
\n\n"

#if 0	// input files:
        SDR_4CB.pttrns  // all patterns (FILE *fp1)
        SDR_4CB.hpt     // hpt          call as "SDR_4CB", assumes "SDR_4CB.hpt"
        SDR_4CB_new.mma // subgroup alignments; open all of these...  (FILE *fp2)
           // pull out cma file for vsi.
        pdb_infile      // pdb file paths
#endif

void	hsc_typ::Init(FILE *mmafp, FILE *hptfp,FILE *ptrnfp)
// mmafp and hptfp are passed in by hierview to avoid too many input files!!
{
	// read in hyperpartition...

	// read in subgroup alignments...
	FILE *fp=0;
	AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62); verbose=TRUE;
	if(mmafp){ IN_CMA=MultiReadCMSA(mmafp,&Number,AB); fclose(mmafp); }
        else {
	   sprintf(string,"%s_new",mcBPPS_prefix); CheckString();
	   fp = open_file(string,".mma","r");
	   IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	} assert(IN_CMA != 0);
#if 1
	if(hptfp){ hpt=new hpt_typ(hptfp); fclose(hptfp); } // for hierview... 
	else {
	  sprintf(string,"%s_new.hpt",mcBPPS_prefix); CheckString();
	  if((fp=fopen(string,"r")) == 0){
	     sprintf(string,"%s_chk",mcBPPS_prefix); CheckString();
	     fp = open_file(string,".hpt","r"); 
	  } hpt=new hpt_typ(fp); fclose(fp);
	}
#else
	fp = open_file(string,".hpt","r"); hpt = new hpt_typ(fp); fclose(fp);
#endif
#if 1	// lengthen IN_CMA array 
	if(hpt->NumSets() > Number) {	// then lengthen the IN_CMA array using Null's
	   Int4	x,y;
	   cma_typ *in_cma; NEW(in_cma,hpt->NumSets() +3, cma_typ);
	   for(x=1,y=1; x <= hpt->NumSets(); x++)
	   // for(x=1,y=1; x < hpt->NumSets(); x++)
	   {
#if 1
	        if(IN_CMA[y] == 0 && x == hpt->NumSets()){
	     	  fp = tmpfile(); 
	          PutRandomCMA(fp,tFreqCMSA(IN_CMA[1]),LengthCMSA(1,IN_CMA[1]),1,AB);
	          rewind(fp); in_cma[x]=IN_CMA[y]=ReadCMSA(fp,AB); fclose(fp); Number++;
	          RenameCMSA("Random",IN_CMA[y]); y++; continue;
		}
#endif
		if(strcmp(NameCMSA(IN_CMA[y]),hpt->ElmntSetName(x)) == 0){ // --> same set.
			in_cma[x]=IN_CMA[y];  y++;
		} else if(y == 1 || hpt->TypeOfSet(x) != '?'){
			print_error("*.hpt and *.mma files are fatally inconsistent");
		} 
	   } y--; 
	   if(y != Number) print_error("*.hpt and *.mma files are fatally inconsistent");
	   Number = hpt->NumSets(); free(IN_CMA); IN_CMA=in_cma; 
	}
#endif
	// read in patterns ...
	if(ptrnfp){ ptrn = new pat_typ(ptrnfp); fclose(ptrnfp); } // for hierview
	else { fp = open_file(mcBPPS_prefix,".pttrns","r"); ptrn = new pat_typ(fp); fclose(fp); }

	mpdb = new mps_typ(pdb_paths);
	if(mpdb->NumPDB <= 0) print_error("hsc_typ::Init( ) fatal error: no pdb files created!");

	Validate();
}

void	hsc_typ::Free()
{
	Int4 i;
	delete p2c;	// need to free this before the rest.
	delete esc; 
	delete ptrn;
	delete mpdb;
	delete hpt;
	free(pdb_paths); free(mcBPPS_prefix);
	for(i=1; i <= Number; i++){ if(IN_CMA[i]) TotalNilCMSA(IN_CMA[i]); } free(IN_CMA);
	NilAlpha(AB);
}

void	hsc_typ::Validate()
{
	Int4 Row,Col;
	// hpt vs cma
	if(hpt->NumSets() != Number) print_error("inconsistent *.hpt and *.mma input files");
	for(Row = 1; Row <= hpt->NumSets(); Row++){
	   if(IN_CMA[Row] == 0 && hpt->TypeOfSet(Row) == '?') continue;
	   if(strcmp(NameCMSA(IN_CMA[Row]),hpt->ElmntSetName(Row))){
		print_error("hsc_typ::Validate(): *.hpt and *.mma files are inconsistent");
	   }
	}
	if(IN_CMA[Number] == 0 || strcmp(NameCMSA(IN_CMA[Number]),"Random") != 0){
		print_error("*.mma file lacks a Random sequence set.");
	}
	// hpt vs ptrn.
	if(hpt->NumBPPS() != ptrn->NumPttrns){
		fprintf(stderr,"hpt->NumSets = %d != ptrn->NumPttrns = %d\n",hpt->NumBPPS(),ptrn->NumPttrns);
		print_error("hyperpartition and pttrns input files are inconsistent");
	}
	for(Col = 1; Col <= hpt->NumBPPS(); Col++){
	   Row= ptrn->NumPttrns - ptrn->PttrnCategory[Col] + 1;	// reversed order for *.pttrns file.
	   if(Row > hpt->NumBPPS()) print_error("*.hpt and *.pttrns files are inconsistent");
	}
	// patterns vs cma
	for(Int4 i=1; i <= Number; i++){
	   if(IN_CMA[i] == 0 && hpt->TypeOfSet(i) == '?') continue;
	    if(LengthCMSA(1,IN_CMA[i]) < ptrn->MaxPttrnPos){
		print_error("cma alignment files and pttrns input files are inconsistent");
	    }
	}
	   // char symbol=hpt->Cell(Row,Col);
	   // if(symbol != '+') print_error("*.hpt and *.pttrns files are inconsistent");
}

pat_typ *hsc_typ::SubPatterns(Int4 Row)
{ return ptrn->SubPatterns(hpt->RtnHyperPartition(Row)); }

int	hsc_typ::Run(char mode,char call,FILE *logfp,FILE *vsifp)
{
	// hpt->PutHyperPartition(stdout);
	// hpt->Put(stdout);
	// ptrn->Put(stdout);

	// if(call == 0) p2c->PrintVSI_Files();
	p2c->PrintVSI_Files(vsifp);
	int rtn=p2c->PrintKLST_Files(call);
        {
	   sprintf(string,"%s",mcBPPS_prefix); CheckString();
	   psc_typ psc(p2c,0,string,ptrn,esc); 
	   if(!verbose) psc.VerboseOff();
	   psc.FindBestPatternContacts2(mode,logfp);
	}
#if 0	// no longer using this...
	FILE *scrfp =0; if(verbose) scrfp=open_file(mcBPPS_prefix,".ca_scores","w");
	for(Int4 Row = 1; Row < hpt->NumSets(); Row++){
	   // skip last row == Random.
	   char c=hpt->TypeOfSet(Row);
	   // if(c != '!') continue;
	   if(IN_CMA[Row] == 0) continue;
	   if(c == '=') continue;
	   // pat_typ *sub_ptrn=SubPatterns(Row);
	   pat_typ *sub_ptrn=ptrn->SubPatterns(hpt->RtnHyperPartition(Row));
	   if(0) sub_ptrn->Put(stderr);
	   sprintf(string,"%s_%d",mcBPPS_prefix,Row); CheckString();
	   assert(strcmp(NameCMSA(IN_CMA[Row]),hpt->ElmntSetName(Row)) == 0);
	   // psc_typ psc(p2c,Row,string,sub_ptrn,esc); 
	   // psc.FindBestPatternContacts();
	   p2c->GetStructAlnScores(scrfp,Row,FALSE);
	   delete sub_ptrn;
	}
	// p2c->GetStructAlnScores(0,0,FALSE);
	if(scrfp) fclose(scrfp);
	// EvaluateStructAln();	// Structurally evaluate the entire alignment.
#endif
	return rtn;
}

void	hsc_typ::GetArg(Int4 argc, char *argv[])
// *************** Get arguments for all program options **********************
{
	Int4 arg,i;
	char	str[200];
        if(argc < 3) print_error(USAGE_START);
        for(arg = 0; arg < argc; arg++){
		fprintf(stdout,"%s ",argv[arg]);
	} fprintf(stdout,"\n");
	pdb_paths=AllocString(argv[1]);
	mcBPPS_prefix=AllocString(argv[2]);
        for(arg = 3; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'c':
		if(sscanf(argv[arg],"-color=%s",str)==1){ 
		   if(!isalpha(str[0])) print_error(USAGE_START);
		   TraceColor=AllocString(str);
		} else print_error(USAGE_START);
	     break;
             case 'C':
		if(sscanf(argv[arg],"-Color=%s",str)==1){ 
		   if(!isalpha(str[0])) print_error(USAGE_START);
		   SideColor=AllocString(str);
		} else print_error(USAGE_START);
	     break;
             case 's':
		if(sscanf(argv[arg],"-seed=%d",&seed)==1){ 
			if(seed < 1) print_error(USAGE_START);
		} else print_error(USAGE_START);
		break;
             default : print_error(USAGE_START);
           }
	 }
	}
}

