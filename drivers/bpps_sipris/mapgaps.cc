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

#include "hat_typ.h"
#include "mgs_typ.h"
#include "cma_gmb.h"
#include <time.h>

#define USAGE_START "*** Multiply Aligned Profiles for Gapped Alignment of Protein Sequences ***\n\
  Usage 1(for creating profiles): mapgaps <prefix> \n\
    Input: <prefix>.tpl: Template alignment in cma format.\n\
	   <prefix>.cma: concatenated alignments in cma format.\n\
           (Use the fa2cma program to convert fasta files to cma format.)\n\
     Note: If <prefix>.tpl and <prefix>.cma files are unavailable,\n\
            a fasta format alignment file (<prefix>.fa) may be used.\n\
     Note: If <prefix>.cma is missing, profiles will be created from the template\n\
           with output stored in files <prefix>_X.tpl and <prefix>_X.cma;\n\
	   call 'mapgaps <prefix>_X <database>' for the search step.\n\
    Output: <prefix>.mpa: profiles required for searches.\n\
            <prefix>.dft: Depth first search tree.\n\
    Options: none\n\
  Usage 2(for searches): mapgaps <prefix> <database> [options]\n\
    Input: <prefix>.mpa: profiles.\n\
           <prefix>.tpl: Template alignment in cma format.\n\
           <prefix>.dft: Optional depth first search tree.\n\
    Output: <database>_A.cma: alignment of sequences found.\n\
            <database>_A.seq: fasta file of sequences found.\n\
    Options: type 'mapgaps -' for search options\n\
  References:\n\
    Neuwald, A.F. 2009. Rapid detection, classification and accurate alignment\n\
      of up to a million or more related protein sequences.  Bioinformatics 25: 1869-1875\n\
    Neuwald, A.F., L. Aravind & S.F. Altschul. 2017. Inferring Joint Sequence-Structural\n\
      Determinants of Protein Functional Specificity. eLife In revision.\n\
\n"

const char USAGE[]="Usage: mapgaps <query_prefix> <dbsfile> [options]\n\
    Input: <prefix>.mpa: profiles.\n\
           <prefix>.tpl: Template alignment in cma format.\n\
           <prefix>.dft: Optional depth first search tree.\n\
    Output: <dbsfile>_A.cma: alignment of sequences found.\n\
            <dbsfile>_A.seq: fasta file of sequences found.\n\
    Options:\n\
        -D        Don't mask low complexity regions in template sequences.\n\
        -e=<real> E-value cutoff (after accounting for the database size).\n\
                  default: 0.01/number of profiles\n\
        -H=<int>  heapsize (default 500000)\n\
	-I=<int>:<int> - lengths of left & right flanking regions to be retained (default: 10:10)\n\
        -l        do low complexity masking of database sequences\n\
        -L        do low complexity and coiled coil masking of database sequences\n\
        -O        output sequence hits in separate files for each profile\n\
	-Refine=<int1>:<int2> Refine the output alignment by retaining \n\
                  sequences with no more than <int1> percent identity and \n\
		  with more than <int2> percent aligned residues. \n\
                  (default: 100% identity; 0% aligned residues).\n\
	    Note: All sequences are retained by default & pdb sequences are always retained.\n\
        -sense=<real>  sensitivity parameter for heuristic search \n\
                   range: 0.01...1.0; lower values: > sensitivity (default: 0.05)\n\
        -size=<int>  set effective length of the database (default: actual database size)\n\
        -T=<int>  blast word hit threshold (default: 11)\n\
        -trigger=<real>  Number of bits to trigger gapping = <real> (default: 22.0)\n\
        -X=<int>  X dropoff for BLAST gapxdrop (default: 25)\n\
\n";

static void     PrintError()
{ PrintLicenseStatement("mapgaps v2.0.3"); print_error(USAGE_START); }

cma_typ	MinColPurgeCMSA(Int4 percent_ident, Int4 mincol, cma_typ cma)
{
	a_type AB=AlphabetCMSA(cma);
	assert(nBlksCMSA(cma) == 1);
	assert(mincol >= 0 && mincol <= 100);
	double	fr,MinCol=(double)mincol/100.0;
	cma_typ	xcma=0;
	FILE	*fp=0;
	e_type	qE,sE;

	// -mincol= option
	Int4	i,j,k,s,J,I,n,r,Len,na;
	Int4    N=NumSeqsCMSA(cma);
        BooLean *skip;
      if(mincol > 0){
	h_type HG=0; // HG=Histogram("fraction aligned",0,1,0.025);
	NEW(skip,N+3,BooLean);
        for(J=1; J <= N; J++){
	   sE=TrueSeqCMSA(J,cma); if(!PdbSeq(sE)) skip[J]=TRUE; 
	}
        for(Len=LengthCMSA(1,cma),n=0,J=1; J <= N; J++){
		for(na=0,s=1 ; s <= Len;s++){ 
		   r=ResidueCMSA(1,J,s,cma); if(r != UndefAlpha(A)) na++; 
		} fr=(double)na/(double)Len;
		if(fr >= MinCol){ skip[J]=FALSE; n++; }
		if(HG) IncdHist(fr, HG);
	} 
	if(n <= 0) print_error("Alignment is too fragmented to be used");
	fp=tmpfile(); PutSelectCMSA(fp,skip,cma); 
	rewind(fp); xcma=ReadCMSA(fp,AB); fclose(fp);
	if(HG){ PutHist(stderr,60,HG); NilHist(HG); } free(skip);
      } else xcma=cma;

	// -U= option...
        set_typ InSet=MakeSet(NumSeqsCMSA(xcma)+4); FillSet(InSet);
        set_typ Set=RtnFastRepSetCMSA(0, percent_ident,InSet,xcma);
	N=NumSeqsCMSA(xcma); 
	for(i=1; i <= N; i++) if(PdbSeq(TrueSeqCMSA(i,xcma))) AddSet(i,Set);	// keep pdb seqs.
	fp=tmpfile(); PutInSetCMSA(fp,Set,xcma);
	NilSet(Set); NilSet(InSet); 
	if(mincol > 0) TotalNilCMSA(xcma);
	rewind(fp); xcma=ReadCMSA(fp,AB); fclose(fp);
	return xcma;
}

BooLean	IsValidConsensusCMSA(Int4 sq, cma_typ cma)
{
	if(sq < 1 || sq > NumSeqsCMSA(cma)) return FALSE;
	gsq_typ *gsq = gsqCMSA(sq,cma);
	if(gsq->Gapped( )) return FALSE;
	if(LenSeq(gsq->FakeSeq()) != LenSeq(gsq->TrueSeq())) return FALSE;
	return TRUE;
}

extern int     fa2cma_main(int argc, char *argv[],FILE *ofp);

int	main(int argc, char *argv[])
// this is the same routine as the public version, except for the usage...
// 
{
   time_t  time1=time(NULL);
   BooLean success,seqs_out=FALSE,FromFA=FALSE;
   Int4    mincol=0,percent_ident=100;
   FILE    *fp=0,*ifp=0,*tfp=0;
   char    str[300];

   if(argc < 2) PrintError(); 
   if(argc == 2 && argv[1][0]=='-'){
         TurnOffLicenseStatement();
	 mgs_typ mgs(argc,argv,USAGE);
   }
   //=========== check for input file ============
   sprintf(str,"%s.tpl",argv[1]); fp=fopen(str,"r"); 
   if(fp == NULL){
        sprintf(str,"%s.fa",argv[1]); fp=fopen(str,"r"); 
        if(fp == NULL) print_error(USAGE_START);  
	else {
	   int Argc=2; FromFA=TRUE;
	   char *Argv[5]; Argv[0]=argv[0]; 
	   sprintf(str,"%s.fa",argv[1]); Argv[1]=AllocString(str);
	   FILE *ofp=open_file(argv[1],".tpl","w"); 
	   fa2cma_main(Argc,Argv,ofp); 	fclose(ofp); free(Argv[1]);
	}
   } else fclose(fp);

   if(argc == 2){	// create profiles.
      TurnOffLicenseStatement();
      sprintf(str,"%s.cma",argv[1]);
      a_type AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
      if(!FromFA && (fp=fopen(str,"r")) != 0){   // tpl & cma files are present...
        Int4 Number,I;
        cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
        cma_typ *OUT_CMA; NEW(OUT_CMA,Number +3, cma_typ);
        for(I=1; I <= Number; I++){ OUT_CMA[I]=RmWrinklesCMSA(IN_CMA[I]); } ifp = tmpfile();
        for(I=1; I <= Number; I++) PutCMSA(ifp,OUT_CMA[I]); rewind(ifp);
        for(I=1; I <= Number; I++){ NilCMSA(OUT_CMA[I]); TotalNilCMSA(IN_CMA[I]); }
	main_setup_gapmaps(argc, argv,USAGE_START,ifp); if(ifp) fclose(ifp);
      } else {			     // create a cma file from tpl; add concensus.
	fprintf(stderr,"WARNING: %s.cma file does not exist;",argv[1]);
        fprintf(stderr," creating from %s.tpl.\n\n",argv[1]);
	sprintf(str,"%s.tpl",argv[1]);
	cma_typ cma=ReadCMSA2(str,AB);
#if 1
	if(nBlksCMSA(cma) > 1) print_error("Template alignment input error");
        Int4 i,j,x,z,r,N=NumSeqsCMSA(cma);
        e_type sE=0,fE=0;
	z=LengthCMSA(1,cma);
        for(Int4 sq=1; sq <= N; sq++){
            fE=FakeSeqCMSA(sq,cma);
	    for(i=1,x=0; i <= z; i++){
	       r=ResidueCMSA(1,sq,i,cma);
	       if(r > 0) x++;
	    }
            double d=(double)x/(double)z;
            if(d < 0.66){
                  fprintf(stderr,"FATAL: Sequence %d in %s.tpl contains %d/%d deletions...\n",
			sq,argv[1],z-x,z);
                  print_error("    ..edit template alignment to eliminate poor quality seqs.\n");
            }
	}
#endif
	cma_typ cma0 = AddConsensusCMSA(cma); TotalNilCMSA(cma);
	BooLean AddX2Ends=TRUE;
	double	fractDeleted=0.75;
	cma=RemoveOverhangsCMSA(cma0,AddX2Ends); TotalNilCMSA(cma0);
	Int4 col_del=RmGappyColumnsCMSA(fractDeleted, cma);
        if(col_del > 0) fprintf(stderr,"\n %d columns removed",col_del);
	sprintf(str,"%s_X",argv[1]);
	fp=open_file(str,".tpl","w"); 
	PutCMSA(fp,cma); TotalNilCMSA(cma); fclose(fp);
	char *tmp=argv[1]; argv[1]=AllocString(str);
	main_setup_gapmaps(argc, argv,USAGE_START,0); free(argv[1]); argv[1]=tmp;
      } NilAlpha(AB);
   } else if(argc > 2 && argv[2][0] == '-'){
      print_error(USAGE_START); 
   } else if(argc >= 3){	// perform gapmaps search.
     sprintf(str,"%s.mpa",argv[1]);
     if((fp=fopen(str,"r")) == 0){   // running without an *.mpa file...
        fprintf(stderr,"%s\n",USAGE_START);
        fprintf(stderr,"**********************************************************************\n");
        fprintf(stderr,"** No %s.mpa file provided. Run 'mapgaps %s' (Usage 1) **\n",
                        argv[1],argv[1]);
        fprintf(stderr,"**********************************************************************\n");
	exit(1);
     } else {
        fclose(fp); TurnOffLicenseStatement();
	Int4	x,z;
	for(Int4 arg = 3; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE);
	   switch(argv[arg][1]) {
		case 'R': if(sscanf(argv[arg],"-Refine=%d:%d",&x,&z) == 2){
                      if(x < 10 || x > 100) print_error(USAGE);
                      if(z < 25 || z > 100) print_error(USAGE);
                      mincol=z; percent_ident=x;  argv[arg][1]=' '; } break;
		default: break;
	   }
	}
        mgs_typ mgs(argc,argv,USAGE);
        a_type AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	tfp=tmpfile(); success=mgs.Run(seqs_out,tfp); rewind(tfp);
	cma_typ xcma=0,cma=ReadCMSA(tfp,AB); fclose(tfp);
	if(mincol > 0 && percent_ident < 100){
	  xcma=MinColPurgeCMSA(percent_ident,mincol,cma);
	  TotalNilCMSA(cma); cma=0;
	} if(cma){ xcma=cma; cma=0; }
	if(xcma != 0){
	   fp=open_file(argv[2],"_A.mma","w"); PutCMSA(fp,xcma); fclose(fp); TotalNilCMSA(xcma); 
	} NilAlpha(AB);
     }
  } else { print_error(USAGE_START); } // will print out usage
  double runtime = difftime(time(NULL),time1);
  fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
  if(success) return 0; else return 1;
}

