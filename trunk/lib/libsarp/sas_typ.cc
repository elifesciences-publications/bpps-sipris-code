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

#include "sas_typ.h"

//*************************** sas_typ ****************************************
//*************************** sas_typ ****************************************
//*************************** sas_typ ****************************************

void	sas_typ::Init( )
{
	// read in subgroup alignments...
	AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	Mode=' '; // Single model, Single value.
	// sprintf(string,"%s_new",mcBPPS_prefix); CheckString();
} 

void	sas_typ::ReadInput( )
{
	sprintf(string,"%s",mcBPPS_prefix); CheckString();
	FILE *fp = open_file(string,".cma","r");
	IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	assert(Number==1);
	assert(IN_CMA != 0);
	const char hpt_str[]="HyperParTition:\n!\n+ 1.Set1?\n- 2.Random=20000.\n\n";
        fp=tmpfile(); fprintf(fp,"%s",hpt_str); rewind(fp); hpt = new hpt_typ(fp); fclose(fp);
        // hpt->Put(stderr);
	IN_CMA[2]=MakeConsensusCMSA(IN_CMA[1]); Number++;	// WARNING: assumes array is longer!!!
	RenameCMSA("Set1",IN_CMA[1]); RenameCMSA("Random",IN_CMA[2]);

// fprintf(stderr,"IN_CMA[1]= '%s'\n",NameCMSA(IN_CMA[1]));
// fprintf(stderr,"IN_CMA[2]= '%s'\n",NameCMSA(IN_CMA[2]));
#if 1	// lengthen IN_CMA array 
	if(hpt->NumSets() > Number) {	// then lengthen the IN_CMA array using Null's
	   Int4	x,y;
	   cma_typ *in_cma; NEW(in_cma,hpt->NumSets() +3, cma_typ);
	   for(x=1,y=1; x <= hpt->NumSets(); x++)
	   // for(x=1,y=1; x < hpt->NumSets(); x++)
	   {
		if(strcmp(NameCMSA(IN_CMA[y]),hpt->ElmntSetName(x)) == 0){ // --> same set.
			in_cma[x]=IN_CMA[y];  y++;
		} else if(y == 1 || hpt->TypeOfSet(x) != '?'){
			print_error("*.hpt and *.cma files are fatally inconsistent");
		} 
	   } y--; 
	   if(y != Number) print_error("*.hpt and *.cma files are fatally inconsistent");
	   Number = hpt->NumSets(); free(IN_CMA); IN_CMA=in_cma; 
	}
#endif
	// read in patterns ...
	char	r,str[200]; 
	e_type  E=TrueSeqCMSA(1,IN_CMA[2]);
	fp=tmpfile();
	for(Int4 i=1; i <= LenSeq(E); i++){
		r=ResSeq(i,E);
		if(AlphaChar(r,AB) != 0){ fprintf(fp,"1: %c%d\n",AlphaChar(r,AB),i); break; }
	}  rewind(fp); ptrn = new pat_typ(fp); fclose(fp); 

	mpdb = new mps_typ(pdb_paths);
	if(mpdb->NumPDB <= 0) print_error("sas_typ::Init( ) fatal error: no pdb files created!");

	Validate();
}

void	sas_typ::Free()
{
	Int4 i;
	delete p2c;	// need to free this before the rest.
	delete esc; 
	delete ptrn;
	delete mpdb;
	delete hpt;
	for(i=1; i <= Number; i++){ if(IN_CMA[i]) TotalNilCMSA(IN_CMA[i]); } free(IN_CMA);
	NilAlpha(AB);
}

void	sas_typ::Validate()
{
	Int4 Row,Col;
	// hpt vs cma
	if(hpt->NumSets() != Number) print_error("inconsistent *.hpt and *.mma input files");
	for(Row = 1; Row <= hpt->NumSets(); Row++){
	   if(IN_CMA[Row] == 0 && hpt->TypeOfSet(Row) == '?') continue;
	   if(strcmp(NameCMSA(IN_CMA[Row]),hpt->ElmntSetName(Row))){
		hpt->Put(stderr);
		fprintf(stderr,"Row %d: IN_CMA name = '%s'",Row,NameCMSA(IN_CMA[Row]));
		fprintf(stderr," != hpt name = '%s'\n",hpt->ElmntSetName(Row));
		// assert(strcmp(NameCMSA(IN_CMA[Row]),hpt->ElmntSetName(Row)) ==0);
		print_error("sas_typ::Validate(): *.hpt and *.mma files are inconsistent");
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

double	sas_typ::Run()
{
	BooLean	PerformControl=FALSE;	PerformControl=TRUE;
	FILE *scrfp = open_file(mcBPPS_prefix,".ca_scores","w");
	p2c->ShowTheBest();
	double	rtn_value=0;
	for(Int4 Row = 1; Row < hpt->NumSets(); Row++){
	   // skip last row == Random.
	   char c=hpt->TypeOfSet(Row);
	   // if(c != '!') continue;
	   if(IN_CMA[Row] == 0) continue;
	   if(c == '=') continue;
	   sprintf(string,"%s_%d",mcBPPS_prefix,Row); CheckString();
	   assert(strcmp(NameCMSA(IN_CMA[Row]),hpt->ElmntSetName(Row)) == 0);
	   p2c->GetStructAlnScores(scrfp,Row,PerformControl);
	} rtn_value=p2c->FinalARSD;
	// p2c->GetStructAlnScores(0,0,PerformControl);
	fclose(scrfp);
	// EvaluateStructAln();	// Structurally evaluate the entire alignment.
	return rtn_value;
}

#define USAGE_START \
"USAGE: sas <pdb_name_file> <cmafile_prefix> [-options]\n\
    -target=<int>  - specify the target # seqs (>= 5) to use (default: 5)\n\
    -max=<int>	   - specify the maximum # of data points to use (default: 2000)\n\
    -mode=<char>   - run mode (default: ' ')\n\
                     'S' - sequence ARSD contributions\n\
                     'C' - column ARSD contributions\n\
    -seq=<int>     - find column contributions for pdb seq <int>\n\
    -dummy \n\
\n\n"

void	sas_typ::GetArg(Int4 argc, char *argv[])
// *************** Get arguments for all program options **********************
{
	Int4 arg,i;
	char	str[200];
        if(argc < 3) print_error(USAGE_START);
	FILE* fp=open_file(argv[2],"_sas.cmd","w");
        for(arg=0; arg < argc; arg++){ fprintf(fp,"%s ",argv[arg]); }
	fprintf(fp,"\n"); fclose(fp);
	pdb_paths=AllocString(argv[1]);
	mcBPPS_prefix=AllocString(argv[2]);
        for(arg = 3; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'm': 
		if(sscanf(argv[arg],"-max=%d",&MaxDataPoints)==1){ 
			if(MaxDataPoints < 10) print_error(USAGE_START);
                } else if(sscanf(argv[arg],"-mode=%c",&Mode)==1){ 
			if(Mode != 'C' && Mode != 'S') print_error(USAGE_START);
		} else print_error(USAGE_START); break;
             case 's': if(sscanf(argv[arg],"-seed=%d",&seed)==1){ 
			if(seed < 1) print_error(USAGE_START);
		} else if(sscanf(argv[arg],"-seq=%d",&KeySeq)==1){ 
			if(KeySeq < 1) print_error(USAGE_START);
			Mode='C';
		} else print_error(USAGE_START); break;
             case 't': if(sscanf(argv[arg],"-target=%d",&Target)==1){ 
			if(Target < 5) print_error(USAGE_START);
		} else print_error(USAGE_START); break;
             default : print_error(USAGE_START);
           }
	 }
	}
}

double	*sas_typ::ColumnARSD_Scores(FILE *fp,Int4 *&NN)
	{
	   Int4 ncols=LengthCMSA(1,IN_CMA[1]);
	   h_type HG=Histogram("change in ARSD scores (x10000)",-50,50,1.0);
	   double *A; NEW(A,ncols+9,double);
	   Int4	  x,*D; NEW(D,ncols+9,Int4);
	   dh_type dH=dheap(ncols+2,4);
	   for(x=0; x <= ncols; x++){
               	// fprintf(fp,"---- column %d -----\n",x);
		adh_typ *adh=0;
		double ARSD,arsd,d;
		Int4 N=1,DataPts,data_pts;
		if(x==0){
			adh=p2c->RoundOne(0,N); 
			ARSD=p2c->LastARSD; DataPts=p2c->LastNumData;
			fprintf(fp,"%d: ARSD=%.6f (%d pts)\n",
				x,ARSD,DataPts);
			d=ARSD*10000;
			A[x]=d; D[x]=DataPts;
		} else {
			adh=p2c->RoundOne(0,N,x,KeySeq); 
			arsd=p2c->LastARSD; data_pts=p2c->LastNumData;
			d=ARSD-arsd;  d*=10000;
			IncdHist(d,HG);
			insrtHeap(x,(keytyp)d,dH); A[x]=d; D[x]=data_pts;
			// fprintf(stdout,"%d: delta ARSD=%.6f (%d pts)\n",x,A[x],D[x]);
		}
               	double  key;
               	Int4    item,i;
               	adv_typ *adv;
               	for(i=0; (adv=adh->DelMin(&key,&item)) != 0; ){
                     		i++; delete adv;
               	} delete adh;
       	   } PutHist(fp,60,HG); NilHist(HG);
	   if(KeySeq > 0 && KeySeq <= p2c->NumPDB( )){
	     fprintf(fp," -------- %d.",KeySeq);
	     this->PutID(fp,KeySeq);
	     fprintf(fp," --------\n\ncol:   sqARSD   sqPairs\n");
	     for(x=1; x <= ncols; x++){
	        fprintf(fp,"%d\t%.6f\t%d\n",x,A[x],D[0]-D[x]);
	     } fprintf(fp,"\n\ncol:   sqARSD   sqPairs\n");
	     while(!emptyHeap(dH)){
		assert((x=delminHeap(dH)) != 0);
	        fprintf(fp,"%d\t%.6f\t%d\n",x,A[x],D[0]-D[x]);
	     }
	   } else { 
	     while(!emptyHeap(dH)){
		assert((x=delminHeap(dH)) != 0);
		fprintf(fp,"%d: delta ARSD=%.6f (%d pts)\n",x,A[x],D[x]);
	     } 
	   } 
	   Nildheap(dH); NN=D;
	   return A; 
}

double	*sas_typ::SeqARSD_Scores(FILE *fp,Int4 *&NN)
	{
	   Int4 ncols=LengthCMSA(1,IN_CMA[1]);
	   p2c->SetDataPoints(MaxDataPoints);
	   h_type HG=Histogram("change in ARSD scores (x10000)",-50,50,1.0);
	   Int4 M=p2c->NumPDB();
	   double *A; NEW(A,M+9,double);
	   Int4	  j,x,*D; NEW(D,M+9,Int4);
	   dh_type dH=dheap(M+2,4);
	   for(x=0; x <= M; x++){
               	// fprintf(fp,"---- column %d -----\n",x);
		adh_typ *adh=0;
		double ARSD,arsd,d;
		Int4 N=1,DataPts,data_pts;
		if(x==0){
			adh=p2c->RoundOne(0,N); 
			ARSD=p2c->LastARSD; // DataPts=p2c->LastNumData;
			DataPts=p2c->LastTotalData; 
			fprintf(fp,"%d: ARSD=%.6f (%d pts)\n",
				x,ARSD,DataPts);
			d=ARSD*10000;
			A[x]=d; D[x]=DataPts;
		} else {
			adh=p2c->RoundOne(0,N,-1,x); 
			arsd=p2c->LastARSD; // data_pts=p2c->LastNumData;
			data_pts=p2c->LastTotalData;
			d=ARSD-arsd;  d*=10000;
			IncdHist(d,HG);
			insrtHeap(x,(keytyp)d,dH); A[x]=d; D[x]=data_pts;
			// fprintf(stdout,"%d: delta ARSD=%.6f (%d pts)\n",x,A[x],D[x]);
		}
               	double  key;
               	Int4    item,i;
               	adv_typ *adv;
               	for(i=0; (adv=adh->DelMin(&key,&item)) != 0; ){ i++; delete adv; } delete adh;
       	   } PutHist(fp,60,HG); NilHist(HG);
	   fprintf(fp,"seq:   Sq_#  sqARSD   sqPairs\n");
	   for(j=0; !emptyHeap(dH); ){
		assert((x=delminHeap(dH)) != 0);
		if(D[0]==D[x]) continue; // This pdb sequence is not in the alignment.
		j++; this->PutID(fp,x);
		fprintf(fp,"\t%d\t%d\t%.6f\t%d\n",j,x,A[x],-(D[x]-D[0]));
	   } Nildheap(dH); NN=D;
	   return A; 
}

