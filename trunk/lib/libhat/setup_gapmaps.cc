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
#include "dft_typ.h"

#define MAX_NUM_CMA_FILES 2500

#if 0	// Moved to cmsa.h, cmsa.cc
void    PutMinimalSeqCMA(FILE *fp, Int4 sq, cma_typ cma)
// print out a concensus sequence cma file
{
 {      
	if(sq < 1 || sq > NumSeqsCMSA(cma)) return;
	if(nBlksCMSA(cma) > 1) print_error("PutMinimalSeqCMA( ) requires only 1 blk");

        gss_typ *gss = gssCMSA(cma);
        char	*aln;
        Int4    len=LengthCMSA(1,cma),pos[5],max_id_len=30;
	Int4	x = len + gss->NumIns(sq) + gss->NumDel(sq);
	NEW(aln,x+55,char);
        e_type sqE= TrueSeqCMSA(sq,cma);
        StrSeqID(aln, max_id_len, sqE);
        fprintf(fp, "[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",aln);
        fprintf(fp,"(%d)",len);
        for(Int4 i=1; i <= len; i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",len + gss->NumIns(sq) - gss->NumDel(sq),len);
        StrSeqID(aln, 30, sqE);
        fprintf(fp,">%s ",aln);
        StrSeqDescript(aln,50, sqE);
        fprintf(fp,"%s\n",aln);
        PosSiteCMSA(1, sq, pos, cma);
        Int4 LenStr=gss->Region(sq,aln,pos[1],len);
        fprintf(fp,"{()%s()}*\n\n_0].\n",aln);
	free(aln);
 }
}

cma_typ MinimizeFirstSeqCMSA(cma_typ cma)
{
        Int4    Number,n;
        a_type  A=AlphabetCMSA(cma);
        if(nBlksCMSA(cma) > 1) print_error("MinimizeFirstSeqCMSA( ) input error");

        // create concensus cma
        FILE *fp=tmpfile();
	PutMinimalSeqCMA(fp,1,cma);
	BooLean	*skip; NEW(skip, NumSeqsCMSA(cma) + 3, BooLean); skip[1]=TRUE;
	PutSelectOneCMSA(fp,skip,cma); free(skip);
        rewind(fp);
        cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,A);
        fclose(fp);

        // merge csq.cma and original cma.
        fp=tmpfile(); PutMergedCMSA(fp,Number,IN_CMA); rewind(fp);
        cma_typ CMA=ReadCMSA(fp,A); fclose(fp);
        for(n = 1; n<= Number; n++) TotalNilCMSA(IN_CMA[n]); free(IN_CMA);
        return CMA;
}
#endif

BooLean	CreateProfilesMAPGAPS(cma_typ *IN_CMA, Int4 num_cma_files, char *file_name, a_type A,
	double	Ethresh,double Ecutoff,double x_parameter, BooLean verbose)
// creates either a 
// If IN_CMA[0] == 0 then input is a false positive file (without a template).
{
	BooLean	success=FALSE,IsFalsPos=FALSE;
	e_type	qE,QueryE[MAX_NUM_CMA_FILES],trueE,fakeE;
	ss_type	data,Data[MAX_NUM_CMA_FILES];
	sap_typ	sap;
	cma_typ	cma;
	Int4    T=11;
	char 	str[100],*Checkpoint;
	Int4	II,maxrounds,Begin=0;
	char	*checkin=0,*checkout=0;

	//***************************** 2. Initialize search *******************************
  assert(num_cma_files < MAX_NUM_CMA_FILES);
  for(II=0; II <= num_cma_files; II++){ QueryE[II]=0; Data[II]=0; }
  if(IN_CMA[0] == 0){	// input file 
    Begin=1;
    FILE *fp=open_file(file_name,".xpq","w");
    for(II=1; II <= num_cma_files; II++){
	cma=IN_CMA[II];
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		// need to renumber sequence ids in alignment for blast.
		EqSeqI(sq,TrueSeqCMSA(sq,cma)); EqSeqI(sq,FakeSeqCMSA(sq,cma));
	}
	QueryE[II]=CopySeq(TrueSeqCMSA(1,cma));
	ChangeInfoSeq(NameCMSA(cma),QueryE[II]); 
	PutSeq(fp,QueryE[II],A);
	Data[II]=TrueDataCMSA(cma);
    } fclose(fp);
    sprintf(str,"%s.xup",file_name);	// checkpoint file name
  } else {
    Begin=0;
    for(II=0; II <= num_cma_files; II++){
	cma=IN_CMA[II];
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		// need to renumber sequence ids in alignment for blast.
		trueE = TrueSeqCMSA(sq,cma);
		fakeE = FakeSeqCMSA(sq,cma);
		EqSeqI(sq,trueE); EqSeqI(sq,fakeE);
	}
	Data[II]=TrueDataCMSA(cma);
	trueE=TrueSeqCMSA(1,cma); fakeE=FakeSeqCMSA(1,cma);
	if(!IdentSeqs(trueE,fakeE) || LenSeq(trueE) != LengthCMSA(1,cma) 
			|| nBlksCMSA(cma) != 1){
		   fprintf(stderr,"file %s_%d.cma: \n\n",file_name,II);
		   fprintf(stderr,"trueE = %d aa; cma = %d aa \n\n",
					LenSeq(trueE),LengthCMSA(1,cma));
		   PutSeq(stderr,trueE,A); // PutSeq(stderr,fakeE,A);
		   print_error("incompatible template consensus sequence");
	}
	QueryE[II]=CopySeq(trueE);
	ChangeInfoSeq(NameCMSA(cma),QueryE[II]); 
    }
    // sprintf(str,"%s.map",file_name);
    sprintf(str,"%s.mpa",file_name);
  }
  Checkpoint = AllocString(str);

  //******************************** 3. Main search loop  ******************************
  FILE *chkfp=open_file(Checkpoint,"","w");
  for(II=Begin; II <= num_cma_files; II++){ // start main loop over all *cma files....
	checkin=0; 
	checkout=Checkpoint; 
	cma=IN_CMA[II]; 
	data=Data[II]; 
	qE=CopySeq(QueryE[II]); 
	EqSeqI(0,qE);
	maxrounds=2;
	UInt4	hpsz=NumSeqsCMSA(cma)+3;
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,maxrounds,checkin);
	gpsi.SpeakUp( );
	if(!verbose) gpsi.KeepQuite();
	sap = gpsi.FakeSearchCMSA(T,cma); 
	if(sap) success=TRUE; else break;
	gpsi.ComputeMatrix(chkfp, checkout);
	NilSeq(QueryE[II]);
   } // end loop over cma files
   fclose(chkfp);
   return success;
}

#define USAGE3 "Usage: mapgaps <prefix> [options]\n\
      NOTE: Put concatenated cma alignments in <prefix>.cma \n\
            with the last being the Template alignment.\n\
      options:\n\
        -exclude  program will look for a *.xpa (excluded profile alignments) input file \n\
                  which will be used to create a *.xpq (query) and a *.xup (unaligned profiles) files.\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
\n"

int	public_setup_mapgaps(int argc, char *argv[])
{ return main_setup_gapmaps(argc, argv,USAGE3); }

#define USAGE2 "usage: run_map in_file [options]\n\
      NOTE: Put concatenated cma alignments in in_file.cma \n\
            with the last being the Template alignment.\n\
      options:\n\
        -exclude  program will look for a *.xpa (excluded profile alignments) input file \n\
                  which will be used to create a *.xpq (query) and a *.xup (unaligned profiles) files.\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
\n"

int	public_setup_gapmaps(int argc, char *argv[])
{ return main_setup_gapmaps(argc, argv,USAGE2); }

#define USAGE "usage: mkmaps in_file [options]\n\
      NOTE: Put concatenated cma alignments in in_file.cma \n\
	    Put the Template alignment in <in_file>.tpl.\n\
	    Outfiles are: <in_file>.mpa (the multiple profile alignment)\n\
      options:\n\
        -exclude  program will look for a *.xpa (excluded profile alignments) input file \n\
                  which will be used to create a *.xpq (query) and a *.xup (unaligned profiles) files.\n\
        -m<int>   print option m=0,1,4 (default 0)\n\
        -v        verbose mode\n\
        -X<int>   X dropoff for gapxdrop\n\
        -Z        Print cma input names without running search\n\
        -z        see posMatrix\n\
\n"

// 	Don't need the following...
//         -False=<str> input the 'false positive' cma file so that number can be computed\n 


#if 0	//**********************************************************************
#endif  //**********************************************************************
int	setup_gapmaps(int argc, char *argv[])
{ return main_setup_gapmaps(argc, argv,USAGE); }

int	main_setup_gapmaps(int argc, char *argv[],const char *usage,FILE *ifp)
{
	Int4    T=11;
	time_t	time1=time(NULL);
        double	x_parameter=25.0;
	Int4	x;
	UInt4	hpsz=100000;
	Int4	printopt=0;
	double	Ethresh=0.001,Ecutoff=0.001;
	BooLean	see_smatrix=FALSE;
	BooLean	success=FALSE;
	BooLean print_names_only=FALSE,verbose=FALSE;
	Int4	num_cma_files,start_arg=3,num_false_pos=0;
	char	FalsePosFile[200];
	char 	str[300];
	char	FileName[300];

	BooLean Open_fpa_file=FALSE;
	Int4    NumFalsePos=0,NumTrueCMAs=0,NumFalseCMAs=0,TotalNumCMA_Files=0;
	cma_typ *FALSE_CMA;

	if(argc < 2) print_error(usage);
	//***************************** 1. Read input files *******************************
	start_arg=2;
	a_type	A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // IMPORTANT!!
	// use distinct *.tpl  input file...
	sprintf(str,"%s.tpl",argv[1]);
	cma_typ tpl_cma=ReadCMSA2(str,A);
	FILE    *fp=0;
	cma_typ *IN_CMA=0;

	// check for input cma file and if doesn't exist create it from the tpl input file.
	if(ifp==0){
	   sprintf(str,"%s.cma",argv[1]);
	   if((fp=fopen(str,"r")) == NULL){
	        // fprintf(stderr,"WARNING: %s.cma file does not exist;",argv[1]);
	        // fprintf(stderr," creating from template alignment\n\n",argv[1]);
	        fp=open_file(argv[1],".cma","w");	// then create the file.
		Int4 x,y,sq,N=NumSeqsCMSA(tpl_cma),MaxLen=0;
		e_type sE=0,fE=0;
		for(sq=1; sq <= N; sq++){
		   fE=FakeSeqCMSA(sq,tpl_cma);
		   x=NonXResSeq(fE,A); y=LenSeq(fE);
		   double d=(double)x/(double)y;
		   if(d < 0.66){
			fprintf(stderr,"Sequence %d in template file contains %d gaps\n",sq,x);
			print_error("Fatal: edit template alignment to eliminate poor quality seqs.");
		   } sE=TrueSeqCMSA(sq,tpl_cma); MaxLen=MAXIMUM(Int4,LenSeq(sE),MaxLen);
		} char *buffer; NEW(buffer,MaxLen+9,char);
		for(sq=2; sq <= N; sq++){
		        sE=TrueSeqCMSA(sq,tpl_cma); 
			char Str[202];
			StrSeqID(Str, 200, sE);
			Int4 len=SeqToString(buffer, sE, A);
			fprintf(fp,"[0_(1)=%s(2){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n(%d)", Str,len);
			for(Int4 x=1; x <= len; x++) fprintf(fp,"*"); fprintf(fp,"\n\n");
			fprintf(fp,"$%d=%d(%d):\n>%s \n{()%s()}*\n\n",sq,len,len,Str,buffer);
			fprintf(fp,"$%d=%d(%d):\n>%s \n{()%s()}*\n\n",sq,len,len,Str,buffer);
			fprintf(fp,"\n_0].\n");
		}
	   } fclose(fp);
	   fp=open_file(argv[1],".cma","r");
	   IN_CMA=MultiReadCMSA(fp,&num_cma_files,A); fclose(fp);
	} else { IN_CMA=MultiReadCMSA(ifp,&num_cma_files,A); }
	if(num_cma_files < 1){
		print_error("FATAL: 0 input files!"); 
		print_error(usage);
	} else if(num_cma_files >= MAX_NUM_CMA_FILES){
		print_error("FATAL: too many input files!");
	} IN_CMA[0]=tpl_cma;
	if(NumSeqsCMSA(IN_CMA[0]) == num_cma_files + 1) {
		sprintf(FileName,"%s",argv[1]);
	} else if(NumSeqsCMSA(tpl_cma) > num_cma_files + 1) {	// salvageable error.
		Int4	i,j,x,NumMissing=0;
		Int4	N=NumSeqsCMSA(tpl_cma),TotalMissing;
		sprintf(FileName,"%s_M",argv[1]);
		BooLean *skip; NEW(skip,N+3,BooLean);
		TotalMissing= N-num_cma_files-1;
		Int4	*OrphanE; NEW(OrphanE,TotalMissing+3,Int4);
		fprintf(stderr,"WARNING: fewer alignments (%d) in %s.cma than template seqs (%d) in %s.tpl!\n",
				num_cma_files,argv[1],N-1,argv[1]);
		for(j=i=1; i < N; j++,i++){
			e_type cE,tE=TrueSeqCMSA(i+1,tpl_cma);
			if(j > num_cma_files){ 	// missing sequences are at the end of tpl file.
				skip[i+1]=TRUE; NumMissing++;
				OrphanE[NumMissing]=i+1;
				if(NumMissing > TotalMissing) break;
		        } else {
			   cE=TrueSeqCMSA(1,IN_CMA[j]);
			   if(!IdentSeqs(tE,cE)){
				skip[i+1]=TRUE; NumMissing++;
				OrphanE[NumMissing]=i+1;
				if(NumMissing > TotalMissing) break;
				j--; // try again with next template seq.
			   }
			}  
		} 
		if(NumMissing != TotalMissing){
			fprintf(stderr,"Fatal inconsistency between %s.cma & %s.tpl files\n",argv[1],argv[1]);
			print_error("exiting...\n");
		} else {
			fprintf(stderr,"===> Attempting to correct the problem...\n");
			for(i=1; i<= NumMissing; i++){
				e_type tE=TrueSeqCMSA(OrphanE[i],tpl_cma);
				fprintf(stderr,"---------- %d.'orphaned' template sequence: ---------\n",OrphanE[i]);
				PutSeq(stderr,tE,A); 
			}
		}
		fprintf(stderr,"   ...(creating %s.tpl file).\n\n",FileName);
		fp=tmpfile(); PutSelectCMSA(fp,skip,tpl_cma); rewind(fp); 
		TotalNilCMSA(tpl_cma); tpl_cma=ReadCMSA(fp,A); fclose(fp);
		sprintf(str,"%s.tpl",FileName);	// new file name.
		WriteCMSA(str,tpl_cma);	IN_CMA[0]=tpl_cma;
		free(skip); free(OrphanE);
	} else if(NumSeqsCMSA(IN_CMA[0]) < num_cma_files + 1) {	// Fatal error.
		fprintf(stderr,"\n*********** Template input error! ***********\n\n");
		fprintf(stderr,"Greater number of alignments (%d) than template sequences (%d)\n\n",
			num_cma_files,NumSeqsCMSA(IN_CMA[0])-1);
		{
		  Int4 i,smaller;
		  smaller=MINIMUM(Int4,num_cma_files-1,NumSeqsCMSA(IN_CMA[0]));
		  for(i=1; i <= smaller; i++){
			e_type tE=TrueSeqCMSA(i+1,tpl_cma);
			e_type cE=TrueSeqCMSA(1,IN_CMA[i]);
			if(!IdentSeqs(tE,cE)){
				fprintf(stderr,"-------------- %d.Sequence mismatch: --------------\n",i);
				PutSeq(stderr,tE,A); PutSeq(stderr,cE,A);
				print_error(usage);
			}
		  }
		}
		print_error(usage);
	}
#if 1	// Check for truncations at either end and modify cma files accordingly.
   TrimToTemplateCMSA(IN_CMA,num_cma_files);
#else
   for(Int4 s=1; s <= num_cma_files; s++){ 
        fprintf(stderr,"================================ %d: %s.\n", s,NameCMSA(IN_CMA[s])); 
	e_type tplSq=TrueSeqCMSA(s+1,tpl_cma);
	if(LenSeq(tplSq) > LengthCMSA(1,IN_CMA[s])) print_error("input error");
	if(LenSeq(tplSq) != LengthCMSA(1,IN_CMA[s])){	// then need to change input cma file.
		Int4 Start,i;
		cma_typ cmaX=IN_CMA[s];
		e_type csqSq=TrueSeqCMSA(1,cmaX);
		char rtn=IsSubSeq(tplSq,csqSq,&Start,FALSE);
		// rtn = 1 if tplSq is a subseq of csqSq.
		if(rtn != 1) print_error("Template and cma files are incompatible");
		if(Start > 0){					// remove N-terminal columns.
			for(i=1; i<=Start; i++){
                           if(LengthCMSA(1,cmaX) <= 3) print_error("input error");
                           RmColumnMSA(1,1,cmaX); // block 1, first column removed.
                        }
		}
		if(LenSeq(tplSq) < LengthCMSA(1,cmaX)) {	// remove C-terminal columns.
			Int4 lenrm = LengthCMSA(1,cmaX) - LenSeq(tplSq);
                        for(i=1; i<=lenrm; i++){
                           Int4 lemon = LengthCMSA(1,cmaX);
                           if(lemon <= 3) print_error("input error");
                           RmColumnMSA(1, lemon, cmaX);
                        }
		}
		IN_CMA[s] = MinimizeFirstSeqCMSA(cmaX); TotalNilCMSA(cmaX);
	}
   }
#endif

	for(Int4 arg=0; arg < argc; arg++) { fprintf(stderr,"%s ",argv[arg]); } fprintf(stderr,"\n\n");
	// fprintf(stderr,"num_cma_files = %d\n",num_cma_files);
	//***************************** 2. Get input options *******************************
        for(Int4 arg = start_arg; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(usage);
           switch(argv[arg][1]) {
	     case 'F': {	
		cma_typ fcma=0;
		  if(sscanf(argv[arg],"-False=%s",FalsePosFile) == 1){
			fcma=ReadCMSA2(FalsePosFile,A);
			num_false_pos = NumSeqsCMSA(fcma);
			if(num_false_pos < 1 || num_false_pos >= num_cma_files){
				print_error(usage);
			} else NilCMSA(fcma);
		  } else print_error(usage);
		} break;
	     case 'e': 	
		if(argv[arg][8] == 0 && strcmp(argv[arg],"-exclude") == 0){
			Open_fpa_file=TRUE;
#if 0
		} else if(sscanf(argv[arg],"-false=%d",&num_false_pos) == 1){
			if(num_false_pos < 1 || num_false_pos >= num_cma_files){
				print_error(usage);
			}
#endif
		} else print_error(usage);
		break;
	     case 'm': 	
		{
		  if(argv[arg][2] != '='){
			printopt=IntOption(argv[arg],'m',0,6,usage); 
		  } else print_error(usage);
		} break;
             case 'v': verbose=TRUE; break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,usage); break;
             case 'Z': print_names_only=TRUE; break;
             case 'z': see_smatrix= TRUE; break;
             default: print_error(usage);
           }
        }

  if(print_names_only){
   for(Int4 s=1; s <= num_cma_files; s++){ 
     fprintf(stderr,"================================ %d: %s.\n", s,NameCMSA(IN_CMA[s])); 
   } 
  }
  if(print_names_only) exit(0);
  if(Open_fpa_file){
	fp=open_file(FileName,".xpa","r");
	FALSE_CMA=MultiReadCMSA(fp,&NumFalseCMAs,A);
	fclose(fp);
  } else FALSE_CMA=0;

#if 1	// create .dft file here for tree...
  // fp=open_file(FileName,".dft","r");
  sprintf(str,"%s.dft",FileName);	// new file name.
  fp = fopen(str,"r");
  if(fp == NULL) {	// then create a dft file.
	Int4 i,j,k,N,*Level;
	N=num_cma_files + NumFalseCMAs;
	NEW(Level,N+4,Int4);
	for(i=0; i <= num_cma_files; i++){
		Level[i]=LevelCMSA(IN_CMA[i]); 
	}
	for(j=1; j <= NumFalseCMAs; j++,i++){
		Level[i]=1; 
		// Level[i]=LevelCMSA(FALSE_CMA[j]);
	}
	dft_typ dft(i, Level);
	fp=open_file(FileName,".dft","w");
	dft.Write(fp); fclose(fp);
  } else fclose(fp);
#endif
  // check to make sure that input file is okay.
  // WARNING: the next line appears to be unnecessary and can probably be deleted at some point.
  if(num_false_pos > 0) SetLevelCMSA(num_false_pos,IN_CMA[0]);
  // else use level given in input file.
//   sprintf(str,"%s.tpl",FileName);
//   WriteCMSA(str,IN_CMA[0]);	// write out template file...

  if(num_false_pos > 0){	// then create a template profile lacking the false positives.
	BooLean *skip=0;
	Int4 N = NumSeqsCMSA(IN_CMA[0]);
	NEW(skip,N+3,BooLean);
	for(Int4 i=N - num_false_pos + 1; i <= N; i++) skip[i]=TRUE;
	FILE *tmpfp=tmpfile();
	PutSelectCMSA(tmpfp,skip,IN_CMA[0]); free(skip);
	rewind(tmpfp);
	TotalNilCMSA(IN_CMA[0]);
	IN_CMA[0]=ReadCMSA(tmpfp,A);
	fclose(tmpfp);
  }
  // sprintf(str,"%s.map",FileName);
  // Checkpoint = AllocString(str);
  success=CreateProfilesMAPGAPS(IN_CMA,num_cma_files,FileName,A,Ethresh,Ecutoff,x_parameter,verbose);
  if(FALSE_CMA){
     success=CreateProfilesMAPGAPS(FALSE_CMA,NumFalseCMAs,FileName,A,Ethresh,Ecutoff,x_parameter,verbose);
     for(Int4 II=1; II <= NumFalseCMAs; II++) TotalNilCMSA(FALSE_CMA[II]);
  }
   //************************* 11. output Main alignment file ************************
   for(Int4 II=0; II <= num_cma_files; II++) TotalNilCMSA(IN_CMA[II]);
   NilAlpha(A);
   double runtime = difftime(time(NULL),time1);
   fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
   if(success) return 0;
   else return 1;
}


