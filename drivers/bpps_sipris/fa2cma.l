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

/** fa2cma **/
%START DEF SEQ STRT MSF MFASTA
A       [a-zA-Z]
U       [A-Z]
L       [a-z]
D       [-]
P       [.]
R       ({D}|{P}|{A})+
S       [ \t]+
%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <sequence.h>
#include <stdinc.h>
#include <afnio.h>
#include "cmsa.h"
#include "cma_gmb.h"
#include "dms_typ.h"

#define YYLMAX  1000000
#define print_error(str) (fprintf(stderr,"%s\n",str)? exit(1): exit(1)) 

unsigned long    num=0,TotalSeq=0,LenAln=0;
unsigned long	 LenCMA=0,LenTrue=0,LenFake=0;
unsigned long    NumP=0,NumD=0,NumU=0,NumL=0,NumAln=0,NumGap=0,NumRF=0;
char	SeqCMA[YYLMAX];
char	SeqDef[YYLMAX];
unsigned long Begin=ULONG_MAX,End=0,Length;
e_type	DummySeq=0;
long SeqOffSet=0,SeqExtend=0;
char Kingdom,Phylum[100];
FILE	*ofp=stdout;

Int4	RmQueryGapsCMSA(cma_typ &cma)
// remove columns that are deleted in the first sequence.
{
    Int4	sq,hits,sq_hits,s,pos[4],len=LengthCMSA(1,cma),N=NumSeqsCMSA(cma);
    Int4	i,j,total=0,cycle;
    // a_type	A=AlphabetCMSA(cma);
    BooLean	*Delete;
    cma_typ	rcma=0;

#if 0
    if((i=Check4NtermExtendCMSA(cma)) > 0){
	fprintf(stderr,"%d out of %d sequences have N-term extensions\n",i,N);
	print_error("Fatal: input file must not have N-term extensions");
    }
#endif
#if 1	// fixes problem with adding insertions!
    ExtendFakeToRealCMSA(cma);
#endif
#if 0	// 
    RmUnAlignedSeqsCMSA(cma);
#endif
#if 0	// 
    rcma=RmWrinklesCMSA(cma);
    if(rcma){ NilCMSA(cma); cma=rcma; }
#endif
    NEW(Delete, LengthCMSA(1,cma)+3, BooLean);
    // Don't remove columns on ends!!
    for(s=2 ; s < LengthCMSA(1,cma); s++){ Delete[s]=IsDeletedCMSA(1,1,s,cma); }
    cycle=0;
    Int4 Len=LengthCMSA(1,cma);
    for(s = Len; s > 0; s--){
	if(Delete[s]){
	   cycle++;
	   Int4 end=s; 
	   while(Delete[s]){ total++; s--; }
	   s++;
#if 0	// ConvertColsToInsertsCMSA is not working!!!
	   rcma=0;
	   rcma=ConvertColsToInsertsCMSA(cma,1,s,end);
	   if(rcma){ NilCMSA(cma); cma=rcma; }
	   fprintf(stderr,"\n********************** cycle %d **************************\n",cycle);
	   char str[200];
           sprintf(str,"%s.cycle%d",NameCMSA(cma),cycle);
           // FILE *fp = open_file(str,".cma","w"); PutCMSA(fp,cma); fclose(fp);
#else	// use this instead...
	   // fprintf(stderr,"removing column %d to %d --> ",s,end);
	   // ColumnsToInsertCMSA(cma,s,end);
	   ColumnsToInsertCMSA2(cma,s,end);
	   // fprintf(stderr,"Len=%d.\n",LengthCMSA(1,cma));
#endif
	}
    } free(Delete); rcma=0;
#if 1
    if(IsDeletedCMSA(1,1,1,cma)) rcma=TrimBlkCMSA(cma,1,1,0,2);
    if(rcma){ NilCMSA(cma); cma=rcma; rcma=0; } 
    Len=LengthCMSA(1,cma); 
    if(IsDeletedCMSA(1,1,Len,cma)) rcma=TrimBlkCMSA(cma,1,0,1,2);
    if(rcma){ NilCMSA(cma); cma=rcma; rcma=0; } 
#endif
    return total;
}

FILE    *open_file2(const char *fstring,const char *subfile,const char *cmnd)
{
        FILE    *fptr;
        char    s[100];

        while(fstring[0] == ' ') fstring++;
        strcpy(s,fstring);
        strcat(s,subfile);
        if((fptr = fopen(s,cmnd)) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",s);
                print_error("File does not exist!\n");
        }
        return(fptr);
}

%}
%%
^>.+$		{ 
			long i,os,x;
			char *str1,*str2,*str3,*str4;
			unsigned char seq[5]={1,1,1,1,1,};
			long len=5;
			if(num==0){
			   fprintf(ofp,"[0_(1)=fa2cma(%d){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",
				TotalSeq);
			   fprintf(ofp,"(%d)",LenAln);
			   for(i=1; i<=LenAln; i++) fprintf(ofp,"*"); fprintf(ofp,"\n\n");
			} else {
			   if(MFASTA){ SeqCMA[LenCMA]='('; LenCMA++; }
			   fprintf(ofp,"$%d=%d(%d):\n",num,LenTrue,LenFake);
			   PutSeqInfo(ofp,DummySeq); SeqCMA[LenCMA]=0;
			   if(NumAln==LenAln && NumRF==0) fprintf(ofp,"{(%s()}*\n\n",SeqCMA);
			   else fprintf(ofp,"{(%s)}*\n\n",SeqCMA);
			} num++;
			SeqOffSet=0; SeqExtend=0;
			e_type  E=MkSeq(yytext, len, seq);
			if(DummySeq) NilSeq(DummySeq);
			DummySeq=E; Kingdom=KingdomSeq(E);
			if(PhylumSeq(E)) strncpy(Phylum,PhylumSeq(E),50);
			else strncpy(Phylum,"unknown",50);
			SeqOffSet=OffSetSeq(E);
			NumP=NumD=NumU=NumL=NumAln=NumGap=NumRF=0;
			LenTrue=LenFake=LenCMA=0;
		}

<MSF>{P}		{ 	NumP++; }

<MFASTA>{L}		{	// lower case...
			   char c=yytext[0];
			   if(NumAln==0){			// left flank
				SeqCMA[LenCMA]=')'; LenCMA++;
			   } SeqCMA[LenCMA]=toupper(c); LenCMA++; LenTrue++; LenFake++;
			   NumAln++;
			}

<MSF>{L}		{	// lower case...
			char c=yytext[0];
			if(NumAln==0){			// left flank
			   SeqCMA[LenCMA]=toupper(c); LenCMA++; LenTrue++; LenFake++; 
			}else if(NumAln==LenAln){	// right flank;
			   if(NumRF==0){ SeqCMA[LenCMA]='('; LenCMA++; } NumRF++;
			   SeqCMA[LenCMA]=toupper(c); LenCMA++; LenTrue++; LenFake++;
			} else {			// within alignment.
			   SeqCMA[LenCMA]=c; LenCMA++; LenTrue++;
			} NumL++; NumGap++;
		}

{D}		{	// delete character '-';
			if(NumAln==0){ SeqCMA[LenCMA]=')'; LenCMA++; }
			SeqCMA[LenCMA]='-'; LenCMA++; LenFake++;
#if 0
			if(NumAln==LenAln){	// start right flank;
			   if(NumRF==0){ SeqCMA[LenCMA]='('; LenCMA++; } 
			}
#endif
			NumD++; NumAln++;
		}

{U}		{	// upper
			if(NumAln==0){ SeqCMA[LenCMA]=')'; LenCMA++; }
			SeqCMA[LenCMA]=yytext[0]; LenCMA++; LenTrue++; LenFake++;
			NumU++; NumAln++;
		}

<MFASTA>{P}	{ print_error("invalid fasta input");  // util/FA2CMA/faToCMA.l
		}

[ \t\n]                  ;
%%
int	yywrap()
{
	fprintf(ofp,"$%d=%d(%d):\n",num,LenTrue,LenFake);
	PutSeqInfo(ofp,DummySeq); SeqCMA[LenCMA]=0;
	if(NumAln==LenAln && NumRF==0) fprintf(ofp,"{(%s()}*\n",SeqCMA);
	else fprintf(ofp,"{(%s)}*\n",SeqCMA);
	fprintf(ofp,"\n_0].\n\n");
	return(1);
}


#define USAGE_START "Usage: fa2cma <infile> [options]\n\
       Options: \n\
	 -cma2fa  convert <infile> from cma to fasta format\n\
          Note: this option eliminates extensions on ends\n\
       \n\n"

int     fa2cma_main(int argc, char *argv[],FILE *outfptr=0)
{
        char    last_c,c,str[300],cma2fa=0;
        long    i,j,len,arg,rlen,rLenAln;
	char	format='m'; // 'm' == MSF; 'f' == mfasta (NCBI CDD).
	cma_typ	cma=0;
	FILE	*tfp;
	char	debug=0; 

	BEGIN MSF;
	ofp=stdout;
        if(argc < 2) print_error(USAGE_START);
        for(arg = 1; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'x': break;
             case 'c': if(strcmp("-cma2fa",argv[arg]) == 0) cma2fa=1; break;
             default : print_error(USAGE_START);
           }
          } else if(arg > 1) print_error(USAGE_START);
        }
    TurnOffLicenseStatement();
    a_type AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
    if(cma2fa){
	cma=ReadCMSA2(argv[1],AB);
	if(!cma) print_error("cma file read error");
	PutFastaCMSA(stdout,cma);
	// PutFastaAlnCMSA(stdout,cma);
	if(cma) TotalNilCMSA(cma);
    } else {
	tfp=open_file2(argv[1],"","r");
	TotalSeq=0; LenAln=0; rLenAln=0;
	last_c='\n'; len=0; rlen=0;
	while((c=fgetc(tfp)) != EOF){
	   if(c=='>'){
		if(last_c != '\n') print_error("fasta file input error");
		if(TotalSeq==1){ LenAln=len; rLenAln=rlen; }
		else if(TotalSeq > 0 && rlen != rLenAln) print_error("fasta file input error");
		if(TotalSeq > 0 && len != LenAln){ BEGIN MFASTA; format='f'; }
		TotalSeq++; len=0; rlen=0;
		ungetc(c,tfp); assert(fgets(SeqCMA,YYLMAX,tfp) != NULL); 
	   } else if(isupper(c) || c == '-'){ len++; rlen++; }
	   else if(islower(c)){ rlen++; }
	   else if(last_c=='\n' && (c == '[' || c == '{')) print_error("fasta file input error");
	   last_c=c;
	} fclose(tfp);
#if 1	// 1. obtain data for consensus sequence...
	// if(format=='f') LenAln=rLenAln;
	UInt4	*ResCounts[256]; 
	for(c='A'; c <='Z'; c++) NEW(ResCounts[c], rLenAln +9, UInt4);
	tfp=open_file2(argv[1],"","r"); rlen=0;
	while((c=fgetc(tfp)) != EOF){
	   if(c=='>'){ rlen=0; ungetc(c,tfp); assert(fgets(SeqCMA,YYLMAX,tfp) != NULL); }
	   else if(isupper(c)){ rlen++; ResCounts[c][rlen]++; }
	   else if(c == '-'){ rlen++; ResCounts['X'][rlen]++; }
	   else if(islower(c)){ rlen++; c=toupper(c); ResCounts[c][rlen]++; }
	} fclose(tfp);

	// 2. compute consensus 
	UInt4	X,Z,r,rx,rz,max=0;
	char	dms_mode='F',wt_factor=1;
	Int4 	pernats=1000;
	UInt4	kount[30];
	dms_typ *dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	double	*BILD=0; NEW(BILD,rLenAln+9,double);
	char	maxC,*Csq; NEW(Csq, rLenAln+9, char);
	double	d,dd,D,RE=0,maxRE=0,sum,p,q=0.05;
	for(i=1; i<=rLenAln; i++) Csq[i]='-';
	for(j=0,i=1; i<=rLenAln; i++){
	  for(r=0; r <= nAlpha(AB); r++) kount[r]=0;
	  for(c='A'; c < 'Z'; c++){
	  	if((r=AlphaCode(c,AB))==0) continue;
	  	// fprintf(stderr,"%d.%c(%d) = %d.\n",i,c,r,ResCounts[c][i]);
		kount[r] += ResCounts[c][i];
	  } D=dms->bild(kount);
	  // if(D > 0) fprintf(stderr,"%d. bild=%.3f\n",i,D);
	  BILD[i]=D;
	  for(RE=0,max=0,maxC=0,c='A'; c < 'Z'; c++){
	      X=ResCounts[c][i];
	      p=(double)X/(double)TotalSeq;
	      // if(X > 0) RE += p*log(p/q);
	      if(c != 'X' && X > 0) RE += p*log(p/q);
	      if(X > max){ max=X; maxC=c; }
	  }
	  d=(double)ResCounts['X'][i]/(double)TotalSeq;
	  if(debug && maxC != 'X'){
		j++; fprintf(stderr,"%d(%c=%d): re=%.3f; bild=%.3f\n",j,maxC,ResCounts[maxC][i],RE,BILD[i]); 
	  }
	  if(d < 0.75 && maxC != 'X') Csq[i]=maxC;
	}
	for(c='A'; c <='Z'; c++) free(ResCounts[c]);
	delete dms; 

	double	*RunScore=0,AveScore=0;	// score in nats.
	NEW(RunScore,rLenAln+9,double);
	Int4 x,strt,end;
	h_type HG = Histogram("Running BILD scores (nats)",-50,200,10);
	for(i=1; !isalpha(Csq[i]); ) i++; strt=i;
	for(i=rLenAln; !isalpha(Csq[i]); ) i--; end=i;
	for(i=strt; i<=end; i++){
	   Int4 s=MAXIMUM(Int4,i-3,1),e=MINIMUM(Int4,i+3,rLenAln);
	   for(D=0,x=0,j=s; j <=e; j++){ D += BILD[j]; x++; }
	   RunScore[i] = D/(double)x;
	   if(debug) fprintf(stderr,"%d(%c): %.1f (%.2f)\n",i,Csq[i],RunScore[i],BILD[i]);
	   IncdHist(RunScore[i], HG);
	} 
	for(i=strt; RunScore[i] < 5.0 && i <= end; i++) Csq[i]='-';
	for(i=end; RunScore[i] < 5.0 && i >= strt;  i--) Csq[i]='-';
	// MeanHist(HG); VarianceHist(HG);
	if(debug) PutHist(stderr,60,HG); 
	NilHist(HG);
	free(BILD); free(RunScore);

	Int4 NumCsqAln=0;
	for(i=strt; i<=end; i++){ 
	   if(isalpha(Csq[i])) NumCsqAln++;
	   if(debug) fprintf(stderr,"%c",Csq[i]);
	} if(debug) fprintf(stderr,"\n num col=%d\n",NumCsqAln);

   if(NumCsqAln >= 50){
	yyin=tmpfile(); 
	fprintf(yyin,">consensus seq\n%s\n",Csq+1); TotalSeq++;
	FILE *fp=open_file2(argv[1],"","r");
	while((c=fgetc(fp)) != EOF) fprintf(yyin,"%c",c); 
	fprintf(yyin,"\n"); fclose(fp); 
	rewind(yyin);
	// fprintf(stderr,">consensus seq\n%s\n",Csq+1); 
	free(Csq); 
   } else {
	free(Csq); yyin=open_file2(argv[1],"","r");
   }
	ofp=tmpfile();
        while(yylex()); fclose(yyin);
#else
	ofp=tmpfile();
	yyin=open_file2(argv[1],"","r");
        while(yylex()); fclose(yyin);
#endif
	rewind(ofp); cma=ReadCMSA(ofp,AB);  fclose(ofp); 
	if(!cma) print_error("cma file read error");
	RmQueryGapsCMSA(cma); 
	if(outfptr) PutCMSA(outfptr,cma); else PutCMSA(stdout,cma); 
	TotalNilCMSA(cma);
    }
	NilAlpha(AB);
}

