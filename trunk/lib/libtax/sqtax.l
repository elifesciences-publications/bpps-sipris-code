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

%p 2000
%p 2000
%e 2000
%n 2000
%k 2000
%a 2000
%o 2000
%{
#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "seqset.h"
#include "residues.h"
#include "txs_typ.h"
#include "n2a_typ.h"

#define MAX_NUM_SEQS	100000
#define MAX_NUM_GROUPS	500
#define MAX_SEQ_LENG	100000
#define NIL 		-1
#define PID		1
#define LOCUS		2
#define SOURCE		3
#define DEFINITION	4
#define ORIGIN		5

#if defined (YYLMAX)
#undef YYLMAX
#endif

#define YYLMAX  5000

a_type	aaAB,dnaAB,AB;
int SECTION=NIL;

char    *inputfile;
e_type  QuerySeq=0;

Int4	SeqId[MAX_NUM_SEQS];
Int4	GroupId[MAX_NUM_SEQS];
e_type	Seq[MAX_NUM_SEQS];
UInt4	Entry=0;
char	locus[YYLMAX];
char	descript[YYLMAX];
char	pid[YYLMAX];
char	defline[YYLMAX];
unsigned char seq[MAX_SEQ_LENG];
Int4	no_res=0;
n2a_typ	*n2a=0;
BooLean	IsDNA=FALSE;
BooLean	GotResidues=FALSE,GotTaxa=FALSE,GotGI=FALSE;
txc_typ	*txc=0;

%}
%%
^[\\]+			{ SECTION = NIL; printf("\n"); }

^LOCUS[ \t]+.*$ 	{ 
			  double score;
			  e_type dnaE;
			  if(Entry > 0){	// Process previous entry.
			  	assert(GotResidues); GotResidues=FALSE;
			  	assert(GotTaxa); GotTaxa=FALSE;
			  	assert(GotGI); GotGI=FALSE;
				// sprintf(defline,"gi|%s|%s %s",pid,locus,descript);
				sprintf(defline,"gi|%d|%s %s",SeqId[Entry],locus,descript);
				if(IsDNA){
				   dnaE= MkSeq(defline,no_res,seq);
				   Seq[Entry] = n2a->BestTranslate(QuerySeq,&score,dnaE);
				   if(!Seq[Entry]){ // then delete this hit.
					PutSeq(stderr,dnaE,dnaAB);
					txc->GroupHits[GroupId[Entry]]--; 
					GroupId[Entry]=0; Entry--; 
					fprintf(stderr,"Translation failure (score = %.2f)\n",
						score); 
				   } NilSeq(dnaE);
				} else {
				   Seq[Entry] = MkSeq(defline,no_res,seq);
				   assert(Seq[Entry]);
				}
			  }
			  no_res=0; 
			  IsDNA=FALSE; AB=aaAB;
			  SECTION = LOCUS; Entry++; // sprintf(pid,"0");
			  sscanf(yytext,"LOCUS %[^ ]",locus);
			  // printf("%s\n", locus);
			}

^SOURCE[ \t]+.*$ 	{
				assert(!GotTaxa); SECTION = SOURCE; GotTaxa=TRUE;
			}

^PID[ \t]+.*$ 		{ 	// Present in proteins files only.
			   {
				Int4	sqid;
				assert(!GotGI);
				// sscanf(yytext,"PID g%[0-9]",pid);
				if(sscanf(yytext,"PID         g%d",&sqid) == 1){
					fprintf(stderr,"seq id = %d\n",sqid);
					SeqId[Entry]=sqid;
				} else {
					fprintf(stderr,"WARNING: sscanf failed!!!");
					SeqId[Entry]=0;
				} GroupId[Entry]=0;
				GotGI=TRUE;
			   }
			}

^VERSION[ \t]+.*$	{
				Int4	sqid,oid=0;
				if(GotGI) { oid=SeqId[Entry]; }
				else GotGI=TRUE;
				char *giptr = strstr(yytext,"GI:");
				// assert(giptr);
				if(giptr && sscanf(giptr,"GI:%d",&sqid) == 1){
					fprintf(stderr,"seq id = %d\n",sqid);
					SeqId[Entry]=sqid;
					if(oid) assert(oid==sqid);
				} else {
					fprintf(stderr,"WARNING: sscanf failed!!!");
					SeqId[Entry]=0;
				} GroupId[Entry]=0;
			}

^DEFINITION[ \t]+.*$ 	{
			  strcpy(descript,yytext+12);
			  // sscanf(yytext,"DEFINITION %s",descript);
				fprintf(stderr,"%s\n",descript);
			}

^ACCESSION[ \t]+.*$	{ } 

^REFERENCE.*$		{ SECTION = NIL; }

^COMMENT.*$		{ SECTION = NIL; }

^FEATURES.*$		{ SECTION = NIL; }

^BASE[ ]COUNT.*$	{	// Nucleotide sequence... // 141 a    107 c    146 g    101 t
				Int4	a,c,g,t;
				Int4	min_len=25,code=1;
				if(sscanf(yytext,
				   "BASE COUNT %d a %d c %d g %d t",&a,&c,&g,&t) == 4){
					fprintf(stderr,"Nucleotide sequence: %d a %d c %d g %d t\n",
						a,c,g,t);
				}
				if(!n2a) { n2a = new n2a_typ(code,min_len,aaAB,dnaAB); }
				IsDNA=TRUE; AB=dnaAB;
			}

^ORIGIN.*$		{
				SECTION = ORIGIN; no_res=0; 
			  	assert(GotTaxa);
			  	assert(!GotResidues);
				if(GroupId[Entry] ==0){
					// fprintf(stderr,"gi|%s not found!\n",pid);
					fprintf(stderr,"gi|%d not found!\n",SeqId[Entry]);
				}
			}

[a-zA-Z]                { if(SECTION == ORIGIN) {
                                no_res++;
				seq[no_res]=AlphaCode(yytext[0],AB);
				// ECHO;
			  }
                        }

[ \t\n]			{ if(SECTION == SOURCE) { } }

^[/][/]			{ SECTION = NIL; }

.*$			{ 
			  Int4 k,s,g;
			  const char *name;
			  if(SECTION == ORIGIN) {
				for(s = 0; s <= yyleng; s++){
				  char c=yytext[s];
				  if(isalpha(c)){
					no_res++;
					seq[no_res]=AlphaCode(c,AB);
				  }
				}
				GotResidues=TRUE;
			  } else if(SECTION == SOURCE) { 
			  	assert(!GotResidues);
				for(g=1; g <= txc->GroupNum; g++){
				  if(!txc->Group[g]) continue;  // deleted group
				  for(s=0; txc->Group[g][s]; s++){
					// Check all names for a match...
					name = txc->Group[g][s];
					if(strstr(yytext,name)){  // match found!
					   fprintf(stderr,"GROUP: %s\n",name);
					   if(GroupId[Entry]==0){ 
						GroupId[Entry]=g; 
						txc->GroupHits[g]++;
					   } else if(g < GroupId[Entry]){
						// if find a higher number group 
						// then use it instead.
						fprintf(stderr,"g=%d < GroupId[%d]=%d\n",
							g,Entry,GroupId[Entry]);
						txc->GroupHits[GroupId[Entry]]--; 
						GroupId[Entry]=g; 
						txc->GroupHits[g]++;
					   } break;
					}
				  }
				} 
			  } 
			}

%%

int	yywrap()
{
	Int4	g,e;
	e_type	dnaE;
	double	score;

	if(Entry > 0){
	     assert(GotTaxa); assert(GotResidues);
	     sprintf(defline,"gi|%s|%s %s",pid,locus,descript);
	     if(IsDNA){
		dnaE=MkSeq(defline,no_res,seq);
		Seq[Entry] = n2a->BestTranslate(QuerySeq,&score,dnaE);
		if(!Seq[Entry]){ 
			PutSeq(stderr,dnaE,dnaAB);
			txc->GroupHits[GroupId[Entry]]--; 
			GroupId[Entry]=0; Entry--; 
			fprintf(stderr,"Translation failure (score = %.2f)\n",
				score); 
		} NilSeq(dnaE);
	     } else Seq[Entry] = MkSeq(defline,no_res,seq);
	}

	e_type		*Seq2;
	Int4		s,G,*GrpSize,QueryGroupPurged=0;
	const char      **GrpName;
	char		*KingSymbol;

	NEW(Seq2,Entry+3,e_type);
	NEW(GrpSize,txc->GroupNum+3,Int4);
	NEWP(GrpName,txc->GroupNum+3,const char);
	NEW(KingSymbol,txc->GroupNum+3,char);

	Int4	QueryGroup=0;
	for(e=1; e <= Entry; e++){
		if(GroupId[e] == 0){
			fprintf(stderr,"Taxonomic Group unknown:\n");
			PutSeq(stderr,Seq[e],aaAB);
		}
	}
	// remove subseqs from Query group.
	s=1; Seq2[s]=QuerySeq; // first sequence in query group.
	for(G=0,e=1; e <= Entry; e++){
	   if(IsSubSeq(QuerySeq,Seq[e])){ // Found query group!
		PutSeq(stderr,Seq[e],aaAB); 
		fprintf(stderr,"GroupId[%d]=%d!\n",e,GroupId[e]);
		if(GroupId[e] == 0){ continue; } // query taxonomic group not found.
		QueryGroup=g=GroupId[e];
		G++; GrpName[G]=txc->Group[g][0];  // G == 1 is query group.
		KingSymbol[G]=txc->Kingdom[g];
		GrpSize[G]=txc->GroupHits[g];      // get query group size.
		// GrpSize[G] = 1;	// replace e with query itself...
		for(Int4 f=1; f <= Entry; f++){
		   if(f == e) continue; // replaced e with query.
		   if(GroupId[f]==QueryGroup){  // Find other sequences in query group.
	   		if(!IsSubSeq(QuerySeq,Seq[f])){
			   s++; Seq2[s]=Seq[f]; 
			   // GrpSize[G]++;  
			} else {	// if subseq of query or vice versa...
			   QueryGroupPurged++;
			   GrpSize[G]--;  
			   fprintf(stderr,"Purged a query-related sequence.\n");
			}
		   }
		} break;
	   }
	} 
	// assert(GrpSize[G] > 0);
	if(QueryGroup==0){
#if 1	// Fix problem when QueryGroup not found...
	    G++; GrpName[G]=txc->Group[0][0];
	    KingSymbol[G]=txc->Kingdom[0];
	    GrpSize[G]=1;
	    // Seq2[1] is the query group sequence.
#endif
	    fprintf(stderr,"WARNING: QueryGroup not found!\n");
	} else assert(GrpSize[G] == s);

	for(g=1; g <= txc->GroupNum; g++){
	   if(g == QueryGroup) continue; // took care of this already.
	   if(txc->Group[g] && txc->GroupHits[g] > 0){
		G++; GrpName[G]=txc->Group[g][0];
		KingSymbol[G]=txc->Kingdom[g];
		GrpSize[G]=txc->GroupHits[g];
		Int4 gs=0;
		for(e=1; e <= Entry; e++){
		   if(GroupId[e]==g){ gs++; s++; Seq2[s]=Seq[e]; }
		} 
		// GrpSize[G]=gs;
		assert(gs == txc->GroupHits[g]);
	   }
	}
	Int4 NumEntries=s;
for(s=0,g=1; g <= G; g++){
	fprintf(stdout,"%s (%c%d)\n",GrpName[g],KingSymbol[g],GrpSize[g]);
	s += GrpSize[g];
} assert(s==NumEntries);
	{ 
// assert(!"debug\n");
	  txs_typ txs(inputfile,GrpSize,G,GrpName,KingSymbol,NumEntries,Seq2,aaAB);
	  // Purge subseqs from each group before outputting...
	  // subSeq: Let X in a subseq match with anything.
	  txs.RmSubseq( ); // remove redundant entries. 
	  // txs.WriteFile("junk");
	  // txs_typ txs2("junk.txs",aaAB);
	  txs.WriteFile(inputfile);
	  FILE *sfp=open_file(inputfile,".tsq","w");
	  txs.PutSeqs(sfp); fclose(sfp);
	} free(Seq2); free(GrpSize); free(GrpName);
	if(n2a) delete n2a;
	if(Entry > NumEntries) 
	   fprintf(stderr,"Unable to classify %d sequences\n",
		(Entry-QueryGroupPurged)-NumEntries);
	return(1);
}

#define DNA_MTRX "-4  -4  -4  -4  -4 \
                  -4   5  -4  -4  -4 \
                  -4  -4   5  -4  -4 \
                  -4  -4  -4   5  -4 \
                  -4  -4  -4  -4   5 "

int	main(int argc, char *argv[])
{

	if(argc != 3) print_error("USAGE: sqtax keyfafile input ");
	aaAB = MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	dnaAB = MkAlpha("NACGT",DNA_MTRX);
	inputfile=argv[1];
	QuerySeq=ReadSeqFA(argv[1],1,aaAB);
	txc = new txc_typ(); // initialize(  );
	yyin=open_file(argv[2],"","r");
        while(yylex()); 
	fclose(yyin);
	NilAlpha(aaAB);
	NilAlpha(dnaAB);
	delete txc;
	return 0;
}

