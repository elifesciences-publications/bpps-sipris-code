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
%START NODES NAMES IDENT SEQS 
I       (0|[1-9][0-9]*)
L       [A-Za-z]
C       {L}+[ ]+{L}+
H       {L}+[/]+{L}+
T	({C}|{H}|{L}+)+
S       [ \t]+
%{
#include "stdinc.h"
#include "alphabet.h"
#include "sequence.h"
#include "seqset.h"
#include "residues.h"
#include "txs_typ.h"
#include "txn_typ.h"
#include "n2a_typ.h"

#define MAX_NUM_GROUPS	500
#define MAX_SEQ_LENG	100000
#define NIL 		-1

#if defined (YYLMAX)
#undef YYLMAX
#endif

#define YYLMAX  100000

#define MAX_SEQ_ID	100000000
// 100,000,000 x 4 = 400 Mb
#define MAX_TAX_NODES	500000
#define MAX_TAX_NAMES	10

typedef struct {
	UInt4   tax_id;         // node id in GenBank taxonomy database
        char            *rank;          // superkingdom, kingdom, ...
        UInt4   parent;         // parent node id in GenBank taxonomy database.
        char            *name;         // names for node.
        unsigned char   gc_code;        // genetic code id.
} taxonomic_node_type;
typedef taxonomic_node_type *tnd_typ;

tnd_typ	tax_node;


UInt4	SeqId2Node[MAX_SEQ_ID];

BooLean	assign_pseudophyla_tax_node(UInt4 n,taxonomic_node_type *tnd)
{
	UInt4 p;
	// if(!strstr(tax_node[n].rank,"no rank")) return FALSE;
	if(tax_node[n].rank && 
		strstr(tax_node[n].rank,"kingdom") ||
		strstr(tax_node[n].rank,"phylum")) return FALSE;
	p=tax_node[n].parent;
	if(tax_node[p].rank && strstr(tax_node[p].rank,"kingdom")){
		if(tax_node[n].rank) free(tax_node[n].rank);
		tax_node[n].rank = AllocString("pseudophylum");
		return TRUE;
	} return FALSE;
}

char	find_kingdom_tax_node(UInt4 n,taxonomic_node_type *tnd)
{
	UInt4 p;
	if(tax_node[n].rank && strstr(tax_node[n].rank,"kingdom") &&
		strstr(tax_node[n].rank,"subkingdom") == NULL) 
		return tax_node[n].name[0];
	while((p=tax_node[n].parent) > 1){
		n = p;
		if(tax_node[n].rank && 
			strstr(tax_node[n].rank,"kingdom")) 
				return tax_node[n].name[0];
	} return 0;
}

char	*find_phylum_tax_node(UInt4 n,taxonomic_node_type *tnd)
{
	UInt4 p;
	if(tax_node[n].rank && strstr(tax_node[n].rank,"phylum")) 
		return tax_node[n].name;
	while((p=tax_node[n].parent) > 1){
		n = p;
		if(tax_node[n].rank && 
			strstr(tax_node[n].rank,"phylum")) return tax_node[n].name;
	} return 0;
}

void	print_tax_node(FILE *fp,taxonomic_node_type tnd)
{
	fprintf(fp,"node: %d | %s | %d | %d |",
		tnd.tax_id,tnd.rank,tnd.parent,tnd.gc_code);
	if(tnd.name){
		fprintf(fp," %s ",tnd.name);
	} else fprintf(fp," no name!");
	fprintf(fp,"\n");
}

void	prune_othertax_nodes(taxonomic_node_type *tnd)
// These seem to be viruses and are not of interest in CHAIN analysis.
{
	for(UInt4 n=2; n <= MAX_TAX_NODES; n++){
	   if(tax_node[n].rank){
		UInt4 p = tax_node[n].parent;
		if(p==1 && strcmp(tax_node[n].name,"other")==0 ||
		     strstr(tax_node[n].name,"environmental samples") ||
		     strstr(tax_node[n].name,"uncultured bacteria") ||
		     strstr(tax_node[n].name,"candidate division")){
			free(tax_node[n].name); tax_node[n].name=0;
			free(tax_node[n].rank); tax_node[n].rank=0;
		} else if(strstr(tax_node[n].rank,"no rank")){
			if(strstr(tax_node[n].name,"unclassified ")){
				free(tax_node[n].name); tax_node[n].name=0;
				free(tax_node[n].rank); tax_node[n].rank=0;
			} else if(tax_node[p].rank && 
					strstr(tax_node[p].rank,"phylum")){
			   free(tax_node[n].name); tax_node[n].name=0;
			   free(tax_node[n].rank); tax_node[n].rank=0;
			} else if(p ==1 && 	// root node as parent...
			       (strstr(tax_node[n].name,"Viruses") ||
				strstr(tax_node[n].name,"Viroids") ||
				strstr(tax_node[n].name,"unclassified"))){
			   free(tax_node[n].name); tax_node[n].name=0;
			   free(tax_node[n].rank); tax_node[n].rank=0;
			   // print_tax_node(stderr,tax_node[n]);
			} else if(find_phylum_tax_node(n,tnd)){
			   free(tax_node[n].name); tax_node[n].name=0;
			   free(tax_node[n].rank); tax_node[n].rank=0;
			}
		} else if(strstr(tax_node[n].rank,"subphylum") &&
			tax_node[p].rank && strstr(tax_node[p].rank,"phylum")){
			   // subphylum of a known phylum.
			   free(tax_node[n].name); tax_node[n].name=0;
			   free(tax_node[n].rank); tax_node[n].rank=0;
		}
	   }
	}
}

Int4	prune_tax_nodes(taxonomic_node_type *tnd)
// Remove names from all nodes pointing to nodes with rank=0;
{
	// n == 1 is root node.
	Int4	total=0;
	for(UInt4 n=2; n <= MAX_TAX_NODES; n++){
	  if(tax_node[n].rank){
		UInt4 p = tax_node[n].parent;
		if(tax_node[p].rank==0){
			total++;
			free(tax_node[n].name); tax_node[n].name=0;
			free(tax_node[n].rank); tax_node[n].rank=0;
		}
	  }
	} return total;
}

// NEW structures for direct fasta indexing without 'gi_taxid_prot.dmp' file.
UInt4	*NumNodesHash[256];	// bos taurus == NumNodesHash['b']['o'];
tnd_typ		*NodesHash[256];	// has table for nodes
// End structures for direct fasta indexing.
Int4 GetTaxNodeID(char *string)

{
	char	*ptr;
	Int4	i;
	ptr=strstr(string,"["); 
	if(ptr==0) return 0;
	else {
		if(NumNodesHash[ptr[1]] ==0){
			NEW(NumNodesHash[ptr[1]],257,UInt4);
		}
		NumNodesHash[ptr[1]][ptr[2]]++;
	}
	return 1;
}

a_type	aaAB,dnaAB,AB;
BooLean	ISDNA=FALSE;
BooLean	PRINTSEQ=FALSE;
BooLean	SKIPCHORDATES=FALSE;
BooLean	SKIPUNKNOWN=TRUE;

Int4	Field=0;
Int4	File=0;
UInt4	Node=0;
UInt4	LineNum=1;
UInt4	SeqId=0;


%}
%%

^">".*$		{
		   BEGIN SEQS;
		   UInt4 i,c,n;
		   if(sscanf(yytext,">gi|%u|",&i) != 1){
				print_error("fasta file lacks gi numbers!");
		   }
		   assert(i < MAX_SEQ_ID);

		   n = SeqId2Node[i];
		   assert(n < MAX_TAX_NODES);

		   if(GetTaxNodeID(yytext) == 0) { // ONly a stub for now.
			fprintf(stderr,"organism not found for gi=%d!\n",i);
		   }

		   char *phylum=find_phylum_tax_node(n,tax_node);
		   if(SKIPUNKNOWN && phylum==0 || 
			 phylum && SKIPCHORDATES && strstr(phylum,"Chordata")){
			PRINTSEQ=FALSE;
		   } else {
			PRINTSEQ=TRUE;
			char kingdom=find_kingdom_tax_node(n,tax_node);
			for(c=0; !isspace(yytext[c]); c++){
				fprintf(stdout,"%c",yytext[c]);		
			} c++;
			if(kingdom==0) kingdom='X';
			if(phylum==0) { 
			  if(ISDNA){
			    fprintf(stdout," {<unknown(%c1)>}%s",kingdom,yytext+c);
			  } else {
			    fprintf(stdout," {<unknown(%c)>}%s",kingdom,yytext+c);
			  }
			} else { 
			  if(ISDNA){
			    fprintf(stdout," {<%s(%c%d)>}%s",
					phylum,kingdom,tax_node[n].gc_code,
					yytext+c);
			  } else {
			    fprintf(stdout," {<%s(%c)>}%s",phylum,kingdom,yytext+c);
			  }
			} 
		   }
		}

<SEQS>.		{ if(PRINTSEQ) ECHO; }

<SEQS>[\n]	{ if(PRINTSEQ) ECHO; }

{I}		{
		   UInt4 i;
		   sscanf(yytext,"%u",&i);
		   if(File==0){
			 BEGIN NODES;
			 if(Field == 0){ 
				Node=i;
				tax_node[Node].tax_id=i;
			   	tax_node[Node].name=0;
				tax_node[Node].rank=0;
				assert(Node > 0 && Node < MAX_TAX_NODES);
				// printf("node %d",Node);
			 } else {
			  assert(Node > 0 && Node < MAX_TAX_NODES);
			  if(Field == 1){
				tax_node[Node].parent =i;
				// printf("--> %d",i);
			  } else if(Field == 6){
				tax_node[Node].gc_code=i;
				// printf("gc: %d",i);
			  } else ; // ECHO;
			 }
		   } else if(File==1){
			   BEGIN NAMES;
			   if(Field == 0){ 
				Node=i;
				assert(Node < MAX_TAX_NODES);
				// printf("Node %d",Node);
			   }
		   } else if(File==2){
			 BEGIN IDENT;
			 if(Field == 0){
				SeqId = i;
				assert(SeqId < MAX_SEQ_ID);
				Field++;
			 } else if(Field == 1){
				assert(i < MAX_TAX_NODES);
				SeqId2Node[SeqId] = i;
			 }
		   } else print_error("fatax input error!");
		}

<NODES>[|]		{ 
				Field++; // ECHO; 
			}

<NODES>{L}+	|
<NODES>{C}	{
			if(Field != 2) ; //ECHO;
			else if(strstr(yytext,"phylum") || 
					strstr(yytext,"class") || 
					strstr(yytext,"order") || 
					strstr(yytext,"no rank") || 
					strstr(yytext,"kingdom")){
				tax_node[Node].rank = AllocString(yytext);
				// print_tax_node(stdout,tax_node[Node]);
			} else if(Node==136) fprintf(stderr,"-->%s\n",yytext);
		}

<NAMES>[|]		{
				Field++; // ECHO; 
			}

<NAMES>{T}+	{
		  if(Field == 1){
		    if(tax_node[Node].rank){	// phylum/kingdom/etc.
			   if(tax_node[Node].name == 0)
			     tax_node[Node].name = AllocString(yytext);
			   // print_tax_node(stderr,tax_node[Node]);
		    }
			else {
			   if(tax_node[Node].name == 0)
				tax_node[Node].name = AllocString(yytext);
			}
		  }
		}

<NAMES>.	{ }

<IDENT>[\n]	|
<NODES>[\n]	|
<NAMES>[\n]	{ 
			Field=0; 
			LineNum++;
			// ECHO;
		}

{S}		{
		}

.		{
			fprintf(stderr,"Node = %d; line = %d\n",Node,LineNum);
			fprintf(stderr,"PARSE ERROR!: %s\n",yytext);
		}

%%

int	yywrap()
{
    if(File==1){
	prune_othertax_nodes(tax_node);
	while(prune_tax_nodes(tax_node) > 0);
	for(Int4 n=2; n <= MAX_TAX_NODES; n++){
		if(tax_node[n].rank){
			assign_pseudophyla_tax_node(n,tax_node);
			// print_tax_node(stdout,tax_node[n]);
			if(strstr(tax_node[n].rank,"phylum") ||
				strstr(tax_node[n].rank,"kingdom")){
				print_tax_node(stderr,tax_node[n]);
			}
		}
	}
    } return 1;
}

#define DNA_MTRX "-4  -4  -4  -4  -4 \
                  -4   5  -4  -4  -4 \
                  -4  -4   5  -4  -4 \
                  -4  -4  -4   5  -4 \
                  -4  -4  -4  -4   5 "

#if 0	//*************************************************************
names.dmp
nodes.dmp
gi_taxid_nucl.dmp
gi_taxid_prot.dmp

delnodes.dmp
division.dmp
gc.prt
gencode.dmp
gi_taxid.readme
merged.dmp
#endif	//*************************************************************

#define USAGE_START "Usage: nrtax fafile [options]\n\
       NOTE: for very large files use 'cat file | fatax [options]' \n\
       -n          nucleotide sequences in input\n\
       -C          Skip Chordate sequence in input\n\
       -u          Don't skip sequences from unknown phylums\n\
       -x          dummy\n\n"

int	main(int argc, char *argv[])
{
	char	str[300];
	Int4	i,j,arg;

	if(argc < 1) print_error(USAGE_START);
	for(arg = 1; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'n': ISDNA=TRUE; break;
             case 'C': SKIPCHORDATES=TRUE; break;
             case 'u': SKIPUNKNOWN=FALSE; break;
             case 'x': break;
             default : print_error(USAGE_START);
           }
          } else if(arg > 1) print_error(USAGE_START);
	}
	if(getenv("TAXDUMPDIR") == 0){
		fprintf(stderr,"Define the path to NCBI taxdump directory ");
		fprintf(stderr,"via the environmental variable 'TAXDUMPDIR'.\n");
		print_error("   fatal error!");
	}
	print_error("This program not finished (problem with nr organism names)!");
	aaAB = MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	dnaAB = MkAlpha("NACGT",DNA_MTRX);
	for(i=0; i<256; i++) NumNodesHash[i]=0;

	NEW(tax_node,MAX_TAX_NODES+2,taxonomic_node_type);

	File=0; LineNum=0;
	strcpy(str,getenv("TAXDUMPDIR"));
	strcat(str,"/nodes.dmp");
	yyin=open_file(str,"","r");
        while(yylex()); fclose(yyin);

	File=1; LineNum=0;
	strcpy(str,getenv("TAXDUMPDIR"));
	// strcat(str,"/names.dmp");
	strcat(str,"/scientific.dmp");	// == grep 'scientific name' names.dmp
	yyin=open_file(str,"","r");
        while(yylex()); fclose(yyin);

	File=2; LineNum=0;
	strcpy(str,getenv("TAXDUMPDIR"));
	if(ISDNA) strcat(str,"/gi_taxid_nucl.dmp");
	else strcat(str,"/gi_taxid_prot.dmp");
	for(Int4 s=0; s < MAX_SEQ_ID; s++) SeqId2Node[s]=0;
	yyin=open_file(str,"","r");
        while(yylex()); fclose(yyin);

	File=3;
	if(argc == 1 || argv[1][0] == '-') yyin=stdin;
	else yyin=open_file(argv[1],"","r");
        while(yylex()); fclose(yyin);
	
	for(i=0; i<256; i++) {
	  if(NumNodesHash[i]){
		for(j=0; j<256; j++) {
	  	  if(NumNodesHash[i][j]){
			fprintf(stderr,"%c%c: %d\n",(char)i,(char)j,NumNodesHash[i][j]);
		  }
		}
	  }
	}

	NilAlpha(aaAB); NilAlpha(dnaAB);
	return 0;
}

