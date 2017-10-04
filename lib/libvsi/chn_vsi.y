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

/* CHAIN analysis 'View Structural Interaction' yacc program */
%{
#include "chn_vsi.h"
/**************************** Global Variables ******************************/
int	yyparse();
int     yylex();
char	*FileName=0;
char	*PDB_File=0;
Int4	NUMBER_OF_PDB_FILES;
char	*PDB_File_Comment=0;
Int4	KEY_FILE=0;
char	KEY_CHAIN=0;
UInt4	PARSELINE=1;

Int4    TRACE_WIDTH=60;
Int4    WIRE_WIDTH=60;
Int4    THIN_WIRE_WIDTH=50;
Int4    SPACEFILL=100;
Int4    BIG_SPACEFILL=160;
BooLean	CIS_TRANS_PRO=FALSE;
BooLean	USE_CURRENT_DIR=FALSE;
BooLean	VERBOSE=FALSE;
BooLean	debug=FALSE;

extern int yyerror(char *s);

/****************************************************************************/
vsi_typ	HEAD_NODE=0;
vsi_typ	TAIL_NODE=0;
/********************* Declarations of Action Symbols ***********************/
%}
%start start
%token <cval>	RESIDUE			/*  residue */
%token <pcval>	MOLECULE		/*  residue */
%token <pcval>	FILENAME		/*  pdbfile */
%token <dummy>	PDBFILE			/*  pdbfile */
%token <cval>	LETTER			/*  letter */
%token <cval>	COLOR			/*  color */
%token <cval>	CHAIN			/*  chain */
%token <pcval>	ECHO_CMD		/*  echo into rasmol file */
%token <pcval>	COMMENT			/*  comment in file */
%token <aval>	AN_ATOM			/*  atom */
%token <aval>	HYDROGEN		/*  hydrogen atom */
%token <ival>	INTEGER			/*  integers */
%token <Vval>	VIEW			/*  file rotation, translation */
%token <dummy>	CARBONYL		// c=o group
%token <cval>	JUNK			/*  syntax error in lex program */
%token <dummy>	NEWLINE
%token <dummy>	END_OF_FILE		/* require program to come to parse entire file.*/
%token <dummy>	CLEAR
%union	{
		int	dummy;		/* no (null) attributes */
		char 	*pcval;		/* pointer to character string */
		char 	*mval;		/* pointer to molecule character string */
		char 	cval;		/* character */
		Int4	*Vval;		/* view of protein */
		char 	aval[9];	/* atom type*/
		Int4	ival;		/* integer type attribute */
		double	rval;		/* double (real) type attribute */
		vsi_typ	Lval;	// item list 
		char	Cval[3];	// color [0]='Y'; [1]=dots? (boolean)
		char	Sval[200];	// color [0]='Y'; [1]=dots? (boolean)
	}
%type	<Lval>	item_list trace item residue molecule echo view cmd
%type	<Cval>	color 
%type	<Sval>	junk fatal_error
%type	<pcval>	atom hydrogen 
%type	<cval>	backbone
%type	<dummy>	pdbfile end_file 
%%
start		: nl item_list end_file {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); YYACCEPT; 
		}
		| COMMENT nl item_list end_file {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); free($1); YYACCEPT;
		}
		| item_list end_file {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); YYACCEPT; 
		}
		| COMMENT nl fatal_error nl { yyerror($3); YYERROR; }
		| nl fatal_error nl { yyerror($2); YYERROR; }
		| fatal_error nl { yyerror($1); YYERROR; }
		| error nl end_file { yyerror("parse error"); YYERROR; }
		;

fatal_error	: item_list junk { strcpy($$,$2); }
		| junk { strcpy($$,$1); }
		;

item_list	: item_list item nl { 
			$$=$2; TAIL_NODE->next=$2; TAIL_NODE=$2; TAIL_NODE->comment=0;
		}
		| item_list item ',' { 
			$$=$2; TAIL_NODE->next=$2; TAIL_NODE=$2; TAIL_NODE->comment=0;
		}
		| item_list item COMMENT nl { 
			$$=$2; TAIL_NODE->next=$2; TAIL_NODE=$2; TAIL_NODE->comment=$3;
		}	// BELOW CONSTITUTE THE START OF THE LIST...
		| pdbfile nl { TAIL_NODE=HEAD_NODE=0; }
		| pdbfile nl item nl { $$=$3; TAIL_NODE=HEAD_NODE=$3; $$->comment=0; }
		| pdbfile nl item ',' { $$=$3; TAIL_NODE=HEAD_NODE=$3; $$->comment=0; }
		| pdbfile nl item COMMENT nl { $$=$3; TAIL_NODE=HEAD_NODE=$3; $$->comment=$4; }
		;

item		: residue { $$=$1; }
		| residue INTEGER { $$=$1; $1->item.residue.width=$2; }  // width == $2...
		| molecule { $$=$1; }
		| molecule INTEGER { $$=$1; $1->item.molecule.width=$2; }
		| cmd { $$=$1; }
		| trace { $$=$1; }
		| trace INTEGER { $$=$1; $1->item.trace.width=$2; }
		| view { $$=$1; }
		; 

cmd		: CLEAR { $$=MakeItemNode("clear.",'C'); }
		| echo { $$=$1; };

echo		: ECHO_CMD { $$=MakeItemNode($1,'E'); };

end_file	: END_OF_FILE ;

view		: VIEW { $$=MakeItemNodeView($1); } ;

trace		: INTEGER '-' INTEGER backbone COLOR 
		{ $$=MakeItemNodeTrace($1,$3,0,$4,$5,TRACE_WIDTH); }
		| INTEGER '-' INTEGER CHAIN backbone COLOR // 304-324A.M or 346-370A.R
		{ $$=MakeItemNodeTrace($1,$3,$4,$5,$6,TRACE_WIDTH); }
		;

backbone	: '^' { $$='^'; } 	// ribbon
		| '.' { $$='.'; } 	// trace
		| '+' { $$='+'; } 	// sticks
		| '~' { $$='~'; } 	// cartoon
		| ':' { $$=':'; } 	// spacefill
		| '@' { $$='@'; };	// dot clouds

residue		: RESIDUE INTEGER '.' color			// R127.Y or R127.{Y}
		{ $$=MakeItemNodeRes($1,$2,KEY_CHAIN,0,0,$4[0],$4[1],THIN_WIRE_WIDTH); } 
		| RESIDUE INTEGER CHAIN '.' color		// R127A.Y or R127A.{Y} 
		{ $$=MakeItemNodeRes($1,$2,$3,0,0,$5[0],$5[1],THIN_WIRE_WIDTH); } 
		| RESIDUE INTEGER '_' atom '.' color		// R127_nh2.X
		{ $$=MakeItemNodeRes($1,$2,KEY_CHAIN,$4,0,$6[0],$6[1],THIN_WIRE_WIDTH); } 
		| RESIDUE INTEGER CHAIN  '_' atom '.' color	// R127A_nh2.X
		{ $$=MakeItemNodeRes($1,$2,$3,$5,0,$7[0],$7[1],THIN_WIRE_WIDTH); } 

		| RESIDUE INTEGER '_' CARBONYL '.' color  
		{ $$=MakeItemNodeRes($1,$2,KEY_CHAIN,AllocString("c"),AllocString("o"),
						$6[0],$6[1],THIN_WIRE_WIDTH); } 
		| RESIDUE INTEGER CHAIN '_' CARBONYL '.' color  // c=o group
		{ $$=MakeItemNodeRes($1,$2,$3,AllocString("c"), AllocString("o"),
						$7[0],$7[1],THIN_WIRE_WIDTH); } 

		| RESIDUE INTEGER '_' atom '-' hydrogen '.' color  
						// R127A_nh2-2hh2.X or 126D295A_c-o.X
		{ $$=MakeItemNodeRes($1,$2,KEY_CHAIN,$4,$6,$8[0],$8[1],THIN_WIRE_WIDTH); } 
		| RESIDUE INTEGER CHAIN  '_' atom '-' hydrogen '.' color  
						// R127A_nh2-2hh2.X or 126D295A_c-o.X
		{ $$=MakeItemNodeRes($1,$2,$3,$5,$7,$9[0],$9[1],THIN_WIRE_WIDTH); } 
		; 

molecule	: MOLECULE INTEGER '.' color 		// !ATG801.C or !ATG801.{C} or ![SO4]1250.X
		{ $$=MakeItemNodeMol($1,$2,0,0,0,$4[0],$4[1],THIN_WIRE_WIDTH); }
		| MOLECULE INTEGER CHAIN '.' color 		// !ATG801A.C
		{ $$=MakeItemNodeMol($1,$2,$3,0,0,$5[0],$5[1],THIN_WIRE_WIDTH); }
		| MOLECULE INTEGER '_' atom '.' color		// !ATG801X_o2b.X
		{ $$=MakeItemNodeMol($1,$2,0,$4,0,$6[0],$6[1],THIN_WIRE_WIDTH); }
		| MOLECULE INTEGER CHAIN  '_' atom '.' color	// !ATG802X_o1g.X !ATG803X_o2*.X 
		{ $$=MakeItemNodeMol($1,$2,$3,$5,0,$7[0],$7[1],THIN_WIRE_WIDTH); }
		| MOLECULE INTEGER '_' atom '-' hydrogen '.' color  // 18HOH301_o-h1.X
		{ $$=MakeItemNodeMol($1,$2,0,$4,$6,$8[0],$8[1],THIN_WIRE_WIDTH); }
		| MOLECULE INTEGER CHAIN  '_' atom '-' hydrogen '.' color  
		{ $$=MakeItemNodeMol($1,$2,$3,$5,$7,$9[0],$9[1],THIN_WIRE_WIDTH); }
						// !ATG803X_n1-h.X 
		; 

atom		: AN_ATOM { $$=AllocString($1); }

hydrogen	: HYDROGEN { $$=AllocString($1); }

color		: '{' COLOR '}' { $$[0]=$2; $$[1]='T'; }
		| '(' COLOR ')' { $$[0]=$2; $$[1]='S'; }
		| COLOR { $$[0]=$1; $$[1]='F'; }
		;

junk		: junk JUNK { char s[3]; s[0]=$2; s[1]=0; strcpy($$,$1); strcat($$,s); 
			if( strlen($$) > 90) yyerror($$); }
		| JUNK { $$[0]=$1; $$[1]=0; $$[2]=0; }

pdbfile		: PDBFILE { PDB_File_Comment=0; }
		| PDBFILE COMMENT { PDB_File_Comment=$2; }
		;

nl		: nl COMMENT NEWLINE { free($2); }
		| nl NEWLINE { }
		| NEWLINE  { }
		;

%%

int	yyerror(char *s)
{ 
	char	str[102];
	if(FileName == 0){
		fprintf(stderr,"File number = %d\n",KEY_FILE);
		print_error("Invalid file number!");
	} else {
		fprintf(stderr,"\nfatal error (line %d): '%s'.\n",PARSELINE,s);
		fprintf(stderr,"%s --> \"%s\"\n",s,fgets(str,100,yyin));
		print_error("Input rejected.\n");
	}
	return 1;
}

void    Close(FILE *fptr) {if(fptr != stderr && fptr != stdout) fclose(fptr); }

#define USAGE_START "Usage: chn_vsi infile.vsi <int> [options]\n\
Usage: chn_vsi infile.vsi <int> <int> [options] (for multiple input files).\n\
       -A<int_array>  sort and print the files designated in <int_array>\n\
       -a<int>     print a copy of the current file adding <int> to residue sites\n\
       -append     append short rasmol add-on to standard output file\n\
       -Append     append full rasmol add-on to standard output file\n\
       -f<int>     print a copy of the current file adding file <int> to items\n\
       -c          show cis vs trans proline\n\
       -d<real>    dmax in Angstroms for classical H-bonds (default: 2.7 Angstoms)\n\
       -C          sort by chain and print a copy of the current file\n\
       -D          look for pdb input file within current directory\n\
       -B          Burley method for aromatic-aromatic interactions\n\
       -pml        generate a PyMol script instead\n\
       -pml=<file> generate a PyMol script instead\n\
       -s          sort and print a copy of the current file\n\
       -S          sort, purge, and print a copy of the current file\n\
       -skip=<char> if skip='W' the skip weak H-bonds; if skip == 'A' skip all H-bonds\n\
       -T          Test mode (fill in H-bonds)\n\
       -out=<name> Specify name of the output file\n\
       -t          Print protein regions w/o backbone traces (default: skip such regions)\n\
       -v          verbose mode\n\
       -w          print a copy of the current file\n\
       -x          dummy\n\n"


// int	main(int argc,char *argv[])
int	ChainVSI(int argc,char *argv[],FILE *ifp, FILE *ofp)
{ 
	Int4	i,n,arg,add=0,file=0,*value;
	char	c,mode=' ',filename[200]; filename[0]=0;
	vsi_typ	node,*headnode=0,*tailnode=0,tmp_node;
	item_type itm;
	char	**FileNames=0,*KeyChains=0,**FileComment=0;
	FILE	*fp;
	BooLean	only_if_trace=TRUE;
	Int4	append_output=0;
	char	modeHB=0,*pml_file=0,str[200];
        float   HA_dmax = 2.7,dmax=3.6;
	int	rtn_number=0;
	Int4	SARP_ID=0;	// the identifier for key file within a SARP multiple VSI input files.
	Int4	arg_start=3;

	if(argc < 3) print_error(USAGE_START);
	if(sscanf(argv[2],"%d",&i)!=1) print_error(USAGE_START);
	KEY_FILE=i;
	if(i <1) print_error(USAGE_START);
#if 0
        for(arg = 0; arg < argc; arg++) fprintf(stderr,"%s ",argv[arg]); 
	fprintf(stderr,"\n");
#endif
	if(argc >= 4 && argv[3][0] != '-'){
		if(sscanf(argv[3],"%d",&SARP_ID)!=1) print_error(USAGE_START); arg_start=4; 
		// fprintf(stderr,"sarp id = %d\n",SARP_ID);
	}
        for(arg = arg_start; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'c': CIS_TRANS_PRO=TRUE; break;
             case 'D': USE_CURRENT_DIR=TRUE; break;
             case 'a': 
		if(strcmp("-append",argv[arg]) == 0){
			append_output=1; 
			// print_error(" -append option not yet implemented (problems)");
		} else {
		  mode='w'; 
		  add=IntOption(argv[arg],'a',-5000,5000,USAGE_START); 
		} break;
             case 'f': mode='w'; file=IntOption(argv[arg],'f',1,5000,USAGE_START); break;
             case 'w': mode='w'; break;
             case 'd': HA_dmax=RealOption(argv[arg],'d',1.0,50,USAGE_START); break;
             case 's': 
		if(sscanf(argv[arg],"-skip=%c",&modeHB) != 1){
			if(argv[arg][2]==0) mode='s'; else print_error(USAGE_START);
		} break;
             case 'S': mode='S'; break;
             case 'C': mode='C'; break;
             case 'B': mode='B'; break;
             case 'T': mode='T'; break;
             case 't': only_if_trace=FALSE; break;
             case 'v': VERBOSE=TRUE; break;
             case 'p':  mode='p';
		if(sscanf(argv[arg],"-pml=%s",str) == 1){ 
			pml_file=AllocString(str);
		} else if(strcmp("-pml",argv[arg]) == 0){ 
			pml_file=0;
		} else print_error(USAGE_START);
		break;
             case 'o': mode=' ';
		if(sscanf(argv[arg],"-out=%s",&filename) != 1) print_error(USAGE_START);
		break;
             case 'A': 
	      if(strcmp("-Append",argv[arg]) == 0){
			append_output=2; 
	      } else {
		mode='A'; 
		Int4 tmp_value[1000];
		n=ParseIntegers(argv[arg]+2,tmp_value,USAGE_START);
		NEW(value,n+3, Int4);
		for(i=0; i<= n; i++) value[i]=tmp_value[i];
		NEW(headnode,n+3, vsi_typ);
		NEW(tailnode,n+3, vsi_typ);
		NEWP(FileNames,n+3, char);
		NEWP(FileComment,n+3, char);
		NEW(KeyChains,n+3, char);
	      } break;
             default : print_error(USAGE_START);
           }
	  }else print_error(USAGE_START);
	}
	if(mode != 'A'){
	   NUMBER_OF_PDB_FILES=0;
	   if(debug) fprintf(stderr,"=========== Parsing file %s ============\n",argv[1]);
	   if(ifp) yyin=ifp; else yyin = open_file(argv[1],"","r"); yyout=stderr; 
	   if(SARP_ID >0){	// then input contains multiple vsi files; go to starting file..
		char str[200];
		char	c; 
		Int4	i;
		do {
		  while((c=fgetc(yyin)) != '~'){
			if(c == EOF){
				fprintf(stderr,"Input does not contain a vsi file with ID %d\n",SARP_ID);
				print_error("Fatal error");
			}
		  } fgets(str,198,yyin);
		  if(sscanf(str,"$=%d.\n",&i)!=1) print_error(USAGE_START);
		} while(i != SARP_ID);
	   }
	   yyparse(); 
	   if(ifp==0) fclose(yyin); 
	   rtn_number=NUMBER_OF_PDB_FILES;
	   if(debug) fprintf(stderr,"=========== %d files parsed. ============\n",NUMBER_OF_PDB_FILES);
	}
	switch (mode){
	  case 'A':
	     file=0;
	     for(i=1; i <=n; i++){
		KEY_FILE=value[i]; 
	   	NUMBER_OF_PDB_FILES=0;
		yyin = open_file(argv[1],"","r"); yyparse(); fclose(yyin);
	   	rtn_number=NUMBER_OF_PDB_FILES;
		// fprintf(stderr,"value[%d]=%d (%c)\n",i,value[i],KEY_CHAIN);
		headnode[i]=HEAD_NODE;
		tailnode[i]=TAIL_NODE;
		FileNames[i]=FileName;
		FileComment[i]=PDB_File_Comment;
		KeyChains[i]=KEY_CHAIN;
		for(node=HEAD_NODE; node; node=node->next){
		    node->file=value[i];
		    itm=node->item;
		    if(node->type == 'R'){
			if(itm.residue.chain == 0) node->item.residue.chain=KEY_CHAIN; 
			// else fprintf(stderr,"res chain: %c\n",itm.residue.chain);
		    } else if(node->type == 'T'){
			if(itm.trace.chain == 0) node->item.trace.chain=KEY_CHAIN;
			// else fprintf(stderr,"trace chain: %c\n",itm.residue.chain);
		    }
		}
		HEAD_NODE=TAIL_NODE=0; FileName=0; PARSELINE=1;
	     }
	     HEAD_NODE=headnode[1];
	     for(i=1; i < n; i++){ tailnode[i]->next=headnode[i+1]; }
	     KEY_FILE=0;
	     HEAD_NODE=SortItems(HEAD_NODE);
	     mode='C';
	  case 'C':
	  case 'S': case 's':
	     if(mode=='S') HEAD_NODE=SortItems(HEAD_NODE,TRUE);
	     else if(mode=='C') HEAD_NODE=SortItemsByChains(HEAD_NODE);
	     else HEAD_NODE=SortItems(HEAD_NODE);
	  case 'w':
	     fp=open_file(argv[1],".crs","w");
	     if(FileNames){
	     	PrintItems(fp,HEAD_NODE,n,FileNames,FileComment,KeyChains,value,add,0,FALSE);
	     } else {
	     	PrintItems(fp,HEAD_NODE,FileName,PDB_File_Comment,KEY_CHAIN,KEY_FILE,
				add,file,FALSE);
	     }
	     Close(fp);
	     break;
	  case 'B':
	     {
//***************** Begin Burley Routine ****************************************
		char str[200],*str2,*DATADIR=0;
		BooLean *skip=0;

	    if(!USE_CURRENT_DIR){
        	DATADIR=getenv("CHN_RAS_DIR");
        	if(debug) fprintf(stderr,"CHN_RAS_DIR = %s\n",DATADIR);
		// if(DATADIR == 0) sprintf(str,"%s",FileName); else 
		sprintf(str,"%s/%s",DATADIR,FileName);
	    } else { sprintf(str,"%s",FileName); }
		str2=strstr(str,".pdb");
		*str2 = 0;
        	// fprintf(stderr,"pdb path = %s\n",str);
		pdb_typ P=MkPDB(str);
	        ErrorToDevNullPDB(P);

        Int4    C=GetChainNumberPDB(P,KEY_CHAIN);
        Int4    resStart=MinResPDB(C,P);
        Int4    resEnd=MaxResPDB(C,P),r;

		NEW(skip,resEnd+5,BooLean);
		for(r=0; r <= resEnd; r++) skip[r]=TRUE;
		for(node=HEAD_NODE; node; node=node->next){
		    node->file=KEY_FILE;
		    itm=node->item;
		    if(node->type == 'R'){
			if(itm.residue.chain == 0) node->item.residue.chain=KEY_CHAIN; 
			if(itm.residue.chain != KEY_CHAIN) continue;
			r=itm.residue.site;
			if(!(r >= resStart && r <= resEnd)){
			  fprintf(stderr,
			  "WARNING: Inconsistent input r=%d; resStart=%d; resEnd=%d\n",
					r,resStart,resEnd);
			  fprintf(stderr," (ignored)\n");
			} else skip[r]=FALSE;
			// else fprintf(stderr,"res chain: %c\n",itm.residue.chain);
		    }
		}
   {
        Int4    file_id=0;
        char    color=0;

	file_id=KEY_FILE;
        if(KEY_CHAIN==0) KEY_CHAIN='X';
#if 0
        PutAroAroInteractions(stdout,file_id,resStart,resEnd,HA_dmax,dmax,C,
                C,KEY_CHAIN,skip,color,P);
#endif
	free(skip);
   }
		NilPDB(P);
  }
//***************** End Burley Routine ****************************************
	     break;
//***************** Begin AutoFormating ***************************************
	  case 'T':
	{
	    char    str[200],*str2,*DATADIR=0;
	    Int4    N;
	    FILE    *fp;

	    if(!USE_CURRENT_DIR){
	      // if(FileName[0] == '/') DATADIR = AllocString(""); else DATADIR=getenv("CHN_RAS_DIR");
	      DATADIR=getenv("CHN_RAS_DIR");
	      if(debug) fprintf(stderr,"CHN_RAS_DIR = %s\n",DATADIR);
	      sprintf(str,"%s/%s",DATADIR,FileName);
	    } else sprintf(str,"%s",FileName);
#if 0
	    str2=strstr(str,".pdb"); *str2 = 0; 
	    // fprintf(stderr,"pdb path = %s\n",str);
	    pdb_typ P=MkPDB(str);
#else
	    pdb_typ P=MakePDB(str);
#endif
	    ErrorToDevNullPDB(P);
		
	    if(ofp) fp=tmpfile(); else
	    if(!USE_CURRENT_DIR){ fp=open_file(argv[1],"_tmp.vsi","w"); }
	    else { sprintf(str,"_tmp%d.vsi",KEY_FILE); fp=open_file(argv[1],str,"w"); }

	    if(KEY_CHAIN==0) fprintf(fp,"File%d=%s:X   \n\n",KEY_FILE,FileName);
	    else fprintf(fp,"File%d=%s:%c   \n\n",KEY_FILE,FileName,KEY_CHAIN);
	    vsi_typ first_head,first_tail=0;
	    first_head=PutAutoFormatVSI(fp,KEY_FILE,HEAD_NODE,KEY_CHAIN,
						HA_dmax,dmax, only_if_trace,P,modeHB);
	    if(ofp==0) fclose(fp); else { rewind(fp); }

	    for(N=0,node=first_head; node; node=node->next){ first_tail=node; N++; }
	    if(N <= 1) fprintf(stderr,"WARNING: ");
	    // fprintf(stderr,"%d items parsed\n",N);

	    if(PDB_File_Comment) free(PDB_File_Comment); free(FileName); 
	    FileName=0; PDB_File=0; PDB_File_Comment=0; KEY_CHAIN=0; PARSELINE=1; 
	    HEAD_NODE=TAIL_NODE=0; 
	    if(ofp){ yyin=fp; } else if(!USE_CURRENT_DIR){ yyin = open_file(argv[1],"_tmp.vsi","r"); }
	    else {
	        if(debug) fprintf(stderr,"=========== Parsing file %s_tmp%d.vsi ============\n",argv[1],KEY_FILE);
		sprintf(str,"_tmp%d.vsi",KEY_FILE); yyin = open_file(argv[1],str,"r"); 
	    }
	    NUMBER_OF_PDB_FILES=0;
	    yyout=stderr; yyparse(); fclose(yyin);	// don't close input file here.. 
	    if(debug) fprintf(stderr,"=========== %d files parsed. ============\n",NUMBER_OF_PDB_FILES);

	    for(N=0,node=HEAD_NODE; node; node=node->next){ node->file=KEY_FILE; N++; }
	    first_tail->next = HEAD_NODE; 
	    // HEAD_NODE=SortItems(first_head);		// sort nodes only...
	    HEAD_NODE=SortItems(first_head,TRUE);  // both sort & purge the nodes
	    if(N <= 1) fprintf(stderr,"WARNING: ");

	    // fprintf(stderr,"%d items parsed\n",N);
	    if(ofp) fp=ofp;
	    else if(!USE_CURRENT_DIR){
		fp=open_file(argv[1],".crs","w"); 
	        if(debug) fprintf(stderr,"=========== Creating file %s.crs ...============\n",argv[1]);
	    // } else { sprintf(str,"%d.vsi",KEY_FILE); fp=open_file(argv[1],str,"w"); }
	    } else { sprintf(str,".crs"); fp=open_file(argv[1],str,"w"); }
	    PrintItems(fp,HEAD_NODE,FileName,PDB_File_Comment,KEY_CHAIN,KEY_FILE,add,file,FALSE);
	    if(ofp == 0) Close(fp);

	    if(USE_CURRENT_DIR){
	        if(debug) fprintf(stderr,"======= Interaction network file '%s%d.vsi' created. =======\n",
			argv[1],KEY_FILE);
	    }
	    for(node=HEAD_NODE; node; ){ tmp_node=node->next; FreeNode(node);  node=tmp_node; }
	    if(PDB_File_Comment) free(PDB_File_Comment); free(FileName); NilPDB(P);
	}
	     break;
//******************* END AutoFormating ***************************************
//*****************************************************************************
	  case 'r':
	   fp=open_file(filename,"","w");
	   PrintItemsRasmol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
#if 1
	   if(append_output){
          	fprintf(fp,"select hetero and within(4.0,*%c)\n",KEY_CHAIN);
          	fprintf(fp,"wireframe 50\n");
          	fprintf(fp,"select water\ncolor cpk\n");
          	fprintf(fp,"select (ions,Mg,Zn,Mn,Fe,Co) and within(4.0,*%c)\n",
                                                                KEY_CHAIN);
          	fprintf(fp,"color cpk\nspacefill 250\n");
		if(append_output == 2){
          	  fprintf(fp,"select (protein or water or (ions,Mg,Zn,Mn,Fe,Co)) ");
          	  fprintf(fp,"and not within(4.0,*%c)\n",KEY_CHAIN);
          	  fprintf(fp,"spacefill off\nwireframe off\n");
		}
	   }
#endif
	   Close(fp);
	   fprintf(stderr,"======= Rasmol script file %s created. =======\n",filename);
	   break;
	  case ' ':
	  case 'p':
	  default:
	   if(filename[0]==0) sprintf(filename,"%s",argv[1]);
	   if(mode=='p'){
	     if(pml_file){
		fp=open_file(pml_file,".pml","w"); 
		fprintf(stderr,"=========== Creating file %s.pml ============\n",pml_file);
		free(pml_file); pml_file=0;
	     } else {
		fp=stdout;
		fprintf(stderr,"=========== Creating file %s.pml ============\n",argv[1]);
	     }
	     PrintItemsPyMol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
	   } else {	// rasmol script instead...
	     fp=open_file(filename,".ras","w");
	     fprintf(stderr,"=========== Creating file %s.ras ============\n",argv[1]);
	     PrintItemsRasmol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
	   } Close(fp); 
	   if(PDB_File_Comment) free(PDB_File_Comment); free(FileName);
	   for(node=HEAD_NODE; node; ){ tmp_node=node->next; FreeNode(node);  node=tmp_node; }
	   break;
	}
	return rtn_number;
}

