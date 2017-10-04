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

#include "hpt_typ.h"
#include "blosum62.h"

#define	HPT_PUBLIC_USAGE2 "\n\
\n\
************************** hyperpartition file syntax: **************************\n\
Hyperpartition:\n\
<partitions> 1.<name><type>\n\
<partitions> 2.<name><type>\n\
   :    :      :   :   :\n\
   :    :      :   :   :\n\
<partitions> <int1>.Random=<int2>.   (<int1> = # sets. <int2> = # random sequences to be generated; default=10000)\n\
\n\
Settings:                            (this section is optional)\n\
1.<category_name> <options>\n\
2.<category_name> <options>\n\
   :    :      :   :\n\
   :    :      :   :\n\
<int>.<category_name> <options>         (where <int> corresponds to the total # categories)\n\
\n\
*************************** non-terminal definitions: ***************************\n\
  <partitions> = a string of '+','-' & 'o' characters representing the partitions to which the set belongs\n\
	'+'  --> the set belongs to the foreground partition\n\
	'-'  --> the set belongs to the background partition\n\
     	'o'  --> the set belongs to neither partition (i.e., set omitted)\n\
                 (Note that there must be one such character for each category)\n\
                 For example, <partitions> = \"+----++++++--oo-o\" (17 categories)\n\
  <name>  =  name of the corresponding seed alignment within the *.dma file\n\
  <type>  =  '.'  --> set is specific (e.g., <name><type> = 'Ran.') or\n\
             '?'  --> set is generic (e.g., <name><type> = 'MiscRasLike?')\n\
\n\
  <category_name> = an arbitrary character string that describes the category\n\
  <options> =\n\
   -quality=<int>       - Quality of the category's foreground sequences (range=0-9; default=2)\n\
				higher --> sequences more stringently match the category's foreground pattern\n\
   -noise=<int>         - Degree to which the average pattern position is contaminated (range=0-50\%; default=10\%)\n\
				higher --> more noise (pattern is less stringently conserved)\n\
   -weight=<int>        - Weight of the prior versus the data when computing the noise (range=2..1000; default=10)\n\
				higher --> more empirical evidence required to alter the noise parameter\n\
   -contrast=<int>      - Number of pattern positions highlighted in output alignments (value > 1; default=10)\n\
				higher --> more positions highlighted (i.e., LOWER CONTRAST)\n\
				lower  --> fewer positions highlighted (i.e., HIGHER CONTRAST)\n\
\n\
\n"

void	MultiPutHpt(FILE *fp, hpt_typ *head)
{
	for(hpt_typ *hpt = head; hpt; hpt=hpt->RtnNext()){
	    hpt->Put(fp); fprintf(fp,"\nEnd.\n");
	}
}

hpt_typ *MultiReadHpt(Int4 &N, char *filename)
{
        FILE *fp=open_file(filename,".hpt","r");
        hpt_typ *hpt,*Head=0,*Tail=0;
        Int4	i=0;
        char    c;
        do {
           do { c=fgetc(fp); if(c == EOF) break; } while(c != 'H');
           if(c==EOF) break; else ungetc(c,fp);
           if(Head==0){ Head=new hpt_typ(fp); Tail=Head; }
	   else { hpt=new hpt_typ(fp); Tail->SetNext(hpt); Tail=hpt; }
	   // i++; if(hpt) fprintf(stderr,"%d. read okay.\n",i);
        } while(TRUE); fclose(fp);
        // for(hpt=Head; hpt; hpt=hpt->RtnNext()){ hpt->Put(stderr); } 
	N=i; return Head;
}

void    hpt_typ::PrintError(){ print_error(HPT_PUBLIC_USAGE); }


void    hpt_typ::Init(char *filename)
{
	FileName=0; NumRandom=0; Mode='I';
	NumberBPPS=0; NumElementarySets=0; Next=0;
	FILE *fp=open_file(filename,".hpt","r"); Read(fp); fclose(fp);
	CreateGroups(); ValidityCheck();
}

void    hpt_typ::Init(FILE *fp)
{
	NumRandom=0; FileName=0; Mode='I';
	NumberBPPS=0; NumElementarySets=0; Next=0;
	Read(fp); CreateGroups(); ValidityCheck();
}

void    hpt_typ::ValidityCheck()
{
	Int4	i,j,m,n,s,t;
	BooLean	found;
	for(n=1; n <= NumberBPPS; n++){
	  for(m=n+1; m <= NumberBPPS; m++){
	   found=FALSE;
	   for(s=1; s <= NumElementarySets; s++){
		if(HyperPartition[s][m] != HyperPartition[s][n]) found=TRUE;
	   }
	   if(!found){
		fprintf(stderr,"Tripartitions %d and %d are identical\n",n,m);
		print_error("Fatal hyperpartition syntax error");
	   }
	  }
	}
	for(s=1; s <= NumElementarySets; s++){
	  for(t=s+1; t <= NumElementarySets; t++){
	     found=FALSE;
	     for(n=1; n <= NumberBPPS; n++){
		if(HyperPartition[s][n] != HyperPartition[t][n]) found=TRUE;
	     }
	     if(!found){
		fprintf(stderr,"Sequence sets %d and %d have identical assignments\n",s,t);
		print_error("Fatal hyperpartition syntax error");
	     }
	  }
	}
}

static BooLean IgnorableLine(char *str)
{
	Int4 j;
	if(str[0] == '#') return TRUE;  	// line commented out. 
	if(str[0] == '\n') return TRUE; 	// blank line. 
	for(j=0; isspace(str[j]); ){ if(str[j] == '\n') break; else j++; }
	if(str[j] == '\n') return TRUE; 	
	return FALSE;
}

static void	print_hpt_error(Int4 line,char *str, char *msg)
{
	fprintf(stderr,"Hyperpartition syntax error on line %d:\n\t(\"%s\")\n",line,str);
	fprintf(stderr,"\n\t(\"%s\")\n",msg);
	assert(!"Fatal error");
	print_error(HPT_PUBLIC_USAGE);
}

void	hpt_typ::initialize()
{
	NEWP(HyperPartition, MAX_NUM_ELMENTARY_SETS +2,char);
	for(Int4 i=0; i < MAX_NUM_ELMENTARY_SETS; i++){
		GroupName[i]=0;
		nGroupsFG[i]=0; nGroupsBG[i]=0;
		GroupsFG[i]=0; GroupsBG[i]=0;
		NameElementarySet[i]=0;
		nArgmnt[i]=0; Argmntv[i]=0;
		sst_string[i]=0; ArgV[i]=0; ArgC[i]=0;
		HyperPartition[i]=0;
		OutputColumns[i]=0;
		ArgmntStr[i]=0; PttrnStr[i]=0; QryGroup[i]=0;
	}
}

void	hpt_typ::Read(FILE *fp,BooLean ignore)
// modified this to allow deleted columns and rows to be ignored.
{
	Int4	b,c,s,i,j,g,n,N,M,y,line;
	int	x;
	char	str2[300];

	initialize();

	NumberBPPS=0; OutputCols=0;
	char      str[2003],whatIs,tmp_str[2003],name[2003],*rtn,tmp_arg[2003];
	BooLean	Start=FALSE;
	tmp_str[2001]=0; str[1995]=0;
	for(line=1; fgets(str,2000,fp) != NULL; line++){
	   if(str[1995]!=0) print_hpt_error(line,str,"input file line length over limit");
	   if(IgnorableLine(str)){ line--; continue; } 
	   if(strcmp(str,"Hyperpartition:\n") == 0){ Mode='I'; Start=TRUE; break; }
	   if(sscanf(str,"Hyperpartition(%[^)]):\n",str2) == 1){
		if(FileName) free(FileName); FileName=AllocString(str2); 
		Mode='I'; Start=TRUE; break; 
	   }
	   if(strcmp(str,"HyperParTition:\n") == 0){ Mode='F'; Start=TRUE; break; }
	   if(sscanf(str,"HyperParTition(%[^)]):\n",str2) == 1){
		if(FileName) free(FileName); FileName=AllocString(str2); 
		Mode='F'; Start=TRUE; break; 
	   }
	}
	if(!Start) print_hpt_error(line,str,"no Start site found");
	for(N=M=0; (rtn=fgets(str,2000,fp)) != NULL; line++){
	   if(str[1995]!=0) print_hpt_error(line,str,"input file line length over limit");
	   if(IgnorableLine(str)){ line--; continue; } 
	   if(N==0 && str[0] == '!' || str[0] == '_' || str[0] == '*'){
	     if(sscanf(str,"%[!_*]\n",tmp_str) == 1){
		NumberBPPS = strlen(tmp_str);
// fprintf(stderr,"\"%s\" == len %d\n",tmp_str,NumberBPPS);
		for(c=0,g=1; g <= NumberBPPS; c++,g++){	
		     if(tmp_str[c]=='!') OutputColumns[g]=TRUE;
		     else if(tmp_str[c]=='_') OutputColumns[g]=FALSE;
		     else if(ignore && tmp_str[c]=='*') ;	// ignore...
		     else if(tmp_str[c]=='*') OutputColumns[g]=FALSE;	// ignore...
		     else print_hpt_error(line,str,"column designator syntax error 1");
		}
		OutputCols=AllocString(tmp_str);
		continue;
	     } else print_hpt_error(line,str,"column designator syntax error 2");
	   } else  if(str[0] == '+' || str[0] == '-' ||str[0] == 'o'){
		BooLean AreArgs=FALSE;
		char CC=0,X[10];
	   	if(sscanf(str,"%[+-o] Random=%d.",tmp_str,&x) == 2){
		   NumRandom=x; whatIs='='; sprintf(name,"Random");
		   if(NumRandom < 1) print_hpt_error(line,str,"hpt_typ::Read( ) Random set syntax error 1");
	   	} else if(sscanf(str,"%[+-o] %d.Random=%d.",tmp_str,&y,&x) == 3){
		   NumRandom=x; whatIs='='; sprintf(name,"Random");
		   if(NumRandom < 1) print_hpt_error(line,str,"hpt_typ::Read( ) Random set syntax error 2");
		   // if((N+1) != y) print_hpt_error(line,str,"*.hpt file syntax (set numbering) error 1");
		   if((M+1) != y) print_hpt_error(line,str,"*.hpt file syntax (set numbering) error 1");
		//===========  below: gcc-4.7.2 and gcc-4.9.0 bug fixes. =============
	   	} else if(sscanf(str,"%[+-o] %d%c%[^.?!=*]%c -%[^\n]\n",
				tmp_str,&x,&CC,name,&whatIs,tmp_arg) == 6 
					&& CC == '.' && strchr("[?.!=*]",whatIs) != 0){
			AreArgs=TRUE;
			// fprintf(stderr,"%s\n",tmp_arg);
			// if((N+1) != x) print_hpt_error(line,str,"*.hpt file syntax (set numbering) error");
			if((M+1) != x) print_hpt_error(line,str,"*.hpt file syntax (set numbering) error");
		} else if(sscanf(str,"%[+-o] %d%c%[^.?!=*]%c",
				tmp_str,&x,&CC,name,&whatIs) == 5 && CC=='.' && strchr("[?.!=*]",whatIs) != 0){
			// fprintf(stderr,"tmp_str = \"%s\"; N=%d; M=%d; x=%d; name = \"%s\"; whatIs='%c'.\n",
			//		tmp_str,N,M,x,name,whatIs);
			if((M+1) != x){
			    // WARNING: sscanf on the gcc-4.7.2 compiler is not reading in the x integer!!.
			    // converting ' %d.' to ' %d[.]' in sscanf argument fixed the problem. 
			    fprintf(stderr,
			       "tmp_str = \"%s\"; N=%d; M=%d; x=%d; name = \"%s\"; whatIs='%c'.\n",
					tmp_str,N,M,x,name,whatIs);
			    print_hpt_error(line,str,"*.hpt file syntax (set numbering) error");
			}
	   	} else if(sscanf(str,"%[+-o] %[^.?!=*]%c -%[^\n]\n",
				tmp_str,name,&whatIs,tmp_arg) == 4 && strchr("[?.!=*]",whatIs) != 0){
			AreArgs=TRUE;
			// fprintf(stderr," ********** %s ***********\n",tmp_arg);
		   // then okay 
	   	} else if(sscanf(str,"%[+-o] %[^.?!=*]%c",
				tmp_str,name,&whatIs) == 3 && strchr("[?.!=*]",whatIs) != 0){
			// fprintf(stderr," ********** last option ***********\n");
		   // then okay 
		} else print_hpt_error(line,str,"line fails to match syntax");

		assert(tmp_str[2001]==0);	// make sure no input overflow...
		M++; // if(whatIs=='*') continue;   // skip these rows.

		N++;
	   	if(N >= MAX_NUM_ELMENTARY_SETS) print_hpt_error(line,str,"too many elmentary sets defined");
		if(N==1 && NumberBPPS==0) NumberBPPS = strlen(tmp_str); 
		else if(NumberBPPS != strlen(tmp_str)){
			fprintf(stderr,"\nN = %d; NumberBPPS = %d; string length = %d\n",
								N,NumberBPPS,strlen(tmp_str));
			fprintf(stderr,"\nOuput columns = %s\n", OutputCols);
			print_hpt_error(line,tmp_str,"NumberBPPS != string length");
		} NameElementarySet[N]=AllocString(name);
		sprintf(name," %s",tmp_str);   // want HP to start at one, not zero.
		HyperPartition[N]=AllocString(name);
#if 1	//**************** New set arguments... ********************
		if(AreArgs){	// arguments passed for this set...
			char *tmp_argv[1003];
			sprintf(tmp_str,"-%s",tmp_arg);
  			ArgC[N]=string2argv(tmp_argv,tmp_str);
			if(ArgC[N] < 1){
				fprintf(stderr,"  %s  ",tmp_str);
				print_hpt_error(line,str,"hpt file set argument input error");
  			} else {
  	   			NEWP(ArgV[N],ArgC[N]+3, char);
  	   			for(j=0; j < ArgC[N]; j++){
					ArgV[N][j] = AllocString(tmp_argv[j]);
					// fprintf(stderr,"%s ",ArgV[N][j]);
  	   			} // fprintf(stderr,"\n");
  	   			ArgV[N][j] = 0;
			} 
		} else ArgC[N]=0;
#endif	//**********************************************************************************
		if(ignore){
		   if(strchr(".?=!",whatIs) == NULL){
			 print_hpt_error(line,str,"Terminal punctuation not .?=!");
		   }
		} else if(strchr(".?=!*",whatIs) == NULL){
			print_hpt_error(line,str,"Terminal punctuation not .?=!*");
		} SetType[N]=whatIs;
	   } else break;
	}
	if(N < 1) print_error("hpt file syntax error (no elementary sets).");
	else NumElementarySets=N; 
     char *OldCols=0;
     if(ignore){	// remove deleted columns from each category.
      if(M < N || strchr(OutputCols,'*')!= 0){	// 
	for(x=y=0; x < NumberBPPS; x++){
		if(OutputCols[x] != '*'){ 
		    tmp_str[y] = OutputCols[x];
		    for(n=1; n <= NumElementarySets; n++){
			HyperPartition[n][y+1]=HyperPartition[n][x+1];
		    } y++;
		}
	} NumberBPPS=y;  OldCols=OutputCols; tmp_str[y]=0;  OutputCols=AllocString(tmp_str);
	y++; 
	for(n=1; n <= NumElementarySets; n++){ HyperPartition[n][y]=0; }
      }
     }
	for(Start=FALSE,N=0; rtn != NULL; line++){
	   if(str[1995]!=0) print_error("input file line length over limit");
	   if(strcmp(str,"Settings:\n") == 0){ Start=TRUE; break; }
	   else {
	   	if(strcmp(str,"End.\n") == 0){ break; }	// allows more than one hpt per file
	        if(IgnorableLine(str)){ line--; } 	// line blank or commented out. 
	   	else {
		   print_hpt_error(line,str,"hpt file syntax error (extraneous text).");
		   // print_error("hpt file syntax error (extraneous text).");
		}
	   }
	   if((rtn=fgets(str,2000,fp)) == NULL) break; 
	   if(strcmp(str,"End.\n") == 0){ break; }	// allows more than one hpt per file
	}
	if(Start){
	   for(N=M=0; fgets(str,2000,fp) != NULL; line++){
	     if(str[1995]!=0) print_error("hpt input file line length over limit");
	     if(IgnorableLine(str)){ line--; continue; } 	// line blank or commented out. 
	     if(isdigit(str[0])){
		if(sscanf(str,"%d.%[^\n]\n",&x,tmp_str) != 2){
		   print_error("hpt_typ::Read( ) initialization file syntax error 3");
		} 
		if(OldCols != 0 && OldCols[N]=='*'){ 
		  N++;
		} else {
		  M++; N++;
		  assert(tmp_str[2001]==0);	// make sure no input overflow...
	   	  if(M > NumberBPPS) print_error("hpt_typ::Read(): inconsistent number of groups");
	   	  if(N != x) print_error("Hyperpartition error: inconsistent numbering of groups");
		  assert(N==x);
	          GetSettings(M, tmp_str);
	        } 
	     } else if(strcmp(str,"End.\n") == 0){ break; }	// allows more than one hpt per file
	     else print_error("*.hpt initialization file syntax error 4");
	  } if(OldCols) free(OldCols);
	} else {	// set up default settings...
		for(n=1; n <= NumberBPPS; n++){
			nArgmnt[n]=0; Argmntv[n]=0;
			sprintf(str,"Group%d",n);
			GroupName[n]= AllocString(str);
			ArgmntStr[n]= AllocString(str);
		}
	}
	for(n=1; n <= NumberBPPS; n++){
		if(OutputCols==0) OutputColumns[n]=TRUE;
		PttrnStr[n]=0; 
		if(GroupName[n]==0){ sprintf(str,"Group%d",n); GroupName[n]=AllocString(str); }
	}
	if(OutputCols==0){
		NEW(OutputCols, NumberBPPS+3, char);
		for(n=0; n < NumberBPPS; n++) OutputCols[n]='!';
		OutputCols[n]=0;
	}
	// this->Put(stdout);
}

void    hpt_typ::GetSettings(Int4 N, char *tmp_str)
{
	Int4	j,tmp_argc;
	char *tmp_argv[1003];

	if(N >= MAX_NUM_ELMENTARY_SETS) print_error("too many groups defined");
	ArgmntStr[N]=AllocString(tmp_str);
  	//fprintf(stderr,"%s\n",ArgmntStr[N]);
  	// nArgmnt[N]=string2argv(tmp_argv,tmp_str);
	tmp_argc=string2argv(tmp_argv,tmp_str); nArgmnt[N]=tmp_argc;
	if(nArgmnt[N] < 1){
		fprintf(stderr,"  %s  ",tmp_str);
		print_error("hpt file Settings input error");
  	}
	GroupName[N]=AllocString(tmp_argv[0]);
	if(nArgmnt[N] > 1){
  	   NEWP(Argmntv[N],nArgmnt[N]+3, char);
  	   for(j=1; j < nArgmnt[N]; j++){
		Argmntv[N][j-1] = AllocString(tmp_argv[j]);
		// fprintf(stderr,"%s ",Argmntv[N][j-1]);
  	   } // fprintf(stderr,"\n");
  	   Argmntv[N][j-1] = 0;
	} nArgmnt[N]--;
	for(j=0; j < tmp_argc; j++) free(tmp_argv[j]);
}

void    hpt_typ::CreateGroups( )
//**************** Create HyperPartition. *******************
{
	Int4	n,s,nFG,nBG,nMG;
	BooLean *IsFG,*IsBG;
	NEW(IsFG,NumElementarySets +3, BooLean); NEW(IsBG,NumElementarySets +3, BooLean);
	for(n=1; n<= NumberBPPS; n++){
	  // NEW(GroupsFG[n], NumberBPPS +3,Int4); NEW(GroupsBG[n], NumberBPPS +3,Int4);
	  NEW(GroupsFG[n], NumElementarySets +3,Int4); NEW(GroupsBG[n], NumElementarySets +3,Int4);
	  sst_string[n]=0;
          for(nFG=nBG=nMG=0,s=1; s<= NumElementarySets; s++){
	    switch (HyperPartition[s][n]){
		case '+': nFG++; GroupsFG[n][nFG]=s; IsFG[s]=TRUE; break;
		case '-': nBG++; GroupsBG[n][nBG]=s; IsBG[s]=TRUE; break;
		case 0: case 'o': nMG++; break;
		default: 
		  fprintf(stderr,"HP[%d][%d] = %d\n",s,n,HyperPartition[s][n]);
		  print_error("hpt_typ::CreateGroups( ) this should not happen."); 
		break;
	    }
          } nGroupsFG[n] = nFG; nGroupsBG[n] = nBG;
	  if(nFG == 0){
		fprintf(stderr,"Tripartition %d category lacks a foreground set\n",n);
		this->Put(stderr);
		print_error("Fatal hyperpartition syntax error.");
	  }
	  if(nBG == 0){
		fprintf(stderr,"Tripartition %d category lacks a background set\n",n);
		this->Put(stderr);
		print_error("Fatal hyperpartition syntax error.");
	  }
	}
        for(s=1; s<= NumElementarySets; s++){
	   if(strcmp("Random",NameElementarySet[s]) == 0){
		if(IsFG[s]){
			fprintf(stderr,"Random set %d contains a foreground ('+') assignment\n",s);
			print_error("FATAL: syntax error in *.hpt file");
		}
	   } else {
		if(!IsFG[s]){
			fprintf(stderr,"Set %d ('%s') lacks a foreground ('+') assignment\n",
								s,NameElementarySet[s]);
			print_error("FATAL: syntax error in *.hpt file");
			// print_error(HPT_PUBLIC_USAGE);
		}
		if(FALSE && !IsBG[s]){
			fprintf(stderr,"Set %d ('%s') lacks a background ('-') assignment\n",
								s,NameElementarySet[s]);
			print_error("FATAL: syntax error in *.hpt file");
		}
	   }
	} free(IsFG); free(IsBG);
}

