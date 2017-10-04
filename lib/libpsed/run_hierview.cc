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

#if 0
#include "stdinc.h"
#include "hsc_typ.h"
#include "lha_typ.h"
#include "omc_typ.h"
#else

#include "hierview.h"

#endif

extern	int ChainVSI(int argc,char *argv[],FILE *ifp, FILE *ofp);

#define USAGE_RTF "\n           === Information regarding output rtf contrast alignments: ===\n\n\
     The rtf output file shows, for each leaf node, an alignment (of a few representative sequences\n\
   taken from that node) multiple times, one time for each node on the lineage from (and including)\n\
   the leaf node back to the root node.  Each node's alignment highlights the conserved residues \n\
   most distinctive of the sequences assigned to the subtree rooted at that node.  \n\n\
     Above these alignments is another alignment highlighting the residues conserved solely within\n\
   the representative (leaf node) set.  If phylum information is provided in the input fasta \n\
   sequence file, this first alignment names the phylum from which each sequence was obtained.\n\n\
     Below each of the other alignments is a summary of the most conserved amino acid residues \n\
   at each position; the number of sequences (assigned to the foreground) is given in parentheses \n\
   on the first line.  The 1st to 3rd lines show up to three residues at each position that occur \n\
   both most frequently and in ≥10% of the sequences.  \n\n\
     Directly below this, the frequencies of the designated residues are given in integer tenths; \n\
   for example, an ‘8’ indicates that 80-90% of the sequences in the foreground alignment match \n\
   the corresponding pattern residues.  To highlight larger integers ‘5’ and ‘6’ are shown in black \n\
   and ‘7’-‘9’ in red. The first of these lines (labeled as 'wt_res_freqs' for 'weighted residue \n\
   frequencies') reports the effective number of aligned sequences. In all of these cases, reported \n\
   frequencies have been down-weighting for redundancy.\n\n\
     The black dots above the alignment indicate the pattern positions that were identified by the\n\
   program.   Pattern-matching (correlated) residues are highlighted in color, with biochemically\n\
   similar residues colored similarly.  For example, acidic residues are shown in red, basic residue\n\
   in cyan and hydrophobic residues in yellow; histidine, glycine and proline are each assigned a \n\
   unique color.\n\n\
     The heights of the red bars above the alignment quantify (using a semi-logarithmic scale) the\n\
   degree to which residue frequencies in the foreground diverge at each position from the \n\
   corresponding positions in the background.\n\n\
     Note that for the second alignment, the foreground corresponds to the root node, that is, to\n\
   the entire tree and thus to the entire superfamily, and the background corresponds to all proteins \n\
   unrelated to the superfamily, which is represented by standard amino acid residue frequencies.\n\n\
        ========================================================================\n\n"

#define USAGE_START "Usage: hierview <file_prefix> <int> [options] \n\
      requires <prefix>.cma, <prefix>.tpl and <prefix>.hpt <prefix>.sets input files (from hieraln)\n\
	<int> = the line within the hpt file specifying the subgroup of interest.\n\
        -pdb=<file> - file with paths to pdb coordinates e.g., /tmp/pdb_nat/4ag9_H.pdb).\n\
        		    list one pdb file per line\n\
        -h   helpful information regarding the *rtf output file.\n\n"

#define PUBLIC_START "\
	<int> = the line within the hpt file specifying the subgroup of interest.\n\
   Input: himsa files created by BPPS 2 (*.cma, *.tpl, *.hpt, *.sets)\n\
   Output:\n\
        <prefix>.tbl: hpt lineage from root to <int>th node.\n\
        <prefix>_aln.rtf: linage constrast alignments for the <int>th node.\n\
      if -pdb option is used and matching pdb sequences are found:\n\
        <prefix>.sprs: input file for SIPRIS analyses.\n\
        <prefix>_<int>_<int1>_<int2>.pml: <int>th node PyMol script with <int1>_<int2> identifier.\n\
   Options:\n\
        -pdb=<file> - file with paths to pdb coordinates e.g., /tmp/pdb_nat/4ag9_H.pdb).\n\
        		    list one pdb file per line\n\
        -h   Information regarding the <prefix>_aln.rtf file.\n\n"


static void PrintError(char SecretCode,char *prgm_name)
{
	if(SecretCode > 0){
                fprintf(stderr,"Usage: %s %d <prefix> <int> [options] \n",prgm_name,(int)SecretCode);
                print_error(PUBLIC_START);
        } else { print_error(USAGE_START); }
}


int	run_hierview(int argc,char *argv[]) { return run_hierview(0,argc,argv); }

int	run_hierview(char SecretCode,int argc,char *argv[])
{ 
	Int4	rtn,N,i,j,arg,Argc,focus=0;
	char	c,*Argv[20],str[502],*pdb_names=0,str0[200],str1[200];

	if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') print_error(USAGE_RTF);
	if(argc == 3 && argv[2][0] == '-' && argv[2][1] == 'h') print_error(USAGE_RTF);
	if(argc < 3) PrintError(SecretCode,argv[0]);
        if(sscanf(argv[2],"%d",&N)!=1) PrintError(SecretCode,argv[0]);
        for(arg=3; arg < argc; arg++){
           if(argv[arg][0] != '-') PrintError(SecretCode,argv[0]);
           switch(argv[arg][1]) {
              case 'p':
		if(sscanf(argv[arg],"-pdb=%s",str)==1){
			pdb_names=AllocString(str);
		} else PrintError(SecretCode,argv[0]); break;
              case 'h': print_error(USAGE_RTF); break;
              default: PrintError(SecretCode,argv[0]); break;
           }
        }
	fprintf(stderr,"======== 1. Creating lineage MSA (liMSA) ========\n");
	// sprintf(str,"hybrid %s %s ",argv[1],argv[2]);
	sprintf(str,"%s %s %s -x ",argv[0],argv[1],argv[2]);
	Argc=string2argv(Argv,str);
	for(i=0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n");
	// 1. run hybrid.cc : name 3 
	// needs <file_name>.hpt, <file_name>.tpl, <file_name>.cma files from hieraln.
	FILE *ofp,*mfp=0,*smfp=0,*stfp=0,*hfp=0,*xfp,*tmpfp,*vsifp;
	mfp=tmpfile(); hfp=tmpfile(); stfp=tmpfile(); smfp=tmpfile(); 
	run_hybrid(Argc,Argv,mfp,hfp,stfp,smfp);
	rewind(mfp); rewind(smfp); rewind(stfp); rewind(hfp);
#if 0
char Str[9002];
while(fgets(Str,9000,mfp) != NULL) fprintf(stderr,"%s",Str);
rewind(mfp);
#endif
	// create an sma file with key sequences in it!!
	// 2. run omcBPPS to create rtf contrast alignment file
	// looks for : name_3 -print=P -rtf
	for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);

	fprintf(stderr,"======== 2. BPPS column optimization; printing contrast alignment(s). ========\n");
#if 0	// -del == don't treat BG as random residues!!!
	{ hpt_typ hpt(argv[1]); focus=hpt.NodeDepth(N); }
	sprintf(str,"omcBPPS %s_%d -print=o -rtf -focus=%d -maxcol=25 -minnats=5",argv[1],N,focus);
#else
	sprintf(str,"omcBPPS %s_%d -run=c -rtf -maxcol=30 -minnats=5",argv[1],N);
#endif
	FILE *mmafp=0,*hptfp=0,*ptrnfp=0;
	mmafp=tmpfile(); hptfp=tmpfile(); ptrnfp=tmpfile();
	Argc=string2argv(Argv,str);
	{  // don't use SecretCode=3; prevents hptfp output ...
	   omc_typ omc(SecretCode,Argc,Argv,mfp,hfp,stfp,smfp); omc.VerboseOff(); 
	   rtn=omc.Run(mmafp,hptfp,ptrnfp); omc.PrintTime(stderr); 
	} // fclose(mfp); fclose(smfp); fclose(stfp); fclose(hfp);
	// WARNING: mfp, smfp, stfp, hfp files are closed by omc_typ!
	rewind(mmafp); rewind(hptfp); rewind(ptrnfp);
	// while((c=fgetc(mmafp)) != EOF){ fprintf(stderr,"%c",c); }  exit(1);
	for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);

	if(pdb_names == 0) return 0;

	// 3. run sarp if passed in pdb_names: (e.g., sarp name_3) 
	// needs <mcBPPS_prefix>_new.hpt <mcBPPS_prefix>.pttrns <mcBPPS_prefix>_new.mma

	fprintf(stderr,"======== 3. Mapping covariant patterns to structures (if available). ========\n");
	Argc=0; Argv[Argc]=AllocString("sarp");  Argc++;
	Argv[Argc]=pdb_names;  Argc++;
	sprintf(str,"%s_%d",argv[1],N);
	Argv[Argc]=AllocString(str);  Argc++;
	for(i = 0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n");
#if 1
	xfp=tmpfile(); vsifp=tmpfile(); fprintf(vsifp,"\n");
	BooLean	verbose=FALSE;
#else
	xfp=vsifp=0;
#endif
	{ hsc_typ hsc(Argc, Argv,mmafp,hptfp,ptrnfp); hsc.VerboseOff(); rtn=hsc.Run(0,3,xfp,vsifp); }
	if(xfp) rewind(xfp); if(vsifp) rewind(vsifp);
	// mmafp and hptfp files are closed by hsc_typ.
	for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);

	// 4. run chn_vsi if passed in pdb names.
	fprintf(stderr,"======== 4. Creating pymol files. ========\n");
	Int4	file,id;
	if(rtn == 0){ fprintf(stderr,"No structures found\n"); }
	else {
	  if(xfp==0){ xfp=open_file(str,"_pdb.log","r"); }
	  for(j=1; fgets(str,500,xfp) != NULL; ){
if(verbose) fprintf(stderr,"%s\n",str);
	   if(sscanf(str,"    ChnSARP %s %d %d",str0,&id,&file) == 3){
	     sprintf(str1,"chn_vsi %s %d %d -T -skip=W -d2.5 -D",str0,id,file);
// fprintf(stderr,"%s\n",str1);
	     Argc=string2argv(Argv,str1);	// mode == 'T'
if(verbose){ for(i=0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n"); }

	     FILE *ifp,*ofp=tmpfile();
	     ChainVSI(Argc,Argv,vsifp,ofp); rewind(ofp); ifp=ofp; if(vsifp) rewind(vsifp);
	     for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);

	     sprintf(str1,"chn_vsi %s.crs %d -d2.5 -c -D -pml=%s_%d_%d_%d",str0,id,argv[1],N,file,id); j++;
// fprintf(stderr,"%s\n",str1);
	     Argc=string2argv(Argv,str1);	// mode == 'p'
if(verbose){ for(i = 0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n"); }

	     ChainVSI(Argc,Argv,ifp,0); fclose(ifp);
	     for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);
	   }
	  } fclose(xfp); if(vsifp) fclose(vsifp);
	} return 0;
}

