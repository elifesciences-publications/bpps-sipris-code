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

#include "ras_typ.h"

ras_typ::ras_typ(char *filename, BooLean cis_trans_pro)
{

	AB=MkAlpha("XCGASTNDEQKRHWYFVILMP",NULL);
	pdbfile=AllocString(filename);
	CIS_TRANS_PRO=cis_trans_pro;
	Init( );
}

void    ras_typ::Init( )
{ 
	static char AARes0[22][4]= {"unk","cys","gly","ala","ser","thr","asn","asp", 			"glu","gln","lys","arg","his","trp","tyr","phe","val","ile","leu","met","pro"};
	for(Int4 i=0;i <=20; i++) strncpy(AARes[i],AARes0[i],4);
}

void    ras_typ::Free( )
{
	NilAlpha(AB);
	free(pdbfile);
}

void	ras_typ::PrintHEADER(FILE *fp,BooLean use_path,char c,Int4 wire_width)
{
        char    *DATADIR=0;
	if(use_path){
          DATADIR=getenv("CHN_RAS_DIR");
          fprintf(stderr,"CHN_RAS_DIR = %s\n",DATADIR);
          if(DATADIR==0) fprintf(fp,"load pdb /usr/molbio/pdb/%s\n",pdbfile);
          else fprintf(fp,"load pdb %s/%s\n",DATADIR,pdbfile);
	} else { fprintf(fp,"load pdb %s\n",pdbfile); }
        // CHANGE TO ENVIRONMENTAL VARIBLE: /usr/molbio/pdb

	fprintf(fp,"set background white\n");
	fprintf(fp,"wireframe off\n");
	// fprintf(fp,"color white\nset strands 1\n");
	// if(USE_TRACE) fprintf(fp,"trace\n");
	// else fprintf(fp,"cartoon\n");
	if(c){ 
	  fprintf(fp,"select not *%c\nstrands\ncolor gray\n",c);
	  // if(USE_TRACE) fprintf(fp,"select not *%c\nstrands\ncolor gray\n",c);
	  // else fprintf(fp,"select not *%c\nbackbone 50\ncolor gray\n",c);
	}
        fprintf(fp,"select ligand and not (hydrogen,tpo,ptr)\n");
        // fprintf(fp,"color [200,150,0]\n");
        // fprintf(fp,"color violet\n");
        fprintf(fp,"color cyan\n");
        // fprintf(fp,"color [150,255,150]\n");	// pale green
	if(c){
          fprintf(fp,"wireframe %d\n",wire_width);
          // fprintf(fp,"select not ligand and not (carbon,oxygen,nitrogen,phosphorus,hydrogen,sulfur)\n");
          // fprintf(fp,"select not (carbon,oxygen,nitrogen,phosphorus,hydrogen,sulfur)\n");
          fprintf(fp,"select (iron,zinc,mg,mn,cobalt,nickel,SODIUM,POTASSIUM,CHLORINE,CHROMIUM,COPPER,MOLYBDENUM,SILVER,GOLD)\n");
          // fprintf(fp,"spacefill\ncolor cpk\n");
          fprintf(fp,"spacefill 250\ncolor [0,144,0]\n");	// dark green.
	}
}

void	ras_typ::PrintTAIL(FILE *fp)
{
        fprintf(fp,"select all\n");
        fprintf(fp,"set specular on\nset specpower 100\n");
	fputc('\n',fp);
	fflush(fp);
}

void	ras_typ::PrintTrueColor(FILE *fp, char c)
{
	short	R,G,B;
	// fprintf(stderr,"c = %c\n",c);
   if(c == 'X' || c == 'x') fprintf(fp,"color cpk\n");
   else {
	// c=toupper(c);
	switch(c){
	  case 'M': G=0; R=B=255; break; // Old majenta...
	  case 'm': G=0; R=B=200; break; // Old majenta...
	  case 'R': R=255; G=B=0; break; 
	  case 'r': R=200; G=B=0; break;
	  case 'Y': R=G=255; B=0; break;
	  case 'y': R=G=200; B=0; break;
	  case 'G': R=0; G=255; B=0; break;
	  case 'g': R=0; G=200; B=0; break;
	  case 'C': R=0; G=B=255; break;
	  // case 'C': R=100; G=200; B=200; break;
	  case 'c': R=0; G=B=200; break;
	  case 'O': R=220; G=140; B=60; break;
	  case 'o': R=200; G=110; B=0; break;
	  case 'B': R=G=0; B=255; break;
	  default: R=G=B=0; break;
	} 
	fprintf(fp,"color [%d,%d,%d]\n",R,G,B);
   }
}

void	ras_typ::PrintColor(FILE *fp, char c)
{
	short	R,G,B;
	// fprintf(stderr,"c = %c\n",c);
   if(c == 'X' || c == 'x') fprintf(fp,"color cpk\n");
   else {
	switch(c){
	  case 'T': R=206; G=167; B=128; break;	// 'amber == lt. brown
	  case 't': R=206; G=167; B=128; break;	// 'amber == lt. brown
	  case 'A': R=190; G=B=115; break;	// tan...
	  case 'a': R=190; G=B=115; break;
	  case 'M': G=100; R=B=180; break; // Old majenta...
	  // case 'm': G=200; R=B=255; break;
	  case 'm': G=180; R=B=255; break;
	  // case 'R': R=255; G=B=85; break;
	  // case 'R': R=220; G=B=120; break; // for orginal Ran CHAIN analysis
	  // case 'R': R=210; G=B=110; break; // for next two papers.
	  // case 'R': R=200; G=B=105; break; 
	  case 'R': R=190; G=B=100; break; 
	  // case 'r': R=255; G=B=170; break;
	  case 'r': R=255; G=B=135; break;
	  case 'S': R=255; G=156; B=0; break; // gold == special residues.
	  case 's': R=255; G=190; B=65; break; // gold == special residues.
	  case 'O': R=200; G=120; B=40; break;
	  // case 'o': R=255; G=215; B=130; break;
	  case 'o': R=255; G=180; B=100; break; // new
	  case 'Y': R=G=160; B=50; break;
	  case 'y': R=255; G=255; B=160; break;
	  case 'G': R=110; G=190; B=75; break;
	  case 'g': R=B=200; G=255; break;
	  case 'C': R=100; G=200; B=200; break;
	  case 'c': R=180; G=B=255; break;
	  case 'B': R=G=110; B=240; break;
	  case 'b': R=G=190; B=255; break;
	  case 'P': R=130; G=50; B=180; break;
	  case 'p': R=200; G=150; B=255; break;
	  case 'W': G=R=B=255; break;
	  case 'w': G=R=B=230; break;
	  case 'L': G=R=B=150; break;
	  case 'l': G=R=B=200; break;
	  case 'D': G=R=B=100; break;
	  case 'd': G=R=B=125; break;
	  default: R=G=B=0; break;
	} 
	fprintf(fp,"color [%d,%d,%d]\n",R,G,B);
   }
}

void	ras_typ::FixAtom(char *atom)
// bug in rasmol does not allow strings like: "ANP1413B.*ho3"
// This routine fixes this...
{ 
  Int4	i;
  for(i=0; atom[i]; i++) if(atom[i]=='*') atom[i]='?'; 
  i--; if(atom[i] == '?') atom[i]='*'; // asteric on end seems to be okay.
}

void	ras_typ::ColorAtom(FILE *fp, char *mol, char c, char color, Int4 i, char *atom)
{
	char	atm;

	if(atom){
	  for(Int4 j=0; !isalpha(atm=atom[j]); j++){
		if(atm==0) print_error("ColorAtom( ) syntax error");
	  }
	} else atm='x';
	if(c){
	 	if(atom) fprintf(fp,"select %s%d%c.%s\n",mol,i,c,atom);
	 	else fprintf(fp,"select %s%d%c\n",mol,i,c);
	} else {
		if(atom) fprintf(fp,"select %s%d.%s\n",mol,i,atom);
		else fprintf(fp,"select %s%d\n",mol,i);
	}
	if(color == 'X'){
	   switch(atm){
	    // case 'n': fprintf(fp,"color [143,143,255]\n"); break;
	    case 'n': fprintf(fp,"color blue\n"); break;
	    case 'o': fprintf(fp,"color red\n"); break;
	    case 'c': fprintf(fp,"color [200,200,200]\n"); break;
	    case 's': fprintf(fp,"color [255,200,50]\n"); break;
	    case 'h': fprintf(fp,"color white\n"); break;
	    // case 'c': fprintf(fp,"color lightgrey\n"); break;
	    default: fprintf(fp,"color cpk\n");
	   }
	}
	else if(color != 'x') PrintColor(fp,color);
	// else if(color != 'x') PrintColor(fp,tolower(color));
}

void	ras_typ::print_bond_only(FILE *fp,char res[4],Int4 i,char c,Int4 width,char color,
		const char *atom1,const char *atom2)
{
#if 0	// Fix problem with rasmol not recognizing arg.hh21 instead of h21, etc.
	if(strcmp("arg",atom1) == 0 && strncmp("hh",atom2,2) == 0){
	  if(c){ fprintf(fp,"select %s%d%c.%s,%s%d%c.%s\n",res,i,c,atom1,res,i,c,atom2+1); }
	  else { fprintf(fp,"select %s%d.%s,%s%d.%s\n",res,i,atom1,res,i,atom2+1); }
	}
#endif
	if(c){
	  fprintf(fp,"select %s%d%c.%s,%s%d%c.%s\n",res,i,c,atom1,res,i,c,atom2);
	} else {
	  fprintf(fp,"select %s%d.%s,%s%d.%s\n",res,i,atom1,res,i,atom2);
	}
	if(color == 0 || isupper(color)) fprintf(fp,"wireframe %d\n",width);
	if(color) PrintColor(fp,tolower(color));
}

void	ras_typ::PrintMolecule(FILE *fp, char *mol, char c, char color, Int4 i, 
	Int4 Width,Int4 spacefill)
{
	if(c != 0) fprintf(fp,"select %s%d%c and not hydrogen\n",mol,i,c);
	else fprintf(fp,"select %s%d and not hydrogen\n",mol,i);
        if(color=='X'){ 
		if(spacefill == 9999) fprintf(fp, "spacefill\n");
		else fprintf(fp, "spacefill %d\n",spacefill);
		fprintf(fp, "wireframe %d\n",Width);
	} else {
		if(spacefill == 9999) fprintf(fp, "spacefill\n");
		else fprintf(fp, "spacefill %d\n",spacefill);
		fprintf(fp, "wireframe %d\n",Width);
	   	// else fprintf(fp,"wireframe %d\n",Width);
	}
        PrintTrueColor(fp,color);
}

void	ras_typ::PrintResidue(FILE *fp, char aa, char c, char color, Int4 i, 
	BooLean sideonly,Int4 Width, Int4 big_spacefill)
{
	Int4	r;
#if 0
	Int4	width=THIN_WIRE_WIDTH;
    	if(isupper(color)) width=WIRE_WIDTH; 
	if(Width > 0) width=Width;
#else
	Int4    width=Width;
#endif
	if(aa=='G'){
	  if(c){ fprintf(fp,"select gly%d%c.ca\n",i,c); }
	  else { fprintf(fp,"select gly%d.ca\n",i); }
	  if(isupper(color)) {
		if(big_spacefill == 9999) fprintf(fp, "spacefill\n");
		else fprintf(fp,"spacefill %d\n",big_spacefill);
	  }
	} else {
	  r = AlphaCode(aa,AB);
	  if(c){ 
	    if(aa=='P'){
	      if(CIS_TRANS_PRO && i > 0){
	       	fprintf(fp,
		"select %s%d%c.ca or %s%d%c.n or (%s%d%c and sidechain and not hydrogen)",
		   AARes[r],i,c,AARes[r],i,c,AARes[r],i,c);
	       	fprintf(fp," or (%d%c and (*.c or *.o))\n", i-1,c);
	      } else {
	       	fprintf(fp,
		"select %s%d%c.ca or %s%d%c.n or (%s%d%c and sidechain and not hydrogen)\n",
		   AARes[r],i,c,AARes[r],i,c,AARes[r],i,c);
	      }
	    } else 
		fprintf(fp,"select %s%d%c.ca or (%s%d%c and sidechain and not hydrogen)\n",
		AARes[r],i,c,AARes[r],i,c);
	    if(isupper(color)) fprintf(fp,"wireframe %d\n",width);
	    else fprintf(fp,"wireframe %d\n",width);
	    if(sideonly && aa != 'P') 
		fprintf(fp,"select %s%d%c and sidechain and not hydrogen\n",AARes[r],i,c);
	  } else {
	    if(aa=='P'){
	      if(CIS_TRANS_PRO && i > 0){
	    	fprintf(fp,
			"select %s%d.ca or %s%d.n or (%s%d and sidechain and not hydrogen)",
			AARes[r],i,AARes[r],i,AARes[r],i);
	       	fprintf(fp," or (%d and (*.c or *.o))\n", i-1);
	      } else {
	    	fprintf(fp,
			"select %s%d.ca or %s%d.n or (%s%d and sidechain and not hydrogen)\n",
			AARes[r],i,AARes[r],i,AARes[r],i);
	      }
	    } else fprintf(fp,"select %s%d.ca or (%s%d and sidechain and not hydrogen)\n",
		AARes[r],i,AARes[r],i);
	    if(isupper(color)) fprintf(fp,"wireframe %d\n",width);
	    else fprintf(fp,"wireframe %d\n",width);
	    if(sideonly && aa != 'P') 
		fprintf(fp,"select %s%d and sidechain and not hydrogen\n",AARes[r],i);
	  }
	} 
	if(aa=='P'){
	   if(c){
	       	fprintf(fp,"select %s%d%c.ca or %s%d%c.n\n",AARes[r],i,c,AARes[r],i,c);
		fprintf(fp,"wireframe off\n");
	       	fprintf(fp,"select (%s%d%c and sidechain and not hydrogen)\n",AARes[r],i,c);
	   } else {
	       	fprintf(fp,"select %s%d.ca or %s%d.n\n",AARes[r],i,AARes[r],i);
		fprintf(fp,"wireframe off\n");
	       	fprintf(fp,"select (%s%d and sidechain and not hydrogen)\n",AARes[r],i);
	   }
	}
	if(isupper(color)) PrintColor(fp,tolower(color));
	else PrintColor(fp,toupper(color));
}

void    ras_typ::PutResidueItem(FILE *fp,char a,Int4 i,char *atom1,char *atom2,char chn,char color,
		Int4 wire_width,Int4 spacefill)
{
    if(atom1){
        if(atom2) PrintResidueAtoms(fp,a,i,atom1,atom2,chn,color,wire_width,spacefill);
        else PrintResidueAtom(fp,a,chn,color,i,atom1,TRUE,spacefill);
    } else PrintResidue(fp,a,chn,color,i,TRUE,wire_width,spacefill);
}

void    ras_typ::PrintResidueAtoms(FILE *fp, char a, Int4 i, char *atom1, char *atom2,
	char chn, char color,Int4 wire_width,Int4 spacefill)
{
        char *mol= AARes[AlphaCode(a,AB)];
	PrintMoleculeAtoms(fp, mol, i, atom1, atom2,chn, color,wire_width,spacefill);
}

void    ras_typ::PutMoleculeItem(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
                        char chn, char color,Int4 wire_width,Int4 spacefill)
{
    if(atom1){
	if(atom2) PrintMoleculeAtoms(fp,mol, i,atom1,atom2,chn,color,wire_width,spacefill);
        else PrintMoleculeAtom(fp,mol,chn,color,i,atom1,TRUE,spacefill);
    } else PrintMolecule(fp,mol,chn,color,i,wire_width,spacefill);
}

void	ras_typ::PrintTrace(FILE *fp, Int4 start, Int4 end, char c,char color,char mode,Int4 diameter)
{
	// fprintf(stderr,"i=%d; j=%d; c=%c\n",i,j,color);
	if(mode == '@'){
	  if(c) fprintf(fp,"select (%d-%d) and *%c\n",start,end,c);
	  else fprintf(fp,"select (%d-%d)\n",start,end);
	  PrintColor(fp,color);
	  return;
	}
	if(mode == '_'){
	  if(c){
		fprintf(fp,"select (%d-%d) and *%c and backbone and (carbon,nitrogen)\n",
			start,end,c);
		// and not (%d%c.n,%d%c.c) // except proline!
	  } else {
		fprintf(fp,"select (%d-%d) and backbone and (carbon,nitrogen)\n",start,end); 
	  }
	} else if(c){
	   if(mode=='^') fprintf(fp,"select %d-%d%c\n",start,end,c);
	   else if(mode==':') fprintf(fp,"select %d-%d%c\n",start,end,c);
	   else if(mode=='+') fprintf(fp,"select %d-%d%c and not water\n",start,end,c);
	   else fprintf(fp,"select %d-%d and *%c and backbone\n",start,end,c);
	} else { 
	   if(mode==':') fprintf(fp,"select %d-%d\n",start,end);
	   if(mode=='+') fprintf(fp,"select %d-%d\n",start,end);
	   else fprintf(fp,"select %d-%d and backbone\n",start,end); 
	}
	// if(mode == '.'){ fprintf(fp,"cartoon off\n"); fprintf(fp,"trace\n"); }
	if(mode == '.'){
		if(diameter > 0) fprintf(fp,"trace %d\n",diameter); 
		else fprintf(fp,"trace\n"); 
	} else if(mode == '~'){ 
		if(diameter > 0) fprintf(fp,"cartoon %d\n",diameter); 
		else fprintf(fp,"cartoon on\n"); 
	} else if(mode == ':'){
		fprintf(fp,"spacefill\n"); 
	} else if(mode == '^'){
		if(diameter > 0) fprintf(fp,"ribbons %d\n",diameter); 
		else fprintf(fp,"ribbons\n"); 
	} else if(mode == '+'){
		if(diameter > 0) fprintf(fp,"wireframe %d\n",diameter); 
		else fprintf(fp,"wireframe \n"); 
	} else { fprintf(fp,"backbone on\n"); }
	PrintColor(fp,color);
}

void	ras_typ::PrintClear(FILE *fp)
{ fprintf(fp,"select all\ntrace off\nstrands off\nwireframe off\nspacefill off\n"); }

void	ras_typ::PrintView(FILE *fp, Int4 rx,Int4 ry,Int4 rz,Int4 tx,Int4 ty,Int4 tz)
{
        fprintf(fp,"reset\n");
        fprintf(fp,"rotate x %d\n",rx);
        fprintf(fp,"rotate y %d\n",ry);
        fprintf(fp,"rotate z %d\n",rz);
        // fprintf(fp,"translate x %d\ntranslate y %d\n",tx,ty);
        fprintf(fp,"translate x %.2f\ntranslate y %.2f\n",(double)tx/100.0,(double)ty/100.0);
        fprintf(fp,"zoom %d\n",tz);
}

char	*ras_typ::DefineMolCloud(FILE *fp, char *mol, Int4 site, char *atom1, char *atom2, char chain)
{
    char str[50];

    if(atom1){
      if(atom2){
	if(chain==0) sprintf(str,"%s%d.%s",mol,site,atom2);
	else sprintf(str,"%s%d%c.%s",mol,site,chain,atom2);
      } else {
	if(chain==0) sprintf(str,"%s%d.%s",mol,site,atom1);
	else sprintf(str,"%s%d%c.%s",mol,site,chain,atom1);
      }
    } else {
	if(chain==0) sprintf(str,"%s%d",mol,site);
	else sprintf(str,"%s%d%c",mol,site,chain);
    } return AllocString(str);
}

char	*ras_typ::DefineResCloud(FILE *fp, char aa, Int4 site, char *atom1, char *atom2, char chain)
{
	char str[50];
	Int4  r= AlphaCode(aa,AB);

    if(atom1){
      if(atom2){
	if(chain==0) sprintf(str,"%s%d.%s",AARes[r],site,atom2);
	else sprintf(str,"%s%d%c.%s",AARes[r],site,chain,atom2);
      } else {
	if(chain==0) sprintf(str,"%s%d.%s",AARes[r],site,atom1);
	else sprintf(str,"%s%d%c.%s",AARes[r],site,chain,atom1);
      }
    } else {
	return DefineResCloud(fp,aa,site,chain);
    } return AllocString(str);
}

char	*ras_typ::DefineResCloud(FILE *fp, char aa, Int4 site, char chain)
{
	char	str[50],str2[50];
	Int4	r= AlphaCode(aa,AB);

	if(chain==0) { sprintf(str2,"%c%dX",aa,site); sprintf(str,"%d",site); }
	else { sprintf(str2,"%c%d%c",aa,site,chain); sprintf(str,"%d%c",site,chain); }
	fprintf(fp,"define %s ",str2);
	switch (aa){
	  case 'A':
	 	fprintf(fp,"%s%s.CB\n",AARes[r],str);
	  break;
	  case 'C':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.SG\n",AARes[r],str);
	  break;
	  case 'I':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG2\n",AARes[r],str);
	  break;
	  case 'L':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD2\n",AARes[r],str);
	  break;
	  case 'V':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG2\n",AARes[r],str);
	  break;
	  case 'Y': case 'F': 
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE2,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD2,",AARes[r],str);
	 	fprintf(fp,"%s%s.CZ\n",AARes[r],str);
	  break;
	  case 'R':
	 	fprintf(fp,"%s%s.NE,",AARes[r],str);
	 	fprintf(fp,"%s%s.NH1,",AARes[r],str);
	 	fprintf(fp,"%s%s.NH2,",AARes[r],str);
	 	fprintf(fp,"%s%s.CZ\n",AARes[r],str);
	  break;
	  case 'P':
	 	// fprintf(fp,"%s%s.CA,",AARes[r],str);
	 	// fprintf(fp,"%s%s.N,",AARes[r],str);
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD\n",AARes[r],str);
	  break;
	  case 'D':
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.OD2,",AARes[r],str);
	 	fprintf(fp,"%s%s.OD1\n",AARes[r],str);
	  break;
	  case 'E':
	 	fprintf(fp,"%s%s.CD,",AARes[r],str);
	 	fprintf(fp,"%s%s.OE2,",AARes[r],str);
	 	fprintf(fp,"%s%s.OE1\n",AARes[r],str);
	  break;
	  case 'K':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE,",AARes[r],str);
	 	fprintf(fp,"%s%s.NZ\n",AARes[r],str);
	  break;
	  case 'M':
	 	fprintf(fp,"%s%s.CB,",AARes[r],str);
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.SD,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE\n",AARes[r],str);
	  break;
	  case 'N':
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.ND2,",AARes[r],str);
	 	fprintf(fp,"%s%s.OD1\n",AARes[r],str);
	  break;
	  case 'H':
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE1,",AARes[r],str);
	 	fprintf(fp,"%s%s.NE2,",AARes[r],str);
	 	fprintf(fp,"%s%s.ND1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD2\n",AARes[r],str);
	  break;
	  case 'W':
	 	fprintf(fp,"%s%s.CG,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CD2,",AARes[r],str);
	 	fprintf(fp,"%s%s.NE1,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE2,",AARes[r],str);
	 	fprintf(fp,"%s%s.CE3,",AARes[r],str);
	 	fprintf(fp,"%s%s.CZ2,",AARes[r],str);
	 	fprintf(fp,"%s%s.CZ3,",AARes[r],str);
	 	fprintf(fp,"%s%s.CH2\n",AARes[r],str);
	  break;
	  default: 
	 	fprintf(fp,"%s%s\n",AARes[r],str);
	  break;
	}
	return AllocString(str2);
}

void	ras_typ::PrintResidueAtom(FILE *fp, char aa, char c, char color, Int4 i, 
		char *atom, BooLean sideonly,Int4 big_spacefill)
{
	char *mol = AARes[AlphaCode(aa,AB)];
	PrintMoleculeAtom(fp,mol,c,color,i,atom,sideonly,big_spacefill);
}

void	ras_typ::PrintMoleculeAtom(FILE *fp,char *mol, char chn, char color, Int4 i, 
		char *atom, BooLean sideonly,Int4 spacefill)
{
        FixAtom(atom);
        char chnX;
        if(chn==' ') chnX = 0;
        else chnX = chn;
        ColorAtom(fp,mol,chnX,color,i,atom);
        if(spacefill == 9999) fprintf(fp,"spacefill\n");
        else fprintf(fp,"spacefill %d\n",spacefill);
}

void    ras_typ::PrintMoleculeAtoms(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
	char chn, char color,Int4 wire_width,Int4 spacefill)
{
        Int4    r,k;
        char chnX,atm;

        FixAtom(atom1); FixAtom(atom2);
        if(chn==' ') chnX = 0;   
        else chnX = chn;              
        if(chnX) fprintf(fp,"select %s%d%c.%s,%s%d%c.%s\n",mol,i,chnX,atom1,mol,i,chnX,atom2);
        else fprintf(fp,"select %s%d.%s,%s%d.%s\n", mol,i,atom1,mol,i,atom2);      
        fprintf(fp,"wireframe %d\n",wire_width);
        {       

          for(k=0; !(isalpha(atm=atom1[k]) || atm=='?'); k++){
                if(atm==0) print_error("atom1 syntax error");
          }     
          if(atm != 'c'){               
                ColorAtom(fp,mol,chnX,color,i,atom1);
                if(isupper(color)){
		    if(spacefill == 9999) fprintf(fp,"spacefill\n");
		    else fprintf(fp,"spacefill %d\n",spacefill);
		}
          }     
          for(k=0; !(isalpha(atm=atom2[k]) || atm=='?'); k++){
                if(atm==0) print_error("atom1 syntax error");
          }     
          if(atm != 'c'){               
                ColorAtom(fp,mol,chnX,color,i,atom2);
                if(isupper(color)){
		    if(spacefill == 9999) fprintf(fp,"spacefill\n");
		    else fprintf(fp,"spacefill %d\n",spacefill);
		}
          }     
        }                             
}       

void	ras_typ::PrintCommand(FILE *fp,char *cmd,char type)
{
        switch (type){
          case 'C': PrintClear(fp); break;
          case 'E': fprintf(fp,"%s\n",cmd); break;
        }
}


