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

#include "pml_typ.h"

pml_typ::pml_typ(char *filename, BooLean cis_trans_pro)
{

	AB=MkAlpha("XCGASTNDEQKRHWYFVILMP",NULL);
	pdbfile=AllocString(filename);
	CIS_TRANS_PRO=cis_trans_pro;
	NEWP(Object,130,set_typ);	// this is to keep track of defined sets...
	for(char c='A'; c <= 'z' ; c++) NEW(Object[c],130,set_typ);
	SetSize=100000;
	NEWP(head,300,obj_typ); NEWP(last,300,obj_typ);
	static char AARes0[22][4]= {"unk","cys","gly","ala","ser","thr","asn","asp",
		"glu","gln","lys","arg","his","trp","tyr","phe","val","ile","leu","met","pro"};
	for(Int4 i=0;i <=20; i++) strncpy(AARes[i],AARes0[i],4);
}

void    pml_typ::Free( )
{
	NilAlpha(AB);
	for(char i='A'; i <= 'z'; i++){
	  if(Object[i] != 0){
	    for(char j='A'; j <= 'Z'; j++){ if(Object[i][j]) NilSet(Object[i][j]); }
	  } free(Object[i]); 
	} free(Object);
	for(char c='A'; c <= 'Z'; c++){ if(head[c] != 0) delete head[c]; } 
	free(head); free(last); free(pdbfile);
}

void	pml_typ::PrintHEADER(FILE *fp,BooLean use_path,char c,Int4 wire_width)
{
        char    *DATADIR=0;
	if(use_path){
          DATADIR=getenv("CHN_PML_DIR");
          fprintf(stderr,"CHN_PML_DIR = %s\n",DATADIR);
          if(DATADIR==0) fprintf(fp,"load pdb /usr/molbio/pdb/%s\n",pdbfile);
          else fprintf(fp,"cmd.load(\"%s/%s\")\n",DATADIR,pdbfile);
	} else { fprintf(fp,"cmd.load(\"%s\")\n",pdbfile); }
        // CHANGE TO ENVIRONMENTAL VARIBLE: /usr/molbio/pdb

	// fprintf(fp,"cmd.bg_color(\"white\")\n"); // keep as black...
	fprintf(fp,"cmd.set(\"seq_view\",1)\n");
#if 0
	fprintf(fp,"cmd.select(\"backbone\",\"name N or name CA or name C or name O\")\n");
	fprintf(fp,"cmd.select(\"sidechain\",\"polymer and not backbone and not hydro\")\n");
#else
	fprintf(fp,"cmd.select(\"backbone\",\"name N or name CA or name C or name O ");
	   fprintf(fp,"or name H or name HA2 or name HA3 or name H1 or name H2 or name H3\")\n");
	fprintf(fp,"cmd.select(\"sidechain\",\"polymer and not backbone\")\n");
#endif

	fprintf(fp,"cmd.color(\"gray\",\"all\")\n");
	fprintf(fp,"cmd.button(\"L\",\"Ctrl\",\"RotZ\")\n");
	fprintf(fp,"cmd.button(\"R\",\"Ctrl\",\"Move\")\n");
	fprintf(fp,"cmd.set_color(\"gray0\",\"[170,170,170]\")\n");
	fprintf(fp,"cmd.set_color(\"blue0\",\"[200,200,255]\")\n");
	fprintf(fp,"cmd.hide(\"everything\")\n");
	fprintf(fp,"cmd.select(\"metals\",\"inorganic and not hydro and not solvent\")\n");
	fprintf(fp,"cmd.select(\"ligand\",\"hetatm and not hydro and not solvent\")\n");

	if(c){
	   fprintf(fp,"cmd.select(\"tmp1\",\"hetatm and not (hydro or solvent)\")\n");
	   fprintf(fp,"cmd.select(\"near_chn\",\"tmp1 and (chain %c around 2.5)\")\n",c);
	   // fprintf(fp,"cmd.select(\"ligand\",\"bymolecule near_chn\")\n",c);
	   // fprintf(fp,"cmd.select(\"ligand\",\"byfragment near_chn\")\n",c);
	   fprintf(fp,"cmd.select(\"ligand\",\"hetatm and byobject near_chn\")\n",c);
	   fprintf(fp,"cmd.delete(\"tmp1\")\n");
	   fprintf(fp,"cmd.delete(\"near_chn\")\n");
	}


	fprintf(fp,"cmd.show(\"sticks\",\"ligand and not hydro\")\n");
	fprintf(fp,"cmd.color(\"deepteal\",\"ligand\")\n");

	fprintf(fp,"cmd.set_color(\"magenta1\",\"[255,100,255]\")\n");
	fprintf(fp,"cmd.set_color(\"red1\",\"[255,100,100]\")\n");
	fprintf(fp,"cmd.set_color(\"orange1\",\"[255,170,90]\")\n");
	fprintf(fp,"cmd.set_color(\"dkorange\",\"[210,110,10]\")\n");
	fprintf(fp,"cmd.set_color(\"yellow1\",\"[255,250,80]\")\n");
	fprintf(fp,"cmd.set_color(\"dkyellow\",\"[180,180,0]\")\n");
	fprintf(fp,"cmd.set_color(\"green1\",\"[120,255,110]\")\n");
	fprintf(fp,"cmd.set_color(\"cyan1\",\"[150,255,255]\")\n");
	fprintf(fp,"cmd.set_color(\"blue1\",\"[170,170,255]\")\n");
	fprintf(fp,"cmd.set_color(\"purple1\",\"[190,120,255]\")\n");

	if(c){ 
	   fprintf(fp,"cmd.show(\"cartoon\",\"chain %c\")\n",c);
	   fprintf(fp,"cmd.color(\"blue0\",\"chain %c\")\n",c);
	   // fprintf(fp,"cmd.cartoon(\"loop\",\"chain %c\")\n",c);
	   fprintf(fp,"cmd.set(\"cartoon_oval_length\",\"1.1\",\"chain %c\")\n",c);
	   fprintf(fp,"cmd.set(\"cartoon_oval_width\",\"0.25\",\"chain %c\")\n",c);
	   fprintf(fp,"cmd.center(\"chain %c\")\n",c);
	} else {
	   fprintf(fp,"cmd.show(\"cartoon\",\"polymer\")\n");
	   fprintf(fp,"cmd.color(\"blue0\",\"polymer\")\n");
	   fprintf(fp,"cmd.set(\"cartoon_oval_length\",\"1.1\",\"all\")\n");
	   fprintf(fp,"cmd.set(\"cartoon_oval_width\",\"0.25\",\"all\")\n");
	   // fprintf(fp,"cmd.cartoon(\"loop\",\"polymer\")\n");
	} // fprintf(fp,"cmd.set(\"seq_view\")\n");
	
}

void	pml_typ::PrintTAIL(FILE *fp)
{ 
	char	c,i;
	const char clss[] = " YROMGCBPADEFHIJKLNQSTUVWXZ ";
	for(i=1,c=clss[i]; isupper(c); i++,c=clss[i]){
// fp=stderr; fprintf(stderr,"c=%c; i=%i\n",c,i);
	    if(head[c]){
		assert(isprint(c));
		fprintf(fp,"cmd.create(\"Class_%c\",\"",c);
		head[c]->Put(fp," or "); fprintf(fp,"\")\n");
		fprintf(fp,"cmd.set(\"seq_view\",0,\"Class_%c\")\n",c);
		fprintf(fp,"cmd.order(\"");
		head[c]->Put(fp," "); fprintf(fp,"\",\"no\",\"bottom\")\n");
	    }
	}
	fprintf(fp,"cmd.set(\"dot_radius\",\"0\",\"all\")\n");
	fprintf(fp,"cmd.set(\"dot_color\",\"gray0\",\"all\")\n");
#if 1
	fprintf(fp,"cmd.set(\"sphere_scale\",\"0.22\",\"all\")\n");
	fprintf(fp,"cmd.set(\"stick_radius\",\"0.20\",\"all\")\n");
#else
	fprintf(fp,"cmd.set(\"sphere_scale\",\"0.33\",\"all\")\n");
	fprintf(fp,"cmd.set(\"stick_radius\",\"0.30\",\"all\")\n");
#endif
	fprintf(fp,"cmd.set(\"dot_width\",\"1\",\"all\")\n");
	fprintf(fp,"cmd.set(\"dot_density\",\"4\",\"all\")\n");
	// fprintf(fp,"cmd.set_bond(\"stick_transparency\",\"0.3\",\"all\")\n");
	fprintf(fp,"cmd.set_bond(\"stick_transparency\",\"0.0\",\"all\")\n");

	fprintf(fp,"cmd.deselect( )\n"); 
	fprintf(fp,"cmd.set(\"ray_shadows\",\"off\")\n"); 
	fprintf(fp,"cmd.zoom( )\n\n"); fflush(fp); 
}

void	pml_typ::PrintColor(FILE *fp, char c)
{
	short	R,G,B;
	// fprintf(stderr,"c = %c\n",c);
   // if(c == 'X' || c == 'x') fprintf(fp,"color cpk\n");
   {
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
	} fprintf(fp,"\"[%d,%d,%d]\"",R,G,B);
   }
}

void	pml_typ::FixAtom(char *atom)
// bug in rasmol does not allow strings like: "ANP1413B.*ho3"
// This routine fixes this...
{ 
  Int4	i;
  for(i=0; atom[i]; i++) if(atom[i]=='*') atom[i]='?'; 
  i--; if(atom[i] == '?') atom[i]='*'; // asteric on end seems to be okay.
}

void	pml_typ::ColorAtom(FILE *fp, char *mol, char c, char color, Int4 i, char *atom)
{
	char	atm,str[200];
	char	*clr,cpk[][10]= {" ","blue","red","white","yellow" };

	if(atom){
	  for(Int4 j=0; !isalpha(atm=atom[j]); j++){
		if(atm==0) print_error("ColorAtom( ) syntax error");
	  }
	} else atm='x';
#if 0
	if(c){
	 	if(atom) fprintf(fp,"cmd.color(\"select %s%d%c.%s\n",mol,i,c,atom);
	 	else fprintf(fp,"select %s%d%c\n",mol,i,c);
	} else {
		if(atom) fprintf(fp,"select %s%d.%s\n",mol,i,atom);
		else fprintf(fp,"select %s%d\n",mol,i);
	}
#endif
	// if(color == 'X')
	{
	   switch(atm){
	    // case 'n': fprintf(fp,"color [143,143,255]\n"); break;
	    case 'n': clr=cpk[1]; break;
	    case 'o': clr=cpk[2]; break;
	    case 'h': clr=cpk[3]; break;
	    case 's': clr=cpk[4]; break;
	    case 'c': clr=0; break;
	    default: clr=0; break;
	   }
	   if(clr==0) return;	// don't color carbon atoms...
	   BooLean found=FALSE;
	   if(IsAARes(mol) && !IsBackBoneAtom(atom) && this->IsObject(i,c)){
	     for(char cls='A'; cls <= 'Z'; cls++){
	       if(Object[c][cls] && MemberSet(i,Object[c][cls])){ 
	         if(c) sprintf(str,"%s%i%c%c",mol,i,c,cls); else sprintf(str,"%s%i",mol,i);
	     	 fprintf(fp,"cmd.color(\"%s\",\"%s and name %s\")\n",clr,str,atom);
	       }
	     }
	   } else if(c) {
	 	fprintf(fp,"cmd.color(\"%s\",\"resn %s and resi %d and chain %c and name %s\")\n",
			clr,mol,i,c,atom);
	   } else {
	 	fprintf(fp,"cmd.color(\"%s\",\"resn %s and resi %d and name %s\")\n",clr,mol,i,atom);
	   }
	}
}

void	pml_typ::print_bond_only(FILE *fp,char res[4],Int4 i,char c,Int4 width,char color,
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

void    pml_typ::PrintResColor(FILE *fp,char *selection,char color)
{
	char	*clr;
#if 0
	char	Colors[][12]={" ", "magenta1","red1","orange1","yellow1","green1","cyan1","blue1","purple1",
					"white","gray" };
	char	DkColors[][16]={" ", "magenta","tv_red","tv_orange","tv_yellow","tv_green","cyan","tv_blue",
					"violetpurple","white","gray" };
#elif 0
	char	Colors[][16]={" ", "magenta","tv_red","tv_orange","tv_yellow","tv_green","cyan","tv_blue",
					"violetpurple","gray80","white" };
	char	DkColors[][16]={" ", "magenta","firebrick","dkorange","dkyellow","splitpea","dealteal","deepblue",
					"violetpurple","gray90","gray50" };
#else
	char	Colors[][16]={" ", "magenta1","red1","orange1","tv_yellow","green1","cyan","blue1",
					"purple1","gray80","white" };
	char	DkColors[][16]={" ", "magenta","firebrick","dkorange","dkyellow","splitpea","dealteal","deepblue",
					"violetpurple","gray90","gray50" };
#endif
	switch (color){
	   case 'm': clr=DkColors[1]; break;
	   case 'r': clr=DkColors[2]; break;
	   case 'o': clr=DkColors[3]; break;
	   case 'y': clr=DkColors[4]; break;
	   case 'g': clr=DkColors[5]; break;
	   case 'c': clr=DkColors[6]; break;
	   case 'b': clr=DkColors[7]; break;
	   case 'p': clr=DkColors[8]; break;
	   case 'w': clr=DkColors[9]; break;

	   case 'M': clr=Colors[1]; break;
	   case 'R': clr=Colors[2]; break;
	   case 'O': clr=Colors[3]; break;
	   case 'Y': clr=Colors[4]; break;
	   case 'G': clr=Colors[5]; break;
	   case 'C': clr=Colors[6]; break;
	   case 'B': clr=Colors[7]; break;
	   case 'P': clr=Colors[8]; break;
	   case 'L': clr=Colors[9]; break;
	   case 'W': clr=Colors[10]; break;
	   default: clr=Colors[9]; break;
	} fprintf(fp,"cmd.color(\"%s\",\"%s\")\n",clr,selection);
}

void	pml_typ::DefineResidueItem(FILE *fp,char aa,Int4 i,char c,char color)
{
	Int4	r = AlphaCode(aa,AB);
	char	str[200],clr=toupper(color);
	// if(c){ sprintf(str,"%s%d%c",AARes[r],i,c); } else { sprintf(str,"%s%d",AARes[r],i); } 
	if(c){ sprintf(str,"%s%d%c%c",AARes[r],i,c,clr); } else { sprintf(str,"%s%d",AARes[r],i); } 
        if(Object[c][clr]==0){ Object[c][clr]=MakeSet(SetSize);  ClearSet(Object[c][clr]); }
#if 1	// avoid duplicate definitions....
	if(MemberSet(i,Object[c][clr])) return; // don't create object again...(due to redundancy in VSI file).
#endif
// if(clr=='Y'){ fprintf(stderr,"--> %s\n",str); }
	AddSet(i,Object[c][clr]);
	if(head[clr]==0){ head[clr] = last[clr] = new obj_typ(str); }
	else { last[clr]->nxt = new obj_typ(str); last[clr] = last[clr]->nxt; }

	fprintf(fp,"cmd.create(\"%s\",\"resn %s and resi %d and chain %c\")\n",str,AARes[r],i,c); 
	if(aa == 'G'){
	   fprintf(fp,"cmd.show(\"spheres\",\"%s and name ca\")\n",str); 
	} else if(aa == 'P'){
	   fprintf(fp,"cmd.show(\"sticks\",\"%s and not (name c or name o or hydro)\")\n",str); 
	} else fprintf(fp,"cmd.show(\"sticks\",\"%s and not (name c or name n or name o or hydro)\")\n",str);
	PrintResColor(fp,str,color); 
	fprintf(fp,"cmd.set(\"seq_view\",0,\"%s\")\n",str);
	fprintf(fp,"cmd.disable(\"%s\")\n",str);
#if 0
	if(aa=='G'){
	  if(c){
		fprintf(fp,"cmd.select(\"gly%d%c\",\"resn gly and resi %d and chain %c",i,c,i,c); 
		sprintf(str,"gly%d%c",i,c); 
		fprintf(fp," and  not (name O or name C or name N or hydro\")\n"); 
	  } else {
		fprintf(fp,"cmd.select(\"gly%d\",\"resn gly and resi %d",i,i); sprintf(str,"gly%d",i);
		fprintf(fp," and  not (name O or name C or name N or hydro\")\n"); 
	  }
		sprintf(str,"gly%d%c",i,c); 
		fprintf(fp," and  not (name O or name C or name N or hydro\")\n"); 
	  PrintResColor(fp,str,color);
	} else {
	  r = AlphaCode(aa,AB);
	  if(aa=='P'){
	    if(c){ 
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
	    }
	    if(isupper(color)) fprintf(fp,"wireframe %d\n",width);
	    else fprintf(fp,"wireframe %d\n",width);
	    if(sideonly && aa != 'P') 
		fprintf(fp,"select %s%d%c and sidechain and not hydrogen\n",AARes[r],i,c);
	    } else { // c == 0..
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
#endif
}

void    pml_typ::PutResidueItem(FILE *fp,char a,Int4 i,char *atom1,char *atom2,char chn,char color,
		Int4 wire_width,Int4 spacefill)
{
    if(atom1){
        if(atom2) PrintResidueAtoms(fp,a,i,atom1,atom2,chn,color,wire_width,spacefill);
        else PrintResidueAtom(fp,a,chn,color,i,atom1,TRUE,spacefill);
    } 
}

void    pml_typ::PrintResidueAtoms(FILE *fp, char a, Int4 i, char *atom1, char *atom2,
	char chn, char color,Int4 wire_width,Int4 spacefill)
{
        char *mol= AARes[AlphaCode(a,AB)];
	PrintMoleculeAtoms(fp, mol, i, atom1, atom2,chn, color,wire_width,spacefill);
}

void    pml_typ::PutMoleculeItem(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
                        char chn, char color,Int4 wire_width,Int4 spacefill)
{
    if(atom1){
	if(atom2) PrintMoleculeAtoms(fp,mol, i,atom1,atom2,chn,color,wire_width,spacefill);
        else PrintMoleculeAtom(fp,mol,chn,color,i,atom1,TRUE,spacefill);
    } 
}

void	pml_typ::PrintTrace(FILE *fp, Int4 start, Int4 end, char c,char color,char mode,Int4 diameter)
{
	static Int4	color_ID=0;
	char		str[100];
	sprintf(str,"bb_color%d",color_ID); color_ID++;
	fprintf(fp,"cmd.set_color(\"%s\",",str); PrintColor(fp,color); fprintf(fp,")\n");
	// fprintf(stderr,"i=%d; j=%d; c=%c\n",i,j,color);
#if 0
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
	} else 
#endif
	if(c){
#if 0
	   if(mode=='^') fprintf(fp,"select %d-%d%c\n",start,end,c);
	   else if(mode==':') fprintf(fp,"select %d-%d%c\n",start,end,c);
	   else if(mode=='+') fprintf(fp,"select %d-%d%c and not water\n",start,end,c);
	   else 
#endif
	   fprintf(fp,"cmd.show(\"cartoon\",\"resi %d-%d and chain %c and backbone\")\n",start,end,c);
	   // fprintf(fp,"cmd.cartoon(\"loop\",\"chain %c and backbone\")\n",c);
	   fprintf(fp,"cmd.set(\"cartoon_oval_length\",\"1.1\",\"resi %d-%d and chain %c and backbone\")\n",
				start,end,c);
	   fprintf(fp,"cmd.set(\"cartoon_oval_width\",\"0.25\",\"resi %d-%d and chain %c and backbone\")\n",
				start,end,c);
	   fprintf(fp,"cmd.color(\"%s\",\"resi %d-%d and chain %c and backbone\")\n",str,start,end,c);
	} else { 
#if 0
	   if(mode==':') fprintf(fp,"select %d-%d\n",start,end);
	   if(mode=='+') fprintf(fp,"select %d-%d\n",start,end);
	   else 
#endif
	   fprintf(fp,"cmd.show(\"cartoon\",\"resi %d-%d and backbone\")\n",start,end);
	   // fprintf(fp,"cmd.cartoon(\"loop\",\"backbone\")\n");
	   fprintf(fp,"cmd.set(\"cartoon_oval_length\",\"1.1\",\"resi %d-%d and backbone\")\n",start,end);
	   fprintf(fp,"cmd.set(\"cartoon_oval_width\",\"0.25\",\"resi %d-%d and backbone\")\n",start,end);
	   fprintf(fp,"cmd.color(\"%s\",\"resi %d-%d and backbone\")\n",str,start,end);
	}
#if 0
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
#endif
}

void	pml_typ::PrintView(FILE *fp, Int4 rx,Int4 ry,Int4 rz,Int4 tx,Int4 ty,Int4 tz)
{
        fprintf(fp,"reset\n");
        fprintf(fp,"rotate x %d\n",rx);
        fprintf(fp,"rotate y %d\n",ry);
        fprintf(fp,"rotate z %d\n",rz);
        // fprintf(fp,"translate x %d\ntranslate y %d\n",tx,ty);
        fprintf(fp,"translate x %.2f\ntranslate y %.2f\n",(double)tx/100.0,(double)ty/100.0);
        fprintf(fp,"zoom %d\n",tz);
}

char	*pml_typ::DefineResCloud(FILE *fp, char aa, Int4 site, char *atom1, char *atom2, char chain)
{
	char str[50];
	Int4  r= AlphaCode(aa,AB);

    if(atom1){
      if(atom2){
	if(chain==0) sprintf(str,"resn %s and resi %d and name %s",AARes[r],site,atom2);
	else sprintf(str,"resn %s and resi %d and chain %c and name %s",AARes[r],site,chain,atom2);
      } else {
	if(chain==0) sprintf(str,"%s%d.%s",AARes[r],site,atom1);
	else sprintf(str,"%s%d%c.%s",AARes[r],site,chain,atom1);
      }
    } else {
	return DefineResCloud(fp,aa,site,chain);
    } return AllocString(str);
}

char	*pml_typ::DefineResCloud(FILE *fp, char aa, Int4 site, char chain)
{
	char	str[50],str2[50];
	Int4	r= AlphaCode(aa,AB);

return 0;
	fprintf(fp,"cmd.show(\"dots\","); 
	if(chain==0) {
	   fprintf(fp,"\"resn %s and resi %d and (",AARes[r],site); 
	} else {
	   fprintf(fp,"\"resn %s and resi %d and chain %c and (",AARes[r],site,chain); 
	}
	switch (aa){
	  case 'A': fprintf(fp,"name CB"); break;
	  case 'C': fprintf(fp,"name CB or name SG"); break;
	  case 'I':
	 	fprintf(fp,"name CB or name CD1 or name CG1 or name CG2");
	  break;
	  case 'L':
	 	fprintf(fp,"name CB or name CD1 or name CD2 or name CG");
	  break;
	  case 'V':
	 	fprintf(fp,"name CB or name CG1 or name CG2");
	  break;
	  case 'Y': case 'F': 
	 	fprintf(fp,"name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ");
	  break;
	  case 'R':
	 	fprintf(fp,"name NE or name NH1 or name NH2 or name CZ");
	  break;
	  case 'P': fprintf(fp,"name CB or name CG or name CD"); break;
	  case 'D': fprintf(fp,"name CG or name OD1 or name OD2"); break;
	  case 'E': fprintf(fp,"name CD or name OE1 or name OE2"); break;
	  case 'K':
	 	fprintf(fp,"name CB or name CG or name CD or name CE or name NZ");
	  break;
	  case 'M': fprintf(fp,"name CB or name CG or name SD or name CE"); break;
	  case 'N': fprintf(fp,"name CG or name ND2 or name OD1"); break;
	  case 'Q': fprintf(fp,"name CD or name NE2 or name OE1"); break;
	  case 'H':
	 	fprintf(fp,"name CG or name CE1 or name NE2 or name ND1 or name CD2");
	  break;
	  case 'W':
	 	fprintf(fp,"name CG or name CD1 or name CD2 or name NE1 or name CE2 or name CE3");
	 	fprintf(fp,"or name CZ2 or name CZ3 or name CH2");
	  break;
	  default: break;
	} fprintf(fp,")\")\n"); 
	return 0;
}

void	pml_typ::PrintResidueAtom(FILE *fp, char aa, char c, char color, Int4 i, 
		char *atom, BooLean sideonly,Int4 big_spacefill)
{
	char *mol = AARes[AlphaCode(aa,AB)];
	PrintMoleculeAtom(fp,mol,c,color,i,atom,sideonly,big_spacefill);
}

BooLean pml_typ::IsAARes(char *mol)
{
	static char AARes0[22][4]= {"unk","cys","gly","ala","ser","thr","asn","asp",
		"glu","gln","lys","arg","his","trp","tyr","phe","val","ile","leu","met","pro"};
	for(Int4 i=1; i <= 20; i++){ 
		if(strcmp(mol,AARes0[i]) == 0) return TRUE;
	} return FALSE;
}

BooLean pml_typ::IsBackBoneAtom(char *atm)
{
	static char Atm0[12][4]= {"x","n","o","c","h","ha2","ha3","h1","h2","h3"};
	for(Int4 i=1; i <= 9; i++){
		if(strcmp(atm,Atm0[i]) == 0) return TRUE;
	} return FALSE;
}

void	pml_typ::PrintMoleculeAtom(FILE *fp,char *mol, char chn, char color, Int4 i, 
		char *atom, BooLean sideonly,Int4 spacefill)
{
        // FixAtom(atom);
        char chnX,str[200];
        if(chn==' ') chnX = 0; else chnX = chn;
	// if(IsAARes(mol) && !IsBackBoneAtom(atom) && IsObject(i,chnX))
	if(IsAARes(mol) && this->IsObject(i,chnX))
	{
	  for(char cls='A'; cls <= 'Z'; cls++){
	    if(Object[chnX][cls] && MemberSet(i,Object[chnX][cls])){
              if(chnX) sprintf(str,"%s%d%c%c",mol,i,chnX,cls); else sprintf(str,"%s%d",mol,i);
	      fprintf(fp,"cmd.show(\"spheres\",\"%s and name %s\")\n",str,atom);
	    }
	  }
	} else if(chnX){
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and chain %c and name %s\")\n",
                        mol,i,chnX,atom);
        } else {
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and name %s\")\n",mol,i,atom);
        } ColorAtom(fp,mol,chnX,color,i,atom);
}

void    pml_typ::PrintMoleculeAtoms(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
	char chn, char color,Int4 wire_width,Int4 spacefill)
{
        Int4    r,k;
        char chnX,atm=0,atm2[10],str[200];

        // FixAtom(atom1); FixAtom(atom2);
        for(k=0; !(isalpha(atm=atom1[k]) || atm=='*'); k++){
                if(atm==0) print_error("atom1 syntax error");
        }     
#if 1
	if(strcmp("arg",mol) == 0){
	  if(atom2[0]=='h' && isdigit(atom2[1]) && isdigit(atom2[2])){ sprintf(atm2,"h%s",atom2); }
	  else sprintf(atm2,"%s",atom2);
	} else sprintf(atm2,"%s",atom2);
#endif
        if(chn==' ') chnX = 0;   else chnX = chn;              

#if 0
	if(IsAARes(mol) && !IsBackBoneAtom(atom1) && !IsBackBoneAtom(atm2) && this->IsObject(i,chnX))
#else
	if(IsAARes(mol) && this->IsObject(i,chnX)) 
#endif
	{
	  for(char cls='A'; cls <= 'Z'; cls++){
	    if(Object[chnX][cls] && MemberSet(i,Object[chnX][cls])){
               if(chnX) sprintf(str,"%s%d%c%c",mol,i,chnX,cls); else sprintf(str,"%s%d",mol,i);
	       fprintf(fp,"cmd.show(\"sticks\",\"%s and (name %s or name %s)\")\n",str,atom1,atm2);
	       if(atm=='c'){
		fprintf(fp,"cmd.show(\"spheres\",\"%s and name %s\")\n",str,atm2);
	       } else {
		fprintf(fp,"cmd.show(\"spheres\",\"%s and (name %s or name %s)\")\n",str,atom1,atm2);
	       }
	    }
	  }
	} else if(chnX){
           fprintf(fp,"cmd.show(\"sticks\",\"resn %s and resi %d and chain %c and (name %s or name %s)\")\n",
                                mol,i,chnX,atom1,atm2);
           if(atm=='c'){
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and chain %c and name %s\")\n",
                        mol,i,chnX,atm2);
           } else {
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and chain %c and (name %s or name %s)\")\n",
                        mol,i,chnX,atom1,atm2);
           }
        } else {
           fprintf(fp,"cmd.show(\"sticks\",\"resn %s and resi %d and (name %s or name %s)\")\n",
                                mol,i,atom1,atm2);
           if(atm=='c'){
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and name %s\")\n",
                                mol,i,atm2);
           } else {
                fprintf(fp,"cmd.show(\"spheres\",\"resn %s and resi %d and (name %s or name %s)\")\n",
                                mol,i,atom1,atm2);
           }
        }
	if(atm != 'c') ColorAtom(fp,mol,chnX,color,i,atom1); 
	ColorAtom(fp,mol,chnX,color,i,atm2);
}       

void	pml_typ::PrintCommand(FILE *fp,char *cmd,char type)
{
        switch (type){
          case 'E': fprintf(fp,"%s\n",cmd); break;
        }
}


