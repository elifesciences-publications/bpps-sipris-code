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

#include "pda_typ.h"

short	GetShapeCMSA2RTF(char shape) 
// return the rtf code for the input shape.
{
	switch (shape){
	   case 'R': return 1;	// pointed rectangle*
	   case 'r': return 2;  // rounded rectangle.
	   case 'E': return 3;  // elipse*
	   case 'D': return 4;  // diamond*.
	   case 'T': return 5;	// triangle*
	   case ' ': return 6;	// ???
	   case 'I': return 7;  // italic rectangle*
	   case 't': return 8;	// trapezoid*
	   case 'H': return 9;	// hexagon*
	   case 'O': return 10;	// octagon*
	   case 'p': return 12;	// 5-pointed star
	   case 'P': return 56;	// pentagon*
	   case 'N': return 57;	// No symbol
	   case 's': return 60;	// 32star*
	   case 'e': return 71;	// explosion*
	   case 'C': return 106;	// cloud*
	   case 'S': return 183;	// Sun*
	   case 'f': return 187;	// four-pointed start
	   default: print_error("GetShapeCMSA2RTF( ) input error");
	}
}

void	GetColorAndShadeCMSA2RTF(char FillColor, char MainColor,
	UInt4 *fillcolor, UInt4 *lineshade, double Evalue)
/**************************** Codes used: *****************************
  FillColors:
         'M' = magenta    'R' = red    'Y' = yellow         'W' = White
         'G' = green      'C' = cyan   'V' = violet (blue)  'B' = Black
 **********************************************************************/
{
	const Int4 	num_shades=5;
	enum shade_type { ghost=0, light, moderate, dark, black};
	Int4		shade;
	UInt4 	grayshade[num_shades]={13158600,10526880,7895160,5592405,0};
	       // [200,200,200]; [160,160,160]; [120,120,120]; [85,85,85]; [0,0,0].
	UInt4 color[num_shades];
	unsigned char set_dark[num_shades] = {215, 175, 125, 75, 0};
	unsigned char set_light[num_shades] = {255, 255, 255, 230, 180};
	BooLean color_on[4]; 

    // 0. Set up color scheme.
	if(islower(FillColor)) FillColor = toupper(FillColor);
	if(FillColor == MainColor) FillColor='W'; // use white for conflicts.
	enum primary_color { blue=1, green, red };
	color_on[0] = FALSE;
	switch (FillColor){	// color_on[4] = [ none, Blue, Green, Red ]; 
	  case 'M': color_on[red] = color_on[blue] = TRUE;
		color_on[green] = FALSE; break;
	  case 'R': color_on[green] = color_on[blue] = FALSE;
		color_on[red] = TRUE; break;
	  case 'Y': color_on[green] = color_on[red] = TRUE;
		color_on[blue] = FALSE; break;
	  case 'G': color_on[red] = color_on[blue] = FALSE;
		color_on[green] = TRUE; break;
	  case 'C': color_on[green] = color_on[blue] = TRUE;
		color_on[red] = FALSE; break;
	  case 'V': color_on[green] = color_on[red] = FALSE;
		color_on[blue] = TRUE; break;
	  case 'W': color_on[green] = color_on[red] = color_on[blue] = TRUE;
		break;
	  case 'B': color_on[green] = color_on[red] = color_on[blue] = FALSE;
		break;
	  default: print_error("PutSchematicPDA( ) FillColor input error");
	}
	unsigned char 	*ptr2color;
	for(shade=ghost; shade < num_shades; shade++){
	   ptr2color = (unsigned char*) &color[shade]; ptr2color[0] = 0;
	   for(Int4 pc=blue; pc <= red; pc++){
		if(color_on[pc]) ptr2color[pc] = set_light[shade]; 
		else ptr2color[pc] = set_dark[shade]; 
	   }
	}
	// 1. Get shade based on Evalue.
	double e = Evalue;
	if(e <= 0.0) shade=ghost; else if(e < 1.5) shade=light;
	else if(e < 3.0) shade=moderate; else if(e < 4.5) shade=dark;
	else shade = black;
	if(FillColor=='W') *fillcolor = 16777215; // = [255,255,255]
	else *fillcolor= color[shade];
	*lineshade= grayshade[shade];
}

static void put_rtf_shape_pda(FILE *fp,BooLean group,Int4 start,Int4 end,
	short shape, double scale, UInt4 fillcolor, 
	UInt4 lineshade, UInt4 linewidth, char mode,
	char *TheName)
{
	UInt4	top,bottom;
	char		*Name,*str,name[20];
	Int4		max_str_len=9;
	if(TheName && mode != 'n'){
	  Name=TheName;
	  strncpy(name,Name,max_str_len); name[max_str_len]=0;
	  for(str=name; isalnum(str[0]); str++) ;
	  if(str[0]) str[0]=0;
	} else Name=0;
	switch (mode){
#if 0
	  case 'n': top=75; bottom=175; break;
	  case 'N': top=50; bottom=200; break;
	  default: top=0; bottom=250; break;
#else
	  case 'n': top=50; bottom=200; break;
	  case 'N': top=0; bottom=250; break;
	  default: top=0; bottom=250; break;
#endif
	}
	if(shape < 60 && Name) fprintf(fp,"{\\shp{"); 
	else fprintf(fp,"{\\shp");
	if(group){
	     fprintf(fp,"{\\*\\shpinst\n");
	     fprintf(fp,"{\\sp{\\sn relLeft}{\\sv %d}}\n",
		(Int4)(start*scale)+100);
	     fprintf(fp,"{\\sp{\\sn relTop}{\\sv %u}}\n",top);
	     fprintf(fp,"{\\sp{\\sn relRight}{\\sv %d}}\n",
		(Int4)(end*scale) +100);
	     fprintf(fp,"{\\sp{\\sn relBottom}{\\sv %u}}\n",bottom);
	     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",shape);
	} else {
	     fprintf(fp,"{\\*\\shpinst\\shpleft%d\\shptop%d\\shpright%d\\shpbottom%d",
		(Int4)(start*scale)+100,top,(Int4)(end*scale) +100,bottom);
	     fprintf(fp,"\\shpfhdr0\\shpbxcolumn\\shpbypara\\shpwr3\\shpwrk0");
	     fprintf(fp,"\\shpfblwtxt0\\shpz0\\shplid1026\n");
	     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",shape);
	}
   if(shape < 60){	// special effects
	// fprintf(fp,"{\\sp{\\sn lTxid}{\\sv 196608}}\n"); // txt???
	if(Name){
	  fprintf(fp,"{\\sp{\\sn dxTextLeft}{\\sv 0}}\n");
	  fprintf(fp,"{\\sp{\\sn dxTextRight}{\\sv 0}}\n");
	  fprintf(fp,"{\\sp{\\sn dyTextTop}{\\sv 0.05}}\n");
	  fprintf(fp,"{\\sp{\\sn dyTextBottom}{\\sv 0}}\n");
	}
      // if(mode=='n' && shape==1) // use for coiled coils...
      if(mode=='n'){
	if(strcmp(TheName,"TM_helix") == 0){ 	// transmembrane helices...
	   lineshade=32896; 
	} else {				// coiled coils
	   fillcolor=lineshade=7895160;	// gray shade (see above).
	}
	fprintf(fp,"{\\sp{\\sn fillColor}{\\sv %u}}\n",fillcolor);
	fprintf(fp,"{\\sp{\\sn fillOpacity}{\\sv %u}}\n",32768);
	fprintf(fp,"{\\sp{\\sn fFilled}{\\sv 1}}\n");
      } else {
	fprintf(fp,"{\\sp{\\sn fillBackColor}{\\sv %u}}\n",271450864); 
	fprintf(fp,"{\\sp{\\sn fillType}{\\sv 7}}\n");
	fprintf(fp,"{\\sp{\\sn fillColor}{\\sv %u}}\n",fillcolor);
	fprintf(fp,"{\\sp{\\sn fillBackColor}{\\sv %u}}\n",271450864);
	fprintf(fp,"{\\sp{\\sn fillFocus}{\\sv %d}}\n",50);
	fprintf(fp,"{\\sp{\\sn fillShadeType}{\\sv 1073741835}}\n");
	fprintf(fp,"{\\sp{\\sn fFilled}{\\sv 1}}\n");
	fprintf(fp,"{\\sp{\\sn fillShape}{\\sv 1}}\n");
      }
   } else { fprintf(fp,"{\\sp{\\sn fillColor}{\\sv %u}}\n",fillcolor); }
	fprintf(fp,"{\\sp{\\sn lineColor}{\\sv %u}}\n",lineshade); // dark gray
	fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv %u}}\n",linewidth);
	if(shape < 60 && Name){
// 1 black; 2 blue; 3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white.
// 9 dk blue; 10 teal; 11 dk green; 12 violet; 13 dk red; 14 dk yellow;
// 15 dk grey; 16 lt grey...
	  fprintf(fp,"{\\sp{\\sn fLayoutInCell}{\\sv 1}}\n");
#if 0
	  fprintf(fp,"{\\shptxt \\pard\\plain \\s2\\qc {\\b\\impr\\f2\\fs16\\cf18 %s \\par }}\n",
		name);
#else 
	  fprintf(fp,"{\\shptxt \\pard\\plain \\s2\\qc {\\b\\impr0\\f2\\fs16\\cf1 %s \\par }}\n",
		name);
#endif
	}
	if(shape < 60 && Name) fprintf(fp,"}\n");
#if 0	// NOT SURE WHY THESE LINES ARE GENERATED BY MSWord; they seem unnecessary.
	fprintf(fp,"{\\shprslt{\\*\\do\\dobxcolumn\\dobypara\\dodhgt8195\\dptxbx\\dptxlrtb\n");
	fprintf(fp,"{\\dptxbxtext\\pard\\plain \\s2\\qc \\li0\\ri0\\keepn\\widctlpar\\aspalpha\\aspnum\\faauto\\outlinelevel2\\adjustright\\rin0\\lin0\\itap0 \\b\\impr\\f4\\fs24\\cf3\\lang1033\\langfe1033\\cgrid\\langnp1033\\langfenp1033 {\\cf13 ARM \\par }}\n");
	fprintf(fp,"\\dpx4500\\dpy264\\dpxsize1620\\dpysize720\\dpfillfgcr240\\dpfillfgcg2\\dpfillfgcb31\\dpfillbgcr255\\dpfillbgcg0\\dpfillbgcb0\\dpfillpat8\\dplinew15\\dplinecor0\\dplinecog0\\dplinecob0}}\n");
#endif
	fprintf(fp,"}}\n");
}

static void put_rtf_shape_pda(FILE *fp,BooLean group,Int4 start,Int4 end,
	short shape, double scale, UInt4 fillcolor, 
	UInt4 lineshade, UInt4 linewidth,char *name)
{ put_rtf_shape_pda(fp,group,start,end,shape,scale,fillcolor,lineshade,linewidth,' ',name); }

static void put_rtf_line_cmsa(FILE *fp, BooLean group, Int4 full_len,
	double scale)
{
   fprintf(fp,"{\\lang1024\n");
   if(group){
     fprintf(fp,"{\\shpgrp{\\*\\shpinst\\shpleft100\\shptop125");
     fprintf(fp,"\\shpright%d\\shpbottom125\\shpbypara\n",
           (Int4)(full_len*scale)+100);

     fprintf(fp,"{\\sp{\\sn groupLeft}{\\sv 100}}\n");
     fprintf(fp,"{\\sp{\\sn groupTop}{\\sv 125}}\n");
     fprintf(fp,"{\\sp{\\sn groupRight}{\\sv %d}}\n",
           (Int4)(full_len*scale)+100);
     fprintf(fp,"{\\sp{\\sn groupBottom}{\\sv 125}}\n");

     fprintf(fp,"{\\shp{\\*\\shpinst\n");
     fprintf(fp,"{\\sp{\\sn relLeft}{\\sv 100}}\n");
     fprintf(fp,"{\\sp{\\sn relTop}{\\sv 125}}\n");
     fprintf(fp,"{\\sp{\\sn relRight}{\\sv %d}}\n",
           (Int4)(full_len*scale)+100);
     fprintf(fp,"{\\sp{\\sn relBottom}{\\sv 125}}\n");

     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",1);
     fprintf(fp,"{\\sp{\\sn shapePath}{\\sv 4}}\n");
#if 1
     fprintf(fp,"{\\sp{\\sn lineColor}{\\sv %u}}\n",5592405); // dark gray
#endif
     // fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 31750}}");
     // fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 12700}}");
     fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 25400}}");
     fprintf(fp,"{\\sp{\\sn fLine}{\\sv 1}}\n}}\n");
   } else {   // old ungrouped method (Two lines for Wilcoxon).
     Int4 center = (Int4)((full_len*scale)/2);
     fprintf(fp,"{\\shp{\\*\\shpinst\\shpleft100\\shptop125");
     fprintf(fp,"\\shpright%d\\shpbottom125\\shpbypara\n",
           center+100);
     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",1);
     fprintf(fp,"{\\sp{\\sn shapePath}{\\sv 4}}\n");
     fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 40000}}");
     fprintf(fp,"{\\sp{\\sn fLine}{\\sv 1}}\n}}\n");

     fprintf(fp,"{\\shp{\\*\\shpinst\\shpleft%d\\shptop125",
           center+100);
     fprintf(fp,"\\shpright%d\\shpbottom125\\shpbypara\n",
           (Int4)(full_len*scale)+100);
     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",1);
     fprintf(fp,"{\\sp{\\sn shapePath}{\\sv 4}}\n");
     fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 10000}}");
     fprintf(fp,"{\\sp{\\sn fLine}{\\sv 1}}\n}}\n");
   }
}

static UInt4 get_explosion_color_pda(Int4 Sq, Int4 d, dom_typ *dom,
	ss_type FullSeq)
{
	UInt4 white = 16777215; // = [255,255,255]
	if(!FullSeq) return white;

	Int4	s,e,len,pos,neg,polar,nonpolar;
	char	r,str[10];
	a_type	A=SeqSetA(FullSeq);

	e_type E = SeqSetE(Sq,FullSeq);
	s = dom->SqStart(Sq,d); e = dom->SqEnd(Sq,d); len = e - s + 1;
	assert(s <= e); pos=neg=polar=nonpolar=0;
	for(Int4 i=s; i <= e; i++){
	   sprintf(str,"%c",AlphaChar(ResSeq(i,E),A));
	   if(strstr("KRH",str) != NULL) pos++;
	   else if(strstr("STGNQWY",str) != NULL) polar++;
	   else if(strstr("DE",str) != NULL) neg++;
	   else if(strstr("AILVFMCP",str) != NULL) nonpolar++;
	}
	enum primary_color { blue=1, green, red };
	unsigned char 	*ptr2color;
	unsigned short c;
	UInt4 color;
	ptr2color = (unsigned char*) &color; ptr2color[0] = 0;

	c=(unsigned short)(150.0*(double)(pos+nonpolar)/(double)len);
	c+=100; if(c > 255) c=255; ptr2color[blue] = (unsigned char) c;
	
	c=(unsigned short)(150.0*(double)(polar+nonpolar)/(double)len);
	c+=100; if(c > 255) c=255; ptr2color[green] = (unsigned char) c;
	
	c=(unsigned short)(150.0*(double)(neg+nonpolar)/(double)len);
	c+=100; if(c > 255) c=255; ptr2color[red] = (unsigned char) c;
	
	return color;
}

void	PutSchematicPDA(FILE *fp, char FillColor, char PageSetUp, 
			BooLean group, ss_type FullSeq, dom_typ *dom)
/*************** Codes used: ************************
  FillColor:
         'M' = magenta    'R' = red    'Y' = yellow
         'G' = green      'C' = cyan   'V' = blue

  ShapeType:
            1 = pointed rectangle* ('R')
            2 = rounded rectangle ('r')
            3 = elipse* ('E')
            4 = diamond* ('D')
            5 = triangle* ('T')
            7 = italic rectangle* ('i')
            8 = trapezoid* ('t')
            9 = hexagon* ('H')
           10 = octagon* ('O')
           20 = line
	   56 = pentagon ('P')
           60 = 32star* ('s')
           71 = explosion* ('e')
           106 = cloud* ('C')
           183 = Sun* ('S')
   maximum points = 8640
  ***************************************************/
{
	a_type	A=SeqSetA(FullSeq);
	e_type	E,E0;
	Int4	J,m,x,*s,pos[5],start,end,repeats;
	Int4	number=NSeqsSeqSet(FullSeq),full_len,maxfull;
	char	id0[202],id[202];
	double	scale,page_width=7000.0;
	char	mode,DomShape;
	short	shape=10,domshape;
	UInt4 white = 16777215; // = [255,255,255]
	UInt4 orange = 3381759; // = [255,153,51]

    // 0. Set up shape & color scheme.
	if(islower(FillColor)) FillColor = toupper(FillColor);
	UInt4 thickline = 19050,mediumline=9525,thinline = 3175;
	UInt4 fillcolor,lineshade,linewidth=thickline;
    // 1. Put header file:
	fprintf(fp,"{\\rtf1\\ansi \\deff4\\deflang1033\n");
	fprintf(fp,"{\\fonttbl {\\f2\\fswiss\\fcharset0\\fprq2 Arial Narrow;}}\n");
	switch(PageSetUp){
	  case 'l': // 8.5" x 11" landscape.
	   fprintf(fp,"\\paperw15840\\paperh12240\\margl720\\margr720");
    	   fprintf(fp,"\\margt1080\\margb1080 ");
    	   fprintf(fp,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   fprintf(fp,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
    	   page_width=12000.0;
	   break;
	  case 'L': // 11" x 17" landscape
	   print_error("PutSchematicPDA( ) 11\" x 17\" landscape not yet implemented");
	   break;
	  case 'p': // 8.5" x 11" portrait.
	   page_width=7000.0;
	   break;
	  case 'P': // 11" x 17" portrait.
	    fprintf(fp,"\\paperw15840\\paperh24480\\margl720\\margr720");
            fprintf(fp,"\\margt1080\\margb1080 ");
            fprintf(fp,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
            fprintf(fp,"\\linex0\\endnhere\\sectdefaultcl \n");
            page_width=12000.0;
	   break;
	  default: print_error("PutSchematicPDA( ) PageSetUp input error");
	}
	fprintf(fp,"\\pard \\nowidctlpar\\tx%d\\adjustright \\sl200\\slmult1 ",
		(Int4)page_width+500); // compressed...

    // 2. Find scaling factor
    	for(maxfull=0,J=1; J <= number; J++) {
	   full_len = SqLenSeqSet(J,FullSeq);
	   maxfull = MAXIMUM(Int4,maxfull,full_len);
	} scale = page_width/(double)maxfull;
    // 3. Loop over all sequences in aligment.
	Int4	lastJ=0;
    	for(id0[0]=0,repeats=0,J=1; J <= number; J++) {
	   E = SeqSetE(J,FullSeq);
	   if(lastJ==0) E0 = E;
	   full_len = SqLenSeqSet(J,FullSeq);
	   StrSeqID(id,200,E);
	   if(strncmp(id,id0,200) != 0){ // New sequence; output line.
		if(lastJ > 0){	// then need to output end of last seq.
	           if(dom){
	             Int4 Sq = lastJ; 
	             for(Int4 d=1; d <= dom->NumDomains(Sq); d++){
			DomShape = dom->Shape(Sq,d);
			if(DomShape == 'e'){ // i.e., biased region.
			   fillcolor = get_explosion_color_pda(Sq,d,dom,FullSeq);
		   	   linewidth = thinline; 
			   lineshade=0;
			} else {
			 GetColorAndShadeCMSA2RTF(dom->Color(Sq,d),FillColor,
				&fillcolor, &lineshade,dom->Evalue(Sq,d));
			 linewidth=mediumline;
			}
			domshape = GetShapeCMSA2RTF(DomShape);
			if(DomShape == 'r') mode = 'n'; else mode = 'N';
// fprintf(stderr,"DomainName(%d,%d)=%s\n",Sq,d,dom->DomainName(Sq,d));
		        put_rtf_shape_pda(fp,group,dom->SqStart(Sq,d),
				dom->SqEnd(Sq,d), domshape,
				scale,fillcolor,lineshade,linewidth,mode,
	             	   	dom->DomainName(Sq,d));
	             }
	           }
		   if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		   else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		   // fprintf(fp,"%s",id0); // old
		   PutShortSeqID(fp,E0);
		   // PutShortSeqID(stderr,E0);
		   // fprintf(fp," (%d)\\par\\par}\n",repeats);
		   fprintf(fp," \\par\\par}\n");
		   E0 = E; repeats=0;
		} lastJ=J;
		// Find length of full sequence
		put_rtf_line_cmsa(fp,group,full_len,scale);
		strcpy(id0,id);
	   }
	}
	if(lastJ > 0){	// then need to output end of last seq.
	    if(dom){
	         Int4 Sq = lastJ;
	         for(Int4 d=1; d <= dom->NumDomains(Sq); d++){
			lineshade=0;
			DomShape = dom->Shape(Sq,d);
			if(DomShape == 'e'){ // i.e., biased region.
			  fillcolor = get_explosion_color_pda(Sq,d,dom,FullSeq);
		   	  linewidth = thinline; 
			} else {
			  GetColorAndShadeCMSA2RTF(dom->Color(Sq,d),FillColor,
				&fillcolor, &lineshade,dom->Evalue(Sq,d));
			  linewidth=mediumline;
			}
			domshape = GetShapeCMSA2RTF(DomShape);
			if(DomShape == 'r') mode = 'n'; else mode = 'N';
// fprintf(stderr,"DomainName(%d,%d)=%s\n",Sq,d,dom->DomainName(Sq,d));
		        put_rtf_shape_pda(fp,group,dom->SqStart(Sq,d),
				dom->SqEnd(Sq,d), domshape,
				scale,fillcolor,lineshade,linewidth,mode,
	             	   	dom->DomainName(Sq,d));
	         }
	    }
	    if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	    else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	    // fprintf(fp,"%s",id0); // old
            PutShortSeqID(fp,E0);
	    // fprintf(fp," (%d)\\par\\par}\n",repeats); E0 = E;

	    if(dom){ // output names of domains...
	     Int4 t;
             // fprintf(fp," (%d)\\par\\par\\par\\par\\par\\par}\n",repeats);
             // fprintf(fp," \\par\\par\\par\\par\\par\\par}\n");
             fprintf(fp," \\par\\par\\par \\page \\par}\n");
	     Int4 maxfull2=0; double scale2,page_width2;
	     for(t=1; t <= dom->NumTypes( ); t++){
	        maxfull2 = MAXIMUM(Int4,maxfull2,dom->AveLen(t)+40);
	     } 
	     page_width2 = page_width *((double)maxfull2/(double) maxfull);
	     scale2 = page_width2/(double)maxfull2;
	     fprintf(fp,"\\pard\\plain \\nowidctlpar\\tx%d\\adjustright \\sl190\\slmult1 ",
		(Int4)page_width2+500); // compressed...
	     for(t=1; t <= dom->NumTypes( ); t++){
		DomShape = dom->Shape(t);
		domshape = GetShapeCMSA2RTF(DomShape);
		if(domshape >= 60) continue;
	        start=1+(maxfull2-dom->AveLen(t))/2;
	        end = start + dom->AveLen(t) -1;
	        put_rtf_line_cmsa(fp,group,maxfull2,scale2);
	        GetColorAndShadeCMSA2RTF(dom->Color(t),FillColor,
			&fillcolor,&lineshade,10.0);
		linewidth=mediumline;
		if(DomShape == 'r') mode = 'n'; else mode = 'N';
// fprintf(stderr,"NameThisType(%d)=%s\n",t,dom->NameThisType(t));
	        put_rtf_shape_pda(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode,
	             	   	dom->NameThisType(t));
	        if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        fprintf(fp,"%s",dom->NameThisType(t)); 
                fprintf(fp," (%d)\\par\\par}\n",dom->NumThisType(t));
	     }
	     Int4 ave=0,total=0;
	     for(t=1; t <= dom->NumTypes( ); t++){
		DomShape = dom->Shape(t);
		domshape = GetShapeCMSA2RTF(DomShape);
		if(domshape < 60) continue;
		else { total++; ave+=dom->AveLen(t); }
	     }
	     if(total > 0 && maxfull2 > 50){	// output biased regions
		ave = (ave/total) + 1;
		domshape = GetShapeCMSA2RTF('e');
		Int4 inc = maxfull2/13;
	   	linewidth = thinline; mode = 'N';
		unsigned char *pc; lineshade=0;
	        enum primary_color { wht=0, blue, green, red };

	        put_rtf_line_cmsa(fp,group,maxfull2,scale2);
		pc = (unsigned char*) &fillcolor; 
		fillcolor = white;

	        start=inc; end = start + 2*inc;
	        put_rtf_shape_pda(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode,0); 
		fillcolor=0; pc[1]=pc[2]=pc[3]=100;
		pc[1] += 150; start=end+inc; end = start + 2*inc;
	        put_rtf_shape_pda(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode,0); 
		fillcolor=0; pc[1]=pc[2]=pc[3]=100;
		pc[2] += 150; start=end+inc; end = start + 2*inc;
	        put_rtf_shape_pda(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode,0);
		fillcolor=0; pc[1]=pc[2]=pc[3]=100;
		pc[3] += 150; start=end+inc; end = start + 2*inc;
	        put_rtf_shape_pda(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode,0);
	        if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		fprintf(fp,"hydrophobic, basic, polar, acidic");
		fprintf(fp," \\par\\par}\n");
	     }
	    } else fprintf(fp," \\par\\par}\n"); E0 = E;
	} fprintf(fp,"}\n");
}

