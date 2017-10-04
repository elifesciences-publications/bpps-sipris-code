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

#include "cmsa.h"

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
         'G' = green      'C' = cyan   'V' = violed (blue)  'B' = Black
 **********************************************************************/
{
	BooLean color_on[4]; 

    // 0. Set up color scheme.
	// enum primary_color { blue=1, green, red };
	enum primary_color { red=0, green, blue};
	// [200,200,200]; [160,160,160]; [120,120,120]; [85,85,85]; [0,0,0].
	if(islower(FillColor)) FillColor = toupper(FillColor);
	if(FillColor == MainColor) FillColor='W'; // use white for conflicts.
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
	  default: print_error("GetColorAndShadeCMSA2RTF( ) FillColor input error");
	}
	const Int4 	num_shades=5;
	Int4		shade;
	UInt4 color[num_shades];
	enum shade_type { ghost=0, light, moderate, dark, black};
	// UInt4 	grayshade[num_shades]={13158600,10526880,7895160,5592405,0};
	// UInt4 	grayshade[num_shades]={0xc8c8c8,0x969696,0x646464,0x323232,0};
	UInt4 	grayshade[num_shades]={0x646464,0x4b4b4b,0x323232,0x191919,0};
#if 1
	unsigned char set_dark[num_shades] = {200, 150, 100, 50, 0};
	// unsigned char set_dark[num_shades] = {215, 175, 125, 75, 0};
	unsigned char set_light[num_shades] = {255, 240, 230, 220, 180};
	//unsigned char set_light[num_shades] = {255, 255, 255, 230, 180};
#else
	unsigned char set_dark[num_shades] = {0xd7, 0xaf, 0x7d, 0x4b, 0};
	unsigned char set_light[num_shades] = {0xff, 0xff, 0xff, 0xe6, 0xb4};
#endif
	for(shade=ghost; shade < num_shades; shade++){
	   color[shade] = 0;
	   UInt4	factor=1;
	   for(Int4 pc=red; pc <= blue; pc++){
		if(color_on[pc]) color[shade] += factor * set_light[shade];
		else color[shade] += factor * set_dark[shade];
		factor = factor * 0x100;	// shift by 
	   }
	}
	// 1. Get shade based on Evalue.
	double e = Evalue;
	if(e <= 0.0) shade=ghost; else if(e < 1.0) shade=light;
	else if(e < 2.0) shade=moderate; else if(e < 3.0) shade=dark;
	else shade = black;
	if(FillColor=='W') *fillcolor = 16777215; // = [0,255,255,255]
	else *fillcolor= color[shade];
	*lineshade= grayshade[shade];
}

static void put_rtf_shape_cmsa(FILE *fp,BooLean group,Int4 start,Int4 end,
	short shape, double scale, UInt4 fillcolor, 
	UInt4 lineshade, UInt4 linewidth, char mode)
{
	UInt4 top,bottom;
	switch (mode){
	  case 'n': top=75; bottom=175; break;
	  case 'N': top=50; bottom=200; break;
	  default: top=0; bottom=250; break;
	}
	if(group){
	     fprintf(fp,"{\\shp");
	     fprintf(fp,"{\\*\\shpinst\n");
	     fprintf(fp,"{\\sp{\\sn relLeft}{\\sv %d}}\n",
		(Int4)(start*scale)+100);
	     fprintf(fp,"{\\sp{\\sn relTop}{\\sv %u}}\n",top);
	     fprintf(fp,"{\\sp{\\sn relRight}{\\sv %d}}\n",
		(Int4)(end*scale) +100);
	     fprintf(fp,"{\\sp{\\sn relBottom}{\\sv %u}}\n",bottom);
	     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",shape);
	} else {
	     fprintf(fp,"{\\shp");
	     fprintf(fp,"{\\*\\shpinst\\shpleft%d\\shptop%d\\shpright%d\\shpbottom%d",
		(Int4)(start*scale)+100,top,(Int4)(end*scale) +100,bottom);
	     fprintf(fp,"\\shpfhdr0\\shpbxcolumn\\shpbypara\\shpwr3\\shpwrk0");
	     fprintf(fp,"\\shpfblwtxt0\\shpz0\\shplid1026\n");
	     fprintf(fp,"{\\sp{\\sn shapeType}{\\sv %d}}",shape);
	}
	fprintf(fp,"{\\sp{\\sn fillColor}{\\sv %u}}\n",fillcolor);
	fprintf(fp,"{\\sp{\\sn lineColor}{\\sv %u}}\n",lineshade); // dark gray
	fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv %u}}}}\n",linewidth);
}

static void put_rtf_shape_cmsa(FILE *fp,BooLean group,Int4 start,Int4 end,
	short shape, double scale, UInt4 fillcolor, 
	UInt4 lineshade, UInt4 linewidth)
{ put_rtf_shape_cmsa(fp,group,start,end,shape,scale,fillcolor,lineshade,linewidth,' '); }

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
     fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 31750}}");
     // fprintf(fp,"{\\sp{\\sn lineWidth}{\\sv 12700}}");
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

static UInt4 get_explosion_color_cmsa_rtf(Int4 Sq, Int4 d, cma_typ cma)
{
	UInt4 white = 16777215; // = [255,255,255]
	if(!cma->FullSeq) return white;

	dom_typ	*dom=DomainsCMSA(cma);
	Int4	s,e,len,pos,neg,polar,nonpolar;
	char	r,str[10];
	a_type	A=SeqSetA(cma->FullSeq);

	e_type E = SeqSetE(Sq,cma->FullSeq);
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

void	PutSchematicCMSA(FILE *fp, char FillColor, unsigned char ShapeType, 
	char PageSetUp, BooLean group, cma_typ cma, ss_type keydata)
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
	gss_typ	*gss=gssCMSA(cma);
	st_type	S=SitesCMSA(cma);
	ss_type	data=DataCMSA(cma);
	a_type	A=SeqSetA(data);
	e_type	E,E0;
	fm_type	*M=ModelsCMSA(cma);
	Int4	J,m,x,*s,pos[5],start,end,repeats;
	Int4	number=gss->NumSeq(),full_len,maxfull;
	char	id0[202],id[202];
	double	scale,page_width=7000.0;
	char	mode,DomShape;
	short	shape=10,domshape;
	dom_typ	*dom=DomainsCMSA(cma);
	UInt4 white = 16777215; // = [255,255,255]
	UInt4 orange = 3381759; // = [255,153,51]

	// if(data != DataCMSA(cma)) print_error("PutGMSA( ) input error");
	if(nBlksCMSA(cma) != 1) print_error("PutSchematicCMSA( ) input error");
	h_type HG=Histogram("numbers of repeats",0,50,1); // added later (11-11-09 by afn)

	BooLean	*arein=0;
	if(keydata) arein=AreInSeqSet(data, keydata);

    // 0. Set up shape & color scheme.
	// if(ShapeType > 0 && (ShapeType <= 10 || ShapeType==56)) shape=ShapeType;
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
	   print_error("PutSchematicCMSA( ) 11\" x 17\" landscape not yet implemented");
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
	  default: print_error("PutSchematicCMSA( ) PageSetUp input error");
	}
#if 0
	fprintf(fp,"\\pard\\plain \\nowidctlpar\\tx%d\\adjustright\n",
		(Int4)page_width+500);
#endif
	fprintf(fp,"\\pard \\nowidctlpar\\tx%d\\adjustright \\sl200\\slmult1 ",
		(Int4)page_width+500); // compressed...

    // 2. Find scaling factor
    	for(maxfull=0,J=1; J <= number; J++) {
	   if(arein && !arein[J]) continue;
	   full_len = gss->TrueFullLength(J);
	   maxfull = MAXIMUM(Int4,maxfull,full_len);
	} scale = page_width/(double)maxfull;
    // 3. Loop over all sequences in aligment.
	Int4	lastJ=0;
    	for(id0[0]=0,repeats=0,J=1; J <= number; J++) {
	   if(arein && !arein[J]){ continue; }
	   E = gss->TrueSeq(J);
	   if(lastJ==0) E0 = E;
	   full_len = gss->TrueFullLength(J);
	   StrSeqID(id,200,E);
	   if(strncmp(id,id0,200) != 0){ // New sequence; output line.
		if(lastJ > 0){	// then need to output end of last seq.
	           if(dom){
	             Int4 Sq = cma->SubToFull[lastJ];
	             for(Int4 d=1; d <= dom->NumDomains(Sq); d++){
			DomShape = dom->Shape(Sq,d);
			if(DomShape == 'e'){ // i.e., biased region.
			   fillcolor = get_explosion_color_cmsa_rtf(Sq,d,cma);
		   	   linewidth = thinline; // lineshade=orange;
			} else {
			 GetColorAndShadeCMSA2RTF(dom->Color(Sq,d),FillColor,
				&fillcolor, &lineshade,dom->Evalue(Sq,d));
			 linewidth=mediumline;
			}
			domshape = GetShapeCMSA2RTF(DomShape);
			if(DomShape == 'r') mode = 'n'; else mode = 'N';
		        put_rtf_shape_cmsa(fp,group,dom->SqStart(Sq,d),
				dom->SqEnd(Sq,d), domshape,
				scale,fillcolor,lineshade,linewidth,mode);
	             	   //             dom->DomainName(Sq,d);
	             }
	           }
		   if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		   else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		   // fprintf(fp,"%s",id0); // old
		   PutShortSeqID(fp,E0);
		   fprintf(fp," (%d)\\par\\par}\n",repeats);
	     	   IncdHist(repeats, HG);
		   E0 = E; repeats=0;
		} lastJ=J;
		// Find length of full sequence
		put_rtf_line_cmsa(fp,group,full_len,scale);
		strcpy(id0,id);
	   }
	   // output block
	   PosSiteCMSA(1,J,pos,cma);
	   start = gss->TrueSite(J,pos[1]);  // offset is added here...
	   end = pos[1]+LengthCMSA(1,cma)-1;
	   end = gss->TrueSite(J,end);
	   // shade according to model probabilities.
	   // double p = GetProbCMSA(1,J,cma);
	   double p = GetGappedProbCMSA(1,J,cma);
	   GetColorAndShadeCMSA2RTF(FillColor,'B',&fillcolor,&lineshade,p);
	   put_rtf_shape_cmsa(fp,group,start,end,shape,scale,
					fillcolor,lineshade,thickline);
	   repeats++;
	}
	if(lastJ > 0){	// then need to output end of last seq.
	    if(dom){
	         Int4 Sq = cma->SubToFull[lastJ];
	         for(Int4 d=1; d <= dom->NumDomains(Sq); d++){
			DomShape = dom->Shape(Sq,d);
			if(DomShape == 'e'){ // i.e., biased region.
			  fillcolor = get_explosion_color_cmsa_rtf(Sq,d,cma);
		   	  linewidth = thinline; // lineshade=orange;
			} else {
			  GetColorAndShadeCMSA2RTF(dom->Color(Sq,d),FillColor,
				&fillcolor, &lineshade,dom->Evalue(Sq,d));
			  linewidth=mediumline;
			}
			domshape = GetShapeCMSA2RTF(DomShape);
			if(DomShape == 'r') mode = 'n'; else mode = 'N';
		        put_rtf_shape_cmsa(fp,group,dom->SqStart(Sq,d),
				dom->SqEnd(Sq,d), domshape,
				scale,fillcolor,lineshade,linewidth,mode);
	         }
	    }
	    if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	    else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	    // fprintf(fp,"%s",id0); // old
            PutShortSeqID(fp,E0);
	    // fprintf(fp," (%d)\\par\\par}\n",repeats); E0 = E;

	    if(dom){ // output names of domains...
	     Int4 t;
             fprintf(fp," (%d)\\par\\par\\par\\par\\par\\par}\n",repeats);
	     IncdHist(repeats, HG);
	     Int4 maxfull2=0; double scale2,page_width2;
	     for(t=1; t <= dom->NumTypes( ); t++){
	        maxfull2 = MAXIMUM(Int4,maxfull2,dom->AveLen(t)+40);
	     } 
	     page_width2 = page_width *((double)maxfull2/(double) maxfull);
	     scale2 = page_width2/(double)maxfull2;
#if 0
	     fprintf(fp,"\\pard\\plain \\nowidctlpar\\tx%d\\adjustright\n",
			(Int4)page_width2+500);
#endif
	     fprintf(fp,"\\pard\\plain \\nowidctlpar\\tx%d\\adjustright \\sl190\\slmult1 ",
		(Int4)page_width2+500); // compressed...
	     // put default domain or repeat...
		start = 1 + ((maxfull2 - LengthCMSA(1,cma))/2);
		end = start + LengthCMSA(1,cma) -1;
	        put_rtf_line_cmsa(fp,group,maxfull2,scale2);
	        GetColorAndShadeCMSA2RTF(FillColor,'B',&fillcolor,&lineshade,10.0);
	        put_rtf_shape_cmsa(fp,group,start,end,shape,scale2,
			fillcolor,lineshade,thickline);
	        if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
                fprintf(fp,"default (%d)\\par\\par}\n",TotSitesFModel(M[1]));
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
	        put_rtf_shape_cmsa(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode);
	        if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        fprintf(fp,"%s",dom->NameThisType(t)); 
                fprintf(fp," (%d)\\par\\par}\n",dom->NumThisType(t));
	     }
#if 1 // NEW
	     Int4 ave=0,total=0;
	     for(t=1; t <= dom->NumTypes( ); t++){
		DomShape = dom->Shape(t);
		domshape = GetShapeCMSA2RTF(DomShape);
		if(domshape < 60) continue;
		else { total++; ave+=dom->AveLen(t); }
	     }
	     if(total > 0){
		ave = (ave/total) + 1;
		domshape = GetShapeCMSA2RTF('e');
	        start=1+(maxfull2-ave)/2; end = start + ave-1;
	   	linewidth = thinline; mode = 'N';
		unsigned char *pc; lineshade=0;
	        enum primary_color { wht=0, blue, green, red };
		for(Int4 x = wht; x <= red;  x++){
	          put_rtf_line_cmsa(fp,group,maxfull2,scale2);
		  pc = (unsigned char*) &fillcolor; fillcolor=0;
		  pc[1]=pc[2]=pc[3]=100;
		  if(x==wht) fillcolor = white;
		  else { pc[x] += 150; }
	          put_rtf_shape_cmsa(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode);
	          if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	          else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
		  switch(x){
		   case wht: fprintf(fp,"hydrophobic_rich");  break;
	           case blue: fprintf(fp,"basic_rich");  break;
	           case green: fprintf(fp,"polar_rich");  break;
	           case red: fprintf(fp,"acidic_rich");  break;
		   default: assert("This should not happen!\n");
		  } fprintf(fp," \\par\\par}\n");
		}
	     }
#endif
#if 0 // old
	     for(t=1; t <= dom->NumTypes( ); t++){
		DomShape = dom->Shape(t);
		domshape = GetShapeCMSA2RTF(DomShape);
		if(domshape < 60) continue;
	        start=1+(maxfull2-dom->AveLen(t))/2;
	        end = start + dom->AveLen(t) -1;
	        put_rtf_line_cmsa(fp,group,maxfull2,scale2);
	        GetColorAndShadeCMSA2RTF(dom->Color(t),FillColor,
			&fillcolor,&lineshade,10.0);
		linewidth=mediumline;
		if(DomShape == 'r') mode = 'n'; else  mode = 'N';
	        put_rtf_shape_cmsa(fp,group,start,end,domshape,
			scale2,fillcolor,lineshade,linewidth,mode);
	        if(!group) fprintf(fp,"}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        else fprintf(fp,"}}}{\\f2\\fs20\\lang1024\\cgrid0\\b \\tab ");
	        fprintf(fp,"%s",dom->NameThisType(t)); 
                fprintf(fp," (%d)\\par\\par}\n",dom->NumThisType(t));
	     }
#endif
	    } else { fprintf(fp," (%d)\\par\\par}\n",repeats); IncdHist(repeats, HG);} E0 = E;
	} fprintf(fp,"}\n");
	if(arein) free(arein);
	PutHist(stderr,60,HG); NilHist(HG);
}

//************************* output rich text format alignment *********************

BooLean	*ConservedCMA(Int4 *observed, double cutoff, double *freq, a_type A)
/* return the conserved residues in column "observed". */
{
	Int4	i,r,r2,nres,n;
	double	total,p,q,f,low_cut;
	BooLean	*conserved,flag,debug=0;

	assert(cutoff <= 0.0);
	low_cut=cutoff + 0.7;
	if(observed != NULL){
	    nres=0;
	    NEW(conserved,nAlpha(A)+2,BooLean);
	    for(total=0.0, r=1; r<=nAlpha(A); r++){
		total += (double)observed[r];
	    }
	    if(debug){
		fprintf(stderr,"cutoff = %.2f; low_cut = %.2f\n",cutoff,low_cut);
		for(r = 1; r <= nAlpha(A); r++){
		   if(observed[r] > 0){
			q = freq[r];
			fprintf(stderr,"%c: %d/%d (p=%.1f)\n", 
				AlphaChar(r,A), observed[r], (Int4)total,
				Log10CBP(observed[r], total, q));
   		   }
		} fprintf(stderr,"\n");
	    }
	    /** 1. find clearly conserved residues **/
	    for(f=0.0, n=0, r=1; r <= nAlpha(A); r++){
		if(observed[r] > 2){
		   q = freq[r];
		   p = Log10CBP(observed[r], total, q);
		   if(p <= cutoff) {
			if(debug) fprintf(stderr," %c: (%.1f)[%d/%d]", 
				AlphaChar(r,A), 
				p,observed[r], (Int4)total);
			n += observed[r]; f += q; 
			conserved[r] = TRUE; nres++;
		    }
		}
	    }
            /** 2. find moderately conserved related residue pairs **/
            if(nres==0){
              for(n=0, f=0.0, r=1; r < nAlpha(A); r++){
                if(!conserved[r] && observed[r] > 2){
                  q = freq[r];
		  p = Log10CBP(observed[r], total, q);
                  if(p <= low_cut){
                    for(r2 = 1; r2 <= nAlpha(A); r2++){
                        if(!conserved[r2] && r2 != r && observed[r2] > 2 
                                        && valAlphaR(r,r2,A) > 0){
                           q = freq[r2];
		  	   p = Log10CBP(observed[r2], total, q);
                           if(p <= low_cut) {
				if(debug) fprintf(stderr," %c %c (%.1f)", 
					AlphaChar(r,A), AlphaChar(r2,A), p);
                                nres+=2;
				n += observed[r] + observed[r2];
				f += freq[r] + freq[r2];
                                conserved[r]=conserved[r2]=TRUE;
                           }
                        }
                    }
                  }
                }
              }
            } 
	    if(nres==0){ free(conserved); return NULL; }
	    /** 3. find weakly conserved residues related to clear ones **/
            do {
	      flag = FALSE;
	      for(r=1; r <= nAlpha(A); r++){
		if(conserved[r]){
	          for(r2 = 1; r2 <= nAlpha(A); r2++){
		   if(!conserved[r2] && observed[r2] >= 2){
			if(valAlphaR(r,r2,A) >= 0){
		   	   q = freq[r2];
		           p = Log10CBP(observed[r2], (total-n), q/(1.0 - f));
		   	   if(p <= low_cut){ 
				if(debug) fprintf(stderr," %c (%.1f)", 
						AlphaChar(r2,A), p);
				flag = TRUE; nres++;
				conserved[r2] = TRUE;
			   }
			}
		   }
		  }
		}
	      }
	    } while(flag);
	    /** 3. reset to off if total conserved residues too low **/
	    for(q=0.0,r2=0, r = 1; r <= nAlpha(A); r++){
	      if(conserved[r]){ r2 += observed[r]; q += freq[r]; }
	    }
	    p = Log10CBP(r2, total, q);
	    if(debug) fprintf(stderr,"  (%.0f%c conserved)(p=%.1f)\n",
		100.0*(double)r2/total,'%',p);
	    // if(p > -3.0) { free(conserved); return NULL; }
	    if(p > -2.0) { free(conserved); return NULL; }
	    p = (double)r2/total;
	    if(p < 0.20) { free(conserved); return NULL; }
	    return conserved;
	} else return NULL;
}

void    CMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, char PageSetUp, cma_typ cma,
	BooLean LabeledOnly)
/********************  from driver program ******************
        FILE fptr=open_file(argv[1],".rtf","w");
        if(key_name) key_seq = SeqSet(key_name,AlphaSMA(MA));
        CMA2RTF(fp, log10(cbp_cut), infoLO,infoHI,key_seq,MA);
        fclose(fp);
        if(key_name) NilSeqSet(key_seq); NilSMA(MA);
/************************************************************/
{
	FILE	*fp=tmpfile();
	double	***freq;
	Int4	b;
	
	PutAlnCMSA(fp,cma); rewind(fp);
        sma_typ MA=ReadSMA(fp); fclose(fp);

	NEWPP(freq,nBlksCMSA(cma)+3,double);
	for(b=1; b <= nBlksCMSA(cma); b++){
	   NEWP(freq[b],LengthCMSA(b,cma)+3,double);
	   // for(Int4 s=1; s <=LengthCMSA(b,cma); s++){ freq[b][s]=FreqSMA(MA);	}
	}
	CMA2RTF(fptr,cbp_cut,purge_cut,infoLO,infoHI,key_seq,PageSetUp,freq,MA,cma,LabeledOnly);
	for(b=1; b <= nBlksCMSA(cma); b++){ free(freq[b]); } free(freq);
	NilSMA(MA);
}

void    CMA2RTF(FILE *fptr, double cbp_cut, double purge_cut, double infoLO, 
	double infoHI, ss_type key_seq, char PageSetUp, double ***freq, 
	sma_typ MA,cma_typ cma,BooLean LowerOnly)
/** PutSites to create a figure for publication **/
/** 1 black;2 blue;3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white **/
/** 9 lt blue;10 lt cyan; 11 lt green;12 ; 13; 14 lt magenta; 15 dk grey; 16 lt grey **/
{
	Int4	t,i,j,n,s,e,end,N,start,length,*observed,gap,ntyp;
	Int4	h,p,aliphatic,polar,consH,consP,consT,consMax,obs;
	float	percentH,percentP,percentT,percentDiff,percentA,percentMax;
	double	info,total,d,f,q;
	a_type	A=AlphabetCMSA(cma);
	char	r,c,str[50];
	char	*cons_state,state,laststate;  
	BooLean	**conserved,flag,*use;
// char PageSetUp='L'; // versus p, P, l.

    assert(cbp_cut < 0.0 && infoLO >= 0.1 && infoLO < infoHI);
    ntyp = ntypSMA(MA);
    fprintf(fptr,"{\\rtf1\\ansi \\deff4\\deflang1033");
    fprintf(fptr,"{\\fonttbl{\\f0\\froman\\fcharset0\\fprq2 Tms Rmn;}\n");
    fprintf(fptr,"{\\f1\\fswiss\\fcharset0\\fprq2 Arial;}\n");
    fprintf(fptr,"{\\f2\\fswiss\\fcharset0\\fprq2 Arial Narrow;}\n");
    fprintf(fptr,"{\\f3\\fmodern\\fcharset0\\fprq1 Courier New;}\n");
    fprintf(fptr,"{\\f4\\froman\\fcharset0\\fprq2 Times New Roman;}\n");
    fprintf(fptr,"{\\f5\\fswiss\\fcharset0\\fprq2 System;}\n");
    fprintf(fptr,"{\\f6\\fmodern\\fcharset0\\fprq1 Courier New;}}\n");
    fprintf(fptr,"{\\colortbl;\\red0\\green0\\blue0;\\red0\\green0\\blue255;");
    fprintf(fptr,"\\red0\\green255\\blue255;\\red0\\green255\\blue0;");
    fprintf(fptr,"\\red255\\green0\\blue255;\\red255\\green0\\blue0;");
    fprintf(fptr,"\\red255\\green255\\blue0;\\red255\\green255\\blue255;");
    fprintf(fptr,"\\red0\\green0\\blue128;\\red0\\green128\\blue128;");
    fprintf(fptr,"\\red0\\green128\\blue0;\\red128\\green0\\blue128;");
    fprintf(fptr,"\\red128\\green0\\blue0;\\red128\\green128\\blue0;");
    /** 15 dk grey; 16 lt grey (original) **/
    fprintf(fptr,"\\red128\\green128\\blue128;\\red192\\green192\\blue192;}");
    fprintf(fptr,"{\\stylesheet{\\widctlpar \\f4\\fs18 \\snext0 Normal;}");
    fprintf(fptr,"{\\*\\cs10 \\additive Default Paragraph Font;}}\n");

    switch(PageSetUp){
	  case 'l': // 8.5" x 11" landscape.
	   fprintf(fptr,"\\paperw15840\\paperh12240\\margl720\\margr720");
    	   fprintf(fptr,"\\margt1080\\margb1080 ");
    	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
    	   // page_width=12000.0;
	   break;
	  case 'L': // 11" x 17" landscape
	   fprintf(fptr,"\\paperw24480\\paperh15840\\margl720\\margr720");
	   fprintf(fptr,"\\margt1080\\margb1080 ");
	   fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    	   fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
	   // print_error("CMA2RTF( ) 11\" x 17\" landscape not yet implemented");
	   break;
	  case 'p': // 8.5" x 11" portrait.
	   // page_width=7000.0;
	   break;
	  case 'P': // 11" x 17" portrait.
	    fprintf(fptr,"\\paperw15840\\paperh24480\\margl720\\margr720");
            fprintf(fptr,"\\margt1080\\margb1080 ");
            fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
            fprintf(fptr,"\\linex0\\endnhere\\sectdefaultcl \n");
            // page_width=12000.0;
	   break;
	  default: print_error("CMA2RTF( ) PageSetUp input error");
    }
    fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs18 \n");
    // fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs18 \\sl180\\slmult1 ");
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
    use=purge_sma_cma(purge_cut,MA,cma);
    NEW(observed,nAlpha(A)+2,Int4);
    h_type HG=Histogram("information content in bits",0,100,0.05);
    for(t=1; t <= ntyp; t++){
        h_type HG2=Histogram("Kyte-Doolittle hydrophobicity",0,100,0.2);
	length = lengthSMA(t,MA);
	NEWP(conserved, length+2, BooLean);
	NEW(cons_state, length+2, char);
	for(s=1; s <= length; s++){
	    for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	    for(n=1; n<=nseqSMA(MA); n++){
		if(use[n]){ r=residueSMA(t,n,s,MA); observed[r]++; }
	    }
	    // use column freq from alternative cma here for subfamily models...
	    if(!freq[t][s]) conserved[s] = ConservedCMA(observed,cbp_cut,FreqSMA(MA),A);
	    else conserved[s] = ConservedCMA(observed,cbp_cut,freq[t][s],A);
	    obs=consMax=consT=consH=consP=aliphatic=polar=0;
	    for(r=1; r<= nAlpha(A); r++){
		n = observed[r];
		sprintf(str,"%c",AlphaChar(r,A));
		if(strstr("AILVFYWMCP",str) != NULL) h=n; else h=0;
		if(strstr("DEKRHSTQN",str) != NULL) p=n; else p = 0;
		if(conserved[s] != NULL && conserved[s][r]){
		    consMax = MAXIMUM(Int4,consMax,n);
		    consT+=n; consH+=h; consP+=p;
		}
		aliphatic+=h; polar+=p; obs+=n;
	    }
/******* compute relative entropy of postions *******/
            for(total=0.0,r=1; r <= nAlpha(A); r++){ total += (double) observed[r]; }
	    double kd_score = 0.0;
            for(info=0.0,r=1; r <= nAlpha(A); r++){
		if(observed[r] > 0){
			kd_score += observed[r]*blsm62kyte_doolittle[r];
                        f = (double) observed[r]/total;
                        if(freq[t][s]) q = freq[t][s][r];
                        else q = freqSMA(r,MA);
                        if(q > 0.0) info += f*log(f/q)/log(2.0);
                } 
            }
// fprintf(stderr,"info[%d][%d] = %.2f\n",t,s,info);
	    kd_score /= total;
            IncdHist(kd_score,HG2); IncdHist(info,HG);
/******* relative entropy of postions *******/
	    percentT = ((float)consT/(float)obs);
	    // percentH = ((float)consH/(float)obs);
	    percentH = ((float)consH/(float)consT);
	    percentP = ((float)consP/(float)obs);
	    percentMax = ((float)consMax/(float)obs);
	    percentA = ((float)aliphatic/(float)obs);
	    percentDiff = (double) fabs(percentH-percentP);
	    if(percentA >= 0.80 && percentH >= 0.70){
	           for(r=1; r<= nAlpha(A); r++){
			sprintf(str,"%c",AlphaChar(r,A));
			if(strstr(" AILVFYWMCP",str) != NULL &&
			    			!conserved[s][r]){
				conserved[s][r]=TRUE;
			}
		   }
	    }
	    if(info >= infoHI){ 
		if(kd_score > 5.0) cons_state[s] = 'H'; else cons_state[s] = 'I'; 
	    } else if(info >= infoLO){ 
		if(kd_score > 5.0) cons_state[s] = 'h'; else cons_state[s] = 'i'; 
	    } else if(percentH >= 0.98 && percentA >= 0.70){ 
		cons_state[s]='h'; 
	    } else if(percentA >= 0.90 && percentH >= 0.70 ){
		cons_state[s] = 'a';
	    } else if(percentT < 0.50){
		if(percentDiff < 0.33) cons_state[s] = 'w';
		else cons_state[s] = 'W';
	    } else if(percentT < 0.90){
		if(percentDiff < 0.66) cons_state[s] = 'm';
		else cons_state[s] = 'M';
	    } else {	// i.e., if(percentT >= 0.90)...
		if(kd_score > 5.0) cons_state[s] = 'a';
		else cons_state[s] = 's'; 
	    }
	}
    	PutHist(stderr,60,HG2);  NilHist(HG2);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
	BooLean	*show_position=0;
	if(key_seq){
	  NEW(show_position,glengthSMA(t,MA) + 2, BooLean);
	  char *gnull = gnullSMA(t,MA);
	  for(i=j=0; j < lengthSMA(t,MA); i++){
	    if(gnull[i] != '_'){ j++; show_position[i]=TRUE; }
	    else for(n=1; n<= nseqSMA(MA); n++){
	      char *gseq = gseqSMA(t,n,MA);
    	      if(IsInSeqSet(key_seq,seqidSMA(n,MA))){
		c = gseq[i];
		if(c != '.') show_position[i] = TRUE;
	      }
	    }
	  }
	}
	for(N=0, n=1; n<= nseqSMA(MA); n++){
    		if(key_seq && !IsInSeqSet(key_seq,seqidSMA(n,MA))) continue;
#if 1   // new option to show lowercase kingdom sequences only.
		if(LowerOnly){
                  e_type trueE=TrueSeqCMSA(n,cma);
                  // fprintf(stderr,"seq n kingdom=%c\n",kingdomSeq(trueE));
                  // PutSeqInfo(stderr,trueE);
                  if(!islower(kingdomSeq(trueE))) continue;
		  else fprintf(stderr,"seq n kingdom=%c\n",kingdomSeq(trueE));
		}
#endif
	   	if(t==1){
		   fprintf(fptr,"{\\f2\\fs16 ");
		   fprintf(fptr,"%s",seqidSMA(n,MA));
		   fprintf(fptr,"\\tab %d\\tab }",startSMA(1,n,MA));
	   	}
		if(t < ntypSMA(MA)) {
			gap = startSMA(t+1,n,MA) - endSMA(t,n,MA) - 1;
		}
		if(IsProbSMA(MA)){
		  if(probSMA(n,t,MA) < 0.0) strcpy(str,"\\b\\i\\strike");
		  else if(probSMA(n,t,MA) < 1.3) strcpy(str,"\\b\\i");
		  else strcpy(str,"\\b");
		} else strcpy(str,"\\b");
// new 
		char *gnull = gnullSMA(t,MA);
	        char *gseq = gseqSMA(t,n,MA);
		unsigned char *seq = seqSMA(t,n,MA);
		for(laststate=0,i=j=0; j < lengthSMA(t,MA); i++){
		   if(gnull[i] == '_'){
			if(!show_position || show_position[i]){
			  state='g'; c = gseq[i];
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    fprintf(fptr,"{%s\\f3\\fs14\\cf16 %c",str,c);  // gap color
			  }
			  laststate = state;
			}
		   } else {
			j++;
			r = seq[j]; 
			if(r==0 && gseq[i] == '-') c = '-';
			else c = AlphaChar(r,A);
			if(conserved[j]==NULL || !conserved[j][r]) state='u';
			else state = cons_state[j];
			if(state == laststate) fprintf(fptr,"%c",c);
			else {
			   if(laststate != 0) fprintf(fptr,"}\n");
			   switch(state){
			     case 'H':
				fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str, c);
				 break;  /** red on yellow **/
			     case 'h':
				fprintf(fptr,"{%s\\f3\\cf2\\highlight7 %c",str, c);
				 break;  /** blue on yellow **/
			     case 'a':
				fprintf(fptr,"{%s\\f3\\cf1\\highlight7 %c",str, c);
				 break;  /** black on yellow **/
			     case 'I': case 'S':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight6 %c",str, c);
				 break;  /** white on red **/
			     case 'i': case 's':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight5 %c",str, c);
				 break;  /** white on magenta **/
			     case 'M':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight1 %c",str, c);
				 break;  /** white on black **/
			     case 'm':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight15 %c",str, c);
				 break;  /** white on dark gray **/
			     case 'W':
				  fprintf(fptr,"{%s\\f3\\cf1 %c",str,c);
				break;  /** black on white **/
			     case 'w':
				  fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
				break;  /** dark gray on white **/
			     case 'u':
			     default : 
				  fprintf(fptr,"{%s\\f3\\cf16 %c",str,c);
				 break;  /** light grey on white **/
			   } 
			}
			laststate = state;
		   }
		}
		fprintf(fptr,"}\n{\\f2\\fs16 ");
		if(t==ntyp) fprintf(fptr,"\\tab %d", endSMA(t,n,MA));
		else if(gap <=0) fprintf(fptr,"\\tab   ");
		else fprintf(fptr,"\\tab (%d)", gap);
#if 1
		if(IsProbSMA(MA)) {
			fprintf(fptr,"\\tab [%.1f]",probSMA(n,t,MA));
		}
#endif
		fprintf(fptr,"\\tab }{\\f3 \n\\par }");
	}  /** end sequences **/
	if(show_position) free(show_position);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3\n\\par\n\\par }\n");
	for(s=1; s <= length; s++)
		if(conserved[s] != NULL) free(conserved[s]);
	free(conserved);
	free(cons_state);
    }
    PutHist(stderr,60,HG);  NilHist(HG);
    free(observed); free(use);
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1  {\\f3\n\\par}}\n");
}

BooLean	*purge_sma_cma(double cutoff, sma_typ MA,cma_typ cma)
/**************************************************************************
  Purge from a list of sequences.
/**************************************************************************/
{
	Int4	totlen,i,k,s,t,entry,N,score,max;
	Int4	item,item1,*hits,num_seqs=0;
	b_type	*edges;
	dh_type	H;
	keytyp	key;
	BooLean	*itemList;
	a_type	A=AlphaSMA(MA);
	unsigned char	**seq1,**seq2;
	double	prcntID;

	assert(cutoff > 0.0);
	N=nseqSMA(MA);
	if(cutoff > 100.){	// return all TRUE if cutoff > 100%
	   NEW(itemList,N+3,BooLean);
	   for(i=1; i <= N; i++) itemList[i]=TRUE;
	   return itemList;
	} 
	for(totlen=0,t=1; t <= ntypSMA(MA); t++) totlen+=lengthSMA(t,MA);
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); NEW(hits,N+2,Int4);
	for(entry =1 ;entry <= N;entry++) edges[entry]=Block(N+1);
	entry = 1;
	for(item =1 ;item <= N;item++) {
	    seq1 = SeqSMA(item,MA);
	    for(item1 =item+1 ;item1 <= N;item1++) {
		seq2 =  SeqSMA(item1,MA);
		for(score=0,t=1; t <= ntypSMA(MA); t++){
                     for(s=1; s <= lengthSMA(t,MA); s++){
			if(seq1[t][s] == seq2[t][s]) score++;
                     }
                } 
		prcntID = 100.*(double)score/(double)totlen;
		if(prcntID >= cutoff){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    } entry++;
	}
	fprintf(stderr,"\n%d items compared; cutoff %.2f\n", entry-1,cutoff); 
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RmHomologs( )");
	    else if(minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=1;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
	num_seqs = ItemsInHeap(H);
	NEW(itemList,N+3,BooLean);
	for(k=0,entry=1; entry <= N; entry++){
	    if(memHeap(entry,H)){ itemList[entry]=TRUE; k++; }
	    else itemList[entry] = FALSE;
	}
	fprintf(stderr,"%d sequences left after purging\n",k);
	for(entry=1; entry <= N; entry++) { NilBlock(edges[entry]); }
	free(hits); free(edges); Nildheap(H);
	return itemList;
}


