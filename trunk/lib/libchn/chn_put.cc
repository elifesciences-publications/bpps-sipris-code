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

#include "chn_typ.h"

void    chn_typ::PutHierarchicalAln(char *filename)
// Create rich text format output file
{
  FILE	*fp;
  if(OutputAln) return;
  assert(isprint(filename[0]));
  assert(filename[0] != '"');
  fp=open_file(filename,".rtf","w");
  // fprintf(stderr,"**** output_name(2)=%s ****\n",filename);
#if 0	// fix this part...core dumping in FixSetStatus() !!!
  if(RTFsCreated == FALSE){ CreateRTF(); RTFsCreated==TRUE; }
#endif
  PutCHA_HEADER(fp);
  // PutCHA(fp," CyPBGRcmbgr                                                        ");
  PutCHA(fp,ColorString);
  PutCHA_TAIL(fp);
  fclose(fp);
}

void	chn_typ::put_marg_prob(FILE *fptr, Int4 gstart, Int4 start, Int4 gend,
	Int4 color,Int4 gapsize,char *gnull,Int4 Analysis)
// Only called by PutCHA above!!! 
{
	char	str[20];
	char	laststate,state;
	Int4	i,j;
	//**************** output marginal probability... *******************
	fprintf(fptr,"{\\b\\f2\\fs%d\\cf%d misalign prob\\tab \\tab }",
		fontsize,color);
	strcpy(str,"\\b");
	char	marg_prob;
	for(laststate=0,i=gstart,j=start; i < gend && j < End; i++){
		   if(gnull[i] == '_'){
			  state=' ';
			  if(state == laststate) fprintf(fptr," ");
                          else {
                            if(laststate != 0) fprintf(fptr,"}\n");
                            fprintf(fptr,"{%s\\f3\\fs%d\\cf1  ",str,gapsize);
                          } laststate = state;
		   } else {
			j++;
			state='X';
			if(MargProb[Analysis][j] == 0.0) marg_prob = '?';
			else {
			  marg_prob=FractionToCharPBC(1.0-MargProb[Analysis][j]);
    			  if(marg_prob == '0') marg_prob = ' ';
			  // if(marg_prob=='!') marg_prob='9';
			}
			if(state == laststate) fprintf(fptr,"%c",marg_prob);
			else { 
			   if(laststate != 0) fprintf(fptr,"}\n");
	        	   fprintf(fptr,"{%s\\f3\\cf1 %c",str,marg_prob);
			} laststate = state;
		   }
	} fprintf(fptr,"}\n{\\f2\\fs%d ",fontsize);
	fprintf(fptr,"\\tab }{\\f3 \n\\par }");
}

void    chn_typ::PutAln(char *filename,char mode)
{
	FILE	*fp,*ofp;
	Int4	s;
	char	str[200];
	assert(isprint(filename[0]));
        fp=open_file(filename,".csq","w"); PutSeq(fp,qE,AB); fclose(fp);
        // fprintf(stderr,"**** output_name(4)=\"%s\" ****\n",filename);
	if(mode=='W'){
#if 1
	  fp=open_file(filename,".chn","r"); 
	  ofp=0;
	  char str2[10003];
          for(Int4 item=0; ; ){
             if(fgets(str2,10000,fp)==NULL) break;
	     // if(strncmp(str2,"[0_(",4) == 0)
	     Int4 level=0;
	     if(sscanf(str2,"[%d_(",&level) == 1){
		if(level < 0) print_error("chn file syntax error");
		if(ofp) fclose(ofp); ofp=0;
		item++;
		if(item==1){
		   sprintf(str,"%s_ds",filename);
	  	   ofp=open_file(str,".cma","w"); 
		} else if(item==2){
		   sprintf(str,"%s_sf",filename);
	  	   ofp=open_file(str,".cma","w"); 
		} else if(item==Number){
		   sprintf(str,"%s_bg",filename);
	  	   ofp=open_file(str,".cma","w"); 
		}
	     } assert(item > 0);
	     if(ofp){
		fprintf(ofp,"%s",str2);
	     }
          } fclose(fp);
	  if(ofp) fclose(ofp);
#else
	  for(s=1; s <= NumAnalysis; s++){
	     if(s == NumAnalysis){
                sprintf(str,"%s_%d.cma",filename,s); WriteCMSA(str,IN_CMA[Number]);
	     } else {
               sprintf(str,"%s_%d.cma",filename,s); WriteCMSA(str,IN_CMA[s]);
               // sprintf(str,"%s%d.cma",filename,s); WriteCMSA(str,IN_CMA[s]);
	     }
	  }
#endif
	} else if(extra_files){
	  for(s=1; s <= NumAnalysis; s++){
             sprintf(str,".%dsq",s-1); 
	     fp=open_file(filename,str,"w");
	     // PutSeqSetEs(fp,TrueDataCMSA(IN_CMA[s])); 
	     PutSeqSetEs(fp,TrueDataCMSA(MCMA[s-1])); 
	     fclose(fp);
	  }
	}
}

