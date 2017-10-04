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

#include "vsi_typ.h"

Int4    PrintDotsItems(FILE *fp, item_list_type *list, pml_typ *pml)
{
	item_list_type *node,*next;
	dots_list_type *dlist=0,*dnode;
	char	*str;
	Int4	N=0,i,j,k,num;

	// 1. get list of things for dot clouds...
	for(num=0,node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  switch(node->type){
#if 0
	   case 'T':	// Trace type
	    if(itm.trace.backbone == '@'){
	      num++;
	      char str2[50];
	      if(itm.trace.chain) {
	         sprintf(str2,"%d-%d%c",itm.trace.start,itm.trace.end,itm.trace.chain);
	      } else { sprintf(str2,"%d-%d",itm.trace.start,itm.trace.end); }
	      str=AllocString(str2);
	      if(dlist==0){ 
		dnode=MakeDotsList(str,itm.trace.color); dlist=dnode;
	      } else {
		dnode->next=MakeDotsList(str,itm.trace.color); dnode=dnode->next;
	      }
	    }
#endif
	   case 'R':
	    if(itm.residue.clouds == 'T'){
	      num++;
#if 0
	      str=pml->DefineResCloud(fp,itm.residue.res,itm.residue.site,
				itm.residue.atom,itm.residue.hydrogen,
				itm.residue.chain);
#else
	      N++; pml->DefineResCloud(fp,itm.residue.res,itm.residue.site,itm.residue.chain);
#endif
#if 0
	      if(dlist==0){ 
		dnode=MakeDotsList(str,itm.residue.color); dlist=dnode;
	      } else {
		dnode->next=MakeDotsList(str,itm.residue.color); dnode=dnode->next;
	      }
#endif
	    } break;
#if 0
	   case 'M':
	    if(itm.molecule.clouds == 'T'){
	      num++;
	      str=pml->DefineMolCloud(fp,itm.molecule.mol, itm.molecule.id,
			itm.molecule.atom,itm.molecule.hydrogen,itm.molecule.chain);
	      if(dlist==0){ 
		dnode=MakeDotsList(str,itm.molecule.color); dlist=dnode;
	      } else {
		dnode->next=MakeDotsList(str,itm.molecule.color);
	     	dnode=dnode->next;
	      }
	    } break;
#endif
	   default: 
	    break;
	  }
	}
#if 0
	// 2. get list of things for dot clouds...
	N=0,i=0;
   if(num > 0){
	for(dnode=dlist; dnode != 0; dnode=dnode->next){
	   if(N % 20 == 0) {
		i++; fprintf(fp,"\ndefine Dots%d %s",i,dnode->str); 
	   } else {
		fprintf(fp,",%s",dnode->str); 
	   } N++;
	} fprintf(fp,"\ndefine DotsAll ");
	for(j=1; j < i; j++) fprintf(fp,"Dots%d,",j); fprintf(fp,"Dots%d\n",i);
	fprintf(fp,"select DotsAll\ndots 1000\n");
   }
#endif
	return N;
}

Int4	DefineItemsPyMol(FILE *fp, item_list_type *list, char key_chain,pml_typ *pml)
{
	item_list_type *node;

	char su=key_chain,chn,c;
        if(islower(su)){ su=toupper(su); }
	double factor=1.5;

        if(su == ' ') chn=0; else chn=su;
	for(node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  Int4 *v;
	  switch(node->type){
	   case 'R':
#if 0	// DEBUG	// last node occurs multiple times; need to debug...
if(itm.residue.color == 'Y'){
	fprintf(stderr,"%c%d%c\n",itm.residue.res, itm.residue.site,itm.residue.color);
}
#endif
	     if(itm.residue.atom == 0){
	        pml->DefineResidueItem(fp,itm.residue.res, itm.residue.site,
				itm.residue.chain, itm.residue.color);
	     }
	    break;
#if 0
	   case 'M': 
	     if(itm.molecule.clouds == 'S') spacefill=9999;
	     else spacefill=(Int4) ceil((double)itm.molecule.width * factor);
	     pml->PutMoleculeItem(fp,itm.molecule.mol, itm.molecule.id,
			itm.molecule.atom, itm.molecule.hydrogen,
			itm.molecule.chain, itm.molecule.color,
			itm.molecule.width, spacefill);
	    break;
#endif
	   case 'T': 
		if(itm.trace.chain==0) c=chn; else c=itm.trace.chain;
		if(c==' ') c=0;
		pml->PrintTrace(fp,itm.trace.start,itm.trace.end,c,
				itm.trace.color,itm.trace.backbone,itm.trace.width);
	    break;
#if 0
	   case 'V': 
		v=itm.view; 
		pml->PrintView(fp,v[0],v[1],v[2],v[3],v[4],v[5]);
	    break;
#endif
	   default: pml->PrintCommand(fp,itm.cmd,node->type); break;
	   // default: break;	// do nothing...
	  }
	}
	return (0);
}


Int4	PrintItemsPyMol(FILE *fp, BooLean look_in_datadir, item_list_type *list,
	char *filename, char key_chain,BooLean cis_trans_pro,
	Int4 wire_width,Int4 spacefill)
{
	item_list_type *node;
	if(filename == 0){
	   fprintf(stderr,"\"File#=pdbfilename:Chain\" not designated!\n"); 
	   print_error("**************** Fatal error! ******************");
	}
	pml_typ *pml = new pml_typ(filename,cis_trans_pro);

	char su=key_chain,chn,c;
        if(islower(su)){ su=toupper(su); }
	double factor=1.5;

        if(su == ' ') chn=0; else chn=su;
        // if(su != '-') pml->PrintHEADER(fp,chn,wire_width); else 
	pml->PrintHEADER(fp,look_in_datadir,key_chain,wire_width);
	DefineItemsPyMol(fp, list, key_chain,pml);	// show sidechains...
	for(node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  Int4 *v,i,j;
	  switch(node->type){
	   case 'R':
	     pml->PutResidueItem(fp,itm.residue.res, itm.residue.site,
				itm.residue.atom, itm.residue.hydrogen,
				itm.residue.chain, itm.residue.color,
				itm.residue.width, spacefill);
	    break;
	   case 'M': 
#if 1	// fixes rasmol formating '[SO4]' to 'SO4'.
	    if(itm.molecule.mol[0]=='['){
	        char *mol=itm.molecule.mol;
		for(i=0,j=1; mol[j] != ']'; i++,j++){ mol[i]=mol[j]; } mol[i]=0;
	    }
#endif
	     pml->PutMoleculeItem(fp,itm.molecule.mol, itm.molecule.id,
			itm.molecule.atom, itm.molecule.hydrogen,
			itm.molecule.chain, itm.molecule.color,
			itm.molecule.width, spacefill);
	    break;
#if 0
	   case 'T': 
	     break;
		if(itm.trace.chain==0) c=chn; else c=itm.trace.chain;
		if(c==' ') c=0;
		pml->PrintTrace(fp,itm.trace.start,itm.trace.end,c,
				itm.trace.color,itm.trace.backbone,itm.trace.width);
	    break;
	   case 'V': 
	     break;
		v=itm.view; 
		pml->PrintView(fp,v[0],v[1],v[2],v[3],v[4],v[5]);
	    break;
	   default: pml->PrintCommand(fp,itm.cmd,node->type); break;
#else
	   default: break;
#endif
	  }
	}
	if(0) PrintDotsItems(fp,list,pml);	// omit for now...
	pml->PrintTAIL(fp);
	delete pml; return (0);
}

