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

item_list_type *MakeItemNodeTrace(UInt4 start, UInt4 end, 
        	char chain, char backbone, char color,Int4 width)
{
	item_list_type *item_node;
	NEW(item_node, 1, item_list_type);
	item_node->type='T';
	item_node->item.trace.start=start;
	item_node->item.trace.end=end;
	item_node->item.trace.chain=chain;
	item_node->item.trace.backbone=backbone;
	item_node->item.trace.color=color;
	item_node->item.trace.width=width;
	return item_node;
}

item_list_type *MakeItemNodeMol(char    *mol,        // molecule (e.g., 'ATG')
        Int4    id,          // molecule id (e.g., 802.
        char    chain,       // chain
        char    *atom,     // atom type (e.g., 'od1')
        char    *hydrogen, // attached hydrogen (e.g., '1hg')
        char    color,       // color
        char    clouds,          // show dots: color=X --> ions; else --> sidechains
        Int4    width)
{
	item_list_type *item_node;
	NEW(item_node, 1, item_list_type);
	item_node->next=0;
	item_node->type='M';
	item_node->item.molecule.mol=mol;
	item_node->item.molecule.id=id;
	item_node->item.molecule.chain=chain;
	item_node->item.molecule.atom=atom;
	item_node->item.molecule.hydrogen=hydrogen;
	item_node->item.molecule.color=color;
	item_node->item.molecule.clouds=clouds;
	item_node->item.molecule.width=width;
	return item_node;
}

item_list_type *MakeItemNodeRes(char    res,  // residue type (e.g., 'T')
        Int4    site,           // residue site.
        char    chain,          // chain
        char    *atom,        // atom type (e.g., 'od1')a	 : may be null
        char    *hydrogen,    // attached hydrogen (e.g., 'od1') : may be null.
        char    color,          // color 
        char    clouds,          // show dots: color=X --> ions; else --> sidechains
        Int4	width          // wireframe width
)
{
	item_list_type *item_node;
	NEW(item_node, 1, item_list_type);
	item_node->next=0;
	item_node->type='R';
	item_node->item.residue.res=res;
	item_node->item.residue.site=site;
	item_node->item.residue.chain=chain;
	item_node->item.residue.atom=atom;
	item_node->item.residue.hydrogen=hydrogen;
	item_node->item.residue.color=color;
	item_node->item.residue.clouds=clouds;
	item_node->item.residue.width=width;
#if 0   // DEBUG
if(item_node->item.residue.color == 'Y'){
	item_list_type *itm=item_node;
        fprintf(stderr,"%c%d%c\n",itm->item.residue.res, itm->item.residue.site,itm->item.residue.color);
}
#endif
	return item_node;
}

item_list_type *MakeItemNode(char *cmd,char type)
// type: 'V' = view; 'C' = clear; 'f' = file; 'F'=File
{
	item_list_type *item_node;
	NEW(item_node, 1, item_list_type);
	item_node->next=0; 
 	item_node->item.cmd=cmd;
	item_node->type=type; 
	return item_node;
}

item_list_type *MakeItemNodeView(Int4 *view)
{
	item_list_type *item_node;
	NEW(item_node, 1, item_list_type);
	item_node->next=0; 
	item_node->type='V'; 
	item_node->item.view=view;
	return item_node;
}

Int4	PutCmdItem(FILE *fp,item_type item,char type,Int4 file)
{
	switch (type){
	  case 'C': 
		if(file) fprintf(fp,"clear%d.",file); 
		else fprintf(fp,"clear."); 
	     break;
	  case 'E': 
		if(file) fprintf(fp,"echo%d=%s",file,item.cmd); 
		else fprintf(fp,"echo=%s",item.cmd); 
	     break;
	}
}

Int4	PutViewItem(FILE *fp,item_type item,char type,Int4 key_file)
{
	Int4 *v=item.view;
	fprintf(fp,"file%d:R=%d,%d,%d;T=%.2f,%.2f,%d.",key_file,
		v[0],v[1],v[2],(double)v[3]/100.0,(double)v[4]/100.0,v[5]);
}

Int4	PutMolItem(FILE *fp,item_type item,Int4 file)
{
	if(file) fprintf(fp,"%d",file);
	if(strcmp("HOH",item.molecule.mol) != 0) fprintf(fp,"!");
#if 0	// this should be unnecessary, as mol="[AF3]" in these cases...
	Int4 end=strlen(item.molecule.mol); end--;
	assert(end >= 0);
	if(isdigit(item.molecule.mol[end])){
	   fprintf(fp,"\[%s\]%d",item.molecule.mol,item.molecule.id);
	} else 
#endif
	fprintf(fp,"%s%d",item.molecule.mol,item.molecule.id);
	if(item.molecule.chain) fprintf(fp,"%c",item.molecule.chain);
	if(item.molecule.atom){ 
	   fprintf(fp,"_%s",item.molecule.atom);
	   if(item.molecule.hydrogen) fprintf(fp,"-%s",item.molecule.hydrogen);
	}
	if(item.molecule.clouds=='T') fprintf(fp,".{%c}",item.molecule.color);
	else if(item.molecule.clouds=='S') fprintf(fp,".(%c)",item.molecule.color);
	else fprintf(fp,".%c",item.molecule.color);
	// item.molecule.width;
}

Int4	PutTraceItem(FILE *fp,item_type item,Int4 file)
{
	if(file) fprintf(fp,"(%d)",file);
	fprintf(fp,"%d-%d",item.trace.start,item.trace.end);
	if(item.trace.chain) fprintf(fp,"%c",item.trace.chain);
	fprintf(fp,"%c",item.trace.backbone);
	fprintf(fp,"%c",item.trace.color);
	fprintf(fp,"%d",item.trace.width);
}

Int4	PutResItem(FILE *fp,item_type item,Int4 add,Int4 file)
{
	if(item.residue.atom) fprintf(fp," ");
	if(file) fprintf(fp,"%d",file);
	fprintf(fp,"%c%d",item.residue.res,item.residue.site+add);
	if(item.residue.chain) fprintf(fp,"%c",item.residue.chain);
	if(item.residue.atom){ fprintf(fp,"_%s",item.residue.atom);
	   if(item.residue.hydrogen){
		char str[9];
		if(item.residue.res == 'R' && strncmp("hh",item.residue.hydrogen,2) == 0){
		   fprintf(stderr,"WARNING: changing ARG%d hydrogen from %s to %s.\n",
			item.residue.site+add,item.residue.hydrogen,item.residue.hydrogen+1);
		   fprintf(fp,"-%s",item.residue.hydrogen+1);
		} else fprintf(fp,"-%s",item.residue.hydrogen);
	   }
	} if(item.residue.clouds == 'T') fprintf(fp,".{%c}",item.residue.color);
	else if(item.residue.clouds == 'S') fprintf(fp,".(%c)",item.residue.color);
	else fprintf(fp,".%c",item.residue.color);
	return (0);
}

Int4	HashCodeAtom(char *atom)
// either atom or 'hydrogen'
// '1'-'4' == 49..52 ACSII code
// n=110, o=111,h=104, c=99
// a,b,g, d,e,z
{
	Int4 factor=0;
	char	*str=atom,c;
	if(str==0) return 0;
	while(str[0] && isalnum(str[0])){
	   if(isdigit(str[0])){
		factor+=(Int4)(str[0] - '0')%10;
	   } else if(isalpha(str[0])){
		factor+=(Int4)(tolower(str[0]) - 'a')%12;
	   } str++;
	} return factor;
}

item_list_typ SortItems(item_list_typ list)
{ return SortItems(list,FALSE); }

item_list_typ SortItems(item_list_typ list,BooLean purge)
{
	item_list_typ node;
	Int4	Number=0,item;
	double	lastkey=0,key=-1.0;
	const double            semitiny = 0.0001;
	const double             tiny = 0.0000001;
	const double       verytiny =0.0000000001;
	const double    moretiny =0.0000000000001;
	const double mosttiny =0.0000000000000001;
	Int4	factor1,factor2;

	for(node=list; node != 0; node=node->next){ Number++; }
	vsih_typ H=MakeVsiHeap(Number);
	// 1. Insert items on heap
	for(node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  switch(node->type){
	   case 'R': key=(double) itm.residue.site; 
		if(itm.residue.atom){
		    factor1=HashCodeAtom(itm.residue.atom);
		    if(itm.residue.hydrogen){
		         factor2=HashCodeAtom(itm.residue.hydrogen);
			 key +=verytiny*(double)factor1 + moretiny*(double)factor2;  
		    } else key +=semitiny*(double)factor1;	// if single atom
		}
		if(key==lastkey) key+= mosttiny;
			break;
	   case 'M': key += mosttiny; break; // keep previous key...
	   case 'T': key=(double)itm.trace.start; key -= 0.1; break;
	   case 'V': key=-1; break;	// Put at top of file.
	   default:  key += mosttiny; // keep right after previous key...
		break; 
	  }
	  if(node->file) key += tiny*(double)node->file;
	  lastkey=key;
	  InsertVsiHeap(node,key,H);
	}
	// 2. Pop items from heap while linking into a list
	item_list_typ head=0,last;
	while(node=DelMinVsiHeap(&key, &item,H)){
	   if(head==0){ head=last=node; }
	   else { last->next=node; last=node; }
	   last->next=0;
	}
	NilVsiHeap(H);
	BooLean rm_something;
	if(purge) do { 
	  rm_something=FALSE;
	  for(last=0,node=head; node != 0; node=node->next){
	   // WARNING: this only eliminates adjacent redundant residue nodes.
	   if(IdenticalNodes(node,last)){
		last->next=node->next; 
		FreeNode(node); node=last;
		rm_something=TRUE;
	   } else last=node;
	  }
        } while(rm_something);
	return head;
}

item_list_typ SortItemsByChains(item_list_typ list)
// maintain the same ordering but sort by chain...
{
	item_list_typ node;
	Int4	Number,J,tmp,item;
	double	lastkey=0,key=-1.0;
	const double             tiny = 0.0001;
	const double       verytiny =0.0000001;
	const double    moretiny =0.0000000001;
	const double mosttiny =0.0000000000001;

	for(Number=0,node=list; node != 0; node=node->next){ Number++; }
	vsih_typ H=MakeVsiHeap(Number);
	// 1. Insert items on heap
	for(J=1,node=list; node != 0; node=node->next,J++){
	  item_type itm = node->item;
	  switch(node->type){
	   case 'R': 
		if(itm.residue.chain){
		    tmp=(Int4)itm.residue.chain;
		    key=(double)tmp + (double)J * verytiny;
		} else key=(double)J * verytiny;
		break;
	   case 'M': 
		if(itm.molecule.chain){
		    tmp=(Int4)itm.molecule.chain;
		    key=(double)tmp + (double)J * verytiny;
		} else key=(double)J*verytiny;
		break; 
#if 1
	   case 'T': 
		if(itm.trace.chain){
		    tmp=(Int4)itm.trace.chain;
		    key=(double)tmp + (double)J * verytiny;
		// 	key=(double)itm.trace.start; key -= 0.1; 
		} else key=(double)J*verytiny;
		break;
#endif
	   case 'V': key=-1; break;	// Put at top of file.
	   default:  key = lastkey + moretiny; // keep right after previous key...
		break; 
	  }
	  lastkey=key;
	  InsertVsiHeap(node,key,H);
	}
	// 2. Pop items from heap while linking into a list
	item_list_typ head=0,last;
	while(node=DelMinVsiHeap(&key, &item,H)){
	   if(head==0){ head=last=node; }
	   else { last->next=node; last=node; }
	   last->next=0;
	}
	NilVsiHeap(H);
	return head;
}

void	FreeNode(item_list_typ node)
{
	if(node==0) return;
	item_type itm = node->item;
	switch(node->type){
	  case 'R':
	    if(itm.residue.atom) free(itm.residue.atom);
	    if(itm.residue.hydrogen) free(itm.residue.hydrogen);
	   break;
	  case 'M':
	    if(itm.molecule.atom) free(itm.molecule.atom); 
	    if(itm.molecule.hydrogen) free(itm.molecule.hydrogen);
	    if(itm.molecule.mol) free(itm.molecule.mol);
	   break;
	  case 'V': if(itm.view) free(itm.view); break;
	  case 'C': case 'E': if(itm.cmd) free(itm.cmd); break;
	  default: break;
	} if(node->comment) free(node->comment);
	free(node);
}

BooLean	IdenticalNodes(item_list_typ node1, item_list_typ node2)
// determine whether or not nodes are identical.
{
	if(node1==0 || node2==0) return FALSE;
	if(node1->type == 'R' && node2->type == 'R'){
	    item_type itm1 = node1->item,itm2=node2->item;
	    if(itm1.residue.clouds != itm2.residue.clouds) return FALSE;
#if 1	// for pymol class buttons...
	    if(itm1.residue.color != itm2.residue.color) return FALSE;
#endif
	    if(itm1.residue.site == itm2.residue.site &&
		itm1.residue.res == itm2.residue.res &&
		itm1.residue.chain == itm2.residue.chain){
#if 0
if(itm1.residue.res == 'T' && itm1.residue.site == 17){
	PutResItem(stderr,itm1,0,0); fprintf(stderr,"\n");
	PutResItem(stderr,itm2,0,0); fprintf(stderr,"\n");
	if(itm1.residue.atom==0 && itm2.residue.atom==0) fprintf(stderr,"... are identical!\n");
	else if(strcmp(itm1.residue.atom,itm2.residue.atom) == 0) fprintf(stderr,"... atoms identical!\n");
	else if(strcmp(itm1.residue.hydrogen,itm2.residue.hydrogen) == 0) fprintf(stderr,"... hydrogens identical!\n");
}
#endif
		  if(itm1.residue.atom==0 && itm2.residue.atom==0) return TRUE;
		  if(itm1.residue.atom == 0) return FALSE;
		  if(itm2.residue.atom == 0) return FALSE;
		  if(strcmp(itm1.residue.atom,itm2.residue.atom) != 0) return FALSE;
		  if(itm1.residue.hydrogen == 0 && itm2.residue.hydrogen==0) return TRUE;
		  if(itm1.residue.hydrogen==0) return FALSE;
		  if(itm2.residue.hydrogen==0) return FALSE;
		  if(strcmp(itm1.residue.hydrogen,itm2.residue.hydrogen) == 0) return TRUE;
		  else return FALSE;
	    }
	     
	} else if(node1->type == 'M' && node2->type == 'M'){
	    item_type itm1 = node1->item,itm2=node2->item;
	    if(itm1.molecule.id == itm2.molecule.id &&
		itm1.molecule.chain == itm2.molecule.chain){
		  assert(itm1.molecule.mol && itm2.molecule.mol);
		  if(strcmp(itm1.molecule.mol,itm2.molecule.mol) != 0) return FALSE;
		  if(itm1.molecule.atom==0 && itm2.molecule.atom==0) return TRUE;
		  if(itm1.molecule.atom==0) return FALSE;
		  if(itm2.molecule.atom==0) return FALSE;
		  if(strcmp(itm1.molecule.atom,itm2.molecule.atom) != 0) return FALSE;
		  if(itm1.molecule.hydrogen==0 && itm2.molecule.hydrogen==0) return TRUE;
		  if(itm1.molecule.hydrogen==0) return FALSE;
		  if(itm2.molecule.hydrogen==0) return FALSE;
		  if(strcmp(itm1.molecule.hydrogen,itm2.molecule.hydrogen) == 0) return TRUE;
		  else return FALSE;
	    } else return FALSE;
	} return FALSE;
}

Int4	PrintItems(FILE *fp, item_list_type *list,char *filename, char *FileComment,
	char key_chain, Int4 key_file,Int4 add,Int4 file,BooLean verbose)
{
	char	*tmpFN[3],tmpKC[3],*tmpFC[3];
	Int4	tmpKF[3];

	tmpFN[1]=filename; tmpKC[1]=key_chain; tmpKF[1]=key_file; tmpFC[1]=FileComment;
	return PrintItems(fp,list,1,tmpFN,tmpFC,tmpKC,tmpKF,add,file,verbose);
}

Int4	PrintItems(FILE *fp, item_list_type *list,Int4 n, char **filename, 
	char **FileComment, char *key_chain, Int4 *key_file,Int4 add,Int4 file,
	BooLean verbose)
// Int4 file term is for changing the file number for a single file.
{
	item_list_type *node,*next,*last;
	char last_type=0;

	for(Int4 i=1; i<=n; i++){
	   assert(filename[i]);
	   if(key_chain[i]==0) key_chain[i]=' ';
	   fprintf(fp,"File%d=%s:%c   ",key_file[i],filename[i],key_chain[i]);
	   if(FileComment[i]) fprintf(fp,"%s\n",FileComment[i]);
	   else fprintf(fp,"\n");
	} fprintf(fp,"\n");

	for(node=list; node != 0; last_type = node->type, node=node->next){
	  item_type itm = node->item;
	  if(last_type && last_type != node->type) fprintf(fp,"\n");
	  switch(node->type){
	   case 'R': 
		if(file) PutResItem(fp,itm,add,file); 
		else PutResItem(fp,itm,add,node->file);
		break;
	   case 'M': 
		if(file) PutMolItem(fp,itm,file); 
		else PutMolItem(fp,itm,node->file);
		break;
	   case 'T': 
		if(file) PutTraceItem(fp,itm,file); 
		else PutTraceItem(fp,itm,node->file); 
		break;
	   case 'V': 
		if(n==1 && key_file[1]) PutViewItem(fp,itm,'V',key_file[1]); 
		else PutViewItem(fp,itm,'V',node->file);
		break;
	   default: 
		if(file) PutCmdItem(fp,itm,node->type,file); 
		else PutCmdItem(fp,itm,node->type,node->file);
		break;
	  }
	  if(node->comment && verbose) fprintf(fp,"  %s\n",node->comment); 
	  else fprintf(fp,"\n");
	} return (0);
}

Int4	NumberAtomsItem(item_type item,char type)
{
	switch(type){
	   case 'R':
		if(item.residue.hydrogen) return 2;
		else if(item.residue.atom) return 1;
		else return 0;
	    break;
	   case 'M': 
		if(item.molecule.hydrogen) return 2;
		else if(item.molecule.atom) return 1;
		else return 0;
	    break;
	   default: return 0; break;
	}
}

dots_list_type *MakeDotsList(char *str,char color)
{
	dots_list_type *dots_node; 
	NEW(dots_node,1,dots_list_type);
	dots_node->str=str;
// fprintf(stderr,"********** %s\n",str);
	dots_node->next=0;
	dots_node->color=color;
	return dots_node;
}

#if 0	// moved to vsi_ras.cc
Int4    PrintDotsItems(FILE *fp, item_list_type *list, ras_typ *ras)
{
	item_list_type *node,*next;
	dots_list_type *dlist=0,*dnode;
	char	*str;
	Int4	N,i,j,k,num;

	// 1. get list of things for dot clouds...
	for(num=0,node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  switch(node->type){
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
	   case 'R':
	    if(itm.residue.clouds == 'T'){
	      num++;
	      str=ras->DefineResCloud(fp,itm.residue.res,itm.residue.site,
				itm.residue.atom,itm.residue.hydrogen,
				itm.residue.chain);
	      if(dlist==0){ 
		dnode=MakeDotsList(str,itm.residue.color); dlist=dnode;
	      } else {
		dnode->next=MakeDotsList(str,itm.residue.color); dnode=dnode->next;
	      }
	    } break;
	   case 'M':
	    if(itm.molecule.clouds == 'T'){
	      num++;
	      str=ras->DefineMolCloud(fp,itm.molecule.mol, itm.molecule.id,
			itm.molecule.atom,itm.molecule.hydrogen,itm.molecule.chain);
	      if(dlist==0){ 
		dnode=MakeDotsList(str,itm.molecule.color); dlist=dnode;
	      } else {
		dnode->next=MakeDotsList(str,itm.molecule.color);
	     	dnode=dnode->next;
	      }
	    } break;
	   default: 
	    break;
	  }
	}
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
	return N;
}

Int4	PrintItemsRasmol(FILE *fp, BooLean look_in_datadir, item_list_type *list,
	char *filename, char key_chain,BooLean cis_trans_pro,
	Int4 wire_width,Int4 spacefill)
{
	item_list_type *node;
	if(filename == 0){
	   fprintf(stderr,"\"File#=pdbfilename:Chain\" not designated!\n"); 
	   print_error("**************** Fatal error! ******************");
	}
	ras_typ *ras=0;
        ras = new ras_typ(filename,cis_trans_pro);

	char su=key_chain,chn,c;
        if(islower(su)){ su=toupper(su); }
	double factor=1.5;

        if(su == ' ') chn=0; else chn=su;
        // if(su != '-') ras->PrintHEADER(fp,chn,wire_width); else 
	ras->PrintHEADER(fp,look_in_datadir,0,wire_width);
	for(node=list; node != 0; node=node->next){
	  item_type itm = node->item;
	  Int4 *v;
	  switch(node->type){
	   case 'R':
	     if(itm.residue.clouds == 'S') spacefill=9999;
	     else spacefill=(Int4) ceil((double)itm.residue.width * factor);
	     ras->PutResidueItem(fp,itm.residue.res, itm.residue.site,
				itm.residue.atom, itm.residue.hydrogen,
				itm.residue.chain, itm.residue.color,
				itm.residue.width, spacefill);
	    break;
	   case 'M': 
	     if(itm.molecule.clouds == 'S') spacefill=9999;
	     else spacefill=(Int4) ceil((double)itm.molecule.width * factor);
	     ras->PutMoleculeItem(fp,itm.molecule.mol, itm.molecule.id,
			itm.molecule.atom, itm.molecule.hydrogen,
			itm.molecule.chain, itm.molecule.color,
			itm.molecule.width, spacefill);
	    break;
	   case 'T': 
		if(itm.trace.chain==0) c=chn; else c=itm.trace.chain;
		if(c==' ') c=0;
		ras->PrintTrace(fp,itm.trace.start,itm.trace.end,c,
				itm.trace.color,itm.trace.backbone,itm.trace.width);
	    break;
	   case 'V': 
		v=itm.view; 
		ras->PrintView(fp,v[0],v[1],v[2],v[3],v[4],v[5]);
	    break;
	   default: ras->PrintCommand(fp,itm.cmd,node->type); break;
	  }
	}
	PrintDotsItems(fp,list,ras);
	ras->PrintTAIL(fp);
	delete ras;
	return (0);
}
#endif

