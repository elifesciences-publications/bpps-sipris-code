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

#include "tree2hpt.h"
#include "random.h"
#include "dheap.h"
#include "set_typ.h"

int	tree2hpt(int argc,char *argv[]){ return tree2hpt(0,argc,argv,20000); }

int	tree2hpt(FILE *fp, int argc,char *argv[]) { return tree2hpt(fp,argc,argv,20000); }

int     tree2hpt(FILE *fp, int argc,char *argv[], Int4 NumRandom) 
	{ return tree2hpt(fp,0,argc,argv,NumRandom); }

// e.g. input string: (((7,((15,16):11,12,(17,18):13,14):8):3,4,(9,10):5,6):1,2):0;

#if 0
//**************************************************************************
   Strategy for a tree based search:
//--------------------------------------------------------------------------
0. Create directories for profiles:
  0a. TreeStructure for Directories mirrors tree structure for Template.
  0b. change ~/sbin/MkGAPMAPS to ~/sbin/MkMAPGAPS to convert this into maps format (define?).
  0c.

1. Read in tree (if file provided as input).
 1a. Define subtrees that need to be converted at each level.
 1b. Everytime a parent node is reached Convert the child alignments to match parent node
       using the SubTreeTemplate
 1c. Also convert the SubTreeTemplate cma file and add to template...How?

2. Create mgs object and pass in tree for search...

3. For tree, call ConvertViaTemplate( ) recursively...

For each lowest level subtree, call:

conversionViaTemplateCMSA(cma_typ TemplateCMA, cma_typ *IN_CMA);

Where TemplateCMA is SubtreeTemplateCMA.

Alsp need to convert SubtreeTemplateCMA to

//**************************************************************************
#endif

#define USAGE_START "Usage: tree2hpt infile [options] \n\
    -b          	print tree in binary format (default: FD-table format)\n\
    -depth=<int> 	print tree as a hpt, but only to depth <int> from root\n\
    -fuse=<int> 	Fuse node <int> into its parent\n\
    -leaf=<int>:<int>	split off leaf node <int1> from node <int2>\n\
    -N          	print tree in Newick format\n\
    -O=<str>    	output filename=<str>\n\
    -parent=<int>:<int> split off parent node <int1> from node <int2>\n\
    -R=<int>          	Run <int> random operations on input tree\n\
    -seed=<int>  	provide seed for random number generator\n\
    -sub=<str>  	print subtree from the node with name=<str>\n\
    -x          	dummy\n\n"

int     tree2hpt(FILE *fp, FILE *iFP, int argc,char *argv[], Int4 NumRandom)
{
	Int4	i,J,n,r,arg,N=0;
	FILE	*outfp=stdout,*ifp=0;
	Int4    DEPTH_FROM_ROOT=0;
	char    *SubTreeName=0;
	char	mode='H',str[200];
	BooLean PrintAsNewick=FALSE;
	Int4	FuseID=-1,targetID=-1,newID=-1,random=0;
	UInt4	seed=18364592;
	Int4	time1=time(NULL);

	if(argc < 2) print_error(USAGE_START);
	if(fp != 0) outfp = fp;
	if(iFP != 0) ifp = iFP;
        for(arg = 2; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'b': mode='B'; break;
             case 'd': 
		if(sscanf(argv[arg],"-depth=%d",&DEPTH_FROM_ROOT)==1){
			if(DEPTH_FROM_ROOT < 1) print_error(USAGE_START);
		} break;
             case 'f': 
		if(sscanf(argv[arg],"-fuse=%d",&FuseID)==1){
			if(FuseID < 0) print_error(USAGE_START);
			mode='f';
		} else print_error(USAGE_START);
		 break;
             case 'l': 
		if(sscanf(argv[arg],"-leaf=%d:%d",&targetID,&newID)==2){
			if(targetID < 0) print_error(USAGE_START);
			if(newID <= 0) print_error(USAGE_START);
		} else print_error(USAGE_START);
		mode='l'; break;
             case 'N': mode='N'; break;
             case 'O':
	        if(argv[arg][2] == '=' && isprint(argv[arg][3])){
			outfp=open_file(argv[arg] + 3,"","w");
		} else print_error(USAGE_START);
		break;
             case 'p': 
		if(sscanf(argv[arg],"-parent=%d:%d",&targetID,&newID)==2){
			if(targetID < 0) print_error(USAGE_START);
			if(newID <= 0) print_error(USAGE_START);
		} else print_error(USAGE_START);
		mode='p'; break;
             case 'R': 
#if 0
		if(sscanf(argv[arg],"-R=%d",&random)==1){
			if(random < 1) print_error(USAGE_START);
		} else print_error(USAGE_START); 
#else
		random = IntOption(argv[arg],'R',1,100000,"tree2hpt -R option input error");
#endif
		mode='R'; break;
             case 's': 
		if(sscanf(argv[arg],"-sub=%s",str)==1){
		   SubTreeName=AllocString(str);
		} else if(sscanf(argv[arg],"-seed=%u",&seed)!=1){
			fprintf(stderr,"input seed = %u\n",seed);
			print_error(USAGE_START);
		} break;
             case 'x': break;
             default : print_error(USAGE_START);
           }
	  }else print_error(USAGE_START);
	}
	if(seed == 18364592)  seed = (UInt4) time(NULL)/2;
        if(seed != 0) sRandom(seed);
	fprintf(stderr,"\n=========== Parsing Newick file %s ============\n",argv[1]);
	btn_typ	*btn=0;
	if(mode != 'R'){
	   if(ifp==0){ ifp = open_file(argv[1],"","r");  }
	   btn=TreeToHpt(ifp); fclose(ifp);
	} // the above uses 'else' statement below...
	if(mode=='R'){
	  dh_type dH=dheap(random+2,4);		// new node heap.
	  dh_type tdH=dheap(random+2,4);	// target node heap.
	  btn_typ *btn= new btn_typ(1,0);	// 1 == the root.
	  insrtHeap(1,(keytyp) Random(),tdH);
	  for(i=2; i <= random; i++){ insrtHeap(i,(keytyp) Random(),dH); }
	  // set_typ Set(); // want sequence set operations as well?
	  do {
		// Int4 N=btn->CountNodes(0);	// number of nodes: starts with 0
		targetID=delminHeap(tdH);
	  	insrtHeap(targetID,(keytyp) Random(),tdH);
		if(targetID == 1) r = 0; else r=random_integer(3);
		switch(r){
		case 0: 	// add a leaf node.
		 {
		   newID=delminHeap(dH);  
		   btn_typ *X= new btn_typ(newID,0);
		   btn_typ *fnd=btn->SplitOffLeaf(X,targetID);
		   if(fnd != 0){
			insrtHeap(newID,(keytyp) Random(),tdH); // ID now in tree.
		   	// btn->PrintNewick(outfp); fprintf(outfp,";\n\n"); 
		   } else { insrtHeap(newID,(keytyp) Random(),dH); delete X; }
		 } break;
		case 1: 	// add a parent node.
		 {
		   newID=delminHeap(dH);  
		   btn_typ *X= new btn_typ(newID,0);
		   btn_typ *fnd=btn->SplitOffLeaf(X,targetID);
		   if(fnd !=0){
			insrtHeap(newID,(keytyp) Random(),tdH); // ID now in tree.
		   	// btn->PrintNewick(outfp); fprintf(outfp,";\n\n"); 
		   } else { insrtHeap(newID,(keytyp) Random(),dH); delete X; }
		 } break;
		case 2: 	// fuse a node to its parent.
		 {
		   btn_typ *X=btn->Fuse(targetID);
		   if(X!=0){
			i=X->Ident(); insrtHeap(i,(keytyp) Random(),dH);
			rmHeap(i,tdH); // remove from target heap.
			delete X;
		   	// btn->PrintNewick(outfp); fprintf(outfp,";\n\n"); 
		   } 
		 } break;
		default: print_error("this should not happen"); break;
		}
	   } while(!emptyHeap(dH)); Nildheap(dH);
	   btn->PrintNewick(outfp); fprintf(outfp,";\n\n"); 
	} else if(mode=='l'){
		btn_typ *X= new btn_typ(newID,0);
		btn_typ *fnd=btn->SplitOffLeaf(X,targetID);
		if(fnd==0) print_error("Target node not found");
		btn->PrintHPT(stderr,DEPTH_FROM_ROOT,NumRandom);
		btn->PrintNewick(outfp); fprintf(outfp,";\n\n");
	} else if(mode=='p'){
		btn_typ *X= new btn_typ(newID,0);
		X->PrintNewick(stderr); fprintf(stderr,";\n\n");
		btn_typ *fnd=btn->SplitOffParent(X,targetID);
		if(fnd==0) print_error("Target node not found");

		fnd->PrintNewick(stderr); fprintf(stderr,";\n\n");
		X->PrintNewick(stderr); fprintf(stderr,";\n\n");

		btn->PrintHPT(stderr,DEPTH_FROM_ROOT,NumRandom);
		btn->PrintNewick(outfp); fprintf(outfp,";\n\n");
	} else if(mode=='f'){
		btn->PrintNewick(outfp); fprintf(outfp,";\n\n");
		btn->PrintHPT(outfp,DEPTH_FROM_ROOT,NumRandom);
		btn_typ *X=btn->Fuse(FuseID);
		if(X){
			X->PrintNewick(stderr); fprintf(stderr,";\n\n");
			delete X;
			btn->PrintNewick(outfp); fprintf(outfp,";\n\n");
			fflush(outfp);
			btn->PrintHPT(outfp,DEPTH_FROM_ROOT,NumRandom);
		} else fprintf(stderr,"Node %d not found\n\n",FuseID);
	} else if(mode=='B'){
		btn->PrintBinary(outfp); fprintf(outfp,";\n\n");
		// N=btn->Print(stderr,0);
	} else if(mode=='N'){
		btn->PrintNewick(outfp); fprintf(outfp,";\n\n");
	} else{
	    if(SubTreeName){ N=btn->PrintSubHPT(outfp,SubTreeName,DEPTH_FROM_ROOT); }
	     else N=btn->PrintHPT(outfp,DEPTH_FROM_ROOT,NumRandom);
	} if(outfp != stdout && outfp != fp) fclose(outfp);
	// fprintf(stderr,"; %d nodes\n",N);
	delete btn;
	fprintf(stderr,"\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	return N;
}

