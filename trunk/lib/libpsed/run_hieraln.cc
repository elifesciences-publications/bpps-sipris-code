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

#include "lha_typ.h"

#define USAGE_START "Usage: hieraln <file_prefix> [options] \n\
      requires <prefix>.mma, <prefix>.sets and <prefix>.hpt input files\n\
        -debug      - output intermediate files for debugging.\n\
        -realign    - realign input alignment using gismo++ (use on poor quality alignments only)\n\
        -Sto        - Output subalignments in a Stockholm formated file\n\
        -seed=<int> - seed random number generator\n\
        -h   	    - help\n\n"

#define PUBLIC_START "\
   Input: <prefix>.mma: multiple sequence alignment (MSA) for BPPS 1\n\
          <prefix>.chk: checkpoint file created in BPPS step 1.\n\
   Output: <prefix>_himsa.cma: hiMSA subgroup alignments in cma format.\n\
           <prefix>_himsa.tpl: subalignment template MSA.\n\
           <prefix>_himsa.sets: sequence subgroup assignments (binary file).\n\
           <prefix>_himsa.dft: depth first search tree used by MAPGAPS.\n\
           <prefix>_himsa.hpt: hierarchy file.\n\
           <prefix>_himsa.ph: hierarchy in Newick format.\n\
       Note: run 'mkmaps <prefix>_himsa' to create profiles for MAPGAPS searches\n\
   Options:\n\
          -Sto        - Output subalignments in Stockholm format\n\
          -seed=<int> - seed for random number generator\n\
   Reference:\n\
          Neuwald, A.F. & S. F. Altschul. 2016. Inference of Functionally-Relevant \n\
	    N-Acetyltransferase Residues Based on Statistical Correlations. \n\
            PLoS Comput. Biol. 12(12):e1005294. PMID: 28002465\n\
        \n"

#include "blosum62.h"

void	Foo(cma_typ cma)
// Test the MkTemplateOperation() routine...
{
	assert(nBlksCMSA(cma)==1);
	// char	dms_mode='S';
	// Int4    aa_per_io=500,aa_per_do=500,exp_ie=1,exp_de=1;
	char	dms_mode='F';
	Int4	aa_per_io=20,aa_per_do=150, exp_ie=1,exp_de=1;
	Int4	end,Len = LengthCMSA(1,cma);

	Int4	start=random_integer(Len)+1;
	Int4	len=random_integer(10) +3;
	if((start+len) >= Len) start=Len-len-1;
	if(start < 2) start=2; end=start+len;
	if(end >= Len) end = Len-1; 
	fprintf(stderr,"start=%d; end=%d; len=%d\n",start,end,Len);
	cma_typ	xcma=ConvertColsToInsertsCMSA(cma,1,start,end);

	Len= LengthCMSA(1,xcma);
	len=random_integer(10)+3; 
	Int4	offset=random_integer(5)-2;
	start=start+offset; if(start < 2) start=2; end=start+len;
	if(end >= Len) end = Len-1; if(start > end) start=end;
	fprintf(stderr,"start=%d; end=%d; len=%d\n",start,end,Len);
	cma_typ	zcma=ConvertColsToInsertsCMSA(xcma,1,start,end);

	Int4	sq=1,*Sites=GetPosSitesCMSA(sq,xcma);
	gsq_typ *gsq=gsqCMSA(sq,xcma);
	char *oper[3]; oper[2]=gsq->Operation(1,Sites,LengthsCMSA(xcma));  //
	Sites=GetPosSitesCMSA(sq,zcma); gsq=gsqCMSA(sq,zcma);
	oper[1]=gsq->Operation(nBlksCMSA(zcma),Sites,LengthsCMSA(zcma));  //
        fprintf(stderr,"oper[p]=%s\n",oper[1]);
        fprintf(stderr,"oper[c]=%s\n",oper[2]);
	gmb_typ	*xgmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,xcma,dms_mode);
	ssx_typ	*ssx=xgmb->RtnSSX();
	gmb_typ *zgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,zcma,dms_mode);
	ssx_typ *SSX=zgmb->RtnSSX();

	char	*operation=MkTemplateOperation(SSX,ssx,oper);
	fprintf(stderr,"opr_out=%s\n",operation);
	delete zgmb; delete xgmb; NilCMSA(xcma); NilCMSA(zcma);
}

static void PrintError(char SecretCode,char *prgm_name)
{
        if(SecretCode > 0){
                fprintf(stderr,"Usage: %s %d <prefix> [options] \n",prgm_name,(int)SecretCode);
                print_error(PUBLIC_START);
        } else { print_error(USAGE_START); }
}

int     run_hieraln(int argc, char *argv[]) { return run_hieraln(0,argc,argv); }

int     run_hieraln(char SecretCode,int argc, char *argv[])
// 0. Creating a MAPGAPS-searchable hierarchical alignment (for CDD).
#if 0   // Collapse Hpt and sets down to lineage only for each internal node!!!
	Adjust patterns!!!
        <file>_Set1.hpt <file>_Set1.sets
        <file>_Set2.hpt <file>_Set2.sets
#endif
{
	Int4	arg,i,j,n,x,y,z,nsets,percent_ident=80,time1=time(NULL),*Parent,parent,s,p;
	Int4	seed=18364592;
	char	str[200]; 
	BooLean	debug=FALSE,realign=FALSE,put_stockholm=FALSE; // realign=TRUE;
	cma_typ	TrueMainCMA,old_cma,new_cma=0,xcma,tcma,rcma,cma,CMA;
	set_typ	*NodeSet=0;

	//============ 1. Scan (and output) argument list. =============
	if(argc < 2) PrintError(SecretCode,argv[0]);
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') PrintError(SecretCode,argv[0]);
	   switch(argv[arg][1]) {
		case 'd':
		   if(SecretCode==0) PrintError(SecretCode,argv[0]);
		   if(strcmp("-debug",argv[arg])==0){ debug=TRUE; }
                   else PrintError(SecretCode,argv[0]); break;
		case 'S': 
		   if(strcmp("-Sto",argv[arg])==0) put_stockholm=TRUE;
		   else PrintError(SecretCode,argv[0]); break;
		case 's': 
		   if(sscanf(argv[arg],"-seed=%d",&seed)!=1) PrintError(SecretCode,argv[0]); break;
		case 'r':
		   if(SecretCode==0) PrintError(SecretCode,argv[0]);
	           if(strcmp("-realign",argv[arg])==0){ realign=TRUE; }
                   else PrintError(SecretCode,argv[0]); break;
		default: PrintError(SecretCode,argv[0]); break;
	   }
	}
	FILE *fp = open_file(argv[1],".cmd","w");
        for(i = 0; i < argc; i++) { fprintf(fp,"%s ",argv[i]); }
        if(seed == 18364592) {  // not provided by user
                seed=(UInt4) time(NULL)/2; // seed = (UInt4) time(NULL);
                fprintf(fp,"-seed=%d",seed);
        } fprintf(fp,"\n"); fclose(fp);
	const Int4 NumIter=1;	// one iteration appears to be enough...

	//============== 2. Read (and check) MSA & hierarchy (BPPS) input files. ===============
	a_type AB=MkAlphabet(AMINO_ACIDS,GBLAST_BLOSUM62,SREL26_BLSM62);
	fp=open_file(argv[1],".mma","r"); TrueMainCMA=ReadCMSA(fp,AB); fclose(fp); fp=0;
// ExtendFakeToRealCMSA(TrueMainCMA); 
	// Foo(TrueMainCMA); exit(1); // test routine...
	if(realign)	// realign main alignment...to get rid of MAPGAPS problems.
          { rcma=PurgeSampleReAlign(argv[1],TrueMainCMA); NilCMSA(TrueMainCMA); TrueMainCMA=rcma; }
     hpt_typ *Hpt=0;
     if(SecretCode > 0){
        fp=open_file(argv[1],".chk","r"); 
        fprintf(stderr,"Reading checkpoint file\n");
        Int4    NumSMA=0,sd=0;
        NodeSet=ReadSets(fp,nsets);
        if(fgets(str,200,fp) == NULL) print_error("*.chk file input error 1.");
        if(sscanf(str,"seed=%d; N=%d.",&sd,&NumSMA) != 2) print_error("*.chk file input error 2.");
        seed=sd;
        cma_typ *chk_sma; NEW(chk_sma, NumSMA+3, cma_typ);
        for(i=1; i<=NumSMA; i++){
               if((chk_sma[i]=ReadCMSA(fp,AB))==NULL) print_error("*.chk file input error 3.");
	       TotalNilCMSA(chk_sma[i]);
        } free(chk_sma);
        Hpt=new hpt_typ(fp); fclose(fp);
        if(nsets != Hpt->NumSets()) print_error("FATAL: *.hpt and *.sets files are inconsistent");
     } else {
	fp = open_file(argv[1],".sets","r"); NodeSet=ReadSets(fp,nsets); fclose(fp);
	fp = open_file(argv[1],".hpt","r"); Hpt = new hpt_typ(fp); fclose(fp);
     }
	sRandom(seed);
	// Hpt->Put(stderr); exit(1);	// test out the gcc-4.9.0 sscanf bug with -O options.
	if(Hpt->NumSets() < 3) print_error("FATAL error: Single node hierarchy");
	assert(nsets == Hpt->NumSets()); assert(Hpt->IsTree(Parent)); // Find parent for each node.
	// assert(Hpt->PutDFTree(stderr)); exit(1);

#if 0	// Debug...
	fp = open_file(argv[1],"_new.mma","w"); 
	///// set_typ setX=CopySet(NodeSet
	for(s=1; s <= Hpt->NumSets(); s++){	// skip root and reject sets!!!
		PutInSetCMSA(fp,NodeSet[s],TrueMainCMA);
	} fclose(fp); exit(1);
#endif
	//============== 3. Create NodeCMA array and root consensus sequence. ===============
	cma_typ	*NodeCMA; NEW(NodeCMA,Hpt->NumSets()+5, cma_typ); 
	e_type	*CsqE; 		NEW(CsqE,Hpt->NumSets()+3,e_type); 
	
	//============= 4. Align sequences within each subgroup. ==================
/******************************************************************
  No need to realign entire set! So no need to AutoPurge!!

  Try a fast purge procedure that uses random sample of a tenth of sequences.
  1. Sample at random say a 1/100 of the sequences.
  2. Find sequence < 95% identity.
  3. 
or
  1. Purge Main set and use the same seqs for other sets?

  * Fix Deletion-to-insertion issue...

  After lengthening each alignment and creating template alignments:
  1. Create a new sub_alignment (cma) for each consisting only of seqs belonging to that node.
  2. Add original sub_cma gsq aligned sequences to new sub_alignment file.
  3. Use Gambit to sampling in previously purged sequences for that node.
  4. For leaf nodes that weren't 
 *****************************************************************/
	//============= 4... Set GISMO parameters. ==================
	gmb_typ	*gmb=0,*xgmb=0;

	//============= 4... For each subtree group... ==================
	dh_type dH=dheap(Hpt->NumSets()+3,4); insrtHeap(1,1,dH); 
	char	**Oper; NEWP(Oper,Hpt->NumSets()+3, char);  // Save indels introduced for each set.
	set_typ	*ChildSet;	NEW(ChildSet,Hpt->NumSets()+3,set_typ);
	set_typ	**RtnSet=0; NEWP(RtnSet,Hpt->NumSets()+3,set_typ);
#if 0	// This is much faster...but doesn't work?
	set_typ subtree=Hpt->MkSubTreeSet(1); PutSet(stderr,subtree);
	// RtnSet[1]=MkSubSets(NodeSet,Hpt->NumSets(),subtree); NilSet(subtree);
	RtnSet[1]=MkSubSets(NodeSet,Hpt->NumSets()-1,subtree); NilSet(subtree);
#if 1
	// fp=tmpfile(); 
	fp=open_file("junk",".cma","w");  // main file.
	PutConsensusCMSA(fp,TrueMainCMA);
	PutInSetCMSA(fp,NodeSet[1],TrueMainCMA);
	fclose(fp);
exit(1);	
#else
	NodeCMA[1]=CopyCMSA(TrueMainCMA); // Don't purge; Main set not realigned!!!
#endif
#else	// this is much slower...
	assert((NodeCMA[1]=MkSubSetCMSA(NodeSet,Hpt,1,FALSE,TrueMainCMA,&RtnSet[1])) != 0);
#endif
	RenameCMSA(Hpt->SetName(1),NodeCMA[1]); CsqE[1]=MkConsensusCMSA(NodeCMA[1]); 
	sprintf(str,"%s consensus",Hpt->SetName(1)); ChangeInfoSeq(str,CsqE[1]);
#if 0	// Run this to fix temporary issue with older output... Use *_himsa.sets and main cma.
   {
	sprintf(str,"%s_xxxx",argv[1]); 
	fp=open_file(str,".cma","w");  // main file.
	PutWithCsqCMSA(fp,CsqE[1],NodeCMA[1]); // skip root & reject sets!
	fclose(fp); 

	sprintf(str,"%s_xxxx",argv[1]); fp=open_file(str,".sets","w"); 
        FILE *xfp=open_file("junk",".sets","r");
	Int4	ii,NN;
     set_typ *xsets=0,**xSets=0; NEWP(xSets,Hpt->NumSets()+3,set_typ);
     for(ii=1; (xsets=ReadSets(xfp,NN)) != 0; ii++){
	if(ii==1) xSets[ii]=RtnSet[1]; else xSets[ii]=xsets;
        Int4 total=0,x;
        for(Int4 i=1; i <= NN; i++){
                x = CardSet(xSets[ii][i]);
                fprintf(stderr,"%d. %d items\n",i,x);
                total += x;
        } WriteSets(fp,NN,xSets[ii]);
     } fclose(fp); fclose(xfp);
   } exit(1);
#endif
	//============= 4... Recursively expand each subgroup alignment. ==================
	//========== modify to start with parent alignment to expand each set. ===========
	do {
	    parent = delminHeap(dH); assert(NodeCMA[parent] != 0);
	    ChildSet[parent]=Hpt->MkSubTreeSet(parent); DeleteSet(parent,ChildSet[parent]);
	    fprintf(stderr,"************** parent = %d **************\n",parent);
	    for(s=2; s < Hpt->NumSets(); s++){	// skip root and reject sets!!!
		if(Parent[s] != parent) continue; assert(Oper[s] == 0);
		insrtHeap(s,(keytyp)s,dH); // Store child on heap as next parent; save its cma.

		//========  4a. Generate subgroup alignment (without concensus). ========
		//========      Returns 0 if set is empty.      ===========
		if((cma=MkSubSetCMSA(NodeSet,Hpt,s,FALSE,TrueMainCMA,&RtnSet[s]))==0) continue;

		//========  4b. Extend subgroup alignment; save as NodeCMS[s]. ========
#if 0
	FILE *dbfp=open_file("debug_in",".cma","w"); PutCMSA(dbfp,cma); fclose(dbfp);
	rcma=QuickOptimizeCMA(cma,s,Oper);	// save add/rm column operation.
	dbfp=open_file("debug_out",".cma","w"); PutCMSA(dbfp,rcma); fclose(dbfp); exit(1);
#else
		rcma=QuickOptimizeCMA(cma,s,Oper);	// save add/rm column operation.
#endif
		assert(rcma != 0); RenameCMSA(Hpt->SetName(s),rcma); NodeCMA[s]=rcma; 

		//========  4c. Make a consensus for template. ========
	        sprintf(str,"%s consensus",Hpt->SetName(s)); 
		CsqE[s]=MkConsensusCMSA(rcma); ChangeInfoSeq(str,CsqE[s]);
		if(debug) PutHieralnDebug("",1,argv[1],Hpt->SetName(s),CsqE[s],rcma);
		NilCMSA(cma); // No longer need input cma; shares data set, so don't TotalNilCMSA();
            } 
	} while(!emptyHeap(dH)); Nildheap(dH); 

	//========  5. Make the template alignments using GAMBIT. ========
	cma_typ	*TmplCMA=MkTemplateCMAs(Hpt,NodeCMA,ChildSet,CsqE,Oper,AB);

	//===========  6. Output MAPGAPS ('root') input files. ===========
	// sprintf(str,"%s_%s",argv[1],Hpt->SetName(1));
	sprintf(str,"%s_himsa",argv[1]);
	fp=open_file(str,".dft","w"); assert(Hpt->PutDFTree(fp)); fclose(fp); 
	// assert(Hpt->PutDFTree(stderr)); exit(1);

	fp=open_file(str,".cma","w");  // main file.
	for(s=1; s < Hpt->NumSets(); s++) PutWithCsqCMSA(fp,CsqE[s],NodeCMA[s]); // skip root & reject sets!
	fclose(fp); 
#if 1
	if(put_stockholm){
	   Int4	aln,Num=0; 
	   FILE *ofp=open_file(str,".sto","w"); PutStockholmCMSA(ofp,TmplCMA[1]);
	   fp=open_file(str,".cma","r");  // main file.
	   cma_typ *IN_CMA=MultiReadCMSA(fp,&Num,AB); fclose(fp);
	   for(aln=1; aln <= Num; aln++){
		PutStockholmCMSA(ofp,IN_CMA[aln]); TotalNilCMSA(IN_CMA[aln]);
	   } fclose(ofp); free(IN_CMA);
	}
#endif
#if 0	// this is not needed...use tpl below.  Delete later.
	sprintf(str,"%s_root",argv[1]);
	fp=open_file(str,".tpl","w"); PutCMSA(fp,TmplCMA[1]); fclose(fp); 
#endif

#if 1	//========= Add an output file containing realigned input sequences... =======
	Int4 Number;
        // fp = open_file(argv[1],".cma","r");
	cma_typ	*IN_CMA; NEW(IN_CMA, Hpt->NumSets() +3, cma_typ);
        for(s=1; s < Hpt->NumSets(); s++){ assert(NodeCMA[s]); IN_CMA[s]=AddThisCsqCMSA(CsqE[s],NodeCMA[s]); }
#if 0
fp=open_file(argv[1],"_debug.cma","w");  // main file.
for(s=1; s < Hpt->NumSets(); s++){ PutCMSA(fp,IN_CMA[s]); } fclose(fp); exit(1);
#endif
        cma_typ *OUT_CMA=ConversionViaTemplateCMSA3(TmplCMA[1],IN_CMA); 
#if 0
fp=open_file(argv[1],"_debug.cma","w");  // main file.
for(s=1; s < Hpt->NumSets(); s++){ PutCMSA(fp,OUT_CMA[s]); } fclose(fp); exit(1);
#endif
	Number = Hpt->NumSets() - 1;
	// for(s=1; s < Hpt->NumSets(); s++){ PutCMSA(fp,OUT_CMA[s]); } // exit(1);
	set_typ	*EachSet; NEW(EachSet, Number +3, set_typ);
	for(s=1; s < Hpt->NumSets(); s++){	// skip root and reject sets!!!
		if(RtnSet[s]==0) continue;
		EachSet[s]=CopySet(RtnSet[s][1]); 
		// set_typ tmpSet=CopySet(RtnSet[s][1]); ClearSet(tmpSet);
		for(i=2; RtnSet[s][i]; i++){	// remove descendent sequences...
			IntersectNotSet(EachSet[s],RtnSet[s][i]);
			// UnionSet(EachSet[s],RtnSet[s][i]); 
		}
	}
#if 0	// don't need this for public version...
	fp=tmpfile(); PutMergedCMSA(fp,Number,EachSet,OUT_CMA,0); 
	rewind(fp); rcma=ReadCMSA(fp,AB); fclose(fp);
	// sprintf(str,"%s_core",argv[1]); fp=open_file(str,".cma","w"); PutCMSA(fp,rcma); fclose(fp);
	//============ Sort the sequences the same as for input file ==============
	sprintf(str,"%s_sort",argv[1]); fp=open_file(str,".cma","w"); 
	xcma=SortByCMA(TrueMainCMA,rcma); PutCMSA(fp,xcma); fclose(fp);
	TotalNilCMSA(rcma); TotalNilCMSA(xcma);
#elif 0
	sprintf(str,"%s_core",argv[1]); fp = open_file(str,".cma","w"); 
	PutMergedCMSA(fp,Number,EachSet,OUT_CMA,0); fclose(fp);
#endif
	// Warning: OUT_CMA shares data with IN_CMA!!!
	for(s=1; s < Hpt->NumSets(); s++){ 
	   NilSet(EachSet[s]); NilCMSA(OUT_CMA[s]); TotalNilCMSA(IN_CMA[s]); 
	} free(IN_CMA); free(OUT_CMA); free(EachSet);
#endif

	//===========  7. Output lineage ('line') hierarchical alignment files. ===========
	sprintf(str,"%s_himsa",argv[1]); fp=open_file(str,".sets","w"); // PutCMSA(fp,TmplCMA[1]); 
	for(s=1; s < Hpt->NumSets(); s++){	// skip root and reject sets!!!
	   if(RtnSet[s]==0) continue;
	   fprintf(stderr,"%d.%s: ",s,Hpt->SetName(s));
	   for(i=1; RtnSet[s][i]; i++){ fprintf(stderr,"%d=%d ",i,CardSet(RtnSet[s][i])); }
	   fprintf(stderr,"\n"); WriteSets(fp,i-1,RtnSet[s]);
	   if(RtnSet[s][0]) NilSet(RtnSet[s][0]);	// this appears to be set == 0.
	   for(i=1; RtnSet[s][i]; i++){ if(RtnSet[s][i]) NilSet(RtnSet[s][i]); }
	   free(RtnSet[s]); RtnSet[s]=0;
	} free(RtnSet); fclose(fp); 

	sprintf(str,"%s_himsa",argv[1]); fp=open_file(str,".tpl","w"); // PutCMSA(fp,TmplCMA[1]); 
	hpt_typ *head=0,*hpt,*tail=0; 
	for(p=1; p < Hpt->NumSets(); p++){	// skip root and reject sets!!!
	   // if(TmplCMA[p] == 0) continue;
	   // sprintf(str,"%s_%s",argv[1],Hpt->SetName(p)); fp=open_file(str,".tpl","w"); 
	   if(TmplCMA[p]) PutCMSA(fp,TmplCMA[p]); // fclose(fp); 
	   if(p==1){ head=tail=Hpt; }	// Create lineage hpt files for each non-root node.
	   else { hpt=Hpt->MkLineageTree(p,n); tail->SetNext(hpt); tail = hpt; } 
	   tail->SetFileName(Hpt->SetName(p));
	} sprintf(str,"%s_root",argv[1]); fclose(fp); 

	//============ ?. Print out hierarchies and patterns. =================
	// hpt_typ *ConversionViaTemplateHpt(cma_typ TmplCMA, hpt_typ *Hpt);
	sprintf(str,"%s_himsa",argv[1]); fp=open_file(str,".hpt","w");
	MultiPutHpt(fp,head); fclose(fp);
	// delete Hpt->RtnNext();	Hpt->SetNext(0);  // frees up rest of linked list...

	//========  9. Output the CDD hierarchy file in Newick format. ========
	fp=open_file(str,".ph","w"); Hpt->PutAsTree(fp); fclose(fp);

	//=========== 10. Free up memory. ===============
	for(s=1; s <= Hpt->NumSets(); s++) if(ChildSet[s]) NilSet(ChildSet[s]); 
	free(ChildSet); free(Parent);
    	for(n=1; n <= nsets; n++) NilSet(NodeSet[n]); free(NodeSet);
	for(s=1; s < Hpt->NumSets(); s++){
		if(Oper[s]) free(Oper[s]); 
		if(CsqE[s]) NilSeq(CsqE[s]);
fprintf(stderr,"%d.%s\n",s,Hpt->SetName(s));
	   	if(TmplCMA[s]) TotalNilCMSA(TmplCMA[s]);
		// NodeCMA[1] is also a copy!
		if(0 && s==1){ NilCMSA(NodeCMA[s]); }  // This shares data with MainCMA!
		else if(NodeCMA[s]) TotalNilCMSA(NodeCMA[s]);
	} TotalNilCMSA(TrueMainCMA); NilAlpha(AB);
	delete Hpt; free(NodeCMA); free(TmplCMA); free(CsqE); free(Oper);
	fprintf(stderr, "\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

