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

#include "omc_typ.h"
// #include "blosum62.h"

void	omc_typ::PutSimulatedAln( )
// Generate a random alignment using the TrueMainCMA model.
// hierachy based on a real/optimized hiearchy.
{
	Int4	rpts=1, Gap_Len[4]; Gap_Len[0]=Gap_Len[1]=Gap_Len[2]=0;
	Int4	Stop,i,j,sq,m,n,s,M,N=NumSeqsCMSA(TrueMainCMA),SizeHpt;
	Int4	*Parent,x,k,leng=LengthCMSA(1,TrueMainCMA);
	FILE	*fp,*ofp,*efp=stderr; efp=0;

	// 1. Don't clobber a symbolically linked Main.cma file...other basic input checks.
	sprintf(str,"%s_sim.mma",infile); fp=fopen(str,"r"); 
	sprintf(str,"%s_sim.hsw",infile); ofp=fopen(str,"r"); 
	if(fp || ofp) print_error("FATAL: Delete existing *_sim.mma and *_sim.hsw output files\n");
	assert(nBlksCMSA(TrueMainCMA) == 1);
	assert(Hpt->IsTree(Parent));

	// 2. Create simulated (emitted) sequences for each node.
	mcs->SampleColumns(); mcs->PutHyperPartition(stderr);
	e_type	*tSeq; NEW(tSeq,N+4,e_type);
	set_typ	*sq_set; NEW(sq_set,Hpt->NumSets() + 5,set_typ); 
#if 1	// compute the set size...
	Int4 reject= mcs->SetCard(Hpt->NumSets()) - mcs->GetNumRandom();
	M = N - reject;	// subtracting out rejected sequences.
	Int4 num_random =mcs->ComputeNumRandom(M);	// Not dependent on mcs configuration!
	Int4 size_set=M + num_random + 1;
	set_typ EqSet=MakeSet(size_set);
	fprintf(stderr,"num reject = %d; set size = %d\n",reject,size_set); // exit(1);
#endif
	for(M=0,i=1; i < Hpt->NumSets(); i++){
	   sq_set[i]=MakeSet(size_set);
	   Stop=mcs->SetCard(i);
	   for(n=1; n <= Stop; n++) {
		e_type E=mcs->EmitRandomSeq(i);	// uses weighted residue frequencies...
                sprintf(str,"Set%d_%d simulated_seq",i,n);
                ChangeInfoSeq(str,E);
		M++; tSeq[M]=E; assert(M <= N);
		AddSet(M,sq_set[i]);
	   }
	} sq_set[i]=MakeSet(size_set); ClearSet(sq_set[i]); // == reject (empty).
	// Add random sequences (non-existant) to reject set.
	for(i=Hpt->NumSets(), sq=M+1; sq < size_set; sq++) AddSet(sq,sq_set[i]); 

	// 3. Obtain pattern sets and create sequence and node sets
	sst_typ **xsst=mcs->RtnCopyOfSSTs( );
	set_typ	RootSet,*node_set; NEW(node_set,Hpt->NumSets() + 3, set_typ);
	RootSet= node_set[1]=MakeSet(Hpt->NumSets()+1);
	set_typ	tmp_set=MakeSet(Hpt->NumSets() + 1); 

	// 4. Shuffle residues across each category within each subgroup...
	//    For each position find each equivalence set and shuffle sequences.
	for(k=1; k <=leng; k++){
	   // 4a. Create and initialize the equivalence sets.
	   for(i=2; i < Hpt->NumSets(); i++){
	       if(xsst[i][k] != 0){
		  node_set[i]=MakeSet(Hpt->NumSets()+1);
		  ClearSet(node_set[i]); AddSet(i,node_set[i]); 
	       } else  node_set[i]=0;
	   } ClearSet(node_set[1]); AddSet(1,node_set[1]); 

	   // 4b. Add associated nodes to each equivalence set i.
	   for(i=1; i <= Hpt->NumBPPS(); i++){
	      if(Hpt->TypeOfSet(i) != '!') continue;
	      if(efp) fprintf(stderr,"%d: ",k);
	      ClearSet(tmp_set);
	      for(x=i; x > 0; x=Parent[x]){
		if(xsst[x][k] == 0){	// then put x into tmp_set;
		   AddSet(x,tmp_set); 
		   if(efp) fprintf(stderr,"-"); 
		   if(Parent[x] == 0){ // i.e., the Root node.
			assert(x==1);
			UnionSet(node_set[x],tmp_set);  // node_set[x]= node_set[x] U tmp_set.
		   }
		} else {	// then add the current tmp_set to this node_set 'x'.
		   assert(node_set[x] != 0); 
	           UnionSet(node_set[x],tmp_set); // node_set[x]= node_set[x] U tmp_set.
		   ClearSet(tmp_set); 
		   if(efp) PutSST(stderr,xsst[x][k],AB);
		} if(efp) fprintf(stderr,"[%d] ",x);
	      } 
	      if(efp) fprintf(stderr,"\n");
	   } 

	   // 4c. Shuffle the sequences residues within each equivalence set.
	   for(x=1; x <= Hpt->NumBPPS(); x++){
	     if(node_set[x]==0) continue;	// these don't correspond to pattern positions.
	     fprintf(stderr,"%d: ",k); PutSST(stderr,xsst[x][k],AB); PutSet(stderr,node_set[x]); 

	     // 4c_i. add all sequences to equivalence set 'x'.
	     for(ClearSet(EqSet), i=1; i <= Hpt->NumBPPS(); i++){
		if(MemberSet(i,node_set[x])) UnionSet(EqSet,sq_set[i]); // EqSet=EqSet U sq_set[i].
	     } 
	     // 4c_ii. Shuffle sequences within each equivalence set.
	     dh_type dH=dheap(M+2,4);
	     unsigned char *Res; NEW(Res,M+3,unsigned char);
	     for(i=0,sq=1; sq <= M; sq++) {
		if(MemberSet(sq,EqSet)){
		  i++; Res[i]=ResSeq(k,tSeq[sq]); 
		  insrtHeap(sq,(keytyp) Random(),dH); 
		}
	     } Int4 numRes=i;
	     for(i=0; (sq=delminHeap(dH)) != 0; ){
		i++; assert(i <= numRes);
		EqSeq(k,Res[i],tSeq[sq]); 
	     } assert(i == numRes);
	     Nildheap(dH); free(Res);
	   }

	   // 4d. Ensure that equivalence sets are disjoint and together contain all nodes.
	   ClearSet(tmp_set);
	   for(i=1; i < Hpt->NumSets(); i++){
		if(node_set[i]){  
		   assert(CardInterSet(tmp_set,node_set[i]) == 0);
	           UnionSet(tmp_set,node_set[i]); // tmp_set = tmp_set U node_set[x].
		}
	   } assert(CardSet(tmp_set) == Hpt->NumBPPS());

	   // 4e. Free up sets.
	   for(i=2; i < Hpt->NumSets(); i++){
		if(node_set[i]) NilSet(node_set[i]); node_set[i] = 0;
	   }

	}

	// 5. Create a cma using simulated sequences; swt & hsw freed when leaving scope.
	cma_typ tcma=SimSeqToCMSA(tSeq,M,AB);
	swt_typ swt=swt_typ(tcma,FALSE);
        hsw_typ hsw=swt.RtnHSW( );
	
	// 6. Print out mma and hsw files.
        fp=open_file(infile,"_sim.hsw","w"); FWriteHSW(fp,hsw); fclose(fp);
	fp=open_file(infile,"_sim.mma","w"); PutCMSA(fp,tcma); fclose(fp);

	// 7. Print out diagnostic files for evalBPPS.
	// 8. Create gold standard hiearchy checkpoint files; the sampler should find this.
	fp=open_file(infile,"_simGld_new.mma","w"); 
	for(i=1; i < Hpt->NumSets(); i++){
	   sprintf(str,"Set%d",i); RenameCMSA(str,tcma);
	   if(CardSet(sq_set[i]) > 0){ PutInSetCMSA(fp,sq_set[i],tcma); }
	} fclose(fp);
	// fp=open_file(infile,"_simGld.pttrns","w"); mcs->PutPttrns(fp); fclose(fp);
        // fp=open_file(infile,"_simGld.hsw","w"); FWriteHSW(fp,hsw); fclose(fp);
	// fp=open_file(infile,"_simGld.mma","w"); PutCMSA(fp,tcma); fclose(fp);
	fp=open_file(infile,"_simGld_new.sets","w"); WriteSets(fp,Hpt->NumSets(),sq_set); fclose(fp);

	char    **pttrn; NEWP(pttrn,MaxNumNodesPlus+3,char);
	hpt_typ *hpt= Hpt->Copy();
	for(i=1; i <= Hpt->NumBPPS(); i++){
                char *tmp[3]; tmp[0]=SST2ArgStrAlpha(xsst[i],leng,AB);
		hpt->ReSetArgv(i,1,tmp); free(tmp[0]);
                sprintf(str,"Set%d",i); hpt->ReNameSet(i,str); hpt->ReNameGroup(i,str);
	}
	fp=open_file(infile,"_simGld_new.hpt","w"); hpt->Put(fp); fclose(fp);
	delete hpt;
	
#if 0	// not needed
	fp=open_file(infile,"_simGld.sma","w"); 
	cma_typ *sma; NEW(sma, Hpt->NumSets() + 3, cma_typ);
	for(i=1; i < Hpt->NumSets(); i++){
		sprintf(str,"Set%d",i);
		// cma_typ sma=mcs->RtnCsqAsCMSA(i,str);
		cma_typ sma=mcs->RtnCsqSstAsCMSA(i,str,xsst[i]);
		PutCMSA(fp,sma); TotalNilCMSA(sma);
	} fclose(fp);
#endif

#if 1	// These won't help here!!
	// fp=open_file(infile,"_simGld.seed","w"); fprintf(fp,"-seed=%d\n",RandomSeed); fclose(fp);
	// fp=open_file(infile,"_simGld_bst.out","w"); mcs->PutHyperPartition(fp); fclose(fp);
#endif

	// 6a. 
	mcs->PutHyperPartition(stderr);
#if 0
	// 7. Open up the new alignment as an omc_typ object.
	Int4	argc=0;
	char	*argv[20];
	argv[argc]=AllocString("omcBPPS"); argc++; // name of program.
	sprintf(str,"%s_sim",infile);
	argv[argc]=AllocString(str); argc++; // prefix of input file name.
	argv[argc]=AllocString("-maxcol=20"); argc++; // maxcolumn setting.
	argv[argc]=0;
	omc_typ omc(argc,argv);
	for(i=0; argv[i] != 0; i++) free(argv[i]);
#endif
	// 4. free up memory.
	for(i=1; i <= Hpt->NumBPPS(); i++){ NilSet(sq_set[i]); free(xsst[i]); } 
	NilSet(node_set[1]); free(node_set); free(sq_set); NilSet(tmp_set); free(xsst);
	NilSet(EqSet); TotalNilCMSA(tcma);
}
	
void	omc_typ::PutAbInitioSimulatedAln( )
// Ab initio simulated hierarchy...
{
	print_error("omc_typ::PutAbInitioSimulatedAln() note yet working");
	Int4	rpts=1, Gap_Len[4]; Gap_Len[0]=Gap_Len[1]=Gap_Len[2]=0;
	Int4	Stop,i,j,sq,m,n,s,M,N=NumSeqsCMSA(TrueMainCMA),SizeHpt;
	Int4	*Parent,x,k,leng=LengthCMSA(1,TrueMainCMA);
	assert(nBlksCMSA(TrueMainCMA) == 1);
	double	d = (double)N/250;
	SizeHpt=(Int4) ceil(d);
	SizeHpt=MINIMUM(Int4,SizeHpt,50);

	e_type *Seq=SimulatedSeqsCMSA(TrueMainCMA,N,rpts,Gap_Len);
	cma_typ tcma=SimSeqToCMSA(Seq,N,AB);

	// 0. Create corresponding hsw file; deleted when out of scope
	swt_typ swt = swt_typ(tcma,FALSE);
        hsw_typ mhsw,hsw=swt.RtnHSW( );

	// 1. create a random hierarchy.
	hpt_typ *hpt=Hpt->Randomize(SizeHpt);

	// 2. Create a MainSet.
	n=1+(NumSeqsCMSA(TrueMainCMA)/3); // sam
	cma_typ rcma,mcma;
	mcma=MkMainFileCMSA(tcma,n,rcma);
	mhsw=AddRandomHSW(hsw,tcma,rcma,mcma);
	TotalNilCMSA(rcma); // destroy the temporary Random sequence alignment.

	// 2. Create sequence sets (includes reject set).
	Int4 SetSize=mcs->GetSetSize();
	set_typ	*set; NEW(set,hpt->NumSets()+3,set_typ);
	for(i = 1; i <= hpt->NumSets(); i++) set[i]=MakeSet(SetSize);
	for(i=1; i <= N; i++){
	   s=random_integer(hpt->NumSets()-1); s++; AddSet(i,set[s]);
	}

	// 2. Add sequences to each set .

	// 3. Add a pattern to each level retaining consistency.

	// 4. output the hpt and sma file
	hpt->Put(stderr);

	// 5. output the hierarchical alignment.

#if 0
	rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,hpt,in_sma,Argc,Argv);
	rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,hpt->NumSets(),set,hpt,in_sma,Argc,Argv);
#endif

        FILE* ofp = open_file(infile,"_sim.hsw","w"); FWriteHSW(ofp,hsw); fclose(ofp);
	

	FILE *fp=open_file(infile,"_sim.mma","w"); PutCMSA(fp,tcma); fclose(fp);
	TotalNilCMSA(tcma);
}

