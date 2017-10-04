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

#include "bsa_typ.h"

cma_typ	bsa_typ::GetSeedAln( )
{
	Int4	I,J,setI,setJ,setIJ,NumPhyla;
	char	*phylumI,*phylumJ,str[200],name[200],old_name[200];
	cma_typ cma=input_cma;
	ss_type data = TrueDataCMSA(cma);
	double	score;

	Iter++; 
	// 3. Find best cross-phylum sets of related sequences...
	Int4 NumPrinted=0;
	sprintf(old_name,"%s",NameCMSA(cma)); 
	ReNameCMSA(old_name,cma); 
	while(!EmptyMheap(mH)){
		score=-MinKeyMheap(mH); current_rank++; 
		Int4 x=DelMinMheap(mH); I=StoreI[x]; J = StoreJ[x];
		assert(I <= N && J <= N);
		if(Rank[I] == 0) Rank[I] = current_rank;   // keep track of highest cross-phyla-scoring seqs.
		if(Rank[J] == 0) Rank[J] = current_rank;
		setI = findDSets(I,sets); setJ = findDSets(J,sets);
#if 1	// try reloading heap here...
		if(EmptyMheap(mH)){
		  fprintf(stderr,"************ reloading heap... ************\n");
		  Int4 sz=ReloadHeap();
		  fprintf(stderr,"***** %d sequence pairs added to the heap *****\n",sz);
		}
#endif
		if(setI == setJ) continue; 	// These already belong to the same set...continue.
	   	if(MemberSet(J,TriedSeeds) && MemberSet(I,TriedSeeds)) LinkDisjointSets(setI,setJ);
		else if(MemberSet(J,TriedSeeds)){	// previously printed...
		   setIJ=LinkDisjointSets(setI,setJ); GetDisplaySet(setIJ,data); UnionSet(TriedSeeds,SeedSet); 
		} else if(MemberSet(I,TriedSeeds)){
		   setIJ=LinkDisjointSets(setI,setJ); GetDisplaySet(setIJ,data); UnionSet(TriedSeeds,SeedSet); 
		} else {
		   setIJ=LinkDisjointSets(setI,setJ);	// link the two sets together and return head of seq-list.
		   GetDisplaySet(setIJ,data);		// Sets SeedSet & DisplaySet for setIJ.
		   NumPhyla=CardSet(DisplaySet);
#if 0
fprintf(stderr,"sets: %d + %d; score = %.1f; phyla=%d; d seqs.\n",
	setI,setJ,score,NumPhyla,I-1);
#endif
		   // 3. Print output files upon first occurrance of reaching the limit...
		   if(NumPhyla >= PrintSize){
			NumPrinted++; NumDisplaySets++;
			fprintf(stderr,"%d: setI=%d; setJ=%d; NumPhyla=%d; score=%.1f.\n", 
				NumPrinted,setI,setJ,NumPhyla,score);
		   	fprintf(stderr,"phyla=%d; %d seqs\n",NumPhyla,I-1);
			UnionSet(TriedSeeds,SeedSet); 	// add all sequences in this set to TriedSeeds.

			fprintf(stderr,"making seed sma file\n");
			cma_typ dcma=CreateDisplayCMA(NumPrinted,cma);
			e_type csq=MkConsensusCMSA(dcma);
		        if(!UsableCsq(csq)) {
			  fprintf(stderr,"Similar to previous seed set...skipping\n");
			  NilSeq(csq); NumPrinted--; NumDisplaySets--;  
			  TotalNilCMSA(dcma);
			} else {
			  fprintf(stderr,"====> printing subgroup seed align to sma file.\n");
			  CSQ[NumDisplaySets]=csq;
			  cma_typ cmaQ=AddConsensusCMSA(dcma); TotalNilCMSA(dcma);
			  return cmaQ;
			} 
		   } 
		}	// else already in the same set, do nothing.
	} 	// end of while loop.
#if 0	// DEBUG...
	FILE  *tmp_fp=open_file(infilename,".tried","w");
	PutInSetCMSA(tmp_fp,TriedSeeds,cma); fclose(tmp_fp);
#endif
	return 0;
}

void	bsa_typ::PutCsqCMSA_File(FILE *fp, char *name, cma_typ cma)
{
	a_type AB=AlphabetCMSA(cma);
        e_type  cE=MkConsensusCMSA(cma);
        // fprintf(stderr,"%d: %s\n",I,NameCMSA(IN_CMA[I]));
        fprintf(fp, "[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n", name);
        fprintf(fp,"(%d)",LenSeq(cE));
        for(Int4 i=1; i <= LenSeq(cE); i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",LenSeq(cE),LenSeq(cE));
        fprintf(fp,">%s consensus\n",name);
        char *seq; NEW(seq,LenSeq(cE) +3, char);
        SeqToString(seq, cE, AB);
        fprintf(fp,"{()%s()}*\n\n_0].\n",seq);
        free(seq); NilSeq(cE);
}

double	bsa_typ::fastest_sq2sq_percent_identCMSA(register Int4 length, register unsigned char *isq, register unsigned char *jsq)
{
	register Int4	score=0,s;
	for(s=length; s >= 1; s--,isq++,jsq++){ if(*isq == *jsq) score++; }
	return (double)score*100.0/(double)length;
}

void	bsa_typ::PutChnFileCMSA(FILE *fp, cma_typ cmaQ, char *main_set_file)
{
	Int4 i;
	PutCMSA(fp,cmaQ);
	FILE *ms_fp=open_file(main_set_file,"","r");
        char ms_c;
        while((ms_c=fgetc(ms_fp)) != EOF){
            if(!isprint(ms_c) && !(ms_c == '\n')){
                fprintf(stderr,"non-printable character found: '%c'\n",ms_c);
                print_error("Fatal error!");
            } fprintf(fp,"%c",ms_c);
        } fclose(ms_fp);
	BooLean *skip;
        NEW(skip,NumSeqsCMSA(cmaQ)+3,BooLean);
        for(i=1; i<=NumSeqsCMSA(cmaQ); i++) skip[i]=TRUE;
        for(i=1; i<=NumSeqsCMSA(cmaQ); i++){
            skip[i]=FALSE; PutSelectCMSA(fp,skip,cmaQ); skip[i]=TRUE;
        }
}

void	bsa_typ::PutChnFileCMSA(FILE *fp, cma_typ cmaQ, cma_typ cmaM)
{
	Int4 i;
	PutCMSA(fp,cmaQ);
	PutCMSA(fp,cmaM);
	BooLean *skip;
        NEW(skip,NumSeqsCMSA(cmaQ)+3,BooLean);
        for(i=1; i<=NumSeqsCMSA(cmaQ); i++) skip[i]=TRUE;
        for(i=1; i<=NumSeqsCMSA(cmaQ); i++){
            skip[i]=FALSE; PutSelectCMSA(fp,skip,cmaQ); skip[i]=TRUE;
        }
}

#if 0	// excluding labeled sequences from seed alignments.
	// 1. read in *.mma file.
	// 2. Create a set to keep track of excluded sequences.
	// 3. Can do anything I like with non-seed sequences, as these are not passed on to cmcBPPS
	// 4. 
#endif

#define	USAGE_START	"USAGE: bsa_typ cmafile [options]\n\
   options:\n\
     -s=<int>    - random seed\n\
     -size=<int> - Minium number of phyla in output seed alignments\n\
     -U=<int>    - Precluster sequences at <int> percent identify (default: 95)\n\
     -Lots       - Print lots of output files\n\
     -F=<int>    - Depth of search setting; higher values = search deeper (default: 5)\n\
     -print=<int> - max number of output alignments to report (default: 10)\n\
     -Diversity=<real> - diversity of seed sequences in fraction of identical residues(default: 0.66)\n\
     -c=<cmafile> - use this cmafile as main set for *.chn output file.\n\
     -random     - Return seed alignments at random (default: best scoring first)\n\
     -x          - dummy\n\n"


void    bsa_typ::InitAsNull( )
{
	main_set_file=0;
	time1=time(NULL); 
	Iter = 0;
	NumDisplaySets=0; 
	NEW(CSQ,MaxNumDisplaySets+3,e_type);
}

void    bsa_typ::InitFlags( )
{
	PrintLots=FALSE;
}

void    bsa_typ::InitDefaults( )
{
	MinSeqIdent=95.0;
	SeedDiversity=0.40;
	PrintSize=4,MaxNumPrint=10;
	RandomlyInsert=FALSE;
	Factor=20;
	MaxNumDisplaySets=10000;	// ten thousand as max.
	// Factor=1;	// for testing reloading of the heap...
}

void	bsa_typ::Init(cma_typ in_cma, Int4 argc,char *argv[])
{ 
	Int4	arg;
	FILE	*fp;
	char	*phylumI,*phylumJ,str[200]; 
	UInt4	seed=18364592;

	InitDefaults( ); // WARNING: need to call this before InitAsNull( )!!!
	InitAsNull( );
	InitFlags( );

	if(argc < 2) print_error(USAGE_START);
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
	     case 'D': 
		if(sscanf(argv[arg],"-Diversity=%lf",&SeedDiversity) == 1){
		  if(SeedDiversity > 1.0 || SeedDiversity < 0.01){
			print_error("-Diversity option input error: values out of range.");
		  }
		} else print_error(USAGE_START);
		break;
	     case 'F': Factor=IntOption(argv[arg],'F',1,500,USAGE_START); break;
	     case 'U': MinSeqIdent = (double) IntOption(argv[arg],'U',1,100,USAGE_START);
		// MinSeqIdent=MinSeqIdent/100.0;
		break;
	     case 'L': 
		if(strcmp("-Lots",argv[arg]) == 0) PrintLots=TRUE;
		else print_error(USAGE_START);
		break;
	     case 'p': 
		if(sscanf(argv[arg],"-print=%d",&MaxNumPrint) == 1){
		   if(MaxNumPrint < 1) print_error(USAGE_START);
		} else print_error(USAGE_START);
		break;
	     case 'c': 
		if(argv[arg][2]=='=' && isprint(argv[arg][3])){
                        main_set_file=argv[arg] + 3;
		}
		break;
	     case 'r': 
		if(strcmp("-random",argv[arg]) == 0) RandomlyInsert=TRUE;
		else print_error(USAGE_START);
		break;
	     case 's': 
		if(sscanf(argv[arg],"-size=%d",&PrintSize) == 1){
		   if(PrintSize < 2) print_error(USAGE_START);
		} else if(sscanf(argv[arg],"-s=%d",&seed)!=1)
                                print_error(USAGE_START);
		break;
	     case 'x': break;
	     case ' ': break; // ignore these...
	     default: print_error(USAGE_START);
	   }
	}
        if(seed == 18364592)  seed = (UInt4) time(NULL)/2;
        if(seed != 0) sRandom(seed);
	current_rank=0;
	input_cma=in_cma;
	cma_typ cma = in_cma; AB=AlphabetCMSA(cma);
	TriedSeeds=MakeSet(NumSeqsCMSA(cma)+1); ClearSet(TriedSeeds);	
	SeedSet=MakeSet(NumSeqsCMSA(cma)+1); ClearSet(SeedSet);	
	DisplaySet=MakeSet(NumSeqsCMSA(cma)+1); ClearSet(DisplaySet);	

	ss_type data = TrueDataCMSA(cma);
	ReNumberSeqSet(data);	// renumber colinearly.
	N=NSeqsSeqSet(data);
	Int4	I,J,S;
	e_type	SeqI,SeqJ;
	h_type	HG=0;
	Int4	setI,setJ;
	char name[200],old_name[200];

	//************************** 1. Initialize data structures. **************************
	NEW(StoreI,Factor*N+3,Int4); NEW(StoreJ,Factor*N+3,Int4);
	NEW(SeqList,N+3,e_type); NEW(Rank,N+3,Int4);
	
	for(I=1; I <= N; I++){		// initialize for below...
	   SeqI=SeqSetE(I,data);
	   assert(I==SeqI(SeqI));	// make sure that identifiers are co-linear.
	   SetNextSeq(SeqI,0); 
	   SeqList[I]=SeqI;	// initialize for below...
	   EqSeqI(I,SeqI);	// set sequence identifier
	}
	// fprintf(stderr,"\ttime (init. rnBPPS): %d sec (%0.2f min)\n", time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

	//*********************** 2a. Put seqs from same phyla with 95% identify into one set. ***********************
	//********************** Keep sequences from unknown kingdom or phyla in their own set. **********************
	//=========== 2b. Load the best pairwise scores on a min/max heap. ===========
	HG=Histogram("sequence percent identities", 0,100,5.0);
	sets=DSets(N);   // puts every sequence in its own set.
	mH = Mheap(Factor*N,4);
 {
	// unsigned char **PtrFakeSq,*fsq;	// globally defined...
	unsigned char *fsq;
	assert(N == NumSeqsCMSA(cma) && nBlksCMSA(cma)==1);
	NEWP(PtrFakeSq,N+3,unsigned char);
	Int4	pos[3];
	for(I=1; I <= N; I++){	//************* For speed, store starting positions in each sequence. ***************
		fsq=SeqPtrCMSA(I,cma);	// pointer to fake seq;
		PosSiteCMSA(1,I,pos,cma); fsq+=pos[1]; PtrFakeSq[I]=fsq;	// starting position in seqI.
	}

	ds_type Xsets=DSets(N);   // puts every sequence in its own set.
	double score_cutoff = 100.0 * SeedDiversity;
	for(I=1; I < N; I++){
	   // if(triedseeds && MemberSet(I,triedseeds)) continue;
	   // set_typ triedseeds passed in; but won't create dsets appropriately...
	   SeqI=SeqSetE(I,data); phylumI=PhylumSeq(SeqI); assert(I == SeqI(SeqI));
	   if(KingdomSeq(SeqI) == 'U' || KingdomSeq(SeqI) == 'X') continue;
	   if(phylumI==0) continue;
	   setI = findDSets(I,sets);	// returns the root.
	   Int4 XsetI = findDSets(I,Xsets);	// returns the root.
	   unsigned char   *jsq,*isq = SeqPtrCMSA(I,cma);  // fake seq...
	   Int4 length=LengthCMSA(1,cma);
	   for(J=I+1; J <= N; J++){
	        // if(triedseeds && MemberSet(J,triedseeds)) continue;
	   	SeqJ=SeqSetE(J,data); phylumJ=PhylumSeq(SeqJ);
	   	if(KingdomSeq(SeqJ) == 'U' || KingdomSeq(SeqJ) == 'X') continue;
	   	if(phylumJ==0) continue;
#if 0		// this may not work right with extensions!!!
		double score = fastest_sq2sq_percent_identCMSA(LengthCMSA(1,cma),
				PtrFakeSq[I], PtrFakeSq[J]);
#else
	   	jsq = SeqPtrCMSA(J,cma);  // fake seq...
		Int4	s,si,sj,pos[3],scr;
		PosSiteCMSA(1,I,pos,cma); si=pos[1];
		PosSiteCMSA(1,J,pos,cma); sj=pos[1];
		for(scr=0,s=1; s <= length; s++,si++,sj++){
			if(isq[si] == jsq[sj]) scr++;
		}
		double score = (double)scr*100.0/(double)length;
#endif
		if(strcmp(phylumI,phylumJ) == 0){
		  if(score >= MinSeqIdent){
		     setJ = findDSets(J,sets);  // returns the root. 
		     setI=linkDSets(setI,setJ,sets);
		     Int4 XsetJ = findDSets(J,Xsets);  // returns the root. 
		     XsetI=linkDSets(XsetI,XsetJ,Xsets);
		  }
		} else {	// if from distinct phyla...
		  if(score >= score_cutoff){
		    IncdHist(score,HG);
#if 1		// try inserting seeds at random.
		    Int4 x;
		    if(RandomlyInsert) x=InsertMheap(Random(),mH); else x=InsertMheap(-score,mH);
		    StoreI[x]=I; StoreJ[x]=J;
#else
		    x=InsertMheap(-score,mH);
		    StoreI[x]=I; StoreJ[x]=J;
#endif
		    Int4 XsetJ = findDSets(J,Xsets);  // returns the root. 
		    XsetI=linkDSets(XsetI,XsetJ,Xsets);
		  }
		}
	   }
	   // if(I % 1000 == 0) fprintf(stderr," %d\n",I); // test run times.
	} PutHist(stdout,50,HG); fflush(stdout);  NilHist(HG); // free(PtrFakeSq);
if(0){	// use this to see how many sets are being clustered...
	BooLean	*skip;
	Int4    Size,*RtnSet=RtnOneDSet(Xsets,22,Size);
	NEW(skip,NumSeqsCMSA(cma)+3,BooLean);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++) skip[sq]=TRUE;
	for(Int4 i=1; i <= Size; i++) skip[RtnSet[i]]=FALSE;
	FILE *fp=open_file("test_junk",".cma","w");
	PutSelectCMSA(fp,skip,cma); fclose(fp);
	
	PutDSets(stderr,Xsets); exit(1); // testing..
} NilDSets(Xsets);
   }
	// fprintf(stderr,"\ttime (link same phyla): %d sec (%0.2f min )\n",time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

	//***************  3. link together sequences in the same set... ***************
	sprintf(str,"the sizes of %.1f%c identity sequence clusters",MinSeqIdent,'%');
	HG=Histogram(str, 0,100,3.0);
	e_type	LastSq;
	for(I=1; I <= N; I++){
	     setI = findDSets(I,sets);	// returns the root.
	     // fprintf(stderr,"DEBUG: %d\n",I);
	     if(setI != I) continue;		// keep going until the canonical sequence is found for the set.
	     // else PutDSet(stderr, I,sets);
	     Int4 n=1; LastSq=SeqList[I]; 	// should 
	     assert(NextSeq(LastSq) == 0);
	     for(J=1; J <= N; J++){
	        if(J == I) continue;		// skip self seq.
	        //fprintf(stderr,"DEBUG: %d %d\n",I,J);
		setJ = findDSets(J,sets);  // returns the root. 
	        if(setJ == I){			// find seqJ in the same set.
	            // fprintf(stderr,"DEBUG: %d %d down\n",I,J);
		    n++; SeqJ = SeqList[J]; 
		    // ConcatSeqList(LastSq,SeqJ);
		    SetNextSeq(LastSq,SeqJ);  // 
		    assert(NextSeq(LastSq) == SeqJ); assert(NextSeq(SeqJ) == 0);
		    LastSq=SeqJ;
	            // fprintf(stderr,"DEBUG: %d %d up\n",I,J);
		} 
		assert(NextSeq(LastSq) == 0);
	     } // if(n > 4) fprintf(stderr,"set %d: %d seq\n",I,n);
	     if(n > 2) IncdHist(n,HG);
	} PutHist(stdout,50,HG); NilHist(HG);
	// fprintf(stderr,"\ttime (seq. lists): %d sec (%0.2f min)\n",time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	fflush(stdout);
}

Int4    bsa_typ::GetDisplaySet(Int4 setI, ss_type data)
// SeqI is the head of a list of sequences contained in disjoint setI.
// want to find all sequences in the seed set and the top ranking sequences for the display set. 
{
	Int4	I,J,LenList;
	char	*phylumI,*phylumJ;
	e_type	sq,*Sq,SeqI;

	// 1. Initialize.
	SeqI=SeqList[setI]; ClearSet(DisplaySet); ClearSet(SeedSet);

	// 2. Arrange the sequences within the new set in an array by Rank.
	dh_type   dH=dheap(N,3);	// N = number of seqs in data.
	NEW(Sq,N+3,e_type);
	for(sq=SeqI; sq; sq=NextSeq(sq)){	
		I=SeqI(sq); assert(I <= N);
		AddSet(I,SeedSet); AddSet(I,DisplaySet);  // Add seqs on list to DisplaySet & SeedSet.
		if(Rank[I] == 0) insrtHeap(I,(keytyp) INT4_MAX,dH);	// want to grab the best from each phylum.
		else insrtHeap(I,(keytyp) Rank[I],dH);	// want to grab the best from each phylum.
	}
	for(J=0;(I=delminHeap(dH)) != 0; ){ J++; Sq[J]=SeqSetE(I,data); } LenList=J; assert(LenList > 0);
	Nildheap(dH);

	// 3. Delete from the DisplaySet all but the highest Ranking seq from each phylum.
	for(I=1; I <= LenList; I++){
		phylumI=PhylumSeq(Sq[I]); assert(phylumI);
		for(J=I+1; J <= LenList; J++){
		   phylumJ=PhylumSeq(Sq[J]); assert(phylumJ);
		   if(strcmp(phylumI,phylumJ)==0){		// two sequences from the same phylum
			   DeleteSet(SeqI(Sq[J]),DisplaySet);	// delete these from display set.
		   }
		} 
	} free(Sq);
	return LenList;
}

cma_typ bsa_typ::CreateDisplayCMA(Int4 NumPrinted, cma_typ cma)
{       
        char str[200],name[200],old_name[200];
        FILE *fp=tmpfile();
        if(Iter > 0) sprintf(name,"Set%d_%d",Iter,NumPrinted-1);
        else sprintf(name,"Set%d",NumPrinted-1);		// later change numbering to NumDisplaySets.
        sprintf(old_name,"%s",NameCMSA(cma)); ReNameCMSA(name,cma);
        PutInSetCMSA(fp,DisplaySet,cma); ReNameCMSA(old_name,cma); rewind(fp);
        cma_typ dcma=ReadCMSA(fp,AB); fclose(fp);
        return dcma;
}


Int4	bsa_typ::FindSeedAlns(Int4 iter,char *infilename,Int4 num_random, cma_typ cma)
{
	Int4	I,J,setI,setJ,setIJ,NumPhyla;
	char	*phylumI,*phylumJ,str[200],name[200],old_name[200];
	ss_type data = TrueDataCMSA(cma);
	e_type	SeqI,SeqJ;
	double	score;
	FILE	*fp;

	Iter++; assert(Iter == iter);

	// 3. Find best cross-phylum sets of related sequences...
	Int4 NumPrinted=0;
	sprintf(old_name,"%s",NameCMSA(cma)); 
	sprintf(str,"MainSet%d",Iter); ReNameCMSA(str,cma);
	FILE	*sma_fp=open_file(infilename,".sma","w");
	PutCsqCMSA_File(sma_fp,cma); ReNameCMSA(old_name,cma); fflush(sma_fp);
	while(!EmptyMheap(mH)){
		score=-MinKeyMheap(mH); current_rank++; 
		Int4 x=DelMinMheap(mH); I=StoreI[x]; J = StoreJ[x];
		assert(I <= N && J <= N);
		if(Rank[I] == 0) Rank[I] = current_rank;   // keep track of highest cross-phyla-scoring seqs.
		if(Rank[J] == 0) Rank[J] = current_rank;
		setI = findDSets(I,sets); setJ = findDSets(J,sets);
#if 1	// try reloading heap here...
		if(EmptyMheap(mH)){
		  fprintf(stderr,"************ reloading heap... ************\n");
		  Int4 sz=ReloadHeap();
		  fprintf(stderr,"***** %d sequence pairs added to the heap *****\n",sz);
		}
#endif
		if(setI == setJ) continue; 	// These already belong to the same set...continue.
	   	if(MemberSet(J,TriedSeeds) && MemberSet(I,TriedSeeds)) LinkDisjointSets(setI,setJ);
		else if(MemberSet(J,TriedSeeds)){	// previously printed...
		   setIJ=LinkDisjointSets(setI,setJ); GetDisplaySet(setIJ,data); UnionSet(TriedSeeds,SeedSet); 
		} else if(MemberSet(I,TriedSeeds)){
		   setIJ=LinkDisjointSets(setI,setJ); GetDisplaySet(setIJ,data); UnionSet(TriedSeeds,SeedSet); 
		} else {
		   setIJ=LinkDisjointSets(setI,setJ);	// link the two sets together and return head of seq-list.
		   GetDisplaySet(setIJ,data);		// Sets SeedSet & DisplaySet for setIJ.
		   NumPhyla=CardSet(DisplaySet);
#if 0
fprintf(stderr,"sets: %d + %d; score = %.1f; phyla=%d; d seqs.\n",
	setI,setJ,score,NumPhyla,I-1);
#endif
		   // 3. Print output files upon first occurrance of reaching the limit...
		   if(NumPhyla >= PrintSize){
			NumPrinted++; NumDisplaySets++;
			fprintf(stderr,"%d: setI=%d; setJ=%d; NumPhyla=%d; score=%.1f.\n", 
				NumPrinted,setI,setJ,NumPhyla,score);
		   	fprintf(stderr,"phyla=%d; %d seqs\n",NumPhyla,I-1);
			UnionSet(TriedSeeds,SeedSet); 	// add all sequences in this set to TriedSeeds.

			fprintf(stderr,"making seed sma file\n");
			cma_typ dcma=CreateDisplayCMA(NumPrinted,cma);
			e_type csq=MkConsensusCMSA(dcma);
		        if(!UsableCsq(csq)) {
			  fprintf(stderr,"Similar to previous seed set...skipping\n");
			  NilSeq(csq); NumPrinted--; NumDisplaySets--;  
			} else {
			  fprintf(stderr,"====> printing subgroup seed align to sma file.\n");
			  CSQ[NumDisplaySets]=csq;
			  cma_typ cmaQ=AddConsensusCMSA(dcma);
			  PutCMSA(sma_fp,cmaQ); NilCMSA(cmaQ); 
			} TotalNilCMSA(dcma);
			if(NumPrinted >= MaxNumPrint) break;
		   } 
		}	// else already in the same set, do nothing.
	} 	// end of while loop.
#if 0	// DEBUG...
	FILE  *tmp_fp=open_file(infilename,".tried","w");
	PutInSetCMSA(tmp_fp,TriedSeeds,cma); fclose(tmp_fp);
#endif
	fclose(sma_fp);
	fprintf(stderr,"\ttime (delmheap): %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

	// X. Print out the corresponding Hyperpartition.
    if(NumPrinted > 0){
      FILE *hpt_fp=open_file(infilename,".hpt","w");
      fprintf(hpt_fp,"HyperParTition:\n");
      {
	Int4 S;
	for(S=0; S <= NumPrinted; S++) fprintf(hpt_fp,"!");	// print one extra.
	for(S=0; S <= NumPrinted; S++){
	   fprintf(hpt_fp,"\n+");
	   for(Int4 s=1; s <= NumPrinted; s++){
		if(s == S) fprintf(hpt_fp,"+");
		else fprintf(hpt_fp,"-");
	   } 
	   if(S==0){
		if(Iter > 0) fprintf(hpt_fp," %d.MainSet%d?",S+1,Iter);
		else fprintf(hpt_fp," %d.MainSet?",S+1);
	   } else {
		if(Iter > 0) fprintf(hpt_fp," %d.Set%d_%d!",S+1,Iter,S-1);
		else fprintf(hpt_fp," %d.Set%d!",S+1,S-1);
	   }
	}
	fprintf(hpt_fp,"\n-");
	for(Int4 s=1; s <= NumPrinted; s++) fprintf(hpt_fp,"o");
	fprintf(hpt_fp," %d.Random=%d.\n\n",S+1,num_random);
      } fclose(hpt_fp);
    } return NumPrinted;
}

Int4	bsa_typ::ReloadHeap( )
// Reload the heap if it is empty and return the number of items on the heap.
// If 0 is returned, this means that all of the sequences have been analyzed.
{ 
	char	*phylumI,*phylumJ; 
	ss_type data = TrueDataCMSA(input_cma);
	Int4	I,J,S,Length=LengthCMSA(1,input_cma);
	e_type	SeqI,SeqJ;
	Int4	setI,setJ;

	//=========== Reload the best pairwise scores on a min/max heap. ===========
	h_type	HG=Histogram("sequence percent identities", 0,100,5.0);
	double score_cutoff = 100.0 * SeedDiversity;
	for(I=1; I < N; I++){			// N=NSeqsSeqSet(data); N is global variable...
	   if(MemberSet(I,TriedSeeds)) continue;
	   // set_typ triedseeds passed in; but won't create dsets appropriately...
	   SeqI=SeqSetE(I,data); phylumI=PhylumSeq(SeqI); assert(I == SeqI(SeqI));
	   if(KingdomSeq(SeqI) == 'U' || KingdomSeq(SeqI) == 'X') continue;
	   if(phylumI==0) continue;
	   setI = findDSets(I,sets);	// returns the root.
	   for(J=I+1; J <= N; J++){
	        if(MemberSet(J,TriedSeeds)) continue;
	   	SeqJ=SeqSetE(J,data); phylumJ=PhylumSeq(SeqJ);
	   	if(KingdomSeq(SeqJ) == 'U' || KingdomSeq(SeqJ) == 'X') continue;
	   	if(phylumJ==0) continue;
		setJ = findDSets(J,sets);  // returns the root. 
		if(setI == setJ) continue;	// already in the same set...
		if(strcmp(phylumI,phylumJ) != 0){	// if from distinct phyla...
		  double score = fastest_sq2sq_percent_identCMSA(Length, PtrFakeSq[I], PtrFakeSq[J]);
		  if(score >= score_cutoff){
		    IncdHist(score,HG);
		    Int4 x = InsertMheap(-score,mH); StoreI[x]=I; StoreJ[x]=J;
		  }
		}
	   }
	} PutHist(stdout,50,HG); NilHist(HG); fflush(stdout);
	return ItemsInMheap(mH);
}

void	bsa_typ::Free()
{
	Int4	i,j;
	// X. Free up memory...
	if(sets) NilDSets(sets);
	free(StoreI); free(StoreJ); free(Rank); free(SeqList); free(PtrFakeSq);
	NilMheap(mH);
	for(i = 1; i <= NumDisplaySets; i++) NilSeq(CSQ[i]); free(CSQ);
	if(TriedSeeds) NilSet(TriedSeeds);
	if(DisplaySet) NilSet(DisplaySet);
	if(SeedSet) NilSet(SeedSet);
	// fprintf(stderr,"\ttime (bsa_typ): %d sec (%0.2f min)\n", time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

