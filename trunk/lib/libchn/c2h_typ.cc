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

#include "c2h_typ.h"

Int4	c2h_typ::CheckInput(Int4 N_cma, Int4 N_rows, cma_typ *CMSA,char **Name_row)
{
        Int4    f,r,n=0;
        if(N_cma > N_rows){
                fprintf(stderr,"Number of cmafiles (%d) inconsistent with Hpt rows (%d)\n",N_cma,N_rows);
                exit(1);
        } else if(N_cma <= N_rows){     // Make sure all CMAs and Hpt rows are in the right order.
             for(r=f=1; f <= N_cma; r++,f++){
                if(r > N_rows) {
                        fprintf(stderr,"cmafile names (%d) inconsistent with Hpt (%d)\n",N_cma,N_rows);
                        exit(1);
                }
                if(strcmp(NameCMSA(CMSA[f]),Name_row[r]) !=0){
			n++;
                        fprintf(stderr,"WARNING: %3d: NameRow=%s missing\n",r,Name_row[r]);
                        f--;
                } // else fprintf(stderr,"%3d: NameCMSA=%s; NameRow=%s\n",r,NameCMSA(CMSA[f]),Name_row[r]);
             }
        } return n;

}

Int4    c2h_typ::Intersection(ss_type P1, ss_type P2)
{
        Int4    n,m,Num=0;
        for(n = 1; n <= NSeqsSeqSet(P1); n++){
            e_type E1 = SeqSetE(n,P1);
            BooLean found=FALSE;
            for(m=1; m <= NSeqsSeqSet(P2); m++){
                e_type E2 = SeqSetE(m,P2);
                if(IdentSeqs(E1,E2)){ found=TRUE; break; }
            } if(found) Num++;
        } return Num;
}

BooLean c2h_typ::CompareCMSAs(FILE *fp)
// BooLean	success = CompareCMSAs(NumIRows,irow,NameIRow,ihpt->TypeOfSet(),
//						NumORows,orow,NameORow,ohpt->TypeOfSet());
{
        Int4 f,g,i,o,x,y;

	if(isets[0] && osets[0]){
		x = CardInterSet(isets[0],osets[0]); y = CardSet(isets[0]);
		fprintf(stdout,"\n Input sets share %d/%d sequences in common.\n",x,y);
		// SetN(isets[0]);
	}
        fdt_typ *FDT= new fdt_typ(NumIRows,irow,NameIRow,ihpt->TypeOfSet(),
					NumORows,orow,NameORow,ohpt->TypeOfSet());
        for(o=f=1; f <= NumOn; o++,f++){
           ss_type odata=TrueDataCMSA(ON_CMA[f]);
           // fprintf(stderr,"NameOut=%s; NameRow=%s\n",NameSeqSet(idata),NameIRow[f]);
           if(strcmp(NameCMSA(ON_CMA[f]),NameORow[o]) !=0){ f--; continue; }
           assert(NSeqsSeqSet(odata) > 0);
           // fprintf(stderr,"DEBUG: seqs = %d\n",NSeqsSeqSet(odata));
	   if(isets[0]){
		assert(o <= NumORows);
		x = CardInterSet(isets[0],osets[o]);
           	FDT->AssignHits(o,0,x);
	   } else FDT->AssignHits(o,0,NSeqsSeqSet(odata));
           for(i=g=1; g <= NumIn; i++,g++){
              ss_type idata=TrueDataCMSA(IN_CMA[g]);
              if(strcmp(NameCMSA(IN_CMA[g]),NameIRow[i]) !=0){ g--; continue; }
              Int4 hits=Intersection(odata,idata);
              FDT->AssignHits(o,i,hits);
           }
        } FDT->Put(stdout); delete FDT;

	fprintf(stdout,"\n Inverse comparison:\n");
        FDT= new fdt_typ(NumORows,orow,NameORow,ohpt->TypeOfSet(),
					NumIRows,irow,NameIRow,ihpt->TypeOfSet());
        for(o=f=1; f <= NumIn; o++,f++){
           ss_type odata=TrueDataCMSA(IN_CMA[f]);
           // fprintf(stderr,"NameOut=%s; NameRow=%s\n",NameSeqSet(idata),NameIRow[f]);
           if(strcmp(NameCMSA(IN_CMA[f]),NameIRow[o]) !=0){ f--; continue; }
           assert(NSeqsSeqSet(odata) > 0);
           // fprintf(stderr,"DEBUG: seqs = %d\n",NSeqsSeqSet(odata));
	   if(osets[0]){
	   	x = CardInterSet(osets[0],isets[o]); FDT->AssignHits(o,0,x);
	   } else FDT->AssignHits(o,0,NSeqsSeqSet(odata));
           for(i=g=1; g <= NumOn; i++,g++){
              ss_type idata=TrueDataCMSA(ON_CMA[g]);
              if(strcmp(NameCMSA(ON_CMA[g]),NameORow[i]) !=0){ g--; continue; }
              Int4 hits=Intersection(odata,idata);
              FDT->AssignHits(o,i,hits);
           }
        } FDT->Put(stdout); fprintf(stdout,"\n\n"); delete FDT;
	return TRUE;
}

void	c2h_typ::Free()
{
        Int4 f,s;
        for(f=1; f <= NumIn; f++) TotalNilCMSA(IN_CMA[f]); free(IN_CMA);
        for(f=1; f <= NumOn; f++) TotalNilCMSA(ON_CMA[f]); free(ON_CMA);
	for(s=1; s <= NumOSets; s++){
	   	if(oSets[s]) NilSet(oSets[s]);
		if(osets[s]) NilSet(osets[s]); 
	   	if(misc_osets[s]) NilSet(misc_osets[s]);
	} free(osets); free(oSets); free(misc_osets);
	for(s=1; s <= NumISets; s++){ 
	   	if(iSets[s]) NilSet(iSets[s]);
		if(isets[s]) NilSet(isets[s]); 
	   	if(misc_isets[s]) NilSet(misc_isets[s]);
	} free(isets); free(iSets); free(misc_isets);
	NilSet(leaf_osets); NilSet(leaf_isets); NilSet(tmp_isets); NilSet(tmp_osets);
	NilSet(inter_isets); NilSet(inter_osets);
	delete ohpt; delete ihpt; 
	if(filename1) free(filename1);
	if(filename2) free(filename2);
	free(orow); free(irow); free(NameORow); free(NameIRow);
	free(onfile); free(infile);
	if(TreeI) NilWdgraph(TreeI);
	if(TreeO) NilWdgraph(TreeO); 
	NilAlpha(AB);
}

set_typ	*c2h_typ::GetJackknifeSets(hpt_typ *Hpt,cma_typ *cma,Int4 num_cma)
// Get the sequence sets for Jackknife tests.
{
	assert(MainCMA); 
        Int4	f,s,sq,x,M,N=NumSeqsCMSA(MainCMA);
	set_typ	*rtnSet=0; 
	Int4 nsets=Hpt->NumSets();
	NEW(rtnSet,nsets+3,set_typ);
	for(s=1; s <= nsets; s++){ rtnSet[s]=MakeSet(N+1); ClearSet(rtnSet[s]); }
	for(f=s=1; s <= nsets; f++,s++){
	   char *NameRow=Hpt->ElmntSetName(s);
           if(strcmp(NameCMSA(cma[f]),NameRow) !=0){ f--; continue; }
	   M=NumSeqsCMSA(cma[f]);
           for(x = 1; x <= M; x++){
	   	e_type Ex=TrueSeqCMSA(x,cma[f]);
           	BooLean found=FALSE;
           	for(sq = 1; sq <= N; sq++){
                   if(IdentSeqs(Ex,TrueSeqCMSA(sq,MainCMA))){
			AddSet(sq,rtnSet[s]); found=TRUE; break; 
		   }
		} assert(found);
	   }
	} f--; assert(f==num_cma); // make sure that all sets were found...
	// make Union of sets for each.
	rtnSet[0]=MakeSet(N+1); ClearSet(rtnSet[0]);
	for(s=1; s <= nsets; s++) UnionSet(rtnSet[0],rtnSet[s]); // Set1 = Set1 U Set2.
	return rtnSet;
}

const char USAGE00[]="usage: evalBPPS run1_prefix run2_prefix [options]\n\
    NOTE: Input hyperpartition files = query.hpt & dbsfile.hpt\n\
            Input set files = query.sets & dbsfile.sets\n\
	    Compares _new.mma output files\n\
    Options:\n\
        -M=<name> Main CMA file (needed for Jackknife test).\n\
        -file1=<name> first file identifier.\n\
        -file2=<name> second file identifier.\n\
    Reference:\n\
	Neuwald, A.F. 2014. Evaluating, comparing and interpreting protein domain hierarchies.\n\
	Journal of Computational Biology 21(4): 287-302.\n\n";

const char USAGE11[]="usage: evalBPPS <prefix1> <prefix2> [options]\n\
    Input: <prefix1>_new.hpt & <prefix2>_new.hpt\n\
           <prefix1>_new.mma & <prefix2>_new.mma\n\
           <prefix1>.chk & <prefix2>.chk\n\
    Output: stdout: primary analysis results\n\
	   <prefix1>_<labe2>_trees.rtf: file for viewing hierachies in MS Excel\n\
    Reference:\n\
	Neuwald, A.F. 2014. Evaluating, comparing and interpreting protein domain hierarchies.\n\
	Journal of Computational Biology 21(4): 287-302.\n\n";

const char USAGE[]="Usage: evalBPPS <prefix1> <prefix2> [options]\n\
    Input: <prefix1>_new.hpt & <prefix2>_new.hpt\n\
           <prefix1>_new.mma & <prefix2>_new.mma\n\
           <prefix1>.chk & <prefix2>.chk\n\
    Output: stdout: primary analysis results\n\
    Options:\n\
        -trees  Output <prefix1>_<prefix2>_trees.rtf file for viewing hierachies in MS Excel\n\
    Reference:\n\
	Neuwald, A.F. 2014. Evaluating, comparing and interpreting protein domain hierarchies.\n\
	Journal of Computational Biology 21(4): 287-302.\n\n";

const char PUBLIC_USAGE[]="<prefix1> <prefix2> [options]\n\
    Input: <prefix1>.mma & <prefix2>.mma: identical input MSAs in cma format.\n\
           <prefix1>.chk & <prefix2>.chk: BPPS 1 checkpoint files.\n\
    Output: stdout: comparisons between hierarchies.\n\
    Options:\n\
        -trees  Output <prefix1>_<prefix2>_trees.rtf file for viewing hierachies in MS Excel\n\
    Reference:\n\
	Neuwald, A.F. 2014. Evaluating, comparing and interpreting protein domain hierarchies.\n\
	Journal of Computational Biology 21(4): 287-302.\n\n";

void	c2h_typ::PrintError(char code,char *program_name)
{
	if(code=='E'){ fprintf(stderr,"Usage: %s E ",program_name); print_error(PUBLIC_USAGE); }
	else print_error(USAGE);
}

void	c2h_typ::Init(char code, int argc, char *argv[])
// this is the same routine as the public version, except for the usage...
{
	FILE		*fp,*tfp;
	char		str[305],mode=' ';
	Int4		c,i,j,o,s,r,arg;
	MainCMA=0;
	AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	filename1=0; filename2=0;
	TreeI=TreeO=0;

	iHpt=oHpt=0;
	NodeMap=0; BackMap=0;
	time1=time(NULL);
	if(argc < 3) PrintError(code,argv[0]);
	infile=argv[2]; onfile=argv[1];
	print_tree=FALSE;

	for(arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') PrintError(code,argv[0]);
           switch(argv[arg][1]) {
             case 'f': 
		if(sscanf(argv[arg],"-file1=%s",str) == 1){
			if(filename1) free(filename1);
			filename1=AllocString(str);
		} else if(sscanf(argv[arg],"-file2=%s",str) == 1){
			if(filename2) free(filename2);
			filename2=AllocString(str);
		} else PrintError(code,argv[0]);
		break;
             case 't': 
		if(strcmp("-trees",argv[arg])==0) print_tree=TRUE;
		break;
             case 'M': 
		if(sscanf(argv[arg],"-M=%s",str) != 1) PrintError(code,argv[0]);
		MainCMA=ReadCMSA2(str,AB);
		break;
             case 'g': mode='g'; break;
             default: PrintError(code,argv[0]);
           }
        }
	buffer[50000]=0;	// set for overflow.
        fprintf(stderr,"\n%s vs %s:\n",onfile,infile);
	//============================== *.chk Input ================================
	BooLean okay=TRUE;
	sprintf(str,"%s.chk",infile); fp=fopen(str,"r");
	if(fp == NULL){ okay=FALSE;  if(code=='E') PrintError(code,argv[0]); }
	else {
	   fclose(fp); sprintf(str,"%s.chk",onfile); fp=fopen(str,"r");
	   if(fp == NULL){ okay=FALSE; PrintError(code,argv[0]); } else fclose(fp);
	} 
	if(okay){
	  fp=open_file(infile,".chk","r");
          {
	    bpcp_typ bpcp; bpcp.Read(fp,AB); fclose(fp); isets=bpcp.RtnSets(NumISets); ihpt=bpcp.RtnHpt();
	    cma_typ *cma=bpcp.RtnSMA(i); for(j=1; j <= i; j++) TotalNilCMSA(cma[j]); free(cma);
	  }
	  fp=open_file(onfile,".chk","r");
          {
	    bpcp_typ bpcp; bpcp.Read(fp,AB); fclose(fp); osets=bpcp.RtnSets(NumOSets); ohpt=bpcp.RtnHpt();
	    cma_typ *cma=bpcp.RtnSMA(i); for(j=1; j <= i; j++) TotalNilCMSA(cma[j]); free(cma);
	  }

	  cma_typ  in_cma,on_cma;
	  fp=open_file(infile,".mma","r"); in_cma=ReadCMSA(fp,AB); fclose(fp);
	  fp=open_file(onfile,".mma","r"); on_cma=ReadCMSA(fp,AB); fclose(fp);

	  tfp=tmpfile();
	  set_typ Rset=CopySet(isets[ihpt->NumSets()]); // random sequence set;
	  for(s=SetN(Rset),s--; s > NumSeqsCMSA(in_cma); s--) DeleteSet(s,Rset); 
          for(s=1; s < ihpt->NumSets(); s++) if(CardSet(isets[s]) > 0){
		ReNameCMSA(ihpt->SetName(s),in_cma); PutInSetCMSA(tfp,isets[s],in_cma); 
	  }
          if(CardSet(Rset) > 0){
		ReNameCMSA(ihpt->SetName(s),in_cma); PutInSetCMSA(tfp,Rset,in_cma); 
	  } rewind(tfp);
          IN_CMA=MultiReadCMSA(tfp,&NumIn,AB); fclose(tfp);

	  tfp=tmpfile();
	  Rset=CopySet(isets[ohpt->NumSets()]); // random sequence set;
	  for(s=SetN(Rset),s--; s > NumSeqsCMSA(on_cma); s--) DeleteSet(s,Rset); 
          for(s=1; s < ohpt->NumSets(); s++) if(CardSet(osets[s]) > 0){
		ReNameCMSA(ohpt->SetName(s),on_cma); PutInSetCMSA(tfp,osets[s],on_cma); 
	  }
          if(CardSet(Rset) > 0){
		ReNameCMSA(ohpt->SetName(s),on_cma); PutInSetCMSA(tfp,Rset,on_cma); 
	  } rewind(tfp);
          ON_CMA=MultiReadCMSA(tfp,&NumOn,AB); fclose(tfp);
	  TotalNilCMSA(in_cma); TotalNilCMSA(on_cma);
	} else {  //=========================== Old Input ================================
	  sprintf(str,"%s_new",infile); ihpt = new hpt_typ(str);
	  sprintf(str,"%s_new",onfile); ohpt = new hpt_typ(str);
          fp=open_file(infile,"_new.mma","r");
          IN_CMA=MultiReadCMSA(fp,&NumIn,AB); fclose(fp);
          fp=open_file(onfile,"_new.mma","r");
          ON_CMA=MultiReadCMSA(fp,&NumOn,AB); fclose(fp);
	  if(MainCMA){
	    isets=GetJackknifeSets(ihpt,IN_CMA,NumIn); NumISets=ihpt->NumSets();
	    osets=GetJackknifeSets(ohpt,ON_CMA,NumOn); NumOSets=ohpt->NumSets();
	  } else {
	    sprintf(str,"%s_chk.sets",infile);
	    fp=fopen(str,"r");
	    if(fp == NULL) fp=open_file(infile,"_new.sets","r"); 
	    isets=ReadSets(fp,NumISets); fclose(fp);
	    sprintf(str,"%s_chk.sets",onfile);
	    fp=fopen(str,"r");
	    if(fp == NULL) fp=open_file(onfile,"_new.sets","r"); 
	    osets=ReadSets(fp,NumOSets); fclose(fp);
	  }
	} //===========================================================================

	// fprintf(stderr,"NumOSets = %d; ohpt.NumSets() = %d\n",NumOSets,ohpt->NumSets());
	// fprintf(stderr,"NumISets = %d; ihpt.NumSets() = %d\n",NumISets,ihpt->NumSets());
	assert(NumOSets == ohpt->NumSets()); 
	assert(NumISets == ihpt->NumSets());

	// ihpt.Put(stderr); ohpt.Put(stderr);
	NumIRows=ihpt->NumSets();
	NumORows=ohpt->NumSets();
	NEWP(orow,NumORows+3,char); NEWP(irow,NumIRows+3,char);
	NEWP(NameORow,NumORows+3,char); NEWP(NameIRow,NumIRows+3,char);
	for(r=1; r <= ihpt->NumSets(); r++){
		irow[r]=ihpt->RtnHyperPartition(r);
		NameIRow[r]=ihpt->ElmntSetName(r);
	}
LenEqOI[0]=ohpt->NumSets(); LenEqOI[1]=ihpt->NumSets();
NEWP(EqualOI,ohpt->NumSets() +3, char);
	for(r=1; r <= ohpt->NumSets(); r++){
NEW(EqualOI[r],ihpt->NumSets() +3, char);
		orow[r]=ohpt->RtnHyperPartition(r);
		NameORow[r]=ohpt->ElmntSetName(r);
	}

        Mempty=CheckInput(NumIn,NumIRows,IN_CMA,NameIRow);

        Nempty=CheckInput(NumOn,NumORows,ON_CMA,NameORow);

	// fprintf(stderr,"NumORows = %d; NumOn = %d; NumOSets = %d\n", NumORows,NumOn,NumOSets);
	// fprintf(stderr,"NumIRows = %d; NumIn = %d; NumISets = %d\n", NumIRows,NumIn,NumISets);

	Int4 size=0,x,N;
	for(o=0,s=1; s <= NumOSets; s++){
		if(osets[s]){
		    if(size == 0) size = SetN(osets[s]);
		    else assert(SetN(osets[s]) == size);
		    x = CardSet(osets[s]);
		    o++; // assert(x == NumSeqsCMSA(ON_CMA[o]));
           	    if(o > NumOn || strcmp(NameCMSA(ON_CMA[o]),NameORow[s]) !=0){ 
			o--; // fprintf(stderr,"%d.%s (%d); [missing]\n",s,NameORow[s],x);
		    } else {
			N=NumSeqsCMSA(ON_CMA[o]);
			// fprintf(stderr,"%d.%s (%d); [%d]\n",s,NameORow[s],x,N); 
			assert(s==NumOSets || x == N);
		    }
		} // else fprintf(stderr,"Set %d 0\n",s);

	} // fprintf(stderr,"\n\n");
	for(i=0,s=1; s <= NumISets; s++){
		if(isets[s]){
		    assert(SetN(isets[s]) == size);
		    x = CardSet(isets[s]); i++;
           	    if(i > NumIn || strcmp(NameCMSA(IN_CMA[i]),NameIRow[s]) !=0){ 
			i--; // fprintf(stderr,"%d.%s (%d); [missing]\n",s,NameIRow[s],x);
		    } else {
			N=NumSeqsCMSA(IN_CMA[i]);
			// fprintf(stderr,"%d.%s (%d); [%d]\n",s,NameIRow[s],x,N); 
			assert(s==NumISets || x == N);
		    }
		} // else fprintf(stderr,"Set %d 0\n",s);
	} // fprintf(stderr,"\n\n");
	// fprintf(stderr,"SetSizeO = %d; SetSizeI = %d\n\n", size,size);
	// make misc sets...
	leaf_isets=MakeSet(NumISets+1); ClearSet(leaf_isets);
	tmp_isets=MakeSet(NumISets+1); ClearSet(tmp_isets);
	inter_isets=MakeSet(NumISets+1); ClearSet(inter_isets);

	leaf_osets=MakeSet(NumOSets+1); ClearSet(leaf_osets);
	tmp_osets=MakeSet(NumOSets+1); ClearSet(tmp_osets);
	inter_osets=MakeSet(NumOSets+1); ClearSet(inter_osets);

	NEW(misc_osets,NumOSets+3,set_typ);
	NEW(oSets,NumOSets+3,set_typ);
	// fprintf(stderr,"\n\n=============== Intermediate Sets =============\n");
	for(s=2; s < NumOSets; s++){
	   if(!osets[s]) continue;
	   // oSets[s]=CopySet(osets[s]); ClearSet(oSets[s]);
	   oSets[s]=CopySet(osets[NumOSets]); // include rejected nodes.
	   for(c=1; c < s; c++){ 	// include root node: c == 1.
		if(orow[s][c] == '+'){
		    if(osets[c]) UnionSet(oSets[s],osets[c]); 
		    // WARNING: Assumes that the input corresponds to a tree!
		}
	   }
	   if(ohpt->TypeOfSet(s) != '?'){ AddSet(s,leaf_osets); continue; }
	   else AddSet(s,inter_osets);
	   misc_osets[s] = CopySet(osets[s]);
	   for(r=1; r <= NumORows; r++){
		if(r == s) continue;
		if(!osets[r]) continue;
		if(orow[r][s] == '+'){
			UnionSet(misc_osets[s], osets[r]);	// modifies first set
		}
	   } // fprintf(stderr,"%d.%s (%d)\n",s,NameORow[s],CardSet(misc_osets[s])); 
	}
	// fprintf(stderr,"\n\n=============== Intermediate Sets =============\n");
	NEW(misc_isets,NumISets+3,set_typ);
	NEW(iSets,NumISets+3,set_typ);
	for(s=2; s < NumISets; s++){
	   if(!isets[s]) continue;
	   // iSets[s]=CopySet(isets[s]); ClearSet(iSets[s]);
	   iSets[s]=CopySet(isets[NumISets]); 	// include rejected nodes.
	   for(c=1; c < s; c++){ 
		if(irow[s][c] == '+'){   // row 's' column 'c'.
		   if(isets[c]) UnionSet(iSets[s],isets[c]); 
		   // WARNING: Assumes that the input corresponds to a tree!
		}
	   } 
	   // if(ihpt->TypeOfSet(s) != '?') continue;
	   if(ihpt->TypeOfSet(s) != '?'){ AddSet(s,leaf_isets); continue; }
	   else AddSet(s,inter_isets);
	   misc_isets[s] = CopySet(isets[s]);
	   for(r=1; r <= NumIRows; r++){
		if(r == s) continue;
		if(!isets[r]) continue;
		if(irow[r][s] == '+'){    // row 'r', column 's'.
			UnionSet(misc_isets[s], isets[r]);	// modifies misc_isets[s].
		}
	   } // fprintf(stderr,"%d.%s (%d).\n",s,NameIRow[s],CardSet(misc_isets[s])); 
	} // fprintf(stderr,"\n\n");
	char *ptr=argv[1]; for(i=0; ptr[i] !=0; i++) ;
	while(i >= 0 && ptr[i] != '/') i--; if(i > 0) ptr += i+1; 
	onfile=AllocString(ptr);

	ptr=argv[2]; for(i=0; ptr[i] !=0; i++) ;
	while(i >= 0 && ptr[i] != '/') i--; if(i > 0) ptr += i+1;
	infile=AllocString(ptr);

}

Int4	c2h_typ::IdenticalFuzzySets(Int4 Num[2], set_typ *Sets[2],set_typ *PSet[2], char **NameRow[2],hpt_typ *hpt,
			set_typ &RtnSet)
//PSet[2] is the parental leafover set for Sets[2].
// Can eventually combine this with IdenticalSets() routine...
{
	Int4	i,j,o,s,r,arg,NumSubTrees=0;
	set_typ	Set0o=0,Set1i=0;
	int	bos=0; bptr = &buffer[0]; bptr[0]=0;	// buffer pointer and offset;
	double	*FO,*FI; NEW(FO,Num[0]+3,double); NEW(FI,Num[1]+3,double);
	
	set_typ match_sets=MakeSet(Num[0]+1);
	for(o=2; o < Num[0]; o++){
	   if(!Sets[0][o]) continue;
	   if(Set0o==0) Set0o=CopySet(Sets[0][o]); 
	   Int4 co = CardSet(Sets[0][o]);
	   for(i=2; i < Num[1]; i++){
	   	if(!Sets[1][i]) continue;
	   	if(Set1i==0) Set1i=CopySet(Sets[1][i]); 

		Int4 ci = CardSet(Sets[1][i]);
		Int4 inter = CardInterSet(Sets[0][o],Sets[1][i]);
		// Int4 min=MINIMUM(Int4,co,ci);
		double Fo,fo=(double)inter/(double)co;
		double Fi,fi=(double)inter/(double)ci;
		if(fo >= 0.90 && fi >= 0.90) continue; // found previously
		else if(fo >= 0.50 && fi >= 0.50){
		   // then check for overlap between core sets...
		   // ignoring sequences placed into opposite parent nodes.
		   assert(PSet[1][i]); assert(PSet[0][o]);
		   IntersectNotSet(Sets[0][o],PSet[1][i],Set0o);
		   IntersectNotSet(Sets[1][i],PSet[0][o],Set1i);
		   // IntersectNotSet(set_typ B1, set_typ B2, set_typ notIB);
		   // modifies notIB to equal B1 intersect not B2
		   // UnionSet3(Set1i,Sets[0][o],PSet[1][i]);  //  
		   Int4 Co = CardSet(Set0o),Ci = CardSet(Set1i);
		   inter = CardInterSet(Set0o,Set1i);
		   Fo=(double)inter/(double)Co;
		   Fi=(double)inter/(double)Ci;
#if 0
		   Int4 Union = CardUnionSet(Set0o,Set1i);
		   // Int4 Union = CardUnionSet(Sets[0][o],Sets[1][i]);
		   double D= (double)inter/(double)Union;
#endif
#if 0
		   fprintf(stdout,"%d.%s (%d) == %d.%s (%d)\n",
					o,NameRow[0][o],Co,i,NameRow[1][i],Ci);
		   fprintf(stdout,"Fo=%g; Fi=%g; inter=%d\n",Fo,Fi,inter);
#endif
		   if(Fo >= 0.90 && Fi >= 0.90){
#if 0	// Debug...
if(FO[o] > 0.0) fprintf(stderr,"%d.%s %.3f vs %.3f\n",o,HPT[0]->ElmntSetName(o),D,FO[o]);
if(FI[i] > 0.0) fprintf(stderr,"%d.%s %.3f vs %.3f\n\n",i,HPT[1]->ElmntSetName(i),D,FI[i]);
#endif
#if 0
			if(D < FO[o] || D < FI[i]) continue;  // already found a better match.
			else { FO[o]=D; FI[i]=D; }
#endif
			sprintf(tmpstr,"~%s",HPT[1]->ElmntSetName(i)); oHpt->ReNameSet(o,tmpstr);
			sprintf(tmpstr,"~%s",HPT[0]->ElmntSetName(o)); iHpt->ReNameSet(i,tmpstr);
			bos=sprintf(bptr,"%d.%s (%d) == %d.%s (%d)\n",
					o,NameRow[0][o],co,i,NameRow[1][i],ci);
#if 1	// consensus hierarchy..
// double  dd=MAXIMUM(double,fo,fi);
double  dd=MINIMUM(double,fo,fi);
char	chr=(char) floor(100*dd);
Int4	ii,oo;
if(NameRow[0]==NameORow){
    if(EqualOI[o][0] == 0){ EqualOI[o][i]=chr; EqualOI[o][0] = chr; }
    else if(chr > EqualOI[o][0]){
for(ii = 1; ii <= LenEqOI[1]; ii++) if(EqualOI[o][ii] > 0) break;
fprintf(stderr,"~i? %d(%d).%s = %d.%s ",o,BackMap[o],NameRow[0][o],ii,NameRow[1][ii]);
fprintf(stderr,"(%d) replaced by %d.%s (%d)\n",EqualOI[o][0],i,NameRow[1][i],chr);
// fprintf(stderr,"nodes %d = %d (%d) replaced by %d (%d)\n",o,ii,EqualOI[o][0],i,chr);
	EqualOI[o][0]=chr; 
	for(Int4 ii = 1; ii <= LenEqOI[1]; ii++) EqualOI[o][ii]=0;
	EqualOI[o][i]=chr; EqualOI[o][0] = chr; 
    } 
} else {	// i[1] is really o[0] and vice versa
    assert(NameRow[0]==NameIRow);
    if(EqualOI[i][0] == 0){ EqualOI[i][o]=chr; EqualOI[i][0] = chr; }
    else if(chr > EqualOI[i][0]){
for(oo = 1; oo <= LenEqOI[0]; oo++) if(EqualOI[i][oo] > 0) break;
fprintf(stderr,"~o? %d.%s = %d.%s ",i,NameRow[1][i],oo,NameRow[0][oo]);
fprintf(stderr,"(%d) replaced by %d.%s (%d)\n",EqualOI[i][0],o,NameRow[0][o],chr);
// fprintf(stderr,"nodes %d = %d (%d) replaced by %d (%d)\n",i,oo,EqualOI[i][0],o,chr);
	EqualOI[i][0]=chr; 
	for(oo = 1; oo <= LenEqOI[0]; oo++) EqualOI[i][oo]=0;
	EqualOI[i][o]=chr; EqualOI[i][0] = chr; 
    } 
}
#else
if(NameRow[0]==NameORow){
	EqualOI[o][i]=1;
} else {
	EqualOI[i][o]=1;
}
#endif
			bptr += bos;
			if(hpt){	// condensed = subtree to leaf.
			   for(r=1; r <= Num[0]; r++){
				if(r==o) continue; // column o and row r == '+'.
				if(hpt->RtnHyperPartition(r,o) == '+'){
				    if(Sets[0][r]){
					bos=sprintf(bptr," child: %d.%s (%d).\n",
					      r,NameRow[0][r],CardSet(Sets[0][r]));
				    } else {
					bos=sprintf(bptr," child: %d.%s.\n",
					      r,NameRow[0][r]);
				    } bptr += bos; AddSet(r,match_sets);
				}
			   } 
		        }
			AddSet(o,match_sets);  NumSubTrees++;
		   }
		}
	   }
	} if(Set0o) NilSet(Set0o); if(Set1i) NilSet(Set1i);
	free(FO); free(FI);
	RtnSet=match_sets;
	return NumSubTrees;
}

set_typ c2h_typ::SplitSets(Int4 Num[2], set_typ *Sets[2],char **NameRow[2],hpt_typ *hpt[2],set_typ &RtnISet)
// Sets[1] is checked set, not Sets[0].
{
        Int4    i,j,o,s,r,arg;
	int	bos=0; bptr = &buffer[0]; bptr[0]=0; // buffer pointer and offset;

        set_typ split_osets=MakeSet(Num[0]+1);
        set_typ tmpset=0;
	set_typ	tmpiset=MakeSet(Num[1]+1),rtniset=MakeSet(Num[1]+1);
	ClearSet(rtniset); ClearSet(split_osets);

        for(o=1; o < Num[0]; o++){
           if(!Sets[0][o]) continue;
	   if(hpt[0]->TypeOfSet(o) == '?') continue;
           if(tmpset==0){ tmpset=CopySet(Sets[0][o]); ClearSet(tmpset); }
	   else ClearSet(tmpset); ClearSet(tmpiset);
           Int4 co = CardSet(Sets[0][o]);
           for(i=1; i < Num[1]; i++){
                if(!Sets[1][i]) continue;
	        if(hpt[1]->TypeOfSet(i) == '?') continue;
                Int4 ci = CardSet(Sets[1][i]);
                Int4 inter = CardInterSet(Sets[0][o],Sets[1][i]);
                double fo=(double)inter/(double)co;
                double fi=(double)inter/(double)ci;
                if(fo >= 0.90 && fi >= 0.90){ continue;	// not checked here...
                } else if(fi >= 0.90){	// o set = to 2 or more i sets?
		   UnionSet(tmpset,Sets[1][i]); AddSet(i,tmpiset);
                   ci = CardSet(tmpset);
                   inter = CardInterSet(Sets[0][o],tmpset);
                   fo=(double)inter/(double)co;
                   fi=(double)inter/(double)ci;
                   if(fo >= 0.90 && fi >= 0.90){ // found a split match...
                        AddSet(o,split_osets);
                        UnionSet(rtniset,tmpiset);

                        bos=sprintf(bptr,"%d.%s (%d) == ",o,NameRow[0][o],co); bptr += bos;
			for(j=1; j < Num[1]; j++){
			    if(MemberSet(j,tmpiset)){
#if 0	// need to work on this... need an example..
			sprintf(tmpstr,"s_%s",HPT[1]->ElmntSetName(i)); oHpt->ReNameSet(o,tmpstr);
			sprintf(tmpstr,"s_%s",HPT[0]->ElmntSetName(o)); iHpt->ReNameSet(i,tmpstr);
#endif
                               bos=sprintf(bptr,"%d.%s (%d)  ",j,NameRow[1][j],CardSet(Sets[1][j]));
			       bptr += bos;
			    }
			} bos=sprintf(bptr,"\n"); bptr += bos;
			// break;  // look no further?
		   }
                }
           }
        }
	if(tmpset) NilSet(tmpset); NilSet(tmpiset);
	RtnISet=rtniset;
	return split_osets;
}

Int4	c2h_typ::IdenticalSets(Int4 Num[2], set_typ *Sets[2],char **NameRow[2],hpt_typ *hpt,set_typ &RtnSet)
{
	Int4	i,j,o,s,r,arg;
	set_typ	match_sets=MakeSet(Num[0]+1);
	int	bos=0; bptr = &buffer[0]; bptr[0]=0; // buffer pointer and offset;
	Int4	NumSubTrees=0;
	double	*FO,*FI; NEW(FO,Num[0]+3,double); NEW(FI,Num[1]+3,double);
	// store scores to see whether they should be replaced...
	// set_typ	tmpiset=0,tmposet=0;
	for(o=2; o < Num[0]; o++){
	   if(!Sets[0][o]) continue;
	   // if(tmpset=0){ tmpset=CopySet(Sets[0][o]); ClearSet(tmpset); }
	   Int4 co = CardSet(Sets[0][o]);
	   for(i=2; i < Num[1]; i++){
	   	if(!Sets[1][i]) continue;
		Int4 ci = CardSet(Sets[1][i]);
		Int4 inter = CardInterSet(Sets[0][o],Sets[1][i]);
		// Int4 min=MINIMUM(Int4,co,ci);
		double fo=(double)inter/(double)co;
		double fi=(double)inter/(double)ci;
		Int4 Union = CardUnionSet(Sets[0][o],Sets[1][i]);
		double D= (double)inter/(double)Union;
		if(fo >= 0.90 && fi >= 0.90){
#if 0	// Debug...
if(FO[o] > 0.0 || o >= 41 && o <= 43) fprintf(stderr,"%d.%s %.3f vs %.3f\n",o,HPT[0]->ElmntSetName(o),D,FO[o]);
if(FI[i] > 0.0 || i >=41 && i <= 43) fprintf(stderr,"%d.%s %.3f vs %.3f\n\n",i,HPT[1]->ElmntSetName(i),D,FI[i]);
#endif
			if(D < FO[o] || D < FI[i]) continue;  // already found a better match.
			else { FO[o]=D; FI[i]=D; }
			sprintf(tmpstr,"x_%s",HPT[1]->ElmntSetName(i)); oHpt->ReNameSet(o,tmpstr);
			sprintf(tmpstr,"x_%s",HPT[0]->ElmntSetName(o)); iHpt->ReNameSet(i,tmpstr);
			bos=sprintf(bptr,"%d.%s (%d) == %d.%s (%d)\n",
					o,NameRow[0][o],co,i,NameRow[1][i],ci);
#if 1	// consensus hierarchy..
double  dd=MAXIMUM(double,fo,fi);
char	chr=(char) floor(100*dd);
Int4	ii,oo,OO,II;
if(NameRow[0]==NameORow){
    if(EqualOI[o][0] == 0){ EqualOI[o][i]=chr; EqualOI[o][0] = chr; }
    else if(chr > EqualOI[o][0]){
for(ii = 1; ii <= LenEqOI[1]; ii++) if(EqualOI[o][ii] > 0) break;;
fprintf(stderr,"i? %d(%d).%s = %d.%s ",o,BackMap[o],NameRow[0][o],ii,NameRow[1][ii]);
fprintf(stderr,"(%d) replaced by %d.%s (%d)\n",EqualOI[o][0],i,NameRow[1][i],chr);
	EqualOI[o][0]=chr; 
	for(Int4 ii = 1; ii <= LenEqOI[1]; ii++) EqualOI[o][ii]=0;
	EqualOI[o][i]=chr; EqualOI[o][0] = chr; 
    }
} else {	
    assert(NameRow[0]==NameIRow);
    if(EqualOI[i][0] == 0){ EqualOI[i][o]=chr; EqualOI[i][0] = chr; }
    else if(chr > EqualOI[i][0]){
for(oo = 1; oo <= LenEqOI[0]; oo++) if(EqualOI[i][oo] > 0) break;
fprintf(stderr,"o? %d.%s = %d.%s ",i,NameRow[1][i],oo,NameRow[0][oo]);
fprintf(stderr,"(%d) replaced by %d.%s (%d)\n",EqualOI[i][0],o,NameRow[0][o],chr);
	EqualOI[i][0]=chr; 
	for(oo = 1; oo <= LenEqOI[0]; oo++) EqualOI[i][oo]=0;
	EqualOI[i][o]=chr; EqualOI[i][0] = chr; 
    }
}
#else
if(NameRow[0]==NameORow) EqualOI[o][i]=1; else EqualOI[i][o]=1;
#endif
			bptr += bos;
#if 0
			fprintf(stdout,"%d.%s(%d)==%d.%s(%d)\n",o,NameRow[0][o],co,i,NameRow[1][i],ci);
#endif
			if(hpt){	// condensed = subtree to leaf.
			   for(r=1; r <= Num[0]; r++){
				if(r==o) continue; // column o and row r == '+'.
				if(hpt->RtnHyperPartition(r,o) == '+'){
				    if(Sets[0][r]){
					bos=sprintf(bptr," child: %d.%s (%d).\n",
					      r,NameRow[0][r],CardSet(Sets[0][r]));
				    } else {
					bos=sprintf(bptr," child: %d.%s.\n",
					      r,NameRow[0][r]);
				    } bptr += bos; AddSet(r,match_sets);
				}
			   } 
		        } AddSet(o,match_sets);  NumSubTrees++;
		}
	   }
	} free(FI); free(FO);
	RtnSet=match_sets;
	return NumSubTrees;
}

set_typ	c2h_typ::NewSets(Int4 Num[2], set_typ *Sets[2],char **NameRow[2],hpt_typ *hpt[2])
// find sets in one hierarchy that are missing from the other hierarchy.
{
	Int4		i,j,o,s,r,arg;
	int	bos=0; bptr = &buffer[0]; bptr[0]=0; // buffer pointer and offset;

	set_typ new_sets=MakeSet(Num[0]+1);
	for(o=1; o < Num[0]; o++){
	   if(!Sets[0][o]) continue;
	   if(hpt[0]->TypeOfSet(o) == '?') continue;
	   BooLean IsNew=TRUE;
	   Int4 co = CardSet(Sets[0][o]);
	   for(i=1; i < Num[1]; i++){
	   	if(!Sets[1][i]) continue;
	   	if(hpt[1]->TypeOfSet(i) == '?') continue;
		Int4 ci = CardSet(Sets[1][i]);
		Int4 inter = CardInterSet(Sets[0][o],Sets[1][i]);
		double fo=(double)inter/(double)co;
		double fi=(double)inter/(double)ci;
		if(fo >= 0.10 && fi >= 0.10){ IsNew=FALSE; break; }
	   } if(IsNew){
		AddSet(o,new_sets); 
		bos=sprintf(bptr,"%d.%s (%d)\n",o,NameRow[0][o],co); bptr += bos;
	   }
	} return new_sets;
}

Int4    *c2h_typ::GetTransitiveParentO(FILE *fp,Int4 Root,Int4 *ParentO)
{
	Int4	*Po,p,gp,ggp,N,i; NEW(Po,NumOSets+3,Int4);

	for(i=0; i <= NumOSets; i++) Po[i]=ParentO[i];
	do {
	  for(N=0,i=2; i < NumOSets; i++){  // set all grandchild or lower nodes to related child of root.
		p=Po[i]; if(p == Root || p == 1) continue;	// root is currently direct parent of 'o'.
		gp=Po[p]; if(gp == Root || gp == 1) continue; // root is a parent of 'p'; p is parent of 'o'.
		ggp=Po[gp];
		N++; Po[i]=Po[Po[i]];	// node 'o' is further down the tree; point to its parent's parent.
		if(fp) fprintf(fp,"o(%d) --> %d --> %d --> %d; now --> %d ...\n",i,p,gp,ggp,Po[i]);
	  }	// continue doing this until all nodes pointing to root or to level 1 ancestor.
	} while(N > 0); 
	if(fp){
	   for(i=2; i < NumOSets; i++){	// set all grand children or lower nodes to ancestral child of root node.
		fprintf(fp,"%d.%s --> %d\n",i,NameORow[i],Po[i]);
	   } fprintf(fp,"\n");
	} return Po;
}

Int4	*c2h_typ::MapO2I( )
// rearrange one of the hierarchies to better coincide with the other (for comparison).
{
	assert(CardSet(deleted_osets)==0);
	assert(ohpt->NumSets() == NumOSets);
	FILE	*fp=stderr; fp=0;
	Int4	*Pi=0,*PrntO; ihpt->IsTree(Pi); ohpt->IsTree(PrntO);
	if(fp) ohpt->Put(fp);
	Int4 *Row=SortByOverlapSetO(fp,1,PrntO,Pi); 
	free(PrntO); free(Pi); 
	return Row;
}

Int4	*c2h_typ::SortByOverlapSetO(FILE *fp,Int4 RootO,Int4 *PrntO, Int4 *Pi)
{
	Int4	i,j,o,s,r,bstI,bstCI,bstC,bstCU,co,ci,x;
	double	D,bstD;
	set_typ	SetI,SetO;
	keytyp	key,max_key=DBL_MAX;	// keytyp == double 

	Int4 *Po=GetTransitiveParentO(fp,RootO,PrntO); 

	dh_type dH=dheap(NumOSets+2,4); 
	if(fp) fprintf(fp,"NumOSets=%d\n",NumOSets);
	for(o=2; o < NumOSets; o++){
	   if(misc_osets[o]) SetO=misc_osets[o]; else SetO=osets[o];
	   assert((co=CardSet(SetO)) > 0);
	   if(Po[o] != RootO) continue;	// order only the direct children of RootO.
	   for(bstI=0,bstD=0.0,i=2; i < NumISets; i++){
	   	// if(Pi[i] != 1) continue;	// one of the other nodes will intersect with this.
	   	if(misc_isets[i]) SetI=misc_isets[i]; else SetI=isets[i];
		assert(SetI); ci=CardSet(SetI);
	   	if(ci == 0) print_error("FATAL! c2h_typ::SortByOverlapSetO(): empty leaf nodes disallowed"); 
		Int4 Inter = CardInterSet(SetO,SetI);
		Int4 Union = CardUnionSet(SetO,SetI);
		double	D=(double)Inter/(double)Union;
		if(D > bstD){ bstC=ci; bstD=D; bstI=i; bstCI=Inter; bstCU=Union; }
	   } 
	   if(bstI > 0){
		if(fp) fprintf(fp,"%d.%s(%d) & %d.%s(%d) = %d/%d (%.3f)\n",
			o,NameORow[o],co,bstI,NameIRow[bstI],bstC,bstCI,bstCU,bstD);
		if(bstD > 0.0){ key=(keytyp)bstI+(bstD-0.0001); } else key=max_key;
	   } else {
		if(fp) fprintf(fp,"%d.%s(%d) no matching set found\n",o,NameORow[o],co);
		key=max_key;
	   } insrtHeap(o,key,dH); 
	}

	Int4 *Row; NEW(Row,NumOSets+3,Int4); Row[1]=1; Row[NumOSets]=NumOSets;
	for(j=2; (o=delminHeap(dH)) != 0; j++){
	  Row[j]=o; 
#if 0
	  for(x=2; x < NumOSets; x++){ if(Po[x]==o){ j++; Row[j]=x; } }
#else
	  Int4 *row=0;
	  for(x=2; x < NumOSets; x++){
		if(Po[x]==o){ row=SortByOverlapSetO(fp,o,PrntO,Pi); break;}
	  } if(row){
		for(x=2; x < NumOSets; x++) if(row[x]){ j++; Row[j]=row[x]; }  free(row);
	  }
	 
#endif
	} Nildheap(dH); 
	for(Int4 x=1; x <= NumOSets; x++){
		if(fp) fprintf(fp,"%d.%s --> %d: Row/Col[%d]=%d\n",x,NameORow[x],Po[x],x,Row[x]);
	} free(Po); 
	return Row;
}

void	c2h_typ::ReInit(Int4 *Row)
// this is the same routine as the public version, except for the usage...
{
	Int4		c,f,o,s,x,z,*P;
	
	assert(CardSet(deleted_osets)==0);
	Int4	*s2f; NEW(s2f,NumOSets+10,Int4);  // mapping of set to file!!!
	for(f=s=1; s <= NumOSets; f++,s++){	// WARNING: for empty sets ON_CMA[f] is missing
                if(f > NumOn || strcmp(NameCMSA(ON_CMA[f]),NameORow[s]) !=0) // s == empty set.
			{ s2f[s]=0; f--;  } else s2f[s]=f;
	} 
	hpt_typ	*xhpt=ohpt->Rearrange(FALSE,Row,Row); 
	// rearrange ohpt using the ordering given by Row and Col (==Row) arrays.
	assert(xhpt->IsTree(P));  free(P);
	// xhpt->Put(stderr); // xhpt->PutSorted(stderr);
	delete ohpt; ohpt=xhpt;
	assert(NumORows==ohpt->NumSets()); assert(ohpt->NumSets() == NumOSets);
	for(s=1; s <= ohpt->NumSets(); s++){
		orow[s]=ohpt->RtnHyperPartition(s); NameORow[s]=ohpt->ElmntSetName(s); 
	}

	set_typ	*RtnSets=0,*RtnSetsM=0;
	NEW(RtnSets,NumOSets +4, set_typ); NEW(RtnSetsM,NumOSets +4, set_typ);
	RtnSets[0]=osets[0]; RtnSets[NumOSets]=osets[NumOSets];
 	RtnSetsM[0]=misc_osets[0]; RtnSetsM[NumOSets]=misc_osets[NumOSets];
	cma_typ	*ocma; NEW(ocma,NumOSets+10,cma_typ); 
	for(z=x=1; x <= NumOSets; x++){
		o=Row[x];
		RtnSets[x]=osets[o]; RtnSetsM[x]=misc_osets[o]; 
		if(s2f[o]){ ocma[z]=ON_CMA[s2f[o]]; ON_CMA[s2f[o]]=0; z++; }
		// i.e., If no file corresponding to set 'x' then skip this position in the array.
	} // ocma[NumOSets]=ON_CMA[NumOSets];
	for(s=1; s <= NumOSets; s++){
		osets[s]=RtnSets[s]; misc_osets[s]=RtnSetsM[s];
	   	// fprintf(stderr,"%d.okay?\n",s); // Debug...
		if(s <= NumOn) { assert(ON_CMA[s]==0); ON_CMA[s]=ocma[s]; }
	} free(RtnSets); free(RtnSetsM); free(ocma); free(s2f);

	for(s=2; s < NumOSets; s++){ if(oSets[s]) NilSet(oSets[s]); oSets[s]=0; }
	for(s=2; s < NumOSets; s++){
	   if(!osets[s]) continue;
	   // fprintf(stderr,"Card osets[NumOSets] = %d\n",CardSet(osets[NumOSets]));
	   oSets[s]=CopySet(osets[NumOSets]); // include rejected nodes.
	   for(c=1; c < s; c++){ 	// include root node: c == 1.
		if(orow[s][c] == '+'){
		    if(osets[c]) UnionSet(oSets[s],osets[c]); 
		    // WARNING: Assumes that the input corresponds to a tree!
		}
	   }
	}
}

Int4	c2h_typ::Run()
// note: mode == 's', 'n', 'e','c', 'm'
{
	Int4	c,i,j,o,O,s,r,arg,card,kard,x,z;
	Int4	Num[2];
	set_typ	*Sets[2],*PSet[2];
	char	**NameRow[2];
	hpt_typ	*Hpt[2];
	FILE	*sofp=stdout;
	FILE	*ofp=stdout;

	deleted_osets=MakeSet(NumOSets+1);
	for(s=2; s < NumOSets; s++) if(!osets[s]) AddSet(s,deleted_osets);
	deleted_isets=MakeSet(NumISets+1);
	for(s=2; s < NumISets; s++) if(!isets[s]) AddSet(s,deleted_isets);
	
	wdg_typ treeO=ohpt->RtnAsTree();
	// ReMap the second hierarchy to the first hierarchy...
	// if(CardSet(deleted_osets)==0){ Int4 *Row=MapO2I(); ReInit(Row); free(Row); }
	if(CardSet(deleted_osets)==0){ BackMap=MapO2I(); ReInit(BackMap); }
	else print_error("FATAL: c2h_typ::Run( ) - one or more input sets are missing.");;
NEW(NodeMap,NumOSets+3,Int4);
for(O=1; O <= NumOSets; O++){
	o = BackMap[O]; NodeMap[o]=O;
	// fprintf(stderr," NodeMap[%d] = %d\n",o,O);
} NodeMap[0]=NumOSets; 

	fprintf(sofp,"==========================================================================\n");
	fprintf(sofp,  "============ Evaluating the consistency between hierarchies : ============\n");
	fprintf(sofp,  "==========================================================================\n\n");

	fprintf(sofp,"**********************************************************\n");
	fprintf(sofp,"************** Hierarchy Rooted at %s (%s): **************\n",NameORow[1],onfile);
	fprintf(sofp,"**********************************************************\n\n");
	Num[0]=NumOSets; Num[1]=NumISets; NameRow[0]=NameORow; NameRow[1]=NameIRow;
	PSet[0]=oSets; PSet[1]=iSets;
	Hpt[0]=ohpt; Hpt[1]=ihpt; HPT[0]=ohpt; HPT[1]=ihpt; 

	iHpt=ihpt->Copy(); oHpt=ohpt->Copy();

Int4 Nextra,Nequal,Nclose,N2tree,N2leaf,Nsplit,Nexpand,Ncondense,NsameTree,Nconflict,Nsubtrees,Nsameleaf;
Int4 Mextra,Mequal,Mclose,M2tree,M2leaf,Msplit,Mexpand,Mcondense,MsameTree,Mconflict,Msubtrees,Msameleaf;
	Sets[0]=osets; Sets[1]=isets; 
	new_osets=NewSets(Num,Sets,NameRow,Hpt);
	PrintBuffer(sofp,CardSet(new_osets),"sets are more or less unique to this hierarchy");
Nextra=CardSet(new_osets);

	Sets[0]=osets; Sets[1]=isets; 
	match_osets=IdenticalSets(Num,Sets,NameRow);
	PrintBuffer(sofp,CardSet(match_osets),"nodes are equivalent");
Nequal=CardSet(match_osets);
Nsameleaf= CardInterSet(match_osets,leaf_osets);

	Sets[0]=osets; Sets[1]=isets; 
	fuzzy_osets=IdenticalFuzzySets(Num,Sets,PSet,NameRow);
	PrintBuffer(sofp,CardSet(fuzzy_osets),"nodes are roughly equivalent");
Nclose=CardSet(fuzzy_osets);
Nsameleaf += CardInterSet(fuzzy_osets,leaf_osets);

	Sets[0]=osets; Sets[1]=misc_isets; 
	Fuzzy_osets=IdenticalFuzzySets(Num,Sets,PSet,NameRow);
	x=CardSet(Fuzzy_osets);
	if(x == 1) PrintBuffer(sofp,x,"'fuzzy' node corresponds to a subtree");
	else PrintBuffer(sofp,x,"'fuzzy' nodes each correspond to a subtree");
N2tree=x;

	Sets[0]=misc_osets; Sets[1]=isets; 
	card=IdenticalFuzzySets(Num,Sets,PSet,NameRow,Hpt[0],Fuzzy1_osets); kard=CardSet(Fuzzy1_osets);
	if(card == 1) PrintBuffer(sofp,card,"'fuzzy' subtree (",kard,"nodes) corresponds to a leaf");
	else PrintBuffer(sofp,card,"'fuzzy' subtrees (",kard,"nodes) correspond to leaves");
N2leaf=card;

	Sets[0]=osets; Sets[1]=isets; 
	source_osets=SplitSets(Num,Sets,NameRow,Hpt,split_isets);
	card=CardSet(source_osets); kard=CardSet(split_isets);
	PrintBuffer(sofp,card,"single leaves are split into",kard,"unrelated leaves");
Nsplit=kard;


	Sets[0]=osets; Sets[1]=misc_isets; 
	expanded_osets=IdenticalSets(Num,Sets,NameRow);
	x=CardSet(expanded_osets);
	if(x == 1) PrintBuffer(sofp,x,"leaf expanded into a subtree");
	else PrintBuffer(sofp,x,"leaves each expanded into a subtree");
Nexpand=x;
	
	Sets[0]=misc_osets; Sets[1]=isets; 
	card=IdenticalSets(Num,Sets,NameRow,Hpt[0],condensed_osets); kard=CardSet(condensed_osets);
	PrintBuffer(sofp,card,"subtrees (",kard,"nodes) each correspond to a single leaf");
Ncondense=kard;

	Sets[0]=misc_osets; Sets[1]=misc_isets; 
	m2m_osets=IdenticalSets(Num,Sets,NameRow);
	x=CardSet(m2m_osets); CountSets(Num,Sets,i,j);
	sprintf(tmpstr,"nearly identical subtrees (%d others)",i-x);
	PrintBuffer(sofp,x,tmpstr);
	// PrintBuffer(sofp,CardSet(m2m_osets),"nearly identical subtrees");
Nsubtrees=i;
NsameTree=x;

#if 0
	Sets[0]=misc_osets; Sets[1]=misc_isets; 
	// m2m_osets=IdenticalSets(Num,Sets,NameRow);
	fuzzy_m2m_osets=IdenticalFuzzySets(Num,Sets,PSet,NameRow);
	x=CardSet(m2m_osets); CountSets(Num,Sets,i,j);
	sprintf(tmpstr,"nearly identical subtrees (%d others)",i-x);
	PrintBuffer(sofp,x,tmpstr);
#endif

      FILE *afp=0;
      if(print_tree){
	sprintf(tmpstr,"%s_%s",onfile,infile);
	afp=open_file(tmpstr,"_trees.rtf","w");
	fprintf(afp,"Open as a Windows(default) file, copy and paste as unformated text\n\n");
	fprintf(afp,"  Select a SmartArt horizontal hierarchy\n\n");

	fprintf(afp,"Original (%s) hierarcy:\n",onfile);
	ohpt->PutAsSmartArt(afp); 
	ohpt->PutAsTree(afp); fprintf(afp,"\n\n");

	fprintf(afp,"Remapped (%s) hierarcy:\n",onfile);
      }
#if 1
	// change all that don't start with '~' nor with x_ to "Set86*".
	// and change all x_ to "Set86" only; don't modify the root node.
	for(o=2; o < oHpt->NumSets(); o++){
		if(sscanf(oHpt->ElmntSetName(o),"x_%s",tmpstr) == 1){
			sprintf(tmp_str,"%s",tmpstr); oHpt->ReNameSet(o,tmp_str);
		} else if(sscanf(oHpt->ElmntSetName(o),"~%s",tmpstr) != 1){
			sprintf(tmpstr,"%s*",oHpt->ElmntSetName(o)); 
			oHpt->ReNameSet(o,tmpstr);
		}
	}
#endif
      if(print_tree){
	oHpt->PutAsSmartArt(afp); 
	oHpt->PutAsTree(afp); fprintf(afp,"\n\n");

	fprintf(afp,"Original (%s) hierarcy:\n",infile);
	ihpt->PutAsSmartArt(afp); 
	ihpt->PutAsTree(afp); fprintf(afp,"\n\n");

	fprintf(afp,"Remapped (%s) hierarcy:\n",infile);
      }
#if 1
	for(i=2; i < iHpt->NumSets(); i++){
		if(sscanf(iHpt->ElmntSetName(i),"x_%s",tmpstr) == 1){
			sprintf(tmp_str,"%s",tmpstr); iHpt->ReNameSet(i,tmp_str);
		} else if(sscanf(iHpt->ElmntSetName(i),"~%s",tmpstr) != 1){
			sprintf(tmpstr,"%s*",iHpt->ElmntSetName(i)); 
			iHpt->ReNameSet(i,tmpstr);
		}
	}
#endif
      if(print_tree){
	iHpt->PutAsSmartArt(afp); 
	iHpt->PutAsTree(afp); fprintf(afp,"\n\n");
	fclose(afp);
      }

	fprintf(sofp,"\n**********************************************************\n");
	fprintf(sofp,"************** Hierarchy Rooted at %s (%s): **************\n",NameIRow[1],infile);
	fprintf(sofp,"**********************************************************\n\n");
	Num[1]=NumOSets; Num[0]=NumISets; NameRow[1]=NameORow; NameRow[0]=NameIRow;
	Hpt[1]=ohpt; Hpt[0]=ihpt; PSet[0]=iSets; PSet[1]=oSets;
	HPT[1]=ohpt; HPT[0]=ihpt; delete iHpt; delete oHpt; oHpt=ihpt->Copy(); iHpt=ohpt->Copy();

	Sets[1]=osets; Sets[0]=isets;
	new_isets=NewSets(Num,Sets,NameRow,Hpt);
	PrintBuffer(sofp,CardSet(new_isets),"sets are more or less unique to this hierarchy");
Mextra=CardSet(new_isets);

	Sets[1]=osets; Sets[0]=isets;
	match_isets=IdenticalSets(Num,Sets,NameRow);
	PrintBuffer(sofp,CardSet(match_isets),"nodes are equivalent");
Mequal=CardSet(match_isets);
Msameleaf= CardInterSet(match_isets,leaf_isets);

	Sets[1]=osets; Sets[0]=isets; 
	fuzzy_isets=IdenticalFuzzySets(Num,Sets,PSet,NameRow);
	x=CardSet(fuzzy_isets);
	PrintBuffer(sofp,x,"nodes are roughly equivalent");
Mclose=x;
Msameleaf += CardInterSet(fuzzy_isets,leaf_isets);

	Sets[1]=misc_osets; Sets[0]=isets; 
	Fuzzy1_isets=IdenticalFuzzySets(Num,Sets,PSet,NameRow);
	x=CardSet(Fuzzy1_isets);
	if(x==1) PrintBuffer(sofp,x,"'fuzzy' node corresponds to a subtree");
	else PrintBuffer(sofp,x,"'fuzzy' nodes each correspond to a subtree");
M2tree=x;

	Sets[1]=osets; Sets[0]=misc_isets;
	card=IdenticalFuzzySets(Num,Sets,PSet,NameRow,Hpt[0],Fuzzy_isets); kard=CardSet(Fuzzy_isets);
	if(card==1) PrintBuffer(sofp,card,"fuzzy subtree (",kard,"nodes) corresponds to a leaf");
	else PrintBuffer(sofp,card,"fuzzy subtrees (",kard,"nodes) correspond to leaves");
M2leaf=card;

	Sets[1]=osets; Sets[0]=isets;
	source_isets=SplitSets(Num,Sets,NameRow,Hpt,split_osets);
	card=CardSet(source_isets); kard=CardSet(split_osets);
	if(card==1) PrintBuffer(sofp,card,"single leaf is split into",kard,"unrelated leaves");
	else PrintBuffer(sofp,card,"single leaves are split into",kard,"unrelated leaves");
Msplit=kard;

	Sets[1]=misc_osets; Sets[0]=isets; 
	expanded_isets=IdenticalSets(Num,Sets,NameRow);
	x=CardSet(expanded_isets);
	if(x==1) PrintBuffer(sofp,x,"leaf expanded into a subtree");
	else PrintBuffer(sofp,x,"leaves each expanded into a subtree");
Mexpand=x;

	Sets[1]=osets; Sets[0]=misc_isets;
	card =IdenticalSets(Num,Sets,NameRow,Hpt[0],condensed_isets); kard=CardSet(condensed_isets);
	if(card==1) PrintBuffer(sofp,card,"subtree (",kard,"nodes) corresponds to a single leaf");
	else PrintBuffer(sofp,card,"subtrees (",kard,"nodes) each correspond to a single leaf");
Mcondense=kard;

	Sets[1]=misc_osets; Sets[0]=misc_isets;
	m2m_isets=IdenticalSets(Num,Sets,NameRow);
	x=CardSet(m2m_isets); CountSets(Num,Sets,i,j);
	sprintf(tmpstr,"nearly identical subtrees (%d others)",i-x);
	PrintBuffer(sofp,x,tmpstr);
	// PrintBuffer(sofp,CardSet(m2m_isets),"nearly identical subtrees");
Msubtrees=i;
MsameTree=x;

	// fprintf(stderr,"\nOriginal: "); iHpt->PutAsTree(stderr); 
	// fprintf(stderr,"\nMapped: "); ohpt->PutAsTree(stderr);
	// fprintf(stderr,"\nOriginal: "); ihpt->PutAsTree(stderr);
	// fprintf(stderr,"\nRemapped: "); oHpt->PutAsTree(stderr); fprintf(stderr,"\n");
	delete iHpt; delete oHpt; 
	fprintf(sofp,"\n*********** Summary for %s (%s) versus %s (%s) hierarchies: ***********\n",
			NameORow[1],onfile,NameIRow[1],infile);
	set_typ okay_osets = CopySet(match_osets);
	UnionSet(okay_osets,fuzzy_osets);
	UnionSet(okay_osets,Fuzzy_osets);
	UnionSet(okay_osets,Fuzzy1_osets);
	UnionSet(okay_osets,new_osets);
	UnionSet(okay_osets,expanded_osets);
	UnionSet(okay_osets,condensed_osets);
	UnionSet(okay_osets,m2m_osets);
	UnionSet(okay_osets,source_osets);
	UnionSet(okay_osets,split_osets);
	UnionSet(okay_osets,deleted_osets);
	Int4 okay=CardSet(okay_osets);
	Int4 inconsistent=NumOSets-okay-2;	// subtract  2 for the root and reject nodes.
	fprintf(sofp," %s (%s): %d nodes okay; %d nodes inconsistent.\n",NameORow[1],onfile,okay,inconsistent);
Nconflict=inconsistent;
	if(inconsistent > 0){
	  fprintf(sofp,"  inconsistent sets:\n");
	  for(s=2; s < NumOSets; s++){
		if(!MemberSet(s,okay_osets) && osets[s]){
	  	   char t='.'; if(ohpt->TypeOfSet(s) == '?') t='?';
		   fprintf(sofp,"  %d.%s (%d)%c\n",s,NameORow[s],CardSet(osets[s]),t);
		}
	  } fprintf(sofp,"\n"); fflush(sofp);
	}

	set_typ okay_isets = CopySet(match_isets);
	UnionSet(okay_isets,fuzzy_isets);
	UnionSet(okay_isets,Fuzzy_isets);
	UnionSet(okay_isets,Fuzzy1_isets);
	UnionSet(okay_isets,new_isets);
	UnionSet(okay_isets,expanded_isets);
	UnionSet(okay_isets,condensed_isets);
	UnionSet(okay_isets,m2m_isets);
	UnionSet(okay_isets,source_isets);
	UnionSet(okay_isets,split_isets);
	UnionSet(okay_isets,deleted_isets);
	okay=CardSet(okay_isets);
	inconsistent=NumISets-okay-2;	// subtract  2 for the root and reject nodes.
	fprintf(sofp," %s (%s): %d nodes okay; %d nodes inconsistent.\n",NameIRow[1],infile,okay,inconsistent);
Mconflict=inconsistent;
	if(inconsistent > 0){
	  fprintf(sofp,"  inconsistent sets:\n");
	  for(s=2; s < NumISets; s++){
		if(!MemberSet(s,okay_isets) && isets[s]){
	  	   char t='.'; if(ihpt->TypeOfSet(s) == '?') t='?';
		   fprintf(sofp,"  %d.%s (%d)%c\n",s,NameIRow[s],CardSet(isets[s]),t);
		}
	  } fprintf(sofp,"\n"); fflush(sofp);
	}

Int4  MsameAll, NsameAll,Xclose;
	UnionSet3(match_osets,fuzzy_osets,tmp_osets);
	j=CardInterSetINotJ(inter_osets,tmp_osets);

	UnionSet3(m2m_osets,tmp_osets,tmp_osets);  // tmp_osets = m2m_osets U tmp_osets
	NsameAll=CardSet(tmp_osets);

	UnionSet3(m2m_osets,match_osets,tmp_osets);  Nequal=CardSet(tmp_osets);
// PutSet(stderr,tmp_osets);

	// m2m_osets == internal nodes for equivalent trees; fuzzy_osets = 'close' nodes (all types)
	IntersectNotSet(fuzzy_osets,m2m_osets, tmp_osets); // tmp_osets = fuzzy_osets & ! m2m_osets
	IntersectNotSet(tmp_osets,match_osets); // tmp_osets = tmp_osets & ! match_osets
	Nclose=CardSet(tmp_osets);

	UnionSet3(match_isets,fuzzy_isets,tmp_isets);  // tmp_isets = match_isets U fuzzy_isets
	i=CardInterSetINotJ(inter_isets,tmp_isets); // i = internal nodes and not the same.

	UnionSet3(m2m_isets,tmp_isets,tmp_isets);  // tmp_isets = m2m_isets U tmp_isets
	MsameAll=CardSet(tmp_isets);

	UnionSet3(m2m_isets,match_isets,tmp_isets);  Mequal=CardSet(tmp_isets);
// PutSet(stderr,tmp_isets);

	IntersectNotSet(fuzzy_isets,m2m_isets, tmp_isets); // tmp_isets = fuzzy_isets & ! m2m_isets
	IntersectNotSet(tmp_isets,match_isets); // tmp_isets = tmp_isets & ! match_isets
	Mclose=CardSet(tmp_isets);
	Xclose= MINIMUM(Int4,Nclose,Mclose);

// M = I = real;  N = O = Sim/Rnd.

fprintf(sofp,
"\n%c cdd_id   nodes   equal   close  subtrees   equalST   extra   2tree   2leaf   split   expand condense conflict diff_intern\n",'$');
fprintf(sofp,
"%c %7s	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d  %d\n",
	'$',filename1,NumOSets-1,Nequal+1,Nclose,Nsubtrees,NsameTree,Nextra,N2tree,N2leaf,
		Nsplit,Nexpand,Ncondense,Nconflict,j,NsameAll);
fprintf(sofp,
"%c %7s	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d  %d\n",
	'$',filename2,NumISets-1,Mequal+1,Mclose,Msubtrees,MsameTree,Mextra,M2tree,M2leaf,
		Msplit,Mexpand,Mcondense,Mconflict,i,MsameAll);


#if 0
fprintf(sofp,
"\n$$        ________nodes_______  ___subtrees__  diff_node _leaf2tree_  _tree2leaf_  __split_  __simulated__\n");
fprintf(sofp,
  "$$  cdd_id real sim same close  real sim same  real sim   real sim     real sim    real sim  expand shrink\n",
	'$');

fprintf(sofp,
      "$$ %7s  %3d %3d  %3d   %3d   %3d %3d  %3d   %3d %3d    %3d %3d      %3d %3d     %3d %3d    %3d    %3d\n",
	filename1,NumISets-1, NumOSets-1, Mequal+1,Mclose,Msubtrees, Nsubtrees,MsameTree,Mextra,Nextra,M2tree,
		N2tree,M2leaf,N2leaf,Msplit,Nsplit,Nexpand,Ncondense); 
#else

fprintf(sofp,"$\n$\n");

fprintf(sofp,
"$$        ________nodes_______  ___subtrees__  diff_node _leaf2tree_  __split_  _expanded_  _condensed_  _empty_\n");
fprintf(sofp,
"$$  cdd_id real sim same close  real sim same  real sim   real sim    real sim   real sim    real sim    real sim\n");

fprintf(sofp,
    "$$ %7s  %3d %3d  %3d   %3d   %3d %3d  %3d   %3d %3d    %3d %3d     %3d %3d   %3d  %3d     %3d %3d     %3d %3d\n",
	filename1,NumISets-1,NumOSets-1,Mequal+1,Xclose,Msubtrees,Nsubtrees,MsameTree,Mextra,Nextra,M2tree,
		N2tree,Msplit,Nsplit,Mexpand,Nexpand,Mcondense,Ncondense,Mempty,Nempty); 
#endif

#if 1	// consensus hierarchy...
	FILE *efp=0; // efp=stderr;

	EqualOI[1][1]=1;
	Int4 n,*No; NEW(No,LenEqOI[0] +3,Int4);
	Int4 e,vi,m,*Ni; NEW(Ni,LenEqOI[1] +3,Int4);
	wdg_typ otree=ohpt->RtnAsTree();
	wdg_typ itree=ihpt->RtnAsTree();
	wdg_typ treeI=ihpt->RtnAsTree(); // treeI should use same node numbering as in hpt.
	if(efp){
	   PutSet(efp,match_osets);
	   PutSet(efp,match_isets);
	   fprintf(efp,"     ");
	   for(i=1; i <= LenEqOI[1]; i++){ if(i %5 == 0) fprintf(efp," %4d",i); } fprintf(efp,"\n");
	}
	for(o=1; o <= LenEqOI[0]; o++){
	   n=0;
	   if(efp) fprintf(efp,"%3d: ",o);
	   for(i=1; i <= LenEqOI[1]; i++){
	   	if(EqualOI[o][i] > 0){
		   if(efp) fprintf(efp,"*"); n++;
	   	   e=FirstInWdgraph(NodeMap[o],treeO); SetWeightWdgraph(i,e,treeO);
	   	   e=FirstInWdgraph(o,otree); SetWeightWdgraph(i,e,otree);
	   	   e=FirstInWdgraph(i,itree); SetWeightWdgraph(o,e,itree);
	   	   e=FirstInWdgraph(i,itree); SetWeightWdgraph(BackMap[o],e,treeI);
		   Ni[i]++; No[o]++; 
		}
		else if(efp) fprintf(efp," ");
	   } if(efp) fprintf(efp," (%d)\n",n);
	} if(efp) fprintf(efp,"\n     ");
	for(n=0,i=1; i <= LenEqOI[1]; i++){
		if(Ni[i] > 0){ if(efp) fprintf(efp,"%d",Ni[i]); n++; } else if(efp) fprintf(efp,"0");
	} if(efp) fprintf(efp," (%d)\n",n);
	EqualArray[1]=Ni; EqualArray[1][0]=LenEqOI[1];
	EqualArray[0]=No; EqualArray[0][0]=LenEqOI[0];
#if 0	// for consensus hierarchy...
	for(o=1; o <= NumOSets; o++){ 
		fprintf(stderr,"%d -> %d\n",o,NodeMap[o]);
	}
	fprintf(stderr,"Original %s tree: # nodes = ",onfile);
	PutWdgraph(stderr,treeO); 
	fprintf(stderr,"remapped %s tree: (parent, child, \"%s\" node == child); # nodes = ",onfile,infile);
	PutWdgraph(stderr,otree); 
	fprintf(stderr,"%s tree: # nodes = ",infile);
	PutWdgraph(stderr,itree); 
	fprintf(stderr,"%s tree: (parent, child, \"%s\" node == child); # nodes = ",infile,onfile);
	PutWdgraph(stderr,treeI); 
#endif
	NilWdgraph(otree); NilWdgraph(itree); 
	// NilWdgraph(treeO); NilWdgraph(treeI);
	TreeO=treeO; TreeI=treeI; 
	// exit(1);
#endif


	fprintf(sofp,"\n==============================================================================\n");
	fprintf(sofp,"============ Mapping of %s (%s) hierarchy to %s (%s) hierarchy: ============\n\n",
			NameORow[1],onfile,NameIRow[1],infile);
	NilSet(match_isets); NilSet(okay_isets);
	NilSet(fuzzy_isets); NilSet(Fuzzy_isets); NilSet(Fuzzy1_isets);
	NilSet(new_isets); NilSet(expanded_isets); NilSet(condensed_isets);
	NilSet(m2m_isets); NilSet(source_isets); NilSet(split_isets); NilSet(deleted_isets);

	NilSet(match_osets); NilSet(okay_osets);
	NilSet(fuzzy_osets); NilSet(Fuzzy_osets); NilSet(Fuzzy1_osets);
	NilSet(new_osets); NilSet(expanded_osets); NilSet(condensed_osets);
	NilSet(m2m_osets); NilSet(source_osets); NilSet(split_osets); NilSet(deleted_osets);

	BooLean	success = CompareCMSAs(ofp);
	if(print_tree){
	  double runtime = difftime(time(NULL),time1);
	  fprintf(stderr,"time = %0.2f minutes (%0.1f seconds)\n",runtime/60.0,runtime);
	}
	if(success) return 0; else return 1;
}

