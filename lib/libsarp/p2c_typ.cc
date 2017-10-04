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

#include "psc_typ.h"

void	p2c_typ::Init()
{
	LastARSD=0.0; LastNumData=0;
	FindBest=FALSE; K1=0,K2=0; MaxMeanDist=20; MaxSqDist=1000; MinDistInSeq=8;
        MinDist=3;      // what about inserts in some proteins?
        begin=0,end=0; Begin=0,End=0; KeyCol=0; MinVar=0.0; bin_size=1.0;
        HA_dmax=0.0, dmax=0.0;
	Target=40; MaxDataPoints=2000;
}

BooLean	p2c_typ::AddColToSeq(Int4 S, Int4 A, Int4 *ColToSeq)
// returns true if new repeat was found and added; false if already present.
{
	Int4 num_cols=NumCol( ),col;
	Int4 R;
	for(R=1; R <= NumFullRpts[S]; R++){
	 	Int4 *Tmp=Col2FullSeq[S][R];
	  	BooLean different=FALSE;
		for(Int4 col=1; col <= num_cols; col++){
		   if(ColToSeq[col] && Tmp[col] && ColToSeq[col] != Tmp[col]){
			different=TRUE; break;
		   }
		} if(!different){ return FALSE; }
	} NumFullRpts[S] = R;
	if(R > MAX_NUMBER_INTERNAL_RPTS) print_error("Too many internal repeats"); 
	RptCategory[S][R]=A;
	Col2FullSeq[S][R]=ColToSeq;
	AddSet(A,RelevantSet[S]);
	return TRUE;
}

BooLean	p2c_typ::AddColToSeq(Int4 i, Int4 S, Int4 A, Int4 *ColToSeq)
// returns true if new repeat was found and added; false if already present.
{
	Int4 num_cols=NumCol( ),col;
	assert(i > 0 && i <= esc->NumPDB_Set[S]);
	Int4 R,I=esc->PDB_SetI[S][i],C=esc->PDB_SetC[S][i];
	for(R=1; R <= NumRpts[I][C]; R++){
	  	Int4 *Tmp=Col2pdbSeq[I][C][R];
		BooLean	different=FALSE;
		for(Int4 col=1; col <= num_cols; col++){
		   if(ColToSeq[col] && Tmp[col] && ColToSeq[col] != Tmp[col]){
			different=TRUE; break;
		   } // ignores 'X' residue in column positions.
		} if(!different){ return FALSE; } // if 
	} NumRpts[I][C] = R;
	if(R > MAX_NUMBER_INTERNAL_RPTS) print_error("Too many internal repeats"); 
	Col2pdbSeq[I][C][R]=ColToSeq;
	// AddSet(A,RelevantSet[S]);
	RelevantSeq[I][C][R]=A;
	return TRUE;
}

void    p2c_typ::MapSqAlnToStruct( )
// Call only after calling esc_typ( )
//************************ 3. Find pdb sequences in cma file ***************************
//************************ and get mapping between seqs ***************************
// Note that 
{
   Int4	I,C,i,j,c1,c2,col,real,R,pdb_real,os,os_cma,os_pdb,N,S,s,sq,A;
   Int4	NumCol=LengthCMSA(1,IN_CMA[1]);
   e_type pdbE,cmaE,*csq;
   char str[58];
   Int4 Score;

   NEW(RelevantSet,esc->NumPDB_Sets + 5,set_typ); 
   NEWP(RptCategory,esc->NumPDB_Sets + 5,Int4);
   NEWPP(Col2FullSeq,esc->NumPDB_Sets + 5,Int4);
   NEW(NumFullRpts,esc->NumPDB_Sets + 5,Int4);
   for(S=1; S <=esc->NumPDB_Sets; S++){
        NEW(RptCategory[S],MAX_NUMBER_INTERNAL_RPTS + 5,Int4);
   	RelevantSet[S] = MakeSet(Number + 5); ClearSet(RelevantSet[S]);
   	NEWP(Col2FullSeq[S],MAX_NUMBER_INTERNAL_RPTS+ 5,Int4);
   }
   NEWPP(RelevantSeq,mpdb->NumPDB + 5,Int4); // is there a corresponding seq. in cma file
   NEWP3(Col2pdbSeq,MAX_NUM_PDB_INPUT_FILES + 5,Int4); // mapping of pdb residue # to cma column #
   NEWP(NumRpts,mpdb->NumPDB +5, Int4);
   for(I=1; I <=mpdb->NumPDB; I++){
	NEWP(RelevantSeq[I],nChainsPDB(mpdb->pdb[I]) + 3, Int4);
	NEWPP(Col2pdbSeq[I],nChainsPDB(mpdb->pdb[I]) + 3, Int4);
   	NEW(NumRpts[I],nChainsPDB(mpdb->pdb[I]) +5, Int4);	// NumRpts[I][C] = 0;
   	for(C=1; C <=nChainsPDB(mpdb->pdb[I]); C++){
		pdbE=mpdb->pdbSeq[I][C];
		if(pdbE && LenSeq(pdbE) >= esc->MinSeqOverlap){
		    NEW(RelevantSeq[I][C],MAX_NUMBER_INTERNAL_RPTS + 3, Int4);
		    NEWP(Col2pdbSeq[I][C],MAX_NUMBER_INTERNAL_RPTS +3, Int4);
		    // NEWP(Col2pdbSeq[I][C],NumCol + 3, Int4);
		}
	}
   }
#if 1	// Diagnostics...
   enum fault { deleted=0, disjoint=1, missing=2, x_over=3, gaps=4, off=5, okay=6, highX=7 };
   Int4	Fault[10]; for(i=0; i <= 9; i++) Fault[i]=0;
#endif
   NEW(csq,Number +3, e_type);
   for(A=1; A < Number; A++){ if(IN_CMA[A]) csq[A]=MkConsensusCMSA(IN_CMA[A]); }
   for(S=1; S <=esc->NumPDB_Sets; S++){
	//======================== Find column residue positions in Full pdb seq =============
	pdbE=esc->FullSeq[S];
	if(pdbE==0 || LenSeq(pdbE) < esc->MinSeqOverlap){ continue; }
	// fprintf(stderr,"******************** Set %d (pdb[%d][%d]) ********************\n",S,I,C);
// PutSeq(stderr,pdbE,AB);
	//*************** check over all CMA files except last, reject subgroup. ****************
#if 0
	for(A=1; A < Number; A++)
#else
	for(A=Number-1; A > 0; A--)	// search backwards to favor subset assignments...
#endif
	{
	     cma=IN_CMA[A]; 
	     if(cma==0){ continue; } // there were no sequences in this Misc set; no need to check.
	     assert(LengthCMSA(1,cma) == NumCol); N = NumSeqsCMSA(cma);
	     // if(FastAlnSeqSW(12,4,csq[A],pdbE,AB) < 10) continue;	// if 
	     for(sq=1; sq <= N; sq++){	//========== looking through sequences in alignment. ======
	  	cmaE = TrueSeqCMSA(sq,cma);	
		if(LenSeq(cmaE) < esc->MinSeqOverlap){ continue; }
		Int4 NumX,adjust=0,MaxNumX;
		// char rtn=IsSameSeqFast(pdbE,cmaE,&os,&NumX,esc->MinSeqOverlap); // ignores 'X' residues...
		char rtn=IsSameSeqFastX(pdbE,cmaE,&os,&NumX,esc->MinSeqOverlap); // ignores 'X' residues...
		if(rtn==0) continue; 
		MaxNumX=(Int4)floor(((double) esc->MinSeqOverlap*0.33));
		// fprintf(stderr," !!!!!!! IsSameSeq(): Set = %d; sq=%d; os=%d; NumX = %d (%d) !!!!!!!\n",
		//	A,sq,os,NumX,MaxNumX);
		if(NumX > MaxNumX) continue;			// Ignore if 2/3rds of residues are 'X's
#if 1	
	if(ColPairSet){
		assert(Number== 2);
		if(MapPdb2cmaSq==0){ NEWP(MapPdb2cmaSq,mpdb->NumPDB + 5,Int4); }
		for(Int4 id=1; id <= esc->NumPDB_Set[S]; id++){	// label the rest of sequences.
		   I=esc->PDB_SetI[S][id]; C=esc->PDB_SetC[S][id]; // C == pdb chain; I = pdb file ID.
		   assert(I <= mpdb->NumPDB && I > 0); 
		   pdb_typ P=mpdb->pdb[I]; assert(C <= nChainsPDB(P) && C > 0);
		   if(MapPdb2cmaSq[I]==0) NEW(MapPdb2cmaSq[I],nChainsPDB(P)+5,Int4);
	   	   e_type pdbIC=mpdb->pdbSeq[I][C]; assert(pdbIC);
		   Int4 os2;
		   char rtn2=IsSameSeqFastX(pdbIC,cmaE,&os2,&NumX,esc->MinSeqOverlap); // ignoring 'X' residues...
		   if(rtn2 && NumX <= MaxNumX){
			MapPdb2cmaSq[I][C]=sq;
		   }
	        }
	}
#endif
		StrSeqID(str,50,pdbE); 
		// fprintf(stderr," a match: %s --> \"%s\"(%d)\n",str,NameCMSA(cma),A);
		// Fixed this so that overhanging regions within sequences are allowed.
		if(rtn == 1) adjust = -os;	// pdbE N-terminus starts within cmaE.
		else if(rtn == 2) adjust = os;	// cmaE N-terminus starts within pdbE.
		else print_error("p2c_typ::MapSqAlnToStruct: this should not happen!");
		os_cma=OffSetSeq(cmaE); Score=PseudoAlnScoreSqToCMSA(csq[A],sq,cma);
		Int4 *TmpColToSeq; NEW(TmpColToSeq ,NumCol + 3, Int4); // temporary array for pdbE.
		for(col = 1; col <= NumCol; col++){
			if(IsDeletedCMSA(1,sq,col,cma)) { Fault[deleted]++; continue; } // assumes single block
// WARNING!!!: need to fix IsDeletedCMSA(sq,col,cma) within cmsa.cc !!!!  afn: 12_27_2010.
			// if(IsDeletedCMSA(sq,col,cma)) continue; // this function has problems!!!
			// NOTE: TruePosCMSA() returns the position in real seq w/o offset
			// col is position in block 1 (assumes only one block) within fake seq
			i=TruePosCMSA(sq,col,cma);	// ignores offset... 
			if(i < 1 || RealToFakeCMSA(sq,os_cma+i,cma) == 0){ Fault[disjoint]++; continue; }
			          // ^ this == 0 if no corresponding position in fake seq.
			j = i+adjust;         // <-- corresponding position in pdbE;
			if(!(j > 0 && j <= LenSeq(pdbE))){ Fault[missing]++; continue; }
					// ^ implies that pdbE lacks these positions
			c1=AlphaChar(ResSeq(i,cmaE),AB); c2=AlphaChar(ResSeq(j,pdbE),AB);
			assert(c1 == 'X' || c2 == 'X' || c1 == c2);
			TmpColToSeq[col]=j;	// position within full pdb sequence.
			Fault[okay]++;
			// assert(TmpColToSeq[col]==j); // run this with old method turned on as a check.
		}
		if(AddColToSeq(S,A,TmpColToSeq)){	// adds columns to Col2FullSeq[S][R]=TmpColToSeq;
			// Then A added to RelevantSet[S] & RptCategory[S][R]=A; increments NumFullRpts[S].
			// ^ this returns TRUE iff TmpCols are unique (not found for this pdb before).
			// TRUE implies a repeat has been added.
			// inserts next to deletions are a problem! Run "tweakcma -iron"
#if 0	// for debugging...can all be turned off.
			if(TmpColToSeq[NumCol]==0){
		 	   for(col = 1; col <= NumCol; col++){
				c1=AlphaChar(ResSeq(TmpColToSeq[col],pdbE),AB);
				fprintf(stderr,"%d: %c%d\n",col,c1,TmpColToSeq[col]);
			   }
			}
			c1=AlphaChar(ResSeq(TmpColToSeq[1],pdbE),AB);
			c2=AlphaChar(ResSeq(TmpColToSeq[NumCol],pdbE),AB);
			// PutSubSeq(stderr,TmpColToSeq[1],TmpColToSeq[NumCol],pdbE,AB);
			PutSeq(stderr,pdbE,AB);
			PutSeq(stderr,cmaE,AB);
			fprintf(stderr,
			    "!!!!!!!! IsSameSeq(): Set=%d; score=%d; os=%d; col(1)=%c%d; col(%d)=%c%d) !!!!!!\n",
				A,Score,os,c1,TmpColToSeq[1],NumCol,c2,TmpColToSeq[NumCol]);
			if(rtn==2) PutDiagonalSeq(stderr,os,pdbE,cmaE,AB);
        		else PutDiagonalSeq(stderr,os,cmaE,pdbE,AB);
#endif
			// break; // ignore any duplicate sequences that might be in cma
		} else free(TmpColToSeq); 
	     } 
        }
	//=================== Find column residue positions in each pdb structure =============
	if(CardSet(RelevantSet[S]) > 0){	// Add the rest of the sequences in the set.
	   Int4		id,NumX,Rj,Start,End,x,y;
	   e_type	pdbIC,pdbS;
	   for(R=1; R <= NumFullRpts[S]; R++){
		A=RptCategory[S][R];
		//*************** Find Start and End for repeat R within pdb S ********************
		for(x=1; x <= NumCol; x++){ Start=Col2FullSeq[S][R][x];  if(Start !=0) break; }
 // if(x > NumCol) break;	// Temporary fix; need to see why this is occurring; run purify!!
		if(x > NumCol){ 
#if 0
			PutSeq(stderr,pdbE,AB); 
			fprintf(stderr,"x=%d; NumCol=%d; S=%d; R=%d\n",x,NumCol,S,R);
#endif
			Fault[x_over]++; continue;
		} assert(x <= NumCol); 
		for(x=NumCol; x > 0; x--){
		   End=Col2FullSeq[S][R][x];  if(End!=0) break; 
		}
		if(!(Start < End && End > 0)){
			fprintf(stderr,"x=%d; NumCol=%d; S=%d; R=%d; Start=%d; End=%d\n",x,NumCol,S,R);
			Fault[off]++; continue;
		} assert(Start < End && End > 0);
		if(!(Start > 0 && End <= LenSeq(FullSeq[S]))){
			PutSeq(stderr,FullSeq[S],AB);
			fprintf(stderr,"x=%d; NumCol=%d; S=%d; R=%d; Start=%d; End=%d; LenSeq(S)=%d\n",
				x,NumCol,S,R,LenSeq(FullSeq[S]));
			assert(Start > 0 && End <= LenSeq(FullSeq[S]));
	 	} pdbS=MkSubSeq(Start,End,FullSeq[S]);	// pdbS = subseq of aligned pdb seq.
		for(id=1; id <= esc->NumPDB_Set[S]; id++){	// label the rest of sequences.
		   I=esc->PDB_SetI[S][id]; C=esc->PDB_SetC[S][id]; // C == pdb chain; I = pdb file ID.
		   pdbIC=mpdb->pdbSeq[I][C]; assert(pdbIC);
		   // char rtn=IsSameSeqFast(pdbS,pdbIC,&os,&NumX,esc->MinSeqOverlap); // ignoring 'X' residues...
		   char rtn=IsSameSeqFastX(pdbS,pdbIC,&os,&NumX,esc->MinSeqOverlap); // ignoring 'X' residues...
		   Int4 MaxNumX = (Int4) floor(((double) esc->MinSeqOverlap*0.33));
		   // Int4 MaxNumX = (Int4) floor(((double) LenSeq(pdbS)*0.25));
		   if(rtn && NumX <= MaxNumX){
			Int4 *TmpColToSeq; NEW(TmpColToSeq ,NumCol + 3, Int4); // temporary array for pdbIC
			for(col=1; col <= NumCol; col++){	//
				x = Col2FullSeq[S][R][col];
				if(x){
				    y = x - OffSetSeq(pdbIC); 
				    if(y < 1 || y > LenSeq(pdbIC)){ Fault[gaps]++; continue; }  // gaps at end...
				    if(ResSeq(y,pdbIC)){ TmpColToSeq[col]=y; }
				    unsigned char r_x,r_y;
				    r_y=ResSeq(y,pdbIC); r_x=ResSeq(x,FullSeq[S]);
				    if(r_x && r_y && r_x != r_y){
					PutSeq(stderr,FullSeq[S],AB);
					fprintf(stderr,
					  "x=%d;y=%d;NumCol=%d;S=%d;R=%d;Start=%d;End=%d;LenSeq(S)=%d\n",
                                		x,y,NumCol,S,R,LenSeq(FullSeq[S]));
					this->DebugMapSqAln2Strct(stderr,S,C,R,I,TmpColToSeq,pdbS,pdbIC,pdbE);
					assert(ResSeq(y,pdbIC) == ResSeq(x,FullSeq[S]));
				    }
				}
		        }
			// assert(AddColToSeq(id,S,A,TmpColToSeq)); 
			// id is passed in in order to retrieve I, C to check for repeats.
			if(!AddColToSeq(id,S,A,TmpColToSeq)){ // sets Col2pdbSeq[I][C][R]= TmpColToSeq;
			   free(TmpColToSeq); continue; // just skip this; not sure that it matters?
			   // In 'id' mode: if(TRUE) TmpCol
			   fprintf(stderr,"!!!!!!!!!!!!!!!!!! ERROR (%d:%s) !!!!!!!!!!!!!!!!!!\n",A,NameCMSA(cma));
			   this->DebugMapSqAln2Strct(stderr,S,C,R,I,TmpColToSeq,pdbS,pdbIC,pdbE);
			   print_error("p2c_typ::MapSqAlnToStruct( ) error: Redundant internal repeats?");
			} 
			// assert(AddColToSeq(id,S,A,Col2FullSeq[S][R])); 
		   } else if(rtn && NumX > MaxNumX) Fault[highX]++; 
		} NilSeq(pdbS);
	   }
	}
   }
   char Str[300];
   sprintf(Str,"deleted=%d, disjoint=%d, missing=%d, x_over=%d, gaps=%d, off=%d; highX=%d;  okay=%d.\n",
			Fault[0],Fault[1],Fault[2],Fault[3],Fault[4],Fault[5],Fault[7],Fault[6]);
   Diagnostic=AllocString(Str);
   for(A=1; A < Number; A++) if(csq[A]) NilSeq(csq[A]); free(csq);
}

#define P2C_USAGE_START \
"USAGE: p2c_typ(argc, *argv[], <pdb_name_file>,<cmafile>)\n\
   argument [-options]\n\
    -best                       Find the atom pairs that deviate least (default: deviate most)\n\
    -beta                       Use beta carbons for non-glycine residues (default: alpha carbons only)\n\
    -col=<int1>                 compute all pairs that include columns <int1>\n\
    -D<float>                   dmax for non-classical H-bonds (default: case dependent)\n\
    -d<float>                   dmax in Angstroms for classical H-bonds (default: 2.5 Angstoms)\n\
    -show=<int>:<int>           Show residues in columns <int1> and <int2>\n\
    -srch=<int1>:<int2>         Search for residues in sequence <int1> that are a better 3D-fit to \n\
                                average distances of other residues within 10 Angstroms of \n\
                                column <int2> residues \n\
    -range=<int1>:<int2>        Only look at columns <int1> to <int2>\n\
    -Range=<int1>:<int2>        Focus on columns <int1> to <int2>\n\
    -bin=<real>                 Set the bin size for histogram (default: 1.0)\n\
    -P=<filename>               Pattern file corresponding to sequence subgroup.\n\
    -maxdist=<int>              Set the Maximum Mean distance to consider (default: 40)\n\
    -seqdist=<int>              Set the minimum distance between residues to consider (default: 4)\n\
    -maxsqdist=<int>            Set the maximum distance between residues to consider (default: 200) \n\
    -mindist=<int>              Set the minimum distance between aligned columns to compare (default: 3)\n\
    -minvar=<int>               Set the minimum variance to consider (default: 0)\n\
\n\n"

void	p2c_typ::GetArg(Int4 argc, char *argv[])
// *************** Get arguments for all program options **********************
{
	Int4 arg,i;
        // if(argc < 3) print_error(P2C_USAGE_START);
        for(arg = 0; arg < argc; arg++){
		fprintf(stdout,"%s ",argv[arg]);
	} fprintf(stdout,"\n");
        for(arg = 1; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'b':
		if(sscanf(argv[arg],"-bin=%lf",&bin_size) == 1){
		     // fprintf(stderr,"binsize = %f\n",bin_size);
		     if(bin_size < 0.001) print_error(P2C_USAGE_START);
		     else if(bin_size > 4.0) print_error(P2C_USAGE_START);
		} else if(strcmp("-best",argv[arg])==0) {
			FindBest=TRUE;
		} else if(strcmp("-beta",argv[arg])==0) {
			UseBeta=TRUE;
		} else print_error(P2C_USAGE_START);
		break;
             case 'c':
		if(sscanf(argv[arg],"-col=%d",&KeyCol) == 1){
		   if(KeyCol < 1) print_error(P2C_USAGE_START);
		} else print_error(P2C_USAGE_START);
		break;
	     case 'D':
                if(sscanf(argv[arg],"-D%f",&dmax) != 1 || dmax > 100 || dmax < 0){
                        print_error(P2C_USAGE_START);
                } break;
             case 'd':
                if(sscanf(argv[arg],"-d%f",&HA_dmax) != 1 || HA_dmax > 100 || HA_dmax < 0)
                        print_error(P2C_USAGE_START);
                break;
             case 'm':
		// if(sscanf(argv[arg],"-maxdist=%lf",&MaxMeanDist) != 1)
		if(sscanf(argv[arg],"-minvar=%lf",&MinVar) != 1){
		 if(sscanf(argv[arg],"-maxdist=%d",&MaxMeanDist) != 1){
		   if(sscanf(argv[arg],"-mindist=%d",&MinDist) != 1) print_error(P2C_USAGE_START);
		   else if(MinDist < 1) print_error(P2C_USAGE_START);
		 } else if(MaxMeanDist < 1) print_error(P2C_USAGE_START);
		} else if(MinVar < 0.0) print_error(P2C_USAGE_START);
		break;
             case 'R':
		if(sscanf(argv[arg],"-Range=%d:%d",&Begin,&End) == 2){
			if(Begin > End || Begin <= 0) print_error(P2C_USAGE_START);
		} else print_error(P2C_USAGE_START);
		// print_error("-R option not yet implemented");
		break;
             case 'r':
		if(sscanf(argv[arg],"-range=%d:%d",&begin,&end) == 2){
			if(begin >= end || begin <= 0) print_error(P2C_USAGE_START);
		} else print_error(P2C_USAGE_START);
		break;
             case 's':
		if(sscanf(argv[arg],"-seqdist=%d",&MinDistInSeq) == 1){
			if(MinDistInSeq < 1) print_error(P2C_USAGE_START);
		} else if(sscanf(argv[arg],"-show=%d:%d",&K1,&K2) == 2){
			if(K1 > K2){ i = K1; K1 = K2; K2 = i; }
			else if(K1 == K2) print_error(P2C_USAGE_START);
		} else print_error(P2C_USAGE_START);
		break;
             default : print_error(P2C_USAGE_START);
           }
	 }
	}
}

int	p2c_typ::PrintKLST_Files(char call )
{
   //************* Print out vsi files for each set.
   Int4 A,x,p,n,i,j,I,C,R,S;
   Int4 vsi_number=0;
   char side_color[]="WYROMGCBDLWWWWWWWWWWWWWWWWWWDDDDDDDDDDDDDDDDDDD";
   Int4 MaxColor=40;
   double maxdist=4.0;
   set_typ SetP=0;

   if(SideColors){
	for(j=0,i=1; isalpha(SideColors[j]); j++,i++){
		side_color[i]=SideColors[j]; if(i >= MaxColor) break;
	}
   }
   FILE *vfp=0;
   for(S=1; S <= esc->NumPDB_Sets; S++){
 	char vsifile[200];
	if(CardSet(RelevantSet[S]) == 0) continue;	// skip irrelevant files.
	assert(esc->PDB_SetI[S]);
	vsi_number++;
	if(vfp==0){
	    if(call > 0) vfp=open_file(OutFile,".sprs","w");
	    else vfp=open_file(OutFile,"_pdb.klst","w");
	}
	// fprintf(vfp,"~$=%d.\n",vsi_number); fflush(vfp);

	for(R=1; R <= NumFullRpts[S]; R++){
#if 1
          for(j=1; esc->PDB_SetI[S][j]; j++){
		I=esc->PDB_SetI[S][j]; C=esc->PDB_SetC[S][j];
		fprintf(vfp,"file: %s\n",mpdb->pdb_file[I]);
		fprintf(vfp,"chain: %c.\n",ChainCharPDB(C,mpdb->pdb[I]));
	  }
#endif
	  Int4 Row = RptCategory[S][R]; // hpt row...
	  assert(MemberSet(Row,RelevantSet[S]));
	  // if(!MemberSet(A,RelevantSet[S])) continue;	// skip irrelevant files.
	  Int4 X;
          // for(X=0,x=1; x <= ptrn->NumPttrns; x++) 
          for(X=0,x=ptrn->NumPttrns;  x > 0; x--)	// stored backwards...
	  {
	    Int4 Col=ptrn->PttrnCategory[x];
	    if(hpt->Cell(Row, Col) != '+') continue; else X++;
	    // print out subgroups before supergroups to ensure proper color.
	    if(X > 40) continue;  // ran out of colors.
	    if(X==1){
	      Int4 start=1,TheEnd=NumCol(),num_ins,num_del;
	      for(j=1; Col2FullSeq[S][R][j] == 0 && j <= TheEnd; j++){ } start=j;
	      for(j=TheEnd; Col2FullSeq[S][R][j] == 0 && j > 0; j--){ } TheEnd=j;
	      for(j=1,num_del=0; j <= NumCol(); j++){ if(Col2FullSeq[S][R][j] == 0) num_del++; }
	      Int4 nAln,strt=Col2FullSeq[S][R][start],end=Col2FullSeq[S][R][TheEnd];
	      nAln=NumCol() - num_del; num_ins=(end-strt+1) - nAln;
#if 1	// Find residue Set...
	      SetP=MakeSet(end +9);
	      for(j=start; j <= TheEnd; j++){ 
		Int4 res_j= Col2FullSeq[S][R][j];
		if(res_j != 0){ assert(res_j >= strt && res_j <= end);  AddSet(res_j,SetP); }
	      }
	      Int4 low,high;
	      char *rtn=RtnStrSet(SetP,low,high);
	      fprintf(stderr,"Set string = \"%s\"; range: %d-%d\n",rtn,low,high);
	      NilSet(SetP); SetP=0;
	      fprintf(vfp,"range: %s(%d;%d).\n",rtn,nAln,num_ins); free(rtn);
#else
	      fprintf(vfp,"range: %d_%d(%d;%d).\n",strt,end,nAln,num_ins);
#endif
	    }
            BooLean first=TRUE;
	    for(p=1; p <= ptrn->NumPttrnRes[x]; p++){
	        // if(p > 10) continue;  // skip less significant.
		Int4 col=ptrn->PosPttrn[x][p];
		Int4 site=Col2FullSeq[S][R][col];
		if(site == 0) continue;  // not visible within structures.
		assert(site <= LenSeq(FullSeq[S]));
		Int4 r=ResSeq(site,FullSeq[S]);
		char Res=AlphaChar(r,AB);
		if(strchr(ptrn->PttrnRes[x][p],Res)){
		   // Int4 X=ptrn->NumPttrns-x+1; 
		   if(first){ fprintf(vfp,"%c=%d",side_color[X],site); first=FALSE; }
		   else fprintf(vfp,",%d",site);
		} else {	// print mismatches in white.
#if 0
		   // See whether shows up at a higher level:
		   BooLean IsHigher=FALSE;
		   for(Int4 x2=x+1;  x2 > 0; x2--){
	    	      Int4 Col2=ptrn->PttrnCategory[x2];
	    	      if(hpt->Cell(Row, Col2) != '+') continue;
	    	      for(Int4 p2=1; p2 <= ptrn->NumPttrnRes[x2]; p2++){
			Int4 col2=ptrn->PosPttrn[x2][p2];
			if(col == col2){ IsHigher=TRUE; break; }
		      } if(IsHigher) break;
		   } if(IsHigher) continue;
#endif
		   // if no higher level pattern, then print out.
		   if(first){ fprintf(vfp,"%c=%d",side_color[X],site); first=FALSE; }
		   else fprintf(vfp,",%d",site);
	       	} 
	  } fprintf(vfp,"\n");
	} fprintf(vfp,"\n");
      }
   } if(vfp){ fclose(vfp); return 1; } else return 0;
}

void	p2c_typ::PrintVSI_Files(FILE *vsifp)
{
   //************* Print out vsi files for each set.
   Int4 A,x,p,n,i,j,I,C,R,S;
   Int4 vsi_number=0;
   char trace_color[]="WYROMGCBDLWWWWWWWWWWWWWWWWWWDDDDDDDDDDDDDDDDDDD";
   char side_color[]="WYROMGCBDLWWWWWWWWWWWWWWWWWWDDDDDDDDDDDDDDDDDDD";
   Int4 MaxColor=40;
   double maxdist=4.0;

   if(TraceColors){
	for(j=0,i=1; isalpha(TraceColors[j]); j++,i++){
		trace_color[i]=TraceColors[j]; if(i >= MaxColor) break;
	}
   }
   if(SideColors){
	for(j=0,i=1; isalpha(SideColors[j]); j++,i++){
		side_color[i]=SideColors[j]; if(i >= MaxColor) break;
	}
   }
   FILE *vfp=0; 
   for(S=1; S <= esc->NumPDB_Sets; S++){
 	char vsifile[200];
// if(mpdb->pdb[S] == 0) continue;	// need to know why this is null for some!!! --> WRONG INDEX!!!
// fprintf(stderr,"%da: %s \n",S,FilenamePDB(mpdb->pdb[S]));
	if(CardSet(RelevantSet[S]) == 0) continue;	// skip irrelevant files.
	assert(esc->PDB_SetI[S]);
	vsi_number++;
	if(vsifp == 0 && vfp==0) vfp=open_file(OutFile,"_pdb.VSI","w");
#if 0
	sprintf(vsifile,"_pdb%d.vsi",vsi_number);
	// sprintf(vsifile,"_pdb%d.vsi",S);
	vfp=open_file(OutFile,vsifile,"w");
#else
	if(vsifp){ fprintf(vsifp,"~$=%d.\n",vsi_number); fflush(vsifp); }
	else { fprintf(vfp,"~$=%d.\n",vsi_number); fflush(vfp); }
#endif

        for(j=1; esc->PDB_SetI[S][j]; j++){
		I=esc->PDB_SetI[S][j]; C=esc->PDB_SetC[S][j];
		if(vsifp){
		    fprintf(vsifp,"File%d=%s:%c  // \n",j,mpdb->pdb_file[I],ChainCharPDB(C,mpdb->pdb[I]));
		} else {
		  fprintf(vfp,"File%d=%s:%c  // \n",j,mpdb->pdb_file[I],ChainCharPDB(C,mpdb->pdb[I]));
		}
	}
	if(vsifp) fprintf(vsifp,"\n1-10000.W15\n"); else fprintf(vfp,"\n1-10000.W15\n");
        for(n=1; esc->PDB_SetI[S][n]; n++){
	  I=esc->PDB_SetI[S][n]; C=esc->PDB_SetC[S][n];
#if 0	// DEBUG
	if(n==1){
	   fprintf(vfp,"\n# NumRpts[%d][%d] = %d\n",I,C,NumRpts[I][C]);
	   fprintf(vfp,"\n# NumFullRpts[%d] = %d\n",S,NumFullRpts[S]);
	}
#endif
	  for(R=1; R<=NumRpts[I][C]; R++){
	     // Find adjacent, hetero subunits...
	     char **AdjMolecule=AdjacentHeteroMolecule(R,C,I,maxdist);
	     for(j=1; AdjMolecule[j]; j++){
		if(AdjMolecule[j][0] == '!'){	// indicates a single atom (ion).
		   if(vsifp) fprintf(vsifp,"\n%d%s.{X}\n",n,AdjMolecule[j]); 
		   else fprintf(vfp,"\n%d%s.{X}\n",n,AdjMolecule[j]); 
		} else {
		   if(vsifp) fprintf(vsifp,"\n%d!%s.C\n",n,AdjMolecule[j]); 
		   else fprintf(vfp,"\n%d!%s.C\n",n,AdjMolecule[j]); 
		} free(AdjMolecule[j]);
	     } free(AdjMolecule);
	  }
	} 
	if(vsifp) fprintf(vsifp,"\n"); else fprintf(vfp,"\n");
        // for(x=ptrn->NumPttrns; x > 0; x--)
#if 1
	BooLean	HptIsTree=FALSE;
	Int4 *Parent=0;
	if(hpt->IsTree(Parent)) HptIsTree=TRUE; 
	if(Parent) free(Parent);
#endif
	for(R=1; R <= NumFullRpts[S]; R++){
	  Int4 Row = RptCategory[S][R]; // hpt row...
	  assert(MemberSet(Row,RelevantSet[S]));
	  // if(!MemberSet(A,RelevantSet[S])) continue;	// skip irrelevant files.
	  Int4 X;
          // for(X=0,x=1; x <= ptrn->NumPttrns; x++) 
          for(X=0,x=ptrn->NumPttrns;  x > 0; x--)	// stored backwards...
	  {
	    Int4 Col=ptrn->PttrnCategory[x];
	    if(hpt->Cell(Row, Col) != '+') continue;
	    else X++;
	    // print out subgroups before supergroups to ensure proper color.
	    if(X > 40) continue;  // ran out of colors.
            BooLean first=TRUE;
#if 1	// add more info to vsi files; need to fix this so that only one trace is printed
	  if(R < MaxColor){
	    if(HptIsTree){
	       // fprintf(vfp,"\n# %d.%s:\n",Row,hpt->ElmntSetName(Row));
	       if(vsifp) fprintf(vsifp,"\n# %d.%s:\n",Col,hpt->ElmntSetName(Col)); // 
	       else fprintf(vfp,"\n# %d.%s:\n",Col,hpt->ElmntSetName(Col)); // 
	    } else if(vsifp) fprintf(vsifp,"\n# %d.%s:\n",Col,hpt->GrpName(Col));
	    else fprintf(vfp,"\n# %d.%s:\n",Col,hpt->GrpName(Col));
	    Int4 start=1,TheEnd=NumCol();
	    for(j=1; Col2FullSeq[S][R][j] == 0 && j <= TheEnd; j++){ } start=j;
	    for(j=TheEnd; Col2FullSeq[S][R][j] == 0 && j > 0; j--){ } TheEnd=j;
	    if(vsifp) fprintf(vsifp,"%d-%d.%c80\n",Col2FullSeq[S][R][start],
			Col2FullSeq[S][R][TheEnd],trace_color[R]);
	    else fprintf(vfp,"%d-%d.%c80\n",Col2FullSeq[S][R][start],
			Col2FullSeq[S][R][TheEnd],trace_color[R]);
	  }
#endif
	    for(p=1; p <= ptrn->NumPttrnRes[x]; p++){
	        // if(p > 10) continue;  // skip less significant.
		Int4 col=ptrn->PosPttrn[x][p];
		Int4 site=Col2FullSeq[S][R][col];
		if(site == 0) continue;  // not visible within structures.
		assert(site <= LenSeq(FullSeq[S]));
		Int4 r=ResSeq(site,FullSeq[S]);
		char Res=AlphaChar(r,AB);
// PutSeq(stderr,FullSeq[S],AB); 
// site+=OffSetSeq(FullSeq[S]);
		if(strchr(ptrn->PttrnRes[x][p],Res)){
		   // Int4 X=ptrn->NumPttrns-x+1; 
		   if(first){ 
			first=FALSE; 
			if(vsifp) fprintf(vsifp,"%c%d.%c",Res,site,side_color[X]); 
			else fprintf(vfp,"%c%d.%c",Res,site,side_color[X]);
		   } else {
			if(vsifp) fprintf(vsifp,",%c%d.%c",Res,site,side_color[X]);
			else fprintf(vfp,",%c%d.%c",Res,site,side_color[X]);
		   }
		} else {	// print mismatches in white.
		   // See whether shows up at a higher level:
		   BooLean IsHigher=FALSE;
		   for(Int4 x2=x+1;  x2 > 0; x2--){
	    	      Int4 Col2=ptrn->PttrnCategory[x2];
	    	      if(hpt->Cell(Row, Col2) != '+') continue;
	    	      for(Int4 p2=1; p2 <= ptrn->NumPttrnRes[x2]; p2++){
			Int4 col2=ptrn->PosPttrn[x2][p2];
			if(col == col2){ IsHigher=TRUE; break; }
		      } if(IsHigher) break;
		   } if(IsHigher) continue;
		   // if no higher level pattern, then print out.
		   if(first){
		      first=FALSE; 
		      if(vsifp) fprintf(vsifp,"%c%d.W",Res,site); 
		      else fprintf(vfp,"%c%d.W",Res,site); 
		   } else if(vsifp) fprintf(vsifp,",%c%d.W",Res,site);
		   else fprintf(vfp,",%c%d.W",Res,site);
	       	} 
	    } if(vsifp) fprintf(vsifp,"\n"); else fprintf(vfp,"\n");
#if 1	// print column positions as well...
	    for(first=TRUE,p=1; p <= ptrn->NumPttrnRes[x]; p++){
	        // if(p > 10) continue;  // skip less significant.
		Int4 col=ptrn->PosPttrn[x][p];
		Int4 site=Col2FullSeq[S][R][col];
		if(site == 0) continue;  // not visible within structures.
		assert(site <= LenSeq(FullSeq[S]));
		Int4 r=ResSeq(site,FullSeq[S]);
		char Res=AlphaChar(r,AB);
		if(strchr(ptrn->PttrnRes[x][p],Res)){
		   if(first){
			first=FALSE; 
			if(vsifp) fprintf(vsifp,"#%s%d",ptrn->PttrnRes[x][p],col); 
			else fprintf(vfp,"#%s%d",ptrn->PttrnRes[x][p],col);
		   } else if(vsifp) fprintf(vsifp,",%s%d",ptrn->PttrnRes[x][p],col);
		   else fprintf(vfp,",%s%d",ptrn->PttrnRes[x][p],col);
		} else {	// print mismatches in white.
		   if(first){
			first=FALSE; 
			if(vsifp) fprintf(vsifp,"#(%s%d)",ptrn->PttrnRes[x][p],col);
			else fprintf(vfp,"#(%s%d)",ptrn->PttrnRes[x][p],col); 
		   } else if(vsifp) fprintf(vsifp,",(%s%d)",ptrn->PttrnRes[x][p],col);
		   else fprintf(vfp,",(%s%d)",ptrn->PttrnRes[x][p],col);
	       	} 
	    } if(vsifp) fprintf(vsifp,"\n"); else fprintf(vfp,"\n");
#endif
	  }
	} if(vsifp) fprintf(vsifp,"\n"); else fprintf(vfp,"\n"); // fclose(vfp);
   } if(vfp) fclose(vfp); if(vsifp) fflush(vsifp);
}

char	**p2c_typ::AdjacentHeteroMolecule(Int4 RR, Int4 CC, Int4 II, double maxdist)
// ,res_typ **ResALL_I,Int4 *num_resALL_I,pdb_typ P)
// Return string array indicating which molecules are are adjacent to chain C0
// sprintf(str,"[%s]%d%c",res_name,res,chain);
{
        Int4    C,a,aa,i,j,jj,x,p;
        atm_typ atm,aatm;
        char	**OutPut,str[100];
	Int4	NumMol=0,MaxNumMol=1000;

	assert(II <= mpdb->NumPDB && II > 0); pdb_typ P=mpdb->pdb[II]; assert(CC <= nChainsPDB(P) && CC > 0);
        NEWP(OutPut,MaxNumMol+3,char);

	for(x=ptrn->NumPttrns; x > 0; x--){
	  for(p=1; p <= ptrn->NumPttrnRes[x]; p++){
	    Int4 pos=ptrn->PosPttrn[x][p];
	    Int4 Res=Col2pdbSeq[II][CC][RR][pos] + OffSetSeq(mpdb->pdbSeq[II][CC]);
	    res_typ *ResP=mpdb->ResALL[II][CC];
	    for(jj=1; jj <= mpdb->num_resALL[II][CC]; jj++){
	      if(ResidueID(ResP[jj]) != Res) continue;  // find the pattern residue (Res).
	      for(aa=1; aa <= ResidueAtomNumber(ResP[jj]); aa++){
		aatm=AtomResidue(aa,ResP[jj]);
	        if(IsWaterAtom(aatm)) continue;
                for(C=1; C <= nChainsPDB(P); C++){
		  if(C == CC) continue; 
		  if(IsProteinPDB(C,P)) continue;
		  res_typ *ResH=mpdb->ResALL[II][C];
		  for(j=1; j <= mpdb->num_resALL[II][C]; j++){
		      for(a=1; a <= ResidueAtomNumber(ResH[j]); a++){
			atm = AtomResidue(a,ResH[j]);
	     		if(IsWaterAtom(atm)) continue;
             		if(!IsHeteroAtom(atm)) break;	// all residue atoms are either hetero or not.
                        if(DistanceAtoms(atm, aatm) <= maxdist){
			   char *mole_name=AtomResName(atm);
			   while(isspace(mole_name[0])){
			      mole_name++;
			      if(!isprint(mole_name[0])) print_error("AdjacentHeteroMolecule() atom error");
			   }
			   if(AtomsInResidue(ResH[j]) == 1) {
			     sprintf(str,"![%s]%d%c",mole_name,ResAtom(atm),AtomChain(atm));
			   } else {
			     sprintf(str,"[%s]%d%c",mole_name,ResAtom(atm),AtomChain(atm));
			   }
			   // sprintf(str,"[%s]%d%c",AtomResName(atm),ResAtom(atm),AtomChain(atm));
			   // ResidueName(R); ResidueID(R); chn= ResidueChain(R); ChainCharPDB(chn,P);
			   BooLean is_new=TRUE;
			   for(i=1; i <= NumMol; i++){
				if(strcmp(str,OutPut[i]) == 0){ is_new=FALSE; break; }
			   }
			   if(is_new){
				NumMol++; OutPut[NumMol]=AllocString(str);  
				if(NumMol >= MaxNumMol) return OutPut;
			   } break;	// No need to look further...
			}
		      }
		  }
		}
	      }
	    }
	  }
        } return OutPut;
}

