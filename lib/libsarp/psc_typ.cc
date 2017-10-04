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

#define USAGE_START \
"USAGE: asepsia <pdb_name_file> <cmafile> [-options]\n\
   (\"automated structural evaluation of protein sequences in alignments\")\n\
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
    -RE=<float>                 Set the minimum column RelEntropy to consider (default: 0.0)\n\
    -seed=<int1>                seed for random number generator\n\
    -seqdist=<int>              Set the minimum distance between residues to consider (default: 4)\n\
    -maxsqdist=<int>            Set the maximum distance between residues to consider (default: 200) \n\
    -mindist=<int>              Set the minimum distance between aligned columns to compare (default: 3)\n\
    -minvar=<int>               Set the minimum variance to consider (default: 0)\n\
\n\n"

void	psc_typ::initialize()
// set parameters to default values.
{
	verbose=TRUE;
	RelevantSeq=0;	// This is set within MapSqAlnToStruct( );
	Int4	i;
	time1=time(NULL);
	AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	nAB = MkAlpha("NACGTU",DNA_MTRX);
      	HA_dmax=0.0, dmax=0.0;
	pttrn_filename=0; ptrn=0;
	NumWeakHBonds=NumModHBonds=0;
	cmafilename = 0;
	mpdb0=mpdb=0;
	ptrn0=ptrn=0;
}

void	psc_typ::Free()
{
        Int4 S,R,I,C;
	for(I=1; I <=mpdb->NumPDB; I++){
     	   for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
     	      for(R=1; R <= NumRpts[I][C]; R++) free(matches[I][C][R]);
	      free(matches[I][C]);
	   } free(matches[I]);
	} free(matches);
        for(I=1; I <=mpdb->NumPDB; I++){
             for(C=1; C <=nChainsPDB(mpdb->pdb[I]); C++){
                if(PatternResHbonds[I][C]) free(PatternResHbonds[I][C]);
                if(WeakHBonds[I][C]) free(WeakHBonds[I][C]);
                if(ModHBonds[I][C]) free(ModHBonds[I][C]);
	     } free(PatternResHbonds[I]); free(WeakHBonds[I]); free(ModHBonds[I]);
        } free(PatternResHbonds); free(WeakHBonds); free(ModHBonds);
	if(ptrn0) delete ptrn0;
	if(mpdb0) delete mpdb0;
	double runtime=difftime(time(NULL),time1);
	free(cmafilename); NilAlpha(AB); NilAlpha(nAB);
	fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
}

void	psc_typ::GetArg(Int4 argc, char *argv[])
// *************** Get arguments for all program options **********************
{
	Int4 arg,i;
        if(argc < 3) print_error(USAGE_START);
        for(arg = 0; arg < argc; arg++){
		fprintf(stdout,"%s ",argv[arg]);
	} fprintf(stdout,"\n");
        for(arg = 3; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
	     case 'D':
                if(sscanf(argv[arg],"-D%f",&dmax) != 1 || dmax > 100 || dmax < 0){
                        print_error(USAGE_START);
                } break;
             case 'd':
                if(sscanf(argv[arg],"-d%f",&HA_dmax) != 1 || HA_dmax > 100 || HA_dmax < 0)
                        print_error(USAGE_START);
                break;
             default : print_error(USAGE_START);
           }
	 }
	}
}

#if 0  // moved to p2c_typ.cc 
char	**psc_typ::AdjacentHeteroMolecule(Int4 RR, Int4 CC, Int4 II, double maxdist)
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
#endif

void	psc_typ::FindMatches( )
{
   Int4 x,p,n,i,j,I,C,R;
   Int4 site,res,pos;
   NEWP3(matches,mpdb->NumPDB +3, BooLean);
   for(I=1; I <=mpdb->NumPDB; I++){
     	NEWPP(matches[I],nChainsPDB(mpdb->pdb[I]) + 3, BooLean);
     	for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
	   NEWP(matches[I][C],NumRpts[I][C]+ 3, BooLean);
     	   for(R=1; R <= NumRpts[I][C]; R++){
		NEW(matches[I][C][R],num_cols + 3, BooLean);
	   }
	}
   }
   for(i=1; i <= NumPDB_Sets; i++){
	if(CardSet(RelevantSet[i]) == 0) continue;	// skip irrelevant files.
        for(n=1; PDB_SetI[i][n]; n++){
	   I=PDB_SetI[i][n]; C=PDB_SetC[i][n];
           if(mpdb->pdbSeq[I][C] == 0) continue;
     	   for(R=1; R <= NumRpts[I][C]; R++){
            for(x=1; x <= ptrn->NumPttrns; x++){
	      for(p=1; p <= ptrn->NumPttrnRes[x]; p++){
		pos=ptrn->PosPttrn[x][p];
		assert(pos > 0 && pos <= num_cols);
		if(Col2pdbSeq[I][C][R][pos] == 0) continue;
		site = Col2pdbSeq[I][C][R][pos];
		res=ResSeq(site,mpdb->pdbSeq[I][C]);
		if(res == 0) continue;
		char Res=AlphaChar(res,AB);
		if(strchr(ptrn->PttrnRes[x][p],Res)){ matches[I][C][R][pos]=TRUE; }
	      }
	    }
	  }
	}
   }
}

#if 0	// moved to p2c_typ.cc
void	psc_typ::PrintVSI_Files( )
{
   //************* Print out vsi files for each set.
   Int4 A,x,p,n,j,I,C,R,S;
   Int4 vsi_number=0;
   char color[]="WMROYGCBDLMROYGCBDLMROYGCBDLWWWWWWWWWWWWWWWWWWW";
   double maxdist=4.0;
   FILE *vfp;

   for(S=1; S <= NumPDB_Sets; S++){
	char vsifile[200];
	if(!MemberSet(cmaFileID,RelevantSet[S])) continue;	// skip irrelevant files.
	assert(PDB_SetI[S]);
	assert(cmafilename);
	vsi_number++;
	sprintf(vsifile,"_pdb%d.vsi",vsi_number);
	vfp=open_file(cmafilename,vsifile,"w");

        for(j=1; PDB_SetI[S][j]; j++){
		I=PDB_SetI[S][j]; C=PDB_SetC[S][j];
		fprintf(vfp,"File%d=%s:%c  // \n",j,mpdb->pdb_file[I],ChainCharPDB(C,mpdb->pdb[I]));
	}
	fprintf(vfp,"\n1-10000.W30\n");
        for(n=1; PDB_SetI[S][n]; n++){
	  I=PDB_SetI[S][n]; C=PDB_SetC[S][n];
	  for(Int4 R=1; R<=NumRpts[I][C]; R++){
	     // Find adjacent, hetero subunits...
	     char **AdjMolecule=AdjacentHeteroMolecule(R,C,I,maxdist);
	     for(j=1; AdjMolecule[j]; j++){
		if(AdjMolecule[j][0] == '!'){	// indicates a single atom (ion).
		   fprintf(vfp,"\n%d%s.{X}\n",n,AdjMolecule[j]); 
		} else {
		   fprintf(vfp,"\n%d!%s.C\n",n,AdjMolecule[j]); 
		}
	     }
	  }
	} fprintf(vfp,"\n");
        // for(x=ptrn->NumPttrns; x > 0; x--)
	for(R=1; R <= NumFullRpts[S]; R++){
	  // A = RptCategory[S][R]; 
	  // if(!MemberSet(A,RelevantSet[S])) continue;	// skip irrelevant files.
          for(x=1; x <= ptrn->NumPttrns; x++) {
	    // print out subgroups before supergroups to ensure proper color.
	    if(x > 40) continue;  // ran out of colors.
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
		   Int4 X=ptrn->NumPttrns-x+1; 
		   if(first){
			fprintf(vfp,"%c%d.%c",Res,site,color[X]); first=FALSE; }
		   else fprintf(vfp,",%c%d.%c",Res,site,color[X]);
		} else {	// print mismatches in white.
		   if(first){ fprintf(vfp,"%c%d.W",Res,site); first=FALSE; }
		   else fprintf(vfp,",%c%d.W",Res,site);
	       	} 
	     } fprintf(vfp,"\n");
	  }
	} fprintf(vfp,"\n"); fclose(vfp);
   }
}

#endif

