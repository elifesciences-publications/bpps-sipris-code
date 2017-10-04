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

#include "mps_typ.h"

//*************************** mps_typ ****************************************
//*************************** mps_typ ****************************************
//*************************** mps_typ ****************************************

void	mps_typ::GetFileNamesPDB( )
//********************** 0. Read pdb input filenames (full path) ***************************
{
        char    str[505];
        Int4    s;

	NumPDB=0;
	NEWP(pdb_file,MAX_NUM_PDB_INPUT_FILES + 5,char);
        FILE *fp = open_file(pdb_infile,"","r");
        while(fgets(str,500,fp)){
                for(s=0; !isspace(str[s]); s++) ;
                if(s == 0){
                        if(NumPDB== 0) print_error("mps_typ GetFileNamesPDB() input error");
                        else break;     // stop here...
                }
                str[s]=0;
                NumPDB++;
                if(NumPDB > MAX_NUM_PDB_INPUT_FILES) print_error("too many pdb input files");
                pdb_file[NumPDB] = AllocString(str);
        } fclose(fp);
	if(NumPDB < 1) print_error("mps_typ: FATAL ERROR - no input pdb files found");
}

void	mps_typ::Free( )
{
	Int4 i,j,c;
	for(i=1; i <= NumPDB; i++){
	   for(c=1; c <= nChainsPDB(pdb[i]); c++){
	      if(ResALL[i][c]){
	      	for(j=1; j <= num_resALL[i][c]; j++){
	   	    if(ResALL[i][c][j]) NilRes(ResALL[i][c][j]);
	      	} free(ResALL[i][c]);
	      }
	      if(Calpha[i][c]) free(Calpha[i][c]);
	      if(pdbSeq[i][c]) NilSeq(pdbSeq[i][c]);
	   } free(pdbSeq[i]); free(ResALL[i]); NilPDB(pdb[i]); free(pdb_file[i]);
	   free(Calpha[i]); free(num_resALL[i]);
	} free(pdb_file); free(pdbSeq); free(ResALL); free(Calpha); free(num_resALL); free(pdb);
	free(pdb_infile);
}

void	mps_typ::Read( )
//********************** 1. Read pdb and cma input files ***************************
{
	FILE *fp = 0;
	
	NEW(pdb,MAX_NUM_PDB_INPUT_FILES + 5,pdb_typ);
        NEWP(pdbSeq,MAX_NUM_PDB_INPUT_FILES + 5,e_type);
	for(Int4 I=1; I <=NumPDB; I++){	
	   fp = open_file(pdb_file[I],"","r");
	   pdb[I]=MakePDB(fp,pdb_file[I]); fclose(fp);
	   if(pdb[I] == 0) print_error("pdb file read error");
	   else {	// clean up files
		   BooLean okay=TRUE;
                   long natms=ChangePO3toAtmPDB(pdb[I]);
		   // fprintf(stderr,"------- %s ---------.\n",pdb_file[I]);
                   if(natms > 0){
			okay=FALSE;
			fprintf(stderr,"%d phosphorylated amino acids changed to ATOM.\\n",natms);
		   }
		   natms=ChangeMSE2MET_PDB(pdb[I]);
		   if(natms > 0){
			okay=FALSE;
			fprintf(stderr,"%d selenomethionine residues changed to Met.\n",natms);
		   }
		   if(!okay){	// then redo pdb[I]...
			fp=tmpfile(); PutPDB(fp, pdb[I]); NilPDB(pdb[I]); 
			rewind(fp); pdb[I]=MakePDB(fp,pdb_file[I]); fclose(fp);
		   }
	   }
	   NEW(pdbSeq[I], nChainsPDB(pdb[I]) + 3, e_type);
	   for(Int4 ch=1; ch <= nChainsPDB(pdb[I]); ch++){
		if(!IsFullProteinChainPDB(ch,pdb[I])){ pdbSeq[I][ch]=0; continue; }
		e_type E=GetPDBSeq(ch,pdb[I]);
		ShortenSeqInfo(E);	// eliminate path information.
		if(E != NULL){
			Int4 numX = NumXSeq(E);
			double frqX = (double) numX / (double) LenSeq(E);
			if(frqX >= 0.25){
			  fprintf(stderr,"%s: frqX = %.3f (skipped).\n",FilenamePDB(pdb[I]),frqX);
			  NilSeq(E); pdbSeq[I][ch]=0;
			} else {
			  pdbSeq[I][ch]=E; // PutSeq(stdout,pdbSeq[I][ch],AB); 
			}
		} else { pdbSeq[I][ch]=NULL; }
	   } // PutPDB(stdout, pdb[I]);
	}
}

void	mps_typ::GetResidues()
//********************** 2. Get residues and Calphas for pdb seqs *********************
{
   atm_typ atm=0;
   res_typ R;
   Int4    n,r,r1,r2,a,i;
   char c,c1,c2;
   a_type  AB=AminoAcidAlphabetPDB(pdb[1]);

   NEWPP(Calpha,MAX_NUM_PDB_INPUT_FILES + 5,atm_typ);
   NEWPP(ResALL,MAX_NUM_PDB_INPUT_FILES + 5,res_typ);
   NEWP(num_resALL,MAX_NUM_PDB_INPUT_FILES + 5,Int4);
   for(Int4 I=1; I <=NumPDB; I++){
        NEWP(ResALL[I],nChainsPDB(pdb[I])+3,res_typ);
        NEW(num_resALL[I],nChainsPDB(pdb[I])+3,Int4);
        NEWP(Calpha[I],nChainsPDB(pdb[I])+3,atm_typ);
        for(Int4 C=1; C <= nChainsPDB(pdb[I]); C++){
          ResALL[I][C] = MakeResPDB(C,&num_resALL[I][C],pdb[I]);
	  e_type E=pdbSeq[I][C];
	  if(E){
	   Int4 array_len =MAXIMUM(Int4,LenSeq(E),num_resALL[I][C]);
           NEW(Calpha[I][C],array_len+3,atm_typ); // may have more residues than LenSeq!
	   for(r=1; r <= num_resALL[I][C]; r++) {
	     R = ResALL[I][C][r];
	     i=ResidueID(ResALL[I][C][r])-OffSetSeq(E);
	     if(i <= 0 || i > LenSeq(E)){
#if 0
		PutSeq(stderr,E,AB);
	     	PrintResidueAtoms(stderr,R);
	        assert(!(i <= 0 || i > LenSeq(E)));
#else
		continue;	// skip over these...
#endif
	     }
	     r1 = ResSeq(i,E); c1 = AlphaChar(r1,AB);
	     c2 = GetCharResidue(R,AB);
	     // printf("%c%d (pdb)  %c%d (Res)\n",c1,i+OffSetSeq(E),c2,ResidueID(R));
	     if(c1 != c2){
#if 0
		PutSeq(stderr,E,AB);
	     	PrintResidueAtoms(stderr,R);
	     	assert(c1 == c2);
#else
		continue;	// skip over these... alternative positions...
#endif
	     }
	     for(a = 1; a <= ResidueAtomNumber(R); a++){
		atm = AtomResidue(a,R);
		if(UseBeta){
		  if(c1 == 'G'){
		    if(AlphaCarbonAtom(atm)){
			Calpha[I][C][i] = atm; // PutAtom(stdout, atm);
			break;
		    }
		  } else if(BetaCarbonAtom(atm)){	// use beta carbons for non-glycine residues
		    Calpha[I][C][i] = atm; // PutAtom(stdout, atm);
		    break;
		  }
		} else if(AlphaCarbonAtom(atm)){
		  Calpha[I][C][i] = atm; // PutAtom(stdout, atm);
		  break;
		}
	     }
	   } 
	  }
 	}
   }
}

