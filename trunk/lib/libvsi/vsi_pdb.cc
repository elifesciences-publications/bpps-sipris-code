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

#if 0
#include "pdb.h"
#include "vsi_typ.h"
#include "histogram.h"
#include "set_typ.h"
#else
#include "vsi_pdb.h"
#endif

BooLean BurleyAromatic(FILE *fp,Int4 file_id, double *Energy,res_typ ResJ, 
		res_typ ResK, char chainJ,char chainK, char colorJ, char colorK,
		pdb_typ P);

#define MAX_angleDHA 90

static BooLean	SkipBB(register Int4 i,register char chain,register set_typ **skip)
{ if(MemberSet(i,skip[0][chain])) return TRUE; else return FALSE; }

static BooLean	SkipSC(register Int4 i,register char chain,register set_typ **skip)
{ if(MemberSet(i,skip[1][chain])) return TRUE; else return FALSE; }

Int4	FindHBondsPDB2(FILE *fp,res_typ Res,Int4 C,float dmax, 
		unsigned short file_id, set_typ **skip, char **Color,pdb_typ P,
		char mode)
// For all C-H, O-H, N-H, etc. bonds in Res search the entire pdb file for atoms
// that can accept the hydrogen.  Ignore any atoms in skip[BB][chainK][rJ]. 
{
	Int4	n,r,r2,rK,gly,rJ,C2,atmJ,atmK;
	Int4	BB=0,SC=1;	// backbone vs sidechain...
	atm_typ	H,Donor,Accept;
	char	colorJ,colorK,chainJ,chainK,*cp;
	BooLean	IsResD,IsResA,skipDSC,skipASC,skipDBB,skipABB;
	
	if(C > P->nchains || C<1) pdb_error("FindBondsPDB( ): input error"); 
        chainJ=ChainCharPDB(C,P);
	gly = AlphaCode('G',P->A);
	for(n=0,atmJ = 1; atmJ <= ResidueAtomNumber(Res); atmJ++){
	     if(!ResidueAtomHydrogens(atmJ,Res)) continue;
	     Donor = AtomResidue(atmJ,Res); rJ=ResAtom(Donor);
#if 0	// DEBUG:
	     // if(ResAtom(Donor) == 169 && AtomChain(Donor) == 'B' && NitrogenAtom(Donor)){
	     if(AtomID(Donor) == 481 || AtomID(Donor) == 482){
		PutAtom(stderr, Donor);
	     }
#endif
	     IsResD=IsAminoAcidAtom(Donor);
	     skipDSC=SkipSC(rJ,chainJ,skip); skipDBB=SkipBB(rJ,chainJ,skip);
	     if(IsResD && SideAtom(Donor) && skipDSC) continue;
	     if(skipDSC && skipDBB) continue;
	     // if made it this far assume backbone to be examined.
	     if(P->maxres[C] < rJ) r = 0;
	     else { if(P->seq[C]) r = P->seq[C][rJ];  else r = 0; }
	     // ******************** Find all Donor atoms: *****************************
	     for(Int4 h=1; (H=ResidueAtomHydrogen(atmJ,h,Res)); h++){
	      for(C2=1; C2 <= P->nchains; C2++){
		chainK=ChainCharPDB(C2,P);
		// ******************** Srch all Acceptor atoms ***********************
		for(atmK = 1; atmK <= P->natoms[C2]; atmK++){
	          Accept = P->atom[C2][atmK];
		  if(Accept == Donor) continue;
#if 1		// Ignore Nucleic acid internal H-bonds
	          if(IsNucleicAcidAtom(Donor) && IsNucleicAcidAtom(Accept) 
					&& AtomChain(Donor) == AtomChain(Accept)) continue;
#endif
		  if(HydrogenAtom(Accept)) continue;
		  IsResA=IsAminoAcidAtom(Accept);
		  if(IsResD && IsResA){
		        // Ignore backbone-only interactions.
			if(!SideAtom(Donor) && !SideAtom(Accept)) continue;
		  }
		  rK= ResAtom(Accept);
		  skipASC=SkipSC(rK,chainK,skip); skipABB=SkipBB(rK,chainK,skip);
		  if(skipASC && skipABB) continue;
		  if(skipASC && IsResA && SideAtom(Accept)) continue;
		  if(skipABB && IsResA && !SideAtom(Accept)) continue;
		  if(CarbonAtom(Donor) && CarbonAtom(Accept)) continue;
		  if(mode == 'W' && (CarbonAtom(Donor) || CarbonAtom(Accept))) continue;
	     	  if(P->seq[C2] && IsResA) r2 = P->seq[C2][rK];  else r2=0;
		  if(C==C2 && rJ == rK) continue;
		  // if(C==C2 && abs(rJ - rK) < 2) continue;
		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MAX_angleDHA){
		    atm_typ	A,B;
		    Int4 rA,rB;
		    unsigned char rcA;
		    Int4 isaccept=0;
		    A=Donor;B=Accept; rcA=r; 
		    BooLean NeitherCarbons=TRUE;
		    if(CarbonAtom(A) || CarbonAtom(B)) NeitherCarbons=FALSE;
		    // Loop for printing out both donor and acceptor atoms:
		    for(rcA=r; isaccept < 2; A=Accept,B=Donor,rcA=r2,isaccept++){
			rA=ResAtom(A); rB=ResAtom(B);
#if 1	// fix problem with negative residues numbers.  AFN: 8_26_2015.
	// if(rA < 0 || rB < 0) continue;
	if(rA < 0) continue;
#endif
			char str[20],*str_ptr = AtomName(A);
			if(IsWaterAtom(A) && A==Donor) str_ptr = AtomName(H);
			if(isspace(str_ptr[0])) str_ptr++;
			Int4 z,zz;
			for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
			} str[z]=0;
			if(IsWaterAtom(A)){
			  if(A==Donor){
			    if(AtomChain(A)) fprintf(fp,"HOH%d%c_o-%s.X\t// to ",rA,AtomChain(A),str);
			    else fprintf(fp,"HOH%d_o-%s.X\t// to ",rA,str);
			  } else {
			    if(AtomChain(A)) fprintf(fp,"HOH%d%c_o.X\t// to ",rA,AtomChain(A));
			    else fprintf(fp,"HOH%d_o.X\t// to ",rA);
			  }
			} else {
			  // add hydrogen too...
			  if(!isaccept){	// i.e. donor atom...
			    str_ptr = AtomName(H);
			    if(isspace(str_ptr[0])) str_ptr++; str[z]='-'; z++;
				
			    for(zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
			    } str[z]=0;
			  }
			 if(strcmp("o",str)==0) strcpy(str,"c-o");
			 else if(str[0]=='c' && strstr(str,"-") == 0){
				 fprintf(fp,"// ");
			 }
			 if(NeitherCarbons){
				if(file_id > 0) fprintf(fp," %d",file_id);
				else fprintf(fp," ");
			 } // else fprintf(fp,"#"); 
			 else fprintf(fp," ");
			 if(IsHeteroAtom(A)){
				fprintf(fp,"!");
			 }
			 if(rcA){
			    if(AtomChain(A) != ' '){
			        fprintf(fp,"%c%d%c",AlphaChar(rcA,P->A),rA,AtomChain(A));
			    } else fprintf(fp,"%c%d", AlphaChar(rcA,P->A),rA);
			 } else {
			    cp=AtomResName(A); while(cp[0]==' ') cp++;
			    Int4 end = strlen(cp); end--;
			    assert(end >= 0);
#if 1
			    BooLean Digit=FALSE;
			    for(char c=0; c <= end; c++) if(isdigit(cp[c])){ Digit=TRUE; break; }
			    if(Digit){
			     if(AtomChain(A) != ' '){
			  	fprintf(fp,"[%s]%d%c",cp,rA,AtomChain(A));
			     } else fprintf(fp,"[%s]%dX",cp,rA);
			    }
#else
			    if(isdigit(cp[end]) || isdigit(cp[0])){
			     if(AtomChain(A) != ' '){
			  	fprintf(fp,"[%s]%d%c",cp,rA,AtomChain(A));
			     } else fprintf(fp,"[%s]%dX",cp,rA);
			    }
#endif
			    else {
			     if(AtomChain(A) != ' '){
			  	fprintf(fp,"%s%d%c",cp,rA,AtomChain(A));
			     } else fprintf(fp,"%s%dX",cp,rA);
			    }
			 } fprintf(fp,"_%s.X\t// to ",str);
			}
			if(AtomChain(B) != ' '){
			  fprintf(fp,"%s %s%d%c (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).",
					AtomName(B),AtomResName(B),
					rB,AtomChain(B),d_DA,d_HA,angle_DHA);
			} else fprintf(fp,"%s %s%d (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).",
					AtomName(B),AtomResName(B),rB,d_DA,d_HA,angle_DHA);
			if(NeitherCarbons) fprintf(fp,"\n"); else fprintf(fp,"@\n");
			// if one of the atoms is a carbon then label the line with '@' symbol.
		    } n++;
		  } 
		 }
		}
	     }
	} return n;
}

void	PutContactsPDB2(FILE *fptr, Int4 file_id, float HA_dmax, float dmax,
		set_typ **skip,char **Color,pdb_typ P, char mode)
{
	Int4	res0,C2,C;
	Int4	rJ,rK,i,j,k,num_resA,num_resD;
	Int4	BB=0,SC=1;	// sidechain and backbone indices...
	char	color=0,chainJ,chainK;
	res_typ *ResCA,*ResOA,*ResCD,*ResOD,**ResALL;
	Int4	num_resC,num_resO,*num_resALL;
	res_typ ResK,ResJ;
	char	colorJ,colorK;

	// if(chain==0) chain='?'; // this will work in rasmol for unlabeled chains!

	NEWP(ResALL,nChainsPDB(P)+3,res_typ);
	NEW(num_resALL,nChainsPDB(P)+3,Int4);
	for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	  ResALL[C2] = MakeResPDB(C2,&num_resALL[C2],P);
	}

	// Move this routine to libpdb.a at some point:
	a_type AB = AminoAcidAlphabetPDB(P);
        // h_type	HG=Histogram("Aro-Aro energy",-100,100,1.0);

for(C=1; C <= nChainsPDB(P); C++){
        chainJ=ChainCharPDB(C,P);
	num_resC=num_resALL[C];
	ResCA=ResCD=ResALL[C];
   	//********************* Srch for ARo-Aro...****************************
        for(j=1; j <= num_resC; j++){
          ResJ=ResALL[C][j];
          rJ=ResidueID(ResJ);
          if(SkipSC(rJ,chainJ,skip)) continue;
	  colorJ=Color[chainJ][rJ];
	  for(C2=C; C2 <= nChainsPDB(P); C2++){
             chainK=ChainCharPDB(C,P);
	     num_resO=num_resALL[C2];
	     if(C2==C) k=j+1; else k=1;
             for( ; k <= num_resO; k++){
                ResK=ResALL[C2][k];
                rK=ResidueID(ResK);
                if(SkipSC(rK,chainK,skip)) continue;
	        colorK=Color[chainK][rK];
                double E=0.0;
                if(BurleyAromatic(fptr,file_id,&E,ResJ,ResK,chainJ,chainK,colorJ,colorK,P)){
                 if(E < 0){
                   // fprintf(stderr,"E = %g\n\n",E);
                 } // else fprintf(stderr,"failed: E = %g\n\n",E);
                 // IncdHist(E,HG);
                }
	     }
          }
        }
   	//********************* End Srch for ARo-Aro...****************************
	//******************* Srch for H-bonds... *************************
        for(j=1; j <= num_resC; j++){	
          ResJ=ResALL[C][j]; rJ=ResidueID(ResJ);
#if 0	// DEBUG
	  if(rJ == 169) PrintResidueAtoms(stderr,ResJ);
#endif
	  if(SkipBB(rJ,chainJ,skip) && SkipSC(rJ,chainJ,skip)) continue;
	  // fprintf(fptr,"\n"); // Note: ResJ == H-bond donor 
	  if(mode != 'A') FindHBondsPDB2(fptr,ResJ,C,HA_dmax,file_id,skip,Color,P,mode); 
	}

	//=========== pi-bonds: ===========
	if(mode != 'A') for(j=1; j <= num_resC; j++){
           ResJ=ResALL[C][j]; rJ=ResidueID(ResJ);  // X-H donor residue...
	   if(SkipBB(rJ,chainJ,skip) && SkipSC(rJ,chainJ,skip)) continue;
	   // Need to also allow donor backbone atoms: don't now...
           for(C2=1; C2 <= nChainsPDB(P); C2++){
              chainK=ChainCharPDB(C2,P);
              ResOD=ResOA=ResALL[C2];   // other than query residue.
              num_resO=num_resALL[C2];
	      for(k=1; k <= num_resO; k++){
           	ResK=ResALL[C2][k]; rK=ResidueID(ResK);  // pi-orbital residue...
		if(C==C2 && rJ==rK) continue;
	        if(SkipBB(rK,chainK,skip) && SkipSC(rK,chainK,skip)) continue;
		if(!SkipSC(rK,chainK,skip)){	// pi-orbitals == K
	   	  if(!(SkipBB(rJ,chainJ,skip) || SkipSC(rJ,chainJ,skip))){
		   FindAromaticHbondsPDB(fptr,ResCD[j],ResOA[k],C,C2,dmax,color,P);
		   FindOtherPiHbondsPDB(fptr,ResCD[j],ResOA[k],C,C2,dmax,color,P);
		  }
		}
                if(SkipSC(rJ,chainJ,skip)){	// pi-orbitals == J
	          if(!(SkipBB(rK,chainK,skip) || SkipSC(rK,chainK,skip))){
		   FindAromaticHbondsPDB(fptr,ResOD[k],ResCA[j],C2,C,dmax,color,P);
		   FindOtherPiHbondsPDB(fptr,ResOD[k],ResCA[j],C2,C,dmax,color,P);
		  }
		}
	      }
	    }
	 } //********** end of res0 loop... *****
   }//************************** End of loop over all chains C *******************
        // PutHist(stderr,60,HG); NilHist(HG);
	 for(C2=1; C2 <= nChainsPDB(P); C2++){ 
		for(i=1; i <= num_resALL[C2]; i++) NilRes(ResALL[C2][i]);
		free(ResALL[C2]);
	 } free(ResALL); free(num_resALL);
} 

vsi_typ	PutAutoFormatVSI(FILE *fp,Int4 file_id, vsi_typ HEAD_NODE,
		char KEY_CHAIN,float HA_dmax,float dmax, BooLean only_if_trace,
		pdb_typ P, char mode)
// Automatically Formating a vsi file by adding hydrogen bonds, Aro-Aro and 
// other interactions.
{
	char	chain,str[200],*str2,*DATADIR=0;
	vsi_typ node,last;
	item_type itm;

       	Int4    C=GetChainNumberPDB(P,KEY_CHAIN);
       	Int4    resStart=MinResPDB(C,P);
       	Int4    resEnd=MaxResPDB(C,P),r,C2;

	char	**Color=0;
	set_typ *skip[3];
	Int4	x,SC=1,BB=0;	// sidechain vs backbone
	NEWP(Color,300,char);	// use printable characters
	NEW(skip[BB],300,set_typ);
	NEW(skip[SC],300,set_typ);
	x=ReLabelPDB(P);
	if(x > 0) fprintf(stderr,"%d chains relabeled\n",x);
	for(C2=1; C2 <= nChainsPDB(P); C2++){
	   chain=ChainCharPDB(C2,P);
	   if(chain==' '){
		fprintf(stderr,"chain='%c'; C2=%d\n",chain,C2);
		fprintf(stderr,"Possibly fatal semantic error: unlabeled chain\n");
		// print_error("Fatal semantic error: unlabeled chain");
	   }
	   if(Color[chain] == 0){	// don't reallocate redundant chain char
	     NEW(Color[chain],MAX_RESIDUES_PDB+2,char);
	     skip[BB][chain]=MakeSet(MAX_RESIDUES_PDB+2); // pdb maximum is 9999
	     skip[SC][chain]=MakeSet(MAX_RESIDUES_PDB+2); // pdb maximum is 9999
	     FillSet(skip[BB][chain]);	// set all members to 'TRUE'
	     FillSet(skip[SC][chain]);	// set all members to 'TRUE'
	   } else {
#if 1
		fprintf(stderr,"pdb filename = '%s'\n",FilenamePDB(P));
		// fprintf(stderr,"Color[%c] = '%c'(%d)\n",chain,Color[chain],Color[chain]);
		for(Int4 c2=1; c2 <= nChainsPDB(P); c2++){
			char chn=ChainCharPDB(c2,P);
			fprintf(stderr,"Chain[%d] = '%c'\n",c2,chn);
		} print_error("Fatal sematic error: redundant chain designations");
#else
		fprintf(stderr,"WARNING: sematic error: redundant chain designations\n");
#endif
	   }
	}
	for(node=HEAD_NODE; node; node=node->next){
	    node->file=file_id;
	    if(node->type == 'T'){
		if(node->item.trace.chain == 0)
			node->item.trace.chain=KEY_CHAIN;
	        itm=node->item; chain=itm.trace.chain;
		if(skip[BB][chain] == 0){
			PutTraceItem(stderr,itm,file_id);
			fprintf(stderr,"File name = %s\n",FilenamePDB(P));
			print_error("\nFatal semantic error: non-existent chain 1.");
		}
		assert(itm.trace.end <= MAX_RESIDUES_PDB);
		for(r=itm.trace.start; r <= itm.trace.end; r++){
		   DeleteSet(r,skip[BB][chain]);
		}
	    } else if(node->type == 'R'){
	      if(node->item.residue.chain==0) node->item.residue.chain=KEY_CHAIN; 
	      if(skip[BB][node->item.residue.chain] == 0){
			PutResItem(stderr,node->item,0,file_id);
			fprintf(stderr,"File name = %s\n",FilenamePDB(P));
			print_error("\nFatal semantic error: non-existent chain 2.");
	      }
	    } else if(node->type == 'M'){	// trace...
	      if(node->item.molecule.chain==0) node->item.molecule.chain=KEY_CHAIN;
	      if(skip[BB][node->item.molecule.chain] == 0){
			PutMolItem(stderr,node->item,file_id);
			fprintf(stderr,"File name = %s\n",FilenamePDB(P));
			print_error("\nFatal semantic error: non-existent chain 3.");
	      }
	    }
	}
	for(last=0,node=HEAD_NODE; node; node=node->next){
	    if(node->type == 'R'){
	      itm=node->item; r=itm.residue.site; chain=itm.residue.chain;
	      // fprintf(stderr,"res=%c%d%c:",itm.residue.res,r,chain);
	      if(only_if_trace){
		if(last==0){
	      	   fprintf(stderr,"residue=%c%d%c:\n",itm.residue.res,r,chain);
		   fprintf(stderr,
		     "FATAL ERROR: a region of chain %c needs to be specified.\n",chain);
		   fprintf(stderr,
		     "The syntax for this is <res_start>-<res_end><chain>.<color><width>.\n");
		   fprintf(stderr,
		     "Example: '1-250A.B50' where 'A' and 'B' imply chain A and color blue.\n");
		   exit(1);
		} assert(last);	// a residue should not be first in this file.
		if(MemberSet(r,skip[BB][chain])){  // then delete this node...
		   last->next=node->next; FreeNode(node); node=last;
		   // fprintf(stderr,"deleted\n");
		} else {
		   DeleteSet(r,skip[SC][chain]);
		   Color[chain][r]=itm.residue.color; 
		   // fprintf(stderr,"retained\n");
		}
	      } else {
		DeleteSet(r,skip[SC][chain]); DeleteSet(r,skip[BB][chain]); 
		Color[chain][r]=itm.residue.color; 
	      }
	    } else if(node->type == 'M'){	// trace...
	        itm=node->item; chain=itm.molecule.chain; r=itm.molecule.id;
		DeleteSet(r,skip[SC][chain]); DeleteSet(r,skip[BB][chain]);
		Color[chain][r]=itm.molecule.color;
	    } last=node;
	}

	PutContactsPDB2(fp,file_id,HA_dmax,dmax,skip,Color,P,mode);
	fprintf(fp,"\n"); 

        for(C2=1; C2 <= nChainsPDB(P); C2++){
            chain=ChainCharPDB(C2,P);
            if(Color[chain]){ free(Color[chain]); Color[chain]=0; }
            if(skip[SC][chain]) NilSet(skip[SC][chain]);
            if(skip[BB][chain]) NilSet(skip[BB][chain]);
        } free(skip[SC]); free(skip[BB]); free(Color);
	return HEAD_NODE;
}


