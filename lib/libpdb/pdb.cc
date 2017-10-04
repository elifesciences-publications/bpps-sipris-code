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

#include "pdb.h"

Int4	ReLabelPDB(pdb_typ P)
// Label all unlabeled or duplicately labeled chains
{
        BooLean *Used;
	Int4	C,NumChngd=0;
        NEW(Used,265,char);
        char ChainLabels[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	if(nChainsPDB(P) > 26) print_error("Fatal: Too many chains in pdb file to relabel");
        for(C=1; C <= nChainsPDB(P); C++){
                char chain=ChainCharPDB(C,P);
                if(chain != ' ' && Used[chain]){
                  fprintf(stderr,"duplicate chain %d: '%c'\n",C,chain);
                  chain = ' ';
                  if(!ReNameChainPDB(C,chain,P)) print_error("Fatal: ReLabelPDB() routine failure");
                  char chain=ChainCharPDB(C,P);
                  fprintf(stderr,"changed %d to: '%c'\n",C,chain);
                } else if(chain != ' '){
                  fprintf(stderr,"chain %d: '%c'\n",C,chain);
                  Used[chain]=TRUE;
                }
        }
        for(C=1; C <= nChainsPDB(P); C++){
                char chain=ChainCharPDB(C,P);
                if(chain == ' '){
		  NumChngd++;
                  Int4 j=0;
                  do { chain=ChainLabels[j]; j++; } while(chain != 0 && Used[chain]);
                  if(chain == 0) print_error("FATAL: too many chains");
                  fprintf(stderr,"chain %d: ' ' --> '%c'\n",C,chain);
                  if(!ReNameChainPDB(C,chain,P)) print_error("Fatal: ReLabelPDB() routine failure");
                  Used[chain]=TRUE;
                }
        } free(Used);
	return (NumChngd);
}

pdb_typ MkPDB(char *infile)
{
	FILE	*fptr = open_file(infile,".pdb","r");
	pdb_typ	P=MakePDB(fptr,infile);
	P->fp_stderr=stderr;
	fclose(fptr); return P;
}

void	ErrorToDevNullPDB(pdb_typ pdb) { pdb->fp_stderr=0; }

pdb_typ  MakePDB(char *infile)
{
	FILE	*fptr;
	fptr = open_file(infile,"","r");
	pdb_typ	P=MakePDB(fptr,infile);
	P->fp_stderr=stderr;
	fclose(fptr); return P;
}

pdb_typ  MakePDB(FILE *fptr) { return MakePDB(fptr,"unknown.pdb"); }

pdb_typ  MakePDB(FILE *fptr,char *infile)
{
	pdb_typ	P;
	Int4	C,r,max;

	NEW(P,1,protdb_type);
	P->filename = AllocString(infile);
	P->nchains = 0;
	P->fp_stderr=stderr;
	count_atoms(fptr,P);
	// fprintf(stderr,"nchains = %d\n",P->nchains);
	for(max = 0, C=1; C <= P->nchains; C++){
#if 0
	   fprintf(stderr,"natoms[%d] = %d; numres[%d] = %d\n",
		C,P->natoms[C],C,P->maxres[C]);
#endif
	   max = MAXIMUM(Int4,max,P->natoms[C]);
	}
	NEW(P->temp, max + 3, Int4);
	P->A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	/** Sander's total contacts: **/
	NEW(P->C_total, nAlpha(P->A) + 3, float);
	r = AlphaCode('X',P->A); P->C_total[r] = 0.0;
	r = AlphaCode('A',P->A); P->C_total[r] = 9.71;
	r = AlphaCode('R',P->A); P->C_total[r] = 70.28;
	r = AlphaCode('N',P->A); P->C_total[r] = 38.77;
	r = AlphaCode('D',P->A); P->C_total[r] = 41.76;
	r = AlphaCode('C',P->A); P->C_total[r] = 19.88;
	r = AlphaCode('Q',P->A); P->C_total[r] = 48.06;
	r = AlphaCode('E',P->A); P->C_total[r] = 53.19;
	r = AlphaCode('G',P->A); P->C_total[r] = 6.24;
	r = AlphaCode('H',P->A); P->C_total[r] = 65.86;
	r = AlphaCode('I',P->A); P->C_total[r] = 36.66;
	r = AlphaCode('L',P->A); P->C_total[r] = 36.73;
	r = AlphaCode('K',P->A); P->C_total[r] = 47.57;
	r = AlphaCode('M',P->A); P->C_total[r] = 40.40;
	r = AlphaCode('F',P->A); P->C_total[r] = 68.84;
	r = AlphaCode('P',P->A); P->C_total[r] = 29.99;
	r = AlphaCode('S',P->A); P->C_total[r] = 19.37;
	r = AlphaCode('T',P->A); P->C_total[r] = 28.74;
	r = AlphaCode('Y',P->A); P->C_total[r] = 77.80;
	r = AlphaCode('V',P->A); P->C_total[r] = 27.98;
	r = AlphaCode('W',P->A); P->C_total[r] = 95.90;
	P->nA = MkAlpha("NACGTU",DNA_MTRX);
	rewind(fptr);
	get_atoms(fptr, P);
	CountResPDB(P); // get the number of residues in each chain.

	P->chain[0]=0; P->chain[nChainsPDB(P)+1]=0;
	for(Int4 Nc=1; Nc <= nChainsPDB(P); Nc++){
	   assert(P->natoms[Nc] > 0);
	   P->chain[Nc]=AtomChain(GetAtomPDB(P,Nc,1));
	}
	return P;
}

void    ReNamePDB(const char *str,pdb_typ P)
{ P->filename=AllocString(str); }

void	LabelNullChainsPDB(pdb_typ P)
{ LabelNullChainsPDB(P,FALSE); }

void	LabelNullChainsPDB(pdb_typ P, BooLean RmX)
// Relabel Chains so that they all have distinct non-null labels...
// Return new chains in order...so can name them...
{
	char	chn,C;
	BooLean	*used,Relabel;
	Int4	Nc,n,chn_ptr;
	atm_typ	atom;

	NEW(used,300,BooLean);
	for(Nc=1; Nc <= nChainsPDB(P); Nc++){
      	   atom = GetAtomPDB(P, Nc, 1);
      	   C=AtomChain(atom);
	   // If blank or has occured before then relabel
	   if(C == ' ' || isdigit(C) || used[C] || (RmX && C == 'X')){
		for(C='A'; C <='Z'; C++){
		   if(!used[C] || (RmX && C == 'X')){ Relabel=TRUE; break; } 
		}
		if(!Relabel){
		  for(C='a'; C <='z'; C++){
		    if(!used[C] || (RmX && C == 'x')){ Relabel=TRUE; break; } 
		  }
		}
		if(!Relabel) print_error("FATAL: pdb file has too many chains");
	   } else Relabel=FALSE;	// Indicates chn should == C
	   used[C]=TRUE;
	   for(n=1; n <= P->natoms[Nc]; n++){
      	   	atom = GetAtomPDB(P,Nc,n); chn=AtomChain(atom);
		if(Relabel){
		    ReNameChainAtom(C,atom);
		} else {
		  if(chn != C){
		    PutAtom(stderr,atom);
		    print_error("FATAL: pdb file mislabeling");
		  }
		}
	   }
	} free(used);
}

void	get_atoms(FILE *fptr, pdb_typ P)
{
	char	str[200],state=' ';
	Int4	C,a,res;
	atm_typ	A;
	BooLean	na;

	for(a=1, C=1; fgets(str,195,fptr) != NULL; ){
	   //  if(strncmp(str,"REMARK",6) == 0) { // skip remarks
	   // } else if(strncmp(str,"HETATM",6) == 0) {
	   if(strncmp(str,"HETATM",6) == 0) {
		// if(state=='A'){ C++; a=1; }	// HETATM right after ATOM == New chain..
		A = MkAtom(str);
		if(A != NULL){ P->atom[C][a] = A; a++; }
		else P->natoms[C]--; // counted these earlier....
		state='H';
	   } else if(strncmp(str,"ATOM",4) == 0) {
		// if(state=='H'){ C++; a=1; }	// ATOM right after HETATM == New chain..
		A = MkAtom(str);
		if(A == NULL){
			fprintf(stderr,"str = \"%s\"",str);
			pdb_error("get_atoms( ) input error 1");;
		}
	        P->atom[C][a] = A;
		if(IsAtom(" P  ",A)){
		  // P->protein[C] = FALSE;
		  if(P->seq[C] != NULL){
		     P->seq[C][ResAtom(A)] = 
			GetResidueAtom(A, P->A, P->nA, &na);
		  }
		  if(na != FALSE) {
			// fprintf(stderr,"\nInput: %s \nna = %d\n",str,na);
			// pdb_error("get_atoms( ) input error 1");
		  }
		} else if(IsAtom(" CA ",A)){
		  if(P->seq[C] != NULL){
		     P->seq[C][ResAtom(A)] = GetResidueAtom(A,P->A,P->nA,&na);
		     res = ResAtom(A);
		  }
		  P->protein[C] = TRUE;
		  if(na != FALSE) {
			fprintf(stderr,"\n**** %s ****\n",str);
			pdb_error("get_atoms( ) input error 2");
		  }
		} a++;
		state='A';
	   } else if(strncmp(str,"TER",3) == 0){	/** TER **/
		C++; a = 1;
		state='T';
	   } else if(strncmp(str,"END",3) == 0){	/** END **/
		if(a > 1) C++; 
		state='E';
	   }				/** IGNORE REST **/
	}
}

BooLean	IsFullProteinChainPDB(Int4 c, pdb_typ P)
// Not just C_alpha carbons....
{
	if(c < 1 || c > nChainsPDB(P)) return FALSE;
	if(!IsProteinPDB(c,P)) return FALSE;
   	for(Int4 a=1; a<= P->natoms[c]; a++){
		atm_typ atom = P->atom[c][a];
		if(NitrogenAtom(atom) && IsAminoAcidAtom(atom)) return TRUE;
	} return FALSE;
}

void	count_atoms(FILE *fptr, pdb_typ P)
{
	char	str[200],tmp[5],state=' ';
	Int4	i,j,C,res,max = MAX_NGRPS_PDB;

	for(C=1; C <= max; C++){
		P->natoms[C] = 0; P->maxres[C] = 0; P->lowres[C] = INT4_MAX;
	}
	for(C=1; fgets(str,195,fptr) != NULL; ){
	   if(strncmp(str,"HETATM",6) == 0) {
		// if(state=='A') C++;	// ATOM right before HETATM == New chain..
		if(str[2] == 'T'){
		   P->natoms[C]++;	
		}
		state='H';
	   } else if(strncmp(str,"ATOM",4) == 0) {
		// if(state=='H') C++;	// HETERO right before ATOM == New chain..
		/** residue number **/
		for(i=22,j=0; i < 26; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
		if(sscanf(tmp,"%d",&res) != 1) {
			fprintf(stderr,"string: %s\n",str);
			pdb_error("input error");
		}
		P->maxres[C] = MAXIMUM(Int4,P->maxres[C],res);
		P->lowres[C] = MINIMUM(Int4,P->lowres[C],res);
		P->natoms[C]++;	
		state='A';
	   } else if(strncmp(str,"TER",3) == 0) {
/*********************************************************
Need to fix bug where TER followed by HETERO is treated like a 
second protein chain?
 *********************************************************/
		P->nchains = C; C++; P->maxres[C] = 0;
		state='T';
	   } else if(strncmp(str,"END",3) == 0){	/** END **/
		if(P->natoms[C] > 0){
		  P->nchains = C; C++; P->lowres[C]=P->maxres[C] = 0;
		}
		state='E';
	   }				/** IGNORE REST **/
	}
	if(P->nchains == 0 && P->natoms[1] > 0) P->nchains = 1;
	NEWP(P->atom,P->nchains +3, atm_typ);
	NEWP(P->seq,P->nchains +3, unsigned char);
	NEWP(P->seq0,P->nchains +3, unsigned char);
	NEW(P->protein,P->nchains +3, BooLean);
	for(C=1; C <= P->nchains; C++){
	   NEW(P->atom[C],P->natoms[C] +3, atm_typ);
	   if(P->maxres[C] > 0){
		if(P->lowres[C] < 0){
		  NEW(P->seq0[C],P->maxres[C]-P->lowres[C]+3, unsigned char);
		  P->seq[C] = P->seq0[C] - P->lowres[C]; 
		} else {
		  NEW(P->seq0[C],P->maxres[C]+3, unsigned char);
		  P->seq[C] = P->seq0[C];
		}
	   } else P->seq[C] = P->seq0[C] = NULL;
	}
}

#if 0
void    PutPDBSeq(FILE *fptr, Int4 Chain, pdb_typ P)
{ PutPDBSubSeq(fptr, Chain, start, end, P); }

void    PutPDBSubSeq(FILE *fptr, Int4 Chain, Int4 start, Int4 end, pdb_typ P)
{
}

Int4	SubSeqPDBSeq(e_type E, Int4 Chain, pdb_typ P)
// return 0 if E is NOT a subseq of Chain in P
// otherwise return the start site of E in P.
{
	Int4	r,i;

	if(Chain < 1 || Chain > P->nchains) pdb_error("PutPDBSeq( ): input error");
	if(P->protein){
	  fprintf(fptr,">%s|%c \n",P->filename,AtomChain(P->atom[Chain][1]));
	  for(i=P->lowres[Chain]; i<=P->maxres[Chain]; i++){
		r = P->seq[Chain][i];
		fprintf(fptr,"%c",AlphaChar(r,P->A));
		if(i%60==0) fprintf(fptr,"\n");
	  }
	  fprintf(fptr,"\n");
	}

}
#endif

e_type	GetPDBSeq(Int4 Chain, pdb_typ P)
{
	FILE *fp=tmpfile();
	assert(fp);
	Int4 len = PutPDBSeq(fp, Chain, P);
	rewind(fp);
	e_type  E = ReadSeq(fp,1,len,P->A);
	fclose(fp);
	return E;
}

Int4	PutPDBSeq(FILE *fptr, Int4 Chain, pdb_typ P)
// output the sequence in fasta format
{
	Int4	r,i,os,start;

	if(Chain < 1 || Chain > P->nchains) pdb_error("PutPDBSeq( ): input error");
	if(P->protein[Chain]){
	  os = P->lowres[Chain] - 1;
	  if(os > 0) fprintf(fptr,">%s|%c {|%d(0)|}\n",
			  P->filename,AtomChain(P->atom[Chain][1]),os);
	  else {
		os=0;
		fprintf(fptr,">%s|%c \n",P->filename,AtomChain(P->atom[Chain][1]));
	  }
#if 0
	  start=P->lowres[Chain];
#else
	  start=MAXIMUM(Int4,1,P->lowres[Chain]);
#endif
	  for(i=start; i<=P->maxres[Chain]; i++){
		r = P->seq[Chain][i];
		fprintf(fptr,"%c",AlphaChar(r,P->A));
		if(i%60==0) fprintf(fptr,"\n");
	  } fprintf(fptr,"\n");
#if 0
	  return (P->maxres[Chain] - P->lowres[Chain] + 1);
#else
	  // return (P->maxres[Chain]);
	  return (P->maxres[Chain] - os);
#endif
	} else return 0;
}

BooLean	ReNameChainPDB(Int4 chain, char Chain, pdb_typ P)
// rename chain number chain to 'Chain'.
{
	Int4    a;
	BooLean success=FALSE;
	atm_typ	atom;

	if(chain > 0 && chain <=  P->nchains){
	     success=TRUE;
   	     for(a=1; a<= P->natoms[chain]; a++){
		atom = P->atom[chain][a];
                ReNameChainAtom(Chain,atom);
             } P->chain[chain]=Chain;
	} return success;
}

BooLean	ShiftResChainPDB(char Chain, Int4 shift, pdb_typ P)
{ return ShiftResChainPDB(Chain, shift, INT4_MIN, INT4_MAX, P); }

BooLean	ShiftResChainPDB(char Chain, Int4 shift, Int4 start, Int4 end, pdb_typ P)
// Shift residue numbering in Chain by 'shift' from residue start to residue end.
{
	Int4    C,a,r;
	BooLean success=FALSE;
	atm_typ	atom;

	// fprintf(stderr,"start=%d; end=%d\n",start,end);
	for(C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) == Chain){
	     success=TRUE;
   	     for(a=1; a<= P->natoms[C]; a++){
		atom = P->atom[C][a];
		r=ResAtom(atom);
		if(r >= start && r <= end) RenumberAtom(r+shift,atom);
             }
	  }
	} return success;
}

BooLean	RenumberChainPDB(char Chain, Int4 start, pdb_typ P)
{
	atm_typ	atom;
#if 1
	Int4    C,a,r,diff;
	C=GetChainNumberPDB(P,Chain);
	if(C==0) return FALSE;
	diff=P->lowres[C]-start;
	if(diff==0) return TRUE;
   	for(a=1; a<= P->natoms[C]; a++){
	   atom = P->atom[C][a]; r=ResAtom(atom); RenumberAtom(r-diff,atom);
	} P->lowres[C] -= diff; P->maxres[C] -= diff;
	return TRUE;
#else
	Int4    C,a,r,res,inc;
	BooLean success=FALSE;
	Int4 r_last=INT4_MIN;
	for(res=start-1,C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) == Chain){
	     success=TRUE;
   	     for(a=1; a<= P->natoms[C]; a++){
		atom = P->atom[C][a];
		r=ResAtom(atom);
		if(r != r_last){ r_last = r; res++; }
                RenumberAtom(res,atom);
             } 
	  }
	} return success;
#endif
}

BooLean	PutNotPrtnChainPDB(FILE *fptr, char chain, pdb_typ P)
// Output all chains other than protein 'Chain' of P.
{
	Int4	C,a;
	BooLean	success=FALSE;

#if 1
	for(C=1; C <= P->nchains; C++){
	   char chn=ChainCharPDB(C,P);
	   if(chain != chn){ PutChainPDB(fptr,chn,P); success=TRUE; }
	} return success;
#else
	atm_typ	atm;
	Int4 Chain=GetChainNumberPDB(P, chain);
	for(C=1; C <= P->nchains; C++){
	     BooLean found = FALSE;
   	     for(a=1; a<= P->natoms[C]; a++){
		atm = P->atom[C][a];
	        if(AtomChain(atm) != Chain || 
		  (AtomChain(atm) == Chain && !IsAminoAcidAtom(atm))){
		  if(!success){
	            fprintf(fptr,"HEADER    %s\n",P->filename);
	            success=TRUE;
	          } if(!found){ found=TRUE; }
                  PutAtom(fptr, P->atom[C][a]);
		}
             } if(found) fprintf(fptr,"TER               \n");
	} return success;
#endif
}

BooLean	PutNotChainPDB(FILE *fptr, char chain, pdb_typ P)
// Output all chains other than protein 'Chain' of P.
{
	Int4	C,a;
	BooLean	success=FALSE;

	Int4 Chain=GetChainNumberPDB(P, chain);
	for(C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) != Chain){
	     if(!success) fprintf(fptr,"HEADER    %s\n",P->filename);
	     success=TRUE;
   	     for(a=1; a<= P->natoms[C]; a++){
                PutAtom(fptr, P->atom[C][a]);
             } fprintf(fptr,"TER               \n");
	  } else if(AtomChain(P->atom[C][1]) == Chain){
	     BooLean found=FALSE;
	     for(a=1; a<= P->natoms[C]; a++){
		if(IsHeteroAtom(P->atom[C][a])) {
		   PutAtom(fptr, P->atom[C][a]);
		   found=TRUE;
		}
             } if(found) fprintf(fptr,"TER               \n");
	  }
	} return success;
}

BooLean	PutPrtnChainPDB(FILE *fptr, char Chain, pdb_typ P)
// Output chain 'Chain' of P.
{
	Int4	C,a;
	BooLean	success=FALSE;

	for(C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) == Chain){
	     success=TRUE;
	     fprintf(fptr,"HEADER    %s\n",P->filename);
   	     for(a=1; a<= P->natoms[C]; a++){
                if(IsAminoAcidAtom(P->atom[C][a])) PutAtom(fptr, P->atom[C][a]);
             } fprintf(fptr,"TER               \n");
	     break;
	  }
	} return success;
}

BooLean	PutChainPDB(FILE *fptr, Int4 StartRes, Int4 EndRes, char Chain, pdb_typ P)
// Output chain 'Chain' of P.
{
	Int4	C,a,r;
	BooLean	success=FALSE;

	for(C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) == Chain){
	     success=TRUE;
	     fprintf(fptr,"HEADER    %s\n",P->filename);
   	     for(a=1; a<= P->natoms[C]; a++){
		r = ResAtom(P->atom[C][a]);
		if(r >= StartRes && r <= EndRes){
                   PutAtom(fptr, P->atom[C][a]);
	        }
             } fprintf(fptr,"TER               \n");
	     break;
	  }
	} return success;
}

BooLean	PutChainPDB(FILE *fptr, char Chain, pdb_typ P)
// Output chain 'Chain' of P.
{
	Int4	C,a;
	BooLean	success=FALSE;

	for(C=1; C <= P->nchains; C++){
	  if(AtomChain(P->atom[C][1]) == Chain){
	     success=TRUE;
	     fprintf(fptr,"HEADER    %s\n",P->filename);
   	     for(a=1; a<= P->natoms[C]; a++){
                PutAtom(fptr, P->atom[C][a]);
             } fprintf(fptr,"TER               \n");
	     break;
	  }
	} return success;
}

static BooLean	NearNonWaterPDB(atm_typ Atm, double mindist,pdb_typ P)
{
	Int4	C,a;
	atm_typ	atm;

	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		atm = P->atom[C][a];
		if(Atm != atm && !IsWaterAtom(atm)){
		   if(DistanceAtoms(Atm,atm) <= mindist) return TRUE;
		}
	   }
	} return FALSE;
}

BooLean *AdjacentHeteroSubunitsPDB(Int4 C0, double maxdist,pdb_typ P)
// Return boolean array indicating which chains are adjacent to chain C0
{
        Int4    C,a,a0;
        atm_typ atm,atm0;
        BooLean *OutPut;

        // find subunits...
        NEW(OutPut,MAX_NGRPS_PDB +3,BooLean);
        for(C=1; C <= P->nchains; C++){
           if(C == C0){ OutPut[C]=FALSE; continue; }
	   if(IsProteinPDB(C,P)) continue;
           for(a=1; a<= P->natoms[C]; a++){
             atm = P->atom[C][a];
	     if(IsWaterAtom(atm)) continue;
             if(IsHeteroAtom(atm)){
               for(a0=1; a0<= P->natoms[C0]; a0++){
                atm0 = P->atom[C0][a0];
                if(DistanceAtoms(atm, atm0) <= maxdist){ OutPut[C]=TRUE; break; }
               } if(OutPut[C]) break;
	     }
           }
        } return OutPut;
}

BooLean	*AdjacentSubunitsPDB(Int4 C0, double maxdist,pdb_typ P)
// Return boolean array indicating which chains are adjacent to chain C
{
	Int4	C,a,a0;
	atm_typ	atm,atm0;
	BooLean	*OutPut;

	// find subunits...
	NEW(OutPut,MAX_NGRPS_PDB +3,BooLean);
	for(C=1; C <= P->nchains; C++){
	   if(C == C0){ OutPut[C]=TRUE; continue; }
	   for(a=1; a<= P->natoms[C]; a++){
	     atm = P->atom[C][a];
	     for(a0=1; a0<= P->natoms[C0]; a0++){
	        atm0 = P->atom[C0][a0];
		if(DistanceAtoms(atm, atm0) <= maxdist){ OutPut[C]=TRUE; break; }
	     } if(OutPut[C]) break;
	   } 
	} return OutPut;
}

void    PutAdjacentSubunitsPDB(FILE *fptr, Int4 C0, double maxdist,pdb_typ P)
{ Int4	c0[4]; c0[0]=C0; c0[1]=0; MultiPutAdjacentSubunitsPDB(fptr, c0, maxdist,P); }

void    MultiPutAdjacentSubunitsPDB(FILE *fptr, Int4 *c0, double maxdist,pdb_typ P)
{
	Int4	C,a,a0,C0,i;
	atm_typ	atm,atm0;
	BooLean	*OutPut;

	// find subunits...
	NEW(OutPut,MAX_NGRPS_PDB +3,BooLean);
	for(C=1; C <= P->nchains; C++){
	  for(i=0; c0[i] > 0; i++){
	   C0=c0[i];
	   if(C == C0){ OutPut[C]=TRUE; continue; }
	   for(a=1; a<= P->natoms[C]; a++){
	     atm = P->atom[C][a];
	     for(a0=1; a0<= P->natoms[C0]; a0++){
	        atm0 = P->atom[C0][a0];
		if(DistanceAtoms(atm, atm0) <= maxdist){ OutPut[C]=TRUE; break; }
	     } if(OutPut[C]) break;
	   } if(OutPut[C]) break;
	  }
	}
	// Output subunits...
	fprintf(fptr,"HEADER    %s\n",P->filename);
	for(C=1; C <= P->nchains; C++){
	   if(!OutPut[C]) continue;
	   for(a=1; a<= P->natoms[C]; a++){
		atm = P->atom[C][a];
		// if(!IsAminoAcidAtom(atm))
		if(IsWaterAtom(atm)){ // then output only if near C0 chain
	         for(i=0; c0[i] > 0; i++){
		   for(C0=c0[i],a0=1; a0<= P->natoms[C0]; a0++){
	        	atm0 = P->atom[C0][a0];
			if(DistanceAtoms(atm, atm0) <= maxdist){ PutAtom(fptr, atm); break; }
	     	   }
		 }
		} else PutAtom(fptr, atm);
			 // if(IsHeteroAtom(A)) { }
	   } fprintf(fptr,"TER               \n");
	} free(OutPut);
}

void    RmDistWaterPDB(FILE *fptr, double maxdist,pdb_typ P)
{
	Int4	C,a;
	atm_typ	atm;

	fprintf(fptr,"HEADER    %s\n",P->filename);
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		atm = P->atom[C][a];
		if(!IsWaterAtom(atm)) PutAtom(fptr,atm);
		else if(NearNonWaterPDB(atm,maxdist,P)) PutAtom(fptr, P->atom[C][a]);
	   } fprintf(fptr,"TER               \n");
	}
}

void    TransformPDB(double xyz[4][5], pdb_typ P)
// Transform pdb file using coordinates as provided, for example, by CE output.
/*******************************************************************************
                  xyz[.][1]       xyz[.][2]        xyz[.][3]        xyz[.][4]
xyz[1]:   X2 = ( 0.983195)*X1 + ( 0.109599)*Y1 + ( 0.146000)*Z1 + (  -13.846273)
xyz[2]:   Y2 = (-0.059883)*X1 + ( 0.949108)*Y1 + (-0.309206)*Z1 + (   14.920112)
xyz[3]:   Z2 = (-0.172458)*X1 + ( 0.295267)*Y1 + ( 0.939721)*Z1 + (  -35.552807)
 *******************************************************************************/
{
	Int4	C,a;
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		TransformAtom(xyz,P->atom[C][a]);
	   }
	}
}

Int4	ChangePO3toAtmPDB(pdb_typ P)
// change SEP (phosphoserines) from HETATM to ATM.
// HETATM 4141  O1P SEP A 659     160.986   6.171  23.825  1.00 79.77
// change PTR (phosphotyrosine) from HETATM to ATM.
// HETATM 2393  O1P PTR B   6     -15.715 -13.488   1.650  1.00  0.00           O
// change THP (phosphothreonine) from HETATM to ATM.
// HETATM 1052  O2P THP   300      22.215  10.973  -0.301  1.00 59.28           O
{
	Int4	C,a,n=0;
	atm_typ	atm;

	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		atm=P->atom[C][a];
		if(IsHeteroAtom(atm)){
		   if(IsResAtom("SEP",atm)){ HetAtm2Atom(atm); n++; }
		   else if(IsResAtom("THP",atm)){ HetAtm2Atom(atm); n++; }
		   else if(IsResAtom("PTR",atm)){ HetAtm2Atom(atm); n++; }
		}
	   }
	}
	// fprintf(stderr,"%d SEP HETATM changed to ATOM.\n",n);
	return n;
	
}

Int4	ChangeMSE2MET_PDB(pdb_typ P)
{
	Int4	C,a,n=0;

	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		if(ChangeMSE2MetAtom(P->atom[C][a])) n++;
	   }
	}
	// fprintf(stderr,"%d MSE atoms changed to Met.\n",n);
	return n;
}

void    PutPDB(FILE *fptr, pdb_typ P)
{
	Int4	C,a;

	fprintf(fptr,"HEADER    %s\n",P->filename);
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		PutAtom(fptr, P->atom[C][a]);
	   }
	   fprintf(fptr,"TER               \n");
	}
}

void    NilPDB(pdb_typ P)
{
	Int4	C,a;

	free(P->filename);
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		if(P->atom[C][a] != NULL) NilAtom(P->atom[C][a]);
	   }
	   free(P->atom[C]);
	   if(P->seq0[C] != NULL) free(P->seq0[C]);
	}
	free(P->seq); free(P->seq0);
	free(P->atom);
	free(P->protein);
	free(P->temp);
	free(P->C_total);
	NilAlpha(P->A); NilAlpha(P->nA);
	free(P);
}

float	*OoiNumberPDB(Int4 C, pdb_typ P)
/** return the percent exposed surface for residues in chain C using 
Ooi numbers **/
{
	Int4	a,a2,r,res2,res;
	float	*CW,X,Y,Z,dmax,dmax2,dx,dy,dz,d;
	double	N,sum,ave,var,sd,x;
	atm_typ	A,A2;
	
	dmax = 8.0;  
	dmax2 = dmax*dmax;
	NEW(CW,P->maxres[C] +3, float);	/** NEED TO ALLOCATE maxres!! */
	for(r=1; r <= P->maxres[C]; r++) P->temp[r] = 0;
	for(sum=N=0, a = 1; a <= P->natoms[C]; a++){
	  A = P->atom[C][a];
	  if(IsAtom(" CA ",A)){
	    res = ResAtom(A);
	    if(P->temp[res] == 0){ N += 1; P->temp[res] = 1; }
	    /*** r = P->seq[C][res]; /****/
	    X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
	    for(a2 = 1; a2 <= P->natoms[C]; a2++){
	          A2 = P->atom[C][a2];
	          if(IsAtom(" CA ",A2)){
		    res2 = ResAtom(A2);
		    if(res2 < (res - 2) || res2 > (res + 2)){
			do{
			  dx = X - AtomX(A2);
			  if((dx *= dx) > dmax2) break;
			  dy = Y - AtomY(A2);
			  if((dy *= dy) > dmax2) break;
			  dz = Z - AtomZ(A2);
			  if((dz *= dz) > dmax2) break;
			  d = dx + dy + dz;
			  if((d=sqrt(d)) <= dmax){
				CW[res] += 1;
			  }
			} while(FALSE);
		    }
		  }
	    }
	    sum += CW[res];
	  }
	}
	/**** find stddev ***/
	ave = sum/N;
	for(var = 0,r = 1; r <=  P->maxres[C]; r++){
		if(P->temp[r] != 0){ x=CW[r]-ave; var += x*x; }
	}
	var = var/(N-1);
	sd = sqrt(var);
	for(r = 1; r <=  P->maxres[C]; r++){
		if(P->temp[r] != 0){ CW[r] = -(CW[r] - ave)/sd; }
	}
	/**** end find deviation ***/
	return CW;
}

void	NeighborsPDB(FILE *fptr, Int4 C, Int4 res, pdb_typ P)
/** return the percent exposed surface for residues in 
  chain C using Sander method (Biophysical J. 1990). **/
{
	Int4	a,a2,r;
	float	X,Y,Z,dmax,dmax2,dmin,dx,dy,dz,d;
	atm_typ	A,A2;
	
	dmax = 2.8 + 3.6;  /* Angst diameter Water + Atom */
	dmin = 3.6; /***/
	/** dmax = 2.1 + 2.7;  /* Angst diameter Water + Atom */
	/* dmin = 2.7; /***/
	dmax2 = dmax*dmax;
	for(a = 1; a <= P->natoms[C]; a++){
		A = P->atom[C][a];
		if(res == ResAtom(A)) break;
	}
	r = P->seq[C][res];  /** skip if r = glycine **/
	if(r != AlphaCode('G',P->A)){
	    do{
	     if(SideAtom(A)){
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(a2 = 1; a2 <= P->natoms[C]; a2++){
	          A2 = P->atom[C][a2];
		  if(ResAtom(A2) != res){
		   do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) < dmax) {
			fprintf(fptr,"select atomno=%d\n",AtomID(A2));
			fprintf(fptr,"color white\nspacefill 50\n");
		     }
		   } while(FALSE);
		  }
		}
	     }
	     a++;
	     if(a > P->natoms[C]) break;
	     A = P->atom[C][a];
	    } while(res == ResAtom(A));
	}
}

Int4	FindPotentialHBondsPDB(Int4 res0, Int4 C, Int4 C2, float dmax, pdb_typ P)
// Too redundant with other code; need to condense this at some point!
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,r,res,gly,res2;
	float	X,Y,Z,dmax2,dx,dy,dz,d;
	atm_typ	A,A2;
	BooLean	*contact;
	BooLean	na;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("ContactsResPDB( ): input error"); 
	}
	NEW(contact,P->maxres[C2]+3,BooLean);
	// dmax *= 3.6;  /* sum of radii of 2 Atoms in Angstroms */
	dmax2 = dmax*dmax;
	gly = AlphaCode('G',P->A);
	for(a = 1; a <= P->natoms[C]; a++){
	  A = P->atom[C][a];
	  if(!(OxygenAtom(A) || NitrogenAtom(A))) continue;
	  res = ResAtom(A);
	  if(res == res0){
	     r = P->seq[C][res0];  
	     // if(r==gly || SideAtom(A))	// requires sidechain!! fixed.
	     // if(C != C2 || SideAtom(A) || abs(res-res0) > 1)
	     if(TRUE)
	     { 
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	          A2 = P->atom[C2][a2];
		  if(A2 == A) continue;
		  if(!(OxygenAtom(A2) || NitrogenAtom(A2))) continue;
		  if(IsAminoAcidAtom(A) && IsAminoAcidAtom(A2) && 
			!SideAtom(A) && !SideAtom(A2)) continue;
		  res2 = ResAtom(A2);
#if 1	// Fix purify array bounds error....
		  char	IsRes;
		  IsRes= GetResidueAtom(A2, P->A, P->nA, &na);  // Is this an amino acid residue?
		// == 0 if not residue...??
#endif
		  // if(!(C==C2 && res2 == res0))
		  if(!(C==C2 && res2 == res0))
		  {
		   do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) <= dmax){

			char str[9],*str_ptr = AtomName(A);
			if(isspace(str_ptr[0])) str_ptr++;
			Int4 z;
			for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
			} str[z]=0;

			if(IsWaterAtom(A)){
			  fprintf(stdout,"#HOH%d\t// to ",res);
			} else {
			 if(AtomChain(A) != ' '){
			  fprintf(stdout,"%c%d%c",AlphaChar(r,P->A),res,AtomChain(A));
			 } else fprintf(stdout,"%c%d",AlphaChar(r,P->A),res);
	  		 // fprintf(stdout,"_%s.X\t// to ", str_ptr);
			 if(strcmp("o",str)==0) strcpy(str,"c=o");
	  		 fprintf(stdout,"_%s.X\t// to ",str);
			}
			if(AtomChain(A2) != ' '){
			  fprintf(stdout,"%s %s%d%c = %.2f A.\n",
					AtomName(A2),AtomResName(A2),res2,AtomChain(A2),d);
			} else fprintf(stdout,"%s %s%d = %.2f A.\n",
					AtomName(A2),AtomResName(A2),res2,d);

			if(IsRes) contact[res2] = TRUE;
		     }
		   } while(FALSE);
		  } 
		}
	     }
	  }
	}
	for(n=0, r=1; r<=P->maxres[C2]; r++) { if(contact[r]) { n++; } }
	free(contact);
	return n;
}

Int4	*ContactsPDB(Int4 C, Int4 N, Int4 *C2, float dmax, pdb_typ P)
/** return the residues in chain C within distance dmax of atoms 
in chain C2 **/
{
	Int4	a,a2,r,i,c,res,gly,*array;
	float	X,Y,Z,dmax2,dx,dy,dz,d;
	atm_typ	A,A2;
	BooLean	flag;
	
	if(C > P->nchains || C<1){ pdb_error("ContactsPDB( ): input error"); }
	dmax *= 3.6;  /* sum of radii of 2 Atoms in Angstroms */
	dmax2 = dmax*dmax;
	gly = AlphaCode('G',P->A);
	for(i=1,a = 1; a <= P->natoms[C]; a2++){
	  A = P->atom[C][a];
	  res = ResAtom(A);
	  r = P->seq[C][res];  /** skip if r = glycine **/
	  flag = FALSE;
	  do{
	     if(!flag && (r==gly || SideAtom(A))){ 
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(c = 1; c <= N; c++){
 		  if(C2[c] == C || C2[c] >P->nchains || C2[c] <1){
			pdb_error("ContactsPDB( ): input error");
		  }
		  for(a2 = 1; a2 <= P->natoms[C2[c]]; a2++){
	           A2 = P->atom[C2[c]][a2];
		   do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) <= dmax) {
	    		flag = TRUE;
		     }
		   } while(FALSE);
		  }
		}
	     }
	     a++;
	     if(a > P->natoms[C]) break;
	     A = P->atom[C][a];
	  } while(res == ResAtom(A));
	  if(flag) P->temp[i++] = res;
	}
	P->temp[i] = 0;
	NEW(array,i+3,Int4);
	for(r=1 ; r <=i; r++) array[r] = P->temp[r];
	return array;
}

float	*ExposedPDB(Int4 C, pdb_typ P)
/** return the percent exposed surface for residues in 
  chain C using Sander method (Biophysical J. 1990). **/
{
	Int4	a,a2,r,res;
	float	*CW,X,Y,Z,dmax,dmax2,dmin,dx,dy,dz,d;
	atm_typ	A,A2;
	
	dmax = 2.8 + 3.6;  /* Angst diameter Water + Atom */
	dmin = 3.6; /***/
	/** dmax = 2.1 + 2.7;  /* Angst diameter Water + Atom */
	/** dmin = 2.7; /***/
	dmax2 = dmax*dmax;
	NEW(CW,P->maxres[C] +3, float);	/** NEED TO ALLOCATE maxres!! */
	A = P->atom[C][1];
	for(a = 1; a <= P->natoms[C]; ){
	   res = ResAtom(A);
	   r = P->seq[C][res];  /** skip if r = glycine **/
	   if(r != AlphaCode('G',P->A)){
	    do{
	     if(SideAtom(A)){
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(a2 = 1; a2 <= P->natoms[C]; a2++){
	          A2 = P->atom[C][a2];
		  if(ResAtom(A2) != res){
		   do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) < dmax){
			if(d <= dmin) CW[res] += 1;
			else CW[res]+=1-(d-dmin)/(dmax-dmin);
		     }
		   } while(FALSE);
		  }
		}
	     }
	     a++;
	     if(a > P->natoms[C]) break;
	     A = P->atom[C][a];
	    } while(res == ResAtom(A));
/*****/
	    CW[res] = 100*(P->C_total[r] - CW[res])/P->C_total[r];
/*****
	    CW[res] = (P->C_total[r] - CW[res])/3.0;
/*****/
	   } else { 
	    a++;
	    if(a <= P->natoms[C]) A = P->atom[C][a]; 
	    CW[res] = 0; 
	   }
	}
	return CW;
}

atm_typ	ContactsAtomPDB(Int4 atom_no, float dmax, pdb_typ P)
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	a,r,res,res2,C;
	float	X,Y,Z,dmax2,dx,dy,dz,d;
	atm_typ	atom,A;
	BooLean	na;
	
	for(C=1; C <= P->nchains; C++){
	   for(atom=0,a=1; a <= MaxAtomsPDB(C,P); a++){
		atom = GetAtomPDB(P,C,a);
		if(AtomID(atom) == atom_no){
			PutAtom(stdout,atom);
			break;
		}
		else atom=0;
	   } if(atom) break;
	}
	// dmax *= 3.6;  /* sum of radii of 2 Atoms in Angstroms */
	if(atom){	// Found atom
	    dmax2 = dmax*dmax;
	    res2 = ResAtom(atom);
	    X = AtomX(atom); Y = AtomY(atom); Z = AtomZ(atom);
	    for(C=1; C <= P->nchains; C++){
		for(a = 1; a <= P->natoms[C]; a++){
	  	  A = P->atom[C][a];
		  res = ResAtom(A);
		  char	aa=GetResidueAtom(A, P->A, P->nA, &na); 
		  // Is this an amino acid residue?
		
		  if(A!= atom){
		     dx = X - AtomX(A);
		     if((dx *= dx) >= dmax2) continue;
		     dy = Y - AtomY(A);
		     if((dy *= dy) >= dmax2) continue;
		     dz = Z - AtomZ(A);
		     if((dz *= dz) >= dmax2) continue;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) <= dmax){
			if(AtomChain(A) != ' '){
			  fprintf(stdout,"#%c%d%c_%s.X\t//",AlphaChar(aa,P->A),res,AtomChain(A),
				AtomLowerName(A));
			} else fprintf(stdout,"#%c%d_%s.X\t//",
				AlphaChar(aa,P->A),res,AtomLowerName(A));
	  		// fprintf(stdout,"%s to ", AtomName(A));
	  		fprintf(stdout," to ");
			if(AtomChain(atom) != ' '){
			  fprintf(stdout,"%s %s%d%c = %.2f A.\n",
				AtomName(atom),AtomResName(atom),res2,AtomChain(atom),d);
			} else fprintf(stdout,"%s %s%d = %.2f A.\n",
				AtomName(atom),AtomResName(atom),res2,d);
		     }
		  }
		}
	    }
	}
	return atom;
}


Int4	ContactsResPDB(Int4 res0, Int4 C, Int4 C2, float dmax, pdb_typ P)
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,r,res,gly,res2;
	float	X,Y,Z,dmax2,dx,dy,dz,d;
	atm_typ	A,A2;
	BooLean	*contact;
	BooLean	na;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		if(C > P->nchains)
		  fprintf(stderr,"C==%d > nchains==%d\n",C,P->nchains);
		if(C2 > P->nchains)
		  fprintf(stderr,"C2==%d > nchains==%d\n",C2,P->nchains);
		if(C<1 || C2<1)
		  fprintf(stderr,"C==%d; C2==%d.\n", C,C2);
		pdb_error("ContactsResPDB( ): input error"); 
	}
	NEW(contact,P->maxres[C2]+3,BooLean);
	// dmax *= 3.6;  /* sum of radii of 2 Atoms in Angstroms */
	dmax2 = dmax*dmax;
	gly = AlphaCode('G',P->A);
	char	IsRes;
	for(a = 1; a <= P->natoms[C]; a++){
	  A = P->atom[C][a];
	  res = ResAtom(A);
	  if(res == res0){
	     if(P->seq[C]) r = P->seq[C][res0];  
	     else r = 0;
	     // if(r==gly || SideAtom(A))	// requires sidechain!! fixed.
	     if(C != C2 || SideAtom(A) || abs(res-res0) > 1)
	     { 
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	          A2 = P->atom[C2][a2];
		  res2 = ResAtom(A2);
#if 1	// Fix purify array bounds error....
		  // if(!IsAminoAcidAtom(A2)) continue;
		  IsRes= GetResidueAtom(A2, P->A, P->nA, &na);  
		  // Is this an amino acid residue?
		// == 0 if not residue...??
#endif
#if 0
		  Int4 r2 = P->seq[C2][res2];  
		  // if(!(C==C2 && res2 == res0) && (r2==gly || SideAtom(A2))){
#endif
		  if(!(C==C2 && res2 == res0)){
		   do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) <= dmax){
			if(AtomChain(A) != ' '){
			  fprintf(stdout,"#%c%d%c.Y\t//",AlphaChar(r,P->A),res,AtomChain(A));
			} else fprintf(stdout,"#%c%d.Y\t//",AlphaChar(r,P->A),res);
	  		fprintf(stdout,"%s to ", AtomName(A));
			if(AtomChain(A2) != ' '){
			  fprintf(stdout,"%s %s%d%c = %.2f A.\n",
					AtomName(A2),AtomResName(A2),res2,AtomChain(A2),d);
			} else fprintf(stdout,"%s %s%d = %.2f A.\n",
					AtomName(A2),AtomResName(A2),res2,d);

			if(IsRes) contact[res2] = TRUE;
		     }
		   }while(FALSE);
		  } 
		}
	     }
	  }
	}
	for(n=0, r=1; r<=P->maxres[C2]; r++) { if(contact[r]) { n++; } }
	free(contact);
	return n;
}

double	DistancePDB(Int4 res1, char *name1, Int4 C1, Int4 res2, 
	char *name2, Int4 C2, pdb_typ P)
// return the distance in Angstroms between res1 in chain C1 and
// res2 in chain C2.
{
	Int4	a,res,leng_name;
	double	dx,dy,dz,d,X1,Y1,Z1,X2,Y2,Z2;
	atm_typ	A,A1,A2;
	char	*name;
	
	if(C1 > P->nchains || C1 < 1 || C2 >P->nchains || C2<1){
		pdb_error("DistancePDB( ): input error"); 
	}
	A1=A2=0;
	for(a=1; a <= P->natoms[C1]; a++){
	  A = P->atom[C1][a];
	  res = ResAtom(A);
	  if(res == res1){
	  	name = AtomName(A);
		leng_name=strlen(name1);
		if(strncmp(name,name1,leng_name)){
		   A1 = A; X1=AtomX(A1); Y1=AtomY(A1); Z1=AtomZ(A1);
		   break;
		}
	  }
	}
	for(a=1; a <= P->natoms[C2]; a++){
	   A = P->atom[C2][a];
	   res = ResAtom(A);
	   if(res == res2){
	  	name = AtomName(A);
		leng_name=strlen(name2);
		if(strncmp(name,name2,leng_name)){
		   A2 = A; X2=AtomX(A2); Y2=AtomY(A2); Z2=AtomZ(A2);
		   break;
		}
	   }
	}
	if(A1 == 0 || A2 == 0) return -1.0; // -1.0 returned if atom nonexistant
	dx = X1 - X2; dx *= dx;
	dy = Y1 - Y2; dy *= dy;
	dz = Z1 - Z2; dz *= dz;
	d = dx + dy + dz; d=sqrt(d); 
#if 1
	fprintf(stdout,"%s%d to %s%d: %.2f A\n",
				AtomResName(A1),res1,AtomResName(A2),res2,d);
#endif
	return d;
}

atm_typ	AtomPDB(Int4 C, Int4 atom_num, pdb_typ P)
{
	if(C > P->nchains || C < 1 ) return NULL;
	if(atom_num < 1 || atom_num > P->natoms[C]) return NULL;
	else return P->atom[C][atom_num];
}

Int4	TriadPDB(Int4 C, pdb_typ P, Int4 *triad[5])
// return the number of residues in chain C within distance dmax of atoms in 
// residue res0 of chain C 
// triad[1] = Oxy; triad[2] = ND1 (his); triad[3] =NE2; triad[4] = Oxy
{
	Int4	n,a,a2,r,r2,res,his,res2;
	float	X,Y,Z,dmax2,dmax,dx,dy,dz,d;
	atm_typ	A,A2;
	char	*atom_name,*atom2_name;
	Int4	num_his,his_id,his_ptr;
	Int4	i,*map;
	BooLean	nd1;
	float	*dist1,*dist2;
	
	if(C > P->nchains || C < 1 ){
		pdb_error("TriadPDB( ): input error"); 
	}
	his = AlphaCode('H',P->A);
	for(num_his=0, r=1; r<=P->maxres[C]; r++) if(P->seq[C][r]==his) num_his++;
	if(num_his == 0) return 0;
	NEW(map,P->maxres[C]+2,Int4);
// fprintf(stdout,"%d his\n",num_his);
	for(i=0; i < 5; i++) NEW(triad[i],num_his+3,Int4);
	NEW(dist1,num_his+3,float); NEW(dist2,num_his+3,float);
	for(i=1; i < num_his; i++) dist1[i]=dist2[i]=-1.0;
	dmax = 1.2*3.6;  /* sum of radii of 2 Atoms in Angstroms */
	dmax2 = dmax*dmax;
	for(his_ptr=0,a=1; a <= P->natoms[C]; a++){
	  A = P->atom[C][a];
	  res = ResAtom(A);
	  r = P->seq[C][res];  
	  if(r==his){ 
	     atom_name=AtomName(A);
	     if(SideAtom(A) && (strncmp(atom_name," ND1",4) == 0
	         			|| strncmp(atom_name," NE2",4) ==0)){ 
	// fprintf(stdout,"atom=%s%d; SideAtom=%d\n",AtomName(A),res,SideAtom(A));
		if(strncmp(atom_name," ND1",4) == 0) nd1 = TRUE;
		else nd1 = FALSE;
		X = AtomX(A); Y = AtomY(A); Z = AtomZ(A);
		for(a2 = 1; a2 <= P->natoms[C]; a2++){
	          A2 = P->atom[C][a2];
		  atom2_name=AtomName(A2);
		  if(SideAtom(A2) && strncmp(atom2_name," O",2)==0){
		   res2 = ResAtom(A2);
		   // fprintf(stdout,"atom=%s%d; SideAtom=%d\n",AtomName(A2),res2,SideAtom(A2));
	     	   r2 = P->seq[C][res2];  
		   if(!(res2 == res)){
		    do{
		     dx = X - AtomX(A2);
		     if((dx *= dx) >= dmax2) break;
		     dy = Y - AtomY(A2);
		     if((dy *= dy) >= dmax2) break;
		     dz = Z - AtomZ(A2);
		     if((dz *= dz) >= dmax2) break;
		     d = dx + dy + dz;
		     if((d=sqrt(d)) <= dmax){
			// found potential catalytic histidine...
			if(map[res] == 0){ his_ptr++; map[res]=his_ptr; } 
			his_id=map[res];
			if(nd1) {
				triad[1][his_id] = a; 
				if(triad[2][his_id] <= 0. || dist1[his_id] > d){
					triad[2][his_id] = a2; 
					dist1[his_id] = d;
				}
			} else {
				triad[3][his_id] = a; 
				if(triad[4][his_id] <= 0. || dist2[his_id] > d){
					triad[4][his_id] =a2; 
					dist2[his_id] = d;
				}
			}
#if 0
			fprintf(stdout,"res=%s%d; atom=%s; %.2f Angstroms from: \n",
				AtomResName(A2),res2,AtomName(A2),d);
	  		fprintf(stdout,"\tres=%s%d; atom=%s.\n",
				AtomResName(A),res,AtomName(A));
#endif
		     }
		    } while(FALSE);
		   } 
		  }
	        }
	     }
	  }
	}
	for(n=0, his_id=1; his_id<=num_his; his_id++) {
		if(triad[1][his_id] != 0 && triad[3][his_id] != 0){
			 triad[0][his_id] = 1; n++;
		} else triad[0][his_id] = 0;
	}
	free(map);
	free(dist1); free(dist2);
	return num_his;
}

void    pdb_error(const char *s){fprintf(stderr,"PDB: %s\n",s);exit(1);}

long	CountHydrogensPDB(pdb_typ P) 
{
  Int4	   Nc,n;
  long	   Hydrogens=0;
  atm_typ  atom;

  for(Nc = 1; Nc <= nChainsPDB(P); Nc++){
    for (n=1; n <= P->natoms[Nc]; n++){
      atom = GetAtomPDB(P, Nc, n);
      if(HydrogenAtom(atom)) Hydrogens++;
    }
  } return Hydrogens;
}

//*********************** From HATA code ******************************
/*------------------------------------------------------------------------*/

void	CountResPDB(pdb_typ P) 
{
  int Nchain, Natom;
  int NrPrevious, Nr, NrCount;

  for (Nchain = 1; Nchain <= P->nchains; Nchain++){
    NrPrevious = INT4_MIN; NrCount = 0;
    for (Natom = 1; Natom <= P->natoms[Nchain]; Natom++){
        Nr = ResAtom(P->atom[Nchain][Natom]);
        if(Nr != NrPrevious){ NrCount++; NrPrevious = Nr; }
    } P->num_res[Nchain] = NrCount;
  }
}

/*------------------------------------------------------------------------*/

Int4 GetChainNumberPDB(pdb_typ P, char Chain)
{
  Int4	Nc,n;
  atm_typ  atom;
  char	chn;

  for(Nc = 1; Nc <= nChainsPDB(P); Nc++){
    for (n=1; n <= P->natoms[Nc]; n++){
      atom = GetAtomPDB(P, Nc, n);
      chn=AtomChain(atom);
      // fprintf(P->fp_stderr,"%d(%d): chn=%c; chain=%c\n",Nc,n,chn,Chain);
      if(isalpha(chn)){
       if(AtomChain(atom) == Chain) return Nc;
       break;
      } 
    }
  }
  return 0;
}

/*------------------------------------------------------------------------*/

Int4 GetResMinAtomNumber(pdb_typ P, Int4 Chain, Int4 Nres)
{
  Int4      i=1;
  atm_typ  atom;

  while (i <= MaxAtomsPDB(Chain, P)){

    atom = GetAtomPDB(P, Chain, i);
    if ( ResAtom(atom) == Nres ) return i; 
    i++;
  }
  return 0;
}

/*------------------------------------------------------------------------*/

Int4 GetResNumberOfAtoms(pdb_typ P, Int4 Chain, Int4 Nres)
{
  Int4 N;

  if (Nres == MaxResPDB(Chain, P))
    N = MaxAtomsPDB(Chain, P) - GetResMinAtomNumber(P, Chain, Nres);
  else 
    N = GetResMinAtomNumber(P, Chain, Nres+1) -
        GetResMinAtomNumber(P, Chain, Nres);

  return N;
}

/*------------------------------------------------------------------------*/

//***************** Hata code: from res_typ ************************
/*------------------------------------------------------------------------*/
static Int4 MaxResidueNumber = 5000;

res_typ *MakeResPDB(Int4 chain,Int4 *num_res,pdb_typ P)
{
  Int4		n,a,NrPrevious,Nr,NrCount,ia;
  atm_typ	atom,atm[5000];
  res_typ	*Res=0;

  NrCount    = 0; NrPrevious = INT4_MIN;
  for(n=1; n<= nChainsPDB(P); n++){
   if(n == chain){
    NEW(Res,MaxResidueNumber+3,res_typ);
    // Res=new res_typ[MaxResidueNumber];
    for(ia=0,a=1; a <= MaxAtomsPDB(n,P); a++){
      atom = GetAtomPDB(P,n,a);
      Nr = ResAtom(atom);
      if(Nr != NrPrevious){
#if 0
	assert(NrCount < MaxResidueNumber);
#else
	if(NrCount >= MaxResidueNumber){
		fprintf(stderr,"%s: NrCount = %d\n",FilenamePDB(P),NrCount);
		assert(NrCount < MaxResidueNumber);
	}
#endif
        if(ia > 0){ Res[NrCount]=MkRes(ia,n,atm); }
	else Res[NrCount]=0;
	NrCount++; NrPrevious=Nr; ia=0; 
      }
      ia++; 
      assert(ia < 5000);
      atm[ia]=atom;
    }
    assert(NrCount < MaxResidueNumber);
    if(ia > 0) Res[NrCount] = MkRes(ia,n,atm);
    *num_res=NrCount;
    break;
   }
  }
  // PutBondLengthHGs(stderr);
  return Res;
}

/*------------------------------------------------------------------------*/


#define MAX_angleDHA 90

Int4	FindWaterHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
		char color, pdb_typ P,Int4 water)
// Too redundant with other code; need to condense this at some point!
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,r,r2,res,gly,res2;
	atm_typ	H,Donor,Accept;
	BooLean	na;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("FindBondsPDB( ): input error"); 
	}
	gly = AlphaCode('G',P->A);
	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;
	     Donor = AtomResidue(a,Res);
	     res=ResAtom(Donor);
// if(res==221) fprintf(stderr,"***************** res 221 **************\n");
	     if(P->maxres[C] < res) r = 0;
	     else { if(P->seq[C]) r = P->seq[C][res];  else r = 0; }
	     // Find all Donor atoms:
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
		for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	          Accept = P->atom[C2][a2];
		  if(HydrogenAtom(Accept)) continue;
		  if(Accept == Donor) continue;
		  res2 = ResAtom(Accept);
		  if(Res2 && Res2 != res2) continue;
	     	  if(P->seq[C2] && IsAminoAcidAtom(Accept)) r2 = P->seq[C2][res2];  
		  else r2=0;
		  // if(C==C2 && abs(res2 - res) < 2) continue;
// if(r2==221) fprintf(stderr,"***************** res 221 **************\n");
		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  // double angle_HAX = CalcAngle(, atm_typ atm2, atm_typ atm3);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  // if(d_HA <= 2.5 && d_DA > d_HA && angle_DHA >= 90){
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MAX_angleDHA){
		    atm_typ	A,B;
		    Int4 rA,rB;
		    Int4 isaccept=0;
		    unsigned char rcA;
		    A=Donor;B=Accept; rcA=r; 
		    BooLean NeitherCarbons=TRUE;
		    if(CarbonAtom(A) || CarbonAtom(B)) NeitherCarbons=FALSE;
		    for(rcA=r; isaccept < 2; A=Accept,B=Donor,rcA=r2,isaccept++){
			rA=ResAtom(A); rB=ResAtom(B);
			if(water == 0) {
			  if(!(IsWaterAtom(A) || IsWaterAtom(B))) continue;
			} else {
			  if(!(IsWaterAtom(A) && rA ==water) &&
				!((IsWaterAtom(B) && rB==water))) continue;
			}
			char str[20],*str_ptr = AtomName(A);
			if(IsWaterAtom(A) && A==Donor) str_ptr = AtomName(H);
			if(isspace(str_ptr[0])) str_ptr++;
			Int4 z,zz;
			for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
			} str[z]=0;
			if(IsWaterAtom(A)){
		  	// if(ResAtom(Accept)==567) PutAtom(stderr,Accept); // debug...
			  if(A==Donor){
			  	fprintf(fp,"#HOH%d_o-%s.X\t// to ",rA,str);
			  } else fprintf(fp,"#HOH%d_o.X\t// to ",rA);
			} else {
			  if(!isaccept){	// i.e. donor atom...
			    str_ptr = AtomName(H);
			    if(isspace(str_ptr[0])) str_ptr++;
			    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
			    } str[z]=0;
			  }
			 if(strcmp("o",str)==0) strcpy(str,"c-o");
			 else if(str[0]=='c' && strstr(str,"-") == 0){
				 fprintf(fp,"// ");
			 }
			 if(NeitherCarbons) fprintf(fp," ");
			 else fprintf(fp," ");
			 // else fprintf(fp,"#");
			 if(IsHeteroAtom(A)) fprintf(fp,"!");
			 if(AtomChain(A) != ' '){
			  if(rcA) fprintf(fp,"%c%d%c",AlphaChar(rcA,P->A),rA,AtomChain(A));
			  else fprintf(fp,"%s%d%c",AtomResName(A),rA,AtomChain(A));
			 } else {
			  if(rcA) fprintf(fp,"%c%d", AlphaChar(rcA,P->A),rA);
			  else fprintf(fp,"%s%dX", AtomResName(A),rA);
			 } fprintf(fp,"_%s.X\t// to ",str);
			}
			if(AtomChain(B) != ' '){
			  fprintf(fp,"%s %s%d%c (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),
					rB,AtomChain(B),d_DA,d_HA,angle_DHA);
			} else fprintf(fp,"%s %s%d (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),rB,d_DA,d_HA,angle_DHA);
		     } n++;
		  } 
		}
	     }
	} return n;
}


Int4	FindHBondsPDB(res_typ Res, Int4 C, Int4 C2, float dmax, pdb_typ P)
{ return FindHBondsPDB(stdout,Res,C,C2,dmax,0,0,P); }

Int4	FindHBondsPDB(FILE *fp,res_typ Res, Int4 C, Int4 C2, float dmax, pdb_typ P)
{ return FindHBondsPDB(fp,Res,C,C2,dmax,0,0,P); }

Int4	FindHBondsPDB(res_typ Res, Int4 C, Int4 C2, float dmax, Int4 Res2, pdb_typ P)
{ return FindHBondsPDB(stdout,Res, C, C2, dmax, Res2, 0,P); }

Int4	FindHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
		char color, pdb_typ P)
{ return FindHBondsPDB(fp,Res,C,C2,dmax,Res2,color,0,P); }

Int4	FindHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
		char color, unsigned short file_id, pdb_typ P)
// Too redundant with other code; need to condense this at some point!
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,r,r2,res,res2;
	atm_typ	H,Donor,Accept;
	BooLean	na;
	FILE	*efp=P->fp_stderr;
	BooLean	comments=FALSE;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("FindBondsPDB( ): input error"); 
	}
	if(efp==stderr) comments=TRUE; else comments=FALSE;
	comments=TRUE; // waiting to debug this option;
	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;
	     Donor = AtomResidue(a,Res);
	     res=ResAtom(Donor);
#if 0	// DEBUG...AFN: 7-31-12; was problem with #define HydrogenAtom(A).
	     if(strstr(AtomName(Donor),"NZ") != 0 && res == 13){
			PutAtom(stderr,Donor); PutAtom(stderr,Accept);
	     }
#endif
	     if(P->maxres[C] < res) r = 0;
	     else if(P->seq[C]) r = P->seq[C][res];  else r = 0; 
	     // Find all Donor atoms:
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
		for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	          Accept = P->atom[C2][a2];
		  if(HydrogenAtom(Accept)) continue;
		  if(Accept == Donor) continue;
		  res2 = ResAtom(Accept);
#if 0	// DEBUG...AFN: 7-31-12.
		  if(strstr(AtomName(Donor),"NZ") != 0 && OxygenAtom(Accept) && !SideAtom(Accept) && res2 > 78 && res2 < 81){
			PutAtom(stderr,Donor); PutAtom(stderr,Accept);
		  }
#endif
#if 1	// Do I want to see backbone-to-backbone interactions?
		  // color --> calling from automatic script...
		  if(color && IsAminoAcidAtom(Donor) && IsAminoAcidAtom(Accept) &&
			!SideAtom(Donor) && !SideAtom(Accept)) continue;
#endif
		  if(Res2 && Res2 != res2) continue; 	// 10-7-09: afn bug fix...???
	     	  if(P->seq[C2] && IsAminoAcidAtom(Accept)) r2 = P->seq[C2][res2];  
		  else r2=0;
		  if(C==C2 && abs(res2 - res) < 2) continue;
#if 1	// Fix purify array bounds error....
		  char	IsRes;
		  IsRes= GetResidueAtom(Accept, P->A, P->nA, &na);
		  // Is this an amino acid residue?
#endif
		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  // double angle_HAX = CalcAngle(, atm_typ atm2, atm_typ atm3);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  // if(d_HA <= 2.5 && d_DA > d_HA && angle_DHA >= 90){
#if 0	// DEBUG...AFN: 7-31-12.
		  if(strstr(AtomName(Donor),"NZ") != 0 && OxygenAtom(Accept) && !SideAtom(Accept) && res2 > 78 && res2 < 81){
			PutAtom(stderr,Donor); PutAtom(stderr,Accept);
			fprintf(stderr,"NZ: d_HA =%.3f; angle_DHA = %.3f; d_DA = %.3f\n", d_HA,angle_DHA,d_DA);
		  }
#endif
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MAX_angleDHA){
		    atm_typ	A,B;
		    Int4 rA,rB;
		    Int4 isaccept=0;
		    unsigned char rcA;
		    A=Donor;B=Accept; rcA=r; 
		    BooLean NeitherCarbons=TRUE;
		    if(CarbonAtom(A) || CarbonAtom(B)) NeitherCarbons=FALSE;
		    for(rcA=r; isaccept < 2; A=Accept,B=Donor,rcA=r2,isaccept++){
			rA=ResAtom(A); rB=ResAtom(B);
			char str[20],*str_ptr = AtomName(A);
			if(IsWaterAtom(A) && A==Donor) str_ptr = AtomName(H);
			if(isspace(str_ptr[0])) str_ptr++;
			Int4 z,zz;
			for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
			} str[z]=0;
			if(IsWaterAtom(A)){
		  	// if(ResAtom(Accept)==567) PutAtom(stderr,Accept); // debug...
			  if(A==Donor){
			     if(comments) fprintf(fp,"#HOH%d_o-%s.X\t// to ",rA,str);
			     else fprintf(fp,"#HOH%d_o-%s.X\n",rA,str);
			  } else {
			     if(comments) fprintf(fp,"#HOH%d_o.X\t// to ",rA);
			     else fprintf(fp,"#HOH%d_o.X\n",rA);
			  }
			} else {
#if 1	// add hydrogen too...
			  if(!isaccept){	// i.e. donor atom...
			    str_ptr = AtomName(H);
			    if(isspace(str_ptr[0])) str_ptr++;
			    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
			    } str[z]=0;
			  }
#endif
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
			 if(AtomChain(A) != ' '){
			  if(rcA) fprintf(fp,"%c%d%c",AlphaChar(rcA,P->A),rA,AtomChain(A));
			  else fprintf(fp,"%s%d%c",AtomResName(A),rA,AtomChain(A));
			 } else {
			  if(rcA) fprintf(fp,"%c%d", AlphaChar(rcA,P->A),rA);
			  else fprintf(fp,"%s%dX", AtomResName(A),rA);
			 }
			 if(comments) fprintf(fp,"_%s.X\t// to ",str);
			 else fprintf(fp,"_%s.X\n",str);
			}
			if(AtomChain(B) != ' '){
			 if(comments) fprintf(fp,
				"%s %s%d%c (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),
					rB,AtomChain(B),d_DA,d_HA,angle_DHA);
			} else if(comments) fprintf(fp,
				"%s %s%d (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),rB,d_DA,d_HA,angle_DHA);
		     } n++;
		  } 
		}
	     }
	} return n;
}

static double	ProjectOnPlane(atm_typ H, atm_typ X, atm_typ C0, atm_typ C1)
// For atom H above a plane P defined by atoms X, C0 and C1, return the
// distance of the Length of the projection of H onto that plane from X.
{
	// vec_typ	U,A,P,V1,V2;

	vec_typ V1 = MkVector(AtomX(C0)-AtomX(X),AtomY(C0)-AtomY(X),AtomZ(C0)-AtomZ(X));
	vec_typ V2 = MkVector(AtomX(C1)-AtomX(X),AtomY(C1)-AtomY(X),AtomZ(C1)-AtomZ(X));
	vec_typ Vd = GetVectorProduct(V1,V2);
	vec_typ VH = MkVector(AtomX(H),AtomY(H),AtomZ(H));
	vec_typ VX = MkVector(AtomX(X),AtomY(X),AtomZ(X));
	vec_typ Vp = GetVectorProjectionOntoPlane(Vd,VX,VH);
#if 0
	atm_typ Virtual =VirtualAtom(Vp->x + AtomX(X),Vp->y + AtomY(X),Vp->z + AtomZ(X));
	// PutAtom(stderr,Virtual);
	// PutAtom(stderr,X);
#endif
	double d= sqrt(GetVectorLength(Vp));
	delete Vp; delete V1; delete V2; delete Vd; delete VX; delete VH;
	return d;
}

Int4	FindAromaticHbondsPDB(res_typ ResD, res_typ ResA,Int4 C, Int4 C2, 
	float dmax, pdb_typ P)
{ return FindAromaticHbondsPDB(stdout,ResD, ResA,C, C2, dmax, 0,P); }

Int4	FindAromaticHbondsPDB(FILE *fp, res_typ ResD, res_typ ResA,Int4 C, Int4 C2, 
	float dmax, char color, pdb_typ P)
{ double d=-1.0; return FindAromaticHbondsPDB(fp,ResD, ResA,C, C2, dmax,0,P,d); }

Int4	FindAromaticHbondsPDB(FILE *fp, res_typ ResD, res_typ ResA,Int4 C, Int4 C2, 
	float dmax, char color, pdb_typ P,double &minDist)
// find potential aromatic pi-H interactions.
{
	Int4	n,a,d,r,r2,res,res2,numhits;
	atm_typ	H,Donor,Accept,atom[20];
	double	d_HA,d_DA,angle_DHA,Ave_d_DA,Ave_d_HA;
	BooLean	na;
	a_type	AB=P->A;
	FILE	*efp=P->fp_stderr;
	BooLean	comments=FALSE;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("FindAromaticHbondsPDB( ): input error"); 
	}
	if(efp == stderr) comments=TRUE; else comments=FALSE;
	comments=TRUE; // waiting to debug this option;
	if(dmax == 0.0) dmax = 4.5;
	if(ResD == ResA) return 0;
	char    aa = GetResidueAtom(AtomResidue(1,ResA),AB,P->nA,&na);
	aa = AlphaChar(aa,AB);
	char    aa2 = GetResidueAtom(AtomResidue(1,ResD),AB,P->nA,&na);
	aa2 = AlphaChar(aa2,AB);
	if(na || !(aa == 'H' || aa == 'F' || aa == 'Y' || aa == 'W')) return 0;
	atm_typ	Near,Far;
	double	d_near,d_far,angle_HNF;
	mh_type	mH=Mheap(6,3);
	atm_typ	atom6[10];
	Int4	x;

	for(n=0,d = 1; d <= ResidueAtomNumber(ResD); d++){
	     if(!ResidueAtomHydrogens(d,ResD)) continue;
	     Donor = AtomResidue(d,ResD);
// PutAtom(stderr,Donor);
	     res=ResAtom(Donor);
	     if(P->maxres[C] < res) r = 0;
	     else {
	       if(P->seq[C]) r = P->seq[C][res];  
	       else r = 0;
	     }
	     for(Int4 h=1; (H=ResidueAtomHydrogen(d,h,ResD)); h++){
		Ave_d_DA=Ave_d_HA=0.0;
// PutAtom(stderr,H);
		d_near=1000.0; d_far=0.0;
	        for(numhits=0,a = 1; a <= ResidueAtomNumber(ResA); a++){
	     	  Accept = AtomResidue(a,ResA);
		  if(HydrogenAtom(Accept)) continue;
		  if(!SideAtom(Accept)) continue;
		  if(IsAtom(" CB ",Accept) || IsAtom(" OH ",Accept)) continue;
		  res2 = ResAtom(Accept);
		  if(C==C2 && abs(res2 - res) < 1) continue;
		  if(P->seq[C2]) r2 = P->seq[C2][res2];  
		  else r2 = 0;
		  d_HA = DistanceAtoms(H,Accept);
		  if(d_HA > dmax) continue;
		  d_DA = DistanceAtoms(Donor,Accept);
		  if(d_DA <= d_HA) continue;
		  angle_DHA = CalcAngle(Donor,H,Accept);
		  // if(d_HA <= 3.5 && d_DA > d_HA && angle_DHA >= 90)
		  // if(d_HA <= 3.5)
		  // if(d_HA <= dmax)
		  // if(d_HA <= dmax && d_DA > d_HA && angle_DHA >= 90)
		  if(angle_DHA >= 90) {
// PutAtom(stderr,Accept);
		    Ave_d_DA += d_DA;
		    Ave_d_HA += d_HA;
		    if(d_near > d_HA){ d_near = d_HA; Near=Accept; }
		    if(d_far < d_HA){ d_far = d_HA; Far=Accept; }
		    // aveX += ...
		    numhits++;	// heavy atom hit;
		    if(numhits < 10) atom[numhits] = Accept;
		    x=InsertMheap((keytyp)d_HA,mH);	// i = 
		    if(x != 0) atom6[x]=Accept;
		  }
		}
		Int4 z;
		atm_typ Virtual=0;
		double	d_DV,d_HV,angle_DHV,d_HpX;
		if(numhits > 1){
		   angle_HNF = CalcAngle(H,Near,Far);
		   Ave_d_DA = Ave_d_DA/(double) numhits;
		   Ave_d_HA = Ave_d_HA/(double) numhits;
		}
		if(numhits > 4){
#if 1
		  // Add dheap or mheap of size six.
		  // use to compute centroid atom from both five and six member rings
		  // for trp. Currently missing some CH-pi bonds to five member ring...
		  // This will also allow us to use a larger value for dmax.
		  numhits = ItemsInMheap(mH);
		  for(a=0; !EmptyMheap(mH); ){ x=DelMinMheap(mH); a++; atom[a]=atom6[x]; } 
#endif
		   Virtual = AverageAtom(numhits,atom);
		   d_DV=DistanceAtoms(Donor,Virtual);
		   d_HV=DistanceAtoms(H,Virtual);
		   angle_DHV = CalcAngle(Donor,H,Virtual);
		   d_HpX = ProjectOnPlane(H,Virtual,atom[1],atom[2]);
// fprintf(stderr,"d_HpX = %.2f\n",d_HpX);
#if 1	// AFN: 9-30-2016; for SIPRES program...
		  minDist=MINIMUM(double,d_HV,minDist);
#endif
		   NilAtom(Virtual);
		}
		// < 90 --> within plane.
		    atm_typ	A,B;
		    Int4 rA,rB;
		    unsigned char rcA;
#if 0
		if((numhits >= 5 && (aa == 'H' || aa == 'W') ||
			numhits >= 6 && (aa == 'F' || aa == 'Y')) && angle_HNF < 80){
#else
		if((numhits >= 5 && (aa == 'H' || aa == 'W') ||
			numhits >= 6 && (aa == 'F' || aa == 'Y')) 
			&& d_HpX <= 1.2 && d_DV <= 4.5 && angle_DHV >= 120){
#endif

#if 0
		    if(fp){
			fprintf(fp,"Donor:\n"); PutAtom(fp,Donor);
		        fprintf(fp,"Hydrogen:\n"); PutAtom(fp,H);
		    	fprintf(fp,"Accept:\n");
		    	for(z=1; z <= numhits && z < 10; z++) PutAtom(fp,atom[z]);
		    }
#endif

		    A=Donor;B=Accept; rcA=r; 
		    rcA=r;
		    rA=ResAtom(A); rB=ResAtom(B);
#if 1	// fix problem with negative residue positions
	if(rA < 0) continue;
#endif
		    char str[9],*str_ptr = AtomName(A);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
		    } str[z]=0;
#if 1	// ***************** Find hydrogen too.. **********************************
		    Int4 zz;
		    str_ptr = AtomName(H);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
		    } str[z]=0;
#endif	//*************************************************************************
		    if(IsWaterAtom(A)){
			  if(fp) fprintf(fp,"#HOH%d\t// to ",rA);
		    } else {
			if(IsHeteroAtom(A)) if(fp) fprintf(fp,"// ");
			if(AtomChain(A) != ' '){ if(fp) fprintf(fp,"%c%d%c",AlphaChar(rcA,P->A),rA,AtomChain(A)); }
			else if(fp) fprintf(fp,"%c%d",AlphaChar(rcA,P->A),rA);
			if(strcmp("o",str)==0) strcpy(str,"c=o");
			if(str[0] == 'n'){ if(fp) fprintf(fp,"_%s.X\t// NH-pi: ",str); }
			else if(str[0] == 'c'){ if(fp) fprintf(fp,"_%s.X\t// CH-pi: ",str); }
			else { if(fp) fprintf(fp,"_%s.X\t// pi: ",str); }
		    }
		    if(AtomChain(B) != ' '){
			if(fp) fprintf(fp,"%s%d%c.%s %.2f(%.2f)A; HNF=%.1f,DHA=%.1f;DHV=%.f.\n",
				AtomResName(B),rB,AtomChain(B),AtomName(B),
				d_DV,d_HV,	// Ave_d_DA,Ave_d_HA,
				angle_HNF,angle_DHA,angle_DHV);
		    } else if(fp) fprintf(fp,
			    "%s%d.%s %.2f(%.2f) A; HNF=%.1f,DHA=%.1f;DHV=%.f..\n",
				AtomResName(B),rB,AtomName(B),
				d_DV,d_HV,  // Ave_d_DA,Ave_d_HA,
				angle_HNF,angle_DHA,angle_DHV);
		     n++;
		} else if(numhits >= 3){
#if 1	// *********** Find Donor and hydrogen for near hits too.. **************
	// NOTE!! I NEED TO MAKE THIS SECTION INTO A ROUTINE!!!
		    A=Donor;B=Accept; rcA=r; 
		    rcA=r;
		    rA=ResAtom(A); rB=ResAtom(B);
		    char str[9],*str_ptr = AtomName(A);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
		    } str[z]=0;
		    Int4 zz;
		    str_ptr = AtomName(H);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
		    } str[z]=0;
	// NOTE!! I NEED TO MAKE THIS SECTION INTO A ROUTINE!!!
#endif	//*************************************************************************
#if 0
		    fprintf(stderr,"Donor:\n"); PutAtom(stderr,Donor);
		    fprintf(stderr,"Hydrogen:\n"); PutAtom(stderr,H);
		    fprintf(stderr,"Accept:\n");
		    for(z=1; z <= numhits && z < 10; z++) PutAtom(stderr,atom[z]);
#endif
#if 0
			fprintf(fp,"// %c%d(%c%d_%s): near hit (%d)! (%.2f:%.2f A).\n",
				aa,res2,aa2,res,str,numhits,Ave_d_DA,Ave_d_HA);
#else
		if(numhits > 4){
		    if(fp){
			fprintf(fp,"#%c%d_%s.X // CH-pi near hit to %c%d ",aa2,res,str,aa,res2);
			fprintf(fp,"(%d) (%.2f:%.2f A); d_HpX = %.2f\n",numhits,Ave_d_DA,Ave_d_HA,d_HpX);
		    }
		} 
		// else fprintf(fp,"// %c%d(%c%d_%s): near hit (%d)! (%.2f:%.2f A).\n",
		//		aa,res2,aa2,res,str,numhits,Ave_d_DA,Ave_d_HA);
#endif
		}
	     }
	}
	NilMheap(mH);
	return n;
}

Int4	FindOtherPiHbondsPDB(res_typ ResD, res_typ ResA,Int4 C, Int4 C2, 
	float dmax, pdb_typ P)
{ return FindOtherPiHbondsPDB(stdout,ResD, ResA,C, C2, dmax,0,P); }

Int4	FindOtherPiHbondsPDB(FILE *fp,res_typ ResD, res_typ ResA,Int4 C, Int4 C2, 
	float dmax, char color, pdb_typ P)
{ double d=-1; return FindOtherPiHbondsPDB(fp,ResD, ResA,C, C2, dmax,0,P,d); }

#define MAX_NUM_ATOMS_FOR_VIRTUAL 30

Int4    FindOtherPiHbondsPDB(FILE *fp,res_typ ResD,res_typ ResA,Int4 C,Int4 C2,
        float dmax, char color, pdb_typ P,double &minDist)
// find potential aromatic pi-H interactions.
{
	Int4	n,a,d,r,r2,res,res2,numhits;
	atm_typ	H,Donor,Accept,atom[MAX_NUM_ATOMS_FOR_VIRTUAL + 3];
	double	d_HA,d_DA,angle_DHA,Ave_d_DA,Ave_d_HA;
	BooLean	na;
	FILE	*efp=P->fp_stderr;
	a_type	AB=P->A;
	BooLean	comments=FALSE;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("FindOtherPiHbondsPDB( ): input error"); 
	}
	if(efp == stderr) comments=TRUE; else comments=FALSE;
	comments=TRUE; // waiting to debug this option;
	if(ResD == ResA) return 0;
	char    aa = GetResidueAtom(AtomResidue(1,ResA),AB,P->nA,&na);
	aa = AlphaChar(aa,AB);
#if 1	// New 10/15/02: afn.
	char    aa2 = GetResidueAtom(AtomResidue(1,ResD),AB,P->nA,&na);
	aa2 = AlphaChar(aa2,AB);
#endif
	if(na || !(aa == 'N' || aa == 'Q' || aa == 'D' || aa == 'E' || aa == 'R')) return 0;
	atm_typ	Near,Far;
	double	d_near,d_far,angle_HNF;
	for(n=0,d = 1; d <= ResidueAtomNumber(ResD); d++){
	     if(!ResidueAtomHydrogens(d,ResD)) continue;
	     Donor = AtomResidue(d,ResD);
// PutAtom(stderr,Donor);
	     res=ResAtom(Donor);
	     if(P->seq[C]) r = P->seq[C][res]; else r = 0;
	     for(Int4 h=1; (H=ResidueAtomHydrogen(d,h,ResD)); h++){
		Ave_d_DA=Ave_d_HA=0.0;
// PutAtom(stderr,H);
		d_near=1000.0; d_far=0.0;
		Int4 numNO=0;
		atm_typ atomNO[MAX_NUM_ATOMS_FOR_VIRTUAL + 3];
	        for(numhits=0,a = 1; a <= ResidueAtomNumber(ResA); a++){
	     	  Accept = AtomResidue(a,ResA);
		  if(HydrogenAtom(Accept)) continue;
		  if(!SideAtom(Accept)) continue;
		  if(IsAtom(" CB ",Accept)) continue;
		  if((aa == 'E' || aa == 'Q' || aa == 'R') && IsAtom(" CG ",Accept)) continue;
		  if((aa == 'R') && IsAtom(" CD ",Accept)) continue;
		  res2 = ResAtom(Accept);
		  if(C==C2 && abs(res2 - res) < 1) continue;
		  if(P->seq[C2]) r2 = P->seq[C2][res2]; else r2 = 0;
		  d_HA = DistanceAtoms(H,Accept);
		  angle_DHA = CalcAngle(Donor,H,Accept);
		  d_DA = DistanceAtoms(Donor,Accept);
		  if(d_HA <= dmax && d_DA > d_HA) 
		  {
// PutAtom(stderr,Accept);
		    Ave_d_DA += d_DA; Ave_d_HA += d_HA;
		    if(d_near > d_HA){ d_near = d_HA; Near=Accept; }
		    if(d_far < d_HA){ d_far = d_HA; Far=Accept; }
		    // aveX += ...
		    numhits++;	// heavy atom hit;
		    if(numhits < MAX_NUM_ATOMS_FOR_VIRTUAL) atom[numhits] = Accept;
		    else {
			for(Int4 x=1; x <= numhits; x++) PutAtom(stderr,atom[x]);
			print_error("Input error: too many overlapping atoms in coordinate file");
		    } 
		    if(!CarbonAtom(Accept)){
		      numNO++;	// oxygen or nitrogen hit;
		      if(numNO < MAX_NUM_ATOMS_FOR_VIRTUAL) atomNO[numNO] = Accept;
		      else {
			for(Int4 x=1; x <= numNO; x++) PutAtom(stderr,atomNO[x]);
			print_error("Input error: too many overlapping atoms in coordinate file");
		      } 
		    }
		  }
		}
		if(numNO < 2) continue;
		Int4 z;
		atm_typ Virtual=0;
		double	d_DV,d_HV,angle_DHV,d_HpX;
		if(numhits > 1){
		   angle_HNF = CalcAngle(H,Near,Far);
		   Ave_d_DA = Ave_d_DA/(double) numhits;
		   Ave_d_HA = Ave_d_HA/(double) numhits;
		}
		if(numhits > 2){
		   Virtual = AverageAtom(numhits,atom);
		   d_DV=DistanceAtoms(Donor,Virtual);
		   d_HV=DistanceAtoms(H,Virtual);
		   angle_DHV = CalcAngle(Donor,H,Virtual);
		   d_HpX = ProjectOnPlane(H,Virtual,atomNO[1],atomNO[2]);
	// fprintf(stderr,"d_HpX = %.2f\n",d_HpX);
#if 1	// AFN: 9-30-2016; for SIPRES program...
		  minDist=MINIMUM(double,d_HV,minDist);
#endif
		}
		// < 90 --> within plane.
		if((numhits >= 3 && (aa == 'N' || aa == 'Q' || aa == 'D' || aa == 'E') ||
			numhits >= 4 && (aa == 'R')) 
			&& d_HpX <= 1.2 && d_DV <= 4.5 && angle_DHV >= 120){
		    atm_typ	A,B;
		    Int4 rA,rB;
		    unsigned char rcA;

#if 0
		    if(fp){
			fprintf(fp,"found H-bond (%d): (%.2f:%.2f A).\n",numhits,
				Ave_d_DA,Ave_d_HA);
		    	fprintf(fp,"Virtual:\n"); PutAtom(fp,Virtual);
		    	fprintf(fp,"Donor:\n"); PutAtom(fp,Donor);
		    	fprintf(fp,"Hydrogen:\n"); PutAtom(fp,H);
		   	fprintf(fp,"Accept:\n");
		        for(z=1; z <= numhits && z < 10; z++) PutAtom(fp,atom[z]);
		    }
#endif

		    A=Donor;B=Accept; rcA=r; 
		    rcA=r;
		    rA=ResAtom(A); rB=ResAtom(B);
#if 1	// fix negative position problem: AFN 8_26_2015.
	if(rA < 0) continue;
#endif
		    char str[9],*str_ptr = AtomName(A);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
				if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
				else str[z]=str_ptr[z];
		    } str[z]=0;
#if 1	// ***************** Find hydrogen too.. **********************************
		    Int4 zz;
		    str_ptr = AtomName(H);
		    if(isspace(str_ptr[0])) str_ptr++;
		    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
		    } str[z]=0;
#endif	//*************************************************************************
		    if(fp){
		     if(IsWaterAtom(A)){ fprintf(fp,"#HOH%d\t// to ",rA); }
		     else {
			 if(IsHeteroAtom(A)) fprintf(fp,"// ");
			 if(AtomChain(A) != ' '){ fprintf(fp,"%c%d%c",AlphaChar(rcA,P->A),rA,AtomChain(A)); }
			 else fprintf(fp,"%c%d",AlphaChar(rcA,P->A),rA);
			 if(strcmp("o",str)==0) strcpy(str,"c=o");
			 if(comments) fprintf(fp,"_%s.X\t// XH-pi to ",str);
			 else fprintf(fp,"_%s.X\n",str);
		     }
		     if(comments){
		        if(AtomChain(B) != ' '){
			   fprintf(fp,"\t// to %c%d%c (%d);d=%.2f(%.2f)A; HNF = %.1f,DHA = %.1f;DHV = %.f.\n",
				aa,rB,AtomChain(B),numhits, d_DV,d_HV,angle_HNF,angle_DHA,angle_DHV);
		        } else fprintf(fp, "\t// to %c%d (%d);d=%.2f(%.2f) A; HNF = %.1f,DHA = %.1f;DHV = %.f.\n",
				aa,rB,numhits,d_DV,d_HV,angle_HNF,angle_DHA,angle_DHV);
		      } else fprintf(fp,"\n");
		    }
		    n++;
		} else if(numhits >= 3){
		    // if(Virtual) NilAtom(Virtual);
#if 0
		    fprintf(stderr,"Donor:\n"); PutAtom(stderr,Donor);
		    fprintf(stderr,"Hydrogen:\n"); PutAtom(stderr,H);
		    fprintf(stderr,"Accept:\n");
		    for(z=1; z <= numhits && z < 10; z++) PutAtom(stderr,atom[z]);
#endif
#if 0
		    if(numhits > 2){	// PROBLEM WITH str!!
			fprintf(fp,
			 "// %c%d(%c%d_%s): near hit (%d)! (%.2f:%.2f A); d_HpX = %.2f\n",
				aa,res2,aa2,res,str,numhits,Ave_d_DA,Ave_d_HA,d_HpX);
		    }
#else
		    if(numhits > 2){
			if(fp) fprintf(fp,
			  "# %c%d.W // XH-pi near hit to %c%d (%d)! (%.2f:%.2f A); d_HpX=%.2f\n",
				aa,res2,aa2,res,numhits,Ave_d_DA,Ave_d_HA,d_HpX);
		    }
#endif
		} if(Virtual) NilAtom(Virtual);
	     }
	}
	return n;
}

void	TriadSuperImposePDB(FILE* fptr, pdb_typ ToP, pdb_typ FromP, pdb_typ P)
// Use 
{
	if(ToP->nchains != 1 || FromP->nchains != 1 || 
		ToP->natoms[1] != 3 || FromP->natoms[1] != 3){
		print_error("SuperImposeViaTriadPDB( ) input error");
	}
        Int4	i,C,a;
        sip_typ sip(ToP->atom[1],FromP->atom[1]);
	fprintf(fptr,"HEADER    %s\n",P->filename);
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		atm_typ atm = sip.MapTo(P->atom[C][a]);
		PutAtom(fptr, atm);
		NilAtom(atm);
	   }
	   fprintf(fptr,"TER               \n");
	}
	fprintf(fptr,"END\n");
}

void	TriadTwistChnPDB(FILE* fptr, pdb_typ ToP, pdb_typ FromP, Int4 chn, 
	Int4 twistpoint,pdb_typ P)
// Use 
{
	BooLean	reverse=FALSE;
	if(twistpoint < 0){ reverse=TRUE; twistpoint = -twistpoint; }
	if(ToP->nchains != 1 || FromP->nchains != 1 || 
		ToP->natoms[1] != 3 || FromP->natoms[1] != 3){
		print_error("TriadTwistChnPDB( ) triad input error");
	}
	if(chn < 1 || chn > P->nchains) 
		print_error("TriadTwistChnPDB( ) chain input error");
	if(twistpoint < MinResPDB(chn,P) || twistpoint > MaxResPDB(chn,P)) 
		print_error("TriadTwistChnPDB( ) twistpoint input error");
        Int4	i,C,a;
        sip_typ sip(ToP->atom[1],FromP->atom[1]);
	fprintf(fptr,"HEADER    %s\n",P->filename);
	for(C=1; C <= P->nchains; C++){
	   for(a=1; a<= P->natoms[C]; a++){
		atm_typ atm0 = P->atom[C][a];
	     if(C != chn) PutAtom(fptr, atm0);
	     else {
		if(!IsAminoAcidAtom(atm0)) PutAtom(fptr, atm0);
		else if((!reverse && ResAtom(atm0) >= twistpoint) || 
			(reverse && ResAtom(atm0) <= twistpoint)){
		  atm_typ atm = sip.MapTo(P->atom[C][a]);
		  // PutAtom(stderr,atm);
		  PutAtom(fptr, atm);
		  NilAtom(atm);
		} else PutAtom(fptr, atm0);
		// if(IsAminoAcidAtom(atm0)) fprintf(stderr,"res = %d\n",ResAtom(atm0));
	     }
	   }
	   fprintf(fptr,"TER               \n");
	}
	fprintf(fptr,"END\n");
}

char	**AdjacentHeteroMoleculePDB(Int4 C0, double maxdist,pdb_typ P)
// Return string array indicating which molecules are are adjacent to chain C0
// sprintf(str,"[%s]%d%c",res_name,res,chain);
{
        Int4    C,a,a0,i,j;
        atm_typ atm,atm0;
        char	**OutPut,str[100],*mole_name;
	Int4	NumMol=0,MaxNumMol=1000;

        // find subunits...
        NEWP(OutPut,MaxNumMol+3,char);
        for(C=1; C <= nChainsPDB(P); C++){
           if(C == C0){ continue; }
	   if(IsProteinPDB(C,P)) continue;
           for(a=1; a<= NumberAtomsPDB(C,P); a++){
             atm = AtomPDB(C,a,P);
	     if(IsWaterAtom(atm)) continue;
             if(IsHeteroAtom(atm)){
               for(a0=1; a0<= NumberAtomsPDB(C0,P); a0++){
                atm0 = AtomPDB(C0,a0,P);
                if(DistanceAtoms(atm, atm0) <= maxdist){
			mole_name=AtomResName(atm);
			while(isspace(mole_name[0])){
			   mole_name++;
			   if(!isprint(mole_name[0])) print_error("AdjacentHeteroMoleculePDB() error");
			}
			sprintf(str,"[%s]%d%c",mole_name,ResAtom(atm),AtomChain(atm));
			// sprintf(str,"[%s]%d%c",AtomResName(atm),ResAtom(atm),AtomChain(atm));
			BooLean is_new=TRUE;
			for(i=1; i <= NumMol; i++){
				if(strcmp(str,OutPut[i]) == 0){ is_new=FALSE; break; }
			}
			if(is_new){
				NumMol++;
				OutPut[NumMol]=AllocString(str);  
				if(NumMol >= MaxNumMol) return OutPut;
			} break;
		}
               } 
	     }
           }
        } return OutPut;
}

char	**AdjacentHeteroMoleculePDB2(Int4 C0, double maxdist,pdb_typ P)
// Return string array indicating which molecules are are adjacent to chain C0
// sprintf(str,"[%s]%d%c",res_name,res,chain);
{
        Int4    C,a,a0,i,j;
        atm_typ atm,atm0;
        char	**OutPut,str[100];
	Int4	NumMol=0,MaxNumMol=1000;

        // find subunits...
        NEWP(OutPut,MaxNumMol+3,char);
        for(C=1; C <= P->nchains; C++){
           if(C == C0){ continue; }
	   if(IsProteinPDB(C,P)) continue;
           for(a=1; a<= P->natoms[C]; a++){
             atm = P->atom[C][a];
	     if(IsWaterAtom(atm)) continue;
             if(IsHeteroAtom(atm)){
               for(a0=1; a0<= P->natoms[C0]; a0++){
                atm0 = P->atom[C0][a0];
                if(DistanceAtoms(atm, atm0) <= maxdist){
			sprintf(str,"[%s]%d%c",AtomResName(atm),ResAtom(atm),AtomChain(atm));
			BooLean is_new=TRUE;
			for(i=1; i <= NumMol; i++){
				if(strcmp(str,OutPut[i]) == 0){ is_new=FALSE; break; }
			}
			if(is_new){
				NumMol++;
				OutPut[NumMol]=AllocString(str);  
				if(NumMol >= MaxNumMol) return OutPut;
			} break;
		}
               } 
	     }
           }
        } return OutPut;
}

