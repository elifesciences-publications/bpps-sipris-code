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

//***************************** KANTI CODE *********************************
#define MIN_angleDHA 90

Int4	HBondAnglesPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
		unsigned short file_id, pdb_typ P)
// Too redundant with other code; need to condense this at some point!
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,a3,r,r2,res,gly,res2;
	atm_typ	H,Donor,Accept,X;
	BooLean	na;
	FILE	*efp=P->fp_stderr;
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1){
		pdb_error("FindBondsPDB( ): input error"); 
	}
	gly = AlphaCode('G',P->A);
	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;
	     Donor = AtomResidue(a,Res);
	     res=ResAtom(Donor);
	     if(P->maxres[C] < res) r = 0;
	     else {
		if(P->seq[C]) r = P->seq[C][res];  
	     	else r = 0;
	     }
	     // Find all Donor atoms:
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
	      for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	        Accept = P->atom[C2][a2];
		if(HydrogenAtom(Accept)) continue;
		if(CarbonAtom(Accept)) continue;
		if(Accept == Donor) continue;
		res2 = ResAtom(Accept);
		if(Res2 && Res2 != res2) continue;
	     	if(P->seq[C2] && IsAminoAcidAtom(Accept)) r2 = P->seq[C2][res2];  
		else r2=0;
		if(C==C2 && abs(res2 - res) < 2) continue;
#if 1	// Fix purify array bounds error....
		char	IsRes;
		IsRes= GetResidueAtom(Accept, P->A, P->nA, &na);
		// Is this an amino acid residue?
#endif
		for(a3 = 1; a3 <= P->natoms[C2]; a3++){
	          X = P->atom[C2][a3];
		  if(Accept == X) continue;
		  if(HydrogenAtom(X)) continue;
		  if(OxygenAtom(X)) continue;

		  double d_AX = DistanceAtoms(X,Accept);
		  if(d_AX > 1.7) continue;

		  double angle_HAX = CalcAngle(H,Accept,X);
		  if(angle_HAX < MIN_angleDHA) continue;

		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MIN_angleDHA){
			fprintf(fp,"%c-%c..%c(-%c):",
				CharacterAtom(Donor),CharacterAtom(H),
				CharacterAtom(Accept),CharacterAtom(X));
			fprintf(fp," DA=%.2f; HA=%.2f; DHA=%.1f; HAX=%.1f.\n",
				CharacterAtom(Accept),d_DA,d_HA,angle_DHA,angle_HAX);
		      n++;
		  }
		} 
	      }
	     }
	} return n;
	// NOTE: switching acceptor and donor violates the geometry!!! need to fix this...
}


void	PutAngularPDB(FILE *fptr, Int4 file_id, Int4 resStart,Int4 resEnd, float HA_dmax, 
		float dmax,Int4 C, Int4 C2, char chain, BooLean *skip, char color, pdb_typ P)
{
	Int4 res0;
	Int4 r,i,j,num_resA,num_resD;
	if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
	if(C==0){ 
		     fprintf(stderr,"chain = %c\n",chain);
		     print_error("chain not found!");
	}
	if(chain==0) chain='?';
	// Int4 Start,End;
	res_typ *ResCA,*ResOA,*ResCD,*ResOD,**ResALL;
	Int4	num_resC,num_resO,*num_resALL;
	NEWP(ResALL,nChainsPDB(P)+3,res_typ);
	NEW(num_resALL,nChainsPDB(P)+3,Int4);
	for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	  ResALL[C2] = MakeResPDB(C2,&num_resALL[C2],P);
	}
	num_resC=num_resALL[C]; ResCA=ResCD=ResALL[C];

	// Move this routine to libpdb.a at some point:
	a_type AB = AminoAcidAlphabetPDB(P);
	Int4	MaxNeighbors=50; // maximum residue neighbors to look at
	res_typ	*Contact;
	double	contact_surface,*Ratio;
	NEW(Contact,MaxNeighbors+3,res_typ);
	NEW(Ratio,MaxNeighbors+3,double);

	for(res0=resStart; res0 <= resEnd; res0++){
	     if(skip && skip[res0]) continue;	// skip unconserved residues...
	     Int4	hpsz=0;
	     dh_type dH = dheap(MaxNeighbors,4);

	     // fprintf(fptr,"\n");
	     // find H-bonds...
	     for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	      ResOD=ResOA=ResALL[C2];	// other than query residue.
	      num_resO=num_resALL[C2];
	      
	      for(j=1; j <= num_resC; j++){
		if(ResidueID(ResCD[j]) == res0) {
		   if(skip){
		     if(C == C2){
			for(Int4 res2=1; res2 <=num_resO; res2++){
			  if(!skip[res2]) 
				HBondAnglesPDB(fptr,ResCD[j],C,C2,HA_dmax,res2,file_id,P); 
			}
		     }
		   } else {
			HBondAnglesPDB(fptr,ResCD[j],C,C2,HA_dmax,0,file_id,P);
		   }
		}
	      }
	     }
	 } //********** end of res0 loop... *****
	 for(C2=1; C2 <= nChainsPDB(P); C2++){ 
		for(i=1; i <= num_resALL[C2]; i++) NilRes(ResALL[C2][i]);
		free(ResALL[C2]);
	 } free(ResALL); free(num_resALL);
} 
 
double	PrintThetaDist(double &theta,FILE *ofp,atm_typ Xd,atm_typ D,atm_typ A,atm_typ Xa,
		double d_HA,double angle_DHA)
// print the angle theta and the distance between 
{
	static double PI=3.14159265358979323846;

	double	dx=AtomX(D) - AtomX(A);
	double	dy=AtomY(D) - AtomY(A);
	double	dz=AtomZ(D) - AtomZ(A);
	atm_typ V=VirtualAtom(dx+AtomX(Xa),dy+AtomY(Xa),dz+AtomZ(Xa));
	double	angle=(PI/180.0)*CalcAngle(Xd,D,V),aDHA=(PI/180.0)*angle_DHA;
	double	distance=DistanceAtoms(D,A);
	NilAtom(V); theta=angle;
	double X = d_HA*sin(2.0*aDHA),Y=d_HA*cos(2.0*aDHA);
	double x = distance*sin(2.0*theta),y=distance*cos(2.0*theta);
#if 1
	fprintf(ofp,"%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
		distance, angle, x, y, d_HA, aDHA, X, Y); 
#else
	if(angle >= 1.7 && angle <= 2.3) fprintf(ofp,"%.3f\t%.3f\t%.3f\n",angle,distance,d_HA); 
#endif
	return distance;
}

Int4	GetHBondThetaAndDistPDB(FILE *fp,FILE *ofp,res_typ Res,Int4 C,Int4 C2,float dmax,res_typ Res2,
		unsigned short file_id, pdb_typ P, char mode, char *donor, char *accept,
		h_type tHG,h_type dHG,h_type dhaHG,h_type haHG)
// Too redundant with other code; need to condense this at some point!
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	x,n,a,a2,r,r2,res,res2;
	atm_typ	H,Donor,Accept,Xd,Xa;
	BooLean	na;
	FILE	*efp=P->fp_stderr;
	BooLean	comments=FALSE;
	static double PI=3.14159265358979323846;
	a_type AB = AminoAcidAlphabetPDB(P);
	
	if(C > P->nchains || C<1 || C2 >P->nchains || C2<1) pdb_error("FindBondsPDB( ): input error"); 
	if(efp==stderr) comments=TRUE; else comments=FALSE;
	comments=TRUE; // waiting to debug this option;
	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;		// Does atom lack hydrogens? Then skip it.
	     Donor = AtomResidue(a,Res);
	     // if(!IsAminoAcidAtom(Donor)) continue;	// Skip if not a bona fide amino acid atom.
	     res=ResAtom(Donor);
	     if(P->maxres[C] < res) r = 0;
	     else { if(P->seq[C]) r = P->seq[C][res];  else r = 0; }
	     // Find all Donor atoms:
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
	        for(a2 = 1; a2 <= ResidueAtomNumber(Res2); a2++){
	          Accept = AtomResidue(a2,Res2);
		  if(HydrogenAtom(Accept)) continue;
#if 1	// Do I want to see backbone-to-backbone interactions?
		  if(IsAminoAcidAtom(Donor) && IsAminoAcidAtom(Accept) &&
			!SideAtom(Donor) && !SideAtom(Accept)) continue;
#endif
		  if(Accept == Donor) continue;

char AA1=GetCharResidue(Res,AB),AA2=GetCharResidue(Res2,AB);
if(donor && strchr(donor,AA1) == 0) continue;  // inacceptable amino acid.
if(accept && strchr(accept,AA2) == 0) continue;  // inacceptable amino acid.
if(mode == 'C'){
	if(!CarbonAtom(Donor) || ! OxygenAtom(Accept)) continue; 
} else if(mode == 'N'){
	if(!NitrogenAtom(Donor) || ! OxygenAtom(Accept)) continue; 
#if 0
#if 0
	if(AA1 != 'R' || AA2 != 'E') continue;
	if(strcmp(AtomName(Donor)," NH2") != 0 && strcmp(AtomName(Donor)," NH1") != 0) continue;
#elif 0
	if(AA1 != 'K' || (AA2 != 'E' && AA2 != 'D')) continue;
#else
	if(AA1 != 'Q' || (AA2 != 'E' && AA2 != 'D')) continue;
#endif
#endif
} else if(mode == 'O'){
	if(!OxygenAtom(Donor) || ! OxygenAtom(Accept)) continue; 
} else if(mode == 'S'){
	if(!SulfurAtom(Donor) || ! OxygenAtom(Accept)) continue; 
	// if(AA1 != 'S') continue;
}

		  res2 = ResAtom(Accept);
		  // if(Res2 && Res2 != res2) continue;	// 10-7-09: afn bug fix...
	     	  if(P->seq[C2] && IsAminoAcidAtom(Accept)) r2 = P->seq[C2][res2];  
		  else r2=0;
		  if(C==C2 && abs(res2 - res) < 2) continue;
#if 1	// Fix purify array bounds error....
		  char	IsRes;
		  IsRes=GetResidueAtom(Accept,P->A,P->nA,&na);
		  // Is this an amino acid residue?
#endif
		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  // if(d_HA <= 2.5 && d_DA > d_HA && angle_DHA >= 90)
		  // if(d_HA <= (double) dmax && d_DA > d_HA)
		  if(d_HA < 1.2) continue;  // this is an artifact!!!
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MIN_angleDHA)
		  {
		    Int4 rA,rB,isaccept=0;
		    atm_typ	A=Donor,B=Accept; 
		    unsigned char rcA=r;
		    BooLean NeitherCarbons=TRUE;
		    if(CarbonAtom(A) || CarbonAtom(B)) NeitherCarbons=FALSE;
#if 1		    // output for Kanti Mardia.
		    if(TRUE || !CarbonAtom(A) && !CarbonAtom(B)){
			x=ResidueNumBoundHeavyAtoms(a,Res); // a == donor...
			if(x <= 0){
				fprintf(stderr,"---------------------\n");
				PutAtom(stderr,A); PutAtom(stderr,B);
				fprintf(stderr,"x = %d; d_HA=%.3f\n",x,d_HA);
				// PutBondsInRes(stderr,Res);
				PutBondHeavyInRes(stderr,a,Res);
				// PrintResidueAtoms(stderr,Res);
				// PrintResidueAtoms(stderr,Res2);
	continue;	// do this for now...
				assert(x > 0);	// must be covalently bound to something...
			}
			Xd=ResidueBoundHeavyAtom(a,1,Res);
			if(mode != 'C') assert(CarbonAtom(Xd));
			x=ResidueNumBoundHeavyAtoms(a2,Res2);
			if(x > 0)
			// if(x > 0 && d_HA < 1.6)
			{	// must be covalently bound to something...
			  Xa=ResidueBoundHeavyAtom(a2,1,Res2);
			  double	d_Xa_D = DistanceAtoms(Xa,Donor);
			  double	d_Xd_A = DistanceAtoms(Accept,Xd);
		  	  double angle_HAXa = CalcAngle(H,Accept,Xa);
			  if(d_Xa_D > d_DA && d_Xd_A > d_DA && angle_HAXa >= MIN_angleDHA){
			    assert(CarbonAtom(Xa));
			    double theta,d=PrintThetaDist(theta,ofp,Xd,Donor,Accept,Xa,d_HA,angle_DHA);
	if(angle_DHA < MIN_angleDHA){
		fprintf(stderr,"angle_DHA = %.3f\n",angle_DHA);
		assert(angle_DHA >= MIN_angleDHA);
	}
			    IncdHist(d,dHG); IncdHist(theta,tHG); 
			    IncdHist(d_HA,dhaHG); IncdHist((PI/180.0)*angle_DHA,haHG); 
			  }
			}
		    }
#endif
#if 0
#endif
		    if(0) for(rcA=r; isaccept < 2; A=Accept,B=Donor,rcA=r2,isaccept++){
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
			 else if(str[0]=='c' && strstr(str,"-") == 0){ fprintf(fp,"// "); }
			 if(NeitherCarbons){
				if(file_id > 0) fprintf(fp," %d",file_id); else fprintf(fp," ");
			 } else fprintf(fp," ");
			 if(IsHeteroAtom(A)){ fprintf(fp,"!"); }
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
		     } n++; // end of for loop.
		  } 
		}
	     }
	} return n;
}

//***************************** END KANTI CODE *********************************

void	PutContactsPDB(FILE *fptr, Int4 file_id, Int4 resStart,Int4 resEnd, float HA_dmax, 
		float dmax,Int4 C, Int4 C2, char chain, BooLean *skip, char color, pdb_typ P)
{
	Int4 res0;
	Int4 r,i,j,num_resA,num_resD;
	if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
	if(C==0){ 
		     fprintf(stderr,"chain = %c\n",chain);
		     print_error("chain not found!");
	}
	if(chain==0) chain='?';
	// Int4 Start,End;
	res_typ *ResCA,*ResOA,*ResCD,*ResOD,**ResALL;
	Int4	num_resC,num_resO,*num_resALL;
	NEWP(ResALL,nChainsPDB(P)+3,res_typ);
	NEW(num_resALL,nChainsPDB(P)+3,Int4);
	for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	  ResALL[C2] = MakeResPDB(C2,&num_resALL[C2],P);
	}
	num_resC=num_resALL[C]; ResCA=ResCD=ResALL[C];

	// Move this routine to libpdb.a at some point:
	a_type AB = AminoAcidAlphabetPDB(P);
	Int4	MaxNeighbors=50; // maximum residue neighbors to look at
	res_typ	*Contact;
	double	contact_surface,*Ratio;
	NEW(Contact,MaxNeighbors+3,res_typ);
	NEW(Ratio,MaxNeighbors+3,double);

	for(res0=resStart; res0 <= resEnd; res0++){
	     if(skip && skip[res0]) continue;	// skip unconserved residues...
	     Int4	hpsz=0;
	     dh_type dH = dheap(MaxNeighbors,4);

     	    //********** Check out surface VDW contacts of res0 to nearby residues.
	     fprintf(fptr,"\n");
	  if(!skip){ 
	     fprintf(fptr,"\n//*********** ResidueToResidueContacts (%d): ***********\n",res0);
	     for(j=1; j <= num_resALL[C]; j++){
	       if(ResidueID(ResCD[j]) == res0) {
		double total_ratio = 0.0;
		BooLean	SideChainOnly;
		char c_res = GetCharResidue(ResCD[j],AB);
		if(c_res == 'G') SideChainOnly=FALSE;
		else SideChainOnly=TRUE;
		for(C2=1; C2 <= nChainsPDB(P); C2++){ 
		   res_typ *ResC2 = ResALL[C2];
		   if(num_resALL[C2] > 0 && !isspace(ChainResidue(ResC2[1])))
			   fprintf(fptr,"//======== Chain %c ==========\n",
					   	ChainResidue(ResC2[1]));
	           for(hpsz=0,i=1; i <= num_resALL[C2]; i++){
		   	if(ResCD[j] == ResC2[i]) continue;
#if 1	// This removes non conserved residues from consideration:
			if(skip && C == C2){
				Int4 res2=ResidueID(ResC2[i]);
				if(skip[res2]) continue;
			}
#endif
			double ratio;
	     		contact_surface = 
				ResidueToResidueContact(ResCD[j],ResC2[i],&ratio,SideChainOnly);
			if(contact_surface > 0.0){
				hpsz++;
				insrtHeap(hpsz,(keytyp)-contact_surface,dH);
				Ratio[hpsz]=ratio; total_ratio+= ratio;
				Contact[hpsz]=ResC2[i];
			}
		   }
		   Int4 hp_pos=1;
		   while(!emptyHeap(dH)){
			contact_surface = (double) -minkeyHeap(dH);
			assert((i=delminHeap(dH)) != 0);
			char    c2_res = GetCharResidue(Contact[i],AB);
			if(!isspace(ChainResidue(Contact[i]))){
			   fprintf(fptr,
				"#%c%d%c.X  // %c%d%c %.1f A^2 contact (%1.f%%:%d).\n",
					c2_res,ResidueID(Contact[i]),ChainResidue(Contact[i]),
					c_res,ResidueID(ResCD[j]),ChainResidue(ResCD[j]),
					contact_surface,Ratio[i]*100.0,hp_pos);
			} else {
			   fprintf(fptr,"#%c%d.X  // %c%d %.1f A^2 contact (%1.f%%:%d).\n",
					c2_res,ResidueID(Contact[i]),
					c_res,ResidueID(ResCD[j]),
					contact_surface,Ratio[i]*100.0,hp_pos);
			} hp_pos++;
			// PrintResidueSurface(fptr,ResCD[j]);
		   } hpsz=0;
		}
		Nildheap(dH);
		fprintf(fptr,"# total percent surface = %.1f\n",total_ratio*100.0);
		break;
	       }
	     }
	} // end if(!skip)...
	     // find H-bonds...
	     for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	      ResOD=ResOA=ResALL[C2];	// other than query residue.
	      num_resO=num_resALL[C2];
	      
	      if(!skip) fprintf(fptr,
			"//=========== H-bond donor (%d%c): ===========\n",res0,chain);
	      for(j=1; j <= num_resC; j++){
		if(ResidueID(ResCD[j]) == res0) {
		   if(skip){
		     if(C == C2){
			for(Int4 res2=1; res2 <=num_resO; res2++){
			  if(!skip[res2]) 
				FindHBondsPDB(fptr,ResCD[j],C,C2,HA_dmax,res2,color,file_id,P); 
			}
		     }
		   } else {
			// PrintResidueAtoms(stdout,ResCD[j]); // ResCD == query.
			// PutBondsInRes(stdout,ResCD[j]);
			FindHBondsPDB(fptr,ResCD[j],C,C2,HA_dmax,0,color,file_id,P);
		   }
		}
	      }
	      if(!skip) fprintf(fptr,
			"//=========== H-bond acceptor (%d%c): ===========\n",res0,chain);
	      for(j=1; j <= num_resO; j++){
//NOTE: !!!!!!!!!!! NEED TO modify FindHBondsPDB so that it can skip unconserved!!!!
#if 1	// This removes non conserved residues from consideration:
		if(skip){
		  if(C == C2){
			Int4 res2=ResidueID(ResOD[j]); 
			if(!skip[res2]) FindHBondsPDB(fptr,ResOD[j],C2,C,HA_dmax,res0,
						color,file_id,P); 
		  }
		  // 	FindHBondsPDB(fp,Res,C,C2,HA_dmax,Res2,0,P);
		} else FindHBondsPDB(fptr,ResOD[j],C2,C,HA_dmax,res0,color,file_id,P); // res0 is acceptor.
#endif
	      }

	      if(!skip) fprintf(fptr,"//=========== pi-bonds: ===========\n");
	      for(i=1; i <= num_resC; i++){
		if(ResidueID(ResCD[i]) == res0){	// --> ResCA[i] == res0
	         for(j=1; j <= num_resO; j++){
#if 1	// This removes non conserved residues from consideration:
		  if(skip){
		    Int4 res2;
		    if(C == C2){
		      res2=ResidueID(ResOA[j]); 
		      if(!skip[res2]){
		       FindAromaticHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		       FindOtherPiHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		      }
		      res2=ResidueID(ResOD[j]); 
		      if(!skip[res2]){
		       FindAromaticHbondsPDB(fptr,ResOD[j],ResCA[i],C2,C,dmax,color,P);
		       FindOtherPiHbondsPDB(fptr,ResOD[j],ResCA[i],C2,C,dmax,color,P);
		      }
		    }
#endif
		  } else {
		   // PrintResidueAtoms(fptr,ResA[i]);
		   FindAromaticHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		   FindAromaticHbondsPDB(fptr,ResOD[j],ResCA[i],C2,C,dmax,color,P);
		   FindOtherPiHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		   FindOtherPiHbondsPDB(fptr,ResOD[j],ResCA[i],C2,C,dmax,color,P);
		  }
		 }
		}
	      }
	     }
	 } //********** end of res0 loop... *****
	 for(C2=1; C2 <= nChainsPDB(P); C2++){ 
		for(i=1; i <= num_resALL[C2]; i++) NilRes(ResALL[C2][i]);
		free(ResALL[C2]);
	 } free(ResALL); free(num_resALL);
} 

Int4	FindQuadWatersPDB(float dmax, pdb_typ P)
// print out 
{
	Int4	a,b,C,D,N,Total=0,res;
	atm_typ	atmA,atmB,atom[103],water[103];
	float	dist[103];
	BooLean	na;
	
	for(C=1; C <= P->nchains; C++){
	   for(atmB=0,a=1; a <= MaxAtomsPDB(C,P); a++){
	     atmA = GetAtomPDB(P,C,a);
	     if(IsWaterAtom(atmA)){
	       N=0;
	       for(D=1; D <= P->nchains; D++){
	        for(atmB=0,b=1; b <= MaxAtomsPDB(D,P); b++){
		  atmB = GetAtomPDB(P,D,b);
		  if(IsAminoAcidAtom(atmB) && !SideAtom(atmB) && OxygenAtom(atmB) 
			&& DistanceAtoms(atmA, atmB) <= dmax){
			N++; atom[N]=atmB; dist[N]=DistanceAtoms(atmA, atmB);
			if(N > 100) print_error("FindQuadWatersPDB( ) error: too many hits");
		  }
		}
	       }
	       if(N >= 4){
		   if(Total==0){
			fprintf(stdout,"%s:\n",FilenamePDB(P));
		   }
		   fprintf(stdout,"%s.%s %d%c: %d hits\n",
                                AtomResName(atmA),AtomName(atmA),ResAtom(atmA),AtomChain(atmA),N);
		   for(Int4 i=1; i <= N; i++){
			res = ResAtom(atom[i]);
			char  aa=GetResidueAtom(atom[i], P->A, P->nA, &na);
			fprintf(stdout,"%c%d%c_%s.X %.2f A.\n",
                                AlphaChar(aa,P->A),res,AtomChain(atom[i]),AtomLowerName(atom[i]),dist[i]);
			// PutAtom(stdout,atom[i]);
		   } Total++;
	       }
	      }
	   }
	}
	return Total;
}


