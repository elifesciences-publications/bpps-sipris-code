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

#include "atom.h"
#include "pdb.h"
#include "alphabet.h"
#include "histogram.h"

// atm_typ AverageAtom(Int4 natom, atm_typ *atm);
// Use this to compute centroid for aromatic residue.

// float   DistanceAtoms(register atm_typ A1, register atm_typ A2);

//*********************************************************************
// Table 1 from Burley & Petsko (1985) Science 229: 23-28.
// 		Element:	a	b	c
static double	abcForH[3]={ 11.66, 108.06, 1.87 };
static double	abcForC[3]={ 49.13, 606.01, 1.80 };
static double	abcForN[3]={ 37.12, 504.51, 1.89 };
static double	abcForO[3]={ 33.52, 479.65, 1.98 };

static char	PheAroAtoms[11][5]={"CG","CD1", "HD1","CD2", "HD2","CE1","HE1",
						"CE2", "HE2", "CZ", "HZ" };
static double	PheDelta_q[11]   = {0.0,-0.153,0.153,-0.153,0.153,-0.153,0.153,
						-0.153,0.153,-0.153,0.153};

static char	TyrAroAtoms[13][5]={"CG","CD1","HD1","CD2","HD2","CE1","HE1",
						"CE2","HE2","CZ","OH","HH"};
static double	TyrDelta_q[13]   = {0.0,-0.153,0.153,-0.153,0.153,-0.153,0.153,
						-0.153,0.153,0.218,-0.655,0.437};

static char	TrpAroAtoms[15][5]={"CG","CD1","HD1","CD2","NE1","HE1","CE2",
				"CE3","HE3","CZ2","HZ2","CZ3","HZ3","CH2","HH2"};
static double	TrpDelta_q[15]={0.0,-0.011,0.153,-0.055,-0.540,0.440,0.142,
				-0.174,0.153,-0.164,0.153,-0.164,0.153,-0.197,0.153};

//*********************************************************************



double	KitaigorodskySum(Int4 Nj, Int4 Nk, double *a_j, double *a_k, double *b_j, 
		double *b_k, double *c_j, double *c_k,double *q_j, 
		double *q_k, double **r)
{
	double	E=0;
	Int4	j,k;

	for(j=1; j <= Nj; j++){
	    for(k=1; k <= Nk; k++){
		E += b_j[j]*b_k[k]*exp(-(c_j[j] + c_k[k])*r[j][k]);
		E -= pow(a_j[j]*a_k[k]*r[j][k],-6.0);
		E += 1.0/(q_j[j]*q_k[k]*r[j][k]);
	    }
	} return E;
}

double	Kitaigorodsky(double a_j, double a_k, double b_j, 
		double b_k, double c_j, double c_k,double q_j, 
		double q_k, double r_jk)
{
	return (b_j*b_k*exp(-(c_j + c_k)*r_jk) - pow(a_j*a_k*r_jk,-6.0) 
				+ 1.0/(q_j*q_k*r_jk));
}

BooLean GetKitaigorodskyParameters(atm_typ atom, double *a, double *b, 
	double *c, double *q, a_type AB, a_type nAB)
{
	BooLean na;
	char aa=GetResidueAtom(atom, AB, nAB, &na); aa = AlphaChar(aa,AB);
	if(na) return FALSE;
	if(CarbonAtom(atom) && AtomIndex(atom)=='B') return FALSE; // CB
	switch (aa){
	  case 'F':
	    if(CarbonAtom(atom)){ *a=49.13; *b=606.01; *c=1.80; 
		if(AtomIndex(atom)=='G') *q=0.0;	// CG
		else *q=-0.153; 			// CD1, CD2, CE1, CE2, CZ
	    } else if(HydrogenAtom(atom)){ *a=11.66; *b=108.06; *c=1.87; *q=0.153; }
	    else return FALSE;
	    break;
	  case 'Y':
	    if(CarbonAtom(atom)){ *a=49.13; *b=606.01; *c=1.80; 
		if(AtomIndex(atom)=='G') *q=0.0; 
		else if(AtomIndex(atom)=='Z') *q=0.218; 
		else *q=-0.153; 
	    } else if(HydrogenAtom(atom)){ *a=11.66; *b=108.06; *c=1.87;
		if(AtomIndex(atom)=='H') *q=0.437; else *q=0.153; 
	    } else if(OxygenAtom(atom)){ *a=33.52; *b=479.65; *c=1.98; *q=-0.655; }
	    else return FALSE;
	    break;
	  case 'W':
	    if(CarbonAtom(atom)){ *a=49.13; *b=606.01; *c=1.80; 
		switch(AtomIndex(atom)){
		  case 'G': *q=0.0; break;
		  case 'D': 
		        if(AtomIndexDigit(atom) == '1') *q=-0.011; 
			else if(AtomIndexDigit(atom) == '2') *q=-0.055; 
			else return FALSE;
		   break;
		  case 'E':
		        if(AtomIndexDigit(atom) == '2') *q=0.142; 
			else if(AtomIndexDigit(atom) == '3') *q=-0.174; 
			else return FALSE;
		   break;
		  case 'Z': *q=-0.164; break;
		  case 'H':
		        if(AtomIndexDigit(atom) == '2') *q=-0.197; 
			else return FALSE;
		   break;
		  default: return FALSE;
		}
	    } else if(NitrogenAtom(atom)){ *a=37.12; *b=504.51; *c=1.89;
		if(AtomIndex(atom)=='E') *q=-0.540; 	// NE1
		else return FALSE;
	    } else if(HydrogenAtom(atom)){ *a=11.66; *b=108.06; *c=1.87;
		if(AtomIndex(atom)=='E'){
		  if(AtomIndexDigit(atom) == '1') *q=0.440; // HE1
		  else if(AtomIndexDigit(atom) == '3') *q=0.153; // HE3
		  else return FALSE;
		} else if(AtomIndex(atom)=='Z') *q=0.153; // HZ2, HZ3
		else if(AtomIndex(atom)=='H') *q=0.153; // HH2
		else return FALSE;
	    } else return FALSE;
	    break;
	  default:
	    return FALSE;
	}
	return TRUE;
}

BooLean	GetKitaigorodsky(double *E, atm_typ atomJ, atm_typ atomK,a_type AB,
	a_type nAB)
{
	double a_j,a_k,b_j,b_k,c_j,c_k,q_j,q_k,r_jk;
	double term1,term2,term3;
	
	if(!GetKitaigorodskyParameters(atomJ,&a_j,&b_j,&c_j,&q_j,AB,nAB)) 
		return FALSE;
	if(!GetKitaigorodskyParameters(atomK,&a_k,&b_k,&c_k,&q_k,AB,nAB)) 
		return FALSE;
	r_jk= DistanceAtoms(atomK,atomJ);
	term1 = b_j*b_k*exp(-(c_j + c_k)*r_jk);
	term2 = -a_j*a_k*pow(r_jk,-6.0);
	term3 = q_j*q_k/r_jk;
	*E= term1 + term2 + term3;
	// fprintf(stderr,"E = %g + %g + %g = %g\n",term1,term2,term3,*E);
	// fprintf(stderr,"q_j = %g; q_k = %g; r_jk = %g\n",q_j,q_k,r_jk);
	return TRUE;
}

BooLean	BurleyAromatic(FILE *fp,Int4 file_id, double *Energy,res_typ ResJ,
	res_typ ResK, char chainJ, char chainK, char colorJ, char colorK, pdb_typ P)
{
	BooLean na;
	double	E=0.0,E_tot=0.0;
	atm_typ atomJ[25], atomK[25],H,atm;
	a_type  AB=AminoAcidAlphabetPDB(P),nAB=NucleotideAlphabetPDB(P);
	Int4	a,d,n,h,Nj,Nk,j,k;
	char	aa,aaJ,aaK;

        if(ResJ == ResK) return FALSE;
        aa = GetResidueAtom(AtomResidue(1,ResJ),AB,P->nA,&na); aa=AlphaChar(aa,AB);
        if(na || !(aa == 'F' || aa == 'Y' || aa == 'W')) return FALSE;
	aaJ=aa;

        aa = GetResidueAtom(AtomResidue(1,ResK),AB,P->nA,&na); aa=AlphaChar(aa,AB);
        if(na || !(aa == 'F' || aa == 'Y' || aa == 'W')) return FALSE;
	aaK=aa;

	for(Nj=0,j=1; j <= ResidueAtomNumber(ResJ); j++){
		atm=AtomResidue(j,ResJ);
		if(!SideAtom(atm)) continue;
// atm_typ AverageAtom(Int4 natom, atm_typ *atm);
		Nj++; atomJ[Nj]=atm; 
		if(!ResidueAtomHydrogens(j,ResJ)) continue;
		for(Int4 hj=1; (H=ResidueAtomHydrogen(j,hj,ResJ)); hj++){
			Nj++; atomJ[Nj]=H; 
		}
	}
	for(Nk=0,k=1; k <= ResidueAtomNumber(ResK); k++){
		atm=AtomResidue(k,ResK);
		if(!SideAtom(atm)) continue;
		Nk++; atomK[Nk]=atm;
		if(!ResidueAtomHydrogens(k,ResK)) continue;
		for(Int4 hk=1; (H=ResidueAtomHydrogen(k,hk,ResK)); hk++){
			Nk++; atomK[Nk]=H; 
		}
	}
	// Calculate nonbonded potential energy

	atm_typ	centroid_J = AverageAtom(Nj,atomJ);
	atm_typ	centroid_K = AverageAtom(Nk,atomK);
	float dist=DistanceAtoms(centroid_J,centroid_K);
	NilAtom(centroid_J); NilAtom(centroid_K);
	if(dist <= 7.0){
	  for(j=1; j < Nj; j++){
	    for(k=j+1; k < Nk; k++){
		if(GetKitaigorodsky(&E,atomJ[j],atomK[k],AB,nAB)) E_tot+=E;
	    }
	  }
	} else return FALSE;
   if(E_tot <= -0.60){	// Burley indicates E <= -0.6 as significant...
	if(chainJ && chainK){
	   fprintf(fp,"%d%c%d%c.{%c},%c%d%c.{%c}  ",file_id,
		aaJ,ResidueID(ResJ),chainJ,colorJ,aaK,ResidueID(ResK),chainK,colorK);
	} else if(chainJ){
	   fprintf(fp,"%d%c%d%c.{%c},%c%d.{%c}  ",file_id,
		aaJ,ResidueID(ResJ),chainJ,colorJ,aaK,ResidueID(ResK),colorK);
	} else if(chainK){
	   fprintf(fp,"%d%c%d.{%c},%c%d%c.{%c}  ",file_id,
		aaJ,ResidueID(ResJ),colorJ,aaK,ResidueID(ResK),chainK,colorK);
	} else fprintf(fp,"%d%c%d.{%c},%c%d.{%c}  ",file_id,
		aaJ,ResidueID(ResJ),colorJ,aaK,ResidueID(ResK),colorK);
	fprintf(fp,"// E = %.2f d_centroids = %.2f \n",E_tot,dist);
   }
#if 0
	fprintf(stderr,"// E = %.2f (%c%d vs %c%d) = %.2f Angstroms\n",E_tot,
		aaJ,ResidueID(ResJ),aaK,ResidueID(ResK),dist);
#endif
	*Energy=E_tot;
	return TRUE;
}

void	PutAroAroInteractions(FILE *fptr, Int4 file_id, Int4 resStart, Int4 resEnd,
                float HA_dmax, float dmax,Int4 C, Int4 C2, char chain,
                BooLean *skip, char color, pdb_typ P)
{
	Int4 res0;
        Int4 r,i,j,k,num_resA,num_resD;
        res_typ ResK,ResJ,**ResALL;
        Int4    num_resC,num_resO,*num_resALL;

        assert(skip);
        if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
        if(C==0){ fprintf(stderr,"chain = %c\n",chain); print_error("chain not found!"); }
        if(chain==0) chain='?';
        NEWP(ResALL,nChainsPDB(P)+3,res_typ);
        NEW(num_resALL,nChainsPDB(P)+3,Int4);
        for(C2=1; C2 <= nChainsPDB(P); C2++){
          ResALL[C2] = MakeResPDB(C2,&num_resALL[C2],P);
        }
        num_resC=num_resALL[C]; 
	// h_type HG=Histogram("gblast -log(p)",-100,100,1.0);
	for(j=1; j < num_resC; j++){
	  ResJ=ResALL[C][j];
	  r=ResidueID(ResJ);
	  if(r < resStart || r > resEnd || skip[r]) continue;
	  for(k=j+1; k <= num_resC; k++){
	  	ResK=ResALL[C][k];
	        r=ResidueID(ResK);
	        if(r < resStart || r > resEnd || skip[r]) continue;
		double E=0.0;
		if(BurleyAromatic(fptr,file_id,&E,ResJ,ResK,chain,chain,color,color,P)){
		 if(E < 0){ 
		   // PrintResidueAtoms(fptr,ResJ);
		   // PrintResidueAtoms(fptr,ResK);
		   // fprintf(stderr,"E = %g\n\n",E);
		 } // else fprintf(fptr,"failed: E = %g\n\n",E);
		 // IncdHist(E,HG);
		} 
	  }
	}
	// PutHist(stderr,60,HG); NilHist(HG);
}





