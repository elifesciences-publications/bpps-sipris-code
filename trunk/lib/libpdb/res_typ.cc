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

#include "res_typ.h"

#if 0

H  1.008         0.161               H bonded to nitrogen atoms
HC 1.008         0.135               H aliph. bond. to C without electrwd.group
H1 1.008         0.135               H aliph. bond. to C with 1 electrwd. group
H2 1.008         0.135               H aliph. bond. to C with 2 electrwd.groups
H3 1.008         0.135               H aliph. bond. to C with 3 eletrwd.groups
HA 1.008         0.167               H arom. bond. to C without elctrwd. groups
H4 1.008         0.167               H arom. bond. to C with 1 electrwd. group
H5 1.008         0.167               H arom.at C with 2 elctrwd. gr,+HCOO group
HO 1.008         0.135               hydroxyl group
HS 1.008         0.135               hydrogen bonded to sulphur (pol?)
HW 1.008         0.135               H in TIP3P water
HP 1.008         0.135               H bonded to C next to positively charged gr
HZ 1.008         0.161               H bond sp C (Howard et al.JCC,16,243,1995)
SH 32.06         2.900               S in cystine

bond		distance
OW-HW  553.0    0.9572    ! TIP3P water
HW-HW  553.0    1.5136      TIP3P waterC -OH  450.0    1.364       JCC,7,(1986),230; (not used any more for TYR) 
C -H4  367.0    1.080       Junmei et al, 1999
C -H5  367.0    1.080       Junmei et al, 1999
CA-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; PHE,TRP,TYR
CA-H4  367.0    1.080       changed from 340. bsd on C6H6 nmodes; no assigned
CD-HA  367.0    1.080       Junmei et al, 1999 
CK-H5  367.0    1.080       changed from 340. bsd on C6H6 nmodes; ADE,GUA
CM-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; CYT,URA
CM-H4  367.0    1.080       changed from 340. bsd on C6H6 nmodes; CYT,URA
CM-H5  367.0    1.080       changed from 340. bsd on C6H6 nmodes; not assigned
CQ-H5  367.0    1.080       changed from 340. bsd on C6H6 nmodes; ADE
CR-H5  367.0    1.080       changed from 340. bsd on C6H6 nmodes;HIS
CV-H4  367.0    1.080       changed from 340. bsd on C6H6 nmodes; HIS
CW-H4  367.0    1.080       changed from 340. bsd on C6H6 nmodes;HIS(epsilon,+)

CT-HC  340.0    1.090       changed from 331 bsd on NMA nmodes; AA, SUGARS
CT-H1  340.0    1.090       changed from 331 bsd on NMA nmodes; AA, RIBOSE
CT-H2  340.0    1.090       changed from 331 bsd on NMA nmodes; SUGARS
CT-H3  340.0    1.090       changed from 331 bsd on NMA nmodes; not assigned
CT-HP  340.0    1.090       changed from 331; AA-lysine, methyl ammonium cation

CZ-HZ  400.0    1.056       Howard et al JCC,16,243,1995

HO-OH  553.0    0.960       JCC,7,(1986),230; SUGARS,SER,TYR
HO-OS  553.0    0.960       JCC,7,(1986),230; NUCLEOTIDE ENDS

H -N2  434.0    1.010       JCC,7,(1986),230; ADE,CYT,GUA,ARG
H -N*  434.0    1.010       for plain unmethylated bases ADE,CYT,GUA,ARG
H -NA  434.0    1.010       JCC,7,(1986),230; GUA,URA,HIS
H -N   434.0    1.010       JCC,7,(1986),230; AA
H -N3  434.0    1.010       JCC,7,(1986),230; LYS    
H -NT  434.0    1.010       for neutral amines 
#endif

static h_type HGO=0;
static h_type HGN=0;
static h_type HGS=0;
static h_type HGC=0;
static h_type HGX=0;

void    PutBondLengthHGs(FILE *fp)
{
	if(HGO){ PutHist(fp,60,HGO); NilHist(HGO); HGO=0; }
	if(HGN){ PutHist(fp,60,HGN); NilHist(HGN); HGN=0; }
	if(HGS){ PutHist(fp,60,HGS); NilHist(HGS); HGS=0; }
	if(HGC){ PutHist(fp,60,HGC); NilHist(HGC); HGC=0; }
	if(HGX){ PutHist(fp,60,HGX); NilHist(HGX); HGX=0; }
}

atm_typ MkAveAtomResidue(res_typ R)
{ if(R == 0) return 0; else return AverageAtom(R->num_heavy,R->atom); }

BooLean	CarbonBondRangeAtom(double &min, double &max, atm_typ atm1, atm_typ atm2)
// return range of covalent bond distances in Angstroms from a carbon atom to atm.
{
	atm_typ atm;
	if(CarbonAtom(atm1)) atm=atm2; else if(CarbonAtom(atm2)) atm=atm1; else return FALSE;
#if 0	// modify this to allow some slack...
	if(OxygenAtom(atm)){ max=2.15; min=0.95; }
	else if(NitrogenAtom(atm)){ max=2.10; min=1.05; }
	else if(CarbonAtom(atm)){ max=1.54; min=1.05; }
	else if(SulfurAtom(atm)){ max=2.55; min=1.25; }
#else	// the following are a bit more liberal...
	if(OxygenAtom(atm)){ max=2.25; min=0.95; }
	else if(NitrogenAtom(atm)){ max=2.15; min=1.05; }
	else if(CarbonAtom(atm)){ max=1.64; min=1.05; }
	else if(SulfurAtom(atm)){ max=2.60; min=1.25; }
#endif
	else return FALSE; return TRUE;
}

void    PutBondHeavyInRes(FILE *fp,Int4 a, res_typ Res)
{
	Int4 n=0;
	assert(a > 0 && a <= Res->num_heavy);
	atm_typ atm=Res->atom[a]; PutAtom(fp, atm); 
	if(Res->num_heavy_bonds[a] > 0){	// for each heavy atom bound to a...
	     for(Int4 j=1; j <= Res->num_heavy_bonds[a]; j++){
	       fprintf(fp," "); PutAtom(fp,Res->heavy[a][j]); n++;
	     }
	}
	if(Res->hydrogen[a]){		// for each hydrogen bound to a...
	       for(Int4 h=1; Res->hydrogen[a][h]; h++){
		  fprintf(fp," "); PutAtom(fp, Res->hydrogen[a][h]); n++;
	       }
	} if(n > 0) fprintf(fp,"\n");
}

void    PutBondsInRes(FILE *fp,res_typ Res)
{
	Int4 n=0;
	for(Int4 a=1; a <= Res->num_heavy; a++){
	   if(Res->num_heavy_bonds[a] > 0){	// for each heavy atom bound to a...
	     atm_typ atm=Res->atom[a]; PutAtom(fp, atm); 
	     for(Int4 j=1; j <= Res->num_heavy_bonds[a]; j++){
	       fprintf(fp," "); PutAtom(fp,Res->heavy[a][j]); n++;
	     }
	   }
	   if(Res->hydrogen[a]){		// for each hydrogen bound to a...
	       for(Int4 h=1; Res->hydrogen[a][h]; h++){
		  fprintf(fp," "); PutAtom(fp, Res->hydrogen[a][h]); n++;
	       }
	   } if(n > 0) fprintf(fp,"\n");
	}
}

res_typ	MkRes(Int4 num_atm, Int4 chain, atm_typ *atm)
// AFN: create residue type
// atm is array of 
{
	atm_typ	*hyd,*heavy;
	Int4	i,j,k,h,a;
#if 1	// DEBUG
	if(HGO==0){
		HGO=Histogram("O-H bond lengths",800,1500,0.05);
		HGN=Histogram("N-H bond lengths",800,1500,0.05);
		HGS=Histogram("S-H bond lengths",800,1500,0.05);
		HGC=Histogram("C-H bond lengths",800,1500,0.05);
		HGX=Histogram("X-H bond lengths",800,1500,0.05);
	}
#endif
	res_typ	res = new residue_type;
	NEW(res->atom,num_atm+2,atm_typ);
	NEW(hyd,num_atm+2,atm_typ);
	strcpy(res->name, AtomResName(atm[1]));
	res->sfc_area = 0.0;
	for(h=0,a=0,i=1; i <= num_atm; i++){ 
		if(!HydrogenAtom(atm[i])){ 
		   a++; res->atom[a] = atm[i];
		} else {
		   h++; hyd[h]=atm[i];
		}
		// res->sfc_area += GetAtomSurface(atom);
	}
	res->natom = num_atm;
	res->num_heavy = a;
	res->num_hyd = h;
	NEWP(res->hydrogen,a+2,atm_typ);
	NEWP(res->heavy,a+2,atm_typ);
	NEW(res->num_heavy_bonds,a+2,char);
	heavy=res->atom;
	Int4	h0=0;
	for(i=1; i <= a; i++){
	  double maxbond,minbond; // minbond avoids problems with strange coordinate files.
#if 1	  // find heavy atoms bound to each heavy atom. // afn: 7/12/12.
	  for(k=0,j=1; j <= a; j++){
	  	if(!CarbonBondRangeAtom(minbond,maxbond,heavy[i],heavy[j])) continue;
		//   ^ this sets the values for minbond & maxbond!
		double d=DistanceAtoms(heavy[i],heavy[j]);
		if(d <= maxbond && d >= minbond){
		   if(!res->heavy[i]) NEW(res->heavy[i],res->num_heavy + 3,atm_typ);
		   k++; assert(k < res->num_heavy); 
		   assert(res->heavy[i][k]==0); 
		   res->heavy[i][k]=heavy[j];
		}
	  } res->num_heavy_bonds[i]=k;
#endif
	  // distances for hydrogen atoms...
#if 1	// new reduce settings...
	  if(OxygenAtom(heavy[i])){ maxbond=1.05; minbond=0.80; }
	  else if(NitrogenAtom(heavy[i])){ maxbond=1.05; minbond=0.80; }
	  else if(CarbonAtom(heavy[i])){ maxbond=1.15; minbond=0.92; }
	  else if(SulfurAtom(heavy[i])){ maxbond=1.35; minbond=1.13; }
	  else { maxbond=1.4; minbond=0.8; }
#else	  // more liberal settings...
	  if(OxygenAtom(heavy[i])){ maxbond=1.15; minbond=0.95; }
	  else if(NitrogenAtom(heavy[i])){ maxbond=1.15; minbond=0.95; }
	  else if(CarbonAtom(heavy[i])){ maxbond=1.25; minbond=1.05; }
	  else if(SulfurAtom(heavy[i])){ maxbond=1.45; minbond=1.25; }
	  else { maxbond=1.4; minbond=0.8; }
#endif
	  for(k=0,j=1; j <= h; j++){
		double d=DistanceAtoms(heavy[i],hyd[j]);
		if(d <= maxbond && d >= minbond){
#if 1	// DEBUG
		  if(OxygenAtom(heavy[i])) IncdHist(d*1000, HGO);
		  else if(NitrogenAtom(heavy[i])) IncdHist(d*1000, HGN);
		  else if(CarbonAtom(heavy[i])) IncdHist(d*1000, HGC);
		  else if(SulfurAtom(heavy[i])) IncdHist(d*1000, HGS);
		  else IncdHist(d*1000, HGX);
#endif
		  if(!res->hydrogen[i]) NEW(res->hydrogen[i],10,atm_typ);
#if 0	// DEBUG
if(maxbond==1.35){
	PutAtom(stderr,heavy[i]);
	PutAtom(stderr,hyd[j]);
	fprintf(stderr,"Distance = %.3f\n",DistanceAtoms(heavy[i],hyd[j]));
}
#endif
		   k++; assert(k < 9); 
		   assert(res->hydrogen[i][k]==0); 
		   res->hydrogen[i][k]=hyd[j];
		   h0++; // found a bond to hydrogen
		}
	  }
	}
	res->chain = chain;
	res->id_num = ResAtom(atm[1]);
	if(h0 != h){
		if(h0 < h) fprintf(stderr,"WARNING: %d hydrogens unassigned.\n",h-h0);
		// else fprintf(stderr,"WARNING: multiply assigned %d hydrogens.\n",h0-h);
		// PrintResidueAtoms(stderr,res);
	}
	free(hyd);
	return res;
}

void    NilRes(res_typ res)
{
	for(Int4 i=1; i <= res->num_heavy; i++){
		if(res->hydrogen[i]) free(res->hydrogen[i]);
		if(res->heavy[i]) free(res->heavy[i]);
	} free(res->hydrogen); free(res->heavy); free(res->num_heavy_bonds);
	free(res->atom); delete res; 
}

//************************ AFN: added the following routine: **************
double	ResidueToResidueContact(res_typ R, res_typ R2,double *ratio, 
			BooLean SideChainOnly)
{
  Int4		Na, N, Ni,C;
  double	S = 0.0,E=0.0;
  atm_typ	A;
  float		A_area, Ratio;
  atm_typ 	A2,neighbors[500];
  Int4		Natom = ResidueAtomNumber(R);

  // Calculate residue by itself...
  for(S=0.0, Na = 1; Na <= Natom; Na++){
    A = R->atom[Na];
    if(HydrogenAtom(A)) continue;
    if(SideChainOnly && !SideAtom(A)) continue;
    for(Ni=0,N=1; N <= Natom; N++){
      A2=R->atom[N];
      if(HydrogenAtom(A2)) continue;
      if(IsWaterAtom(A2)) continue;
      if(SideChainOnly && !SideAtom(A2)) continue;
      if(N != Na){	// if not this atom.
	float d  =  DistanceAtoms(A2, A);
	float rr =  EffectiveRadius(A) + EffectiveRadius(A2);
	if(d <= rr) { neighbors[Ni] = A2; Ni++; }
      }
    } 
    CalculateAtomSurface(A, neighbors, Ni, &A_area, &Ratio);
    // PutAtom(stderr,A);
    // fprintf(stderr,"LEU Surface = %.3f; %.3f\n",A_area,Ratio);
    if(SideAtom(A)){ S += A_area; }
    else if(!SideChainOnly) S += A_area;
  }
  R->sfc_area_exclusive = S;

  // Calculate test residue aInt4 with the residue of interest.
  S=E=0.0;
  for(Na = 1; Na <= Natom; Na++){
    A = R->atom[Na]; 
    if(HydrogenAtom(A)) continue;
    if(SideChainOnly && !SideAtom(A)) continue;
    // First find neighbor atoms from query residue.
    for(Ni=0,N=1; N <= Natom; N++){
      A2=R->atom[N];
      if(HydrogenAtom(A2)) continue;
      if(IsWaterAtom(A2)) continue;
      if(SideChainOnly && !SideAtom(A2)) continue;
      if(A != A2){	// if not this atom.
	float d  =  DistanceAtoms(A2, A);
	float rr =  EffectiveRadius(A) + EffectiveRadius(A2);
	if(d <= rr) { neighbors[Ni] = A2; Ni++; }
      }
    } 
    // Next find neighbor atoms from other residue.
    Int4 Natom2 = ResidueAtomNumber(R2);
    for(N = 1; N <= Natom2; N++){
      A2=R2->atom[N];
      if(A2 == NULL) continue;
      if(HydrogenAtom(A2)) continue;
      if(SideChainOnly && !SideAtom(A2)) continue;
      if(IsWaterAtom(A2)) continue;
      if(A != A2){ // if not this atom.
	float d  =  DistanceAtoms(A2,A);
	float rr =  EffectiveRadius(A) + EffectiveRadius(A2);
	if(d <= rr) { neighbors[Ni] = A2; Ni++; 
		// fprintf(stderr,"LEU Neighbor: ");
		// PutAtom(stderr,A2);
	}
      }
    }
    CalculateAtomSurface(A, neighbors, Ni, &A_area, &Ratio);
    // PutAtom(stderr,A);
    // fprintf(stderr,"LEU Surface = %.3f; %.3f\n",A_area,Ratio);
    if(SideAtom(A)){
        S += A_area;
        // E += A_area * Ratio; // AFN modification...
    } else if(!SideChainOnly) S += A_area;
  }
  R->sfc_area = S;	// AFN modification...
  R->sfc_area_ratio = R->sfc_area/R->sfc_area_exclusive;
  // SetResidueSurfaceArea(R, S);
  double contact_area = R->sfc_area_exclusive - R->sfc_area;
  *ratio = contact_area/R->sfc_area_exclusive;
  return contact_area;
}

//*************************************************************************
/*------------------------------------------------------------------------*/

double	CalculateResidueSurface(res_typ R,Int4 NumChains,Int4 *NumInChain, 
		                                        atm_typ **AllAtoms)
{
  Int4		Na, N, Ni,C;
  double	S = 0.0,E=0.0;
  atm_typ	A;
  float		A_area, Ratio;

  Int4		NeighborCount = 0;
  Int4		Natom = ResidueAtomNumber(R);
  atm_typ neighbors[500];


  // Calculate residue by itself...
  for (Na = 1; Na <= Natom; Na++){
    A = R->atom[Na];
    if(HydrogenAtom(A)) continue;
    for(Ni=0,N=1; N <= Natom; N++){
      if(HydrogenAtom(R->atom[N])) continue;
      if(IsWaterAtom(R->atom[N])) continue;
      if(N != Na){	// if not this atom.
	float d  =  DistanceAtoms(R->atom[N], A);
	float rr =  EffectiveRadius(A) + EffectiveRadius(R->atom[N]);
	if(d <= rr) { neighbors[Ni] = R->atom[N]; Ni++; }
      }
    } 
    CalculateAtomSurface(A, neighbors, Ni, &A_area, &Ratio);
    // PutAtom(stderr,A);
    // fprintf(stderr,"LEU Surface = %.3f; %.3f\n",A_area,Ratio);
    if(SideAtom(A)){ S += A_area; }
    else S += A_area;
  }
  R->sfc_area_exclusive = S;
  // Calculate residue aInt4 with all other atoms
  S=E=0.0;
  for (Na = 1; Na <= Natom; Na++){
    A = R->atom[Na]; 
    if(HydrogenAtom(A)) continue;
    for(Ni=0,C = 1; C <= NumChains; C++){
     if(NumInChain[C] == 0) continue;
     for(N = 1; N <= NumInChain[C]; N++){
      atm_typ A2=AllAtoms[C][N];
      if(HydrogenAtom(A2)) continue;
      if(IsWaterAtom(A2)) continue;
      if(A2 == NULL) continue;
      if(A != A2){ // if not this atom.
	float d  =  DistanceAtoms(A2,A);
	float rr =  EffectiveRadius(A) + EffectiveRadius(A2);
	if(d <= rr) { neighbors[Ni] = A2; Ni++; 
		// fprintf(stderr,"LEU Neighbor: ");
		// PutAtom(stderr,A2);
	}
      }
     }
    } NeighborCount = Ni;
    CalculateAtomSurface(A, neighbors, NeighborCount, &A_area, &Ratio);
    // PutAtom(stderr,A);
    // fprintf(stderr,"LEU Surface = %.3f; %.3f\n",A_area,Ratio);
    if(SideAtom(A)){
      S += A_area;
      // E += A_area * Ratio; // AFN modification...
    }
    else S += A_area;
  }
  R->sfc_area = S;	// AFN modification...
  R->sfc_area_ratio = R->sfc_area/R->sfc_area_exclusive;
  // SetResidueSurfaceArea(R, S);
  return R->sfc_area;
}

char    GetCharResidue(res_typ R,a_type AB)
{
	BooLean nucleic_acid;
	atm_typ	Atom=AtomResidue(1,R);
	char r=GetResidueAtom(Atom,AB,AB,&nucleic_acid);
	if(nucleic_acid) return 'X';
	return AlphaChar(r,AB);
}

/*------------------------------------------------------------------------*/

void PrintResidueSurface(FILE *fp, res_typ R)
{  
    fprintf(fp,"%s  ", R->name);
    fprintf(fp,"%4d  ", R->id_num);
    fprintf(fp,"%8.2f", R->sfc_area);
    fprintf(fp,"%8.2f", R->sfc_area_exclusive-R->sfc_area);
    fprintf(fp,"%8.2f", R->sfc_area_exclusive);
    fprintf(fp,"%8.3f\n", R->sfc_area_ratio);
}

/*------------------------------------------------------------------------*/

void PrintResidueAtoms(FILE *fp,res_typ Res, BooLean putHydrogens) 
{ 
	for(Int4 a=1; a <= Res->num_heavy; a++){
	   PutAtom(fp, Res->atom[a]); 
	   if(putHydrogens && Res->hydrogen[a]){
	     for(Int4 h=1; Res->hydrogen[a][h]; h++){
		PutAtom(fp, Res->hydrogen[a][h]);
	     }
	   }
	}
}

/*------------------------------------------------------------------------*/

int GetHBAtomsInResidue(res_typ r, atm_typ *hbatoms)
{
  int i, count = 0;
  atm_typ atom, atomhb;
  res_typ reshb;
  int res_num = 0;

  for (i = 1; i <= ResidueAtomNumber(r); i++){
    atom = GetAtomFromResidue(r, i);
    if (IsAtomHB(atom)) {
      hbatoms[count] = atom;
      count++;
    }
  }
  return count;
}

/*------------------------------------------------------------------------*/
#if 0
atm_typ  GetAtomFromResidueByName(res_typ r, char *atomname)
{
  atm_typ A=NULL;
  for (Int4 n=1; n<= ResidueAtomNumber(r); n++){
    A = GetAtomFromResidue(r, n);
    if(!strcmp(AtomRealName(A),atomname)) return A;
  } return NULL;
}
#endif
		
/*------------------------------------------------------------------------*/
