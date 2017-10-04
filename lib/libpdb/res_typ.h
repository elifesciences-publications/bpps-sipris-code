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

#if !defined(RES_TYP)
#define RES_TYP

#include "atom.h"

//typedef struct residue_type{
typedef struct {
  char		name[4];	// residue name 
  Int4          id_num;         // ID number  
  Int4		chain;		// chain name
  atm_typ       *atom;          // array of atoms in residue
  atm_typ	aromatic[2];	// X center atoms for aromatic compounds.
  atm_typ       **hydrogen;     // array of hydrogens bound to each atom in residue
  Int4          natom;          // number of atoms in residue
  short		num_hyd;        // number of hydrogens in residue
  short		num_heavy;      // number of non-hydrogens in residue
  char          structure;      // secondary sturcture
  BooLean       hbond;          // has hydrogen bond
  double	sfc_area;       // surface accessibility
  double	sfc_area_exclusive;  // surface accessibility w/o other res
  double	sfc_area_ratio;    // sfc_area / sfc_area_exclusive

  atm_typ       **heavy;	// array of heavy atoms bound to each heavy atom in residue
  char		*num_heavy_bonds;// number of covalent bonds to other heavy atoms.
} residue_type;

typedef residue_type *res_typ; 

/******************************** Private **********************************/
double  ResidueToResidueContact(res_typ R, res_typ R2,double *ratio,
		BooLean SideChainOnly);
double	CalculateResidueSurface(res_typ R,Int4 NumChains,Int4 *NumInChain, 
					atm_typ **AllAtoms);
void PrintResidueSurface(FILE *fp,res_typ R);
void PrintResidueAtoms(FILE *fp,res_typ Res,BooLean putHydrogens=TRUE);
inline void SetResidueSurfaceArea(res_typ R, float S){
  R->sfc_area_exclusive = S;  R->sfc_area_ratio = R->sfc_area / S;};
inline void SetResidueHB(res_typ R, BooLean hb){ R->hbond = hb;};
inline void SetSecondaryStructure(res_typ R, char S){R->structure = S;};
inline atm_typ GetAtomFromResidue(res_typ r, int Na){return r->atom[Na];};
atm_typ  GetAtomFromResidueByName(res_typ r, char *atomname);
char	GetCharResidue(res_typ R,a_type AB);
/******************************** Public ***********************************/
atm_typ MkAveAtomResidue(res_typ R);
void	PutBondLengthHGs(FILE *fp);
void	PutBondsInRes(FILE *fp,res_typ Res);
void    PutBondHeavyInRes(FILE *fp,Int4 a, res_typ Res);
res_typ MkRes(Int4 num_atm, Int4 chain, atm_typ *atm);
int     GetHBAtomsInResidue(res_typ r, atm_typ *hbatoms);
void	NilRes(res_typ res);
BooLean CarbonBondRangeAtom(double &min, double &max, atm_typ atm);
/******************************** MACROS ***********************************/
#define ResidueName(R)              ((R)->name)
#define AtomResidue(a,R)            (((a) > 0 && (a) <= (R)->natom) ? (R)->atom[(a)]: 0)
#define ChainResidue(R)             AtomChain((R)->atom[1])
#define AtomsInResidue(R)           ((R)->natom)
#define ResidueID(R)                ((R)->id_num)
#define ResidueChain(R)             ((R)->chain)
#define ResidueAtomNumber(R)        ((R)->num_heavy)
#define ResidueAtomHydrogen(a,h,R)  ((R)->hydrogen[(a)][(h)])
#define ResidueAtomHydrogens(a,R)   ((R)->hydrogen[(a)])
#define ResidueBoundHeavyAtom(a,b,R)((R)->heavy[(a)][(b)])
#define ResidueNumBoundHeavyAtoms(a,R)  ((R)->num_heavy_bonds[(a)])
#define ResidueSurfaceArea(R)       ((R)->sfc_area)
#define ResidueSurfaceExclusive(R)  ((R)->sfc_area_exclusive)
#define ResidueSurfaceAccessibility(R)       ((R)->sfc_area_ratio)
/***************************************************************************/

#endif

