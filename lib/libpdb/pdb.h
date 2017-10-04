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

/****************** pdb.h ***************/ 
#if !defined(_PDB_)
#define _PDB_
#include "afnio.h"
#include "stdinc.h"
#include "alphabet.h"
#include "residues.h"
#include "histogram.h"
#include "probability.h"
#include "random.h"
#include "dheap.h"
#include "mheap.h"
#include "atom.h"
#include "res_typ.h"
#include "geometry.h"
#include "sip_typ.h"
#include "sequence.h"

/************************** pdb datatype *****************************/
#define MAX_NGRPS_PDB		200
#define MAX_RESIDUES_PDB	30000

#if 0
typedef struct bnd_type {
        atm_typ Donor;
        atm_typ Accept;
        atm_typ Hydrogen;
        char    type;   // aromatic or classical
} ;

typedef struct hbd_type {
        atm_typ Donor;
        atm_typ Accept;
        atm_typ Hydrogen;
	atm_typ	Virtual;	// aromatic only.
	atm_typ	Projected;	// aromatic only.
	double	d_DA;		// or d_DV;
	double	d_HA;		// or d_HV;
	double	angleDHA;	// or angleDHV;
	double	d_HpX;		// projected on aromatic plane
        char    type;   // aromatic or classical
} ;

// Compare difference between H-bond lists.

#endif

typedef struct {
	char		*filename;
	Int4		nchains;
	char		chain[MAX_NGRPS_PDB + 3];
	Int4		natoms[MAX_NGRPS_PDB + 3];
	atm_typ		**atom; 	/** atoms[chain][atom] **/
	Int4		maxres[MAX_NGRPS_PDB + 3];  /** #res in chain C **/
	Int4		lowres[MAX_NGRPS_PDB + 3];  /** lowest #res in chain C **/
	unsigned char	**seq;		/** sequence[chain][res] **/
	unsigned char	**seq0;		/** offset sequence pointers **/
	float		*C_total;	/** Sander contacts **/
	BooLean		*protein;	/** protein[chain] = TRUE? **/
	Int4		*temp;		/** temporary residue list **/
	a_type		A,nA;		/** alphabet **/
	Int4            nresidues;      // Number of residues
	FILE		*fp_stderr;	// stderr that can be set to fopen("/dev/null","w);
	Int4		num_res[MAX_NGRPS_PDB + 3];
	res_typ		*residues[MAX_NGRPS_PDB + 3];  /** res in chain C **/
} protdb_type;

typedef	protdb_type	*pdb_typ;

/******************************** private **********************************/
void	count_atoms(FILE *fptr, pdb_typ P);
Int4    bubble_sort_pdb(Int4 *L);
void	get_atoms(FILE *fptr, pdb_typ P);
char    get_residue_pdb(register char *name, register pdb_typ P,
        register BooLean *na);
void    pdb_error(const char *s);
void	CountResPDB(pdb_typ P);

/******************************** PUBLIC ***********************************/
pdb_typ MakePDB(FILE *fptr,char *infile);
pdb_typ MakePDB(FILE *fptr);
pdb_typ MakePDB(char *infile);
pdb_typ	MkPDB(char *infile);

void	ErrorToDevNullPDB(pdb_typ pdb);
void    PutPDB(FILE *fptr, pdb_typ P);
Int4	PutPDBSeq(FILE *fptr, Int4 Chain, pdb_typ P);
BooLean PutChainPDB(FILE *fptr, char Chain, pdb_typ P);
BooLean PutChainPDB(FILE *fptr, Int4 StartRes, Int4 EndRes, char Chain, pdb_typ P);
BooLean PutPrtnChainPDB(FILE *fptr, char Chain, pdb_typ P);
BooLean PutNotChainPDB(FILE *fptr, char Chain, pdb_typ P);
void	PutAdjacentSubunitsPDB(FILE *fptr, Int4 C0, double maxdist,pdb_typ P);
void    MultiPutAdjacentSubunitsPDB(FILE *fptr, Int4 *c0, double maxdist,pdb_typ P);
BooLean *AdjacentSubunitsPDB(Int4 C0, double maxdist,pdb_typ P);
BooLean *AdjacentHeteroSubunitsPDB(Int4 C0, double maxdist,pdb_typ P);
char    **AdjacentHeteroMoleculePDB2(Int4 C0, double maxdist,pdb_typ P);
char    **AdjacentHeteroMoleculePDB(Int4 C0, double maxdist,pdb_typ P);
BooLean PutNotPrtnChainPDB(FILE *fptr, char Chain, pdb_typ P);

e_type  GetPDBSeq(Int4 Chain, pdb_typ P);
BooLean RenumberChainPDB(char Chain, Int4 start, pdb_typ P);
Int4    FindPotentialHBondsPDB(Int4 res0, Int4 C, Int4 C2, float dmax, pdb_typ P);
void    RmDistWaterPDB(FILE *fptr, double maxdist,pdb_typ P);
BooLean ShiftResChainPDB(char Chain, Int4 shift, pdb_typ P);
BooLean ShiftResChainPDB(char Chain, Int4 shift, Int4 start, Int4 end, pdb_typ P);
BooLean ReNameChainPDB(Int4 chain, char Chain, pdb_typ P);
void    TransformPDB(double xyz[4][5], pdb_typ P);
void	LabelNullChainsPDB(pdb_typ P);
void    LabelNullChainsPDB(pdb_typ P, BooLean RmX);
Int4    ReLabelPDB(pdb_typ P);
void	ReNamePDB(const char *str,pdb_typ P);

float   *OoiNumberPDB(Int4 C, pdb_typ P);
float   *ExposedPDB(Int4 C, pdb_typ P);
void	NeighborsPDB(FILE *fptr, Int4 C, Int4 res, pdb_typ P);
Int4    *ContactsPDB(Int4 C, Int4 N, Int4 *C2, float dmax, pdb_typ P);
Int4    ContactsResPDB(Int4 res0, Int4 C, Int4 C2, float dmax, pdb_typ P);
atm_typ ContactsAtomPDB(Int4 atom_no, float dmax, pdb_typ P);
Int4    CntExposedPDB(Int4 *surface, Int4 *buried, pdb_typ P);
double  DistancePDB(Int4 res1, char *name1, Int4 C1, Int4 res2,   
        char *name2, Int4 C2, pdb_typ P);
Int4    TriadPDB(Int4 C, pdb_typ P,Int4 *triad[5]);
void    TriadSuperImposePDB(FILE* fptr, pdb_typ ToP, pdb_typ FromP, pdb_typ P);
void    TriadTwistChnPDB(FILE* fptr, pdb_typ ToP, pdb_typ FromP, Int4 chn,
        Int4 twistpoint,pdb_typ P);

long    CountHydrogensPDB(pdb_typ P);
Int4	ChangeMSE2MET_PDB(pdb_typ P);
Int4    ChangePO3toAtmPDB(pdb_typ P);

atm_typ	AtomPDB(Int4 C, Int4 atom_num, pdb_typ P);
void	NilPDB(pdb_typ P);
BooLean IsFullProteinChainPDB(Int4 c, pdb_typ P);

//======================== pdb_misc.cc ========================
Int4	FindQuadWatersPDB(float dmax, pdb_typ P);
void    PutContactsPDB(FILE *fptr, Int4 file_id, Int4 resStart,Int4 resEnd, float HA_dmax,
                float dmax,Int4 C, Int4 C2, char chain, BooLean *skip, char color, pdb_typ P);
Int4    HBondAnglesPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
                unsigned short file_id, pdb_typ P);
void    PutAngularPDB(FILE *fptr, Int4 file_id, Int4 resStart,Int4 resEnd, float HA_dmax,
                float dmax,Int4 C, Int4 C2, char chain, BooLean *skip, char color, pdb_typ P);
double  PrintThetaDist(double &theta,FILE *ofp,atm_typ Xd,atm_typ D,atm_typ A,atm_typ Xa,
                double d_HA,double angle_DHA);
Int4    GetHBondThetaAndDistPDB(FILE *fp,FILE *ofp,res_typ Res,Int4 C,Int4 C2,float dmax,res_typ Res2,
                unsigned short file_id, pdb_typ P, char mode, char *donor, char *accept,
                h_type tHG,h_type dHG,h_type dhaHG,h_type haHG);

/******************************** MACROS ***********************************/
#define SREL2_DNA "AGTC"
        /********* N   A   C   G   T   U *******/
#define DNA_MTRX "-4  -4  -4  -4  -4  -4\
                  -4   5  -4  -4  -4  -4\
                  -4  -4   5  -4  -4  -4\
                  -4  -4  -4   5  -4  -4\
                  -4  -4  -4  -4   5   5\
                  -4  -4  -4  -4   5   5"

//*********************** AFN: H-bonds *******************************
Int4    FindHBondsPDB(FILE *fp,res_typ Res, Int4 C, Int4 C2, float dmax, pdb_typ P);
Int4    FindHBondsPDB(res_typ Res, Int4 C, Int4 C2, float dmax, pdb_typ P);
Int4    FindHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
		char color, pdb_typ P);
Int4    FindHBondsPDB(res_typ Res, Int4 C, Int4 C2, float dmax, Int4 Res2, pdb_typ P);
Int4    FindHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
                char color, unsigned short file_id, pdb_typ P);

Int4    FindWaterHBondsPDB(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
                char color, pdb_typ P,Int4 water);

Int4    FindAromaticHbondsPDB(FILE *fp,res_typ ResD,res_typ ResA,Int4 C,Int4 C2,
        float dmax, char color, pdb_typ P,double &minDist);
Int4    FindAromaticHbondsPDB(FILE *fp,res_typ ResD,res_typ ResA,Int4 C,Int4 C2,
        float dmax, char color, pdb_typ P);
Int4    FindAromaticHbondsPDB(res_typ ResD, res_typ ResA,Int4 C, Int4 C2,
        float dmax, pdb_typ P);

Int4    FindOtherPiHbondsPDB(FILE *fp,res_typ ResD,res_typ ResA,Int4 C,Int4 C2,
        float dmax, char color, pdb_typ P,double &minDist);
Int4    FindOtherPiHbondsPDB(FILE *fp,res_typ ResD,res_typ ResA,Int4 C,Int4 C2,
        float dmax, char color, pdb_typ P);
Int4    FindOtherPiHbondsPDB(res_typ ResD, res_typ ResA,Int4 C, Int4 C2,
        float dmax, pdb_typ P);

//*********************** From Hata code ********************************
/********************************* private *********************************/
// inline void  SetResidueNumber(pdb_typ P, knt N){P->nresidues = N;};
/******************************** public ***********************************/
inline atm_typ GetAtomPDB(pdb_typ P, Int4 Nc, Int4 Na){ return P->atom[Nc][Na];};
Int4     GetChainNumberPDB(pdb_typ P, char Chain);
Int4     GetResMinAtomNumber(pdb_typ P, Int4 Chain, Int4 Nres);
Int4     GetResNumberOfAtoms(pdb_typ P, Int4 Chain, Int4 Nres);
res_typ *MakeResPDB(Int4 chain,Int4 *num_res,pdb_typ P);
/***************************************************************************/
#define AminoAcidAlphabetPDB(P) ((P)->A)
#define NucleotideAlphabetPDB(P) ((P)->nA)
#define NumAtomsPDB(P)            ((P)->natoms)
#define NumberAtomsPDB(c,P)       (((c) > 0 && (c) <= (P)->nchains)?((P)->natoms[(c)]):0)
#define AtomsPDB(P)            	  ((P)->atom)
#define FilenamePDB(P)   	((P)->filename)
#define MaxAtomsPDB(c,P)          (((c) > 0 && (c) <= (P)->nchains)?((P)->natoms[(c)]):0)
#define LowestResPDB(c,P)         (((c) > 0 && (c) <= (P)->nchains)?((P)->lowres[(c)]):0)
#define MaxResPDB(c,P)         (((c) > 0 && (c) <= (P)->nchains)?((P)->maxres[(c)]):0)
#define MinResPDB(c,P)         (((c) > 0 && (c) <= (P)->nchains)?((P)->lowres[(c)]):0)
#define nChainsPDB(P)		((P)->nchains)
#define IsProteinPDB(c,P)         (((c) > 0 && (c) <= (P)->nchains)?((P)->protein[(c)]):0)
#define ChainCharPDB(c,P)         (((c) > 0 && (c) <= (P)->nchains)?((P)->chain[(c)]):0)


#endif




