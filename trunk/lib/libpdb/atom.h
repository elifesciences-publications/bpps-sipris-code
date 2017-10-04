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

/****************** atom.h ***************/
#if !defined(ATOM)
#define ATOM
#include "afnio.h"
#include <math.h>
#include "alphabet.h"
#include "residues.h"
#include "histogram.h"
#include "probability.h"
#include "random.h"

/************************** atom datatype *****************************/
/**********************************************************************
TER		
HELIX    1     ARG      8  LYS     14 
SHEET    1       GLY    23  GLN    29
CONECT  144  667
 **********************************************************************
ATOM    241  N3    T B   1      -8.580  -3.604  49.586  1.00 14.42     (DNA)
ATOM    449  CA  ARG C   3     -13.582   8.064  60.257  1.00 39.38     (protein)
ATOM     62  NH1 ARG     8      -7.802  -5.691  29.031  1.00 27.67      1AAK
HETATM 1205  O   HOH   151      -5.065  11.390  30.224  1.00  3.66      1AAK

<----><---> <-->^<-> ^<-->^   <---X--><---Y--><---Z--><----><----> <->
....|....|....|....|....|....|....|....|....|....|....|....|....|....|....|
    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75
/**********************************************************************/
/**********************************************************************
PDB Record Format:
COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier, left-justified.
77 - 78        LString(2)      element       Element symbol, right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
 **********************************************************************/

#if 0
#define MAXELEMNO  104
typedef struct {
           char symbol[2];
           int covalrad;
           int vdwrad;
           int cpkcol;
           char *name;
} ElemStruct;

typedef	ElemStruct	*ele_typ;

ElemStruct Element[MAXELEMNO] =  {
    { { ' ', ' ' }, 180, 360, 12, ""             },  /*   0 */
    { { 'H', ' ' },  80, 275,  4, "HYDROGEN"     },  /*   1 */
    { { 'H', 'e' }, 400, 550,  5, "HELIUM"       },  /*   2 */
    { { 'L', 'i' }, 170, 305, 14, "LITHIUM"      },  /*   3 */
    { { 'B', 'e' },  88, 157, 12, "BERYLLIUM"    },  /*   4 */
    { { 'B', ' ' }, 208, 387, 13, "BORON"        },  /*   5 */
    { { 'C', ' ' }, 180, 387,  0, "CARBON"       },  /*   6 */
    { { 'N', ' ' }, 170, 350,  1, "NITROGEN"     },  /*   7 */
    { { 'O', ' ' }, 170, 337,  2, "OXYGEN"       },  /*   8 */
    { { 'F', ' ' }, 160, 325,  6, "FLUORINE"     },  /*   9 */
    { { 'N', 'e' }, 280, 505, 12, "NEON"         },  /*  10 */
    { { 'N', 'a' }, 243, 550,  7, "SODIUM"       },  /*  11 */
    { { 'M', 'g' }, 275, 375, 15, "MAGNESIUM"    },  /*  12 */
    { { 'A', 'l' }, 338, 375,  9, "ALUMINIUM"    },  /*  13 */
    { { 'S', 'i' }, 300, 550,  6, "SILICON"      },  /*  14 */
    { { 'P', ' ' }, 259, 470,  8, "PHOSPHORUS"   },  /*  15 */  /* 262? */
    { { 'S', ' ' }, 255, 452,  3, "SULPHUR"      },  /*  16 */
    { { 'C', 'l' }, 250, 437, 13, "CHLORINE"     },  /*  17 */
    { { 'A', 'r' }, 392, 692, 12, "ARGON"        },  /*  18 */
    { { 'K', ' ' }, 332, 597, 12, "POTASSIUM"    },  /*  19 */
    { { 'C', 'a' }, 248, 487,  9, "CALCIUM"      },  /*  20 */
    { { 'S', 'c' }, 360, 330, 12, "SCANDIUM"     },  /*  21 */
    { { 'T', 'i' }, 368, 487,  9, "TITANIUM"     },  /*  22 */
    { { 'V', ' ' }, 332, 265, 12, "VANADIUM"     },  /*  23 */
    { { 'C', 'r' }, 338, 282,  9, "CHROMIUM"     },  /*  24 */
    { { 'M', 'n' }, 338, 297,  9, "MANGANESE"    },  /*  25 */
    { { 'F', 'e' }, 335, 487,  8, "IRON"         },  /*  26 */
    { { 'C', 'o' }, 332, 282, 12, "COBALT"       },  /*  27 */
    { { 'N', 'i' }, 405, 310, 10, "NICKEL"       },  /*  28 */  /* >375! */
    { { 'C', 'u' }, 380, 287, 10, "COPPER"       },  /*  29 */
    { { 'Z', 'n' }, 362, 287, 10, "ZINC"         },  /*  30 */
    { { 'G', 'a' }, 305, 387, 12, "GALLIUM"      },  /*  31 */
    { { 'G', 'e' }, 292, 999, 12, "GERMANIUM"    },  /*  32 */  /* 1225? */
    { { 'A', 's' }, 302, 207, 12, "ARSENIC"      },  /*  33 */
    { { 'S', 'e' }, 305, 225, 12, "SELENIUM"     },  /*  34 */
    { { 'B', 'r' }, 302, 437, 10, "BROMINE"      },  /*  35 */
    { { 'K', 'r' }, 400, 475, 12, "KRYPTON"      },  /*  36 */
    { { 'R', 'b' }, 368, 662, 12, "RUBIDIUM"     },  /*  37 */
    { { 'S', 'r' }, 280, 505, 12, "STRONTIUM"    },  /*  38 */
    { { 'Y', ' ' }, 445, 402, 12, "YTTRIUM"      },  /*  39 */
    { { 'Z', 'r' }, 390, 355, 12, "ZIRCONIUM"    },  /*  40 */
    { { 'N', 'b' }, 370, 332, 12, "NIOBIUM"      },  /*  41 */
    { { 'M', 'o' }, 368, 437, 12, "MOLYBDENUM"   },  /*  42 */
    { { 'T', 'c' }, 338, 450, 12, "TECHNETIUM"   },  /*  43 */
    { { 'R', 'u' }, 350, 300, 12, "RUTHENIUM"    },  /*  44 */
    { { 'R', 'h' }, 362, 305, 12, "RHODIUM"      },  /*  45 */
    { { 'P', 'd' }, 375, 360, 12, "PALLADIUM"    },  /*  46 */
    { { 'A', 'g' }, 398, 387,  9, "SILVER"       },  /*  47 */
    { { 'C', 'd' }, 422, 437, 12, "CADMIUM"      },  /*  48 */
    { { 'I', 'n' }, 408, 362, 12, "INDIUM"       },  /*  49 */
    { { 'S', 'n' }, 365, 417, 12, "TIN",         },  /*  50 */
    { { 'S', 'b' }, 365, 280, 12, "ANTIMONY"     },  /*  51 */
    { { 'T', 'e' }, 368, 315, 12, "TELLURIUM"    },  /*  52 */
    { { 'I', ' ' }, 350, 437, 11, "IODINE"       },  /*  53 */
    { { 'X', 'e' }, 425, 525, 12, "XENON"        },  /*  54 */
    { { 'C', 's' }, 418, 752, 12, "CAESIUM"      },  /*  55 */
    { { 'B', 'a' }, 335, 602,  8, "BARIUM"       },  /*  56 */
    { { 'L', 'a' }, 468, 457, 12, "LANTHANUM"    },  /*  57 */
    { { 'C', 'e' }, 458, 465, 12, "CERIUM"       },  /*  58 */
    { { 'P', 'r' }, 455, 405, 12, "PRASEODYMIUM" },  /*  59 */
    { { 'N', 'd' }, 452, 447, 12, "NEODYMIUM"    },  /*  60 */
    { { 'P', 'm' }, 450, 440, 12, "PROMETHIUM"   },  /*  61 */
    { { 'S', 'm' }, 450, 435, 12, "SAMARIUM"     },  /*  62 */
    { { 'E', 'u' }, 498, 490, 12, "EUROPIUM"     },  /*  63 */
    { { 'G', 'd' }, 448, 422, 12, "GADOLINIUM"   },  /*  64 */
    { { 'T', 'b' }, 440, 415, 12, "TERBIUM"      },  /*  65 */
    { { 'D', 'y' }, 438, 407, 12, "DYSPROSIUM"   },  /*  66 */
    { { 'H', 'o' }, 435, 402, 12, "HOLMIUM"      },  /*  67 */
    { { 'E', 'r' }, 432, 397, 12, "ERBIUM"       },  /*  68 */
    { { 'T', 'm' }, 430, 392, 12, "THULIUM"      },  /*  69 */
    { { 'Y', 'b' }, 485, 385, 12, "YTTERBIUM"    },  /*  70 */
    { { 'L', 'u' }, 430, 382, 12, "LUTETIUM"     },  /*  71 */
    { { 'H', 'f' }, 392, 350, 12, "HAFNIUM"      },  /*  72 */
    { { 'T', 'a' }, 358, 305, 12, "TANTALUM"     },  /*  73 */
    { { 'W', ' ' }, 342, 315, 12, "TUNGSTEN"     },  /*  74 */
    { { 'R', 'e' }, 338, 325, 12, "RHENIUM"      },  /*  75 */
    { { 'O', 's' }, 342, 395, 12, "OSMIUM"       },  /*  76 */
    { { 'I', 'r' }, 330, 305, 12, "IRIDIUM"      },  /*  77 */
    { { 'P', 't' }, 375, 387, 12, "PLATINUM"     },  /*  78 */
    { { 'A', 'u' }, 375, 362,  6, "GOLD"         },  /*  79 */
    { { 'H', 'g' }, 425, 495, 12, "MERCURY"      },  /*  80 */
    { { 'T', 'l' }, 388, 427, 12, "THALLIUM"     },  /*  81 */
    { { 'P', 'b' }, 385, 540, 12, "LEAD"         },  /*  82 */
    { { 'B', 'i' }, 385, 432, 12, "BISMUTH"      },  /*  83 */
    { { 'P', 'o' }, 420, 302, 12, "POLONIUM"     },  /*  84 */
    { { 'A', 't' }, 302, 280, 12, "ASTATINE"     },  /*  85 */
    { { 'R', 'n' }, 475, 575, 12, "RADON"        },  /*  86 */
    { { 'F', 'r' }, 450, 810, 12, "FRANCIUM"     },  /*  87 */
    { { 'R', 'a' }, 358, 642, 12, "RADIUM"       },  /*  88 */
    { { 'A', 'c' }, 295, 530, 12, "ACTINIUM"     },  /*  89 */
    { { 'T', 'h' }, 255, 460, 12, "THORIUM"      },  /*  90 */
    { { 'P', 'a' }, 222, 400, 12, "PROTACTINIUM" },  /*  91 */
    { { 'U', ' ' }, 242, 437, 12, "URANIUM"      },  /*  92 */
    { { 'N', 'p' }, 238, 427, 12, "NEPTUNIUM"    },  /*  93 */
    { { 'P', 'u' }, 232, 417, 12, "PLUTONIUM"    },  /*  94 */
    { { 'A', 'm' }, 230, 415, 12, "AMERICIUM"    },  /*  95 */
    { { 'C', 'm' }, 228, 412, 12, "CURIUM"       },  /*  96 */
    { { 'B', 'k' }, 225, 410, 12, "BERKELIUM"    },  /*  97 */
    { { 'C', 'f' }, 222, 407, 12, "CALIFORNIUM"  },  /*  98 */
    { { 'E', 's' }, 220, 405, 12, "EINSTEINIUM"  },  /*  99 */
    { { 'F', 'm' }, 218, 402, 12, "FERMIUM"      },  /* 100 */
    { { 'M', 'd' }, 215, 400, 12, "MENDELEVIUM"  },  /* 101 */
    { { 'N', 'o' }, 212, 397, 12, "NOBELIUM"     },  /* 102 */
    { { 'L', 'r' }, 210, 395, 12, "LAWRENCIUM"   }   /* 103 */ /* Lw? */
};

#endif

typedef	struct atom_type  *atm_typ;

typedef struct atom_type {
	BooLean		atom;		/** ATOM or HETERO: 1-6.  **/
	Int4		ser_num;	/** serial number **/
	char		name[5];	/** atom name **/
	char		lowername[5];	/** atom name **/
	char		alt;		/** alternate location indicator **/
	char		res_name[4];	/** residue name **/
	Int4		res;		/** residue number **/
	char		chain;		/** chain identifier **/
	char		insrt;		/** code for inserted residues **/
	float		X;
	float		Y;
	float		Z;
	float		occ;		/** occupancy **/
	float		temp;		/** temperature **/
/*---------------------------------------------------------------------*/
	BooLean		side;		/** is amino acid sidechain atom **/
/*---------------------------------------------------------------------*/
//**** From Hata code 
        float           surface;        /** surface **/
        float           surface_ratio;  /** surface accessability **/
        atm_typ         hbond;
        //float           pradius;        /** VdW radius plus probe radius **/
//**** End from Hata code.
} ;		/** 36 bytes per atom = 360 kbytes for 10,000 atoms ***/


/******************************** private **********************************/
void    atom_error(const char *s);

/******************************** PUBLIC ***********************************/
atm_typ	MkAtom(char *str);
void    PutAtom(FILE *fptr, atm_typ A);
float   DistanceAtom(register float X, register float Y, register float Z,
        register atm_typ A);
float   DistanceAtoms(register atm_typ A1, register atm_typ A2);
char    GetResidueAtom(atm_typ Atom, register a_type A, register a_type nA,
	register BooLean *na);
BooLean IsAtom(const char *str, register atm_typ A);
BooLean	IsWaterAtom(register atm_typ A);
BooLean IsAminoAcidAtom(register atm_typ A);
BooLean IsResAtom(const char *str, register atm_typ A);
BooLean IsBackboneHydrogenAtom(register atm_typ A);
double	CalcAngle(atm_typ atm1, atm_typ atm2, atm_typ atm3);
atm_typ	VirtualAtom(float X, float Y,float Z);
atm_typ	AverageAtom(Int4 natom, atm_typ *atm);
atm_typ CopyCoreAtom(atm_typ atm);
void    NilAtom(atm_typ A);

BooLean ChangeMSE2MetAtom(atm_typ atm);
BooLean IsNucleicAcidAtom(atm_typ atm);

void    TransformAtom(double xyz[4][5],atm_typ atm);

//********************* Surface Calculations *****************


/******************************** MACROS ***********************************/
#define	ResAtom(A)		((A)->res)
#define RenumberAtom(r,A)	((A)->res = (r))
#define ReNameChainAtom(c,A)	((A)->chain= (c))
#define HetAtm2Atom(A)		((A)->atom=TRUE)
#define	IsHeteroAtom(A)		(!(A)->atom)
#define	AminoAcidAtom(A)	((A)->atom)
#define SideAtom(A)		((A)->side)
#define	IsBackboneAtom(A)	(!SideAtom(A))
#define AtomID(A)		((A)->ser_num)
#define AtomChain(A)		((A)->chain)
#define AtomName(A)		((A)->name)
#define AtomLowerName(A)	((A)->lowername)
#define AtomResName(A)		((A)->res_name)
#define SetAtomX(x,A)		((A)->X=(x))
#define SetAtomY(y,A)		((A)->Y=(y))
#define SetAtomZ(z,A)		((A)->Z=(z))
#define AtomX(A)		((A)->X)
#define AtomY(A)		((A)->Y)
#define AtomZ(A)		((A)->Z)
#define AtomIndex(A)		((A)->name[2])
#define AtomIndexDigit(A)	((A)->name[3])
#define CharacterAtom(A)	((A)->name[1])
#define OxygenAtom(A)		((A)->name[1]=='O')
#define NitrogenAtom(A)		((A)->name[1]=='N')
#define CarbonAtom(A)		((A)->name[1]=='C')
#define AlphaCarbonAtom(A)	((A)->name[1]=='C' && (A)->name[2]=='A')
#define BetaCarbonAtom(A)	((A)->name[1]=='C' && (A)->name[2]=='B')
#define HydrogenAtom(A)		(((A)->name[0]==' ' && (A)->name[1]=='H') || (A)->name[0]=='H')
// #define HydrogenAtom(A)		(((A)->name[1]=='H'))
#define SulfurAtom(A)		((A)->name[1]=='S')

//***************************** From Hata code ******************************

#include "radius.h"

/***  private functions ************************************************/
void  SetupPlaneData (atm_typ atom, float *xydists, float *xydists_sq,
                     float *betas, atm_typ *neighbors, int Ncount);
void  CalculateAtomSurface (atm_typ atom, atm_typ *neighbors, int Ncount,
                            float *area, float *ratio);
int	insert_arc (float *arcs, float start, float end, int arcsnum);
float sum_of_arcs(float *arcs, int arcsnum);
#if 0
inline void SetAtomSurface(atm_typ atom, float s, float s_ratio){
            atom->surface = s;    atom->surface_ratio = s_ratio; };
inline void SetAtomHB(atm_typ atom, atm_typ hbatom){
  atom->hbond = new atom_type; atom->hbond = hbatom;};
//inline void SetAtomRadius(atm_typ atom, float r){ atom->pradius = r; }
#endif
/******************************** Public ***********************************/
BooLean IsAtomType(const char name, register atm_typ A);
inline float GetAtomSurface(atm_typ A) { return A->surface; };
inline float GetAtomSurfaceRatio(atm_typ A){ return A->surface_ratio; };
float  EffectiveRadius(atm_typ A);
float  GetVDWRadius(atm_typ A);
void	AtomRealName(char *realname,atm_typ A);
inline atm_typ GetAtomHB(atm_typ atom){return atom->hbond;};
inline BooLean IsAtomHB(atm_typ atom){if (GetAtomHB(atom) != NULL) return TRUE;
                                                                   return FALSE;};
float  GetCosineOfAtomVectors(atm_typ A, atm_typ B, atm_typ C);
/***********************************************************************/
/*
static float Probe      = 1.40;
static float RadiusCA   = 1.70;
static float RadiusO    = 1.52;
static float RadiusN    = 1.55;
static float RadiusC    = 1.80;
static float RadiusFE   = 0.64;
static float RadiusOthers = 1.80;
*/
/***********************************************************************/


#endif

