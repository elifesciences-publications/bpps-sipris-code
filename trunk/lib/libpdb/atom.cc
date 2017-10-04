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

atm_typ MkAtom(char *str)
/**********************************************************************
TER		
HELIX    1     ARG      8  LYS     14 
SHEET    1       GLY    23  GLN    29
CONECT  144  667
 **********************************************************************
ATOM     62  NH1 ARG A   8      -7.802  -5.691  29.031  1.00 27.67      1AAK N
HETATM 1205  O   HOH W 151      -5.065  11.390  30.224  1.00  3.66      1AAK O

<----><---> <-->^<-> ^<-->^   <---X--><---Y--><---Z--><----><----> <->  <--><><>
....|....|....|....|....|....|....|....|....|....|....|....|....|....|....|....|
    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75   80
{
typedef struct {
	BooLean		atom;		= ATOM or HETERO: 1-6.
	Int4		ser_num;	= serial number
	char		name[4];	= atom name
	char		alt;		= alternate location indicator
	char		res_name[4];	= residue number 
	char		chain;		= chain identifier
	char		insrt;		= code for inserted residues
	float		X;
	float		Y;
	float		Z;
	float		occ;		= occupancy
	float		temp;		= temperature 
} atom_type;		36 bytes per atom = 360 kbytes for 10,000 atoms
/**********************************************************************/
{
	atm_typ	A=NULL;
	char	tmp[50];
	Int4	i,j,k,lenstr;

	if(strncmp(str,"HETATM",6) == 0) {
		NEW(A,1,atom_type); A->atom = FALSE;
	} else if(strncmp(str,"ATOM",4) == 0) {
		NEW(A,1,atom_type); A->atom = TRUE;
	}
	if(A == NULL) return NULL;
	  /** serial number **/
	lenstr=strlen(str);
	if(strlen(str) < 54) atom_error("input error1");
	for(i=6,j=0; i < 11; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	if(sscanf(tmp,"%d",&A->ser_num) != 1){
	    fprintf(stderr,"Input error: %s...ignored.\n",str);
#if 0
	    atom_error("input error2");
#else
	    return NULL;	// skip these...
#endif
	}
	  /** atom name **/
	for(i=12,j=k=0; i < 16; i++,j++) {
	   A->name[j]= str[i]; 
	   if(!isspace(str[i])){ A->lowername[k] = tolower(str[i]); k++; }
	} A->name[j]=0; A->lowername[j]=0;
	  /** alternate location **/
	A->alt = str[16];
	  /** residue name **/
	for(i=17,j=0; i < 20; i++,j++) A->res_name[j]= str[i]; A->res_name[j]=0;
	  /** chain designator **/
	A->chain = str[21];
	  /** residue number **/
	for(i=22,j=0; i < 26; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	if(sscanf(tmp,"%d",&A->res) != 1){
		fprintf(stderr,"%s\n",str);
		atom_error("input error3");
	}
	  /** insert designator **/
	A->insrt= str[26];

	  /** X coordinate **/
	for(i=30,j=0; i < 38; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	if(sscanf(tmp,"%f",&A->X) != 1) atom_error("input error4");

	  /** Y coordinate **/
	for(i=38,j=0; i < 46; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	if(sscanf(tmp,"%f",&A->Y) != 1) atom_error("input error5");

	  /** Z coordinate **/
	for(i=46,j=0; i < 54; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	if(sscanf(tmp,"%f",&A->Z) != 1) atom_error("input error6");

	if(lenstr >= 60){	// occupancy (Not required).
	  for(i=54,j=0; i < 60; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	  if(sscanf(tmp,"%f",&A->occ) != 1) atom_error("input error7");
	} else A->occ = 0;

	if(lenstr >= 66){	// temperature (not required).
	  for(i=60,j=0; i < 66; i++,j++) tmp[j]= str[i]; tmp[j] = 0;
	  if(sscanf(tmp,"%f",&A->temp) != 1){
		fprintf(stderr,"Input: %s\n",str);
		atom_error("input error8");
	  }
	} else A->temp = 0;

/*** essential info ***/
	if(strncmp(str,"ATOM",4) == 0) {
	  if(A->res_name[0] != ' ') {	/** not DNA or RNA **/
            if(IsAtom(" CA ",A) || IsAtom(" N  ",A) || IsAtom(" C  ",A)
			|| IsAtom(" O  ",A) || IsAtom(" OXT",A)
			|| IsAtom(" H  ",A) || IsAtom(" HA ",A)	// include sidechain hydrogens...
		     ){
                 		 A->side = FALSE;
            } else A->side = TRUE;
	  } else A->side = FALSE;
	}
/*** essential info ***/
	return A;
}

BooLean IsBackboneHydrogenAtom(register atm_typ A)
{ if(IsAtom(" H  ",A) || IsAtom(" HA ",A)) return TRUE; else return FALSE; }

BooLean IsWaterAtom(register atm_typ A)
{
	if(strncmp(A->res_name,"HOH",3)==0) return TRUE;
	else return FALSE;
}

BooLean	IsAminoAcidAtom(register atm_typ Atom)
// Return true if the residue is an amino acid; otherwise return false.
{
	register char *name = Atom->res_name;
	switch(name[0]){
	  case 'A':		/** ALA ARG ASN ASP **/
	    if(strncmp("ALA",name,3) == 0) return TRUE;
	    if(strncmp("ARG",name,3) == 0) return TRUE;
	    if(strncmp("ASN",name,3) == 0) return TRUE;
	    if(strncmp("ASP",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'C':		/** CYS SELENOCYSTEINE **/
	    if(strncmp("CYS",name,3) == 0) return TRUE;
	    if(strncmp("CSE",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'G':		/** GLN GLU GLY **/
	    if(strncmp("GLY",name,3) == 0) return TRUE;
	    if(strncmp("GLU",name,3) == 0) return TRUE;
	    if(strncmp("GLN",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'H':		/** HIS **/
	    if(strncmp("HIS",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'I':		/** ILE **/
	    if(strncmp("ILE",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'L':		/** LEU LYS **/
	    if(strncmp("LEU",name,3) == 0) return TRUE;
	    if(strncmp("LYS",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'M':		/** MET **/
	    if(strncmp("MET",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'P':		/** PRO PHE PTR (phosphotyrosine)**/
	    if(strncmp("PRO",name,3) == 0) return TRUE;
	    if(strncmp("PHE",name,3) == 0) return TRUE;
	    if(strncmp("PTR",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'S':		/** SER SEP=Phosphoserine **/
	    if(strncmp("SER",name,3) == 0) return TRUE;
	    if(strncmp("SEP",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'T':		/** THR TRP TYR THP=phosphothreonine**/
	    if(strncmp("THR",name,3) == 0) return TRUE;
	    if(strncmp("THP",name,3) == 0) return TRUE;
	    if(strncmp("TRP",name,3) == 0) return TRUE;
	    if(strncmp("TYR",name,3) == 0) return TRUE;
	    else return FALSE;
	  case 'V':		/** VAL **/
	    if(strncmp("VAL",name,3) == 0) return TRUE;
	    else return FALSE;
	  default:  return FALSE;
	}
}

BooLean	IsResAtom(const char *str, register atm_typ A)
{
	if(A->res_name[0] != str[0] || A->res_name[1] != str[1]) return FALSE;
	if(A->res_name[2] != str[2]) return FALSE;
	return TRUE;
}

BooLean	IsAtom(const char *str, register atm_typ A)
/** e.g., IsAtom(" CA ",A) **/
{
	if(A->name[0] != str[0] || A->name[1] != str[1]) return FALSE;
	if(A->name[2] != str[2] || A->name[3] != str[3]) return FALSE;
	return TRUE;
}

float	DistanceAtoms(register atm_typ A1, register atm_typ A2)
{ return DistanceAtom(A1->X, A1->Y, A1->Z, A2); }

float	DistanceAtom(register float X, register float Y, register float Z,
	register atm_typ A)
/** return the distance from atom A to point [X,Y,Z] **/
{
	register float	dx,dy,dz;

	dx = X - A->X;
	dy = Y - A->Y;
	dz = Z - A->Z;
	return  sqrt(dx*dx + dy*dy + dz*dz);
}

BooLean	IsNucleicAcidAtom(atm_typ atm)
{
	register char *name = atm->res_name;
	if(!IsHeteroAtom(atm)) return FALSE;
	if(name[0]!=' ') return FALSE;
	if(name[1]!=' ') return FALSE;
	switch(name[2]){
		case 'A': case 'C': case 'G': case 'T': case 'U':
		 return TRUE;
		default: return FALSE;
	}
}

char	GetResidueAtom(atm_typ Atom, register a_type A, register a_type nA,
	register BooLean *na)
{
	register char *name = Atom->res_name;
	*na = FALSE;
	switch(name[0]){
	  case ' ':  
		if(name[1]!=' ') return AlphaCode('X',A);
		*na = TRUE;
		switch(name[2]){
	  	   case 'A':  return AlphaCode('A',nA);
	  	   case 'C':  return AlphaCode('C',nA);
	  	   case 'G':  return AlphaCode('G',nA);
	  	   case 'T':  return AlphaCode('T',nA);
	  	   case 'U':  return AlphaCode('U',nA);
		   default: return AlphaCode('N',nA);
		}
	  case 'A':		/** ALA ARG ASN ASP **/
		switch(name[1]){
		   case 'L': return AlphaCode('A',A);
		   case 'R': return AlphaCode('R',A);
		   case 'S': 
		     switch(name[2]){
			case 'N': return AlphaCode('N',A);
			case 'P': return AlphaCode('D',A);
			default:  return AlphaCode('X',A);
		     }
		   default:  return AlphaCode('X',A);
		}
	  case 'C':		/** CYS **/
		if(name[1]=='Y' && name[2] =='S') 
			return AlphaCode('C',A);
		else return AlphaCode('X',A);
	  case 'G':		/** GLN GLU GLY **/
	    if(name[1] == 'L'){
		switch(name[2]){
		   case 'N': return AlphaCode('Q',A);
		   case 'U': return AlphaCode('E',A);
		   case 'Y': return AlphaCode('G',A);
		   default:  return AlphaCode('X',A);
		}
	    } else return AlphaCode('X',A);
	  case 'H':		/** HIS **/
		if(name[1]=='I' && name[2] =='S') 
			return AlphaCode('H',A);
	    	else return AlphaCode('X',A);
	  case 'I':		/** ILE **/
		if(name[1]=='L' && name[2] =='E') 
			return AlphaCode('I',A);
	    	else return AlphaCode('X',A);
	  case 'L':		/** LEU LYS **/
		switch(name[1]){
		   case 'E': 
			if(name[2]== 'U') return AlphaCode('L',A);
	    		else return AlphaCode('X',A);
		   case 'Y': 
			if(name[2]== 'S') return AlphaCode('K',A);
		   default:  return AlphaCode('X',A);
		}
	  case 'M':		/** MET **/
		if(name[1]=='E' && name[2] =='T')
                        return AlphaCode('M',A);
	    	else return AlphaCode('X',A);
	  case 'P':		/** PRO PHE PTR **/
		switch(name[1]){
		   case 'T': 	// PTR = phosphotyrosine
			if(name[2]=='R') return AlphaCode('Y',A);
	    		else return AlphaCode('X',A);
		   case 'R': 
			if(name[2]=='O') return AlphaCode('P',A);
	    		else return AlphaCode('X',A);
		   case 'H': 
			if(name[2]=='E') return AlphaCode('F',A);
		   default:  return AlphaCode('X',A);
		}
	  case 'S':		/** SER SEP SEP **/
		if(name[1]=='E' && name[2] =='R')
                        return AlphaCode('S',A);
		if(name[1]=='E' && name[2] =='P') // phosphoserine
                        return AlphaCode('S',A);
	    	else return AlphaCode('X',A);
	  case 'T':		/** THR TRP TYR THP (phosphothreonine) **/
		switch(name[1]){
		   case 'H': 
			if(name[2]== 'R') return AlphaCode('T',A);
			if(name[2]== 'P') return AlphaCode('T',A);
	    		else return AlphaCode('X',A);
		   case 'R': 
			if(name[2]== 'P') return AlphaCode('W',A);
	    		else return AlphaCode('X',A);
		   case 'Y': 
			if(name[2]== 'R') return AlphaCode('Y',A);
		   default:  return AlphaCode('X',A);
		}
	  case 'V':		/** VAL **/
		if(name[1]=='A' && name[2] =='L')
                        return AlphaCode('V',A);
	    	else return AlphaCode('X',A);
	  default:  return AlphaCode('X',A);
	}
}

atm_typ VirtualAtom(float X, float Y,float Z)
// Create and return a virtual atom.
{
	char str[100];
	sprintf(str,"ATOM      1  C   XXX X   1    %8.3f%8.3f%8.3f  1.00  0.00      1XXX",
		X,Y,Z);
	return MkAtom(str);
}

atm_typ AverageAtom(Int4 natom, atm_typ *atm)
{
	double aveX=0.0,aveY=0.0,aveZ=0.0;
	for(Int4 a=1; a <= natom; a++){
		aveX+=AtomX(atm[a]);
		aveY+=AtomY(atm[a]);
		aveZ+=AtomZ(atm[a]);
	}
	return VirtualAtom(aveX/(float)natom,aveY/(float)natom,aveZ/(float)natom);
}

void    TransformAtom(double xyz[4][5],atm_typ atm)
/*******************************************************************************
                  xyz[.][1]       xyz[.][2]        xyz[.][3]        xyz[.][4]
xyz[1]:   X2 = ( 0.983195)*X1 + ( 0.109599)*Y1 + ( 0.146000)*Z1 + (  -13.846273)
xyz[2]:   Y2 = (-0.059883)*X1 + ( 0.949108)*Y1 + (-0.309206)*Z1 + (   14.920112)
xyz[3]:   Z2 = (-0.172458)*X1 + ( 0.295267)*Y1 + ( 0.939721)*Z1 + (  -35.552807)
 *******************************************************************************/
{
        double  X,Y,Z;
        X=atm->X; Y = atm->Y; Z = atm->Z;
        atm->X = xyz[1][1]*X + xyz[1][2]*Y + xyz[1][3]*Z + xyz[1][4];
        atm->Y = xyz[2][1]*X + xyz[2][2]*Y + xyz[2][3]*Z + xyz[2][4];
        atm->Z = xyz[3][1]*X + xyz[3][2]*Y + xyz[3][3]*Z + xyz[3][4];
}

BooLean	ChangeMSE2MetAtom(atm_typ atm)
// Returns True if was MSE and changed to Met.
// If atoms Selenium (= SE ) then change to "  SD" 
#if 0 //************************************************************
HETATM  956  N   MSE A  86      25.636  22.984  36.896  1.00 24.72        
HETATM  957  CA  MSE A  86      25.796  21.712  36.227  1.00 27.66        
HETATM  958  C   MSE A  86      25.447  20.546  37.134  1.00 26.17        
HETATM  959  O   MSE A  86      24.680  19.666  36.748  1.00 26.91        
HETATM  960  CB  MSE A  86      27.213  21.580  35.666  1.00 29.94        
HETATM  961  CG  MSE A  86      27.496  22.571  34.581  1.00 34.98        
HETATM  962 SE   MSE A  86      28.739  21.923  33.274  1.00 44.76        
HETATM  963  CE  MSE A  86      27.494  20.948  32.165  1.00 38.63        
   Change to:
ATOM   1943  N   MET A 130       9.128  17.706  29.550  1.00 23.84           N
ATOM   1944  CA  MET A 130       8.821  18.332  28.266  1.00 24.11           C
ATOM   1945  C   MET A 130      10.035  18.996  27.616  1.00 24.20           C
ATOM   1946  O   MET A 130       9.910  19.672  26.598  1.00 25.67           O
ATOM   1947  CB  MET A 130       7.660  19.329  28.433  1.00 24.78           C
ATOM   1948  CG  MET A 130       6.295  18.639  28.618  1.00 27.09           C
ATOM   1949  SD  MET A 130       4.908  19.726  29.095  1.00 33.88           S
ATOM   1950  CE  MET A 130       4.520  20.519  27.520  1.00 31.23           C
^^^^^^      ^^^^ ^^^  (change these fields only)
#endif //*********************************************************
{
	Int4	i;
	if(!atm) return FALSE;
	if(!IsHeteroAtom(atm)) return FALSE;
	if(atm->res_name[0] != 'M') return FALSE;
	if(atm->res_name[1] != 'S') return FALSE;
	if(atm->res_name[2] != 'E') return FALSE;
	atm->atom = TRUE;
	strncpy(atm->res_name,"MET",4);	// rename residue to MET.
	if(strncmp(atm->name,"SE  ",4) == 0){
		strncpy(atm->name," SD ",5);
	}
	return TRUE;
}

atm_typ CopyCoreAtom(atm_typ atm)
// Copy the Core files in atm to new_atom
{
	Int4	i;
	atm_typ	new_atm;
	if(!atm) return 0;
	NEW(new_atm,1,atom_type); 
	new_atm->atom = atm->atom;
	new_atm->ser_num = atm->ser_num;
	for(i=0; i < 5; i++) new_atm->name[i]=atm->name[i];
	new_atm->alt = atm->alt;
	for(i=0; i < 4; i++) new_atm->res_name[i]=atm->res_name[i];
	new_atm->res= atm->res;
	new_atm->chain= atm->chain;
	new_atm->insrt= atm->insrt;
	new_atm->X= atm->X;
	new_atm->Y= atm->Y;
	new_atm->Z= atm->Z;
	new_atm->occ= atm->occ;
	new_atm->temp= atm->temp;
	return new_atm;
}

void	PutAtom(FILE *fptr, atm_typ A)
{
	Int4	i;
	if(IsAminoAcidAtom(A)) fprintf(fptr,"ATOM  ");
	else fprintf(fptr,"HETATM"); 
#if 0	// above checks to make sure designation is correct.
	if(A->atom) fprintf(fptr,"ATOM  ");
	else fprintf(fptr,"HETATM"); /** serial number **/
#endif
	fprintf(fptr,"%5d ",A->ser_num); /** atom name **/
	for(i=0; i < 4; i++) fprintf(fptr,"%c",A->name[i]);
	  /** alternate location **/
	fprintf(fptr,"%c",A->alt); /** residue name **/
	for(i=0; i < 3; i++) fprintf(fptr,"%c",A->res_name[i]);
	  /** chain designator **/
	fprintf(fptr," %c",A->chain);	// residue number 
	fprintf(fptr,"%4d",A->res);	// insert designator 
	fprintf(fptr,"%c   ",A->insrt); 	
	fprintf(fptr,"%8.3f",A->X);	// X coordinate 
	fprintf(fptr,"%8.3f",A->Y);	// Y coordinate 
	fprintf(fptr,"%8.3f",A->Z);	// Z coordinate 
	// if(A->occ > 0) fprintf(fptr,"%6.2f",A->occ);	// occupancy 
	fprintf(fptr,"%6.2f",A->occ);	// occupancy 
	if(A->temp > 0) fprintf(fptr,"%6.2f        ",A->temp);	// temperature
	fprintf(fptr,"\n");
}

void	NilAtom(atm_typ A) { free(A); }

void    atom_error(const char *s){fprintf(stderr,"atom: %s\n",s);exit(1);}

//***************************** From Hata code ****************************
/*------------------------------------------------------------------------*/
static void SetAtomSurface(atm_typ atom, float s, float s_ratio){
            atom->surface = s;    atom->surface_ratio = s_ratio; };
static void SetAtomHB(atm_typ atom, atm_typ hbatom){
  atom->hbond = new atom_type; atom->hbond = hbatom;};
//inline void SetAtomRadius(atm_typ atom, float r){ atom->pradius = r; }

   
//static float zstep = 0.25 ; 
static float zstep = 0.15 ; 
static float two_pi = M_PI + M_PI ; 
static int   MaxArcNumbers = 500;

void CalculateAtomSurface (atm_typ atom, atm_typ *neighbors, int Ncount,
			   float *surface, float *surface_ratio)
{
  float xx = AtomX(atom), yy = AtomY(atom), zz = AtomZ(atom);
  float xydists[500], betas[500], xydists_sq[500]; 
  float eradius = EffectiveRadius(atom);
  float eradsq = eradius * eradius;
  float total_arc = 0.0;

  int   ii, jj; 
  int   grid_points_number = (int) ( 2.0*eradius/zstep + 0.5 );
  float zgrid = zz - eradius - zstep * 0.5;
  int   arcsnum;
  float arcs[500]; 

  if (Ncount < 1) {
    *surface = 2.0 * two_pi * eradsq;
    *surface_ratio = 1.00;
    return;
  }    

  SetupPlaneData(atom, xydists, xydists_sq, betas, neighbors, Ncount);
  
  for (ii = 0; ii < grid_points_number; ii++) {
    zgrid = zgrid + zstep;
    {
      float aa = zgrid - zz; 
      float plane_radius_sq = eradsq - aa*aa;
      float plane_radius;
      arcsnum = 0;
      
      if (plane_radius_sq < 0.0) plane_radius_sq = 0.000001;
      plane_radius = sqrt(plane_radius_sq);
      
      for (jj = 0; jj < Ncount; jj++) {
	{
	  atm_typ neighbor = neighbors[jj]; 
	  float neighbor_plane_radius;
	  float aa = zgrid - AtomZ(neighbor);
	  float ERn = EffectiveRadius(neighbor);
	  float neighbor_plane_radius_sq = ERn * ERn - aa * aa;
	  if (neighbor_plane_radius_sq > 0.0) {
	    float neighbor_plane_radius = sqrt(neighbor_plane_radius_sq);
	    float xydist = xydists[jj];
	    
	    if (xydist < neighbor_plane_radius + plane_radius) {
	      float temp = plane_radius - neighbor_plane_radius;
	      
	      if (xydist <= fabs(temp)) {
		if (temp > 0)
		  continue;
		//else goto next_plane;
		else goto next_neighbor;
	      }
	      else {
		float alpha;
		double d = (xydists_sq[jj] + plane_radius_sq -
				neighbor_plane_radius_sq)/(2.0 * xydist * plane_radius);
		if(d < -1 || d > 1.){
		    if((fabs(d)-1.0) > 0.000000015){
		       	fprintf(stderr,"CalculateAtomSurface( ) error\n");
		       	fprintf(stderr,"\tacos(d) where d = %.9f\n",d);
			fprintf(stderr,"\trounding to nearest in-range value...\n");
		    } if(d < -1.0) d = -1.0; else if(d > 1.0) d = 1.0;
		} alpha = acos(d);
		float beta = betas[jj];
		arcsnum = insert_arc(arcs, beta-alpha, beta+alpha, arcsnum);
	      }
	    }
	  }
	next_neighbor: ;
	}
      }
      total_arc += two_pi - sum_of_arcs(arcs, arcsnum);
    next_plane: ;
    }
  }
  *surface       = total_arc * zstep *  eradius;
  *surface_ratio = *surface / (2.0 * two_pi * eradsq);

  return;
}

/*------------------------------------------------------------------------*/


/* Insert an arc into the arcs array. each arc is two angles which are
   moved into the range -pi to pi, and arcsnum is the number of
   arcs. If the arc go over pi, it is cut to two arcs.  Before
   actually inserting we check if it overlaps with any existing arc,
   and if it is extend the other arc to include both, remove it from
   arcs and reinsert it retunrn the number of arcs after the addition. */

int insert_arc (float *arcs, float start, float end, int arcsnum)
{
  //  float small = 0.00000001;
  float small = 0.00001;

  if (start < (- M_PI - small) ) start += two_pi;
  if (end >  M_PI + small) end -= two_pi;
  if (start > end)
    {
      arcsnum = insert_arc(arcs, (- M_PI) , end, arcsnum) ;
      return insert_arc(arcs, start, M_PI, arcsnum) ;
    }
  else
    {
      int ii ;
      for(ii=0; ii< arcsnum; ii += 2)
	{  float e_start = arcs[ii] ;
	float e_end  = arcs[ii+1] ;
	
	if(start > e_start)
	  {
	    if (start <= e_end)
	      {
		if (end > e_end)
                  {
		    arcs[ii] = arcs[arcsnum-2];
		    arcs[ii+1] = arcs[arcsnum-1];
		    return insert_arc(arcs, e_start, end, arcsnum-2);
                  }
	        else return arcsnum;
              }
	  }
	else if (end >= e_start)
	  {
	    if (end < e_end)  end = e_end;
	    arcs[ii] = arcs[arcsnum-2];
	    arcs[ii+1] = arcs[arcsnum-1];
	    return insert_arc(arcs, start, end, arcsnum-2);
	  }
	}
      arcs[arcsnum++] = start;
      arcs[arcsnum++] = end;
      return arcsnum;
    }
}

/*------------------------------------------------------------------------*/

float sum_of_arcs(float *arcs, int arcsnum)
{
  int ii = 0;
  float sum = 0.0;

  for (ii = 0; ii < arcsnum; ii += 2)
      sum+= arcs[ii+1] - arcs[ii];

  return sum;
}

/*------------------------------------------------------------------------*/

/* 
    This setup all the relations betwen atom and its neibours in the
    xy plane, in parallel arrays.
*/

void SetupPlaneData (atm_typ atom, float *xydists, float *xydists_sq,
		       float *betas, atm_typ *neighbors, int NeighborCount)
{
   int ii;
   float xx = AtomX(atom);
   float yy = AtomY(atom);
   float zz = AtomZ(atom);
   float deltax, deltay, beta, sq;

   for (ii = 0 ; ii < NeighborCount ; ii++)
     {
	atm_typ   neighbor = neighbors[ii] ;
	deltax =  AtomX(neighbor) - xx;
	deltay =  AtomY(neighbor) - yy;
	beta = atan2(deltay,deltax) ;
        betas[ii] = beta ;
	sq = deltax * deltax + deltay * deltay ;
	xydists[ii] = sqrt(sq) ;
	xydists_sq[ii] = sq;
    }
   return;
}

/*------------------------------------------------------------------------*/

BooLean IsAtomType(const char name, register atm_typ A)
{
  int   i = 0, j = 0;

  if ( A->name[1] == name) return TRUE;

  return FALSE;
}    

/*------------------------------------------------------------------------*/
void	AtomRealName(char *realname, atm_typ A)
{
  int  i = 0, j = 0;

  while(A->name[i] == ' ' ) i++;
  while(A->name[i] != ' '  && A->name[i] != 0){
    realname[j] = A->name[i];
    i++; j++;
  }
  if ( j < 10 ) {
    realname[j] = '\0';
    j++;
  }
  fprintf(stderr,"%s and real name:%s\n", A->name, realname);
}    

/*------------------------------------------------------------------------*/

float EffectiveRadius(atm_typ A)
{
  float ER;
  //char *name, *rname;
  char name[10], rname[10];

  //name  = AtomName(A);
  //rname = AtomResName(A);
  // AtomRealName(name,A);
  strcpy(name, AtomName(A));
  strcpy(rname, AtomResName(A));
  /*
  ER = RadiusOthers;

  if (name[1] == 'C')
    if ( name[2] == 'A')
      ER = RadiusCA + Probe;
    else
      ER = RadiusC + Probe;
  if (name[1] == 'N')
    ER = RadiusN + Probe;
  if (name[1] == 'O')
    ER = RadiusO + Probe;
  if (name[1] == 'F' && name[2] == 'E')
    ER = RadiusFE + Probe;
  */

  ER = GetAtomEffectiveRadius(rname, name);

  return ER;
}
     
/*------------------------------------------------------------------------*/

float GetVDWRadius(atm_typ A)
{
  float ER;
  //char *name, *rname;
  char name[10], rname[4];

  // AtomRealName(name,A);
  strcpy(name, AtomName(A));
  strcpy(rname, AtomResName(A));

  ER = GetAtomVDWRadius(rname, name);

  return ER;
}
     
/*------------------------------------------------------------------------*/

float  GetCosineOfAtomVectors(atm_typ A, atm_typ B, atm_typ C)
{
  // Calculate cosine of vec(AC) and vec (BC).

  float ac, bc, product;

  ac = DistanceAtoms(A, C);
  bc = DistanceAtoms(B, C);

  product =  (AtomX(A) - AtomX(C)) * (AtomX(B) - AtomX(C))
            +(AtomY(A) - AtomY(C)) * (AtomY(B) - AtomY(C))
            +(AtomZ(A) - AtomZ(C)) * (AtomZ(B) - AtomZ(C));

  return product/ac/bc;
}
  
#if 1	
#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif
#define Rad2Deg      (180.0/PI)

double CalcAngle(atm_typ atm1, atm_typ atm2, atm_typ atm3)
// From Rasmol code: Rasmol_16bit/rasmol/RasMol2-sun/abstree.c file
{
    register double ulen2,vlen2;
    register double ux,uy,uz;
    register double vx,vy,vz;
    register double temp;

    ux = atm1->X - atm2->X;
    uy = atm1->Y - atm2->Y;
    uz = atm1->Z - atm2->Z;
    if( !ux && !uy && !uz ) return 0.0;
    ulen2 = ux*ux + uy*uy + uz*uz;

    vx = atm3->X - atm2->X;
    vy = atm3->Y - atm2->Y;
    vz = atm3->Z - atm2->Z;
    if( !vx && !vy && !vz ) return 0.0;
    vlen2 = vx*vx + vy*vy + vz*vz;

    temp = (ux*vx + uy*vy + uz*vz)/sqrt(ulen2*vlen2);
    return Rad2Deg*acos(temp);
}

#endif

