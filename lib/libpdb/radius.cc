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

#include <stdlib.h>
#include "radius.h"

typedef struct char8 { char name[8]; } char8;

static int MAX_ATOM_TYPE = 120;    // used to be #define

static double atom_vdw[120];
static double atom_covalent[120];
static double atom_charge[120];

static char8 atom_type_names[120] = {0};
static char8 atom_name_for_residue ={0};

static float probe_radius = 1.30;       // For water use 1.2-1.4 Angstroms.

/*------------------------------------------------------------------------*/

FILE *stream_from_string (char *string)
{

        FILE * tt = tmpfile() ;

        fwrite(string,1,strlen(string), tt) ;

        rewind (tt);

        return tt ;
}

/*------------------------------------------------------------------------*/

char * atom_type_name(int index)
{
  return atom_type_names[index].name ;
}

/*------------------------------------------------------------------------*/

int find_type_definition (char *name)
{
  int ii ;
  
  /* ii used to start at 0 ! */
  for (ii =1 ; ii < MAX_ATOM_TYPE ; ii++) {
    char *k =  atom_type_names[ii].name ;
    if (k == NULL ) break ;
    if(!strncmp(name, k,8) )
      return ii ;
  }
  //serious_error(1, "Unknown atom type", name);
  printf("Unknown atom type  %d\n", name);
  return 0 ;
}


/*------------------------------------------------------------------------*/

int allocate_atom_type(char *name)
{
  
  int ii ;
  
  /* ii used to start at 0 ! */
  for (ii =1 ; ii < MAX_ATOM_TYPE ; ii++) {
    char *k =  atom_type_names[ii].name ;
    if (k[0] == '\0') {
      strncpy(k,name,8);
      return ii ;
    }
    if(!strncmp(name, k ,8) )
      return ii ;
  }
  //serious_error(1, "Too many atom types","") ; 
  printf("Too many atom types\n"); 
  return 0 ;
}

/*------------------------------------------------------------------------*/

void read_atom_types (FILE *ff)
{
  char string[200] ; 
  
  while (fgets(string,200,ff) != NULL) {
    char name [200] ; 
    double charge,vdw,covalent;
    int ii ;
    
    //ii = sscanf(string, LF_STRING_02, name,&vdw,&covalent,&charge) ;
    ii = sscanf(string, "%s%lf%lf%lf", name, &vdw, &covalent, &charge) ;

    if (ii > 0) {
      if (ii < 3)
        //serious_error(1, "Something wrong in atom type", string);
	printf("Something wrong in atom type  %s\n", string);
      
      ii = allocate_atom_type(name) ;
      atom_charge[ii] = charge;
      atom_covalent[ii] = covalent;
      atom_vdw[ii] = vdw;
      //printf("name %s  charge %f vdw %f\n", name, charge, vdw);
    }
  }
}

/*------------------------------------------------------------------------*/

char *GetAtomNotation(char *residuename, char *atomname)
{

  print_error("This code needs to be fixed: it's broken and not yet used.");
  // Hata got rid of spaces ' ' which is breaking other parts of the code.
  char *name = atom_name_for_residue.name;
  //char name[6];

  if ( !strncmp(residuename, "GLY", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
  if ( !strncmp(residuename, "ALA", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "VAL", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "CG1") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CG2") ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "LEU", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "CD1") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CD2") ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "ILE", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "CG1") ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG2") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CD1") ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "PRO", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CD")  ) strcpy(name,  "CPRO");
  if ( !strncmp(residuename, "MET", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "SD")  ) strcpy(name,  "S2");
    else if ( !strcmp(atomname, "CE")  ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "PHE", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CD1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CD2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CE1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CE2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CZ")  ) strcpy(name,  "C3H");
  if ( !strncmp(residuename, "TYR", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CD1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CD2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CE1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CE2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CZ")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "OH")  ) strcpy(name,  "O2H");
  if ( !strncmp(residuename, "TRP", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CD1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "NE1") ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CD2") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CE2") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CE3") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CZ3") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CZ2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CH2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CEH2")) strcpy(name,  "CH2");
  if ( !strncmp(residuename, "SER", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "OG")  ) strcpy(name,  "O2H");
    else if ( !strcmp(atomname, "OG1") ) strcpy(name,  "=OG");
  if ( !strncmp(residuename, "THR", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "OG1") ) strcpy(name,  "O2H");
    else if ( !strcmp(atomname, "CG2") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "=CG2");
  if ( !strncmp(residuename, "ASN", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "OD1") ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "ND2") ) strcpy(name,  "N3HH");
  if ( !strncmp(residuename, "GLN", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CD")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "OE1") ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "NE2") ) strcpy(name,  "N3HH");
  if ( !strncmp(residuename, "CYS", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "SG")  ) strcpy(name,  "S2");
  if ( !strncmp(residuename, "HIS", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "ND1") ) strcpy(name,  "N3");
    else if ( !strcmp(atomname, "CD2") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CE1") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "NE2") ) strcpy(name,  "N3H");
  if ( !strncmp(residuename, "GLU", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CD")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "OE1") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OE2") ) strcpy(name,  "O1O2H");
  if ( !strncmp(residuename, "ASP", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "OD1") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OD2") ) strcpy(name,  "O1O2H");
  if ( !strncmp(residuename, "ARG", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CD")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "NE")  ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CZ")  ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "NH1") ) strcpy(name,  "N3HH");
    else if ( !strcmp(atomname, "NH2") ) strcpy(name,  "N3HH");
  if ( !strncmp(residuename, "LYS", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "O1");
    else if ( !strcmp(atomname, "C")   ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CA")  ) strcpy(name,  "C4H");
    else if ( !strcmp(atomname, "N")   ) strcpy(name,  "N3H");
    else if ( !strcmp(atomname, "CB")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CG")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CD")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CE")  ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "NZ")  ) strcpy(name,  "N4HHH");
  if ( !strncmp(residuename, "ZN", 2) )
    if      ( !strcmp(atomname, "ZN")  ) strcpy(name,  "Z2ION");
  if ( !strncmp(residuename, "HEM", 3) )
    if      ( !strcmp(atomname, "FE")  ) strcpy(name,  "FE");
    else if ( !strcmp(atomname, "CMA") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CMB") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CMC") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CMD") ) strcpy(name,  "C4HHH");
    else if ( !strcmp(atomname, "CAA") ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CBA") ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CAD") ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CBD") ) strcpy(name,  "C4HH");
    else if ( !strcmp(atomname, "CAB") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CH")  ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CAC") ) strcpy(name,  "C3H");
    else if ( !strcmp(atomname, "CBC") ) strcpy(name,  "C3HH");
    else if ( !strcmp(atomname, "CBB") ) strcpy(name,  "C3HH");
    else if ( !strcmp(atomname, "NA")  ) strcpy(name,  "N3");
    else if ( !strcmp(atomname, "NB")  ) strcpy(name,  "N3");
    else if ( !strcmp(atomname, "NC")  ) strcpy(name,  "N3");
    else if ( !strcmp(atomname, "ND")  ) strcpy(name,  "N3");
    else if ( !strcmp(atomname, "O1A") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "O1B") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "O1C") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "O1D") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "C1A") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C2A") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C3A") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C4A") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "CGA") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C1B") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C2B") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C3B") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C4B") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C1C") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C2C") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C3C") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C4C") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C1D") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C2D") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C3D") ) strcpy(name,  "C3");
    else if ( !strcmp(atomname, "C4D") ) strcpy(name,  "C3");
  if ( !strncmp(residuename, "ACE", 3) )
    if      ( !strcmp(atomname, "CH3") ) strcpy(name,  "C4HHH");
  if ( !strncmp(residuename, "ASH", 3) )
    if      ( !strcmp(atomname, "CH3") ) strcpy(name,  "");
  if ( !strncmp(residuename, "WAT", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "OWAT");
    else if ( !strcmp(atomname, "OW")  ) strcpy(name,  "=O");
  if ( !strncmp(residuename, "HOH", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "OWAT");
    else if ( !strcmp(atomname, "OW")  ) strcpy(name,  "=O");
  if ( !strncmp(residuename, "DOD", 3) )
    if      ( !strcmp(atomname, "O")   ) strcpy(name,  "OWAT");
    else if ( !strcmp(atomname, "OW")  ) strcpy(name,  "=O");
  if ( !strncmp(residuename, "DUM", 3) )
    if      ( !strcmp(atomname, "OT")  ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OXT") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OT1") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OT2") ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "OE")  ) strcpy(name,  "O1O2H");
    else if ( !strcmp(atomname, "NT")  ) strcpy(name,  "N4HHH");
    else if ( !strcmp(atomname, "Z")   ) strcpy(name,  "Z2ION");
    else if ( !strcmp(atomname, "DUM") ) strcpy(name,  "DEFAULT");
    else if ( !strcmp(atomname, "NDUM")) strcpy(name,  "N4HHH");
    else if ( !strcmp(atomname, "ODUM")) strcpy(name,  "O2H");
    else if ( !strcmp(atomname, "CDUM")) strcpy(name,  "C4HH");

  //printf("name %s\n", name);
  return name;
}

/*------------------------------------------------------------------------*/

void InitializeAtomRadii()
{
  char input_string[5000];  
  strcpy(input_string,  RICHARDS_RADII);
  read_atom_types ( stream_from_string(input_string) );

  return;
}

/*------------------------------------------------------------------------*/

float GetAtomVDWRadius(char *resname, char *atomname)
{
  char  *atom_notation;
  float vdw_radius;

#if 1
  vdw_radius = 1.80;	// use default for now (need to implement RICHARDS_RADII!
#else
  atom_notation = GetAtomNotation(resname, atomname); 
  vdw_radius    = atom_vdw[allocate_atom_type(atom_notation)];
#endif

  return vdw_radius;
}

/*------------------------------------------------------------------------*/

float GetAtomEffectiveRadius(char *resname, char *atomname)
{
  char  *atom_notation;
  float vdw_radius;

#if 0
  atom_notation = GetAtomNotation(resname, atomname); 
  vdw_radius    = atom_vdw[allocate_atom_type(atom_notation)];
  fprintf(stderr,"vdw_radius = %.3f\n",vdw_radius);
#else
  vdw_radius = 1.80;    // use default for now (need to implement RICHARDS_RADII!
#endif
  return vdw_radius + probe_radius;
}

/*------------------------------------------------------------------------*/
