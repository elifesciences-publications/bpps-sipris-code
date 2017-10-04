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

/* pml_typ.h - PyMol type output program (for creating PyMol scripts). */
#if !defined (_PML_TYP_)
#define _PML_TYP_

#include <time.h>
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "set_typ.h"

class pml_typ {         // Rasmol type
public:
	pml_typ( ){ assert(!"Illegal constructor"); }
	pml_typ(char *filename,BooLean cis_trans_pro);
	~pml_typ( ){ Free( ); }

	void	PrintHEADER(FILE *fp,BooLean use_path, char c,Int4 wire_width);
	void    PrintTAIL(FILE *fp);
	void    PrintColor(FILE *fp, char c);
	void    ColorAtom(FILE *fp, char *mol, char c, char color, Int4 i, char *atom);
	void	PrintView(FILE *fp, Int4 rx,Int4 ry,Int4 rz,Int4 tx,Int4 ty,Int4 tz);
	char	*DefineResCloud(FILE *fp, char aa, Int4 site, char chain);
	char    *DefineResCloud(FILE *fp, char aa, Int4 site, char *atom1,char *atom2,
					char chain);
	void    PrintTrace(FILE *fp, Int4 start, Int4 end, char c,char color,char mode,
			Int4 diameter);
	void	PutResidueItem(FILE *fp, char a, Int4 i, char *atom1,    
                        char *atom2, char chn, char color,Int4 wire_width,Int4 spacefill);
	void    PutMoleculeItem(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
        		char chn, char color,Int4 wire_width,Int4 spacefill);
	void	PrintCommand(FILE *fp,char *cmd,char type);
private:
	void    PrintResidueAtoms(FILE *fp, char a, Int4 i, char *atom1, 
			char *atom2, char chn, char color,Int4 wire_width,
			Int4 spacefill);
	void    PrintMoleculeAtoms(FILE *fp, char *mol, Int4 i, char *atom1, char *atom2,
        		char chn, char color,Int4 wire_width,Int4 spacefill);

	void    PrintResidueAtom(FILE *fp, char aa, char c, char color, Int4 i,
        		char *atom, BooLean sideonly,Int4 big_spacefill);
	void    PrintMoleculeAtom(FILE *fp, char *mol, char c, char color, Int4 i,
        		char *atom, BooLean sideonly,Int4 big_spacefill);

#if 0	// obsolete functions...
	void    PrintClear(FILE *fp);
	char    *DefineMolCloud(FILE *fp, char *mol, Int4 site, char *atom1,
			char *atom2, char chain);
	void    PrintResidue(FILE *fp, char aa, char c, char color, Int4 i,
        		BooLean sideonly,Int4 Width, Int4 big_spacefill);
	void    PrintMolecule(FILE *fp, char *mol, char c, char color, Int4 i,
			Int4 Width,Int4 spacefill);
	void    PrintTrueColor(FILE *fp, char c);
#endif
//======= NEW ================
public:
	void    DefineResidueItem(FILE *fp,char a,Int4 i,char chn,char color);
private:
	void    PrintResColor(FILE *fp,char *aa,char color);
	set_typ	**Object;
	BooLean	IsObject(Int4 i, char chnX){
		   if(Object[chnX] == 0) return FALSE;
		   for(char cls='A'; cls <= 'Z'; cls++){
			if(Object[chnX][cls] && MemberSet(i,Object[chnX][cls])) return TRUE;
		   } return FALSE;
		}
	Int4	SetSize;
	BooLean	IsAARes(char *mol);
	BooLean	IsBackBoneAtom(char *atm);
	class obj_typ {	// name of an object for display...
public:
		obj_typ(char *str){ name=AllocString(str); nxt=0; }
		~obj_typ(){ if(name) free(name); if(nxt) delete nxt; }
		void Put(FILE *fp, const char *sep=" "){
			fprintf(fp,"%s",name);
			if(nxt != 0){ fprintf(fp,"%s",sep); nxt->Put(fp); }
		}
		obj_typ *nxt;
private:
		char *name;
	};
	obj_typ	**head;	// head of linked lists.
	obj_typ	**last;	// tail of linked lists.
	
//======= NEW ================
	a_type	AB;
	char	*pdbfile;
	BooLean	CIS_TRANS_PRO;
	char	AARes[21][4];
	void    Free();
	void    FixAtom(char *atom);
	void    print_bond_only(FILE *fp,char res[4],Int4 i,char c,Int4 width,char color,
                const char *atom1,const char *atom2);
	void    PrintResidue1(FILE *fp, char aa, char c, char color, Int4 i,
        		BooLean sideonly,Int4 Width,Int4 big_spacefill);
};

#endif

