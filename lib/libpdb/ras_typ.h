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

/* ras_typ.h - rasmol type output program (for creating rasmol scripts). */
#if !defined (_RAS_TYP_)
#define _RAS_TYP_

#include <time.h>
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"

class ras_typ {         // Rasmol type
public:
	ras_typ( ){ assert(!"Illegal constructor"); }
	ras_typ(char *filename,BooLean cis_trans_pro);
	~ras_typ( ){ Free( ); }

	void	PrintHEADER(FILE *fp,BooLean use_path, char c,Int4 wire_width);
	void    PrintTAIL(FILE *fp);
	void    PrintTrueColor(FILE *fp, char c);
	void    PrintColor(FILE *fp, char c);
	void    ColorAtom(FILE *fp, char *mol, char c, char color, Int4 i, char *atom);
	void    PrintClear(FILE *fp);
	void	PrintView(FILE *fp, Int4 rx,Int4 ry,Int4 rz,Int4 tx,Int4 ty,Int4 tz);
	char	*DefineResCloud(FILE *fp, char aa, Int4 site, char chain);
	char    *DefineResCloud(FILE *fp, char aa, Int4 site, char *atom1,char *atom2,
					char chain);
	char    *DefineMolCloud(FILE *fp, char *mol, Int4 site, char *atom1,
			char *atom2, char chain);
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

	void    PrintResidue(FILE *fp, char aa, char c, char color, Int4 i,
        		BooLean sideonly,Int4 Width, Int4 big_spacefill);
	void    PrintMolecule(FILE *fp, char *mol, char c, char color, Int4 i,
			Int4 Width,Int4 spacefill);
	void    Init();
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

