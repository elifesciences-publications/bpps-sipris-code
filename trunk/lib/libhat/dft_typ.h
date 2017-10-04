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

#if !defined (_DFT_TYP_)
#define _DFT_TYP_
#include "stdinc.h"
#include "afnio.h"
// Header file for a depth-first traversal (integer Level) representation of a tree.

class dft_typ {
public:
		dft_typ( ){ assert(!"Illegal constructor"); }
		dft_typ(Int4 N, Int4 *Array){ Init(N,Array); }
		dft_typ(char *file){
			FILE *fp = open_file(file,".dft","r");
			BooLean okay=Read(fp); fclose(fp);
			// if(!okay) print_error("depth-first traversal tree file syntax error");
			if(!okay) {
				fprintf(stderr,"WARNING: depth-first traversal tree file syntax error.\n");
				fprintf(stderr,"   Using simple two-level pruning tree instead.\n");
			}
		}
		~dft_typ( ){ Free(); }
	void	Write(char *file){
			FILE *fp = open_file(file,".dft","w");
			Write(fp); fclose(fp);
		}
	void	Write(FILE *);
	Int4	*RtnLevels( ) { return Level; }
	Int4	NumNodes( ){  return num_nodes; }
	BooLean	IsValid( ){  return Valid; }
private:
	BooLean Read(FILE *);
	Int4	num_nodes;		// number of nodes in tree.
	Int4	*Level;			// Level for levels of each node
	void	Free();
	BooLean	ValidityCheck();
	BooLean	Valid;
	void	Init(Int4 N, Int4 *Array);
};

#endif
