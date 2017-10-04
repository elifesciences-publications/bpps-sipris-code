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

#if !defined(_HAT_TYP_H_)
#define _HAT_TYP_H_
#include "alphabet.h"
#include "sset.h"
#include "probability.h"
#include "histogram.h"
#include "dheap.h"

#include "gpsi_typ.h"
#include "patch_sap.h"
#include "residues.h"
#include "smatrix.h"
#include "cmsa.h"
#include "gpsi2cma.h"

#include "sph_typ.h"
#include "dsets.h"


class hat_typ {		// Hierarchical alignment template type
public:
		hat_typ(FILE *fp);
		~hat_typ(){ Free(); }
	void	Put(FILE *fp);
	void    hat_srch(FILE *fp);
private:
	void	Init(FILE *fp);
	void	Free();

};

int     cma_gblastpgp(int argc, char *argv[]);
int     curated_srch(int argc, char *argv[]);
int     gapmap_srch(int argc, char *argv[]);
int     setup_gapmaps(int argc, char *argv[]);
int     public_setup_mapgaps(int argc, char *argv[]);
int     public_setup_gapmaps(int argc, char *argv[]);
int     main_setup_gapmaps(int argc, char *argv[],const char *usage,FILE *fp=0);
int     matblast_srch(int argc, char *argv[]);
int     matblast_srch(int argc, char *argv[],const char *usage);
void    TrimToTemplateCMSA(cma_typ *IN_CMA, Int4 num_cma_files);
cma_typ *ConversionViaTemplateCMSA(cma_typ TemplateCMA, cma_typ *IN_CMA);

// New version (debugged) afn: 7/18/2013.
char    *AddInsertToOperationArray3(Int4 start, Int4 end, char *operation);
cma_typ *ConversionViaTemplateCMSA3(cma_typ TemplateCMA, cma_typ *IN_CMA);

cma_typ GrowEdgesCMSA(BooLean add2Nterm, Int4 NumGaps, cma_typ cma);
// cma_typ SplitBlkCMSA3(Int4 x, Int4 left_leng, Int4 minlen, cma_typ L);
// cma_typ InsertColumnsCMA3(Int4 start, Int4 N, cma_typ cma);

BooLean IsOkayTemplateCMA(cma_typ TemplateCMA);
void    PutVSIregionsCMSA(FILE *fp,Int4 sq,char *color,Int4 *start,Int4 *end,
                Int4 num_regions,cma_typ cma);
Int4    ParseInputRegionsVSI(char *input_string,Int4 **Start,Int4 **End,
                char **Colors,const char *usage);
int     convert_template_main(Int4 argc,char *argv[]);

#endif
