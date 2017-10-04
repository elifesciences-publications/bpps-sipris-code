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

#if !defined (DMS_TYP)
#define DMS_TYP

#include "stdinc.h"
#include "alphabet.h"
#include "residues.h"
#include <math.h>

class dms_typ {
public:
	dms_typ(char wf, double pn, char mode, a_type ab)
		{ wt_factor=wf; pernats=pn; init(mode,ab); }
	// void	score(double count[],double score[]);
	void	score(UInt4 *kount,double *score);
	double	*BackGrnd(){ return back; }
	double  bild(UInt4 *kount);
	double  sdotfs(UInt4 *kount1,UInt4 *kount2);
	double  ComputeLogLike(double *obs,FILE *efp=0);
	double	Typical(UInt4 *kount);
	~dms_typ(){ Free(); }
private:
	a_type	AB;
	double	*RtnData1(double *data,a_type xAB);
	char	wt_factor;
	double  pernats;
	double	*count,*countX;
	void	init(char mode,a_type ab);
	void	Free( );
	Int4	nDC;
	double	*offset;	// used to speed up the calculation...
	double  *weight;             /*  Weight of Dirichlet component               */
	double  *Alpha;              /*  Concentration parameters                    */
	double  **alpha;          /*  Dirichlet parameters                        */
	double  *back;               /*  Implicit background amino acid frequencies  */

}; 

#endif

