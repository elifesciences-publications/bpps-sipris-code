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

#if !defined (_RAI_TYP_)
#define _RAI_TYP_

#include "histogram.h"
#include "dheap.h"
#include "random.h"
#include "sequence.h"
#include "mheap.h"

#include "residues.h"
#include "atom.h"

#if 0
class rai_typ {         // residue atomic interaction type.
public:
		rai_typ(){ print_error("rai_typ: this constructor not allowed"); }
		rai_typ(Int4 argc, char *argv[]){ Init(argc,argv); }
		~rai_typ( ){ Free(); }
	Int4	Run( );
	Int4	Put(FILE *fptr);
private:
	void	Init(Int4 argc, char *argv[]);
	void	initialize( );
	void	Free( );
	char	type;
	rai_typ	*next;
	a_type  AB;
};
#endif

class rhb_typ {		// residue hydrogen bond type.
public:
		rhb_typ(){ print_error("rhb_typ: this constructor not allowed"); }
		rhb_typ(atm_typ hydrogen, atm_typ donor, atm_typ accept,
			double D_DA, double D_HA, double Angle_DHA){
			Init(hydrogen, donor, accept, D_DA, D_HA, Angle_DHA);
	 	}
		~rhb_typ( ){ Free(); }
	Int4	PutAll(FILE *fp,char color, unsigned short file_id, a_type AB,a_type nAB)
		{
			Put(fp,color, file_id, AB,nAB);
			if(next) next->PutAll(fp,color, file_id, AB,nAB);
		}
	void    Put(FILE *fp,char color, unsigned short file_id, a_type AB,a_type nAB);
	rhb_typ	*Next(){ return next; }
	void	Link(rhb_typ *rhb){ this->next = rhb; }
	rhb_typ	*next;
private:
	void	Init(atm_typ hydrogen, atm_typ donor, atm_typ accept,
                        double D_DA, double D_HA, double Angle_DHA);
	atm_typ	H,Donor,Accept;
	double	d_HA,d_DA,angle_DHA;
	void	Free( ){ if(next) delete next; }
};

/*************************************************
K125.M
 K125A_nz-2hz.X // to  OD1 ASP92A (DA=3.06 A; HA=2.19; DHA=143.7 degrees).
 K125A_nz-1hz.X // to  O   GLY21A (DA=2.86 A; HA=1.90; DHA=161.7 degrees).
 G21A_c-o.X     // to  NZ  LYS125A (DA=2.86 A; HA=1.90; DHA=161.7 degrees).
 K125A_ce-1he.X // to  N9  GSP1174B (DA=3.67 A; HA=2.58; DHA=170.8 degrees).
 10!GSP1174B_n9.X // to  CE  LYS125A (DA=3.67 A; HA=2.58; DHA=170.8 degrees).
/*************************************************/

#endif

