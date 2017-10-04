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

#if !defined (_HSC_TYP_)
#define _HSC_TYP_

#include "psc_typ.h"

class hsc_typ {         // hyperpartition structural comparison type
public:
		hsc_typ(){ print_error("hsc_typ( ) constructor is disallowed"); }
		hsc_typ(Int4 argc, char *argv[], FILE *mmafp=0, FILE *hptfp=0, FILE *ptrnfp=0){
			TraceColor=SideColor=0;
			seed = 18364592;
			time1=time(NULL); StampString();
			GetArg(argc, argv);
        		if(seed == 18364592)  seed=(UInt4) time(NULL)/2;
			Init(mmafp,hptfp,ptrnfp); 
			Int4 MinSeqOverlap = (Int4) ceil(LengthCMSA(1,IN_CMA[1])* 0.75);
   			if(MinSeqOverlap < 10) MinSeqOverlap=10;
			esc=new esc_typ(mpdb,MinSeqOverlap);
			p2c=new p2c_typ(0,0,mpdb,IN_CMA,Number,mcBPPS_prefix,esc,hpt,ptrn,
						TraceColor,SideColor);
		}
		~hsc_typ( ){
			Free(); 
			double runtime=difftime(time(NULL),time1);
			fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
		}
	int	Run(char mode=0,char call=0,FILE *logfp=0,FILE *vsifp=0);
	pat_typ *SubPatterns(Int4 Row);
	void    VerboseOff(){ verbose=FALSE; }
        void    VerboseOn(){ verbose=TRUE; }
private:
	BooLean	verbose;
	char	*TraceColor,*SideColor;
	time_t  time1;
	void	Validate();
	void	Init(FILE *mmafp=0, FILE *hptfp=0,FILE *ptrnfp=0);
	void	GetArg(Int4 argc, char *argv[]);
	void	Free();
        UInt4   seed;
	char	*pdb_paths,*mcBPPS_prefix;

	char	string[300];
	void	StampString() { string[298] = 0; }
	void	CheckString() { if(string[298] != 0) print_error("input argument too long"); }

	// Hyperpartition input files
	a_type	AB;
	hpt_typ	*hpt;
	pat_typ *ptrn;
	mps_typ *mpdb;
	esc_typ	*esc;
	cma_typ	*IN_CMA;
	p2c_typ *p2c;
	Int4	Number;		// number of cma input alignments.
};

#endif

