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

#if !defined (_SAS_TYP_)
#define _SAS_TYP_

#include "psc_typ.h"
#include "swaln.h"

class sas_typ {         // structural alignment-based scoring type
public:
		sas_typ(){ print_error("sas_typ( ) constructor is disallowed"); }
		sas_typ(Int4 argc, char *argv[],set_typ **ColPairSet=0){
			TraceColor=SideColor=0;
			seed = 18364592; Target=5; MaxDataPoints=4000; KeySeq=-1;
			time1=time(NULL); StampString();
			Init( ); GetArg(argc, argv); ReadInput( );
        		if(seed == 18364592)  seed=(UInt4) time(NULL)/2;
			Int4 MinSeqOverlap = (Int4) ceil(LengthCMSA(1,IN_CMA[1])* 0.75);
   			if(MinSeqOverlap < 10) MinSeqOverlap=10;
   			if(MinSeqOverlap > 100) MinSeqOverlap=100;
			esc=new esc_typ(mpdb,MinSeqOverlap);
			p2c=new p2c_typ(0,0,mpdb,IN_CMA,Number,mcBPPS_prefix,esc,hpt,ptrn,
						TraceColor,SideColor,ColPairSet);
			p2c->SetTarget(Target);
			p2c->SetDataPoints(MaxDataPoints);
			p2c->FinalPnts=0; p2c->FinalCols=0; 
		}
		~sas_typ( ){
			Free(); 
			double runtime=difftime(time(NULL),time1);
			fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
		}
	double	Run();
	// void	AssignSets(set_typ **sets){ p2c->AssignSets(sets); }
	Int4	RtnFinalPts() { return p2c->FinalPnts; }
	Int4	RtnFinalCols() { return p2c->FinalCols; }
	void	SetDataPoints(Int4 X)
		   { assert(X > 0); MaxDataPoints=X; p2c->SetDataPoints(MaxDataPoints); }
	e_type	RtnCsq(){ return TrueSeqCMSA(1,IN_CMA[2]); }
	cma_typ	RtnCMSA(){ return IN_CMA[1]; }
	a_type	RtnAlphabet(){ return AB; }
	adv_typ **RunAndRtnADV(FILE *fp,Int4 &N){ return p2c->RunAndRtnADV(fp,N); }
	double	*ColumnARSD_Scores(FILE *fp,Int4 *&N);
	double	*SeqARSD_Scores(FILE *fp,Int4 *&N);
	void	PutAdvList(FILE *fp, Int4 N, adv_typ **ADV,Int4 shown=0,char *message=0)
			{ p2c->PutAdvList(fp,N,ADV,shown,message); }
	char	Mode;
	void	PutID(FILE *fp,Int4 i){ 
			char *p,*s=AllocString(p2c->RtnPDB_ID(i));
			p=strstr(s,"_H.pdb"); assert(p > 0); p[0]=0; p -= 4; 
			fprintf(fp,"%s",p);  free(s);
		}
	Int4	KeySeq;
private:
	Int4	Target;		// target number of seqs to use in histogram.
	Int4	MaxDataPoints;
	char	*TraceColor,*SideColor;
	time_t  time1;
	void	Validate();
	void	Init( );
	void    ReadInput( );
	void	GetArg(Int4 argc, char *argv[]);
	void	Free();
        UInt4   seed;
	char	*pdb_paths,*mcBPPS_prefix;

	char	string[300];
	void	StampString() { string[298] = 0; }
	void	CheckString() { if(string[298] != 0) print_error("input argument too long"); }

	// Hyperpartition input files
	a_type	AB;
	hpt_typ	*hpt;		// hyperpartition.	NOT NEEDED.
	pat_typ *ptrn;		// patterns...		NOT NEEDED.
	mps_typ *mpdb;		// multiple protein structure type
	esc_typ	*esc;		// equivalent (protein) structural chains type.
	cma_typ	*IN_CMA;
	p2c_typ *p2c;		// pdb to cma mapping type.
	Int4	Number;		// number of cma input alignments.
};

#endif

