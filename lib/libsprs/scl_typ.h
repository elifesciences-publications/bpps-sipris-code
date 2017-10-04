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

#if !defined (_SCL_TYP_)
#define _SCL_TYP_

#include "stdinc.h"
#include "pdb.h"
#include "dheap.h"
#include "swaln.h"
#include "cmsa.h"
#include "set_typ.h"

	//================== structural cluster heap ======================
class	sch_typ {	// structural cluster heap.
public:
	   sch_typ(){ print_error("sch_typ( ) constructor is disallowed"); }
	   sch_typ(Int4 n,Int4 mnp, Int4 strt,Int4 end, char cls, set_typ stX, res_typ *RA, a_type ab){
	     resEnd=end; size=n; mH=Mheap(size,3); NEWP(SCI,size +5, sci_typ);Clss=cls; 
	     resStart=strt; MaxNumPttrns=mnp; Bonferroni=size;
	     // Warning: cannot make Mheap > size of n=MaxPrint or else will return erroneous values..
	     ResALL=RA; AB=ab; SetX=stX; // don't own these...
	   }
           ~sch_typ( ){ Free(); }

	sch_typ	*GetSets(Int4 &II,double *&pval, Int4 *&II2ii,set_typ *&SetII);
	set_typ	Put(FILE *fptr, FILE *vsifp, FILE *efp, FILE *tfp, char FullCluster=0, char *aStr=0, 
			double cutoff=0.05, BooLean PutAll=FALSE,BooLean TakeUnion=FALSE);
	Int4	Insert(Int4 id, Int4 bstm,Int4 mxstr,Int4 *strng,double p){
		if(mxstr < 1) return 0;
		Int4 i=InsertMheap((keytyp)p,mH);
		if(i==0) return 0;
		if(SCI[i]){ delete SCI[i]; }
		SCI[i]= new sci_typ(id,bstm,mxstr,strng,p);
		return i;
	}
	Int4	PrintVSI(FILE *fptr, Int4 i){ return PrintSite(fptr,i,'v'); }
	Int4	PrintWhiteVSI(FILE *fptr, Int4 i){ return PrintSite(fptr,i,'W'); }
	Int4	PrintSite(FILE *fptr, Int4 ii, char mode='m',FILE *tfp=0);
	void	PutItem(FILE *fp, Int4 i){ if(SCI[i]!=0){ SCI[i]->Put(fp); } }
	Int4	Remove(Int4 i){
		   if(i < 1 || i > size) return 0;
		   Int4 rtn=RmMheap(i,mH);
		   if(rtn == 0) return 0;
		   assert(SCI[i]); delete SCI[i]; SCI[i]=0;
		   return rtn;
		}
	Int4	DeleteMin(Int4 &id, Int4 &bstm,Int4 &maxstr,Int4 *&strng,double &p){
		Int4 i=DelMinMheap(mH); 
		if(i==0) return 0;
		SCI[i]->Rtn(id,bstm,maxstr,strng,p); delete SCI[i]; SCI[i]=0;
		return i;
	}
private:
	Int4	Bonferroni;
	class sci_typ {	// structural cluster information type.
	public:
	   	  sci_typ(){ print_error("sci_typ( ) constructor is disallowed"); }
		  sci_typ (Int4 id, Int4 bstm,Int4 mxstr,Int4 *strng,double p){
			ID=id; pval=p; bestm=bstm; MaxStr=mxstr; String=strng;
		  }
		  ~sci_typ(){ if(String) free(String); }
		sci_typ *Copy(){
			   sci_typ *sci= new sci_typ(ID,bestm,MaxStr,this->CopyStrng(),pval);
			   return sci;
			}
		Int4	Rtn(Int4 &id, Int4 &bstm,Int4 &maxstr,Int4 *&strng,double &p){
		   		id=ID; bstm=bestm; maxstr=MaxStr; strng=String; p=pval;
		   		String=0; ID=bestm=MaxStr=0; pval=0;
		}
		void	Put(FILE *fp){
			   fprintf(fp,"bestm=%d; MaxStr=%d; ID=%d; pval=%g\n",bestm,MaxStr,ID,pval);
			}
		
		double	pvalue(){ return pval; }
		double	pval;
		Int4	bestm,MaxStr,ID,*String; 
	private:
		Int4	*CopyStrng(){
			   Int4 i,*s; NEW(s,MaxStr+9,Int4);
			   for(i=0; i <= MaxStr; i++) s[i]=String[i]; 
			   return s;
			}
	};
	sch_typ	*Copy(){
		  sch_typ *sch = new sch_typ(size,MaxNumPttrns,resStart,resEnd,Clss,CopySet(SetX),ResALL,AB);
		  for(Int4 i=1; i <= size; i++){
		     if(SCI[i]){ Int4 x=sch->Insert(SCI[i]->Copy()); assert(x != 0); }
		  } return sch;
		}
	Int4	Insert(sci_typ *sci){
		    Int4 i=InsertMheap((keytyp)sci->pvalue(),mH);
		    if(i==0) return 0;
		    if(SCI[i]){ delete SCI[i]; }
		    SCI[i]= sci; return i;
		}
	sci_typ **SCI;
	void	Free( ){
	  Int4 i;
	  // AB=0;
	  while((i=DelMinMheap(mH))){
	    Int4 id; Int4 bstm,maxstr,*strng; double p;
	    if(SCI[i]){ 
		   SCI[i]->Rtn(id, bstm,maxstr,strng,p);
		   if(strng) free(strng); delete SCI[i]; 
	    }
	  } free(SCI); NilMheap(mH);
	}
	mh_type	mH;
	res_typ	*ResALL;
	set_typ	SetX;
	a_type	AB;
	Int4	size,MaxNumPttrns,resStart,resEnd,Clss; // Clss = 'M', 'R', etc.
}; //================== structural cluster heap ======================

class scl_typ {         // protein structure clustering type
		friend class scm_typ;
public:
                scl_typ(){ time1=0; infile=0; program=0; PutPvals=FALSE; } // called by scm_typ...
                scl_typ(int argc,char *argv[],char mode=' '){
			// fprintf(stderr,"DEBUG SCL: '%c'\n",mode);
			if(mode=='D') InitDCA(argc,argv); else Init(argc,argv);
		}
                ~scl_typ( ){ Free(); }
	Int4	Run();
	Int4	RunDCA();
private:
	char	Mode;		// type of clustering ('S' 'C' 'H' 'B'...)
	Int4	OverlapOfSets();
	BooLean	FindOverlap;
	BooLean	PutAllSCH;
// protected:
	Int4	Driver();	// for scm_typ...
	void	Free();
	void	Init(int argc,char *argv[]){
		   if(strcmp(argv[0],"soprise") == 0 || strcmp(argv[0],"Soprise") == 0) InitPrivate(argc,argv);
		   else InitPublic(argc,argv);
		}
	void	InitDefaults();
	void	InitPublic(int argc,char *argv[]);
	void	InitPrivate(int argc,char *argv[]);
	void	InitDCA(int argc,char *argv[]);
	Int4    ParseResSites(char *str, Int4 *values, const char *msg);
	Int4    ReadPttrnRes(FILE *fp, Int4 &Start, Int4 &End, set_typ &RtnSet, Int4 &nAln,
                char **&pdb_file, char *&chain, char *&clss, Int4 **&PttrnRes);
	Int4	ChainExistsPDB(char chain,pdb_typ pdb);	// returns chain id if it exists.
	UInt4   seed;
	time_t	time1;
	char	*program,*infile,PutTheRest;
	BooLean	PutPvals,UnionCons;
	BooLean OverlappingSeqs(Int4 &offset, Int4 MinOverlap, Int4 MaxMisMatch, e_type E1, e_type E2);
	void	PrintError(char *);
private:	// accessed by scm_typ.
	BooLean	shuffle,BestOnly,dca_vs_pdb;	// if TRUE then shuffle, else permute.
	char	*dca_file,Ignore,FindHbonds,*ChainI,FullCluster,OnlyOne;
	e_type	dcaE;
	Int4	dcaOS;	// offset between dca and pdb sequence files.
	BooLean	FindKeyResPDB(char, Int4, pdb_typ);
	Int4	KeyRes;		// select a key (e.g., active site) residue as a starting point.
	res_typ	KeyResX;	// KeyRes within chain of interest.
	char	KeyChain;
	Int4	KeyAtmNum;
	BooLean	ShowKey;
	float	min_buried;
	Int4    MaxPrint,MaxNumPttrns,method,MaxMinAdj,EdgeDrop,Simulate,BackboneHbonds;
        double  cutoff;
	BooLean	vsi_file,UseHydrogens,UseJeffreys,AroPiBonds,OtherPiBonds;
	h_type  xHG,zHG,xxHG;
	FILE    *efptr,*vsifp;
	FILE	*tab_fp;
	a_type	AB;
};

class scm_typ: private scl_typ {	// structural clustering method type
   public:
                scm_typ(){ print_error("scm_typ( ) constructor is disallowed"); }
                scm_typ(Int4 rS, Int4 rE, set_typ rSt, Int4 na,
			   char ch, char *rcl, Int4 **PR, pdb_typ p,scl_typ *S)
		{
        		resStart=rS; resEnd=rE; nAln=na; PttrnRes=PR;
			chain=ch; res_class=rcl; P=p; resSet=rSt; Init(S); 

		}
                ~scm_typ( ){ Free2(); }
	BooLean	Run(FILE *fptr,FILE *vsifptr,FILE *tfp=0);
	Int4	RtnNumResC(){ return num_resC; }
	// float	**RtnDistMtrxPDB(){ return CalcDistMtrx(); }
	long double	DCAvsPDB(FILE *fp, Int4 NumRed, float **MtrxPDB, float **MtrxDCA, BooLean **Used);
	float   **CalcDistMtrx();
	float   **GetMtrxDCA(char *dca_file, BooLean **&Used);
   private:
        void	Init(scl_typ *S);
	Int4    pvcalc(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jf=FALSE,BooLean cf=TRUE,BooLean pf=TRUE){
		  return pvcalcF(L,D,pos,pval,jf,cf,pf);
		}
	Int4    pvcalcD(Int4 L,Int4 D,Int4 *pos,long double &pval, BooLean jf=FALSE,BooLean cf=TRUE,BooLean pf=TRUE);
	Int4    pvcalcF(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jf=FALSE,BooLean cf=TRUE,BooLean pf=TRUE);
	Int4    pvcalc1(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jeff=FALSE,BooLean cflag=TRUE);
	Int4    pvcalc3(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jeff=FALSE,BooLean cflag=TRUE);
	Int4    pvcalc2(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jeff=FALSE);
	Int4    pvcalc0(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jeff=FALSE);
	float   **GetHbondDistMtrx();
	float   *CalcBuried( );
	Int4	Shuffle();
	Int4	Permute();
	Int4	FindIndex(Int4 s){
		   for(Int4 i=1; i <= num_resC; i++){ if(ResidueID(ResALL[i]) == s) return i; }
		   return 0;
		}
	//========= Methods =============
	double  FixedInterface(char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D);
	double  SphericalClust(Int4 II,Int4 SR, char *hits,Int4 *Str, Int4 &L, Int4 &bestm, Int4 &D);
	double  CoreClust(Int4 II,Int4 SR, char *hits,Int4 *Str, Int4 &L, Int4 &bestm, Int4 &D);
	double  NetworkClust(Int4 II,Int4 SR, char *hits,Int4 *Str, Int4 &L, Int4 &bestm, Int4 &D);
	double  SurfaceClust(Int4 II,Int4 jj, char *hits,Int4 *Str, Int4 &L, Int4 &bestm, Int4 &D);
	double  OptInterfaceClust(char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D);

	//=====================================
	Int4	C,resStart,resEnd,nAln,**PttrnRes,num_resC;
	char	chain,*res_class;
	pdb_typ P;
	set_typ resSet,SetX;	// SetX = discriminating residue set.
	set_typ	SetI;		// SetI = intersecting discriminating residue set.
	res_typ	*ResALL;
	atm_typ	KeyAtm;		// Virtual atom from which ordering residues...
	float	*KeyAtmDist,**DistMtrx,*Buried;
	void	Free2(){ 	// free up in calling environment; don't free twice!!!
		   ChainI=0; program=0; infile=dca_file=0; xHG=zHG=xxHG=0; dcaE=0;
		   for(Int4 i=1; i <= num_resC; i++){
			NilRes(ResALL[i]); if(DistMtrx) free(DistMtrx[i]); 
		   } free(ResALL); if(DistMtrx) free(DistMtrx);
		   if(Buried) free(Buried);
		   if(SetX) NilSet(SetX);
		   if(SetI) NilSet(SetI);
		   if(KeyAtm) NilAtom(KeyAtm);
		   if(KeyAtmDist) free(KeyAtmDist);
		   KeyResX=0; // don't want inherited scl_typ to free this...
		   AB=0;	// not needed since not calling scl_typ::Free()???
		}
};

#endif

