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

#include "scl_typ.h"

Int4     scl_typ::ParseResSites(char *str, Int4 *values, const char *msg)
// input string: "3,5-7,9,11-17"
// returns:      12 and sets values=[12,3,5,6,7,9,11,12,13,14,15,16,17]
//                                    0 1 2 3 ...
{
        Int4     n,v,w;

        if(!isdigit(str[0])) print_error(msg);
        for(n=1; str[0] != 0 && str[0] != '\n'; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%d", &v) != 1) print_error(msg); 
                else { // fprintf(stderr,"v=%d\n",v);
                   values[n] = v; while(isdigit(str[0])) str++;
                }
           } else print_error(msg);
        } return n;
}


Int4	scl_typ::ReadPttrnRes(FILE *fp, Int4 &Start, Int4 &End, set_typ &RtnSet, Int4 &nAln,
		char **&pdb_file, char *&chain, char *&clss, Int4 **&PttrnRes)
/********************************************************************
file: /tmp/pdb_sema/1olz_H.pdb
chain: A.
file: /tmp/pdb_sema/1olz_H.pdb
chain: B.
range: 28-351,356-461,464-480(447;6).
M=105,260,305,114,292,236,190,211,234,307,253,106,107,239,269,89,44,252,425,153,303,302,255,258,216
R=372,379,376,349,344,367,342,346,334,366,132,374,345,347,81,318,80,101,66,446,320,179,207,378,209
O=180,299,111,369,298,301,221,325,130,47,245,159,398,164,363,219,154,270,480,161,170,184,328,122,84

file: /tmp/pdb_sema/3ol2_H.pdb
chain: A.
range: 49-372,377-482,485-501(447;6).
M=126,281,326,135,313,257,211,232,255,328,274,127,128,260,290,110,65,273,446,174,324,323,276,279,237
R=393,400,397,370,365,388,363,367,355,387,153,395,366,368,102,339,101,122,87,467,341,200,228,399,230
O=201,320,132,390,319,322,242,346,151,68,266,180,419,185,384,240,175,291,501,182,191,205,349,143,105

 ********************************************************************/
{
	Int4	i,j,n,N,*Ptr[300],Values[2000],Ins;
	char	c,str[509],name[100],Clss[300],state='s',chn[300],*Names[300];
	for(j=N=0,i=1; (fgets(str,500,fp) != NULL); i++){
	   if(str[0]=='#') continue;
	   if(!isalpha(str[0])) break;
// fprintf(stderr,"%s\n",str);
	   switch(str[0]){
	     case 'f': {
		if(state != 'c' && state != 's') print_error("SIPRIS syntax error 1");
		if(sscanf(str,"file: %s\n",name) != 1){
		   print_error("SIPRIS pdb input file name error");
		} N++; Names[N]=AllocString(name); state='f'; 
	       } break;
	     case 'c': {
		if(state != 'f') print_error("SIPRIS syntax error 2");
		if(sscanf(str,"chain: %c.\n",&c) != 1 || !isupper(c)){
		   while((fgets(str,500,fp) != NULL) && isalpha(str[0])) ; 
		   fprintf(stderr,"  (SIPRIS syntax error: ");
		   fprintf(stderr," skipping %s chain %c)\n",Names[N],c);
		   return -1;
		   print_error("SIPRIS chain input error");
		} chn[N]=c; state='c';
	       } break;
	     case 'r': {
		if(state != 'c') print_error("SIPRIS syntax error 3");
		char str2[500];
		if(sscanf(str,"range: %[^(](%d;%d).",str2,&nAln,&Ins) != 3){
		   print_error("SIPRIS range formatting error");
		} RtnSet=ParseSet(str2,Start,End); 
#if 0
	Int4	low, high;
	char *str0=RtnStrSet(RtnSet,low,high);
	fprintf(stderr,"str in=%s\n range: %d..%d.\n",str2,Start,End);
	fprintf(stderr,"CardSet=%d\n",nAln); 
	fprintf(stderr,"str out=%s\n range: %d..%d\n",str0,low,high); free(str0);
	fprintf(stderr,"CardSet=%d\n",CardSet(RtnSet)); 
#endif
		if(Start > End) print_error("SIPRIS range input error");
		else {
		   n=End-Start+1;
		   // fprintf(stderr,"n=%d; nAln=%d; Ins=%d; cardSet=%d\n",n,nAln,Ins,CardSet(RtnSet));
		   if(n != (nAln+Ins)) print_error("SIPRIS inconsistent range input");
		   if(CardSet(RtnSet) != nAln) print_error("SIPRIS residue set inconsistent input");
		} state = 'r';
	       } break;
	     default: {
	        if(!isupper(str[0]) || !(state == 'r' || state == 'i')){
			print_error("SIPRIS syntax error 4");;
		}
		Int4 *values; j++; Clss[j]=str[0]; 
		char *Str = &str[2];
// fprintf(stderr,"%s\n",Str);
		n=ParseResSites(str +2,Values,"SIPRIS_pdb input error");
	  	n= MINIMUM(Int4,n,MaxNumPttrns);  // only use up to MaxNumPttrns.
		NEW(values,n+5,Int4); values[0]=n;
		for(Int4 x=1; x <=n; x++) values[x]=Values[x];
		Ptr[j]=values; state='i';
	       } break;
	   }
	}
	if(N > 0){
	  if(state != 'i'){
	   fprintf(stderr,"state = %c; N=%d; j=%d; i=%d; name=%s\n",state,N,j,i,name);
	   fprintf(stderr,"str = %s\n",str);
	   print_error("SIPRIS syntax error 5");
	  }
	  NEWP(pdb_file,N+5,char); NEW(chain,N+5,char);
	  for(i=1; i <= N; i++){ pdb_file[i]=Names[i]; chain[i]=chn[i]; }
	}
	if(j > 0){
	  NEWP(PttrnRes,j+5,Int4); NEW(clss,j+5,char);
	  for(i=1; i <= j; i++){ PttrnRes[i]=Ptr[i]; clss[i]=Clss[i]; }
	} return j;
}

Int4	scl_typ::ChainExistsPDB(char chain,pdb_typ pdb)
{
	for(Int4 c=1; c <= nChainsPDB(pdb); c++) if(ChainCharPDB(c,pdb)==chain) return c;
	return 0;
}

Int4    scl_typ::OverlapOfSets()
// Find the overlap between category residue sets...
{
	FILE	*infp=open_file(infile,"","r");
	Int4	i,j,n,resStart=0,resEnd=0,**PttrnSites,nAln; 
	char	*chain=0,**pdb_file=0,*res_class;
	set_typ RtnSet=0,SetX[50];
	Int4	OL[50][50],Union[50][50];
	do {
	  n=ReadPttrnRes(infp,resStart,resEnd,RtnSet,nAln,pdb_file,chain,res_class,PttrnSites);
	  if(n > 0){
	   assert(n < 50);
	   Int4 N=CardSet(RtnSet);  
	   ClearSet(RtnSet);
	   for(i=1; i <= n; i++){ 
		SetX[i]=CopySet(RtnSet);
		for(j=1; PttrnSites[i][j]; j++){
		   AddSet(PttrnSites[i][j],SetX[i]);
		   // fprintf(stderr,"%d.%d: %d\n",i,j,PttrnSites[i][j]); 
		} free(PttrnSites[i]);
	   }
	   for(i=1; i < n; i++){ 
	      Int4 red,black,out,out_red,u;
	      red = CardSet(SetX[i]); black = N - red;
	      for(j=i+1; j <= n; j++){ 
		out = CardSet(SetX[j]); out_red = CardInterSet(SetX[i],SetX[j]);
		u = CardUnionSet(SetX[i],SetX[j]);
		double prob = CumHypGeomProb(red,black,out,out_red);
		fprintf(stdout,"%c(%d) & %c(%d) = %d; U = %d; N = %d; prob=%.3g\n",
		   res_class[i], red, res_class[j], out, out_red,u,N,prob);
		set_typ TmpSet=CopySet(SetX[i]); IntersectSet1(SetX[i],SetX[j],TmpSet);
		PutSet(stdout,TmpSet); NilSet(TmpSet); 
		OL[i][j]=OL[j][i]=out_red; Union[i][j]=Union[j][i]=u;
	      } OL[i][i]=Union[i][i]=red;
	   } OL[n][n]=Union[n][n]=CardSet(SetX[n]); 
	   fprintf(stdout,"\n     ");
	   for(i=1; i <= n; i++) fprintf(stdout,"%4c ",res_class[i]); fprintf(stdout,"\n");
	   for(i=1; i <= n; i++){ 
	      fprintf(stdout,"%4c ",res_class[i]);
	      for(j=1; j <= n; j++) fprintf(stdout,"%4d ",OL[i][j]);
	      fprintf(stdout,"\n");
	   }fprintf(stdout,"\n\n");
#if 0
	   fprintf(stdout,"\n     ");
	   for(i=1; i <= n; i++) fprintf(stdout,"%4c ",res_class[i]); fprintf(stdout,"\n");
	   for(i=1; i <= n; i++){ 
	      fprintf(stdout,"%4c ",res_class[i]);
	      for(j=1; j <= n; j++) fprintf(stdout,"%4d ",Union[i][j]);
	      fprintf(stdout,"\n");
	   }fprintf(stdout,"\n\n");
#endif
	   for(i=1; i <= n; i++) NilSet(SetX[i]);
	   for(j=1; pdb_file[j] != 0; j++) free(pdb_file[j]);
	   free(chain); free(pdb_file); free(PttrnSites); free(res_class); NilSet(RtnSet);
	  }
	} while(n); fclose(infp);
	return 0;
}

Int4	scl_typ::Run( )
{ 
	if(FindOverlap){ OverlapOfSets(); return 0; }
	Int4	nruns=1;
	if(!dca_file && (Simulate > 0 || efptr || PutPvals)){
	   nruns=MAXIMUM(Int4,1,Simulate);
	   xHG=Histogram("P-values",0,1,0.025);
	   zHG=Histogram("bestm values",0,500,1.0); 
	   xxHG=Histogram("P-values",0,1,0.001); 
	}
	for(Int4 i=1; i <= nruns; i++) Driver();
	if(zHG){ PutHist(stdout,60,zHG); NilHist(zHG); zHG=0; }
	if(xHG){ 
	   Int4 total=TotalDataHist(xHG); 
	   double d=1.0/(double)total; PutHist(stdout,60,xHG); NilHist(xHG); xHG=0; 
	   d = -expm1(-d); fprintf(stdout,"  expected P-value = %g (N=%d)\n\n",d,total);
	} if(xxHG){ PutHist(stdout,60,xxHG); NilHist(xxHG); xxHG=0; }
}

Int4	scl_typ::Driver()
{
	Int4	i,j,k,n,x;
	Int4	resStart=0,resEnd=0,**PttrnSites,nAln,II=0; 
	char	c=' ',*chain=0,**pdb_file=0,*res_class;
	FILE	*infp=open_file(infile,"","r");
	pdb_typ	PDB;
	set_typ RtnSet=0;
	do {
	  n=ReadPttrnRes(infp,resStart,resEnd,RtnSet,nAln,pdb_file,chain,res_class,PttrnSites);
#if 0
	  for(i=1; i <= n; i++){ 
		for(j=1; PttrnSites[i][j]; j++) fprintf(stderr,"%d.%d: %d\n",i,j,PttrnSites[i][j]); 
	  } fprintf(stderr,"n=%d\n",nAln); exit(1);
#endif
	  if(n < 0) continue;
	  if(n > 0){
	    for(j=1; pdb_file[j] != 0; j++){
	      assert(chain[j] != 0);
	      PDB = MakePDB(pdb_file[j]);	// PutPDB(stdout,P);
	      II++;
	      if(ChainI){
		for(k=0; ChainI[k]; k++){
	          if(!ChainExistsPDB(ChainI[k],PDB)){
		    free(pdb_file[j]); NilPDB(PDB); PDB=0; break;
		  }
		} if(PDB==0) continue;
		// Don't include query chain!
	      	if(strchr(ChainI,chain[j]) != NULL){ free(pdb_file[j]); NilPDB(PDB); continue; }
	      }
	      fprintf(stdout,"\n=========== %d. %s:%c ============\n",II,pdb_file[j],chain[j]);
	      if(KeyChain){	// check key residue and key chain...
		if(!FindKeyResPDB(KeyChain,KeyRes,PDB))
		  { free(pdb_file[j]); NilPDB(PDB); continue; }
		if(chain[j]==KeyChain){ 		// Is key residue in this chain?
		  Int4 nn=n;
		  for(k=1; PttrnSites[k]; k++){
		     Int4 *PSJ=PttrnSites[k];
		     for(i=1; PSJ[i]; i++){
		        if(KeyRes == PSJ[i]){
		           fprintf(stderr,"%d %c.%d: %d (%d)\n",k,res_class[k],i,PSJ[i],KeyRes); 
			   fprintf(stderr,"Class %c: ",res_class[k]); 
			   fprintf(stderr,"-K option start point cannot be a distinguishing residue.\n"); 
			   free(PttrnSites[k]); 
			   for(x=k; x <= nn; x++){
			     PttrnSites[x]=PttrnSites[x+1]; res_class[x]=res_class[x+1];
			   } k--; nn--;
			   PttrnSites[nn+1]=0; res_class[nn+1]=0; break;
			}
		     }
		  }
		  if(nn==0 && PttrnSites[1]==0){ free(pdb_file[j]); NilPDB(PDB); continue; }
		}
		if(KeyAtmNum && AtomResidue(KeyAtmNum,KeyResX) == 0)
		  { free(pdb_file[j]); NilPDB(PDB); continue; }
		if(ShowKey){
	           // fprintf(stdout,"\n=========== %d. %s:%c ============\n",II,pdb_file[j],chain[j]);
        	   for(Int4 a=1; a <= ResidueAtomNumber(KeyResX); a++){
		     fprintf(stdout,"%3d: ",a); PutAtom(stdout,AtomResidue(a,KeyResX));
        	   } // PrintResidueAtoms(stdout,KeyResX, FALSE); 
		   free(pdb_file[j]); NilPDB(PDB); continue; 
		}
	      }
	      if(dcaE){	// check for consistency...
		e_type E=GetPDBSeq(GetChainNumberPDB(PDB,chain[j]),PDB);
        	if(!this->OverlappingSeqs(dcaOS,30,0,dcaE,E)){
           	  fprintf(stderr,"Sequences mismatch (%s)\n",pdb_file[j]);
		  // a_type AB=AminoAcidAlphabetPDB(PDB); AlnSeqSW(stderr,11, 1,dcaE, E, AB);
           	  // PutSeq(stderr,dcaE,AB); PutSeq(stderr,E,AB);
		  free(pdb_file[j]); NilPDB(PDB); continue; 
        	} 
		Int4 pdbOS=OffSetSeq(E),os=OffSetSeq(dcaE); dcaOS=os-pdbOS+dcaOS;
		fprintf(stderr,"pdbOS= %d; os=%d; dcaOS=%d\n",pdbOS,os,dcaOS);
		NilSeq(E);
	      }
#if 1	// Compute significance of DCA-structural overlap.
	      scm_typ scm(resStart,resEnd,RtnSet,nAln,chain[j],res_class,PttrnSites,PDB,this);
	      if(dca_file && dca_vs_pdb){
		assert(dcaE);
		long double	pMin=1.0,pmin=1.0,score,pp;
		Int4	xx,Nr,bstNr=0,BstNr=0,start=LenSeq(dcaE),end=15*start;
		Int4	inc=10,Inc=25,Sim=0;
		float    **MtrxPDB=0,**MtrxDCA=0;
		BooLean	**Used=0;
		MtrxPDB=scm.CalcDistMtrx(); MtrxDCA=scm.GetMtrxDCA(dca_file,Used); 
#if 0
		// start=end=MaxNumPttrns;
		bstNr=MaxNumPttrns;
#else
		// long double ddd=lgammal(20); // some C++ libraries give an Uninitialized memory read error.
		fprintf(stderr,"start=%d; end=%d\n",start,end);
		h_type	HG=0;
		if(Simulate > 0){ HG=Histogram("-log(p)",0,20,0.25); Sim=Simulate; }
		else HG=Histogram("-log(p)",0,500,5.0);
		do {
		 for(xx=0,Nr=start; Nr <= end; Nr+= inc){ 
		   score,pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp > 0){ score=-log10(pp); IncdHist(score,HG); }
		   if(pp < pmin){
			fprintf(stderr,"%d. %.1Lf (%d)\n",Nr,score,xx);
			pmin=pp; bstNr=Nr; xx=0; 
		   } else xx++;
		 } BstNr=bstNr;
		 for(Nr=BstNr+Inc; Nr > BstNr; Nr--){ 
		   pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp > 0){ score=-log10(pp); IncdHist(score,HG); }
		   if(pp < pmin){ pmin=pp; bstNr=Nr; }
		 }
		 for(Nr=BstNr-Inc; Nr < BstNr; Nr++){ 
		   pp=scm.DCAvsPDB(0,Nr,MtrxPDB,MtrxDCA,Used); 
		   if(pp > 0){ score=-log10(pp); IncdHist(score,HG); }
		   if(pp < pmin){ pmin=pp; bstNr=Nr; }
		 }
		 if(Sim > 0){ 
		   fprintf(stderr,"======== %d ==========\n",Sim); Sim--; 
		   fprintf(stderr,"pmin=%.3Lg; bstNr=%d\n",pmin,bstNr);
		 }
		} while(Sim > 0);
		PutHist(stderr,60,HG); NilHist(HG);
#endif
		if(Simulate == 0){
	           fprintf(stdout,"\n=========== %d. %s:%c ============\n",II,pdb_file[j],chain[j]);
		   scm.DCAvsPDB(stdout,bstNr,MtrxPDB,MtrxDCA,Used); 
		} else fprintf(stderr,"pmin=%.3Lg; bstNr=%d\n",pmin,bstNr);
		Int4 num_resC=scm.RtnNumResC();
		for(i=1; i <= num_resC; i++){ free(MtrxPDB[i]); } free(MtrxPDB);
		for(i=1; i <= num_resC; i++){ free(MtrxDCA[i]); } free(MtrxDCA);
		for(i=1; i <= num_resC; i++){ free(Used[i]); } free(Used);
		free(pdb_file[j]); NilPDB(PDB); break; 
	      }
#endif
	      if(vsifp){	// then print out the header for the next structure.
#if 0
	        fprintf(vsifp,"~$=%d.\nFile1=%s:%c  // \n\n%d-%d.B80\n",II,pdb_file[j],chain[j],resStart,resEnd);
#else
	        fprintf(vsifp,"~$=%d.\nFile1=%s:%c  // SIPRIS generated file\n\n",II,pdb_file[j],chain[j]);
		if(resStart < 1) print_error("sipris input error: residue start < 1");
		Int4	r,rm1,strt;	// will get a likely error for MemberSet(1+SetN(RtnSet)!!! 
		for(r=resStart,rm1=r-1; r <= resEnd; r++,rm1++){
		   if(MemberSet(r,RtnSet) && !MemberSet(rm1,RtnSet)){ 
			strt=r; r++; rm1++;
			if(r <= SetN(RtnSet)){
			   while(MemberSet(r,RtnSet) && MemberSet(rm1,RtnSet)){ r++; rm1++; if(r > SetN(RtnSet)) break; }
			} fprintf(vsifp,"%d-%d.B80\n",strt,rm1);
		    }
		}  // fprintf(stderr,"N=%d (%d)\n",SetN(RtnSet),MemberSet(0,RtnSet));
#endif
		// "  // " needed for parsing...fix chn_vsi eventually.
	      }
	      // fprintf(stderr,"scm_typ test: %c\n",chain[j]);
	      if(scm.Run(stdout,vsifp,tab_fp)){ if(vsifp) fprintf(vsifp,"\n\n"); }
	      free(pdb_file[j]); NilPDB(PDB);
	      if(Simulate && !shuffle) break;	// do only one simulation...otherwise may repeat!!!
	    }
	    for(i=1; i <= n; i++){ if(PttrnSites[i]) free(PttrnSites[i]); }
	    free(chain); free(pdb_file); free(PttrnSites); free(res_class); NilSet(RtnSet);
	  }
	} while(n); fclose(infp);
	return 0;
}

BooLean	scl_typ::FindKeyResPDB(char KeyChain,Int4 KeyRes,pdb_typ PDB)
{
	Int4	i,j,num_resX;
	if(KeyResX){ NilRes(KeyResX); KeyResX=0; }
	Int4 KeyChn=ChainExistsPDB(KeyChain,PDB);
	if(KeyChn<=0) return FALSE;
	res_typ *ResAllX=MakeResPDB(KeyChn,&num_resX,PDB);
	for(i=1; i <= num_resX; i++){
		if(ResidueID(ResAllX[i]) == KeyRes){
		   KeyResX=ResAllX[i];
		   for(j=i+1; j <= num_resX; j++){ NilRes(ResAllX[j]); }
		   break;
		} else NilRes(ResAllX[i]);
	} free(ResAllX);
	if(KeyResX==0) return FALSE; 
#if 0	// DEBUG...
	fprintf(stderr,"KeyChn=%d; KeyChain=%c; KeyRes=%d.\n",KeyChn,KeyChain,KeyRes);
	atm_typ KeyAtm=MkAveAtomResidue(KeyResX); PutAtom(stderr,KeyAtm);
#endif
	return TRUE;
}

#include "stdinc.h"

extern  int ChainVSI(int argc,char *argv[],FILE *ifp, FILE *ofp);

void	scl_typ::Free()
{
	if(vsifp){
	  if(strcmp(program,"soprise") != 0 && strcmp(program,"Soprise") != 0){
	     Int4	i,Argc,file=1,id=1;
	     char str[200],str0[100],str1[100],*Argv[20];
	     rewind(vsifp);
             sprintf(str0,"%s_X",infile);
#if 0
             sprintf(str1,"chn_vsi %s %d %d -T -skip=A -d2.5 -D",str0,id,file);
#else
             sprintf(str1,"chn_vsi %s %d %d -T -skip=W -d2.5 -D",str0,id,file);
#endif
             Argc=string2argv(Argv,str1);       // mode == 'T'
             // for(i=0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n");
             FILE *ifp,*ofp = tmpfile();
             // ChainVSI(Argc,Argv,0,ofp); rewind(ofp); ifp=ofp;
             ChainVSI(Argc,Argv,vsifp,ofp); rewind(ofp); ifp=ofp;
             for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);

#if 1	
	     // fprintf(stderr,"infile=%s %c\n",infile,Mode); 
             sprintf(str1,"chn_vsi %s.crs %d -d2.5 -c -D -pml=%s",str0,id,infile);
#else
             sprintf(str1,"chn_vsi %s.crs %d -d2.5 -c -D -pml=%s_%c",str0,id,infile,Mode);
#endif
             Argc=string2argv(Argv,str1);       // mode == 'p'
             // for(i = 0; i < Argc; i++) { fprintf(stderr,"%s ",Argv[i]); } fprintf(stderr,"\n");

             // ChainVSI(Argc,Argv,0,0);
             ChainVSI(Argc,Argv,ifp,0); fclose(ifp);
             for(Argc-- ; Argc >= 0; Argc--) free(Argv[Argc]);
	  } fclose(vsifp); 
	}
	if(tab_fp) fclose(tab_fp); 
	if(dcaE) NilSeq(dcaE);  // don't free up in scm_typ!!!
	if(infile) free(infile);
	if(dca_file) free(dca_file);
	if(program) free(program);
	if(ChainI) free(ChainI);
	if(KeyResX) NilRes(KeyResX);
	if(AB) NilAlpha(AB); AB=0;
	if(time1 != 0){
	   double runtime=difftime(time(NULL),time1);
	   fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	}
}

#define WUSAGE "Usage: SIPRIS klst_file [options]\n\
       -B            Search for H-bond networks in each category plus higher categories.\n\
       -Bb           Same as -B option, but also include backbone-to-backbone H-bonds.\n\
       -b            Search for H-bond networks in each category alone.\n\
       -bb           Same as -b option, but also include backbone-to-backbone H-bonds.\n\
       -best         Output best hits only\n\
       -C=<real>     E-value cutoff (default: 0.05).\n\
       -c=<int>      Minimum adjacency cutoff for method 2 (default: 7; range 2-15).\n\
       -d=<int>      Drop point for minimum # of adjacent residues (default: c+1: range 0-15).\n\
       -EV=<file>    Use EVfold CouplingScores.csv to compute structural contacts (for -m=1-4 only).\n\
       -ev=<file>    Compute p-value for overlap between EVfold scores and structural contacts.\n\
       -e            Print additional information out to stderr.\n\
       -F=<char>     Show all discriminating residues in pymol for class <char> (use with -V option).\n\
       -f=<char>     Set full cluster viewing in pymol for class <char> (used with -V option).\n\
       -H            Ignore hydrogens; base distances on heavy atoms only.\n\
       -h            Print histogram of computed (unadjusted) P-values.\n\
       -I=<string>   Interface with subunits <string> clustering (used with method=6 or 7).\n\
       -i=<char>     Ignore this category of residues.\n\
       -J            Use (i.e., don't correct for) Jeffreys priors.\n\
       -K=<int><char> Use <int> residue in chain <char> as starting point (used with -m=1 or -m=3).\n\
       -K=<int1><char><int2> Use <int2> atom in <int1> residue in chain <char> as starting point \n\
       -k=<int><char> Print <int> residue in chain <char> only (for -K=<int><char><int> option).\n\
       -M=<int>      Maxumum number of discriminating positions (default: 25; maximum 100).\n\
       -map=<file>   Homology modeling: \n\
		       maps residues in 1st seq to 2nd seq based on <file>.cma and <file>.pttrns\n\
       -mb=<real>    Minimum interface contact area (for -I option; default: 1.5 squared Angstroms)\n\
       -m=<int>      Method (default: 2)\n\
                       1: Network clustering (= H-bond networks with -b or -B option).\n\
                       2: Core clustering.\n\
                       3: Spherical clustering.\n\
                       4: Spherical clustering for main class; core clustering for others.\n\
                       5: P-value for surface residues enrichment.\n\
                       6: Fixed interface clustering (resets -P=1; requires using -I= option).\n\
                       7: Optimized interface clustering(resets -P=1; requires using -I= option).\n\
			  (This still needs to be properly adjusted; currently too conservative.)\n\
                       8: Spherical clustering using 1 key residue only.\n\
                       9: Interface + adjacent H-bond network (not yet implemented).\n\
       -o=<char>     Only report this category of residues.\n\
       -overlap      Compute the overlap between input sets.\n\
       -P=<int>      Number of top positions to use as starting points (with -m=1..4; default: 5).\n\
       -p            Permute instead of shuffle (use with -S= option).\n\
       -putall       Show all hits, not just the best.\n\
       -S=<int>      Randomize input to perform <int> simulations.\n\
       -seed=<int>   Provide random seed (for simulations).\n\
       -table        Output a *.tbl file.\n\
       -U            Use the union of signficant hits to construct consensus cluster.\n\
       -V            Create VSI output files to structurally visualize results.\n\
       -X            Exclude CH-pi, NH-pi bonds when using -b option.\n\
       -z            dummy\n\n"

BooLean	scl_typ::OverlappingSeqs(Int4 &offset, Int4 MinOverlap, Int4 MaxMisMatch, e_type E1, e_type E2)
/*************************************************************
return TRUE if E1 & E2 are overlapping fragments of the same sequence.
This sets offset > 0 if query starts before pdbseq else it sets dcaOS <= 0.
It also allows for 'X' residues in pdb files...

Start:					end1=21; end2=24; MinOverlap=5;
       	E1  ----+----+----+----+-     	start1=end1-MinOverlap+1 = 21-5+1=17;
                            |||||	start2 = 1;
	E2                  ----+----+----+----+----           
		:	:	:		(offset= start1 - start2 = 17 -1 = 16)
		:	:	:		( add 16 to second seq)
	E1  ----+----+----+----+-	
	    |||||||||||||||||||||		start1=start2=1 offset = 1 - 1 = 0
	E2  ----+----+----+----+----
		:	:	:	
		:	:	:	
End:					
	E1                     ----+----+----+----+-
                               |||||	
	E2  ----+----+----+----+----    start2=end2 - MinOverlap +1 = 24 -5 + 1 = 20       
					offset = 1 - 20 = -19 (add 19 to first seq)
 (move to sequence.cc eventually)
 *************************************************************/
{
	Int4	start1,start2,end1,end2;
	Int4	s1,s2,e1,e2,stop1,stop2;
	unsigned char	*sq1,*sq2;
	Int4    mismatches=0;

	end1=LenSeq(E1); end2=LenSeq(E2); sq1=SeqPtr(E1); sq2=SeqPtr(E2);
	stop1=end1 - MinOverlap +1; stop2=end2 - MinOverlap +1;
	// for(start1=stop1,start2=1; start1 >= 1 && start2 <= stop2; start1--,start2++)
        for(start1=stop1 ; start1 >= 1; start1--)
	   {
           for(start2=1; start2 <= stop2; start2++){

		mismatches=0;
		for(s1=start1, s2=start2; s1 <= end1 && s2 <= end2; s1++,s2++){
	       		if(sq1[s1]==0 || sq2[s2] == 0) continue;	// one or two == 'X'.
	       		if(sq1[s1] != sq2[s2]){ 
			   mismatches++;
			   if(sq1[s1] != 0 && sq2[s2] != 0 && mismatches > MaxMisMatch) break; 
			}	// allow for missing residues in pdb files...
		}
		if(s1 > end1 || s2 > end2){ offset= start1-start2; return TRUE; }
#if 0
		else {
		   fprintf(stderr,"start1=%d; end1=%d; offset=%d; start2=%d; end2=%d\n",
			start1,end1,start1-start2,start2,end2);
		}
#endif
	   }
	} return FALSE;
}

void	scl_typ::InitDefaults( )
{
	Mode=0;
	FindOverlap=FALSE; dca_vs_pdb=FALSE;
	OnlyOne=0; Ignore=0; UseJeffreys=FALSE; min_buried=1.5; PutTheRest=0; PutAllSCH=FALSE;
	MaxPrint=25; MaxNumPttrns=25; cutoff=0.05; seed=18364592;
	Simulate=0; method=3; xxHG=xHG=zHG=0; efptr=0; vsifp=0;
	AroPiBonds=TRUE; OtherPiBonds=TRUE; BackboneHbonds=FALSE;
	UseHydrogens=TRUE; PutPvals=FALSE; EdgeDrop=-1; MaxMinAdj=7; 
	FullCluster=' ';  FindHbonds=0; ChainI=0; dca_file=0; shuffle=TRUE;
	KeyRes=0; KeyResX=0; KeyChain=0; ShowKey=FALSE; KeyAtmNum=0;
	dcaE=0; dcaOS=0; UnionCons=FALSE; tab_fp=0; BestOnly=FALSE;
}

void	scl_typ::InitPrivate(int argc,char *argv[])
{
	Int4	arg;
	char	str[150],*map_file=0;
	BooLean	openVSI=FALSE,OpenTable=FALSE;

	InitDefaults( );
	if(argc < 2) PrintError(argv[0]);
	TurnOffLicenseStatement();
	AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	program=AllocString(argv[0]);
	infile=AllocString(argv[1]);
        for(arg = 2; arg < argc; arg++){
           if(argv[arg][0] != '-') PrintError(argv[0]);
           switch(argv[arg][1]) {
                case 'B': {
			FindHbonds='B'; 
			if(argv[arg][2] == 'b' && argv[arg][3] == 0) BackboneHbonds=TRUE;
			else if(argv[arg][2] != 0) PrintError(argv[0]);
			} break;
                case 'b': {
		        if(strcmp("-best",argv[arg]) == 0){ BestOnly=TRUE; }
			else {
			  FindHbonds='b'; 
			  if(argv[arg][2] == 'b' && argv[arg][3] == 0) BackboneHbonds=TRUE;
			  else if(argv[arg][2] != 0) PrintError(argv[0]);
			}
			} break;
                case 'C': cutoff=RealOption(argv[arg],'C',0,10,WUSAGE); break;
                case 'c': MaxMinAdj=IntOption(argv[arg],'c',2,15,WUSAGE); break;
                case 'd': EdgeDrop=IntOption(argv[arg],'d',0,16,WUSAGE); break;
		case 'E': {
                        if(sscanf(argv[arg],"-EV=%s",str) != 1) PrintError(argv[0]);
			else dca_file=AllocString(str);
		     } break;
                case 'e': {
                        if(sscanf(argv[arg],"-ev=%s",str) == 1) {
			   dca_file=AllocString(str); dca_vs_pdb=TRUE;
			} else if(argv[arg][2] == 0) efptr=stderr; else PrintError(argv[0]); 
		     } break;
		case 'F': 
                     if(argv[arg][2] == '=' && isupper(argv[arg][3]) && argv[arg][4] == 0){
			PutTheRest=argv[arg][3];
                     } else PrintError(argv[0]); break;
		case 'f': 
                     if(argv[arg][2] == '=' && isupper(argv[arg][3]) && argv[arg][4] == 0){
			FullCluster=argv[arg][3];
                     } else PrintError(argv[0]); break;
                case 'H': { if(argv[arg][2] != 0) PrintError(argv[0]); else UseHydrogens=FALSE; } break;
                case 'h': { if(argv[arg][2] != 0) PrintError(argv[0]); else PutPvals=TRUE; } break;
                case 'I': {
                      if(sscanf(argv[arg],"-I=%s",str) != 1) PrintError(argv[0]);
		      else ChainI=AllocString(str); 
		      for(Int4 i=0; ChainI[i] != 0; i++) if(!isalpha(ChainI[i])) PrintError(argv[0]);
		     } break;
                case 'i': {
                        if(sscanf(argv[arg],"-i=%c",&Ignore) != 1) PrintError(argv[0]);
                        if(!isupper(Ignore)) PrintError(argv[0]);
		     } break;
                case 'J': { if(argv[arg][2] != 0) PrintError(argv[0]); else UseJeffreys=TRUE; } break;
                case 'K': 
		     if(sscanf(argv[arg],"-K=%d%c%d",&KeyRes,&KeyChain,&KeyAtmNum) == 3){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
			if(KeyAtmNum < 1) PrintError(argv[0]);
		     } else if(sscanf(argv[arg],"-K=%d%c",&KeyRes,&KeyChain) == 2){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
		     } else PrintError(argv[0]); break;
                case 'k': 
		     if(sscanf(argv[arg],"-k=%d%c",&KeyRes,&KeyChain) == 2){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
			else ShowKey=TRUE;
		     } else PrintError(argv[0]); break;
                case 'M': MaxNumPttrns=IntOption(argv[arg],'M',5,10000,WUSAGE); break;
                case 'm': 
		   if(sscanf(argv[arg],"-mb=%f",&min_buried) == 1){
			if(min_buried < 0.0) PrintError(argv[0]);
		   } else if(sscanf(argv[arg],"-map=%s",str) == 1){
			map_file=AllocString(str);
		   } else {
			method=IntOption(argv[arg],'m',1,9,WUSAGE); break;
		   } break;
                case 'P': MaxPrint=IntOption(argv[arg],'P',1,100,WUSAGE); break;
                case 'o': {
			if(strcmp("-overlap",argv[arg])==0) FindOverlap=TRUE;
                        else if(sscanf(argv[arg],"-o=%c",&OnlyOne) == 1){
                           if(!isupper(OnlyOne)) PrintError(argv[0]);
			} else PrintError(argv[0]); 
		   } break;
                case 'p': {
		   if(strcmp("-putall",argv[arg]) == 0){ PutAllSCH=TRUE; }
		   else if(argv[arg][2] != 0) PrintError(argv[0]); 
		   else shuffle=FALSE; 
		   } break;
                case 'S': Simulate=IntOption(argv[arg],'S',1,10000,WUSAGE); break;
		case 's':
                     if(sscanf(argv[arg],"-seed=%u",&seed) != 1) PrintError(argv[0]); break;
		case 't':
		     if(strcmp("-table",argv[arg]) == 0) OpenTable=TRUE; 
		     else PrintError(argv[0]); break;
                case 'U': {
		      if(argv[arg][2] == 0){ UnionCons=TRUE; } else PrintError(argv[0]); 
		   } break;
                case 'V': {
		      if(argv[arg][2] == 0){ openVSI=TRUE; } else PrintError(argv[0]); 
		   } break;
                case 'X': { if(argv[arg][2] != 0) PrintError(argv[0]);
			else AroPiBonds=OtherPiBonds=FALSE; } break;
                default: PrintError(argv[0]);
            }
        }
	if(KeyRes!=0 && method != 3 && method != 1)
		{ print_error("-K=<int> option requires -m=1 or -m=3."); }
	else if(KeyRes && method == 1 && !FindHbonds)
		{ print_error("-K=<int> option with -m=1 requres -b or -B option."); }
	else if(KeyRes && method == 3 && FindHbonds)
		{ print_error("-K=<int> option with -m=3 prohibits -b and -B options."); }
	if(ChainI==0 && method==6) print_error("For -m=6 option use -I=<str> to define interface.");
	if(ChainI==0 && method==7) print_error("For -m=7 option use -I=<str> to define interface. ");
	if(ChainI && method < 6) print_error("Require -m=6 option with -I=<str> to define interface. ");
	if(FindHbonds && method >=6) print_error("Options -b and -B are inconsistent with -m=6-7. ");
	if(dca_file && method > 4) print_error("The -EV option can only be used with -m=1-4.");
	if(!FindOverlap && OpenTable) tab_fp=open_file(infile,".tbl","w");
#if 1   // create DCA Seq...move to scl_typ eventually.
	if(dca_file){
          FILE *fp=0,*ifp=0;
	  char	Str[105],cI,cJ;
	  double Dd;
	  // if using GaussDCA output as input...then need cma file to compute true positions.
	  if((fp=fopen(dca_file,"r")) == NULL){	// then look for dca_file.cma & dca_file.dca as input.
	     fp=open_file(dca_file,".cma","r");
	     cma_typ cma=ReadCMSA(fp,AB); fclose(fp);
	     assert(nBlksCMSA(cma) == 1);
	     Int4 ii,jj,i,j,len=LengthCMSA(1,cma);
	     char *Seq; NEW(Seq, len+5, char); 
	     e_type qE=TrueSeqCMSA(1,cma);
	     ifp=open_file(dca_file,".dca","r");
	     fp=open_file(dca_file,"","w");	// dca_file as output...adjusts pos and adds residues.
	     Int4 pc[12],os=OffSetSeq(qE);
	     double *X; NEW(X,LenSeq(qE)+os+5,double);
             while(fgets(Str,100,ifp) != NULL){
              if(sscanf(Str,"%d %d %d %d %lf",&i,&j,&pc[1],&pc[2],&Dd) == 5){	// PSICOV input
		assert(i > 0 && i <= LenSeq(qE)); assert(j > 0 && j <= LenSeq(qE));
		ii=i; jj=j; // ii=ii+os; jj=jj+os; 
		// ii=TruePosCMSA(1,i,cma)+os; jj=TruePosCMSA(1,j,cma)+os;
		cI=AlphaChar(ResSeq(ii,qE),AB); cJ=AlphaChar(ResSeq(jj,qE),AB);
                fprintf(fp,"%d,%d,%.5lf,0,0,0,0,0,0,0,%c,%c\n",ii,jj,Dd,cI,cJ);
		if(Dd > 0) { X[ii] += Dd; X[jj] += Dd; }
              } else if(sscanf(Str,"%d %d %lf",&i,&j,&Dd) == 3){		// gDCA input
		assert(i > 0 && i <= len); assert(j > 0 && j <= len);
		// ii=TruePosCMSA(1,i,cma)+os; jj=TruePosCMSA(1,j,cma)+os;
		ii=TruePosCMSA(1,i,cma); jj=TruePosCMSA(1,j,cma);
		assert(ii > 0 && ii <= LenSeq(qE)); assert(jj > 0 && jj <= LenSeq(qE));
		cI=AlphaChar(ResidueCMSA(1,1,i,cma),AB); cJ=AlphaChar(ResidueCMSA(1,1,j,cma),AB);
                fprintf(fp,"%d,%d,%.5lf,0,0,0,0,0,0,0,%c,%c\n",ii,jj,Dd,cI,cJ);
		if(Dd > 0) { X[ii] += Dd; X[jj] += Dd; }
	      } else print_error("-EV option input error 0.");
	     } fclose(ifp); fclose(fp); 
	     fp=open_file(dca_file,".clst","w");	// Order DCA residues...
	     dh_type dH=dheap(LenSeq(qE)+5,4); 
	     for(i=1; i <= LenSeq(qE); i++) if(X[i] > 0.0) insrtHeap(i,-(keytyp)X[i],dH);;
	     for(j=0; !emptyHeap(dH); j++){
		i=delminHeap(dH); 
		// fprintf(stdout, "%d = %.5lf\n",i,X[i]); 
		if(j > 0) fprintf(fp, ",%d",i); else fprintf(fp, "C=%d",i); 
		if(j > 50) break;
	     } free(X); Nildheap(dH); fprintf(fp,"\n"); fclose(fp);
	     TotalNilCMSA(cma); fp=open_file(dca_file,"","r");
	  }
	  // Read in EVcouplings-formatted input file...
          Int4  max_ij=0,min_ij=9999;
          unsigned char *seq;     NEW(seq,1005,unsigned char); // 30-500 is evfold max...
          while(fgets(Str,100,fp) != NULL){
                Int4	i,j,dm[12];    // dummy variables.
                if(sscanf(Str,"%d,%d,%lf,%d,%d,%d,%d,%d,%d,%d,%c,%c,",
                   &i,&j,&Dd,&dm[1],&dm[2],&dm[3],&dm[4],&dm[5],&dm[6],&dm[7],&cI,&cJ) != 12)
                        print_error("GetMtrxDCA() input error 1.");
                assert(i < 1000 && j < 1000 && i > 0 && j > 0);
                if(i > max_ij) max_ij=i; if(j > max_ij) max_ij=j;
                if(i < min_ij) min_ij=i; if(j < min_ij) min_ij=j;
                seq[i]=AlphaCode(cI,AB); seq[j]=AlphaCode(cJ,AB);
                // EqSeq(i,AlphaCode(cI,AB),dcaE); EqSeq(j,AlphaCode(cJ,AB),dcaE);
          } fclose(fp);
          dcaE=MkSeq("DCA seq",max_ij-min_ij +1,seq + min_ij -1); free(seq);
          SetOffSetSeq(min_ij -1,dcaE);
          PutSeq(stderr,dcaE,AB); 
	}
#endif
#if 1	// test map_cma...
	if(map_file){
	// 3. Read in pttrn file for Node 14 (= column numbers in alignment).
	// 4. Convert position to a new *.klst file for the pdb file named in Node2 cma.
	// 5. Output this file, exit and call program using new *klst file.
	// NOte: create alignment files using hybrid program (not hierview).
	  // 0. Open files...
	  Int4	x,y,z,Z,n,N,Number=0,*Pos[102];
	  Int4	i,j,k,s,s1,s2,r1,r2,pos[4];
	  sprintf(str,"%s.cma",map_file);
	  FILE *fp = open_file(map_file,".cma","r");
	  char	*sp,Str[1004],xstr[1004];
          cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	  if(Number != 2) print_error("-map option input error");
	  cma_typ xcma,map_cma=IN_CMA[1];
	  // multi read == 2 files; Node 22 and Node 2...
	  assert(nBlksCMSA(map_cma)==1);
	  assert(NumSeqsCMSA(map_cma) > 1);

	  xcma=IN_CMA[2];
	  assert(nBlksCMSA(xcma)==1); assert(NumSeqsCMSA(xcma) == 1);
	  e_type E=TrueSeqCMSA(1,xcma),E2=TrueSeqCMSA(1,map_cma); 
	  if(!IdentSeqs(E,E2) || OffSetSeq(E) != OffSetSeq(E2)){
		 print_error("-map option: Seqs inconsistent");
	  }
	  N=LenSeq(E) + OffSetSeq(E);
	  sst_typ *sst[102];
	  fp = open_file(map_file,".pttrns","r");
	  
	  //**********************************************************************
	  //***** 1. Map Node 14 columns to true sema4A sites for all categories ****
	  //**********************************************************************
	  Int4 *PosC2A; NEW(PosC2A, LengthCMSA(1,xcma)+5,Int4);
	  for(s=1; s <= LengthCMSA(1,xcma); s++){
		assert(PosSiteCMSA(1, 1, pos, xcma));
		if(IsDeletedCMSA(1,pos[1]+s-1,xcma)) continue;
		s1=TruePosCMSA(1,s,xcma) + OffSetSeq(TrueSeqCMSA(1,xcma));
		  r1=ResidueCMSA(1,1,s,xcma);
		  fprintf(stderr,"%d: %c%d\n",s,AlphaChar(r1,AB),s1);
		PosC2A[s]=s1;	
	  }

	  //**********************************************************************
	  //************ 2. Map sema4A (Node 2) to sema4D (Node 2). *************
	  //**********************************************************************
	  Int4 *PosA2H; NEW(PosA2H, N+5,Int4);
	  fprintf(stderr,"Source --> Target\n");
	  for(s=1; s <= LengthCMSA(1,map_cma); s++){
		assert(PosSiteCMSA(1, 1, pos, map_cma));
		if(IsDeletedCMSA(1,pos[1]+s-1,map_cma)) continue;
		assert(PosSiteCMSA(1, 2, pos, map_cma));
		if(IsDeletedCMSA(2,pos[1]+s-1,map_cma)) continue;
		s1=TruePosCMSA(1,s,map_cma) + OffSetSeq(TrueSeqCMSA(1,map_cma));
		s2=TruePosCMSA(2,s,map_cma) + OffSetSeq(TrueSeqCMSA(2,map_cma));
		r1=ResidueCMSA(1,1,s,map_cma);
		r2=ResidueCMSA(1,2,s,map_cma);
		fprintf(stderr,"%d.%c%d --> %c%d\n",s,AlphaChar(r1,AB),s1,AlphaChar(r2,AB),s2);
		assert(s1 <= N);
		PosA2H[s1]=s2;
	  }
	  //*************************************************************************
	  //***** 3. Get Node 14 patterns for sema4A and map to sema4D (Node 2) ****
	  //*************************************************************************
	  static char c,level[]=" YROMGCBDDDDDDDDD  ";
	  for(i=0; fgets(Str,1000,fp) != NULL; ){
		i++; if(i >= 100) print_error("-map option: too many patterns");
		if(sscanf(Str,"%d: %s\n",&s,xstr) != 2) print_error("-map option: pattern too long");
		NEW(Pos[i], N+5, Int4); NEW(sst[i], N+5, sst_typ);
		n=ParseResidueSets(xstr,Pos[i],sst[i],AB,"-map option: pattern parsing error");
		Pos[i][0]=n;
	  } fclose(fp);
	  Int4 end=i;
	  //**********************************************************************
	  //*** 3. Map Node 14 patterns to true sema4D sites for each category ***
	  //**********************************************************************
	  for(k=1,i=end; i > 0; i--,k++){
		if(k < 9) c=level[k]; else c = 'W'; 
		fprintf(stdout,"%c=",c); n=Pos[i][0];
		for(z=Z=0,j=1; j <= n; j++){
		    x=Pos[i][j]; y=PosC2A[x]; z=PosA2H[y]; 
		    if(Z > 0 && z > 0) fprintf(stdout,","); 
		    if(z > 0){ fprintf(stdout,"%d",z); Z=1; }
		} fprintf(stdout,"\n");
	  }
	  //**********************************************************************
	  //*** 4. Print out residues for a VSI file.
	  //**********************************************************************
	  fprintf(stdout,"\n\n\n");
	  for(k=1,i=end; i > 0; i--,k++){
		if(k < 9) c=level[k]; else c = 'W'; 
		n=Pos[i][0];
		for(z=Z=0,j=1; j <= n; j++){
		    x=Pos[i][j]; y=PosC2A[x]; z=PosA2H[y]; 
		    if(Z > 0 && z > 0) fprintf(stdout,","); 
		    if(z > 0){
	  	        E=TrueSeqCMSA(2,map_cma);  
			Int4 p=z-OffSetSeq(E);
			char rr=ResSeq(p,E);
			fprintf(stdout,"%c%d.%c",AlphaChar(rr,AB),z,c); Z=1; 
		    }
		} fprintf(stdout,"\n");
	  }
	  // X. Free up input files.
	  TotalNilCMSA(map_cma); free(map_file); free(IN_CMA);
	  // X. free infile and replace with outfile...create...and exit; call again...
	  exit(1);
	}
#endif
	if(ChainI && method == 6){ MaxPrint=1; }
	if(ChainI && method == 7){ MaxPrint=1; }
	if(openVSI){
	    if(ChainI) sprintf(str,"%s_%d%s",argv[1],method,ChainI); 
	    else if(FindHbonds) sprintf(str,"%s_%d%c",argv[1],method,FindHbonds);
	    else sprintf(str,"%s_%d",argv[1],method); 
	    vsifp=open_file(str,".VSI","w");
	}
	if(EdgeDrop < 0) EdgeDrop=MaxMinAdj + 1;	// Don't use this parameter
	time1=time(NULL);
        if(seed == 18364592) seed = (unsigned int) (time(NULL)/2);
	sRandom(seed);
        for(arg=0; arg < argc; arg++) fprintf(stdout,"%s ",argv[arg]); 
	if(Simulate) fprintf(stderr,"seed = %u\n",seed); else fprintf(stdout,"\n");
}

#define PUSAGE "Usage: sipris sprs_file mode [options]\n\
   Note: First add hydrogens to pdb files using the Reduce program by J. Michael Word.\n\
     (http://kinemage.biochem.duke.edu/software/reduce.php).\n\
   Possible modes: \n\
      H        hydrogen bond network (sidechain-to-sidechain and sidechain-to-backbone)\n\
      B        hydrogen bond network (also includes backbone-to-backbone H-bonds)\n\
      C        core clustering\n\
      S        spherical clustering\n\
      P=<str>  predefined clustering (<str> = contacting chain(s); e.g., P=BDF)\n\
   options: \n\
       -C=<real>     p-value cutoff (default: 0.05).\n\
       -c=<int>      Minimum adjacency cutoff for core clustering (default: 7; range 2-15).\n\
       -EV=<file>    Use EVfold CouplingScores.csv to compute structural contacts (C or S modes only).\n\
       -K=<int><char> Use <int>th residue in chain <char> as starting point (H, B or S modes only)).\n\
       -K=<int1><char><int2> Use <int2>th atom in <int1>th 'residue' in chain <char> as starting point \n\
       -k=<int><char> Show atom numbers for <int>th 'residue' in chain <char>.\n\
       -M=<int>      Maxumum number of discriminating positions (default: 25; maximum 100).\n\
       -m=<real>     Redefine the minumum residue contact area (default: 1.5 squared Angstroms).\n\
       -P=<int>      Number of top positions to use as starting points (default: up to 25).\n\
       -pml          Create a PyMOL script (for the 1st pdb file in sprs_file).\n\
       -S=<int>      Randomize input to perform <int> simulations.\n\
       -seed=<int>   Provide random seed (for simulations).\n\
                     \n\n"

void    scl_typ::PrintError(char *program_name)
{
	if(strcmp(program_name,"soprise") == 0) print_error(WUSAGE);
	else if(strcmp(program_name,"Soprise") == 0) print_error(WUSAGE);
	else { PrintLicenseStatement("sipris v1.0.2"); print_error(PUSAGE); }
}

void	scl_typ::InitPublic(int argc,char *argv[])
{
	Int4	arg;
	char	str[150],*map_file=0;
	BooLean	openVSI=FALSE,OpenTable=FALSE;

	InitDefaults( );
	if(argc < 3) PrintError(argv[0]);
	TurnOffLicenseStatement();
	AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);
	program=AllocString(argv[0]);
	infile=AllocString(argv[1]);
	Mode=argv[2][0];
	if(Mode != 'P' && argv[2][1] != 0) PrintError(argv[0]);
	switch (Mode){
                case 'B': BackboneHbonds=TRUE;
                case 'H': { method=1; FindHbonds='b'; } break;
                case 'C': { method=2; } break;
                case 'S': { method=3; } break;
                case 'P': {
                      if(sscanf(argv[2],"P=%s",str) != 1) PrintError(argv[0]);
		      ChainI=AllocString(str); method=6;
		      for(Int4 i=0; ChainI[i] != 0; i++) if(!isalpha(ChainI[i])) PrintError(argv[0]);
		  } break;
		default: PrintError(argv[0]); break;
	}
        for(arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') PrintError(argv[0]);
           switch(argv[arg][1]) {
                case 'C': cutoff=RealOption(argv[arg],'C',0,10,PUSAGE); break;
                case 'c': MaxMinAdj=IntOption(argv[arg],'c',2,15,PUSAGE); break;
		case 'E': {
                        if(sscanf(argv[arg],"-EV=%s",str) != 1) PrintError(argv[0]);
			else dca_file=AllocString(str);
		     } break;
		case 'e': {
                        if(sscanf(argv[arg],"-ev=%s",str) == 1) {
			   dca_file=AllocString(str); dca_vs_pdb=TRUE;
			} else PrintError(argv[0]); 
		     } break;
                case 'K': 
		     if(sscanf(argv[arg],"-K=%d%c%d",&KeyRes,&KeyChain,&KeyAtmNum) == 3){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
			if(KeyAtmNum < 1) PrintError(argv[0]);
		     } else if(sscanf(argv[arg],"-K=%d%c",&KeyRes,&KeyChain) == 2){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
		     } else PrintError(argv[0]); break;
                case 'k': 
		     if(sscanf(argv[arg],"-k=%d%c",&KeyRes,&KeyChain) == 2){
			if(!isalpha(KeyChain)) PrintError(argv[0]);
			else ShowKey=TRUE;
		     } else PrintError(argv[0]); break;
                case 'M': MaxNumPttrns=IntOption(argv[arg],'M',5,10000,PUSAGE); break;
                case 'm': min_buried=RealOption(argv[arg],'m',0.1,10,PUSAGE); break;
                case 'P': MaxPrint=IntOption(argv[arg],'P',1,100,PUSAGE); break;
                case 'p': {
		     if(strcmp("-pml",argv[arg]) == 0){ openVSI=TRUE; } else PrintError(argv[0]); 
		   } break;
                case 'S': Simulate=IntOption(argv[arg],'S',1,10000,PUSAGE); break;
		case 's':
                     if(sscanf(argv[arg],"-seed=%u",&seed) != 1) PrintError(argv[0]); break;
                default: PrintError(argv[0]);
            }
        }
	if(KeyRes && method != 3 && method != 1) { print_error("For the -K=<int> option use mode H, B or S."); }
	if(dca_file && Mode != 'C' && Mode != 'S') print_error("For the -EV option use mode C or S.");
//*********** Below here could be combined with InitPrivate() ***************
	// create DCA Seq...move to separate routine eventually.
	if(dca_file){
          FILE *fp=0,*ifp=0;
	  char	Str[105],cI,cJ;
	  double Dd;
	  fp=open_file(dca_file,"","r");
          Int4  max_ij=0,min_ij=9999;
          unsigned char *seq;     NEW(seq,1005,unsigned char); // 30-500 is evfold max...
          while(fgets(Str,100,fp) != NULL){
                Int4	i,j,dm[12];    // dummy variables.
                if(sscanf(Str,"%d,%d,%lf,%d,%d,%d,%d,%d,%d,%d,%c,%c,",
                   &i,&j,&Dd,&dm[1],&dm[2],&dm[3],&dm[4],&dm[5],&dm[6],&dm[7],&cI,&cJ) != 12)
                        print_error("GetMtrxDCA() input error 1.");
                assert(i < 1000 && j < 1000 && i > 0 && j > 0);
                if(i > max_ij) max_ij=i; if(j > max_ij) max_ij=j;
                if(i < min_ij) min_ij=i; if(j < min_ij) min_ij=j;
                seq[i]=AlphaCode(cI,AB); seq[j]=AlphaCode(cJ,AB);
                // EqSeq(i,AlphaCode(cI,AB),dcaE); EqSeq(j,AlphaCode(cJ,AB),dcaE);
          } fclose(fp);
          dcaE=MkSeq("DCA seq",max_ij-min_ij +1,seq + min_ij -1); free(seq);
          SetOffSetSeq(min_ij -1,dcaE);
          PutSeq(stderr,dcaE,AB); 
	}
	if(ChainI && method == 6){ MaxPrint=1; }
	if(openVSI){
#if 0
	    if(ChainI) sprintf(str,"%s_%d%s",argv[1],method,ChainI); 
	    else if(FindHbonds) sprintf(str,"%s_%d%c",argv[1],method,FindHbonds);
	    else sprintf(str,"%s_%d",argv[1],method); 
	    vsifp=open_file(str,".VSI","w");
#else
	    vsifp=tmpfile();
#endif
	}
	if(EdgeDrop < 0) EdgeDrop=MaxMinAdj + 1;	// Don't use this parameter
	time1=time(NULL);
        if(seed == 18364592) seed = (unsigned int) (time(NULL)/2);
	sRandom(seed);
        for(arg=0; arg < argc; arg++) fprintf(stdout,"%s ",argv[arg]); 
	if(Simulate) fprintf(stderr,"seed = %u\n",seed); else fprintf(stdout,"\n");
}

