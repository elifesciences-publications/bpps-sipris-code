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

#include "gmb_typ.h"

e_type	gmb_typ::DelimitSeqAln(Int4 i, Int4 &strt, Int4 &end, Int4 &len, Int4 &del)
{
	cma_typ CMA=SSX->RtnCMA();
	Int4	blk=nBlksCMSA(CMA),Len=LengthCMSA(blk,CMA);
	e_type	E=TrueSeqCMSA(i,CMA); strt=TruePosCMSA(i,1,1,CMA); end=TruePosCMSA(i,blk,Len,CMA);
	for(del=0; del < LengthCMSA(1,CMA) && IsDeletedCMSA(1,i,del+1,CMA); del++) ;
        if(del >= LengthCMSA(1,CMA)) return 0;  // if entire first block deleted then skip.
        else if(del > 0) strt++;  // deletions start at the preceding residue;
	len = end-strt+1; 
	return E;
}

Int4	*gmb_typ::MultiAlnToLongest(FILE *fp, Int4 *ListE)
{
	sap_typ	sap=0,head=0,tail=0;
        Int4	i,j,ssq,qsq,qst,sst,os,end,x,max,J,sq=ListE[0],Qst,Qend,Sst,Send,len;
	Int4	s1,s2,e1,e2,len1,len2,d1,d2,S1,S2,*OffSet;
        e_type	sE,qE;
	cma_typ CMA=SSX->RtnCMA();

	// 1. find the longest sequence.
	for(max=0,J=0,j=1; j <= sq; j++){
		sE=TrueSeqCMSA(j,CMA);
		if(LenSeq(sE) > max){ max = LenSeq(sE); J=j; }
	} if(J != 1){ qsq=ListE[J]; ListE[J]=ListE[1]; ListE[1]=qsq; } else qsq=ListE[1];
	
	// 2. Align the rest of the sequences to the longest.
	NEW(OffSet,sq+5,Int4);
	assert((qE=DelimitSeqAln(qsq,qst,e1,len1,d1)) != 0);
        for(j=1; j <= sq; j++){
                ssq=ListE[j]; 
	        assert((sE=DelimitSeqAln(ssq,sst,e2,len2,d2)) != 0);
		S1=qst-d1; S2=sst-d2; os=S1-S2; OffSet[j]=os;

		// if(fp) this->PutDiagonalSeq(fp,os,qE,sE); // DEBUG
		len=GetDiagonalEnds(os,qE,sE,Qst,Sst,Qend,Send);
		// if(fp) fprintf(fp,"%d: %d vs %d Qst=%d; Sst=%d; os=%d\n",s,qsq,ssq,Qst,Sst,os);
		if(fp){   // WARNING: sap_typ start at zero while cma_typ starts at one!!!
        	  char	*op; NEW(op,len+9,char); op[0]='E'; op[1]='M';
		  for(x=2; x <= len; x++) op[x]='m'; op[x]='E';
		  sap = ToGSeqAlign(strlen(op),op,qE,sE,Qst-1,Sst-1);
                  if(head==0){ head=tail=sap; } else { tail->next=sap; tail=sap; }
		  free(op);
		}
#if 0	// Error checking...
	     //if(OffSet[j] >= LenSeq(sE))
	     if(ssq==2293)
	     {
		this->PutDiagonalSeq(stderr,os,qE,sE); // DEBUG
		fprintf(stderr,"%d: os=%d; OffSet[%d] = %d; adj=%d; net=%d; len=%d\n",
			ListE[j],os,j,OffSet[j],OffSetSeq(qE),os+OffSetSeq(qE),LenSeq(sE));
		// assert(OffSet[j] < LenSeq(xE));
	     }
#endif
	} 

	// 3. Print out the alignment.
	if(fp){
          fprintf(fp,"cluster (%d seqs):",sq);
          for(j=1; j <= sq && j <=5; j++) fprintf(fp," %d",ListE[j]);
          fprintf(stderr,"\n"); PutSeqInfo2(fp,qE); fprintf(fp,"\n");
          PutMultiGSeqAlign(fp,head,60,AB); FreeGSeqAlignList(head);
	} return OffSet;
}

BooLean	IsOnTheList(Int4 *lst)
{
	BooLean found=FALSE;
	Int4	xx,sq=lst[0];
	for(xx=1; xx <= sq; xx++){
		if(lst[xx] == 1792) return TRUE;
		if(lst[xx] == 1886) return TRUE;
		if(lst[xx] == 1904) return TRUE;
		if(lst[xx] == 14185) return TRUE;
#if 0
		if(lst[xx]==3304) return TRUE;	// 70732608
		if(lst[xx]==3300) return TRUE;	// 1CHM_B
		if(lst[xx]==3444) return TRUE; 	// 226360923 
		if(lst[xx]==1854) return TRUE; 	// 491025417 
		if(lst[xx]==3435) return TRUE; 	// 495820312 
		if(lst[xx]==3310) return TRUE; 	// 89055103 
		if(lst[xx]==3377){ found=TRUE; break; }
		if(lst[xx]==3317){ found=TRUE; break; }
		if(lst[xx]==1885){ found=TRUE; break; }
#endif
	} return FALSE;
}

e_type	gmb_typ::MultiAlnToConsensus(FILE *fp, Int4 *OffSet, Int4 *ListE)
{
        Int4	i,j,r,R,s,os,x,min,J,sq=ListE[0],end;
	Int4	len,lenN,lenC,max,maxN,maxC;
        e_type	sE,cE=0;
	h_type	HG=0;

	HG=Histogram("Offsets for sequence cluster", -30,30,1.0);
	// 1. Find the consensus length.
/***************************************************************************************
                        !   .    |    .    |     len os lenC lenN  os-min  Sst  end
        sq1          ---XXXXXXXXXXXXXXXXXXXX---	  20  0   20    0     3      4    20
        sq2          xxxXXXXXXXXXXXXXXXXXX-----   21 -3   18    3     0      1    21
        sq3          --------XXXXXXXXXXXXXXXXXX   18  5	  23    0     8      9
        consensus    xxxXXXXXXXXXXXXXXXXXXXXXXX   26  =   23 +  3     0      1
	             |                                            min_os = -3
                     os - min = 0 
 ***************************************************************************************/
	cma_typ CMA=SSX->RtnCMA();
	for(min=INT4_MAX,max=maxN=maxC=0,j=1; j <= sq; j++){
		s=ListE[j]; os=OffSet[j];
		sE=TrueSeqCMSA(s,CMA);
		IncdHist(os,HG);	//  OffSet == S1-S2; if < 0 --> sE longer Nterm.
		if(min > os) min = os;
		lenN=MAXIMUM(Int4,-os,0);
		if(maxN < lenN) maxN=lenN;
		lenC=LenSeq(sE) + os;	// length from fiducial mark to end.
		if(maxC < lenC) maxC=lenC;
		if(max < LenSeq(sE)) max = LenSeq(sE);
	} if(fp && HG) PutHist(fp,60,HG); if(HG) NilHist(HG);
	assert(min <= 0); len=maxN+maxC;
	if(fp) fprintf(fp,"range=%d..%d; len=%d; max_sqlen=%d\n",maxN,maxC,len,max);

	cE=MkEmptySeq(1,"consensus seq.",len);
	for(r=0; r <= nAlpha(AB); r++) NEW(CsqCnts[r],len+5,Int4);
	HG=Histogram("maximum residue freqs at each position",1,sq,1.0);
	for(i=1; i <= sq; i++){
	   s=ListE[i];
	   sE=TrueSeqCMSA(s,CMA);
	   os=OffSet[i] - min; 	// reset offsets relative to consensus start position.
	   assert(os >= 0);
	   end=LenSeq(sE);
	   for(j=os+1,x=1; x <= end; j++, x++){
		assert(j <= len);
		r=ResSeq(x,sE); CsqCnts[r][j]++; 
	   }
	   OffSet[i]=os;  // reset relative to consensus starting point.
	}
	for(j=1; j <= len; j++){
	   for(max=0,R=0,r=0; r <= nAlpha(AB); r++){
		if(CsqCnts[r][j] > max){ max=CsqCnts[r][j]; R=r; }
	   } EqSeq(j,R,cE);
	   IncdHist(max,HG);	//  maximum N-terminal extension beyond start.
	}
	if(fp && HG) PutHist(fp,60,HG); if(HG) NilHist(HG);
#if 0	// DEBUG...
	// if(fp)
	{
	   for(i=1; i <= sq; i++){
		s=ListE[i]; sE=TrueSeqCMSA(s,CMA);
		// if(OffSet[i] > LenSeq(sE)) 
	        // if(s==2293)
		{
		  this->PutDiagonalSeq(stderr,OffSet[i],cE,sE);
		  fprintf(stderr,"%d(%d): os=%d\n",i,s,OffSet[i]);
		}
	   }
	}
#endif
	for(r=0; r <= nAlpha(AB); r++) free(CsqCnts[r]);
	return cE;
}

Int4	gmb_typ::ScoreSeqAlnSMatrix(char *operation, e_type sbjE, Int4 strt)
{
        Int4    o,i,j,k,score;
	smx_typ *smx=SSX->RtnSMX();
	smx_typ SMX=smx[1];
	unsigned char	*seq1,*seq=SeqPtr(sbjE)+strt-1;
	// a_type	A = SMatrixA(SMX);

	/** get concensus sequence for smatrix **/
	NEW(seq1, LenSMatrix(SMX) +3, unsigned char); 
	// MaxSegSMatrix(seq1, SMX);	// don't see that this is needed...
	for(score=0,k=o=i=1,j=1; operation[o] != 'E'; o++){ 
		switch(operation[o]){
		    case 'M': case 'm':
			j++; score+=ValSMatrix(k,seq[j],SMX); 
			i++; k++;
			break;
		    case 'D': case 'd': i++;  k++; break;
		    case 'I': j++; break;
		    case 'i': j++; break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		} assert(j <= LenSeq(sbjE));
		if(k > LenSMatrix(SMX)) break; 
	} free(seq1); 
	return score;
}

Int4	gmb_typ::AlignLenOperation(char *Op)
// return the length of the alignment within the sequence from start;
// WARNING: Assumes only one block...
{
        Int4    N=0,i,j,S=0;
        char    state;
	assert(Op[0]=='E');
        for(i=1; Op[i] != 'M' && Op[i]!='D'; i++) ;
	for(S=0,state=Op[i]; state != 'E'; i++,state=Op[i]){
	   switch (state){
		case 'I':  N++; break;
		case 'M':  N++; S++; break;
		case 'm':  N++; break;
		case 'D':  S++; break;
		case 'd':  break;
		case 'i':  // break;
		default: print_error("gmb_typ::AlignLenOperation(): input error"); 
		  break;
	   } 
	   if(S != 1){ fprintf(stderr,"Op=%s\n",Op); assert(S == 1); }
	} return N;
}

Int4	gmb_typ::GetMultipleTraces(e_type csqE,double Temp)
// Find all of the possible alignments against the consensus sequence csqE.
{
	e_type	CsqE=CopySeq(csqE);
	if(OPS != 0) delete OPS;
	OPS= new ops_typ();
	Int4 strt,scr,i,j,len,end; 
	for(i=1; i < OPS->size; i++){
		char *op=SSX->GapAlnTrace(CsqE,Temp,strt,scr);
		assert(op[1]=='M' || op[1]=='D');
		len = this->AlignLenOperation(op);
		if(i > 1 && (scr < 5000 || len < 5)){ free(op); break; }
                OPS->Add(op,strt,scr);
		end=strt+len-1;
		// this->PutGappedAln(stderr,CsqE,op,strt); OPS->Put(stderr,i);
		for(j=strt; j <= end; j++){ EqSeq(j,0,CsqE); }
		// PutSeq(stderr,CsqE,AB);
        } NilSeq(CsqE);
#if 0	// DEBUG...
	OPS->Put(stderr); OPS->Sort(); OPS->Put(stderr);
#else
	OPS->Sort();
#endif
	return i;
}

char	*gmb_typ::ConvertOperation(char *operation, e_type sE, Int4 mst, Int4 offset)
// convert the operation for the master seq to a corresponding operation for the slave seq.
/*******************************************************************
                             mst=9
                             |                    
        master seq;  kinrkispDLKVKPNFItkkrLKEDEINFdkLRSYRLdrvkselkknnle*
        master:             EMmmmmmmmmiiiiMmmmmmmmiiMmmmmmE
        slave (os=-13):     EDddddmmmmiiiiMmmmmmmmiiMmmmddE  
        slave seq;   xxxxxxxxXXXXXPNFItkkrLKEDEINFdkLRSYRLdrvkselkknnle*
                             |   012...
                             sst= mst + os = 9 + -13 = -4.

        slave (os=-4):      EMmmmmmmmmiiiiMmmmmmmmiiMmmmddE  


        master:             EDddmmmmmmiiiiMmmmmmmmiiMmmmmmE
        slave (os=?):       EDddddmmmmiiiiMmmmmmmmiiMmmmddE  

 *******************************************************************/

{
	Int4	i,j,sst,len=strlen(operation);
	Int4	maxlen=len+LenSeq(sE)+2;
	char	c,*sop; NEW(sop,maxlen+5,char);	// slave operation.

	FILE *fp=0; // fp=stderr;
#if 0	// DEBUG...
	StrSeqID(str, 100, sE);
	if(strstr(str,"512922887")){ fp=stderr; }
#endif
	// 1. 
	assert(operation[0] == 'E'); assert(operation[len-1] == 'E');
	// assert(operation[len-2] == 'm');  // disallow deletions on ends of master for now.
	assert(offset >= 0);
	// assert(offset < LenSeq(sE));
	sst = mst - offset;
	sop[0]='E'; i=1; j=1;
	if(1 || operation[1]=='M'){	
	  if(sst < 1){  // slave starts with a deletion.
// fprintf(stderr,"mst=%d; offset=%d; sst=%d\n",mst,offset,sst);
	    for(sop[j]='D',i++,j++,sst++; sst < 1 && operation[i] != 'E'; ){
		if(operation[i] == 'D' || operation[i] == 'M'){
			sop[j]='D'; i++; j++; sst++;
		} else {
			// assert(operation[i] == 'm');
			sop[j]='d'; i++; j++; sst++;
		}
	    }
	  }
	  for( ; i < len && (c=operation[i]) != 'E'; i++,sst++){
	   if(sst > LenSeq(sE)){	// sE sequence ends...
	     switch(c){
		case 'M': sop[j]='D'; j++; break;
		case 'D': sop[j]='D'; j++; break;
		case 'm': sop[j]='d'; j++; break;
		case 'd': sop[j]='d'; j++; break;
		case 'I': case 'i': sst--; break;
		default:  
		  fprintf(stderr,"mst=%d; c='%c'=%d; i=%d; operation =%s\n",mst,c,c,i,operation);
		  fprintf(stderr,"sst=%d; j=%d <= maxlen=%d; LengSeq()=%d; offset=%d; strlen=%d\n",
			sst,j,maxlen,LenSeq(sE),offset,len);
		  PutSeq(stderr,sE,AB);
		  fprintf(stderr,"slave operation =%s\n",sop);
		  assert(sst <= LenSeq(sE));
		  print_error("ConvertOperation() input error"); break;
	     }
	   } else { sop[j]=c; j++; }
	  } sop[j]='E';
	   assert(j <= maxlen);
	} else {	// disallow deletions in start of master for now.
	   fprintf(stderr,"operation =%s\n",operation);
	   assert(operation[1]=='M');	// master starts with a match.
	}
	if(fp){
	  fprintf(stderr,"mst=%d; operation =%s\n",mst,operation);
	  fprintf(stderr,"sst=%d; LengSeq()=%d; offset=%d\n",sst,LenSeq(sE),offset);
	  PutSeq(stderr,sE,AB);
	  fprintf(stderr,"slave operation =%s\n",sop);
	}
	return sop;
}

Int4	**gmb_typ::GetClusters(FILE *fp, Int4 **&OffSet, Int4 &Nset)
{
	cma_typ CMA=SSX->RtnCMA();
        Int4	i,j,offset,*self_score,scr,score,N=NumSeqsCMSA(CMA);
	Int4	len,len1,len2,s1,s2,e1,e2,Len,blk,d1,d2,S1,S2,dnet;
	Int4	seti,setj;
	Int4	dropoff=8,hits=0;
	double	d,score_cutoff=0.70;
	e_type	E1,E2;
	ds_type	sets=DSets(N);
	blk=nBlksCMSA(CMA); Len = LengthCMSA(blk,CMA);
	gss_typ *gss=gssCMSA(CMA);
	h_type	HG=0;
	if(fp) HG=Histogram("Diagonal score to self score percentages", 10,100,2.0);
	NEW(self_score,N+5,Int4);
        for(i=1; i <= N; i++){
	  E1=TrueSeqCMSA(i,CMA);
	  self_score[i]=this->Extend(SeqPtr(E1)+1,SeqPtr(E1)+1,LenSeq(E1),AlphaR(AB),dropoff);
	}
        for(i=1; i <= N; i++){

	  if((E1=DelimitSeqAln(i,s1,e1,len1,d1)) == 0) continue;

	  seti=findDSets(i,sets);
	  // self_score=this->Extend(SeqPtr(E1)+s1,SeqPtr(E1)+s1,len1,AlphaR(AB),dropoff);
// char *qop=gss->Operation(0,i); 	// this does not matter for sampling together!!!
// I only need to save the offsets from the master sequence and may need to add deletions??!!
          for(j=i+1; j <= N; j++){

	        if((E2=DelimitSeqAln(j,s2,e2,len2,d2)) == 0) continue;

		if((len1 + d1) != (len2 + d2)) continue;
		S1=s1-d1; S2=s2-d2; offset=S1-S2;
		if(offset < 0){	// S2 > S1.
		  offset = -offset; len=MINIMUM(Int4,LenSeq(E1),LenSeq(E2)-offset);
		  score=this->Extend(SeqPtr(E1)+1,SeqPtr(E2)+offset+1,len,AlphaR(AB),dropoff);
		  offset = -offset;
		} else {	// S1 >= S2.
/***************************************************************************************
				S1=s1-d1=0.
                                s1   S1=s1+(d2-d1)=6. e1
                                |    |                |  len=e1-S1+1
		sq1            -MMMMMMMMMMMMMMMMMMMMMMM  d1=1.
	        sq2            ------MMMMMMMMMMMMMMMMMM	 d2=6.
                               |-d2-||                |
		                     S2=s2-d2=1-6=-5  e2  len2=e2-S2+1;  offset = S1-S2 = 5
			offset = S1-S2=0-5 = -5.
 ***************************************************************************************/
		  len=MINIMUM(Int4,LenSeq(E1)-offset,LenSeq(E2));
		  score=this->Extend(SeqPtr(E1)+offset+1,SeqPtr(E2)+1,len,AlphaR(AB),dropoff);
		}
#if 0
if((i==64 || j == 64) && score > 0) {
	fprintf(stderr,"i=%d: s1=%d; len1=%d, d1=%d. j=%d: s2=%d; len2=%d; d2=%d; score=%d\n",
		i,s1,len1,d1,j,s2,len2,d2,score);
}
#endif
		if(score > 0){
		   scr= MINIMUM(Int4,self_score[i],self_score[j]);
		   d = (double)score/scr;
	           if(HG) IncdHist(100*d,HG);
		   if(d < score_cutoff) continue;
// char *op=gss->Operation(0,j); if(strcmp(qop,op) != 0){ free(op); continue; } else free(op);
		   setj=findDSets(j,sets); seti=linkDSets(seti,setj,sets);
		   hits++; 
#if 0	//======================= DEBUG: look at diagonal alignment ==============
		   if(fp && (i==63 || j == 63 || i==64 || j == 64)) {
                	// offset=s1-s2;
			fprintf(fp,"*** %d(%d..%d) vs %d(%d..%d): score = %d ***\n",i,s1,e1,j,s2,e2,score);
                	this->PutDiagonalSeq(fp,offset,E1,E2);
		   }
#elif 0
		   if(offset != 0){
                	this->PutDiagonalSeq(stderr,offset,E1,E2);
			fprintf(stderr,"*** %d vs %d: score = %d; scr=%d ; os=%d; len=%d ***\n",
				i,j,score,scr,offset,len);
		   }
#endif	//=====================================================================
		}
          } 
// free(qop);
        } if(fp) fprintf(fp," %d hits out of %d seqs\n\n",hits,N);
        if(HG && fp){ PutHist(fp,60,HG); } if(HG) NilHist(HG); 

        // 2. Within each set pick the sequence with the highest profile score.
	HG=0;
	if(fp) HG=Histogram("sequence cluster sizes", 1,30,1.0);
        Int4    s,**ListE,*tmpList;
        NEWP(ListE, N+2,Int4); NEW(tmpList, N+2,Int4); NEWP(OffSet,N+5,Int4); 
#if 1	// redo CsqQuery
	if(CsqQuery){
           for(Int4 sq=1; sq <= NumSeqsCMSA(CMA); sq++){
                if(CsqQuery[sq]) NilSeq(CsqQuery[sq]);
           } free(CsqQuery);
        } NEW(CsqQuery,N+5,e_type);  // initially created within Init();
#endif
        for(s=0,i=1; i <= N; i++){
          seti=findDSets(i,sets);
          if(i==seti){  // is this the canonical sequence?
             Int4 sq=0; sq++; tmpList[sq] = i;
             for(j=1; j <= N; j++){
                if(j!=i){
                  setj=findDSets(j,sets);
                  if(seti == setj){ sq++; tmpList[sq] = j; }
                }
             } s++;
             NEW(ListE[s],sq+3,Int4);
             ListE[s][0]=sq;
	     if(HG) IncdHist(sq,HG);
             for(j=1; j <= sq; j++) ListE[s][j]=tmpList[j];
#if 0	//==== DEBUG: put out a full alignment of the cluster using sap_typ ====
	     // j=ListE[s][1]; 
	     // if(ListE[s][1] !=8342) continue;
	     // if(s==20 || s==21 || s==481 || s==580 || s==599)
	     // if(s==599)
	     // if(s==3)
	     // if(sq == 10)
	     // for 580 look at sq #10648
	     // if(j == 142 || j==134 || j==144 || j==150) 
	     // if(j == 63 || j == 64) 
	     {
		// if(sq > 1) this->PutAlnToLongest(stderr, ListE, s);
		if(sq > 1) this->MultiAlnToLongest(stderr, ListE, s);
	     }
#else
		// if(sq > 1) OffSet[s]=this->MultiAlnToLongest(stderr,ListE[s]);
		if(sq > 1) OffSet[s]=this->MultiAlnToLongest(0,ListE[s]);
#if 0	// Error checking...
	     if(sq > 1) for(Int4 x=1; x <= ListE[s][0]; x++){
		e_type xE=TrueSeqCMSA(ListE[s][x],CMA);
		fprintf(stderr,"%d: OffSet[%d][%d] = %d; adj=%d; net=%d; len=%d\n",
			ListE[s][x],s,x,OffSet[s][x],OffSetSeq(xE),OffSet[s][x]+OffSetSeq(xE),LenSeq(xE));
		assert(OffSet[s][x] < LenSeq(xE));
	     }
#endif
#endif	//=====================================================================
#if 1
		if(sq > 1){
		   e_type cE=this->MultiAlnToConsensus(0,OffSet[s],ListE[s]);
		   // e_type cE=this->MultiAlnToConsensus(stderr,OffSet[s],ListE[s]);
		   if(CsqQuery[s]){ NilSeq(CsqQuery[s]); }
	           CsqQuery[s]=cE;
#if 1
		if(sq > 1){
		   if(fp && IsOnTheList(ListE[s]))
		   // if(IsOnTheList(ListE[s]))
		   {
		   	this->PutMultiAlign(stderr,ListE[s],OffSet[s],cE); 
		   }
		} else {
		   if(ListE[s][1] == 1792) fprintf(stderr,"sq 1792 is a singleton\n");
		   if(ListE[s][1] == 1886) fprintf(stderr,"sq 1886 is a singleton\n");
		   if(ListE[s][1] == 1904) fprintf(stderr,"sq 1904 is a singleton\n");
		   if(ListE[s][1] == 1854) fprintf(stderr,"sq 1854 is a singleton\n");
		   if(ListE[s][1] == 3377) fprintf(stderr,"sq 3377 is a singleton\n");
		   if(ListE[s][1] == 3300) fprintf(stderr,"sq 3300 is a singleton\n");
		}
#endif
		}
#endif
          }
        } if(HG && fp) PutHist(fp,60,HG); if(HG) NilHist(HG); 
	Nset=s; fprintf(stderr,"Nset=%d\n",s);
	free(tmpList); NilDSets(sets); free(self_score);
// exit(1);
        // 3. return the sequence sets
        return ListE;
} 

Int4    gmb_typ::Extend(register unsigned char  *p1,register unsigned char  *p2,
		register Int4 len, register char **R, Int4 dropoff)
// compute the straight diagonal score without gaps...
/*******************************************
  Watch out for short, in-and-out frame shifts!!!
  Make sure that all sequences work in a cluster (graph and clique instead of dsets??
  Make sure that subject score is 75% of self score.
 *******************************************/
{
        register Int4	i,s,max;
        s = R[p1[0]][p2[0]];
        for(max=s, i=1; i < len; i++) {
                if((s += R[p1[i]][p2[i]]) > max) { max=s; }
                else if(((max-s) >= dropoff) || s < 0) return 0;
                // else if((s <= (max-dropoff)) || s < 0) return 0;
        } return max;
}

