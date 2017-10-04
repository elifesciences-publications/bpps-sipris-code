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

Int4 **gmb_typ::ClusterGPSI_CMSA(char method, double cutoff, Int4 &Nset)
// char	method='E'; double cutoff=0.01;
{
	Int4	N,**seqs,set,s,sq;
	Int4	**List,*tmpList;
	cma_typ cma=SSX->RtnCMA();
        a_type  AB=AlphabetCMSA(cma);

	ss_type data=CopySeqSet(TrueDataCMSA(cma));
	ReNumberSeqSet(data);
	e_type  **ListE=ClusterGPSI(method,cutoff,data,AB,&N);
	Nset=N;
	NEWP(List, NumSeqsCMSA(cma)+2,Int4);
	NEW(tmpList, NumSeqsCMSA(cma)+2,Int4);
	for(s=0,set=1; set <= N; set++){
	   for(sq=1; ListE[set][sq] != 0; sq++){
	   	tmpList[sq] = SeqI(ListE[set][sq]);
	   } s++;
	   NEW(List[s],sq+3,Int4);
	   List[s][0]=sq;
	   fprintf(stderr,"set %d: %d\n",s,sq);
	   for(Int4 j=1; j <= sq; j++) List[s][j]=tmpList[j];
	} assert(Nset==s);
	fprintf(stderr,"cutoff = %g; Nset=%d\n",cutoff,s);
	NilSeqSet(data);
        return List;
}

Int4 **gmb_typ::ClusterCMSA(Int4 percent,Int4 *Nset)
// return a representative set of sequences from cma.
#if 0	//***************************************************************************
#endif	//***************************************************************************
{
	cma_typ cma=SSX->RtnCMA();
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=TrueDataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;
	Int4	cutoff,maxscore;

	assert(percent > 0 && percent < 100);
	// 1. Find the maximum self score for the first seq.
	maxscore=0;
	for(i = 1; i <= N; i++){
	   isq = SeqPtrCMSA(i,cma);
	   for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			score += valAlphaR(isq[si],isq[si],A);
	        }
	   } maxscore+=score;
	} maxscore = maxscore/N;	// use ave self score.
fprintf(stderr,"ave maxscore = %d\n",maxscore);

	// 1. Cluster sequences into sets at the cutoff score.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = 1; i < N; i++){
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
	  for(j=i+1; j <= N; j++){
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			score += valAlphaR(isq[si],jsq[sj],A);
		}
	      } 
	      double d=(double)score/(double)maxscore;
	      Int4   pc=(Int4)floor(100*d+0.5);
// fprintf(stderr,"score/maxscore = %d/%d = %.2f (%d\%)\n",score,maxscore,d,pc);
	      if(percent <= pc) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	Int4	**ListE,*tmpList;
	NEWP(ListE, NumSeqsCMSA(cma)+2,Int4);
	NEW(tmpList, NumSeqsCMSA(cma)+2,Int4);
        for(s=0,i=1; i <= N; i++){
	  seti=findDSets(i,sets);
	  if(i==seti){	// is this the canonical sequence?
	     Int4 sq=0;
	     sq++; tmpList[sq] = i; 
	     for(j=1; j <= N; j++){
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){ sq++; tmpList[sq] = j; }
		}
	     } s++;
	     NEW(ListE[s],sq+3,Int4);
	     ListE[s][0]=sq;
	     for(j=1; j <= sq; j++) ListE[s][j]=tmpList[j];
	  }
	} *Nset=s;
	fprintf(stderr,"Nset=%d\n",s);
	free(tmpList);
	// 3. return the sequence sets 
	NilDSets(sets);
	return ListE;
}

Int4 **gmb_typ::GetRepSetCMSA2(Int4 similarity,Int4 *Nset,cma_typ cma)
// Redundant with e_type **GetRepSetCMSA(FILE *fp, Int4 similarity,Int4 *Nset,cma_typ cma)
// eventually fix this redundancy...
// return a representative set of sequences from cma.
#if 0	//***************************************************************************
#endif	//***************************************************************************
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=TrueDataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;

	assert(similarity > 0 && similarity <= 100);
	// 1. Cluster sequences into sets at the percent identity cutoff.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = 1; i < N; i++){
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
	  for(j=i+1; j <= N; j++){
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			if(valAlphaR(isq[si],jsq[sj],A) > 0) score++; // percent similarity
			// if(isq[si] == jsq[sj]) score++;
		}
	      } score = (Int4) floor(((double)score*100.0/(double)total) +0.5);
	      if(score >= similarity) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	Int4	**ListE,*tmpList;
	NEWP(ListE, NumSeqsCMSA(cma)+2,Int4);
	NEW(tmpList, NumSeqsCMSA(cma)+2,Int4);
        for(s=0,i=1; i <= N; i++){
	  seti=findDSets(i,sets);
	  if(i==seti){	// is this the canonical sequence?
	     Int4 sq=0;
	     sq++; tmpList[sq] = i; 
	     for(j=1; j <= N; j++){
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){ sq++; tmpList[sq] = j; }
		}
	     } s++;
	     NEW(ListE[s],sq+3,Int4);
	     ListE[s][0]=sq;
	     for(j=1; j <= sq; j++) ListE[s][j]=tmpList[j];
	  }
	} *Nset=s;
	fprintf(stderr,"Nset=%d\n",s);
	free(tmpList);
	// 3. return the sequence sets 
	NilDSets(sets);
	return ListE;
}

#if 0
void    gmb_typ::ReplaceSSX(ssx_typ *ssx)
{
        cma_typ CMA=SSX->RtnCMA();
        gss_typ *gss=gssCMSA(CMA);
        cma_typ cma=ssx->RtnCMA();
        // fprintf(stderr,"checking gss in old versus new cma\n");
        assert(gssCMSA(cma)->TrueSqSet() == gssCMSA(CMA)->TrueSqSet());
        assert(NumSeqsCMSA(cma) == NumSeqsCMSA(CMA));
        assert(nBlksCMSA(cma) == nBlksCMSA(CMA));
        for(Int4 b=1; b <= nBlksCMSA(CMA); b++){
            assert(LengthCMSA(b,cma) == LengthCMSA(b,CMA));
        } delete SSX; gss->~gss_typ(); NilCMSA(CMA);    // gss_typ separate from NilCMA!
        SSX=ssx; AB=AlphabetCMSA(cma);
}

cma_typ	gmb_typ::ProgressiveSampling(char *name, Int4 start, Int4 stop,Int4 inc)
{
   Int4		NumClusters,iter=0;
   double	Temperature=300.0;
   double	map,bestmap;
   cma_typ	CMA=0;
   for(Int4 similarity=start; similarity <= stop; similarity+=inc){
       do {
                map=GambitClusters('S',similarity,Temperature);
                if(bestmap < map){
			CMA=SSX->RtnCMA();
                        sprintf(str,"%s.new",name);
                        WriteMtfCMSA(str, CMA, NULL);
                        sprintf(str,"%s.new.cma",name);
                        WriteCMSA(str, CMA);
                        bestmap=map;
                        fprintf(stderr,"########## map %.1f ###########\n",map);
                } else break;
                iter++;
                fprintf(stderr,"********** end (%d\%) iter %d ***********\n",similarity,iter);
       } while(TRUE);
       CMA=SSX->RtnCMA();
       if(NumClusters == NumSeqsCMSA(CMA)) break;
   } CMA=SSX->RtnCMA(); return CMA;
}

BooLean	gmb_typ::SampleSingle0(Int4 sq,double Temperature,double &MAP)
// sample insertions and deletions in sequence sq.
{
	Int4	start,trace_length,score,*newpos,**opos,cluster[5];

	//======== 1. Remove sq from duplicate alignment and sample scoring matrix. =====
	ssx_typ *ssx= new ssx_typ(SSX); cluster[0]=1; cluster[1]=sq; cluster[2]=0;
	cma_typ cma=ssx->RtnCMA(); ssx->RmFromAlign(cluster); 
	smx=ssx->GetSMXs(Temperature);

	//========== 2. Sample a gapped alignment for sequence s. ===========
	gss_typ *gss=gssCMSA(cma);
	e_type	sbjE = gss->TrueSeq(sq); 	
	char	*operation=jlh->GapAlnTrace(sbjE,nBlksCMSA(cma),smx,&start,&score);
	trace_length=strlen(operation);
	// this->PutGappedAln(stderr,sbjE,operation,start,cma);

	//========== 3. Create a gapped alignment. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(cma)+3,Int4);  
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,sbjE,newpos);
	free_smx(); free(operation); 

	//========== 4. If sq alignment unchanged then return. =========
	if(gss->Identical(sq,*gsq)){ 
	    this->DeleteSSX(ssx); free(newpos); delete []gsq; return FALSE; 
	}

	//========== 5. Otherwise sample one of the two alignments. =========
	cma_typ CMA=SSX->RtnCMA(); 
	ReplaceCMSA(sq,gsq,cma); ssx->AddToAlign(sq,newpos); free(newpos); 
	double	oldmap=SSX->RelMap( ) + jlh->IndelPenalty(stderr,'S',CMA);
        double	newmap=ssx->RelMap( ) + jlh->IndelPenalty(stderr,'S',cma);
	if(comp_map_cmsa(newmap, oldmap, 300./Temperature)){
		MAP=newmap; this->ReplaceSSX(ssx); return TRUE;
	} else { MAP=oldmap; this->DeleteSSX(ssx); return FALSE; }
}

BooLean	gmb_typ::SampleCluster0(Int4 *cluster, Int4 *offset,e_type csqE,double Temperature,double &MAP)
// sample a gapped alignment for a cluster of sequences.
{
	Int4	start,strt,trace_length,score,sq,ii,*newpos;
	BooLean	Modified=FALSE,found=FALSE;

	//======= 1. Remove seq cluster from a new alignment and sample matrix.  ==========
	ssx_typ *ssx= new ssx_typ(SSX); assert(cluster[0] > 1);
	cma_typ	cma=ssx->RtnCMA(); ssx->RmFromAlign(cluster);
	gsq_typ **gsq; NEW(gsq,cluster[0]+2,gsq_typ *);
	smx=ssx->GetSMXs(Temperature);

	//========== 2. Sample a gapped alignment for master sequence csqE. ===========
	char *operation=jlh->GapAlnTrace(csqE,nBlksCMSA(cma),smx,&start,&score);
	trace_length=strlen(operation);
	assert(operation[1]=='M' || operation[1]=='D');

	// ========== 3. Align sequences based on the master. ===========
	gss_typ	*gss=gssCMSA(cma);  
	for(ii=1; (sq=cluster[ii]) != 0; ii++) {
	   e_type sbjE=gss->TrueSeq(sq); 	
	   char *op=ConvertOperation(operation,sbjE,start,offset[ii]);
	   strt=start-offset[ii]; trace_length=strlen(op);
	   if(strt < 1) strt=1;
	   // if((offset[ii] != 0 && op[1]=='D') || op[trace_length-2]=='d')
	   if(sq==3304 || sq==3444 || sq==3300 || sq==1854 || sq==3310)
	   // if(0)	// DEBUG...
	   { this->PutMasterSlaveAln(sq,csqE,start,operation,sbjE, strt,op,offset[ii],cma); found=TRUE; }
	   NEW(newpos,nBlksCMSA(cma)+2,Int4);  // for new sites.
	   gsq[ii] = new gsq_typ[1]; // NEEDS TO BE AN ARRAY FOR gss_typ!!!
	   gsq[ii]->initialize(gss->LeftFlank(),gss->RightFlank(),op,trace_length,strt,sbjE,newpos);
	   if(gss->Identical(sq,*gsq[ii])){ delete []gsq[ii]; }
	   else { Modified=TRUE; ReplaceCMSA(sq,gsq[ii],cma); }
	   ssx->AddToAlign(sq,newpos); free(newpos); free(op); 
	} free(gsq); free(operation); free_smx(); 

	//========== 4. Sample one of the two alignments. =================
	cma_typ CMA=SSX->RtnCMA(); 
	double	newmap,oldmap=SSX->RelMap() + jlh->IndelPenalty(stderr,'S',CMA);
	if(Modified){ newmap=ssx->RelMap() + jlh->IndelPenalty(stderr,'S',cma); }
	if(Modified && comp_map_cmsa(newmap, oldmap, 300./Temperature)){
	   MAP=newmap; this->ReplaceSSX(ssx); return TRUE;
	} else {	// reject the new cmsa file.
	   // if(found){ this->PutClusterAln("junk_debug",cluster,cma); } 
	   MAP=oldmap; this->DeleteSSX(ssx); return FALSE; 
	}
}
#endif
