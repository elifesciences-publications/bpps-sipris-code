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

#include "esc_typ.h"

//*************************** esc_typ ****************************************
//*************************** esc_typ ****************************************
//*************************** esc_typ ****************************************

void    esc_typ::Free( )
{
	Int4	s,i,j;
	for(s=1; s <= NumPDB_Sets; s++){
	   free(PDB_SetI[s]); free(PDB_SetC[s]); NilSeq(FullSeq[s]);
	} free(FullSeq); free(BestInSet); free(PDB_SetI); free(PDB_SetC); free(NumPDB_Set);
}

void    esc_typ::CreateDsets( )
// 
{
   Int4		C,D,I,J,i,j,N,setI,setJ,os,osI,osJ;
   ds_type	sets;
   e_type	seqI,seqJ;
   Int4		**MapSeqID;     // MapSeqID[I][C] == i; use this for disjoint set determination.
   Int4		*MapSeqChn;     // MapSeqChn[i] == C;   use this to find chain for i.
   Int4		*MapSeqI;       // MapSeqI[i] == I;     use this to find I.

   // 1. Create disjoint sets.
   NEWP(MapSeqID,NumFiles+3,Int4);
   FILE *ofp=0;
   // FILE *ofp=open_file("test",".seq","w");
   for(N=0,I=1; I <=NumFiles; I++){
   	NEW(MapSeqID[I],nChainsPDB(pdb[I])+3,Int4);
        for(C=1; C <= nChainsPDB(pdb[I]); C++){
	    if(pdbSeq[I][C] && LenSeq(pdbSeq[I][C]) >= MinSeqOverlap){
		N++; MapSeqID[I][C]=N;
		if(ofp) PutSeq(ofp,pdbSeq[I][C],AB);
	    }
	}
   } sets = DSets(N);	// puts every sequence in its own set.
   if(ofp) fclose(ofp);

   // 2. Find mapping from each set back to pdb structure.
   NEW(MapSeqChn,N+3,Int4); NEW(MapSeqI,N+3,Int4);
   for(I=1; I <= NumFiles; I++){
     for(C=1; C <= nChainsPDB(pdb[I]); C++){
	i=MapSeqID[I][C];
	if(i){ MapSeqChn[i]=C; MapSeqI[i]=I; }
	assert(i <= N);
     }
   } 

   // 3. Cluster pdb files into identical sets.
   for(I=1; I <= NumFiles; I++){
     for(C=1; C <= nChainsPDB(pdb[I]); C++){
        seqI=pdbSeq[I][C];	// sequences with > 10% 'X' residues are deleted...
	Int4 x=0;
	if(seqI) x = LenSeq(seqI);
// fprintf(stderr,"pdb_seq[%d][%d]: len = %d; min_overlap=%d\n",I,C,x,MinSeqOverlap);

	if(seqI == 0 || LenSeq(seqI) < MinSeqOverlap) continue; 
// PutSeq(stderr,seqI,AB);
	i=MapSeqID[I][C]; assert(i > 0 && i <= N); setI = findDSets(i,sets);
	osI =OffSetSeq(seqI);
	for(J=I; J <=NumFiles; J++){
           for(D=1; D <= nChainsPDB(pdb[J]); D++){
	     if(J==I && C==D) continue;
             seqJ=pdbSeq[J][D];
	     if(seqJ == 0 || LenSeq(seqJ) < MinSeqOverlap) continue; 
	     j=MapSeqID[J][D]; assert(j > 0 && j <= N); 
	     osJ =OffSetSeq(seqJ);
	     // char rtn=IsSameSeqFast(seqI,seqJ,&os,MinSeqOverlap); // ignoring 'X' residues...
	     char rtn=IsSameSeqFastX(seqI,seqJ,&os,MinSeqOverlap); // ignoring 'X' residues...
	     if(rtn){	// then these belong in the same set.
	       Int4 d_os;
	       switch (rtn){
		case 1: d_os=osI-osJ; break; 	// seqI N-terminus starts within seqJ.
		case 2: d_os=osJ-osI; break; 	// seqJ N-terminus starts within seqI.
	        default: assert(rtn < 3 && rtn > 0);
	       } os = os-d_os;	// adjust for inherent offset in sequence.
	       if(os == 0){	// then these belong in the same set.
			j=MapSeqID[J][D]; assert(j > 0 && j <= N); setJ = findDSets(j,sets);
			if(setI != setJ){
				// make sure that seqJ matches all other seqs within setI.
				BooLean okay=TRUE;	// may differ due to distinct extensions.
				for(Int4 K=1; K <= NumFiles && okay; K++){
           			   for(Int4 E=1; E <= nChainsPDB(pdb[K]) && okay; E++){
					if(I==K && C==E) continue;
					if(J==K && D==E) continue;
             				e_type seqK=pdbSeq[K][E];
	     				if(seqK == 0 || LenSeq(seqK) < MinSeqOverlap) continue; 
					Int4 k=MapSeqID[K][E]; assert(k > 0 && k <= N);
					Int4 setK = findDSets(k,sets);
					if(setK == setI){	// seqK is in setI; compatible with setJ?
	     					// rtn=IsSameSeqFast(seqJ,seqK,&os,MinSeqOverlap);
	     					rtn=IsSameSeqFastX(seqJ,seqK,&os,MinSeqOverlap);
						if(rtn == 0) okay=FALSE;
// else fprintf(stderr," a match\n");
					}
					if(setK == setJ){	// seqK is in setJ; compatible with setI?
	     					// rtn=IsSameSeqFast(seqI,seqK,&os,MinSeqOverlap);
	     					rtn=IsSameSeqFastX(seqI,seqK,&os,MinSeqOverlap);
						if(rtn == 0) okay=FALSE;
// else fprintf(stderr," a match\n");
					}
				   }
				}
				if(okay) setI = linkDSets(setI,setJ,sets); 
// else fprintf(stderr," not okay...don't link sets %d and %d\n",setI,setJ);
			}
	       } else {		// check to see whether these are offset consistently...
#if 0	// for debugging...
		fprintf(stderr,"============= Inconsistent numbering: rtn=%d ==============\n",rtn);
		fprintf(stderr,"seqI(%d) == seqJ(%d); os = %d; osI=%d; osJ=%d\n",i,j,os+d_os,osI,osJ);
		// PutSeq(stderr,seqI,AB); PutSeq(stderr,seqJ,AB);
		// AlnSeqSW(stderr,11,2,seqI,seqJ,AB);
		if(rtn==2) PutDiagonalSeq(stderr, os+d_os, seqI,seqJ,AB);
		else PutDiagonalSeq(stderr, os+d_os, seqJ,seqI,AB);
#endif
	       }
	     }
	  }
       }
     }
   } // PutDSets(stderr,sets);
   // PutDSets(stdout,sets);
   //********************* Create PDB sets **************************
   NumPDB_Sets=NSetsDSets(sets);
   fprintf(stderr,"Initial NumPDB_Sets = %d\n",NumPDB_Sets); // exit(1);
   NEW(BestInSet,NumPDB_Sets+3,Int4);
   NEW(NumPDB_Set,NumPDB_Sets +3,Int4); NEWP(PDB_SetI,NumPDB_Sets +3,Int4);
   NEWP(PDB_SetC,NumPDB_Sets +3,Int4);
   NEW(FullSeq,NumPDB_Sets+3,e_type);
   Int4	num_pdb_sets=0,R,BestR;
   for(i=1; i <= N; i++){
	setI = findDSets(i,sets);
	if(i != setI) continue; 	// i is not the canonical sequence...
	if(MapSeqI[i] == 0) continue;	// i is canonical but was found to be irrelevant.
	//*************** Find sequences in the same set. ****************
	num_pdb_sets++;
	assert(num_pdb_sets <= NumPDB_Sets);
	// fprintf(stderr,"SetI = %d: \n",setI);
	//************ Print out vsi file for setI *************
	Int4  NumList=0;
	NEW(PDB_SetC[num_pdb_sets],N+3,Int4);
	NEW(PDB_SetI[num_pdb_sets],N+3,Int4);
        for(BestR=0,NumList=0,j=1; j <= N; j++){
            setJ = findDSets(j,sets);
            if(setJ == setI){		// sequence j is in setI.
		C = MapSeqChn[j]; I = MapSeqI[j]; assert(I > 0);
		Int4 R=NonXResSeq(pdbSeq[I][C],AB);
		NumList++;
		if(R > BestR){ BestInSet[num_pdb_sets]=NumList; BestR=R; }
		assert(NumList <= N);
		PDB_SetC[num_pdb_sets][NumList]=C;
		PDB_SetI[num_pdb_sets][NumList]=I;
            }
	} assert(NumList > 0);
	NumPDB_Set[num_pdb_sets] = NumList; 

	//****************** Find the start and end of the full sequence
	Int4 End=0,S=num_pdb_sets;	// latest set.
	for(j=1; PDB_SetI[S][j]; j++){
		I=PDB_SetI[S][j]; C=PDB_SetC[S][j]; e_type Sq=pdbSeq[I][C];
		Int4 len=LenSeq(Sq) + OffSetSeq(Sq);
		End=MAXIMUM(Int4,len,End);
	}
	//******************* Obtain a full length sequence spanning all fragments.
	unsigned char *seq; NEW(seq,End +5,unsigned char);
	char defline[200]; 
	for(j=1; PDB_SetI[S][j]; j++){
		I=PDB_SetI[S][j]; C=PDB_SetC[S][j]; e_type Sq=pdbSeq[I][C];
		if(j==1) StrSeqID(defline,30,Sq);
		Int4 os=OffSetSeq(Sq),s,sx;
		for(s=1,sx=s+os; s<= LenSeq(Sq); s++,sx++){ 
			assert(sx <= End);
			unsigned char r=ResSeq(s,Sq);
			if(r == 0) continue;
			if(seq[sx]) assert(seq[sx] == r); else seq[sx]=r;
		}
	} 
#if 1
	FullSeq[S]=MkSeq(defline,End,seq);
#else	// Now do this on the fly within IsSameSeqFastX() routine!!
	e_type tmpE=MkSeq(defline,End,seq);
	Int4 start=0,end=0;
	for(start=1; ResSeq(start,tmpE)==0; start++) ;
	for(end=LenSeq(tmpE); ResSeq(end,tmpE)==0; end--) ;
	if(start > 1 || end < LenSeq(tmpE)){
		FullSeq[S]=MkSubSeq(start,end,tmpE); NilSeq(tmpE);
	} else FullSeq[S]=tmpE;
#endif
#if 0	// for debugging...
	PutSeq(stderr,FullSeq[S],AB);
	PutSeq(stderr,pdbSeq[I][C],AB);
	for(j=1; PDB_SetI[S][j]; j++){
		I=PDB_SetI[S][j]; C=PDB_SetC[S][j]; e_type Sq=pdbSeq[I][C];
		AlnSeqSW(stderr,11,2,Sq,FullSeq[S],AB);
	}
#elif 0
	for(j=1; PDB_SetI[S][j]; j++){
		I=PDB_SetI[S][j]; C=PDB_SetC[S][j]; e_type Sq=pdbSeq[I][C];
		AlnSeqSW(stderr,11,2,Sq,FullE,AB);
	}
        //*************** output seq files ********************
	char str[100]; sprintf(str,"%d.seq",S);
        FILE *fp=open_file("pdbEqSet",str,"w");
        PutSeq(fp,FullSeq[S],AB);
        for(j=1; PDB_SetI[S][j]; j++){
                I=PDB_SetI[S][j]; C=PDB_SetC[S][j];
                if(!pdbSeq[I][C]) continue;
                PutSeq(fp,pdbSeq[I][C],AB);
        } fclose(fp);
#endif
	free(seq);
   } NumPDB_Sets=num_pdb_sets;	// adjust for irrelevant sequences...
   fprintf(stderr,"NumPDB_Sets = %d\n",NumPDB_Sets); // exit(1);
   NilDSets(sets);
   for(N=0,I=1; I <=NumFiles; I++) free(MapSeqID[I]);
   free(MapSeqID); free(MapSeqChn); free(MapSeqI);
}

