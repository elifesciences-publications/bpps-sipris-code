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

/* fsm.c - */
#include "fsm.h"

fsm_typ MakeFSM(Int4 T, Int4 numlet, Int4 leng, Int4 **matrix)
{
	fsm_typ B;
	Int4	n,q,q0,q1,q2,s,t,c,d,e,g,h,i,j,k,x,y,z,hits,nlet;
	Int4	blsm62max[40] = { 2418, 1826, 1540, 1264, 1130, 995, 
			759, 627, 458, 304, 222, 164, 105, 83, 70, 61, 
			58, 58, 58, 49, 28, 15, 10, 7, 4, 2, 1, 1, 1, 
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}; 

	if(T >= 40) print_error("T too large; no hits possible.");
	NEW(B,1,fsm_type);
	B->nhits = 0;
	B->nlet = numlet;
	B->len_mtrx= leng;
	B->matrix= matrix;
	/***** create finite automaton *****/
	/*** create transition function d(r,q). ***/
	n = numlet +1; n = 1 + n + n*n;  
	NEWP(B->d,numlet+2,Int4);
	for(c=0; c <= numlet; c++) NEW(B->d[c],n+2,Int4);
	for(q=0, c=0; c <= numlet; c++){
	    q++; B->d[c][0] = q; q0=q;		/** q0 = state "^c" **/
	    for(d=0; d <= numlet; d++){
	    	q++; B->d[d][q0] = q;        	/** q1 = state "cd" **/
	    }
	}
	B->nQ = q;				/** nQ = #states **/
	for(q=0, c=0; c <= numlet; c++){
	    q0 = B->d[c][0];			/** q = state "^c" **/
	    for(d=0; d <= numlet; d++){
	    	q1 = B->d[d][q0];		/** q = state "cd" **/
	    	q2 = B->d[d][0];		/** q2 = state "^d" **/
		for(e=0; e <= numlet; e++){
			q = B->d[e][q2];	/** q2 = state "de" **/
			B->d[e][q1] = q;	/** "cde" = "^de" = "de" **/
		}
	    }
	}
	/*** create acceptance states & lists[q][r][1..] ***/
	// time_t	time1=time(NULL);
	nlet = numlet +1;
	n = B->nQ;
	B->mlist=MkMList(9261*leng,MaxStateGB(nlet));
        // fprintf(stderr,"\nallocation time: %0.1f seconds\n", difftime(time(NULL),time1));
	// time1=time(NULL);

	n = leng - 1;
	NEW(B->tmp,n+3,Int4); 
	g = n; h = n+1;		/*** use a circular matrix **/
	for(hits=0,i=1; i <= n; g=h, h=i, i++){
	   /** look through related 3-words for "matches". **/
	   for(x=0; x <= numlet; x++){
	        q0 = B->d[x][0];
		s = T - matrix[g][x];
		for(y=0; y <= numlet; y++){	
	           q1 = B->d[y][q0];		/** q1 = state 'yx' ***/
		   t = s - matrix[h][y];
		   for(z=0; z <= numlet; z++){	
			if(t <= matrix[i][z]){
				Add2MList(h,StateGB(q1,z), B->mlist);
				hits++;
				/** add to acceptor states **/
			} 
		   }
		}
	   }
	   /*** fprintf(stderr,"\n");  /****/
	}
        // fprintf(stderr,"\nneighborhood time: %0.1f seconds\n",difftime(time(NULL),time1));
	/****/ fprintf(stderr,"hits = %d\n",hits); 
        // fprintf(stderr,"\tcreate time: %0.1f seconds\n", time(NULL)-time1);
	return B;
}

void	NilFSM(fsm_typ B)
{
	Int4	c;

	for(c=0; c <= B->nlet; c++) { free(B->d[c]); }
	NilMList(B->mlist);
	free(B->tmp); free(B->d); 
	free(B);
}

Int4	*Scan4MatchesFSM(register Int4 len, register unsigned char *seq, 
	register Int4 i, Int4 *pos, Int4 *number, register fsm_typ B)
{
	register Int4	q,n;

	q = B->d[seq[i]][0];		/** q = state "^c" **/
	for(i++; i < len; i++){
	   q = B->d[seq[i]][q];		/** q = state "xd" **/
				/** seq[i+1] = next token (mealy model) **/
	   if((n=GetListMList(B->tmp,StateGB(q,seq[i+1]),B->mlist))){ 
		*pos = i; *number = n; return B->tmp;
	   }
	}
	return NULL;
}

Int4	MatcherFSM(Int4 len, register unsigned char *seq, register fsm_typ B)
/************************************************************************
	(list(z) = 1 -> -2 -> 0 -> z)
    pos  =  z,...,-2,-1, 0, 1, 2,...,len  
/************************************************************************/
{
	register Int4	i,d,q,j;
	Int4		nxt,n,s,e,tot,*tmp=B->tmp,item;
	char	res[] =  "XCGASTNDEQKRHWYFVILMP";

	q = B->d[seq[1]][0];		/** q = state "^c" **/
	for(tot=0, i=2; i < len; i++){
	   q = B->d[seq[i]][q];		/** q = state "xd" **/
				/** seq[i+1] = next token (mealy model) **/
	   if((n=GetListMList(tmp,StateGB(q,seq[i+1]),B->mlist))){ 
					/** then q + t signals acceptance **/
		printf("word: %c%c%c: seq. position %d; motif positions: ",
			res[seq[i-1]],res[seq[i]],res[seq[i+1]],i);
		for(j=0;  j < n; j++){
                   s = tmp[j];
		   printf(" %2d", s);
		}
		printf("\n");
		/*** report this on list??? ***/
		tot += n;
	   }
	}
	return tot;
}

/******************************* private *******************************/

void	fsm_error(char *s) { fprintf(stderr,"fsm: %s\n",s); exit(1); }

