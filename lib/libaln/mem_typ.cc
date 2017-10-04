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

#include "mem_typ.h"

mem_typ::mem_typ(Int4 NumActions,char typ){ init(5,25,NumActions,0.25,typ); }

mem_typ::mem_typ(Int4 NumActions, char typ, unsigned short LenShort,
		unsigned short LenLong, double percentPs)
{ init(LenShort,LenLong,NumActions,percentPs,typ); }

void    mem_typ::init(unsigned short LenShort,unsigned short LenLong,
	Int4 NumActions, double percentPs,char typ)
{
	assert(percentPs >= 0.0 && percentPs <= 100.0);
	assert(LenShort < LenLong);

	ActionType=typ;
	N = NumActions;
	shorterm = new set_typ [N+1];
	longterm = new set_typ [N+1];
	eventL	= new unsigned short [N+1];
	eventS	= new unsigned short [N+1];
	longLen=LenLong; shortLen=LenShort;
	Ps=percentPs;
	
	for(Int4 i=1; i<=N; i++){
	   // FillSet == assumes that each action has always succeeded.
	   shorterm[i] = MakeSet(shortLen); FillSet(shorterm[i]);
	   longterm[i] = MakeSet(longLen); FillSet(longterm[i]);
	   eventL[i] = eventS[i] = 0;
	}
}

void	mem_typ::Free( )
{
	for(Int4 i=1; i<=N; i++){ NilSet(shorterm[i]); NilSet(longterm[i]); }
	delete [] eventS; delete [] eventL;
	delete [] shorterm; delete [] longterm;
}

BooLean	mem_typ::Full(Int4 i)
// is the memory empty (i.e., no actions stored).
{
	assert(i > 0 && i <= longLen);
	if(CardSet(longterm[i])==longLen) return TRUE; else return FALSE;
}

BooLean	mem_typ::Empty(Int4 i)
// Is the memory empty (i.e., no actions stored).
{
	assert(i > 0 && i <= longLen);
	if(CardSet(longterm[i])==0) return TRUE; else return FALSE;
}

void	mem_typ::Clear(Int4 i)
{
	assert(i > 0 && i <= longLen);
	ClearSet(shorterm[i]); ClearSet(longterm[i]); eventL[i]=eventS[i]=0;
}

void	mem_typ::Failure(unsigned short action)
// action = the observation.
{
	assert(action > 0 && action <= N);
	eventL[action]++; eventS[action]++; 
	if(eventL[action] >= longLen) eventL[action]=0;
	DeleteSet(eventL[action],longterm[action]);
	if(eventS[action] >= shortLen) eventS[action]=0;
	DeleteSet(eventS[action],shorterm[action]);
}

void	mem_typ::Success(unsigned short action)
// Record that ith feature has changed (sqalns, blocks, columns, gaps, sites?
{
	assert(action > 0 && action <= N);
	eventL[action]++; eventS[action]++; 
	if(eventL[action] >= longLen) eventL[action]=0;
	AddSet(eventL[action],longterm[action]);
	if(eventS[action] >= shortLen) eventS[action]=0;
	AddSet(eventS[action],shorterm[action]);
}

void	mem_typ::Put(FILE *fptr, unsigned short action)
{
	Int4	i,j;
	assert(action > 0 && action <= N);
	switch(ActionType){
	  case 'S': fprintf(fptr,"realign sequence "); break;
	  case 'B': fprintf(fptr,"move block "); break;
	  case 'G': fprintf(fptr,"realign gapped sequence "); break;
	  case 'C': fprintf(fptr,"move column "); break;
	  default: fprintf(fptr,"action "); break;
	} fprintf(fptr,"%d (p=%.3f):\nS: ",action,Chance(action));
	for(j=eventS[action],i=0; i < shortLen; i++){
		if(MemberSet(j,shorterm[action])) fprintf(fptr,"*");
		else fprintf(fptr,".");
		if(i%60==59 && (i+1) < shortLen) fprintf(fptr,"\n   ");
		if(j<=0) j=shortLen; j--; 
	} fprintf(fptr,"\nL: ");
	for(j=eventL[action],i=0; i < longLen; i++){
		if(MemberSet(j,longterm[action])) fprintf(fptr,"*");
		else fprintf(fptr,".");
		if(i%60==59 && (i+1) < longLen) fprintf(fptr,"\n   ");
		if(j<=0) j=longLen; j--; 
	} fprintf(fptr,"\n");
}

double	mem_typ::Chance(unsigned short action)
// Probability (chance) of attempting a particular action.
{
	assert(action > 0 && action <= N);
	double	L,S,D;
	Int4	l,s;

	L = (double) CardSet(longterm[action])/(double) longLen;
	S = (double) CardSet(shorterm[action])/(double) shortLen;
	// D = S/L;	// derivative: how fast are things changing?
	// D = 1.0 --> no change.  D > 1.0 --> change for better.
	// D < 1.0 --> change for worse.
	return (L + S + 2.0*Ps)/(2.0*(Ps + 1));
}


