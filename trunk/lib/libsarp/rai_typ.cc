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

#include "rai_typ.h"

void	rhb_typ::Init(atm_typ hydrogen, atm_typ donor, atm_typ accept,
              double D_DA, double D_HA, double Angle_DHA)
{
	H = hydrogen;   Donor=donor; Accept=accept;
        d_HA=D_HA; d_DA=D_DA; angle_DHA=Angle_DHA;
	next=0;
}

void	rhb_typ::Put(FILE *fp,char color, unsigned short file_id, a_type AB, a_type nAB)
{
	Int4	a,a2,r,r2,rA,rB;
	atm_typ	A,B;
	Int4	isaccept=0;
	char	rcA;

	assert(fp);
	A=Donor;B=Accept; 
	BooLean NeitherCarbons=TRUE,na;
	if(CarbonAtom(A) || CarbonAtom(B)) NeitherCarbons=FALSE;
	for( ; isaccept < 2; A=Accept,B=Donor,isaccept++){
		rcA = GetResidueAtom(A, AB, nAB, &na);
		rA=ResAtom(A); rB=ResAtom(B);
#if 0	// afn: 8_26_2015
		if(rA < 0) continue;	// leads to syntax error...
		if(rB < 0) continue;	// leads to syntax error...
#endif
		char str[50],*str_ptr = AtomName(A);
		if(IsWaterAtom(A) && A==Donor) str_ptr = AtomName(H);
		if(isspace(str_ptr[0])) str_ptr++;
		Int4 z,zz;
		for(z=0; str_ptr[z] && !isspace(str_ptr[z]); z++){
			if(isalpha(str_ptr[z])) str[z]=tolower(str_ptr[z]);
			else str[z]=str_ptr[z];
		} str[z]=0; assert(z < 50);
		if(IsWaterAtom(A)){
		 	// if(ResAtom(Accept)==567) PutAtom(stderr,Accept); // debug...
			if(A==Donor){
			     fprintf(fp,"#HOH%d_o-%s.X\t// to ",rA,str);
			} else {
			     fprintf(fp,"#HOH%d_o.X\t// to ",rA);
			}
		} else {
#if 1	// add hydrogen too...
			if(!isaccept){	// i.e. donor atom...
			    str_ptr = AtomName(H);
			    if(isspace(str_ptr[0])) str_ptr++;
			    for(str[z]='-',z++,zz=0; str_ptr[zz] && !isspace(str_ptr[zz]); zz++,z++){
				if(isalpha(str_ptr[zz])) str[z]=tolower(str_ptr[zz]);
				else str[z]=str_ptr[zz];
			    } str[z]=0; assert(z < 50);
			}
#endif
			if(strcmp("o",str)==0) strcpy(str,"c-o");
			else if(str[0]=='c' && strstr(str,"-") == 0){ fprintf(fp,"// "); }
			if(NeitherCarbons){
				if(file_id > 0) fprintf(fp," %d",file_id); else fprintf(fp," ");
			} // else fprintf(fp,"#"); 
			else fprintf(fp," ");
			if(IsHeteroAtom(A)){ fprintf(fp,"!"); }
			if(AtomChain(A) != ' '){
			  if(rcA && !na) fprintf(fp,"%c%d%c",AlphaChar(rcA,AB),rA,AtomChain(A));
			  else fprintf(fp,"%s%d%c",AtomResName(A),rA,AtomChain(A));
			} else {
			  if(rcA && !na) fprintf(fp,"%c%d", AlphaChar(rcA,AB),rA);
			  else fprintf(fp,"%s%dX", AtomResName(A),rA);
			} fprintf(fp,"_%s.X\t// to ",str);
		}
		if(AtomChain(B) != ' '){
		    fprintf(fp, "%s %s%d%c (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),
					rB,AtomChain(B),d_DA,d_HA,angle_DHA);
		} else fprintf(fp, "%s %s%d (DA=%.2f A; HA=%.2f; DHA=%.1f degrees).\n",
					AtomName(B),AtomResName(B),rB,d_DA,d_HA,angle_DHA);
	} 
	
}


