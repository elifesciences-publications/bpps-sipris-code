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

#include "rasmol.h"

void    ScriptRasMol(char *filename, char su, Int4 *res, float *score,
	unsigned char *seq)
/** create a RasMol script protein chain 1 >= stdev **/
{
	FILE	*fptr;
	Int4	i,s,red=255,green=0,blue=0,width=200;
	char	r;
	char	*resname;
	const char	*residue[] = { "nil", /** X CGASTNDEQK RHWYFVILMP **/
		"cys","gly","ala","ser","thr","asn","asp","glu","gln","lys",
		"arg","his","trp","tyr","phe","val","ile","leu","met","pro"};

	fptr = open_file(filename,".ras","w");
	fprintf(fptr,"load pdb %s.pdb\n",filename);
	fprintf(fptr,"set background black\n");
	fprintf(fptr,"color [100,100,100]\nwireframe off\n");
	// fprintf(fptr,"color [90,163,218]\nwireframe off\n"); 
	if(su == ' ') {
		fprintf(fptr,"select protein\nset strands 1\nstrands\n");
		fprintf(fptr,"strands\n");
	} else {
		fprintf(fptr,"select *%c\nset strands 1\nstrands\n",su);
		fprintf(fptr,"select not *%c and protein\n",su);
		fprintf(fptr,"backbone 200\ncolor red\n");
	}
	for(i=1; (s=res[i]) > 0; i++){ 
	   get_zcolor_rasmol(score[i], &red, &green, &blue);
	   if(su == ' ') fprintf(fptr,"select %d\n",s);
	   else fprintf(fptr,"select %d%c\n",s,su);
	   fprintf(fptr,"strands off\n");
	   fprintf(fptr,"color [%d,%d,%d]\n",red,green,blue);
	   fprintf(fptr,"cartoon 300\n");
	   /*** fprintf(fptr,"ribbons %d\n",width); /** for PCs **/
	   if(score[i] >= 2.0){
	     resname = (char *) residue[seq[s]];
	     if(su!=' '){
		fprintf(fptr,"select %s%d%c.ca,(%d%c and sidechain)\n",
			resname,s,su,s,su);
	        fprintf(fptr,"wireframe 80\n");
		fprintf(fptr,"select %d%c and sidechain\n", s,su);
	        fprintf(fptr,"color orange\n");
	     } else {   
		fprintf(fptr,"select %s%d.ca,(%d and sidechain)\n",
			resname,s,s);
	        fprintf(fptr,"wireframe 80\n");
		fprintf(fptr,"select %d and sidechain\n", s);
	        fprintf(fptr,"color orange\n");
	     }
	   }
	}
	fprintf(fptr,"select not protein and not water\n");
	fprintf(fptr,"spacefill true\ncolor red\n");
	fprintf(fptr,"select hetero and not water\n");
	fprintf(fptr,"spacefill true\ncolor magenta\n");
	fprintf(fptr,"select nucleic\nwireframe 40\ncolor red\n");
	if(su != ' ') fprintf(fptr,"select *%c\n",su);
	else fprintf(fptr,"select protein\n");
	fclose(fptr);
}

void    ScriptRasMolFile(FILE *fptr, char *filename, char su, Int4 *res, 
	float *score)
/** create a RasMol script protein chain 1 >= stdev **/
{
	Int4	i,s,red=255,green=0,blue=0,width=200;
	/** XCGASTNDEQKRHWYFVILMP **/

	fprintf(fptr,"load pdb %s.pdb\n",filename);
	fprintf(fptr,"color [90,163,218]\nwireframe off\n");
	fprintf(fptr,"set background black\n");
	if(su == ' ') fprintf(fptr,"select protein\nset strands 1\nstrands\n");
	else fprintf(fptr,"select *%c\nset strands 1\nstrands\n",su);
	if(su!=' '){
	   fprintf(fptr,"backbone 50\nselect not *%c and protein\n",su);
	   fprintf(fptr,"backbone 200\ncolor red\n");
	} else fprintf(fptr,"strands\n");
	for(i=1; (s=res[i]) > 0; i++){ 
	   get_zcolor_rasmol(score[i], &red, &green, &blue);
	   if(score[i] >= 1.0){
	     if(su!=' '){
		fprintf(fptr,"select %d%c\n",s,su);
	     } else {   
		fprintf(fptr,"select %d\n",s);
	     }
	   }
	   if(su == ' ') fprintf(fptr,"select %d\n",s);
	   else fprintf(fptr,"select %d%c\n",s,su);
	   fprintf(fptr,"color [%d,%d,%d]\n",red,green,blue);
	   fprintf(fptr,"strands off\n");
	   fprintf(fptr,"cartoon 300\n");
	   /*** fprintf(fptr,"ribbons %d\n",width); /****/
	}
	fprintf(fptr,"select not protein and not water\n");
	fprintf(fptr,"spacefill true\ncolor red\n");
	fprintf(fptr,"select hetero and not water\n");
	fprintf(fptr,"spacefill true\ncolor magenta\n");
	fprintf(fptr,"select nucleic\nwireframe 40\ncolor red\n");
	if(su != ' ') fprintf(fptr,"select *%c\n",su);
	else fprintf(fptr,"select protein\n");
}

void    get_zcolor_rasmol(double Z, Int4 *red, Int4 *green, Int4 *blue)
{
        Int4    z;

	Z *= 4;
        if(Z <= 1.0) Z = 0.0;
        z = (Int4) floor(Z);
        switch(z) {
                case 0: *red = 105;  *green= 150; *blue = 200; break;
                case 1: *red = 120;  *green= 137; *blue = 182; break;
                case 2: *red = 135;  *green= 124; *blue = 164; break;
                case 3: *red = 150;  *green= 111; *blue = 146; break;
                case 4: *red = 165;  *green= 98; *blue = 128; break;
                case 5: *red = 180;  *green= 85; *blue = 110; break;
                case 6: *red = 195;  *green= 72; *blue = 92; break;
                case 7: *red = 210;  *green= 59; *blue = 74; break;
                case 8: *red = 225;  *green= 46; *blue = 56; break;
                case 9: *red = 240;  *green= 33; *blue = 38; break;
                default: *red = 255;  *green= 20; *blue = 20; break;
        }
}

