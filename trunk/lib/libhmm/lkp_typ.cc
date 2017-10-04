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

#include "lkp_typ.h"

lkp_typ::lkp_typ(shm_typ **Shm, Int4 Data_len)
// 3-letter word lookup table: assuming a 20 residue alphabet.
{
	fprintf(stderr,"Calling 2 argument lkp_typ() constructor\n");
	data_len = Data_len;
	Int4 m,sum,i,j,n,nh,**nhits,*nhits_i;
	unsigned short *position,*sh,**ent_i,*ent_i_n;
	unsigned char *word;
	Int4 k;
	NEW(word,4,unsigned char);
	unsigned short ***ent;
	NEWPP(ent,8000,unsigned short);
	NEWP(nhits,8000,Int4);
	NEW(start,8001,Int4);
	max_cnt = 0;
	for(total=0,i=0;i<8000;i++){
		start[i] = total+1;
		MEW(nhits[i],data_len+1,Int4);
		NEWP(ent[i],data_len+1,unsigned short);
		ent_i = ent[i];
		nhits_i = nhits[i];
		k=i;
		word[1] = (k/400)+1;
		k = k%400;
		word[2] = (k/20)+1;
		k = k%20;
		word[3] = k+1;
		for(sum=0,n=1;n<=data_len;n++){
			NEW(position,Shm[n]->GetLength()+1,unsigned short);
			nh = Shm[n]->WordHits(word,3,Shm[n]->GetWthreshold(),position);
			nhits_i[n] = nh;
			MEW(ent[i][n],nh+1,unsigned short);
			ent_i_n = ent_i[n];
			for(j=1;j<=nh;j++){
				ent_i_n[j] = position[j];
			}
			free(position);
			sum += nh;
		}
		if(sum > max_cnt) max_cnt = sum;
		total += sum;
	}
	start[i] = total +1;
	NEW(look,total+2,Int4);
	for(m=1,i=0;i<8000;i++) {
		nhits_i = nhits[i]; ent_i = ent[i];
		for(n=1;n<=data_len;n++) {
			ent_i_n = ent_i[n];
			for(j=1;j<=nhits_i[n];j++) {
				sh = (unsigned short *) &(look[m++]);
				sh[0] = n;
				sh[1] = ent_i_n[j];
			}
		}
	}
	for(i=0;i<8000;i++) {
		ent_i=ent[i];
		for(n=1;n<=data_len;n++) {
			free(ent_i[n]);
		}
		free(ent[i]); free(nhits[i]);
	}
	free(ent); free(nhits);
	free(word);
}

lkp_typ::lkp_typ(Int4 Data_len,Int4 Total,Int4 Max_cnt,Int4 *Start,Int4 *Look)
{
	fprintf(stderr,"Calling 5 argument lkp_typ() constructor\n");
	data_len = Data_len;
	total = Total;
	max_cnt = Max_cnt;
	start = Start;
	look = Look;
}

lkp_typ::lkp_typ(FILE *fp)
{
	fprintf(stderr,"Calling 1 argument lkp_typ() constructor\n");
        fread(&data_len,sizeof(Int4),1,fp);
        fread(&total,sizeof(Int4),1,fp);
        fread(&max_cnt,sizeof(Int4),1,fp);
        NEW(start,8001,Int4);
        fread(start,sizeof(Int4),8001,fp);
        NEW(look,total+2,Int4);
        fread(look,sizeof(Int4),total+1,fp);
}

void lkp_typ::Write(FILE *fp)
{
	fwrite(&data_len,sizeof(Int4),1,fp);
	fwrite(&total,sizeof(Int4),1,fp);
	fwrite(&max_cnt,sizeof(Int4),1,fp);
	fwrite(start,sizeof(Int4),8001,fp);
	fwrite(look,sizeof(Int4),total+1,fp);
}

void lkp_typ::Free()
{
	if(start != NULL) free(start);
	if(look != NULL) {
		free(look);
	}
}


