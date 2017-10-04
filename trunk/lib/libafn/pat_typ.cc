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

#include "pat_typ.h"

//*************************** pat_typ ****************************************
//*************************** pat_typ ****************************************
//*************************** pat_typ ****************************************
double	pat_typ::SimilarPatterns(FILE *fp,double cutoff,Int4 I,Int4 J, pat_typ *pat, a_type AB)
{
	Int4	i,j,ni,nj,pi,pj,r;
	double	N;
	char	ci,cj,ri,rj;
	sst_typ seti,setj,setx;
	char	IsMatch,str[5000],str1[100],str2[100];
	BooLean	found=FALSE;

	for(i=1; i <= this->NumPttrns; i++) if(this->PttrnCategory[i] == I) break;
	if(i > this->NumPttrns) return 0.0;
	for(j=1; j <= pat->NumPttrns; j++) if(pat->PttrnCategory[j] == J) break;
	if(j > pat->NumPttrns) return 0.0;
        ni=this->NumPttrnRes[i];
	str[0]=0; str2[0]=0; str1[0]=0; N=0.0;
        nj=pat->NumPttrnRes[j];
        for(pi=1; pi <= ni; pi++){
          for(pj=1; pj <= nj; pj++){
             if(this->PosPttrn[i][pi] == pat->PosPttrn[j][pj]){
		seti=0;
		for(ci=0; ri=this->PttrnRes[i][pi][ci]; ci++){
		   r=AlphaCode(ri,AB); setx=SsetLet(r); seti=UnionSset(seti,setx);
		}
		setj=0;
		for(cj=0; rj=pat->PttrnRes[j][pj][cj]; cj++){
			r=AlphaCode(rj,AB); setx=SsetLet(r); setj=UnionSset(setj,setx);
		}
		IsMatch=0;
		if(strcmp(this->PttrnRes[i][pi],pat->PttrnRes[j][pj]) == 0) {
		     	IsMatch='I';	// identical.
		} else  if(SubSset(seti,setj) || SubSset(setj,seti)){
		     	IsMatch='S';	// subset.
		} else if(IntersectSset(seti,setj) != 0){
		     	IsMatch='O';	// Overlap.
		}
		if(IsMatch){
                        sprintf(str1,"%d%s%d vs %d%s%d",
				this->PttrnCategory[i], this->PttrnRes[i][pi], this->PosPttrn[i][pi],
				pat->PttrnCategory[j], pat->PttrnRes[j][pj], pat->PosPttrn[j][pj]);
			switch (IsMatch){
			  case 'I': N += 2.0; sprintf(str2," !\n"); break;
			  case 'S': N += 1.0; sprintf(str2,"\n"); break;
			  case 'O': N += 0.5; sprintf(str2," ?\n"); break;
			  default: print_error("pat_typ::SimilarPatterns( ) error!");
			}
			strcat(str1,str2); strcat(str,str1); 
		}
             }
          }
	} if(fp && N >= cutoff) fprintf(fp,"%s  score = %.1f\n\n",str,N);
	return N;
}

void	pat_typ::SimilarPatterns(FILE *fp, pat_typ *pat, a_type AB)
{
	Int4	i,j,ni,nj,pi,pj,r;
	double	N;
	char	ci,cj,ri,rj;
	sst_typ seti,setj,setx;
	char	IsMatch,str[5000],str1[100],str2[100];

        for(i=1; i < this->NumPttrns; i++){
           ni=this->NumPttrnRes[i];
           for(j=i+1; j <= pat->NumPttrns; j++){
	      str[0]=0; str2[0]=0; str1[0]=0; N=0.0;
              nj=pat->NumPttrnRes[j];
              for(pi=1; pi <= ni; pi++){
                for(pj=1; pj <= nj; pj++){
                   if(this->PosPttrn[i][pi] == pat->PosPttrn[j][pj]){

		     seti=0;
		     for(ci=0; ri=this->PttrnRes[i][pi][ci]; ci++){
			r=AlphaCode(ri,AB); setx=SsetLet(r); seti=UnionSset(seti,setx);
		     }
		     setj=0;
		     for(cj=0; rj=pat->PttrnRes[j][pj][cj]; cj++){
			r=AlphaCode(rj,AB); setx=SsetLet(r); setj=UnionSset(setj,setx);
		     }
		     IsMatch=0;
		     if(strcmp(this->PttrnRes[i][pi],pat->PttrnRes[j][pj]) == 0) {
		     	IsMatch='I';	// identical.
		     } else  if(SubSset(seti,setj) || SubSset(setj,seti)){
		     	IsMatch='S';	// subset.
		     } else if(IntersectSset(seti,setj) != 0){
		     	IsMatch='O';	// Overlap.
		     }
		     if(IsMatch){
                        sprintf(str1,"%d%s%d vs %d%s%d",
				this->PttrnCategory[i], this->PttrnRes[i][pi], this->PosPttrn[i][pi],
				pat->PttrnCategory[j], pat->PttrnRes[j][pj], pat->PosPttrn[j][pj]);
			switch (IsMatch){
			  case 'I': N += 2.0; sprintf(str2," !\n"); break;
			  case 'S': N += 1.0; sprintf(str2,"\n"); break;
			  case 'O': N += 0.5; sprintf(str2," ?\n"); break;
			  default: print_error("pat_typ::SimilarPatterns( ) error!");
			}
			strcat(str1,str2); strcat(str,str1); 
		     }
                   }
                }
              }
	      if(N >= 5.0) fprintf(fp,"%s  score = %.1f\n\n",str,N);
           }
        }
}

void	pat_typ::SimilarPatterns(a_type AB)
{
	Int4	i,j,ni,nj,pi,pj,r;
	double	N;
	char	ci,cj,ri,rj;
	sst_typ seti,setj,setx;
	char	IsMatch,str[5000],str1[100],str2[100];

        for(i=1; i < NumPttrns; i++){
           ni=NumPttrnRes[i];
           for(j=i+1; j <= NumPttrns; j++){
	      str[0]=0; str2[0]=0; str1[0]=0; N=0.0;
              nj=NumPttrnRes[j];
              for(pi=1; pi <= ni; pi++){
                for(pj=1; pj <= nj; pj++){
                   if(PosPttrn[i][pi] == PosPttrn[j][pj]){

		     seti=0;
		     for(ci=0; ri=PttrnRes[i][pi][ci]; ci++){
			r=AlphaCode(ri,AB); setx=SsetLet(r); seti=UnionSset(seti,setx);
		     }
		     setj=0;
		     for(cj=0; rj=PttrnRes[j][pj][cj]; cj++){
			r=AlphaCode(rj,AB); setx=SsetLet(r); setj=UnionSset(setj,setx);
		     }
		     IsMatch=0;
		     if(strcmp(PttrnRes[i][pi],PttrnRes[j][pj]) == 0) {
		     	IsMatch='I';	// identical.
		     } else  if(SubSset(seti,setj) || SubSset(setj,seti)){
		     	IsMatch='S';	// subset.
		     } else if(IntersectSset(seti,setj) != 0){
		     	IsMatch='O';	// Overlap.
		     }
		     if(IsMatch){
                        sprintf(str1,"%d%s%d vs %d%s%d",
				PttrnCategory[i], PttrnRes[i][pi], PosPttrn[i][pi],
				PttrnCategory[j], PttrnRes[j][pj], PosPttrn[j][pj]);
			switch (IsMatch){
			  case 'I': N += 2.0; sprintf(str2," !\n"); break;
			  case 'S': N += 1.0; sprintf(str2,"\n"); break;
			  case 'O': N += 0.5; sprintf(str2," ?\n"); break;
			  default: print_error("pat_typ::SimilarPatterns( ) error!");
			}
			strcat(str1,str2); strcat(str,str1); 
		     }
                   }
                }
              }
	      if(N >= 5.0) fprintf(stdout,"%s  score = %.1f\n\n",str,N);
           }
        }
}


pat_typ *pat_typ::SubPatterns(char *set)
{
	pat_typ *sub_ptrn;
	FILE	*fp=tmpfile();
	Int4	i,j,k;
	BooLean	*IsPrinted,first;
	NEW(IsPrinted, MaxPttrnPos + 3,BooLean);
	assert(set);
	for(j=NumPttrns, i=1; i<= NumPttrns; j--,i++){
	   if(set[j] != '+') continue;
	   // if(NumPttrnRes[i] <= 0) continue;
	   first=TRUE;
	   fprintf(fp,"%d: ",PttrnCategory[i]);
	   for(k=1; k<= NumPttrnRes[i]; k++){
		Int4 pos=PosPttrn[i][k];
		if(IsPrinted[pos]) continue; 
		if(first) fprintf(fp,"%s%d",PttrnRes[i][k],pos);
		else fprintf(fp,",%s%d",PttrnRes[i][k],pos);
		first=FALSE; IsPrinted[pos]=TRUE;
	   } fprintf(fp,"\n");
	} // fprintf(fp,"\n");
	free(IsPrinted);
	rewind(fp); sub_ptrn= new pat_typ(fp); fclose(fp); 
	return sub_ptrn;
}

void	pat_typ::Put(FILE *fp)
{
	for(Int4 i=1; i<= NumPttrns; i++){
	   fprintf(fp,"%d: ",PttrnCategory[i]);
	   for(Int4 j=1; j<= NumPttrnRes[i]; j++){
		if(j==1) fprintf(fp,"%s%d",PttrnRes[i][j],PosPttrn[i][j]);
		else fprintf(fp,",%s%d",PttrnRes[i][j],PosPttrn[i][j]);
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
}

void	pat_typ::Free( )
{
	Int4	i,j;
	for(i=1; i <= NumPttrns; i++){
	    for(j=1; j <= NumPttrnRes[i] ; j++){
		if(PttrnRes[i][j]) free(PttrnRes[i][j]); 
	    } free(PosPttrn[i]); free(PttrnRes[i]); 
	} free(PttrnCategory); free(NumPttrnRes); free(PosPttrn); free(PttrnRes);
}

void	pat_typ::OpenPttrnFile(FILE *fp)
// Obtain several categories of pattern residues from pttrn input file...
{
	NumPttrns=0;
        char    *str,str0[1005],residue_str[50];
	Int4	category,position,pos;

        NEW(PttrnCategory,MAX_NUM_PATTERNS+ 5,Int4);	// 
        NEW(NumPttrnRes,MAX_NUM_PATTERNS+ 5,Int4);	// 
        NEWP(PosPttrn,MAX_NUM_PATTERNS+ 5,Int4); 	// 
        NEWPP(PttrnRes,MAX_NUM_PATTERNS+ 5,char);	// 

	while(fgets(str0,1000,fp)){
	    str=str0;
	    BooLean stop=FALSE;
	    NumPttrns++; pos=0;

	    NEW(PosPttrn[NumPttrns],MAX_NUM_PATTERN_POS+3,Int4);
            NEWP(PttrnRes[NumPttrns],MAX_NUM_PATTERN_POS + 3,char);	// 

            if(NumPttrns >= MAX_NUM_PATTERNS) print_error("Too many pattern categories");
	    do{
// fprintf(stderr,"str[%d]=\"%s\"",pos,str);
	       if(pos == 0){
                  if(sscanf(str,"%d: %[A-Z]%d", &category,residue_str, &position) != 3){
			fprintf(stderr,"str=\"%s\"",str);
			assert(!"fatal error!\n");
			print_error("pttrn_file input error 1");
		  }
		  while(!isupper(str[0])) str++;
		  PttrnCategory[NumPttrns]=category;
	       } else if(sscanf(str,"%[A-Z]%d", residue_str, &position) != 2){
			print_error("pttrn_file input error 2");
	       } pos++;
	       PosPttrn[NumPttrns][pos] = position;
	       PttrnRes[NumPttrns][pos] = AllocString(residue_str);
	       while(str[0] != ','){
		  if(str[0] ==EOF || str[0] == '\n'){ stop=TRUE; break; }
		  assert(isdigit(str[0]) || isupper(str[0])); str++;
	       } if(str[0] == ',') str++; 
            } while(!stop);
	    NumPttrnRes[NumPttrns] = pos;
        }

#if 0
        for(category =1; category <= NumPttrns; category++){
	    fprintf(stderr,"%d: ",PttrnCategory[category]);
	    for(pos=1; pos <= NumPttrnRes[NumPttrns]; pos++){
		fprintf(stderr,"%s%d ",PttrnRes[category][pos],PosPttrn[category][pos]);
	    } fprintf(stderr,"\n");
	}
#endif
        //****************** Determine the Maximum value Pattern position *************
	Int4 x,p;
	MaxPttrnPos=0;
	for(x=NumPttrns; x > 0; x--){
	  for(p=1; p <= NumPttrnRes[x]; p++){
	    Int4 pos=PosPttrn[x][p];
	    MaxPttrnPos=MAXIMUM(Int4,MaxPttrnPos,pos);
	  }
	} //*****************************************************************************
}

