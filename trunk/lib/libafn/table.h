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

#if !defined(TAB)
#define TAB
#include "probability.h"
#include "stdinc.h"
#include "alphabet.h"
#include "sset.h"

/************************ tables.h - tables data type *********************

		T = <2,2>

		      | column 1 | column 2 |      
		------+----------+----------+ 
		row 1 |  cell 00 | cell 10  |      
		------+----------+----------+ 
		row 2 |  cell 01 | cell 11  |      
		------+----------+----------+  

***************************************************************************/

/******************************* table type ******************************/
typedef struct {
	a_type	A;			/* alphabet for table */
	Int4	i[2];			/* location of partitions */
	sst_typ	bin[2][2];		/* partitioned sets for bins */
	double	**cell;			/* cells with real values */
} tables_type;

typedef tables_type	*t_type;
/******************************* private ********************************/
Int4	expct_tab(double E[2][2], t_type T);
Int4     get_name_tab(char *name, Int4 r, Int4 c, t_type T);

/******************************* Public *********************************/
/*************************** table operations ***************************/
t_type  Table(Int4 i, sst_typ up, sst_typ down, a_type A);   /* create a 2x2 table */
t_type  AddColumnsTab(Int4 i, sst_typ left, sst_typ right, t_type T);
void	PutTable(FILE *fptr,t_type T);	/* output table */
double  ExactTab(t_type T);
double	CnTabEntropy(t_type T);
t_type	NilTab(t_type T);				/* destroy tables */	
t_type	ClearTab(t_type T); 		/* reset tables = 0 */	
Int4    ExpectTable(double E[2][2], t_type T);
double  RelEntropyTable(t_type T);
double  RelEntropyTableMod(double **Obs, double Exp[2][2]);
void    PutTableMod(FILE *fptr,t_type T);
void    PutTableMod(FILE *fptr,char *Title,t_type T);

/************************** macro operations ***************************
***************************************************************************/
#define TABSTINY 1.0e-30
#define TableA(T)       ((T)->A)
#define table_fatal(m)	(fprintf(stderr,"%s\n",m),exit(1),0)
#define	rTab(T)		((T)->i[0])
#define	cTab(T)		((T)->i[1])
#define	valTab(i,j,T)	((T)->cell[(i)][(j)])
#define	AddTab(i,j,V,T)	((T)->cell[(i)][(j)] += (double) (V))
#define	EqTab(i,j,V,T)	((T)->cell[(i)][(j)] = (double) (V))
#define TabBin(r,c,T)  (((r)==0)?(((c)==0)?((T)->bin[0][0]):\
                              ((T)->bin[0][1]))\
                             :(((c)==0)?((T)->bin[1][0]):\
                              ((T)->bin[1][1])) )
#endif

