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

/* defines.h - generic codes and constants for afn biosequence programs. */
#if !defined (_STDINC_)
#define _STDINC_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

extern double lngamma(double x);

#ifdef __cplusplus
}
#endif

/* VALUES */
#define NIL    		-1

#ifndef FALSE 
#define FALSE		0	
#endif

#ifndef TRUE
#define TRUE		1	
#endif

#define BooLean		char

#define ILLEGAL		-1.0
#define BELL		((char) 7)	

#define FILE_BEGIN	0
#define FILE_CURRENT	1
#define FILE_END	2

/* CONSTANTS */

typedef short		Int2;
typedef unsigned short	UInt2;

typedef int		Int4;
typedef unsigned int	UInt4;

typedef long		Int8;
typedef unsigned long	UInt8;

#define INT2_MIN        SHRT_MIN
#define INT2_MAX        SHRT_MAX

#define INT4_MIN	INT_MIN
#define INT4_MAX	INT_MAX
#define UINT4_MAX	UINT_MAX
#define INT8_MAX	LONG_MAX
#define INT8_MIN	LONG_MIN
#define UINT8_MAX	ULONG_MAX

/* MACROS - standard macro definitions and static types */
#define	MEW(x,n,t)	(( (x=(t*) malloc(((n)*sizeof(t))))==NULL) ? \
			 (assert(!"Out of Memory."),exit(1),(t*) 0):x)
// 			 (fprintf(stderr,"Out of Memory."),exit(1),(t*) 0):x)

#define	NEW(x,n,t)	(( (x=(t*) calloc(n,sizeof(t)))==NULL) ? \
			 (assert(!"Out of Memory."),exit(1),(t*) 0):x)
// 			 (fprintf(stderr,"Out of Memory."),exit(1),(t*) 0):x)

#define	NEWP(x,n,t)	(( (x=(t**) calloc(n,sizeof(t*)))==NULL) ? \
			 (assert(!"Out of Memory."),exit(1),(t**) 0):x)
// 			 (fprintf(stderr,"Out of Memory."),exit(1),(t**)0):x)

#define	NEWPP(x,n,t)	(( (x=(t***) calloc(n,sizeof(t**)))==NULL) ? \
			 (assert(!"Out of Memory."),exit(1),(t***) 0):x)
// 			 (fprintf(stderr,"Out of Memory."),exit(1),(t***) 0):x)

#define	NEWP3(x,n,t)	(( (x=(t****) calloc(n,sizeof(t***)))==NULL) ? \
			 (assert(!"Out of Memory."),exit(1),(t****) 0):x)
//			 (fprintf(stderr,"Out of Memory."),exit(1),(t****) 0):x)

#define	GETCHAR(m,C)	do{ fprintf(stderr,"%s ", m); \
			  if(fscanf(stdin,"%c",(C)) == 1) { \
                	    while(getchar()!='\n') if(feof(stdin)) exit(1);\
			    break;\
 			  } while(getchar()!='\n') if(feof(stdin)) exit(1);\
			} while(TRUE);

#define print_error(str) (fprintf(stderr,"%s\n",str)? exit(1): exit(1)) 
#define DIGIT2INT(c)    ((int)(c - 48))
#define MINIMUM(t,x,y)	(((t)(x) < (t)(y)) ? (t)(x) : (t)(y))
#define MAXIMUM(t,x,y)	(((t)(x) > (t)(y)) ? (t)(x) : (t)(y))
#define SUM(x)		(((x) * (x+1)) / 2)

#endif

