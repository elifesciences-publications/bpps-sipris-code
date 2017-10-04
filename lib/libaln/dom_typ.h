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

#if !defined(_DOM_TYP_)
#define _DOM_TYP_
#include "stdinc.h"
#include "afnio.h"
#include "dsets.h"

/************************* Locations of known Domains **************************
<Domains(#seqs)=
->(sq_num)seq_id={[start-end->domain:start-end(extend);-log10(E-value)]...};
>.

Example:
<Domains(3)=
->(1)ctnD_worm={[132-218->fn3:1-84(0);1.6][240-327->fn3:1-84(0);13.2][353-434->fn3:1-84(0);8.4][453-540->fn3:1-84(0);4.4][649-689->Armadillo_seg:1-42(0);1.6][693-734->Armadillo_seg:1-42(0);1.6][795-836->Armadillo_seg:1-42(0);-0.9][881-926->Armadillo_seg:1-42(0);-1.0][932-972->Armadillo_seg:1-42(0);0.1]};
->(2)STU2_YEAST={};
->(3)plakophilin_mouse={[277-317->Armadillo_seg:1-42(0);4.7][319-360->Armadillo_seg:1-42(0);-0.1][516-557->Armadillo_seg:1-42(0);-0.4]};
>.

 *******************************************************************************/
class pdl_typ {  // protein domain list type
public:
	pdl_typ(){ init(); }
	pdl_typ& operator=(const pdl_typ& pdl);
	pdl_typ(UInt4 , UInt4 , unsigned short , 
		unsigned short, unsigned short, char *, float );
	// pdl_typ(FILE *); // ???
	void	Put(FILE *,char *);
	void	Put(FILE *fp){ Put(fp,0); }
	void	Append(UInt4 , UInt4 , unsigned short , 
		unsigned short, unsigned short, char *, float );
	Int4	SqStart( );
	Int4	SqStart(unsigned short d);
	Int4	SqEnd( );
	Int4	SqEnd(unsigned short d);
	float	Eval(unsigned short d);
	char	*DomainName(unsigned short d);
	pdl_typ	*Domain(unsigned short d);
	Int4	LengDom( ){ return len_dom; }
	char	Color( ){ return color; }
	char	Shape( ){ return shape; }
	void	SetColor(char c){ color = c; }
	void	SetShape(char s){ shape = s; }
	void	Copy(const pdl_typ& pdl){ copy(pdl); }
	~pdl_typ( ){ Free(); }
private:
	void    	copy(const pdl_typ& pdl);
	void		Free();		// free memory...
	void		init( );

	UInt4	strt_sq,end_sq;
	unsigned short	strt_dom,end_dom,len_dom;
	char		*dom_name;
	short		dom_num;  // number of domain on linked list.
	float		Evalue;
	pdl_typ		*next;
	char		color;	  // for display: 'R','G','V','M','Y','C','B','W'.
	char		shape;	  // for display: 'R','E','H','O','D','t','T','P','I'.
};

class dom_typ {
public: 
	dom_typ( ){ init(); }
	dom_typ(FILE *);
	dom_typ(FILE *,const char *,const char *);
	dom_typ& operator=(const dom_typ&);     // assignment operator
	~dom_typ( ){ Free(); }
	void	Put(FILE *,BooLean *,char *,Int4);
	void	Put(FILE *,BooLean *,char *);
	void	Put(FILE *fp,BooLean *skip){ Put(fp,skip,0); }
	void	Put(FILE *fp){ Put(fp,0,0); }
	Int4	SqStart(UInt4 s, unsigned short d);
	Int4	SqEnd(UInt4 s, unsigned short d);
	float	Evalue(UInt4 s, unsigned short d);
	char	*DomainName(UInt4 s, unsigned short d);
	Int4	NumDomains(UInt4 s);
	char	Color(UInt4 s, unsigned short d);
	char	Shape(UInt4 s, unsigned short d);
	UInt4	NumSeq( ){ return NumSeqs; }
	unsigned short	NumTypes( ){ if(!colored) Color( ); return dom_numtypes; }
	char		*NameThisType(unsigned short t){
		  	   if(!colored) Color( );
			   if(t > dom_numtypes) return 0; return dom_name[t];
			}
	unsigned short	NumThisType(unsigned short t){
		  	   if(!colored) Color( );
			   if(t > dom_numtypes) return 0; return dom_num[t];
			}
	char		Color(unsigned short t){
		  	   if(!colored) Color( );
			   if(t > dom_numtypes) return 0; return dom_color[t];
			}
	char		Shape(unsigned short t){
		  	   if(!colored) Color( );
			   if(t > dom_numtypes) return 0; return dom_shape[t];
			}
	UInt4	AveLen(unsigned short t){
		  	   if(!colored) Color( );
			   if(t > dom_numtypes) return 0; return dom_avelen[t];
			}
private:
	Int4		Color( );
	void		copy(const dom_typ& dom);
	void		init();
	void		Init(FILE *fp,const char *,const char *);
	void		Free();		// free memory...
	UInt4	NumSeqs;
	char		**sq_ids;
	pdl_typ		**domains;
	BooLean		colored;
	char		*shapes;
	Int4		num_shapes;
	char		*colors;
	Int4		num_colors;
	unsigned short	*NumDoms;
	// information on the types of domains.
	unsigned short	dom_numtypes;
	char		**dom_name,*dom_color,*dom_shape;
	unsigned short	*dom_num;
	UInt4	*dom_avelen;
};

#endif

