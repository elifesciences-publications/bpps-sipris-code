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

#include "dom_typ.h"

dom_typ::dom_typ(FILE *fp) { Init(fp,"REHODtTPI","YGVMYC"); }

dom_typ::dom_typ(FILE *fp,const char *Shapes,const char *Colors)
{ Init(fp,Shapes,Colors); }

void	dom_typ::Init(FILE *fp,const char *Shapes,const char *Colors)
{
	init( );
	num_colors=strlen(Colors);
	colors = new char[num_colors+3];
	strcpy(colors,Colors);

	num_shapes=strlen(Shapes);
	shapes = new char[num_shapes+3];
	strcpy(shapes,Shapes);

	char	str[1000];
	assert(fscanf(fp,"<Domains(%u)=\n",&NumSeqs) == 1);
	domains = new pdl_typ*[NumSeqs+2];
	sq_ids = new char*[NumSeqs+2];
	NumDoms = new unsigned short[NumSeqs+2];
	for(Int4 s=0; s <= NumSeqs; s++){ domains[s]=0; sq_ids[s]=0; }

	UInt4   strt_sq,end_sq;
        UInt4	strt_dom,end_dom,len_dom,extend;
	float	evalue;
	Int4	i=1,err,SqNum;
	char	c;
	do {
	    err = fscanf(fp,"->(%d)%[^=]={",&SqNum, str); 
	    if(err != 2){
		if(fscanf(fp,">.\n") == 0) break;
		else {
		  fprintf(stderr,"%s\n",fgets(str,190,fp));
		  print_error("dom_typ input error");
		}
	    } 
	    assert(SqNum > 0 && SqNum <= NumSeqs);
	    assert(i == SqNum);
	    sq_ids[i] = NewString(str);
	    Int4 J=0;
	    do {
	      if((c=getc(fp)) == '}'){
		ungetc(c,fp);
		if(fscanf(fp,"};\n") == 0) break; 
		else{
		 fprintf(stderr,"%s\n",fgets(str,190,fp));
		 print_error("dom_typ input error");
	        }
	      } 
	      ungetc(c,fp);
              err=fscanf(fp,"[%u-%u->%[^:]:%u-%u(%u);%f]",
                   &strt_sq,&end_sq,str,&strt_dom,&end_dom,&extend,&evalue);
                // evalue == inf???
	      if(err != 7){
		fprintf(stderr,"%s\n",fgets(str,190,fp));
		print_error("dom_typ input error"); 
	      }
	      J++;
	      len_dom=end_dom+extend;
	      if(domains[i]==0){
		domains[SqNum] = 
		  new pdl_typ(strt_sq,end_sq,strt_dom,end_dom,len_dom,str,evalue);
	      } else {
		domains[SqNum]->Append(strt_sq,end_sq,strt_dom,end_dom,
			len_dom,str,evalue);
	      }
	    } while(TRUE);
	    NumDoms[SqNum] = J;
	    i++;
	} while(TRUE);
}

Int4    dom_typ::SqStart(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->SqStart(d);
	else return 0; // this should not happen.
}

Int4    dom_typ::SqEnd(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->SqEnd(d);
	else return 0; // this should not happen.
}

char    *dom_typ::DomainName(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->DomainName(d);
	else return 0; // this should not happen.
}

float	dom_typ::Evalue(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->Eval(d);
	else return 0; // this should not happen.
}

Int4	dom_typ::NumDomains(UInt4 s)
{ if(s < 1 || s > NumSeqs) return 0; else return NumDoms[s]; }

void    dom_typ::Put(FILE *fp, BooLean *skip, char *skipdom)
{ Put(fp, skip, skipdom, 0); }

void    dom_typ::Put(FILE *fp, BooLean *skip, char *skipdom, Int4 offset)
{
	Int4	s,S,n;
	assert(offset >= 0);
	if(skip){
	   for(n=0,s=1; s <= NumSeqs; s++) if(!skip[s]) n++;
	} else n=NumSeqs;
	fprintf(fp,"<Domains(%u)=\n",n);
	for(s=0,S=1; S <= NumSeqs; S++){
	   if(skip && skip[S]) continue;
	   else s++;
	   fprintf(fp,"->(%u)%s={",s+offset,sq_ids[S]);
	   if(domains[S]) domains[S]->Put(fp,skipdom);
	   fprintf(fp,"};\n");
	} fprintf(fp,">.\n");
}

char	dom_typ::Color(UInt4 s, unsigned short d)
{
	if(!colored) Color();	// initializes colors.
	pdl_typ *pdl=domains[s]->Domain(d);
	return pdl->Color( );
}

char	dom_typ::Shape(UInt4 s, unsigned short d)
{
	if(!colored) Color();	// initializes colors & shapes.
	pdl_typ *pdl=domains[s]->Domain(d);
	return pdl->Shape( );
}


Int4	dom_typ::Color()
{
	Int4	s,i,j,d,n,m,N=0,set1,set2,NSets;
	char	*name1,*name2,**name;
	char	biased_shapes[5] = "esSC",*str;

	colored=TRUE;
	for(s=1; s <= NumSeqs; s++){ N += NumDomains(s); }
	ds_type	Sets=DSets(N); NSets=N;
	pdl_typ **list = new pdl_typ*[N+3];
	// NEWP(name,N+3, char);
	name = new char*[N+3];
	for(n=0;n <= N; n++) name[n] = 0;
	for(n=0,s=1; s <= NumSeqs; s++){
	  for(d=NumDomains(s); d > 0; d--){
		n++; list[n] = domains[s]->Domain(d);
		name[n] = domains[s]->DomainName(d);
	  }
	}
	for(n=1; n < N; n++){
	    name1 = name[n]; set1 = findDSets(n,Sets);
	    for(m=n+1; m <= N; m++){
	    	name2 = name[m];
		set2 = findDSets(m,Sets);
		if(set1 != set2 && strcmp(name1,name2)==0){
			NSets--; set1 = linkDSets(set1,set2,Sets);
		}
	    }
	}
	char *colormap,*shapemap;
	Int4	c,t,x,rs,rc;
	// NEW(colormap,N+3,char); NEW(shapemap,N+3,char);
	colormap = new char[N+3];
	shapemap = new char[N+3];
	// Store domain information.
	dom_numtypes = NSets;
	dom_name = new char*[NSets + 3];
	dom_color = new char[NSets + 3];
	dom_shape = new char[NSets + 3];
	dom_num = new unsigned short[NSets + 3];
	dom_avelen = new UInt4[NSets + 3];
	Int4 *dommap;
	// NEW(dommap,N+3,Int4);
	dommap = new Int4[N+3];
	for(n=1; n <= N; n++){ dommap[n]=0; colormap[n]=shapemap[n]=0; }
	for(t=0; t <= NSets; t++){
		dom_num[t]=0; dom_avelen[t] = 0; 
		dom_name[t]=0;
	}
	for(t=1,rs=rc=s=c=0,n=1; n <= N; n++){
	    set1 = findDSets(n,Sets);
	    if(colormap[set1] == 0) { 
		dommap[set1]=t; 
		dom_name[t] = name[set1];
		if(strcmp(name[set1],"TM_helix") == 0){
		  dom_color[t] = colormap[set1] = 'Y'; // set to yellow.
		  dom_shape[t] = shapemap[set1] = 'r';
		  assert(t <= NSets); t++; 
		} else if(strcmp(name[set1],"coiled_coil") == 0){
		  dom_color[t] = colormap[set1] = 'B'; // set to black.
		  dom_shape[t] = shapemap[set1] = 'r';
		  assert(t <= NSets); t++; 
		} else if((str=strstr(name[set1],"_rich")) 
					&& strcmp("_rich",str)==0){
		  dom_color[t] = colormap[set1] = colors[rc];
		  dom_shape[t] = shapemap[set1] = 'e'; // old: = biased_shapes[rs]; 
		  assert(t <= NSets); rc++; t++; 
		  if(rc >= num_colors){ rc=0; rs++; } // black & white are reserved colors.
		  if(rs >= 4){ rs=0; }
		} else { 
		  dom_shape[t] = shapemap[set1] = shapes[s]; 
		  dom_color[t] = colormap[set1] = colors[c]; 
		  assert(t <= NSets); c++; t++; 
		  if(c >= num_colors){ c=0; s++; } // black & white are reserved colors.
		  if(s >= num_shapes){ s=0; }
		}
	    }
	    x = dommap[set1]; 
	    dom_num[x]++; dom_avelen[x] += list[n]->LengDom( );
	    list[n]->SetColor(colormap[set1]);
	    list[n]->SetShape(shapemap[set1]);
	}
	for(t=1; t <= NSets; t++){ 
	    dom_avelen[t]=(UInt4)((double)dom_avelen[t]/(double)dom_num[t]);
	}
	// Free memory.
	for(n=0;n <= N; n++) { name[n] = 0; list[n] = 0; } 
	for(n=1; n <= N; n++){ dommap[n]=0; colormap[n]=shapemap[n]=0; }
	delete [] list; NilDSets(Sets); 
	delete [] name;
	// free(colormap); free(shapemap); free(dommap);
	delete [] colormap; delete [] shapemap; delete [] dommap;
	return NSets;
}

void    dom_typ::Free()
{
	for(Int4 s=1; s <= NumSeqs; s++){ 
	   if(domains[s]) delete domains[s];
	   if(sq_ids[s]) delete [] sq_ids[s];
	}
	if(shapes) delete [] shapes;
	if(sq_ids) delete [] sq_ids;
	if(domains) delete [] domains;
	if(NumDoms) delete [] NumDoms;
	if(dom_name) delete [] dom_name;
	if(dom_shape) delete [] dom_shape;
	if(dom_color) delete [] dom_color;
	if(dom_num) delete [] dom_num;
	if(dom_avelen) delete [] dom_avelen;
	init();
}

void    dom_typ::init()
{ 
  shapes=0; num_shapes=0;
  colors=0; num_colors=0;
  domains = 0; sq_ids = 0; NumSeqs=0; NumDoms=0; colored=FALSE; 
  dom_name=0; dom_shape=0; dom_color=0; dom_num=0; dom_avelen=0;
}

dom_typ& dom_typ::operator=(const dom_typ& dom)
// called for dom_typ dom2; dom2=dom;
{ if (this != &dom) { Free(); copy(dom); } return *this; }

void    dom_typ::copy(const dom_typ& dom) 
// private function that assumes 'this' has been initialized and
// copies dom to 'this'.
{
	NumSeqs=dom.NumSeqs;
	
	domains = new pdl_typ*[NumSeqs+2];
	sq_ids = new char*[NumSeqs+2];
	NumDoms = new unsigned short[NumSeqs+2];
	for(Int4 s=1; s <= NumSeqs; s++){
	   if(dom.domains[s]){ 
		domains[s] = new pdl_typ;
		domains[s]->Copy(*(dom.domains[s]));
		// domains[s] = dom.domains[s]; 
	   }
	   else domains[s] = 0;
	   if(dom.sq_ids[s]){ sq_ids[s] = NewString(dom.sq_ids[s]); }
	   else sq_ids[s]=0; 
	   NumDoms[s] = dom.NumDoms[s];
	}
	colored = dom.colored;
}

//*************************** pdl_typ *******************************

pdl_typ& pdl_typ::operator=(const pdl_typ& pdl)
// called for pdl_typ pdl2; pdl2=pdl;
{ if (this != &pdl) { Free(); copy(pdl); } 
assert(!"don't use this routine"); return *this; }

void    pdl_typ::copy(const pdl_typ& pdl) 
// private function that assumes 'this' has been initialized and
// copies pdl to 'this'.
{
	strt_sq = pdl.strt_sq;
	end_sq = pdl.end_sq;
        strt_dom = pdl.strt_dom;
	end_dom = pdl.end_dom;
	len_dom = pdl.len_dom;
        dom_name=NewString(pdl.dom_name);
        dom_num = pdl.dom_num;  
        Evalue = pdl.Evalue;
        color = pdl.color;
	if(pdl.next){ next = new pdl_typ; next->Copy(*(pdl.next)); }
	else next=0;
}

pdl_typ::pdl_typ(UInt4 ss, UInt4 es,
        unsigned short sd, unsigned short ed, unsigned short ld, 
	char *name,float e)
{
	init();
	strt_sq =ss; end_sq = es;
        strt_dom = sd; end_dom = ed; len_dom=ld;
        dom_name = NewString(name);
        Evalue = e;
	dom_num=1;
	next = 0;
}

void	pdl_typ::Append(UInt4 ss, UInt4 es,
        unsigned short sd, unsigned short ed, unsigned short ld,
	char *name,float e)
{
	pdl_typ *tail,*nxt;
	short	dn=2;
	if(next){
	  for(nxt=next; nxt; nxt=nxt->next){ tail = nxt; dn++;}
	  tail->next = new pdl_typ(ss,es,sd,ed,ld,name,e);
	  tail->next->dom_num = dn;
	} else {
	  next = new pdl_typ(ss,es,sd,ed,ld,name,e);
	  next->dom_num = dn;
	}
}

Int4    pdl_typ::SqStart( )
{ if(this == 0) return 0; else return strt_sq; }

Int4    pdl_typ::SqEnd( )
{ if(this == 0) return 0; else return end_sq; }

Int4    pdl_typ::SqStart(unsigned short d)
{
	if(dom_num == d) return strt_sq; else if(next == 0) return 0;
	else return next->SqStart(d);
}

Int4    pdl_typ::SqEnd(unsigned short d)
{
	if(dom_num == d) return end_sq; else if(next == 0) return 0;
	else return next->SqEnd(d);
}

pdl_typ	*pdl_typ::Domain(unsigned short d)
{
	if(dom_num == d) return this; 
	else if(next == 0) return 0;
	else return next->Domain(d);
}

char    *pdl_typ::DomainName(unsigned short d)
{
	if(dom_num == d) return dom_name; else if(next == 0) return 0;
	else return next->DomainName(d);
}

float	pdl_typ::Eval(unsigned short d)
{
	if(dom_num == d) return Evalue; else if(next == 0) return -9999.;
	else return next->Eval(d);
}

void	pdl_typ::Put(FILE *fp, char *skip_name)
// [start-end->domain:start-end(extend);-log10(E-value)]
{
	if(!skip_name || (strcmp(skip_name,dom_name) != 0)){
	   fprintf(fp,"[%d-%d->%s:%d-%d(%d);%.1f]",
		strt_sq,end_sq,dom_name,strt_dom,
		end_dom,len_dom-end_dom,Evalue);
	}
	if(next) next->Put(fp,skip_name);
}

void	pdl_typ::Free( )
{
	if(next) delete next;
	if(dom_name) delete [] dom_name;
	init( );
}

void	pdl_typ::init( )
{
	strt_sq=end_sq=0;
	strt_dom=end_dom=0; len_dom=0;
	Evalue=0.0; dom_name=0;
	color = 0; shape=0;
	next=0;
}

