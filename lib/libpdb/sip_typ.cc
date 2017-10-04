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

#include "sip_typ.h"

double	sip_typ::DegreeOfColinearity(double **pt)
//checks for the degree of colinearity between the three pts by computing
//the minimal hight of the triangle spaned by these pts.
{
	double magn,d1,d2,d3,n1,n2,n3;
	double **x = pt;

	d1=sqrt(pow(x[3][1]-x[1][1],2)+pow(x[3][2]-x[1][2],2)+pow(x[3][3]-x[1][3],2));
	d2=sqrt(pow(x[2][1]-x[1][1],2)+pow(x[2][2]-x[1][2],2)+pow(x[2][3]-x[1][3],2));
	d3=sqrt(pow(x[2][1]-x[3][1],2)+pow(x[2][2]-x[3][2],2)+pow(x[2][3]-x[3][3],2));

        n1=(x[2][2]-x[1][2])*(x[3][3]-x[1][3])-
                (x[3][2]-x[1][2])*(x[2][3]-x[1][3]);
        n2=(x[3][1]-x[1][1])*(x[2][3]-x[1][3])-
                (x[2][1]-x[1][1])*(x[3][3]-x[1][3]);
        n3=(x[2][1]-x[1][1])*(x[3][2]-x[1][2])-
                (x[3][1]-x[1][1])*(x[2][2]-x[1][2]);

	magn=sqrt(pow(n1,2)+pow(n2,2)+pow(n3,2));

	if(((magn/d1)<=(magn/d2)) && ((magn/d1)<=(magn/d3)))	
		return (magn/d1);
	else if ((magn/d2)<=(magn/d3)) 
		return (magn/d2);
	else return (magn/d3);
}

double	*sip_typ::AtomToPoint(atm_typ atm)
{
	double	*Pt;
	NEW(Pt,5,double);
	Pt[1] = AtomX(atm); Pt[2] = AtomY(atm); Pt[3] = AtomZ(atm);
	return Pt;
}

void	sip_typ::Init(atm_typ *ToAtm, atm_typ *FromAtm)
// ptA[4][4] where ptA[1] = center, ptA[2] = linept,  ptA[3] = plane pt.
// same for ptB.  And ptA[1][1] = x, ptsA[1][2] = y, ptA[1][3] = z
// returns NULL if operation failed.
{
	Int4	a,j;
	// crr_typ C; NEW(C,1,pdb_crr_type);
	double **x,**y;
        NEWP(FromPts,5,double); NEWP(ToPts,5,double);
	for(a=1; a<= 3; a++){
		FromPts[a] = AtomToPoint(FromAtm[a]);
		ToPts[a] = AtomToPoint(ToAtm[a]);
        } NEW(oldPt,5,double); NEW(newPt,5,double);

	ptA=x=ToPts; ptB=y=FromPts;
	fprintf(stderr,"colinearlity x =%lf\n",DegreeOfColinearity(x));
	fprintf(stderr,"colinearlity y =%lf\n",DegreeOfColinearity(y));
	double detA,t0,t1,t2,t3,n1,n2,n3,m1,m2,m3,u,v,w,s,d,g,h,i;

	NEWP(ptIM,5,double);
	for(j=0;j<=4;j++) NEW(ptIM[j],5,double);

	ptIM[1][1]=x[1][1]; ptIM[1][2]=x[1][2]; ptIM[1][3]=x[1][3];

	t0=sqrt((pow(y[2][1]-y[1][1],2)+pow(y[2][2]-y[1][2],2)+
		pow(y[2][3]-y[1][3],2))/(pow(x[2][1]-x[1][1],2)+
		pow(x[2][2]-x[1][2],2)+pow(x[2][3]-x[1][3],2)));

	ptIM[2][1]=x[1][1]+(x[2][1]-x[1][1])*t0;
	ptIM[2][2]=x[1][2]+(x[2][2]-x[1][2])*t0; 
	ptIM[2][3]=x[1][3]+(x[2][3]-x[1][3])*t0; 

	n1=(x[2][2]-x[1][2])*(x[3][3]-x[1][3])-
		(x[3][2]-x[1][2])*(x[2][3]-x[1][3]);
	n2=(x[3][1]-x[1][1])*(x[2][3]-x[1][3])-
                (x[2][1]-x[1][1])*(x[3][3]-x[1][3]);
	n3=(x[2][1]-x[1][1])*(x[3][2]-x[1][2])-
                (x[3][1]-x[1][1])*(x[2][2]-x[1][2]);

	m1=n2*(x[2][3]-x[1][3])-n3*(x[2][2]-x[1][2]);
	m2=n3*(x[2][1]-x[1][1])-n1*(x[2][3]-x[1][3]);
	m3=n1*(x[2][2]-x[1][2])-n2*(x[2][1]-x[1][1]);

	t1=((y[3][3]-y[1][3])*(y[2][3]-y[1][3])+(y[3][2]-y[1][2])*
		(y[2][2]-y[1][2])+(y[3][1]-y[1][1])*(y[2][1]-y[1][1]))/
		(pow(y[2][1]-y[1][1],2)+pow(y[2][2]-y[1][2],2)+
		pow(y[2][3]-y[1][3],2));

	g=y[1][1]+(y[2][1]-y[1][1])*t1;
	h=y[1][2]+(y[2][2]-y[1][2])*t1;
	i=y[1][3]+(y[2][3]-y[1][3])*t1;

	s=sqrt(pow(t1,2)*(pow(y[2][1]-y[1][1],2)+pow(y[2][2]-y[1][2],2)+
		pow(y[2][3]-y[1][3],2)));

	d=sqrt(pow(g-y[3][1],2)+pow(h-y[3][2],2)+pow(i-y[3][3],2));	

	t2=s/sqrt(pow(x[2][1]-x[1][1],2)+pow(x[2][2]-x[1][2],2)+
		pow(x[2][3]-x[1][3],2));

	u=x[1][1]+(x[2][1]-x[1][1])*t2;
	v=x[1][2]+(x[2][2]-x[1][2])*t2;
	w=x[1][3]+(x[2][3]-x[1][3])*t2;

	t3=d/sqrt(pow(m1,2)+pow(m2,2)+pow(m3,2));
	
	ptIM[3][1]=u+m1*t3;
	ptIM[3][2]=v+m2*t3;
	ptIM[3][3]=w+m3*t3;

        detA=y[1][1]*y[2][2]*y[3][3]+y[2][1]*y[3][2]*y[1][3]+   
                y[3][1]*y[1][2]*y[2][3]-y[3][1]*y[2][2]*y[1][3]-
                y[1][1]*y[3][2]*y[2][3]-y[2][1]*y[1][2]*y[3][3];

        NEWP(M,5,double);
        for(j=0;j<=4;j++) NEW(M[j],5,double);	

        M[1][1]=(y[2][2]*y[3][3]-y[3][2]*y[2][3])/detA;
        M[1][2]=(y[3][1]*y[2][3]-y[2][1]*y[3][3])/detA;
        M[1][3]=(y[2][1]*y[3][2]-y[3][1]*y[2][2])/detA;
        M[2][1]=(y[3][2]*y[1][3]-y[1][2]*y[3][3])/detA;
        M[2][2]=(y[1][1]*y[3][3]-y[3][1]*y[1][3])/detA;
        M[2][3]=(y[3][1]*y[1][2]-y[1][1]*y[3][2])/detA;
        M[3][1]=(y[1][2]*y[2][3]-y[2][2]*y[1][3])/detA;
        M[3][2]=(y[2][1]*y[1][3]-y[1][1]*y[2][3])/detA;
        M[3][3]=(y[1][1]*y[2][2]-y[2][1]*y[1][2])/detA;
}

atm_typ	sip_typ::PointToAtom(double *Pt,atm_typ atm)
{
	atm_typ A = CopyCoreAtom(atm);
	SetAtomX(Pt[1],A); SetAtomY(Pt[2],A); SetAtomZ(Pt[3],A);
	return A;
}

atm_typ	sip_typ::MapTo(atm_typ atm)  
// void	sip_typ::MapTo(double *oldPt, double *newPt)  
// point[4] where point[1] = x, [2] -> y, [3] -> z.
// maps the point in space to new position given three key pairs of points
// that were passed to Init().  Returns 1 if operation succeeded and
// 0 otherwise.
{
	double x1,x2,x3,y1,y2,y3,m1,m2,m3,n1,n2,n3,t0,t1,al,be,ga;
	oldPt[1] = AtomX(atm); oldPt[2] = AtomY(atm); oldPt[3] = AtomZ(atm);

	double *z=oldPt;

	n1=(ptIM[2][2]-ptIM[1][2])*(ptIM[3][3]-ptIM[1][3])-
                (ptIM[3][2]-ptIM[1][2])*(ptIM[2][3]-ptIM[1][3]);
        n2=(ptIM[3][1]-ptIM[1][1])*(ptIM[2][3]-ptIM[1][3])-
                (ptIM[2][1]-ptIM[1][1])*(ptIM[3][3]-ptIM[1][3]);
        n3=(ptIM[2][1]-ptIM[1][1])*(ptIM[3][2]-ptIM[1][2])-
                (ptIM[3][1]-ptIM[1][1])*(ptIM[2][2]-ptIM[1][2]);

	m1=(ptB[2][2]-ptB[1][2])*(ptB[3][3]-ptB[1][3])-
                (ptB[3][2]-ptB[1][2])*(ptB[2][3]-ptB[1][3]);
        m2=(ptB[3][1]-ptB[1][1])*(ptB[2][3]-ptB[1][3])-
                (ptB[2][1]-ptB[1][1])*(ptB[3][3]-ptB[1][3]);
        m3=(ptB[2][1]-ptB[1][1])*(ptB[3][2]-ptB[1][2])-
                (ptB[3][1]-ptB[1][1])*(ptB[2][2]-ptB[1][2]);	

	t0=((ptB[1][1]-z[1])*m1+(ptB[1][2]-z[2])*m2+(ptB[1][3]-z[3])*m3)/
		(pow(m1,2)+pow(m2,2)+pow(m3,2));

	y1=z[1]+m1*t0; y2=z[2]+m2*t0; y3=z[3]+m3*t0;

	al=M[1][1]*y1+M[1][2]*y2+M[1][3]*y3;
	be=M[2][1]*y1+M[2][2]*y2+M[2][3]*y3;
	ga=M[3][1]*y1+M[3][2]*y2+M[3][3]*y3;

	x1=al*ptIM[1][1]+be*ptIM[2][1]+ga*ptIM[3][1];
	x2=al*ptIM[1][2]+be*ptIM[2][2]+ga*ptIM[3][2];
	x3=al*ptIM[1][3]+be*ptIM[2][3]+ga*ptIM[3][3];
	t1=t0*sqrt((pow(n1,2)+pow(n2,2)+pow(n3,2))/
		(pow(m1,2)+pow(m2,2)+pow(m3,2)));
	newPt[1]=x1-n1*t1;
	newPt[2]=x2-n2*t1;
	newPt[3]=x3-n3*t1;
	return PointToAtom(newPt,atm);
}

void	sip_typ::Free( )
{
   Int4	i;	
   free(oldPt); free(newPt);
   if(ptA != NULL){ for(i=1;i<=3;i++) free(ptA[i]); free(ptA); }
   if(ptB != NULL){ for(i=1;i<=3;i++) free(ptB[i]); free(ptB);      }
   if(ptIM != NULL){ for(i=1;i<=3;i++) free(ptIM[i]); free(ptIM);   }
   if(M != NULL){ for(i=1;i<=3;i++) free(M[i]); free(M); }
}
