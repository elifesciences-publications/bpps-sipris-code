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

#include "tpb_typ.h"

tpb_typ::tpb_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length,
	unsigned short Rpts, Int4 NumRes,double *ExpGap,double *ExpRptGap,Int4 Pernats)
// constructor for multiblock tpb_typ...
{
	unsigned short BlkLen[3]; BlkLen[1] = Length;
	Init(IdpArg, Mode,Prob,Length,Rpts,1,BlkLen,NumRes,ExpGap,ExpRptGap,Pernats);
}

tpb_typ::tpb_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length, 
	unsigned short Rpts, Int4 NumRes,unsigned short nBlk,
	unsigned short *BlkLen,double *ExpGap,double *ExpRptGap,Int4 Pernats)
// new constructor for tpb_typ...
{ Init(IdpArg, Mode,Prob,Length,Rpts,nBlk,BlkLen,NumRes,ExpGap,ExpRptGap,Pernats); }

tpb_typ& tpb_typ::operator=(const tpb_typ& tpb)
// called for tpb_typ tpb2; tpb2=tpb; 
{
	if (this != &tpb) { 
		Free(); 
		Init(tpb.tpb_arg, tpb.mode,tpb.prob_c,tpb.length, tpb.rpts,
			tpb.nblk, tpb.blklen, tpb.num_res,tpb.exp_gap,tpb.exp_rpt_gap,tpb.pernats);
	} return *this; 
}

void	tpb_typ::Init(char* IdpArg, char Mode,double *Prob,
	unsigned short Length, unsigned short Rpts, unsigned short nBlk,
	unsigned short *BlkLen,Int4 NumRes,double *ExpGap,double *ExpRptGap,Int4 Pernats)
// input IdpArg as "-P100..1100,20..120:500,40..400"
{
        Int4    iol=100,ioh=1100,ixl=20,ixh=120,dol=500,doh=500,dxl=40,dxh=400;

        if(IdpArg && IdpArg[2] !=0){
	   tpb_arg = AllocString(IdpArg);
           if(sscanf(tpb_arg,"-P%d..%d,%d..%d:%d..%d,%d..%d",
                    &iol,&ioh,&ixl,&ixh,&dol,&doh,&dxl,&dxh) != 8){
                                 print_error("tpb_typ constructor input error");
           }
	   if(!(iol <= ioh && ixl <= ixh) || !(dxl <= dxh) 
			|| (iol < 0 || ixl < 0 || dxl < 0 || 
			  dol < 0 || dol > doh))
                                 print_error("tpb_typ argument input error");
                // std::cerr << tpb_arg; std::cerr << "\n";
        } // else use default settings.

	pernats=Pernats;
	length = Length; rpts = Rpts; mode=Mode; nblk = nBlk;
	inso_low = iol; inso_high = ioh; insx_low = ixl;insx_high=ixh;
        delo_low = dol; delo_high=doh; delx_low = dxl; delx_high = dxh; 
	num_res=0; num_res=NumRes; 

	NEW(blklen,nblk+3,unsigned short);
	Int4	k,n=0;
	for(k=1; k <= nblk; k++){ n+=blklen[k]=BlkLen[k]; }
	assert(n == length);

	NEW(prob_c,Length+3,double);
	for(Int4 i=0; i<=Length; i++){ prob_c[i]=Prob[i]; }
	
	NEW(exp_gap,nblk+3,double);
	for(k=1; k <= nblk; k++){
		exp_gap[k] = ExpGap[k];
		// fprintf(stderr,"exp_gap[%d] = %.2f\n",k,exp_gap[k]); 
	}
	NEW(exp_rpt_gap,Rpts+3,double);
	for(k=0; k<=Rpts; k++) { exp_rpt_gap[k] = ExpRptGap[k]; }
 	init( );
}

void	tpb_typ::Free( )
{
	free(m2m); free(i2m); free(d2m);
	free(m2i); free(i2i); free(m2d); free(d2d); 
	free(tpb_arg); free(blklen);
	free(exp_gap); free(exp_rpt_gap);
	free(prob_c);
}

Int4	tpb_typ::EndBlk(Int4 blk)
{
	assert(blk > 0 && blk <= nblk);
	Int4	end=0;
	for(Int4 b=1; b <= blk; b++) end+=blklen[b];
	return end;
}

Int4	*tpb_typ::StartBlk( )
{
	Int4	s,b,*Start;

	NEW(Start,nblk+3,Int4);
	for(s=1,b=1; b <= nblk; b++){ Start[b]=s; s+=blklen[b]; }
	Start[b]=s;
	return Start;
}

Int4	tpb_typ::StartBlk(Int4 blk)
{
	assert(blk > 0 && blk <= nblk);
	Int4	start=0;
	for(Int4 b=1; b < blk; b++) start+=blklen[b];
	return (start+1);
}

void	tpb_typ::init( )
// ins_open = m2i,ins_extend=i2i,del_open=m2d,
// del_extend=d2d
{
        Int4    i,j,k;
	double	temp,temp2;

        NEW(m2i,rpts*length+3,Int4);
        NEW(i2i,rpts*length+3,Int4);

        NEW(m2d,rpts*length+3,Int4);
        NEW(d2d,rpts*length+3,Int4);

        NEW(i2m,rpts*length+3,Int4);
        NEW(d2m,rpts*length+3,Int4);
        NEW(m2m,rpts*length+3,Int4);

	// i2i[0]=floor(0.5+pernats*log(exp_rpt_gap[0]/(exp_rpt_gap[0]+1.)));
	// i2i[0]=floor(0.5+pernats*log(350./351.));
	i2i[0]=0;
#if 0
	m2m[0]=floor(0.5+pernats*log(1-exp(temp/pernats)));
	i2i[length+1]=floor(0.5+pernats*log(exp_rpt_gap[nblk]/(exp_rpt_gap[nblk]+1.)));		
	i2m[length+1]=floor(0.5+pernats*log(1./(exp_rpt_gap[nblk]+1.)));
	temp=(-inso_low-insx_low)*(prob_c[1]) +(-inso_high-insx_high)*(1.0-prob_c[1]);
	m2i[0]=(Int4) floor(0.5+temp);
#endif
	m2i[0]=0;
	temp = (-delo_low)*prob_c[1]+(-delo_high)*(1.0-prob_c[1]);
	m2d[0]=(Int4) floor(0.5+temp);
	Int4 del_pos;
        for(i=1,k=1;k<=nblk;k++){
		for(j=1;j<blklen[k];j++){
			temp2 = (prob_c[i]+prob_c[i+1])/2.0;
			temp=(-inso_low-insx_low)*(temp2)
					+(-inso_high-insx_high)*(1.0-temp2);
	                m2i[i]=(Int4) floor(0.5+temp);
			temp=(-insx_low)*temp2+(-insx_high)*(1.0-temp2); 
	                i2i[i]=(Int4) floor(0.5+temp);

			temp=(-delo_low-delx_low)*prob_c[i+1]
					+(-delo_high-delx_high)*(1.0-prob_c[i+1]);
        	        m2d[i]=(Int4) floor(0.5+temp);
			temp=(-delx_low)*prob_c[i+1]+(-delx_high)*(1.0-prob_c[i+1]);
                	d2d[i]=(Int4) floor(0.5+temp);

			i++;
		}
		if(k!=nblk){
			if(i%length==0) del_pos=1;
			else del_pos = i+1;
			i2i[i]=floor(0.5+pernats*log(exp_gap[k]/(exp_gap[k]+1.)));

			m2i[i]=floor(0.5+pernats*log(exp_gap[k]/(exp_gap[k]+1.)));
#if 1
			temp=exp(((-delo_low-delx_low)*prob_c[del_pos]
				+(-delo_high-delx_high)*(1.0-prob_c[del_pos]))/pernats);
			m2d[i]=floor(0.5+pernats*log(temp/(exp_gap[k]+1.)));
#endif
#if 0
			temp=((-delo_low-delx_low)*prob_c[del_pos]
				+(-delo_high-delx_high)*(1.0-prob_c[del_pos]));
			m2d[i]=floor(0.5+temp);
#endif
			temp=(-delx_low)*prob_c[del_pos]+(-delx_high)*(1.0-prob_c[del_pos]);
			d2d[i]=(Int4) floor(0.5+temp);

       	         	i++;
		}
        }
#if 0	// This doesn't seem to work for HMMs...
        if(mode == 'S'){
          for(i=1,k=1;k<=nblk;k++){
	    Int4 len = blklen[k];
	    Int4 win = MINIMUM(Int4,8,len/2);
            for(j=1; j < win; j++) {
		Int4 los=i+j-1,ros=i+len-j;		// right and left offsets
		Int4 rj=win-j;
                // lower insert transition probabilities near ends.
                m2i[los]=(j*m2i[los]+(rj)*(-insx_high-inso_high))/win; 
                m2i[ros]=(j*m2i[ros]+(rj)*(-insx_high-inso_high))/win;

                i2i[los]=(j*i2i[los]+(rj)*(-insx_high))/win; 
                i2i[ros]=(j*i2i[ros]+(rj)*(-insx_high))/win;

                // raise delete transition probabilities near ends.
                d2d[los]=(j*d2d[los]+(rj)*(-delx_low))/win;
                d2d[ros]=(j*d2d[ros]+(rj)*(-delx_low))/win;
            } i+= len;
	  } i--; m2i[i]=i2i[i]=d2d[i]=0; // reset last position to zero.
        }
#endif
        for(j=1;j<rpts;j++){
            for(i=1;i<length;i++){
                m2i[j*length+i]=m2i[i];
                m2d[j*length+i]=m2d[i]; 
		i2i[j*length+i]=i2i[i];
		d2d[j*length+i]=d2d[i];
            }
        }
#if 1	// AFN: Add back in my routine...
        // Calculate i2m, m2m,*m2d, d2m, // later i2d, & d2i.
        double tmp_d,tmp_i;
        tmp_d= exp((double)m2d[0]/(double)pernats);
        tmp_i= exp((double)m2i[0]/(double)pernats);
        // m2m[0] = (Int4) floor((pernats*log(1.0 - (tmp_d + tmp_i))+0.5));
	m2m[0] = 0;
        for(i=1;i<=length;i++){
          if(m2i[i] != 0){
                tmp_i= exp((double)i2i[i]/(double)pernats);
                i2m[i] = (Int4) floor((pernats*log(1.0 - tmp_i))+0.5);
// fprintf(stderr,"tmp_i=%g; i2m[%d]=%d\n",tmp_i,i,i2m[i]);

                tmp_d= exp((double)m2d[i]/(double)pernats);
                tmp_i= exp((double)m2i[i]/(double)pernats);
                m2m[i] = (Int4) floor((pernats*log(1.0 - (tmp_d + tmp_i))+0.5));

                tmp_d= exp((double)d2d[i]/(double)pernats);
                d2m[i] = (Int4) floor((pernats*log(1.0 - tmp_d))+0.5);
          }
        }
        for(j=1;j<rpts;j++){
           for(i=1;i<=length;i++){
                m2m[j*length+i]=m2m[i]; 
                i2m[j*length+i]=i2m[i];
                d2m[j*length+i]=d2m[i]; 
           }
        }
#endif	// AFN: end Add back in my routine...

#if 0
	for(j=1;j<rpts;j++){
		m2i[j*length]=i2i[j*length]=floor(0.5+pernats*log(exp_rpt_gap[j]/(exp_rpt_gap[j]+1.)));
		temp=(-delo_low-delx_low)*prob_c[1]+(-delo_high-delx_high)*(1.0-prob_c[1]);
		tm=exp(temp/pernats);
                m2d[j*length]=floor(0.5+pernats*log(tm/(exp_rpt_gap[j]+1.)));
		m2m[j*length]=floor(0.5+pernats*log((1.-tm)/(exp_rpt_gap[j]+1.)));
                i2m[j*length]=floor(0.5+pernats*log(1./(exp_rpt_gap[j]+1.)));
		temp=(-delx_low)*prob_c[1]+(-delx_high)*(1.0-prob_c[1]);
                d2d[j*length]=temp;
		tm=exp(temp/pernats);
                d2m[j*length]=floor(0.5+pernats*log(1.-tm));
	}	
#endif
}

void	tpb_typ::Put(FILE *fp)
{
        fprintf(fp,"POS:  m2i  i2i  m2d  d2d  i2m  m2m  d2m  (COIL REGULAR)\n");
        for(Int4 i=1;i<=length;i++){
            fprintf(fp,"%-3d: %4d %4d %4d %4d %4d %4d %4d (%.2f  %.2f)\n",
                       i,-m2i[i],-i2i[i],-m2d[i],-d2d[i],-i2m[i],-m2m[i],-d2m[i],
                        prob_c[i],1.0-prob_c[i]);
        }
}


