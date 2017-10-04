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

#include "idp_typ.h"

idp_typ::idp_typ(idp_typ *idp1, double *prob2, double wt1)
{
	double *Prob;
	assert(wt1 >= 0.0 && wt1 <= 1.0);
	NEW(Prob,idp1->length+3,double);
	for(Int4 i=1; i<=idp1->length; i++){ 
	  Prob[i]=idp1->prob[i]*wt1 + prob2[i]*(1.0 - wt1); 
	}
	Init(idp1->idp_arg,idp1->mode,Prob,idp1->length,idp1->rpts, 
		idp1->nblk,idp1->blklen,idp1->num_res,idp1->seq_ins_prob);
	free(Prob);
}

idp_typ::idp_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length,
	unsigned short Rpts, Int4 NumRes,float *SeqProb)
// constructor for multiblock idp_typ...
{
	unsigned short BlkLen[3]; BlkLen[1] = Length;
	Init(IdpArg, Mode,Prob,Length,Rpts,1,BlkLen,NumRes,SeqProb);
}

idp_typ::idp_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length, 
	unsigned short Rpts, Int4 NumRes,float *SeqProb,unsigned short nBlk,
	unsigned short *BlkLen,Int4 PerNats)
{ Init(IdpArg, Mode,Prob,Length,Rpts,nBlk,BlkLen,NumRes,SeqProb,PerNats); }

idp_typ::idp_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length, 
	unsigned short Rpts, Int4 NumRes,float *SeqProb,unsigned short nBlk,
	unsigned short *BlkLen)
// new constructor for idp_typ...
{ Init(IdpArg, Mode,Prob,Length,Rpts,nBlk,BlkLen,NumRes,SeqProb); }

idp_typ& idp_typ::operator=(const idp_typ& idp)
// called for idp_typ idp2; idp2=idp; 
{
	if (this != &idp) { 
		Free(); 
		Init(idp.idp_arg,idp.mode,idp.prob,idp.length,idp.rpts,
			idp.nblk,idp.blklen,idp.num_res,idp.seq_ins_prob);
	} return *this; 
}

void	idp_typ::Init(char *IdpArg, char Mode,double *Prob,
	unsigned short Length, unsigned short Rpts, unsigned short nBlk,
	unsigned short *BlkLen,Int4 NumRes,float *SeqProb)
{ Init(IdpArg, Mode,Prob,Length,Rpts,nBlk,BlkLen,NumRes,SeqProb,693); }
// 693 pernats == 1000 perbits.

void	idp_typ::Init(char* IdpArg, char Mode,double *Prob,
	unsigned short Length, unsigned short Rpts, unsigned short nBlk,
	unsigned short *BlkLen,Int4 NumRes,float *SeqProb,Int4 Pernats)
// input IdpArg as "-P100..1100,20..120:500,40..400"
// or as "-P100..1100,20..120/20..200:500,40..400"
// or as "-P100..1100,20..120/20..200:500,40..400+i0..30" == length cutoff
// or as "-P100..1100,20..120/20..200:500,40..400+i20" == interblock penalty
{
        Int4    iol=100,ioh=1100,ixl=20,ixh=120,dolh=500,dxl=40,dxh=400;
        Int4    isl=0,ish=0,ibl=-1,ibh=-1;
	Int4	dol=0,doh=0;
	BooLean	use_interblock_penalty=FALSE;

        if(IdpArg && IdpArg[2] !=0){
	   idp_arg = AllocString(IdpArg);
           if(sscanf(idp_arg,"-P%d..%d,%d..%d:%d,%d..%d",
                    &iol,&ioh,&ixl,&ixh,&dolh,&dxl,&dxh) != 7){
             if(sscanf(idp_arg,"-P%d..%d,%d..%d:%d..%d,%d..%d",
                    &iol,&ioh,&ixl,&ixh,&dol,&doh,&dxl,&dxh) != 8){
               if(sscanf(idp_arg,"-P%d..%d,%d..%d/%d..%d:%d,%d..%d+i%d..%d",
                      &iol,&ioh,&ixl,&ixh,&isl,&ish,&dolh,&dxl,&dxh,&ibl,&ibh) != 11){
                 if(sscanf(idp_arg,"-P%d..%d,%d..%d/%d..%d:%d,%d..%d+i%d",
                      &iol,&ioh,&ixl,&ixh,&isl,&ish,&dolh,&dxl,&dxh,&ibl) != 10){
                    if(sscanf(idp_arg,"-P%d..%d,%d..%d/%d..%d:%d,%d..%d",
                      &iol,&ioh,&ixl,&ixh,&isl,&ish,&dolh,&dxl,&dxh) != 9){
                                 // print_error("idp_typ constructor input error");
                                 assert(!"idp_typ constructor input error");
                    }
	         } else {  // check interblkL & interblkH with gap penalty
		     // if(ibl < 0) print_error("idp_typ argument input error");
		     if(ibl < 0) assert(!"idp_typ argument input error");
		     ibh=-1;
		 }
	       } else {  // check interblkL & interblkH
		     // if(ibl < 0 || ibh < ibl) print_error("idp_typ argument input error");
		     if(ibl < 0 || ibh < ibl) assert(!"idp_typ argument input error");
	       }
	     } else {	// if sscanf( ) == 8.
		     if(dol == doh) dolh=dol;
		     else print_error("dol/doh not yet implemented: idp_typ argument input error");
	     }
           }
	   if(!(isl <= ish && isl >= 0) || !(iol <= ioh && ixl <= ixh) 
			|| !(dxl <= dxh) 
			|| (iol < 0 || ixl < 0 || dxl < 0 || dolh < 0)){
                                 // print_error("idp_typ argument input error");
                                 assert(!"idp_typ argument input error");
	   }
                // std::cerr << idp_arg; std::cerr << "\n";
        } else idp_arg=0; // use default settings.
	if(dol == doh ==0) { dol=doh=dolh; }	// use same value for both...

	pernats=Pernats;
	length = Length; rpts = Rpts; mode=Mode; nblk = nBlk;
	interblkL=ibl; interblkH=ibh;
	inso_low = iol; inso_high = ioh; insx_low = ixl;insx_high=ixh;
        delo = dolh; 
	delo_low=dol; delo_high = doh;
	delx_low = dxl; delx_high = dxh; 
	ins_seq_low=ins_seq_high=0;
	num_res=0; seq_ins_prob=0;
	num_res=NumRes; ins_seq_low=isl; ins_seq_high=ish;

	NEW(blklen,nblk+3,unsigned short);
	Int4	k,n=0;
	for(k=1; k <= nblk; k++){ n+=blklen[k]=BlkLen[k]; }
	assert(n == length);

	NEW(prob,Length+3,double);
	for(Int4 i=0; i<=Length; i++){ prob[i]=Prob[i]; }
 	init( );

	NEW(seq_ins_prob,NumRes+2,float);
	for(k=0; k<=NumRes; k++) seq_ins_prob[k]=SeqProb[k];

	NEWP(gpen,(nBlk*Rpts)+2,Int4);
	if(interblkL != -1){
	   if(interblkH == -1){		//  then use interblock gap penalties...
		for(Int4 leng=0,t=1;t<=nblk;t++) {
		    leng+=blklen[t];
		    if(t < nblk){
		      for(Int4 j=0;j < Rpts;j++) {
			Int4 pos=leng + j*length;
		  	i2i[pos]=-interblkL;
		      }
	    	    }
		}
	   } else {
	     assert(nBlk == 1);
	     for(Int4 r=1; r < Rpts; r++){
		NEW(gpen[r],interblkH+4,Int4);
#if 0
	        for(Int4 low=0; low < interblkL; low++){
		   gpen[r][low]=SHRT_MIN+2;  // forbidden lengths...
		}
#endif
		gpen[r][interblkH+1]=SHRT_MIN;
	     }
	   }
	}
}

Int4	*idp_typ::ComputeGREASE(unsigned char *seq, Int4 len, Int4 wind)
// average of grease prob within (closed) wind- nbrhd of each residue
{
        Int4            *grease_penalty;

        NEW(grease_penalty,len+2,Int4);
        for(Int4 i=1;i<=len;i++){
	   grease_penalty[i] = (Int4) floor((pernats*seq_ins_prob[seq[i]])+0.5);
        } return grease_penalty;
}

void	idp_typ::Free( )
{
	free(m2m); free(i2m); free(d2m);
	free(m2i); free(i2i); 
	free(m2d); free(d2d); 
	if(seq_ins_prob) free(seq_ins_prob);
	if(idp_arg) free(idp_arg);
	free(blklen);
	if(gpen){
	    for(Int4 b=nblk*rpts; b >=0; b--) if(gpen[b]) free(gpen[b]);
	    free(gpen);
	} free(prob);
}

Int4	idp_typ::EndBlk(Int4 blk)
{
	assert(blk > 0 && blk <= nblk);
	Int4	end=0;
	for(Int4 b=1; b <= blk; b++) end+=blklen[b];
	return end;
}

Int4	idp_typ::StartBlk(Int4 blk)
{
	assert(blk > 0 && blk <= nblk);
	Int4	start=0;
	for(Int4 b=1; b < blk; b++) start+=blklen[b];
	return (start+1);
}

void	idp_typ::init( )
// ins_open = m2i,ins_extend=i2i,del_open=m2d,
// del_extend=d2d
{
        Int4    i,j,aoi=0,aei=0,aod=0,aed=0;

        NEW(m2i,rpts*length+2,Int4);
        NEW(m2d,rpts*length+2,Int4);
        NEW(i2i,rpts*length+2,Int4);
        NEW(d2d,rpts*length+2,Int4);

        NEW(m2m,rpts*length+2,Int4);
        NEW(i2m,rpts*length+2,Int4);
        NEW(d2m,rpts*length+2,Int4);

        for(i=1;i<=length;i++){
                m2i[i]=(-inso_low)*prob[i]+(-inso_high)*(1.0-prob[i]);
                i2i[i]=(-insx_low)*prob[i]+(-insx_high)*(1.0-prob[i]);
                m2d[i]=-delo;
                // m2d[i]=(-delo_low)*prob[i]+(-delo_high)*(1.0-prob[i]);
                d2d[i]=(-delx_low)*prob[i]+(-delx_high)*(1.0-prob[i]);
                aoi+=m2i[i]; aei+=i2i[i]; 
		aod+=m2d[i]; aed+=d2d[i];
        }
        // Take average of both adjacent columns for insertions
        Int4 t=0;
        for(i=1; i<=length; i++){
	  m2i[i-1] = t;
	  t= (m2i[i]+m2i[i+1])/2; 
	}
        for(t=0,i=1; i<=length; i++){
	  i2i[i-1] = t; t=(i2i[i]+i2i[i+1])/2; 
	}
        if(mode == 'S'){
          aoi/=length; aei/=length; aod/=length; aed/=length;
          for(i=1; i < 8; i++) {
#if 1
		// soften secondary structure predictions near ends.
                m2i[i]=(i*m2i[i]+(8-i)*(-inso_high))/8; 
                m2i[length-i+1]=(i*m2i[length-i+1]+(8-i)*(-inso_high))/8;
                i2i[i]=(i*i2i[i]+(8-i)*(-insx_high))/8; 
                i2i[length-i+1]=(i*i2i[length-i+1]+(8-i)*(-insx_high))/8;
#endif
#if 0
		// soften secondary structure predictions near ends.
                m2i[i]=(i*m2i[i]+(8-i)*aoi)/8; 
                m2i[length-i+1]=(i*m2i[length-i+1]+(8-i)*aoi)/8;
                i2i[i]=(i*i2i[i]+(8-i)*aei)/8; 
                i2i[length-i+1]=(i*i2i[length-i+1]+(8-i)*aei)/8;
#endif

		// lower deletion extension penalty near ends.
		d2d[i]=(i*d2d[i]+(8-i)*(-delx_low))/8;
                d2d[length-i+1]=(i*d2d[length-i+1]+(8-i)*(-delx_low))/8;
          }
        }
        for(j=1;j<rpts;j++){
            for(i=1;i<=length;i++){
                m2i[j*length+i]=m2i[i]; 
		i2i[j*length+i]=i2i[i];
                m2d[j*length+i]=m2d[i]; 
		d2d[j*length+i]=d2d[i];
            }
        }
	// for(j=1;j<=rpts;j++) m2i[j*length]=i2i[j*length]=0; // OLD
	Int4 leng,pos;
	for(leng=0,t=1;t<=nblk;t++) {
	    leng+=blklen[t];
	    for(j=0;j < rpts;j++) {
		pos=leng + j*length;
	  	m2i[pos]=i2i[pos]=0;	// set to zero between blocks???
	    }
	}
	// Calculate i2m, m2m,*m2d, d2m, // later i2d, & d2i.
	double tmp_d,tmp_i;
	for(i=1;i<=length;i++){
	  if(m2i[i] != 0){
		tmp_i= exp((double)i2i[i]/(double)pernats);
		i2m[i] = (Int4) (pernats*log(1.0 - tmp_i));

		tmp_d= exp((double)m2d[i]/(double)pernats);
		tmp_i= exp((double)m2i[i]/(double)pernats);
		// if((tmp_d + tmp_i) > 1.0) m2m[i] = 0;
		if((tmp_d + tmp_i) >= 1.0){
			fprintf(stderr,"pernats = %d\n",pernats);
			fprintf(stderr,"tmp_d = %g; tmp_i = %g (sum = %g >= 1.0)\n",
					tmp_d,tmp_i,tmp_d + tmp_i);
			print_error("idp_typ argument input overflow error");
		} else m2m[i] = (Int4) (pernats*log(1.0 - (tmp_d + tmp_i)));

		tmp_d= exp((double)d2d[i]/(double)pernats);
		d2m[i] = (Int4) (pernats*log(1.0 - tmp_d));
#if 0
fprintf(stderr,"%-3d: %4d %4d %4d %4d %4d %4d %4d \n",
	i,-m2i[i],-i2i[i],-m2d[i],-d2d[i],-i2m[i],-m2m[i],-d2m[i]);
#endif
	  }
	}
	for(j=1;j<rpts;j++){
	   for(i=1;i<=length;i++){
                m2m[j*length+i]=m2m[i]; 
		i2m[j*length+i]=i2m[i];
                d2m[j*length+i]=d2m[i]; 
	   }
	}

}

void	idp_typ::Put(FILE *fp)
{
        fprintf(fp,"POS:  m2i  i2i  m2d  d2d  i2m  m2m  d2m  (COIL REGULAR)\n");
        for(Int4 i=1;i<=length;i++){
            fprintf(fp,"%-3d: %4d %4d %4d %4d %4d %4d %4d (%.2f  %.2f)\n",
                       i,-m2i[i],-i2i[i],-m2d[i],-d2d[i],-i2m[i],-m2m[i],-d2m[i],
                        prob[i],1.0-prob[i]);
        }
}


