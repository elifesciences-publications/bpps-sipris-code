# makefile for libblst

CC = ${CURRENT_CC}

CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a 
IDIR = ../include 
INCDIRS = -I. -I../include -I../libafn 

LIBS = -lstat -lafn ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ 

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = gpxdrop.o my_blastkar.o my_ncbimath.o my_ncbi.o \
	brm_typ.o sfp_typ.o bsb_typ.o hspheap.o brh_typ.o \
	gblast.o 
INC = gpxdrop.h my_blastkar.h my_ncbimath.h my_ncbi.h \
	brm_typ.h sfp_typ.h bsb_typ.h hspheap.h brh_typ.h \
	gblast.h $(STDINC) $(ARCHIVE)

libblst.a: $(OBJS) $(INC)
	ar rcs libblst.a $(OBJS) 
#	ranlib libblst.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

gpxdrop.o: gpxdrop.h $(INC)
my_blastkar.o: my_blastkar.h $(INC)
my_ncbimath.o: my_ncbimath.h $(INC)
brh_typ.o: brh_typ.h $(INC)
my_ncbi.o: my_ncbi.h $(INC)
#my_posit.o: my_posit.h $(INC)
#gpsiblst.o: gpsiblst.h $(INC)
gblast.o: gblast.h $(INC)
hspheap.o: hspheap.h $(INC)
brm_typ.o: brm_typ.h my_blastkar.h $(INC)
bsb_typ.o: bsb_typ.h $(INC)
sfp_typ.o: sfp_typ.h $(INC)
sbp.o: sbp.h $(INC)

clean: 
	@ \rm *.o *.a

