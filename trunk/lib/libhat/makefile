# makefile for scan
CC = ${CURRENT_CC}
CFLAGS = $(INCDIRS) $(LIBDIRS)
ARCHIVES = ../libafn/libafn.a ../libpdb/libpdb.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libblst/libblst.a ../libtools/libtools.a \
	../libstat/libstat.a ../libalex/libalex.a ../libtax/libtax.a 

INCDIRS = -I. -I../include -I../libafn -I../libblst -I../libgsq -I../libaln \
	-I../libpdb -I../libstat -I../libtools -I../libtax \
	-I../libalex 

LIBS = -ltax -ltools -lalex -laln -lgsq -lblst -lpdb -lafn -lstat ${MATH_LIBRARY}
LIBDIRS = -L../libafn -L../libgsq -L../libaln -L../libpdb -L../libblst \
	-L../libstat -L../libalex -L../libtools -L../libtax 

STDINC = ../include/stdinc.h ../include/residues.h

OBJS = patch_sap.o sph_typ.o hat_typ.o cma_gblastpgp.o curated_srch.o  matblast_srch.o setup_gapmaps.o gapmap_srch.o mgs_typ.o dft_typ.o convert_template.o

INC = mgs_typ.h hat_typ.h sph_typ.h patch_sap.h dft_typ.h $(ARCHIVES) $(STDINC)

libhat.a: $(OBJS) $(INC)
	ar rcs libhat.a $(OBJS)

.cc.o:	$(INC)
	${CC} ${CFLAGS} -c $*.cc

patch_sap.o: $(INC)
sph_typ.o: $(INC)
hat_typ.o: $(INC)
dft_typ.o: $(INC)
cma_gblastpgp.o: $(INC)
curated_srch.o: $(INC)
gapmap_srch.o: $(INC)
mgs_typ.o: $(INC)
setup_gapmaps.o: $(INC)
matblast_srch.o: $(INC)
convert_template.o: $(INC)

clean: 
	@ \rm *.o *.a

