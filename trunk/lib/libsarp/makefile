# makefile for sarp library
CC = ${CURRENT_CC}
CFLAGS = $(INCDIRS) $(LIBDIRS)
ARCHIVES = ../libafn/libafn.a ../libpdb/libpdb.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libchn/libchn.a ../libblst/libblst.a ../libtools/libtools.a \
	../libstat/libstat.a ../libalex/libalex.a ../libtax/libtax.a \
	../libvsi/libvsi.a 

INCDIRS = -I. -I../include -I../libafn -I../libblst -I../libgsq -I../libaln \
	-I../libvsi -I../libpdb -I../libstat -I../libtools \
	-I../libtax -I../libchn -I../libalex

LIBS = -lchn -ltax -ltools -lalex -lvsi -laln -lgsq -lblst -lpdb -lafn -lstat ${MATH_LIBRARY} -lfl

LIBDIRS = -L../libafn -L../libgsq -L../libaln -L../libpdb -L../libblst \
	-L../libstat -L../libalex -L../libtools -L../libtax \
	-L../libvsi -L../libchn

STDINC = ../include/stdinc.h ../include/residues.h 

INC = $(ARCHIVES) $(STDINC) adh_typ.h psc_typ.h rai_typ.h rih_typ.h hsc_typ.h esc_typ.h mps_typ.h sas_typ.h

OBJS = psc_typ.o adh_typ.o psc_pttrn.o rai_typ.o rih_typ.o hsc_typ.o p2c_CaScrs.o p2c_typ.o \
	p2c_debug.o esc_typ.o mps_typ.o sas_typ.o

libsarp.a: $(OBJS) $(INC)
	ar rcs libsarp.a $(OBJS)

.cc.o: $(INC)
	${CC} ${CFLAGS} -c $*.cc

adh_typ.o: adh_typ.cc $(INC)
psc_typ.o: psc_typ.cc $(INC)
p2c_typ.o: p2c_typ.cc $(INC)
p2c_debug.o: p2c_debug.cc $(INC)
hsc_typ.o: hsc_typ.cc $(INC)
sas_typ.o: sas_typ.cc $(INC)
esc_typ.o: esc_typ.cc $(INC)
mps_typ.o: mps_typ.cc $(INC)
psc_pttrn.o: psc_pttrn.cc $(INC)
#p2c_ca_scores.o: p2c_ca_scores.cc $(INC)
p2c_CaScrs.o: p2c_CaScrs.cc $(INC)

clean: 
	@ \rm *.o  *.a


