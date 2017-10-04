# makefile for scan
CC = ${CURRENT_CC}
CFLAGS = $(INCDIRS) $(LIBDIRS)
BIN = ../bin
ARCHIVES = ../libafn/libafn.a ../libpdb/libpdb.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libblst/libblst.a ../libtools/libtools.a ../libstat/libstat.a \
	../libalex/libalex.a ../libtax/libtax.a ../libchn/libchn.a 

#CAUTION: make sure that libvsi is not included as both libvsi and libchn archive lex/yacc files.

INCDIRS = -I. -I../include -I../libafn -I../libblst -I../libgsq -I../libaln \
	-I../libpdb -I../libstat -I../libtools -I../libtax \
	-I../libalex -I../libchn 

LIBS = -lchn -ltax -ltools -lalex -laln -lgsq -lblst -lpdb -lafn -lstat ${MATH_LIBRARY}

LIBDIRS = -L../libafn -L../libgsq -L../libaln -L../libpdb -L../libblst \
	-L../libstat -L../libalex -L../libtools -L../libtax \
	-L../libchn 

STDINC = ../include/stdinc.h ../include/residues.h ../include/blosum62.h

OBJS = mcs_arg.o mcs_typ.o mcs_put.o mcs_cdd.o mcs_sample.o mcs_pttrn.o mcs_init.o mcs_rtrn.o mcs_updown.o mcs_partition.o \
	mcs_debug.o mcsBPPS.o hsi_typ.o omc_operate.o omc_fuse.o omc_addleaf.o omc_up.o omc_down.o omc_hihmm.o \
	omc_debug.o omc_tweak.o omc_typ.o omc_simulate.o omc_init.o omc_run.o omc_junk.o omc_insert.o \
	omc_cmh.o mcs_junk.o mad_typ.o

INC = mcs_typ.h omc_typ.h hsi_typ.h mad_typ.h $(ARCHIVES) $(STDINC)

libomc.a: $(OBJS) $(INC) ${ARCHIVES}
	\rm -f libomc.a
	ar rcs libomc.a $(OBJS)

.cc.o:	$(INC)
	${CC} ${CFLAGS} -c $*.cc

hsi_typ.o: ${INC}
mad_typ.o: ${INC}
mcs_arg.o: ${INC}
mcs_debug.o: ${INC}
mcs_typ.o: ${INC}
mcs_put.o: ${INC}
mcs_cdd.o: ${INC}
mcs_sample.o: ${INC}
mcs_pttrn.o: ${INC}
mcs_junk.o: ${INC}
mcs_rtrn.o: ${INC}
mcs_updown.o: ${INC}
mcs_init.o: ${INC}
mcs_partition.o: ${INC}
omc_typ.o: ${INC}
omc_init.o: ${INC}
omc_simulate.o: ${INC}
omc_run.o: ${INC}
omc_debug.o: ${INC}
omc_cmh.o: ${INC}
omc_hihmm.o: ${INC}
omc_junk.o: ${INC}
omc_insert.o: ${INC}
omc_tweak.o: ${INC}
omc_operate.o: ${INC}
omc_down.o: ${INC}
omc_addleaf.o: ${INC}
omc_up.o: ${INC}
omc_fuse.o: ${INC}

clean: 
	@ \rm *.o *.a

