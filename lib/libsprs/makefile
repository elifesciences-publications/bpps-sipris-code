# makefile for gibbs_prop

BIN = ../../bin
CC = ${CURRENT_CC}
#CC = lint 
#CC = CC -signed 
#CC = CC -signed -ansi -fullwarn 
#CC = CC -sopt -signed -ansi -fullwarn 
#CC = CC -wlint,-pq -signed -ansi -fullwarn 
# for profiles 
#CFLAGS = -p -g $(LIBDIRS) $(INCDIRS) 
CFLAGS = -g $(LIBDIRS) $(INCDIRS) 
#CFLAGS = -O ${INCDIRS} $(LIBDIRS) 
#CFLAGS = -O2 $(LIBDIRS) ${INCDIRS} 
ARCHIVE = ../../lib/libafn/libafn.a ../../lib/libpdb/libpdb.a  ../../lib/libgsq/libgsq.a \
	../../lib/libaln/libaln.a ../../lib/libblst/libblst.a ../../lib/libtools/libtools.a \
	../../lib/libstat/libstat.a ../../lib/libalex/libalex.a ../../lib/libhmm/libhmm.a \
        ../../lib/libtax/libtax.a ../../lib/libhat/libhat.a ../../lib/libgsmo/libgsmo.a \
	../../lib/libomc/libomc.a ../../lib/libvsi/libvsi.a ../../lib/libchn/libchn.a \
	../../lib/libsarp/libsarp.a

#CAUTION: make sure that libvsi is not included as both libvsi and libchn archive lex/yacc files

INCDIRS = -I. -I../../lib/include -I../../lib/libafn -I../../lib/libblst -I../../lib/libgsq \
	-I../../lib/libaln \
        -I../../lib/libstat -I../../lib/libtools -I../../lib/libtax -I../../lib/libalex \
        -I../../lib/libhmm -I../../lib/libhat -I../../lib/libvsi -I../../lib/libpdb  \
	-I../../lib/libomc  -I../../lib/libchn -I../../lib/libsarp/ -I../../lib/libgsmo/ 

LIBS = -lsarp -lgsmo -lomc -lchn -lvsi -lhat -ltools -lhmm -lalex -laln -lgsq -lblst -ltax -lpdb -lafn -lstat ${MATH_LIBRARY}

LIBDIRS = -L../../lib/libafn -L../../lib/libgsq -L../../lib/libaln -L../../lib/libpdb \
	-L../../lib/libblst -L../../lib/libstat -L../../lib/libalex \
	-L../../lib/libtools -L../../lib/libhmm -L../../lib/libtax -L../../lib/libhat \
	-L../../lib/libvsi -L../../lib/libchn -L../../lib/libsarp -L../../lib/libgsmo \
	-L../../lib/libomc

STDINC = ../../lib/include/stdinc.h ../../lib/include/residues.h

OBJS = sch_typ.o scl_dca.o scl_typ.o scm_typ.o $(ARCHIVES) 

INC = scl_typ.h $(ARCHIVES) $(STDINC)

libsprs.a: $(OBJS) $(INC) ${ARCHIVES}
	\rm -f libsprs.a
	ar rcs libsprs.a $(OBJS)

.cc.o: $(INC)
	${CC} $(CFLAGS) -c $*.cc ${LIBS}

sch_typ.o: ${INC}
scl_dca.o: ${INC}
scl_typ.o: ${INC}
scm_typ.o: ${INC}

clean: 
	@ \rm *.[oa] 

