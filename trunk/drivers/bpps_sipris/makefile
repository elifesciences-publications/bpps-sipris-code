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
	../../lib/libstat/libstat.a ../../lib/libalex/libalex.a \
        ../../lib/libtax/libtax.a ../../lib/libhat/libhat.a ../../lib/libgsmo/libgsmo.a \
	../../lib/libomc/libomc.a ../../lib/libvsi/libvsi.a ../../lib/libchn/libchn.a \
	../../lib/libsarp/libsarp.a ../../lib/libpsed/libpsed.a ../../lib/libsprs/libsprs.a

INCDIRS = -I. -I../../lib/include -I../../lib/libafn -I../../lib/libblst -I../../lib/libgsq \
	-I../../lib/libaln -I../../lib/libstat -I../../lib/libtools -I../../lib/libtax \
	-I../../lib/libalex -I../../lib/libhat -I../../lib/libvsi \
	-I../../lib/libpdb  -I../../lib/libomc  -I../../lib/libchn -I../../lib/libsarp/ \
	-I../../lib/libgsmo -I../../lib/libpsed -I../../lib/libsprs

LIBS = -lsprs -lpsed -lsarp -lgsmo -lomc -lchn -lvsi -lhat -ltools -lalex -laln -lgsq -lblst -ltax -lpdb -lafn -lstat ${MATH_LIBRARY}

LIBDIRS = -L../../lib/libafn -L../../lib/libgsq -L../../lib/libaln -L../../lib/libpdb \
	-L../../lib/libblst -L../../lib/libstat -L../../lib/libalex \
	-L../../lib/libtools -L../../lib/libtax -L../../lib/libhat \
	-L../../lib/libvsi -L../../lib/libchn -L../../lib/libsarp -L../../lib/libgsmo \
	-L../../lib/libomc -L../../lib/libpsed -L../../lib/libsprs

STDINC = ../../lib/include/stdinc.h ../../lib/include/residues.h ../../lib/include/sset.h

OBJS = $(ARCHIVE) 

INC = $(ARCHIVE) $(STDINC)

all: mapgaps bpps sipris cma2fa fa2cma

.cc.o:
	${CC} $(CFLAGS) -c $*.cc ${LIBS}

mapgaps: ${INC2} mapgaps.cc fa2cma.o $(ARCHIVE)
	${CC} ${CFLAGS} -o mapgaps mapgaps.cc ${LIBS} fa2cma.o
	@ \mv mapgaps ${BIN}

bpps: $(INC) bpps_main.cc ${OBJS} ${INC} 
	${CC} $(CFLAGS) -o bpps bpps_main.cc ${OBJS} ${LIBS}
	@ \mv bpps ${BIN}

cma2fa: ${INC} cma2fa.o
	${CC} ${CFLAGS} -o cma2fa cma2fa.o ${LIBS}
	@ \mv cma2fa ${BIN}

fa2cma.o: fa2cma.l $(INC) ${OBJS}
	lex fa2cma.l
	\mv -f lex.yy.c fa2cma.c
	${CC} $(CFLAGS) -c fa2cma.c ${OBJS} ${LIBS}

fa2cma: fa2cma.o fa2cma_main.cc $(INC) ${OBJS}
	${CC} $(CFLAGS) -o fa2cma fa2cma_main.cc fa2cma.o ${OBJS} ${LIBS}
	@ \mv fa2cma ${BIN}

sipris: soprise.o $(ARCHIVE) ${OBJS}  ${INC}
	${CC} ${CFLAGS} -o sipris soprise.o ${OBJS} ${LIBS}
	@ \mv sipris ${BIN}

clean: 
	@ \rm *.o 

