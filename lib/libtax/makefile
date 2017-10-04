# makefile for sqtax
CC = ${CURRENT_CC}
CFLAGS = $(INCDIRS) $(LIBDIRS)
BIN = ../../bin
ARCHIVES = ../libafn/libafn.a  ../libgsq/libgsq.a ../libaln/libaln.a \
	../libblst/libblst.a ../libtools/libtools.a \
	../libstat/libstat.a
INCDIRS = -I. -I../include -I../libafn -I../libgsq -I../libaln \
	-I../libtools -I../libblst -I../libstat
LIBS = -ltools -lblst -laln -lgsq -lafn -lstat ${MATH_LIBRARY} -lfl
LIBDIRS = -L. -L../libafn -L../libgsq -L../libaln \
	-L../libblst -L../libtools -L../libstat
STDINC = ../include/stdinc.h ../include/residues.h 

INC = txs_typ.h n2a_typ.h tax_typ.h bsa_typ.h rsq_typ.h $(ARCHIVES) $(STDINC)
OBJS = txs_typ.o n2a_typ.o tax_typ.o bsa_typ.o rsq_typ.o 

.cc.o:
	${CC} ${CFLAGS} -c $*.cc

libtax.a: $(OBJS) $(INC)
	ar rcs libtax.a $(OBJS)

txs_typ.o: txs_typ.cc $(INC)
tax_typ.o: tax_typ.cc $(INC)
n2a_typ.o: n2a_typ.cc $(INC)
bsa_typ.o: bsa_typ.cc $(INC)
rsq_typ.o: rsq_typ.cc $(INC)

nrtax:  libtax.a nrtax.l $(INC)
	lex nrtax.l
	\mv lex.yy.c lex.yy.cc
	${CC} ${CFLAGS} -o nrtax lex.yy.cc -ltax ${LIBS} 
	\rm lex.yy.cc 
	\mv nrtax $BIN/

sqtax:  libtax.a sqtax.l $(INC)
	lex sqtax.l
	\mv lex.yy.c lex.yy.cc
	${CC} ${CFLAGS} -o sqtax lex.yy.cc -ltax ${LIBS} 
	\rm lex.yy.cc 
	\mv sqtax $BIN/

clean: 
	@ \rm *.o *.a


