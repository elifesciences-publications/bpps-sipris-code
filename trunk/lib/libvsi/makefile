# makefile for chn_vsi: Aligned Sequence Statistical Evaluation Tool 

CC = ${CURRENT_CC}
#CC = g++ -static -g -O3 
CFLAGS = $(INCDIRS) $(LIBDIRS)

YACC = yacc -vd
ARCHIVE = ../libafn/libafn.a ../libpdb/libpdb.a ../libstat/libstat.a 
INCDIRS = -I. -I../include -I../libafn -I../libpdb 
LIBS = -lpdb -lafn -lstat ${MATH_LIBRARY} 
LIBDIRS = -L../libafn -L../libpdb -L../libstat

OBJS = y.tab.o lex.yy.o vsi_pdb.o aromatic.o 

INC = chn_vsi.h aromatic.h vsi_pdb.h $(ARCHIVE) 

#all: chn_vsi libvsi.a

libvsi.a: $(OBJS) $(INC)
	ar rcs libvsi.a $(OBJS)

.cc.o:
	${CC} ${CFLAGS} -c $*.cc

#BIN = ../../bin
#chn_vsi: ${ARCHIVE} ${OBJS} ${INC}
#	${CC} ${CFLAGS} -o chn_vsi ${OBJS} ${LIBS}
#	@ mv -f chn_vsi ${BIN}

y.tab.o: y.tab.c ${ARCHIVE}
	${CC} ${CFLAGS} -c y.tab.c 
#	@ \rm -f y.tab.c

y.tab.c: chn_vsi.y ${INC} ${ARCHIVE}
	${YACC} chn_vsi.y

lex.yy.o: y.tab.o chn_vsi.l ${INC} ${ARCHIVE}
	flex chn_vsi.l 
	${CC} ${CFLAGS} -c lex.yy.c -lfl
#	@ \rm -f lex.yy.c

vsi_pdb.o: ${INC} 
aromatic.o: ${INC} 

clean: 
	@ \rm *.o *.a lex.yy.c y.tab.c y.tab.h


