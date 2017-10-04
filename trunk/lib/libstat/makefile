# makefile for scan

STDINC = 

OBJS = libstat.a

INC = $(ARCHIVES) $(STDINC)

libstat.a: lgamma.o 
	ar rcs libstat.a lgamma.o 

lgamma.o: 
	${CURRENT_CC} -c lgamma.c ${MATH_LIBRARY}

#lgamma.o: 
#	gcc -O -c lgamma.c -lm

clean: 
	@ \rm *.[oa] 

