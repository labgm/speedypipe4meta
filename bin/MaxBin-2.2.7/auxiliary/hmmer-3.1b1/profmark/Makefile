top_srcdir = ..
srcdir     = .


CC       = /home/lucas/anaconda3/bin/x86_64-conda_cos6-linux-gnu-cc
CFLAGS   = -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/lucas/anaconda3/include -pthread -fPIC
SIMDFLAGS= -msse2
CPPFLAGS = -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/lucas/anaconda3/include
LDFLAGS  = -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/lucas/anaconda3/lib -Wl,-rpath-link,/home/lucas/anaconda3/lib -L/home/lucas/anaconda3/lib
DEFS     = -DHAVE_CONFIG_H
LIBS     = -lhmmer -leasel    -lm
MPILIBS  = 
IMPLDIR  = impl_sse

ESLDIR    = easel

ESLINC   = -I../${ESLDIR} -I${top_srcdir}/easel
SRCINC   = -I../src   -I${top_srcdir}/src

PROGS    = create-profmark\
	   rocplot

PROGOBJS  =\
	create-profmark.o\
	rocplot.o

# beautification magic stolen from git 
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR $$subdir; \
		        ${MAKE} -C $$subdir
endif

.PHONY: all dev distclean clean

all:    ${PROGS}
dev:    ${PROGS}

${PROGS}: %: %.o  ../${ESLDIR}/libeasel.a ../src/libhmmer.a 
	${QUIET_GEN}${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} -L../${ESLDIR} -L../src -o $@ $@.o ${LIBS}

${PROGOBJS}: ../src/hmmer.h ../src/p7_config.h

.c.o: 
	${QUIET_CC}${CC} ${ESLINC} ${SRCINC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -o $@ -c $<		

distclean: clean
	-rm -f Makefile

clean:
	-rm -f *.o *~ ${PROGS} 
	-rm -f *.gcno

################################################################
# HMMER - Biological sequence analysis with profile HMMs
# Version 3.1b1; May 2013
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# HMMER is distributed under the terms of the GNU General Public License
# (GPLv3). See the LICENSE file for details.
# 
# SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/branches/3.1/profmark/Makefile.in $
# SVN $Id: Makefile.in 3844 2012-02-07 21:51:20Z eddys $
################################################################

