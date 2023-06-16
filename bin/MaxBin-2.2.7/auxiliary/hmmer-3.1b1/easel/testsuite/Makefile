top_srcdir = ../..
srcdir     = .
VPATH      =  ../..
SHELL      = /bin/sh

CC         = /home/lucas/anaconda3/bin/x86_64-conda_cos6-linux-gnu-cc
CFLAGS     = -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/lucas/anaconda3/include -pthread -fPIC
CPPFLAGS   = -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/lucas/anaconda3/include
LDFLAGS    = -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/lucas/anaconda3/lib -Wl,-rpath-link,/home/lucas/anaconda3/lib -L/home/lucas/anaconda3/lib
SIMDFLAGS  = -msse2
DEFS       = -DHAVE_CONFIG_H
LIBS       =   -lm
MPILIBS    = 

SQCLEVEL = 2

QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR  $$subdir; \
		        ${MAKE} -C $$subdir
endif

.PHONY: check dcheck clean distclean

check:
	@echo 
	@echo Running Easel test suite...
	@echo
	@${srcdir}/../devkit/sqc ${SQCLEVEL} ${srcdir}/testsuite.sqc .. ${srcdir}/..

dcheck: 
	@echo    Running developer testsuite...
	@echo -n Running unit tests...
	@${srcdir}/../devkit/sqc ${SQCLEVEL} ${srcdir}/testsuite.sqc .. ${srcdir}/..  >  utest.report
	@echo    done. [See utest.report]
	@echo -n Testing that all drivers compile...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./driver_report.pl       >  driver.report
	@echo    done. [See driver.report]
	@echo -n Measuring test code coverage...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./coverage_report.pl     >  coverage.report
	@echo    done. [See coverage.report]
	@echo -n Looking for memory leaks...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./valgrind_report.pl     >  valgrind.report
	@echo    done. [See valgrind.report]
	@cat driver.report coverage.report valgrind.report utest.report > dcheck.report
	@echo
	@echo    Finished developer testsuite.
	@echo    For a complete report, see dcheck.report.

clean:
	-rm -f *.gcno *.gcda *.gcov
	-rm -f esltmp??????
	-rm -f *.o *~ 
	-rm -f driver.report coverage.report valgrind.report utest.report
	-rm -f ${PROGS} 

distclean: clean
	rm Makefile

################################################################
# Easel - a library of C functions for biological sequence analysis
# Version h3.1b1; May 2013
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# Easel is distributed under the Janelia Farm Software License, a BSD
# license. See the LICENSE file for more details.
################################################################
