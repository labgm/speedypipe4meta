# Makefile for fm-index

# VPATH and shell configuration
top_srcdir = ..
srcdir     = .

SHELL      = /bin/sh

# sources
OBJS				= divsufsort.o 
TARGET				= libdivsufsort.a
MAKEFILE			= Makefile

# options
CC        = /home/lucas/anaconda3/bin/x86_64-conda_cos6-linux-gnu-cc
CFLAGS    = -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/lucas/anaconda3/include -pthread -fPIC
CPPFLAGS  = -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/lucas/anaconda3/include
MPILIBS   = 
AR        = /home/lucas/anaconda3/bin/x86_64-conda_cos6-linux-gnu-ar rc
RANLIB    = /home/lucas/anaconda3/bin/x86_64-conda_cos6-linux-gnu-ranlib
INSTALL   = /usr/bin/install -c

#CC					= gcc 
#CFLAGS  	 		= -g -O3 -fomit-frame-pointer -fstrict-aliasing 
#CPPFLAGS			= -Wall 
#AR       = /usr/bin/ar rc
#RANLIB   = ranlib


# beautification magic stolen from git (added within hmmer source)
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_AR      = @echo '    ' AR $@;
endif


# targets
all: $(TARGET) 
$(TARGET): $(OBJS)

.c.o:
	${QUIET_CC}${CC} -I. ${CFLAGS} ${CPPFLAGS} -o $@ -c $<


libdivsufsort.a: $(OBJS)
	${QUIET_AR}${AR} libdivsufsort.a $(OBJS)
	@${RANLIB} libdivsufsort.a
	@chmod 644 libdivsufsort.a


clean:
	$(RM) $(TARGET) $(OBJS)

distclean:
	$(RM) $(TARGET) $(OBJS) $(MAKEFILE) divsufsort.h

# dependencies
$(OBJS): divsufsort.h $(MAKEFILE)
