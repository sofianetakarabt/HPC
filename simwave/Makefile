SHELL = /bin/bash
VERBOSE  := @
DEBUG    := 0
OPENMP   := 1
MACHINE  := $(shell uname -s)

SRCDIR   := src
OBJDIR   := src
BINDIR   := bin
LIBDIR   := lib
INCDIR   := include
DATADIR  := data

LIB=-lm -L/usr/local/cuda/lib64/ #-lcuda -lcudart -lstdc++

VISU     := ${VISU_BINDIR}/xmovie

CC       := mpicc
NVCC=nvcc
LD       := ${CC}
SRCS_C   := ${wildcard ${SRCDIR}/*.cu}
OBJS_CU	 := ${SRCS_C:%.cu=%.o}

SRCS     := ${wildcard ${SRCDIR}/*.c}
OBJS     := ${SRCS:%.c=%.o}

EXE      := ${BINDIR}/simwave

CPPFLAGS := -I${INCDIR}

LOCAL_PATH=$(HOME)/local/
INC_DIR += -I$(LOCAL_PATH)include/
LIB_DIR += -L$(LOCAL_PATH)lib/

# librt is not defined on Darwin and MacOS
ifeq (${MACHINE},Darwin)
LDLIBS   := -lm
else
LDLIBS   += -lm -lrt
endif

ifeq (${DEBUG},1)
CFLAGS   := -g -D_DEBUG
LDFLAGS  :=
CUDA_FLAGS = --generate-code arch=compute_20,code=sm_20 --generate-code arch=compute_35,code=sm_35 -g -G #-DCUDA_DEBUG
else
CFLAGS   := -Wall -O3 -msse3  -funroll-loops -ftree-vectorize  -mfpmath=sse
CUDA_FLAGS = --generate-code arch=compute_20,code=sm_20 --generate-code arch=compute_35,code=sm_35 -O3 #-DCUDA_DEBUG
ifeq (${MACHINE},Darwin)
CFLAGS   += -std=c99
else
CFLAGS   += -std=c99 -D_POSIX_C_SOURCE=200112L -D_GNU_SOURCE
endif
LDFLAGS  :=
ifeq (${OPENMP},1)
CFLAGS   += -fopenmp
LDLIBS   += -lgomp
else
CFLAGS   += -Wno-unknown-pragmas
endif
endif

all: ${EXE}

${EXE}: ${OBJS} ${OBJS_CU} ${wildcard ${INCDIR}/*.h}
	@echo "**************************************************************************",${OBJS_CU}
	@echo "**************************************************************************",${OBJS}
	${VERBOSE}${CC} ${LDFLAGS} -o ${EXE} ${OBJS} ${LDLIBS} ${LIB}


${OBJDIR}/%.o: ${SRCDIR}/%.c
	${VERBOSE}${CC} ${CPPFLAGS} ${CFLAGS} -o $@ -c $<

#${OBJDIR}/gpu_%.o: ${SRCDIR}/gpu_%.cu
#	${VERBOSE}${NVCC} -c $+  ${CUDA_FLAGS} ${CPPFLAGS} ${LIB} -o $@

exec:
	./bin/simwave -v --grid 100,100,100 --iter 100  --epsilon 1e-17 --opt --bloc 32,4,1
	./bin/simwave -v --grid 400,400,400 --iter 100  --epsilon 1e-17 --opt --bloc 32,4,1
	./bin/simwave -v --grid 400,400,400 --iter 100  --epsilon 1e-17 --opt --bloc 32,4,1
source:
	${VERBOSE}gnuplot ${DATADIR}/source.gp
	${VERBOSE}evince  ${DATADIR}/source.eps &

clean:
	${VERBOSE}rm -rf ${SRCDIR}/*.o
	${VERBOSE}find . -name "*~" -exec rm -rf '{}' \;
	${VERBOSE}find . -name "#*" -exec rm -rf '{}' \;
	${VERBOSE}find . -name "*.log" -exec rm -rf '{}' \;

distclean: clean
	${VERBOSE}rm -rf ${EXE}
	${VERBOSE}rm -rf ${DATADIR}/*.bin

-include make/display.mk
-include make/doc.mk
-include make/release.mk
-include make/tex.mk

.PHONY: all clean distclean
# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .c
