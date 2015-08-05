CC := gcc
CXX := g++
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
objs := pairdist.cpp pairdist.pyx
autowrap:
	autowrap --out pairdist.pyx pairdist_h.pxd

clean:
	@rm -f $(objs)

uninstall:
	GREP = $(pip list)
	ifeq ($(GREP), )
	    @echo "Nothing was found for current build"
	else
	    @echo "***Found string in Findings.txt***"
	endif

install:
	CC=$(CC) CXX=$(CXX) pip install $(ROOT_DIR)

check:
	@echo $(ROOT_DIR)
