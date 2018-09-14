QUEST_DIR = .

include make.inc

#all: example_mkl

all : example_

example_: liblapack libblas libdqmc
	(cd applications; $(MAKE))

example_mkl: libdqmc
	$(MAKE) -C applications	

libblas:
	(cd libs/BLAS; $(MAKE))

liblapack:
	(cd libs/LAPACK; $(MAKE))

libdqmc:
	$(MAKE) -C src

clean:
	(cd libs/BLAS; $(MAKE) clean)
	(cd libs/LAPACK; $(MAKE) clean)
	(cd src; $(MAKE) clean)
	(cd applications; $(MAKE) clean)
	(rm -f $(DQMCLIB))

