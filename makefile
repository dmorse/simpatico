include src/compiler.mk
# ==============================================================================
.PHONY: mcMd ddMd mcMd-mpi \
        clean-serial clean-parallel clean clean-bin \
        html clean-html 

all:
	cd build/serial; make mcMd
	cd build/parallel; make ddMd
	cd build/parallel; make mcMd-mpi

mcMd:
	cd build/serial; make mcMd

mcMd-mpi: 
	cd build/parallel; make mcMd-mpi

ddMd:
	cd build/parallel; make ddMd

# ==============================================================================
# Clean targets

clean-serial:
	cd build/serial; make clean

clean-parallel:
	cd build/parallel; make clean

clean:
	cd src; make clean
	cd build/serial; make clean
	cd build/parallel; make clean

clean-bin:
	cd src; make clean-bin

# Remove all generated files, including those created by the setup script.
veryclean:
	cd html; make clean
	cd src; make veryclean
	cd build/serial; make clean
	cd build/parallel; make clean

# ==============================================================================
# HTML Documentation
 
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
