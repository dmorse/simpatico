# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

mcMd: 
	cd src; make mcMd

mcMd-mpi: 
	cd src; make mcMd-mpi

ddMd:
	cd src; make ddMd

clean-mcMd:
	cd src; make clean-mcMd

clean-ddMd:
	cd src; make clean-ddMd

clean:
	cd src; make clean

clean-bin:
	cd src; make clean-bin
 
# ==============================================================================
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
veryclean:
	cd html; make clean
	cd src; make veryclean

# ==============================================================================
