# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

mcMd: 
	cd src; make mcMd

mcMd_mpi: 
	cd src; make mcMd_mpi

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
 
veryclean:
	cd src; make veryclean

# ==============================================================================
