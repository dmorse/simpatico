# ==========================================================================
.PHONY: help util simp mcMd mcMd-mpi ddMd mdPp clean clean-tests veryclean

# -----------------------------------------------------------------------
# Generate help dialog 

help: ## This help dialog.
	@IFS=$$'\n' ; \
        help_lines=(`fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//'`); \
        for help_line in $${help_lines[@]}; do \
           IFS=$$'#' ; \
           help_split=($$help_line) ; \
           help_command=`echo $${help_split[0]} | sed -e 's/^ *//' -e 's/ *$$//'` ; \
           help_info=`echo $${help_split[2]} | sed -e 's/^ *//' -e 's/ *$$//'` ; \
           printf "%-18s %s\n" $$help_command $$help_info ; \
        done

# -----------------------------------------------------------------------
# Main makefile targets

util:        ## Compile src/util directory (MPI state unchanged)
	cd util; $(MAKE) all

simp:        ## Compile src/util, src/simp directories (MPI state unchanged)
	cd util; $(MAKE) all
	cd simp; $(MAKE) all

mcMd:        ## Build serial versions of mdSim and mcSim (disables MPI)
	./configure -m0
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd mcMd; $(MAKE) all

mcMd-mpi:    ## Build parallel versions of mdSim and mcSim (enables MPI)
	./configure -m1
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd mcMd; $(MAKE) all
 
ddMd:        ## Build ddSim parallel MD program (enables MPI) 
	./configure -m1
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd ddMd; $(MAKE) all

mdPp:        ## Build mdPp analysis program (disables MPI)
	./configure -m0
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd mdPp; $(MAKE) all

# -----------------------------------------------------------------------
# Clean-up targets

clean:       ## Remove all intermediate build files (*.o, *.d, *.a files)
	cd util; $(MAKE) clean
	cd simp; $(MAKE) clean
	cd mcMd; $(MAKE) clean
	cd mdPp; $(MAKE) clean
	cd ddMd; $(MAKE) clean

clean-tests: ## Clean unit test directories (remove files created by tests)
	cd util/tests; $(MAKE) clean
	cd simp/tests; $(MAKE) clean
	cd mcMd/tests; $(MAKE) clean
	cd mdPp/tests; $(MAKE) clean
	cd ddMd/tests; $(MAKE) clean

veryclean:  ## Remove all generated files, revert to as-distributed state
	cd util; $(MAKE) veryclean
	cd simp; $(MAKE) veryclean
	cd mcMd; $(MAKE) veryclean
	cd ddMd; $(MAKE) veryclean
	cd mdPp; $(MAKE) veryclean
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f simp/user/sources.mk
	rm -f mcMd/user/sources.mk
	rm -f ddMd/user/sources.mk
	rm -f mdPp/user/sources.mk
	rm -f mcMd/modules/sources.mk
	rm -f mcMd/mcSim.cpp
	rm -f mcMd/mdSim.cpp
endif
	rm -f util/config.mk
	rm -f simp/config.mk
	rm -f mcMd/config.mk
	rm -f ddMd/config.mk
	rm -f mdPp/config.mk
	rm -f config.mk

# =========================================================================
