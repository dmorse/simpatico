BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/sources.mk

TEST=$(BLD_DIR)/ddMd/tests/Test

.SECONDARY:$(ddMd_tests_OBJS)

all: $(ddMd_tests_EXES) 

run: $(TEST)
	$(MPIRUN) 6 $(TEST) $(SRC_DIR)/ddMd/tests/ 
	cd $(SRC_DIR)/ddMd/tests; $(MAKE) clean-outputs

quiet: $(TEST)
	$(MPIRUN) 6 $(TEST) $(SRC_DIR)/ddMd/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in ddMd/tests/log" > count
	@cat count
	@cd $(SRC_DIR)/ddMd/tests; $(MAKE) clean-outputs

clean:
	rm -f $(ddMd_tests_EXES)
	rm -f $(ddMd_tests_OBJS)
	rm -f $(ddMd_tests_OBJS:.o=.d)
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd chemistry; $(MAKE) clean
	cd communicate; $(MAKE) clean
	cd configIos; $(MAKE) clean
	cd neighbor; $(MAKE) clean
	cd potentials; $(MAKE) clean
	cd simulation; $(MAKE) clean
	cd storage; $(MAKE) clean
else
	cd $(SRC_DIR)/ddMd/tests; $(MAKE) clean-outputs
endif

clean-outputs:
	@cd chemistry; $(MAKE) clean-outputs
	@cd configIos; $(MAKE) clean-outputs
        
-include $(ddMd_tests_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

