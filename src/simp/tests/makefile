BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/tests/sources.mk

all: $(simp_tests_EXES)

TEST=$(BLD_DIR)/simp/tests/Test

run: $(TEST)
	$(TEST) $(SRC_DIR)/simp/tests/ 
	cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs

quiet: $(TEST)
	$(TEST) $(SRC_DIR)/simp/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in simp/tests/log" > count
	@cat count
	@cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs

clean:
	rm -f $(simp_tests_EXES)
	rm -f $(simp_tests_OBJS:.o=.d) 
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd interaction; $(MAKE) clean
	cd species; $(MAKE) clean
else
	cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs
endif

clean-outputs:
	@cd species; $(MAKE) clean-outputs

-include $(simp_tests_OBJS:.o=.d)
