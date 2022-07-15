BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/mdPp/tests/sources.mk

TEST=$(BLD_DIR)/mdPp/tests/Test

.Secondary: $(mdPp_tests_OBJS)

all: $(mdPp_tests_EXES)

run: $(TEST)
	$(TEST) $(SRC_DIR)/mdPp/tests/ 

quiet: $(TEST)
	$(TEST) $(SRC_DIR)/mdPp/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in mdPp/tests/log" > count
	@cat count

clean:
	rm -f $(mdPp_tests_OBJS) 
	rm -f $(mdPp_tests_OBJS:.o=.d)
	rm -f $(mdPp_tests_OBJS:.o=)
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd chemistry; $(MAKE) clean
	cd config; $(MAKE) clean
	cd neighbor; $(MAKE) clean
	cd storage; $(MAKE) clean
	cd processor; $(MAKE) clean
	rm -f */*.o */*.d
#	$(MAKE) clean-outputs
#else
#	cd $(SRC_DIR)/mdPp/tests; $(MAKE) clean outputs
endif

#clean-outputs:
#	cd chemistry; $(MAKE) clean-outputs
        
-include $(mdPp_tests_OBJS:.o=.d)
-include $(mdPp_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

