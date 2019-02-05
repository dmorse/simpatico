BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/chemistry/sources.mk

#ifdef UTIL_MPI
#TEST=$(BLD_DIR)/mcMd/tests/chemistry/MpiTest
#else
#endif
TEST=$(BLD_DIR)/mcMd/tests/chemistry/Test

.Secondary: $(mcMd_tests_chemistry_EXES) 

all: $(mcMd_tests_chemistry_EXES) 

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_chemistry_OBJS)
	rm -f $(mcMd_tests_chemistry_OBJS:.o=.d)
	rm -f $(mcMd_tests_chemistry_OBJS:.o=)
	$(MAKE) clean-outputs

clean-outputs:
	@rm -f tmp/*

-include $(mcMd_tests_chemistry_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

