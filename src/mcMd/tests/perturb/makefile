SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include sources.mk

TEST=$(BLD_DIR)/mcMd/tests/perturb/Test

.Secondary: $(mcMd_tests_perturb_OBJS) 

all: $(mcMd_tests_perturb_EXES) 

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_perturb_EXES) 
	rm -f $(mcMd_tests_perturb_OBJS) 
	rm -f $(mcMd_tests_perturb_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_tests_perturb_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(simp_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(simp_LIB) $(util_LIB)

-include $(mcMd_tests_perturb_OBJS:.o=.d)

