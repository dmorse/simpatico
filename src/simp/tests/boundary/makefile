BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/boundary/sources.mk

all: $(simp_tests_boundary_EXES) 

clean:
	rm -f $(simp_tests_boundary_EXES) 
	rm -f $(simp_tests_boundary_OBJS) 
	rm -f $(simp_tests_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_tests_boundary_OBJS:.o=.d)

-include $(simp_tests_boundary_OBJS:.o=.d)

