SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/simp/config.mk
include $(SRC_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include sources.mk

all: $(mcMd_tests_perturb_OBJS) Test

clean:
	rm -f $(mcMd_tests_perturb_OBJS) Test $(tests_mcMd_perturb_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_tests_perturb_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(simp_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(simp_LIB) $(util_LIB)

-include $(tests_mcMd_perturb_OBJS:.o=.d)

