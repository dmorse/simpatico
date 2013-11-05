SRC_DIR_REL =../../../src
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/inter/config.mk
include $(SRC_DIR)/mcMd/config.mk
include $(TESTS_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include sources.mk

all: $(tests_mcMd_perturb_OBJS) Test

clean:
	rm -f $(tests_mcMd_perturb_OBJS) Test $(tests_mcMd_perturb_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_mcMd_perturb_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include $(tests_mcMd_perturb_OBJS:.o=.d)

