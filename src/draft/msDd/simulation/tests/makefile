SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/ddMd/config.mk
include $(SRC_DIR_REL)/msDd/config.mk

include $(SRC_DIR_REL)/msDd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/simp/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/msDd/sources.mk
include sources.mk

all: $(msDd_simulation_tests_OBJS) Test

clean:
	rm -f $(msDd_simulation_tests_OBJS)
	rm -f $(msDd_simulation_tests_OBJS:.o=.d)
	rm -f Test

clean-deps:
	rm -f $(msDd_simulation_tests_OBJS:.o=.d)

Test: Test.o $(msDd_LIB) $(ddMd_LIB) $(mcMd_LIB) $(simp_LIB) $(util_LIB) 
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(msDd_LIB) $(ddMd_LIB) $(mcMd_LIB) $(simp_LIB) $(util_LIB) 

-include $(msDd_simulation_tests_OBJS:.o=.d)

