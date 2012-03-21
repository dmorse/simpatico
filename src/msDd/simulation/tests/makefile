SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/msDd/defines.mk

include $(SRC_DIR_REL)/msDd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/msDd/sources.mk
include sources.mk

all: $(msDd_simulation_tests_OBJS) Test

clean:
	rm -f $(msDd_simulation_tests_OBJS) Test $(msDd_simulation_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(msDd_simulation_tests_OBJS:.o=.d)

Test: Test.o $(msDd_LIB) $(ddMd_LIB) $(mcMd_LIB) $(inter_LIB) $(util_LIB) 
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(msDd_LIB) $(ddMd_LIB) $(mcMd_LIB) $(inter_LIB) $(util_LIB) 

-include $(msDd_simulation_tests_OBJS:.o=.d)

