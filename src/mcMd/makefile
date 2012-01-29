SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: $(mcMd_OBJS) $(mcMd_LIB)

clean:
	rm -f $(mcMd_OBJS) $(mcMd_LIB)
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f Test Test.o Test.d

clean-deps:
	rm -f $(mcMd_OBJS:.o=.d)

clean-test:
	rm -f Test Test.o Test.d

Test: $(mcMd_LIB) $(SRC_DIR)/mcMd/Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test $(SRC_DIR)/mcMd/Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include $(mcMd_OBJS:.o=.d)
-include Test.d

