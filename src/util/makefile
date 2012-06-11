SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

all: $(util_OBJS) $(util_LIB)

clean:
	rm -f $(util_OBJS) $(util_OBJS:.o=.d) $(util_LIB)
	rm -f $(util_LIB)
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

clean-deps:
	rm -f $(util_OBJS:.o=.d)

-include $(util_OBJS:.o=.d)

