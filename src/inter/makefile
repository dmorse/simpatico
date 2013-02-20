SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/sources.mk

all: $(inter_OBJS) $(inter_LIB)

clean:
	rm -f $(inter_OBJS) $(inter_LIB)
	rm -f $(inter_OBJS:.o=.d)
	rm -f */*.o */*.d

clean-deps:
	rm -f $(inter_OBJS:.o=.d)
	rm -f */*.d

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)


