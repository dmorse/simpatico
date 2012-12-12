SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/math/sources.mk

all: $(util_math_OBJS)

clean:
	rm -f $(util_math_OBJS) $(util_math_OBJS:.o=.d)

clean-deps:
	rm -f $(util_math_OBJS:.o=.d)

-include $(util_math_OBJS:.o=.d)

