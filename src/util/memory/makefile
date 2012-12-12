SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/memory/sources.mk

all: $(util_memory_OBJS)

clean:
	rm -f $(util_memory_OBJS) $(util_memory_OBJS:.o=.d)

clean-deps:
	rm -f $(util_memory_OBJS:.o=.d)

-include $(util_memory_OBJS:.o=.d)

