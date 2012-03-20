SRC_DIR =../..

include $(SRC_DIR)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/python/defines.mk
include $(SRC_DIR)/python/patterns.mk
include $(SRC_DIR)/python/cppmodule/sources.mk

all: $(python_cpp_OBJS)

clean:
	rm -f $(python_cpp_OBJS) $(python_cpp_OBJS:.o=.d)

clean-deps:
	rm -f $(python_cpp_OBJS:.o=.d)

-include $(python_cpp_OBJS:.o=.d)

