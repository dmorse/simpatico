SRC_DIR =../..

include $(SRC_DIR)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/python/config.mk
include $(SRC_DIR)/python/patterns.mk
include $(SRC_DIR)/python/cppmodule/sources.mk

all: $(python_cpp_OBJS)

clean:
	rm -f $(python_cpp_OBJS) $(python_cpp_OBJS:.o=.d)

clean-deps:
	rm -f $(python_cpp_OBJS:.o=.d)

-include $(python_cpp_OBJS:.o=.d)

