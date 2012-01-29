ifdef PYTHON_FLAG

include $(SRC_DIR)/python/cpp/sources.mk

include $(SRC_DIR)/mcMd/defines.mk

python_SRCS=$(python_cpp_SRCS)

python_OBJS=$(python_SRCS:.cpp=.o)

$(python_LIB): $(python_OBJS)
	$(CXX) -shared -Wl,-soname,$(python_LIB) -Wl,--no-undefined -Wl,--whole-archive -o $(python_LIB) -L$(mcMd_LIBDIR) -l$(util_LIBNAME) -l$(mcMd_LIBNAME) -Wl,--no-whole-archive $(PYTHON_LIB) $(HOOMD_LIB) -lboost_python $(python_OBJS)

endif
