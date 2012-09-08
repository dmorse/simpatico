
util_misc_SRCS=$(SRC_DIR)/util/misc/Exception.cpp \
    $(SRC_DIR)/util/misc/initStatic.cpp $(SRC_DIR)/util/misc/ioUtil.cpp \
    $(SRC_DIR)/util/misc/Log.cpp 

util_misc_OBJS=$(util_misc_SRCS:.cpp=.o)

