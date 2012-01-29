
util_util_SRCS=$(SRC_DIR)/util/util/Exception.cpp \
    $(SRC_DIR)/util/util/initStatic.cpp $(SRC_DIR)/util/util/ioUtil.cpp \
    $(SRC_DIR)/util/util/Log.cpp 

util_util_OBJS=$(util_util_SRCS:.cpp=.o)

