
util_space_SRCS=$(SRC_DIR)/util/space/Grid.cpp \
    $(SRC_DIR)/util/space/IntVector.cpp $(SRC_DIR)/util/space/Tensor.cpp \
    $(SRC_DIR)/util/space/Vector.cpp 

util_space_OBJS=$(util_space_SRCS:.cpp=.o)

