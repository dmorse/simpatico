
util_space_=util/space/Grid.cpp \
    util/space/IntVector.cpp util/space/Tensor.cpp \
    util/space/Vector.cpp 

util_space_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_space_))
util_space_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_space_:.cpp=.o))

