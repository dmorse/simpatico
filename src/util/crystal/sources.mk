
util_crystal_SRCS=$(SRC_DIR)/util/crystal/LatticeSystem.cpp \
    $(SRC_DIR)/util/crystal/PointGroup.cpp \
    $(SRC_DIR)/util/crystal/PointSymmetry.cpp 

util_crystal_OBJS=$(util_crystal_SRCS:.cpp=.o)

