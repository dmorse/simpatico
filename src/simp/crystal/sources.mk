
simp_crystal_=\
    simp/crystal/LatticeSystem.cpp \
    simp/crystal/PointGroup.cpp \
    simp/crystal/PointSymmetry.cpp 

simp_crystal_SRCS=$(addprefix $(SRC_DIR)/, $(simp_crystal_))
simp_crystal_OBJS=$(addprefix $(BLD_DIR)/, $(simp_crystal_:.cpp=.o))

