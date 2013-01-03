inter_dihedral_SRCS=\
    $(SRC_DIR)/inter/dihedral/CosineDihedral.cpp \
    $(SRC_DIR)/inter/dihedral/MultiHarmonicDihedral.cpp 

inter_dihedral_OBJS=$(inter_dihedral_SRCS:.cpp=.o)

