inter_dihedral_=\
    inter/dihedral/CosineDihedral.cpp \
    inter/dihedral/MultiHarmonicDihedral.cpp 

inter_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_dihedral_))
inter_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_dihedral_:.cpp=.o))

