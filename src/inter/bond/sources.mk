inter_bond_=\
    inter/bond/FeneBond.cpp \
    inter/bond/HarmonicBond.cpp \
    inter/bond/HarmonicL0Bond.cpp 

inter_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_bond_))
inter_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_bond_:.cpp=.o))

