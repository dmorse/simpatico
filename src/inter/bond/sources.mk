inter_bond_SRCS=\
    $(SRC_DIR)/inter/bond/FeneBond.cpp \
    $(SRC_DIR)/inter/bond/HarmonicBond.cpp \
    $(SRC_DIR)/inter/bond/HarmonicL0Bond.cpp 

inter_bond_OBJS=$(inter_bond_SRCS:.cpp=.o)

