
ddMd_chemistry_SRCS=\
    $(SRC_DIR)/ddMd/chemistry/Atom.cpp \
    $(SRC_DIR)/ddMd/chemistry/AtomArray.cpp \
    $(SRC_DIR)/ddMd/chemistry/AtomType.cpp \
    $(SRC_DIR)/ddMd/chemistry/Mask.cpp \
    $(SRC_DIR)/ddMd/chemistry/MaskPolicy.cpp 

ddMd_chemistry_OBJS=$(ddMd_chemistry_SRCS:.cpp=.o)

