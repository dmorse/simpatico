
mcMd_chemistry_SRCS=$(SRC_DIR)/mcMd/chemistry/Atom.cpp \
    $(SRC_DIR)/mcMd/chemistry/AtomType.cpp $(SRC_DIR)/mcMd/chemistry/Mask.cpp \
    $(SRC_DIR)/mcMd/chemistry/MaskPolicy.cpp 

mcMd_chemistry_OBJS=$(mcMd_chemistry_SRCS:.cpp=.o)

