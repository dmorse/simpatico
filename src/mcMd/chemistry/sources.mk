
mcMd_chemistry_=mcMd/chemistry/Atom.cpp \
    mcMd/chemistry/AtomType.cpp mcMd/chemistry/Mask.cpp \
    mcMd/chemistry/MaskPolicy.cpp 

mcMd_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_chemistry_))
mcMd_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_chemistry_:.cpp=.o))

