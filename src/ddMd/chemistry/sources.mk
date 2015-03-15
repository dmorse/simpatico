ddMd_chemistry_=\
    ddMd/chemistry/Atom.cpp \
    ddMd/chemistry/AtomArray.cpp \
    ddMd/chemistry/AtomType.cpp \
    ddMd/chemistry/Mask.cpp \
    ddMd/chemistry/MaskPolicy.cpp 

ddMd_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_chemistry_))
ddMd_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_chemistry_:.cpp=.o))

