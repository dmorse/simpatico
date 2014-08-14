
spAn_chemistry_= \
    spAn/chemistry/Molecule.cpp \
    spAn/chemistry/Species.cpp 

# Create lists of source (*.cpp) and object (*.o) files
spAn_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_chemistry_))
spAn_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_chemistry_:.cpp=.o))

