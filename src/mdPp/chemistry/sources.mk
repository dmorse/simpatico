
mdPp_chemistry_= \
    mdPp/chemistry/Molecule.cpp \
    mdPp/chemistry/Species.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_chemistry_))
mdPp_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_chemistry_:.cpp=.o))

