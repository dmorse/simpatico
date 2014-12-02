
tools_chemistry_= \
    tools/chemistry/Molecule.cpp \
    tools/chemistry/Species.cpp 

# Create lists of source (*.cpp) and object (*.o) files
tools_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_chemistry_))
tools_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_chemistry_:.cpp=.o))

