
mdCf_chemistry_= \
    mdCf/chemistry/Molecule.cpp \
    mdCf/chemistry/Species.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdCf_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_chemistry_))
mdCf_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_chemistry_:.cpp=.o))

