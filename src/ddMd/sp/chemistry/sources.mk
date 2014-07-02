
ddMd_sp_chemistry_= \
    ddMd/sp/chemistry/SpMolecule.cpp \
    ddMd/sp/chemistry/SpSpecies.cpp 

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_chemistry_))
ddMd_sp_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_sp_chemistry_:.cpp=.o))

