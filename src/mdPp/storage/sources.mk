
mdPp_storage_= \
    mdPp/storage/AtomStorage.cpp \
    mdPp/storage/SpeciesStorage.cpp \
    mdPp/storage/Configuration.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mdPp_storage_))
mdPp_storage_OBJS=\
    $(addprefix $(BLD_DIR)/, $(mdPp_storage_:.cpp=.o))

