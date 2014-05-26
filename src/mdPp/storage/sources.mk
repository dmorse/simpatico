
mdPp_storage_= \
    mdPp/storage/Storage.cpp \
    mdPp/storage/AtomStorage.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mdPp_storage_))
mdPp_storage_OBJS=\
    $(addprefix $(OBJ_DIR)/, $(mdPp_storage_:.cpp=.o))

