
mdCf_storage_= \
    mdCf/storage/AtomStorage.cpp \
    mdCf/storage/System.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdCf_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mdCf_storage_))
mdCf_storage_OBJS=\
    $(addprefix $(OBJ_DIR)/, $(mdCf_storage_:.cpp=.o))

