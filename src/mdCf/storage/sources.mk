
ddMd_sp_storage_= \
    mdCf/storage/SpAtomStorage.cpp \
    mdCf/storage/SpConfiguration.cpp 

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(ddMd_sp_storage_))
ddMd_sp_storage_OBJS=\
    $(addprefix $(OBJ_DIR)/, $(ddMd_sp_storage_:.cpp=.o))

