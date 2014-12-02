
tools_storage_= \
    tools/storage/AtomStorage.cpp \
    tools/storage/Configuration.cpp 

# Create lists of source (*.cpp) and object (*.o) files
tools_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(tools_storage_))
tools_storage_OBJS=\
    $(addprefix $(BLD_DIR)/, $(tools_storage_:.cpp=.o))

