
spAn_storage_= \
    spAn/storage/AtomStorage.cpp \
    spAn/storage/Configuration.cpp 

# Create lists of source (*.cpp) and object (*.o) files
spAn_storage_SRCS=\
    $(addprefix $(SRC_DIR)/, $(spAn_storage_))
spAn_storage_OBJS=\
    $(addprefix $(BLD_DIR)/, $(spAn_storage_:.cpp=.o))

