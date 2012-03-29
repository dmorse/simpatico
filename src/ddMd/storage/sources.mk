
ddMd_storage_SRCS=$(SRC_DIR)/ddMd/storage/AtomStorage.cpp \
    $(SRC_DIR)/ddMd/storage/BondStorage.cpp \
    $(SRC_DIR)/ddMd/storage/AngleStorage.cpp \
    $(SRC_DIR)/ddMd/storage/DihedralStorage.cpp 

ddMd_storage_OBJS=$(ddMd_storage_SRCS:.cpp=.o)

