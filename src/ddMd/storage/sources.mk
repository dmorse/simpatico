
ddMd_storage_SRCS=$(SRC_DIR)/ddMd/storage/AtomStorage.cpp \
                  $(SRC_DIR)/ddMd/storage/GroupExchanger.cpp \
                  $(SRC_DIR)/ddMd/storage/BondStorage.cpp \

ifdef INTER_ANGLE
ddMd_storage_SRCS+=$(SRC_DIR)/ddMd/storage/AngleStorage.cpp
endif

ifdef INTER_DIHEDRAL
ddMd_storage_SRCS+=$(SRC_DIR)/ddMd/storage/DihedralStorage.cpp
endif

ddMd_storage_OBJS=$(ddMd_storage_SRCS:.cpp=.o)
