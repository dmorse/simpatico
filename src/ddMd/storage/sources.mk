ddMd_storage_=\
   ddMd/storage/AtomMap.cpp \
   ddMd/storage/AtomStorage.cpp \
   ddMd/storage/GroupExchanger.cpp \
   ddMd/storage/BondStorage.cpp \

ifdef SIMP_ANGLE
ddMd_storage_+=ddMd/storage/AngleStorage.cpp
endif

ifdef SIMP_DIHEDRAL
ddMd_storage_+=ddMd/storage/DihedralStorage.cpp
endif

ddMd_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_storage_))
ddMd_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_storage_:.cpp=.o))

