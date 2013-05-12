
ddMd_configIos_SRCS=\
   $(SRC_DIR)/ddMd/configIos/ConfigIo.cpp \
   $(SRC_DIR)/ddMd/configIos/DdMdConfigIo.cpp \
   $(SRC_DIR)/ddMd/configIos/DdMdOrderedConfigIo.cpp \
   $(SRC_DIR)/ddMd/configIos/LammpsConfigIo.cpp \
   $(SRC_DIR)/ddMd/configIos/SerializeConfigIo.cpp \
   $(SRC_DIR)/ddMd/configIos/ConfigIoFactory.cpp 

ddMd_configIos_OBJS=$(ddMd_configIos_SRCS:.cpp=.o)

