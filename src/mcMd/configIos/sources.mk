mcMd_configIos_SRCS=$(SRC_DIR)/mcMd/configIos/ConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/ConfigIoFactory.cpp \
    $(SRC_DIR)/mcMd/configIos/DdMdConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/LammpsConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/McConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/McMdConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/MdConfigIo.cpp \
    $(SRC_DIR)/mcMd/configIos/PmcConfigIo.cpp 

mcMd_configIos_OBJS=$(mcMd_configIos_SRCS:.cpp=.o)

