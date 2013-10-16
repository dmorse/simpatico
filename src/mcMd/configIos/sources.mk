mcMd_configIos_=mcMd/configIos/ConfigIo.cpp \
    mcMd/configIos/ConfigIoFactory.cpp \
    mcMd/configIos/DdMdConfigIo.cpp \
    mcMd/configIos/McConfigIo.cpp \
    mcMd/configIos/McMdConfigIo.cpp \
    mcMd/configIos/MdConfigIo.cpp \
    mcMd/configIos/PmcConfigIo.cpp \
    mcMd/configIos/LammpsConfigIo.cpp 

mcMd_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_configIos_))
mcMd_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_configIos_:.cpp=.o))

