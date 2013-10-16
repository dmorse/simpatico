ddMd_configIos_=\
   ddMd/configIos/ConfigIo.cpp \
   ddMd/configIos/DdMdConfigIo.cpp \
   ddMd/configIos/DdMdOrderedConfigIo.cpp \
   ddMd/configIos/LammpsConfigIo.cpp \
   ddMd/configIos/SerializeConfigIo.cpp \
   ddMd/configIos/ConfigIoFactory.cpp 

ddMd_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_configIos_))
ddMd_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_configIos_:.cpp=.o))

