spAn_configIos_=\
   spAn/configIos/ConfigIo.cpp \
   spAn/configIos/ConfigIoFactory.cpp \
   spAn/configIos/DdMdConfigIo.cpp 

spAn_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_configIos_))
spAn_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_configIos_:.cpp=.o))

