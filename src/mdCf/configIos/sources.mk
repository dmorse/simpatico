mdCf_configIos_=\
   mdCf/configIos/ConfigIo.cpp \
   mdCf/configIos/ConfigIoFactory.cpp \
   mdCf/configIos/DdMdConfigIo.cpp 

mdCf_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_configIos_))
mdCf_configIos_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_configIos_:.cpp=.o))

