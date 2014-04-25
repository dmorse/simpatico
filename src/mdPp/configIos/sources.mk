mdPp_configIos_=\
   mdPp/configIos/ConfigIo.cpp \
   mdPp/configIos/ConfigIoFactory.cpp \
   mdPp/configIos/DdMdConfigIo.cpp 

mdPp_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_configIos_))
mdPp_configIos_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_configIos_:.cpp=.o))

