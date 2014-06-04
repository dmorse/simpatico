ddMd_sp_configIos_=\
   mdCf/configIos/SpConfigIo.cpp \
   mdCf/configIos/SpConfigIoFactory.cpp \
   mdCf/configIos/DdMdSpConfigIo.cpp 

ddMd_sp_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_configIos_))
ddMd_sp_configIos_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_sp_configIos_:.cpp=.o))

