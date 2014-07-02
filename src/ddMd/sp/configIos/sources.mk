ddMd_sp_configIos_=\
   ddMd/sp/configIos/SpConfigIo.cpp \
   ddMd/sp/configIos/SpConfigIoFactory.cpp \
   ddMd/sp/configIos/DdMdSpConfigIo.cpp 

ddMd_sp_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_configIos_))
ddMd_sp_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_sp_configIos_:.cpp=.o))

