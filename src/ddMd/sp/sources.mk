include $(SRC_DIR)/ddMd/sp/chemistry/sources.mk
include $(SRC_DIR)/ddMd/sp/storage/sources.mk
include $(SRC_DIR)/ddMd/sp/configIos/sources.mk
include $(SRC_DIR)/ddMd/sp/processor/sources.mk
include $(SRC_DIR)/ddMd/sp/analyzers/sources.mk

ddMd_sp_= \
    $(ddMd_sp_chemistry_) \
    $(ddMd_sp_storage_) \
    $(ddMd_sp_configIos_) \
    $(ddMd_sp_processor_) \
    $(ddMd_sp_analyzers_) 

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_))
ddMd_sp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_sp_:.cpp=.o))

