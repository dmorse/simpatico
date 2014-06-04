include $(SRC_DIR)/mdCf/chemistry/sources.mk
include $(SRC_DIR)/mdCf/storage/sources.mk
include $(SRC_DIR)/mdCf/configIos/sources.mk
include $(SRC_DIR)/mdCf/analyzers/sources.mk
include $(SRC_DIR)/mdCf/processor/sources.mk

ddMd_sp_= \
    $(ddMd_sp_chemistry_) \
    $(ddMd_sp_storage_) \
    $(ddMd_sp_configIos_) \
    $(ddMd_sp_analyzers_) \
    $(ddMd_sp_processor_) 

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_))
ddMd_sp_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_sp_:.cpp=.o))

$(ddMd_sp_LIB): $(ddMd_sp_OBJS)
	$(AR) rcs $(ddMd_sp_LIB) $(ddMd_sp_OBJS)

