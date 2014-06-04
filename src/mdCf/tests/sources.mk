include $(SRC_DIR)/mdCf/tests/chemistry/sources.mk
include $(SRC_DIR)/mdCf/tests/storage/sources.mk
include $(SRC_DIR)/mdCf/tests/configIos/sources.mk
include $(SRC_DIR)/mdCf/tests/processor/sources.mk

ddMd_tests_sp_= \
    $(ddMd_tests_sp_chemistry_) \
    $(ddMd_tests_sp_storage_) \
    $(ddMd_tests_sp_configIos_) \
    $(ddMd_tests_sp_processor_) 

# Create lists of source (*.cc) and object (*.o) files
ddMd_tests_sp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_))
ddMd_tests_sp_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_tests_sp_:.cc=.o))

