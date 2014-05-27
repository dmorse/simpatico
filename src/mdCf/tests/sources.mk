include $(SRC_DIR)/mdCf/tests/chemistry/sources.mk
include $(SRC_DIR)/mdCf/tests/storage/sources.mk
include $(SRC_DIR)/mdCf/tests/configIos/sources.mk
include $(SRC_DIR)/mdCf/tests/processor/sources.mk

mdCf_tests_= \
    $(mdCf_tests_chemistry_) \
    $(mdCf_tests_storage_) \
    $(mdCf_tests_configIos_) \
    $(mdCf_tests_processor_) 

# Create lists of source (*.cc) and object (*.o) files
mdCf_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_tests_))
mdCf_tests_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_tests_:.cc=.o))

