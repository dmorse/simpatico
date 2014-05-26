include $(SRC_DIR)/mdPp/tests/chemistry/sources.mk
include $(SRC_DIR)/mdPp/tests/storage/sources.mk
include $(SRC_DIR)/mdPp/tests/configIos/sources.mk
include $(SRC_DIR)/mdPp/tests/processor/sources.mk

mdPp_tests_= \
    $(mdPp_tests_chemistry_) \
    $(mdPp_tests_storage_) \
    $(mdPp_tests_configIos_) \
    $(mdPp_tests_processor_) 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_))
mdPp_tests_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_tests_:.cc=.o))

