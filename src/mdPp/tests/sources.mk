#include $(SRC_DIR)/mdPp/tests/chemistry/sources.mk
#include $(SRC_DIR)/mdPp/tests/storage/sources.mk
#include $(SRC_DIR)/mdPp/tests/config/sources.mk
#include $(SRC_DIR)/mdPp/tests/processor/sources.mk

mdPp_tests_= \
    mdPp/tests/Test.cc
    #$(mdPp_tests_chemistry_) \
    #$(mdPp_tests_storage_) \
    #$(mdPp_tests_config_) \
    #$(mdPp_tests_processor_) \

# Create lists of source (*.cc) and object (*.o) files
mdPp_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_))
mdPp_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_:.cc=.o))
mdPp_tests_EXES=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_:.cc=))

