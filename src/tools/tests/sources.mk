#include $(SRC_DIR)/tools/tests/chemistry/sources.mk
#include $(SRC_DIR)/tools/tests/storage/sources.mk
#include $(SRC_DIR)/tools/tests/config/sources.mk
#include $(SRC_DIR)/tools/tests/processor/sources.mk

tools_tests_= \
    tools/tests/Test.cc
    #$(tools_tests_chemistry_) \
    #$(tools_tests_storage_) \
    #$(tools_tests_config_) \
    #$(tools_tests_processor_) \

# Create lists of source (*.cc) and object (*.o) files
tools_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_))
tools_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_:.cc=.o))
tools_tests_EXES=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_:.cc=))

