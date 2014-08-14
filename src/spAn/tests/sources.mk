#include $(SRC_DIR)/spAn/tests/chemistry/sources.mk
#include $(SRC_DIR)/spAn/tests/storage/sources.mk
#include $(SRC_DIR)/spAn/tests/configIos/sources.mk
#include $(SRC_DIR)/spAn/tests/processor/sources.mk

spAn_tests_= \
    spAn/tests/Test.cc
    #$(spAn_tests_chemistry_) \
    #$(spAn_tests_storage_) \
    #$(spAn_tests_configIos_) \
    #$(spAn_tests_processor_) \

# Create lists of source (*.cc) and object (*.o) files
spAn_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_))
spAn_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_:.cc=.o))

