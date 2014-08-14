include $(SRC_DIR)/spAn/chemistry/sources.mk
include $(SRC_DIR)/spAn/storage/sources.mk
include $(SRC_DIR)/spAn/configIos/sources.mk
include $(SRC_DIR)/spAn/neighbor/sources.mk
include $(SRC_DIR)/spAn/processor/sources.mk
include $(SRC_DIR)/spAn/analyzers/sources.mk

spAn_= \
    $(spAn_chemistry_) \
    $(spAn_storage_) \
    $(spAn_configIos_) \
    $(spAn_neighbor_) \
    $(spAn_processor_) \
    $(spAn_analyzers_) 

# Create lists of source (*.cpp) and object (*.o) files
spAn_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_))
spAn_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_:.cpp=.o))

$(spAn_LIB): $(spAn_OBJS)
	$(AR) rcs $(spAn_LIB) $(spAn_OBJS)

