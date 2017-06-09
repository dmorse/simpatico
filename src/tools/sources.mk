# Include source files lists from subdirectories
include $(SRC_DIR)/tools/chemistry/sources.mk
include $(SRC_DIR)/tools/storage/sources.mk
include $(SRC_DIR)/tools/config/sources.mk
include $(SRC_DIR)/tools/trajectory/sources.mk
include $(SRC_DIR)/tools/neighbor/sources.mk
include $(SRC_DIR)/tools/processor/sources.mk
include $(SRC_DIR)/tools/analyzers/sources.mk
include $(SRC_DIR)/tools/user/sources.mk

# Concatenate source file lists from subdirectories
tools_= \
    $(tools_chemistry_) \
    $(tools_storage_) \
    $(tools_config_) \
    $(tools_trajectory_) \
    $(tools_neighbor_) \
    $(tools_processor_) \
    $(tools_analyzers_) \
    $(tools_user_) 

# Create lists of src (*.cpp) and object (*.o) files, 
# with absolute paths to each file
tools_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_))
tools_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_:.cpp=.o))

# Target to create library file for Tools namespace
$(tools_LIB): $(tools_OBJS)
	$(AR) rcs $(tools_LIB) $(tools_OBJS)

