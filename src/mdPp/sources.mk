# Include source files lists from subdirectories
include $(SRC_DIR)/mdPp/chemistry/sources.mk
include $(SRC_DIR)/mdPp/storage/sources.mk
include $(SRC_DIR)/mdPp/config/sources.mk
include $(SRC_DIR)/mdPp/trajectory/sources.mk
include $(SRC_DIR)/mdPp/neighbor/sources.mk
include $(SRC_DIR)/mdPp/processor/sources.mk
include $(SRC_DIR)/mdPp/analyzers/sources.mk
include $(SRC_DIR)/mdPp/user/sources.mk

# Concatenate source file lists from subdirectories
mdPp_= \
    $(mdPp_chemistry_) \
    $(mdPp_storage_) \
    $(mdPp_config_) \
    $(mdPp_trajectory_) \
    $(mdPp_neighbor_) \
    $(mdPp_processor_) \
    $(mdPp_analyzers_) \
    $(mdPp_user_) 

# Create lists of src (*.cpp) and object (*.o) files, 
# with absolute paths to each file
mdPp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_))
mdPp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_:.cpp=.o))

# Target to create library file for MdPp namespace
$(mdPp_LIB): $(mdPp_OBJS)
	$(AR) rcs $(mdPp_LIB) $(mdPp_OBJS)

