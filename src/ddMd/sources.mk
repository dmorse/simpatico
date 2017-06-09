# Include source files lists from subdirectories
include $(SRC_DIR)/ddMd/chemistry/sources.mk
include $(SRC_DIR)/ddMd/storage/sources.mk
include $(SRC_DIR)/ddMd/communicate/sources.mk
include $(SRC_DIR)/ddMd/neighbor/sources.mk
include $(SRC_DIR)/ddMd/simulation/sources.mk
include $(SRC_DIR)/ddMd/configIos/sources.mk
include $(SRC_DIR)/ddMd/potentials/sources.mk
include $(SRC_DIR)/ddMd/integrators/sources.mk
include $(SRC_DIR)/ddMd/analyzers/sources.mk
include $(SRC_DIR)/ddMd/misc/sources.mk

# Concatenate source file lists from subdirectories
ddMd_=$(ddMd_chemistry_) $(ddMd_storage_) \
    $(ddMd_communicate_) $(ddMd_neighbor_) \
    $(ddMd_simulation_) $(ddMd_configIos_) \
    $(ddMd_potentials_) $(ddMd_integrators_) \
    $(ddMd_analyzers_) $(ddMd_misc_) 

ifdef DDMD_MODIFIERS
include $(SRC_DIR)/ddMd/modifiers/sources.mk
ddMd_+= $(ddMd_modifiers_) 
endif

# Add user source files, if any
include $(SRC_DIR)/ddMd/user/sources.mk
ddMd_+= $(ddMd_user_) 

# Create lists of src and object files, with absolute paths
ddMd_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_))
ddMd_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_:.cpp=.o))

# Target to create static library file for DdMd namespace
$(ddMd_LIB): $(ddMd_OBJS)
	$(AR) rcs $(ddMd_LIB) $(ddMd_OBJS)

