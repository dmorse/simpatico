include $(SRC_DIR)/ddMd/chemistry/sources.mk
include $(SRC_DIR)/ddMd/storage/sources.mk
include $(SRC_DIR)/ddMd/communicate/sources.mk
include $(SRC_DIR)/ddMd/neighbor/sources.mk
include $(SRC_DIR)/ddMd/simulation/sources.mk
include $(SRC_DIR)/ddMd/configIos/sources.mk
include $(SRC_DIR)/ddMd/potentials/sources.mk
include $(SRC_DIR)/ddMd/integrators/sources.mk
include $(SRC_DIR)/ddMd/diagnostics/sources.mk
include $(SRC_DIR)/ddMd/misc/sources.mk
include $(SRC_DIR)/ddMd/user/sources.mk

ddMd_SRCS=$(ddMd_chemistry_SRCS) $(ddMd_storage_SRCS) \
    $(ddMd_communicate_SRCS) $(ddMd_neighbor_SRCS) \
    $(ddMd_simulation_SRCS) $(ddMd_configIos_SRCS) \
    $(ddMd_potentials_SRCS) $(ddMd_integrators_SRCS) \
    $(ddMd_diagnostics_SRCS) \
    $(ddMd_misc_SRCS) \
    $(ddMd_user_SRCS) 

ddMd_OBJS=$(ddMd_SRCS:.cpp=.o)

$(ddMd_LIB): $(ddMd_OBJS)
	$(AR) rcs $(ddMd_LIB) $(ddMd_OBJS)

