include $(SRC_DIR)/ddMd/chemistry/sources.mk
include $(SRC_DIR)/ddMd/storage/sources.mk
include $(SRC_DIR)/ddMd/communicate/sources.mk
include $(SRC_DIR)/ddMd/neighbor/sources.mk
include $(SRC_DIR)/ddMd/ensembles/sources.mk
include $(SRC_DIR)/ddMd/system/sources.mk
include $(SRC_DIR)/ddMd/configIo/sources.mk
include $(SRC_DIR)/ddMd/potentials/sources.mk
include $(SRC_DIR)/ddMd/integrator/sources.mk
include $(SRC_DIR)/ddMd/util/sources.mk

ddMd_SRCS=$(ddMd_chemistry_SRCS) $(ddMd_storage_SRCS) \
    $(ddMd_communicate_SRCS) $(ddMd_neighbor_SRCS) $(ddMd_ensembles_SRCS) \
    $(ddMd_system_SRCS) $(ddMd_configIo_SRCS) \
    $(ddMd_potentials_SRCS) $(ddMd_interaction_SRCS) \
    $(ddMd_integrator_SRCS) \
    $(ddMd_util_SRCS)

ddMd_OBJS=$(ddMd_SRCS:.cpp=.o)

$(ddMd_LIB): $(ddMd_OBJS)
	$(AR) rcs $(ddMd_LIB) $(ddMd_OBJS)

