include $(SRC_DIR)/mcMd/chemistry/sources.mk
include $(SRC_DIR)/mcMd/species/sources.mk
include $(SRC_DIR)/mcMd/neighbor/sources.mk
include $(SRC_DIR)/mcMd/potentials/sources.mk
include $(SRC_DIR)/mcMd/simulation/sources.mk
include $(SRC_DIR)/mcMd/mcSimulation/sources.mk
include $(SRC_DIR)/mcMd/mdSimulation/sources.mk
include $(SRC_DIR)/mcMd/configIos/sources.mk
include $(SRC_DIR)/mcMd/trajectoryIos/sources.mk
include $(SRC_DIR)/mcMd/mdIntegrators/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/sources.mk
include $(SRC_DIR)/mcMd/user/sources.mk

mcMd_SRCS=\
    $(mcMd_potentials_SRCS) $(mcMd_chemistry_SRCS) \
    $(mcMd_species_SRCS) $(mcMd_neighbor_SRCS) \
    $(mcMd_simulation_SRCS) $(mcMd_mcSimulation_SRCS) \
    $(mcMd_mdSimulation_SRCS) $(mcMd_configIos_SRCS) \
    $(mcMd_trajectoryIos_SRCS) $(mcMd_mdIntegrators_SRCS) \
    $(mcMd_mcMoves_SRCS) $(mcMd_diagnostics_SRCS) \
    $(mcMd_user_SRCS)


ifdef MCMD_PERTURB
include $(SRC_DIR)/mcMd/perturb/sources.mk
mcMd_SRCS+=$(mcMd_perturb_SRCS)
endif

ifdef MCMD_LINK
include $(SRC_DIR)/mcMd/links/sources.mk
mcMd_SRCS+=$(mcMd_links_SRCS)
endif

ifdef INTER_TETHER
include $(SRC_DIR)/mcMd/tethers/sources.mk
mcMd_SRCS+=$(mcMd_tethers_SRCS)
endif

mcMd_OBJS=$(mcMd_SRCS:.cpp=.o)

$(mcMd_LIB): $(mcMd_OBJS)
	$(AR) rcs $(mcMd_LIB) $(mcMd_OBJS)

