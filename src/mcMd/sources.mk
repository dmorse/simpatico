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
include $(SRC_DIR)/mcMd/analyzers/sources.mk
include $(SRC_DIR)/mcMd/user/sources.mk

mcMd_=\
    $(mcMd_chemistry_) \
    $(mcMd_species_) \
    $(mcMd_neighbor_) \
    $(mcMd_potentials_) \
    $(mcMd_simulation_) \
    $(mcMd_mcSimulation_) \
    $(mcMd_mdSimulation_) \
    $(mcMd_configIos_) \
    $(mcMd_trajectoryIos_) \
    $(mcMd_mdIntegrators_) \
    $(mcMd_mcMoves_) \
    $(mcMd_analyzers_) \
    $(mcMd_user_)

ifdef MCMD_PERTURB
include $(SRC_DIR)/mcMd/perturb/sources.mk
mcMd_+=$(mcMd_perturb_)
endif

ifdef MCMD_LINK
include $(SRC_DIR)/mcMd/links/sources.mk
mcMd_+=$(mcMd_links_)
endif

ifdef INTER_TETHER
include $(SRC_DIR)/mcMd/tethers/sources.mk
mcMd_+=$(mcMd_tethers_)
endif

mcMd_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_))
mcMd_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_:.cpp=.o))

$(mcMd_LIB): $(mcMd_OBJS)
	$(AR) rcs $(mcMd_LIB) $(mcMd_OBJS)
