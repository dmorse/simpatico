include $(SRC_DIR)/mcMd/diagnostics/util/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/system/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mcSystem/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mdSystem/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/simulation/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mutable/sources.mk

mcMd_diagnostics_SRCS=$(mcMd_diagnostics_util_SRCS) \
    $(mcMd_diagnostics_system_SRCS) $(mcMd_diagnostics_mcSystem_SRCS) \
    $(mcMd_diagnostics_mdSystem_SRCS) $(mcMd_diagnostics_simulation_SRCS) \
    $(mcMd_diagnostics_mutable_SRCS) \
    $(SRC_DIR)/mcMd/diagnostics/Diagnostic.cpp \
    $(SRC_DIR)/mcMd/diagnostics/DiagnosticManager.cpp 

ifdef MCMD_PERTURB
include $(SRC_DIR)/mcMd/diagnostics/perturb/sources.mk
mcMd_diagnostics_SRCS+=$(mcMd_diagnostics_perturb_SRCS)
endif

mcMd_diagnostics_OBJS=$(mcMd_diagnostics_SRCS:.cpp=.o)

