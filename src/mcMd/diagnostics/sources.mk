include $(SRC_DIR)/mcMd/diagnostics/util/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/system/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mcSystem/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mdSystem/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/simulation/sources.mk
include $(SRC_DIR)/mcMd/diagnostics/mutable/sources.mk

mcMd_diagnostics_=\
    $(mcMd_diagnostics_util_) \
    $(mcMd_diagnostics_system_) \
    $(mcMd_diagnostics_mcSystem_) \
    $(mcMd_diagnostics_mdSystem_) \
    $(mcMd_diagnostics_simulation_) \
    $(mcMd_diagnostics_mutable_) \
    mcMd/diagnostics/Diagnostic.cpp \
    mcMd/diagnostics/DiagnosticManager.cpp 

ifdef MCMD_PERTURB
include $(SRC_DIR)/mcMd/diagnostics/perturb/sources.mk
mcMd_diagnostics_+=$(mcMd_diagnostics_perturb_)
endif

mcMd_diagnostics_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_))
mcMd_diagnostics_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_diagnostics_:.cpp=.o))

