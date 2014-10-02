include $(SRC_DIR)/mcMd/analyzers/util/sources.mk
include $(SRC_DIR)/mcMd/analyzers/system/sources.mk
include $(SRC_DIR)/mcMd/analyzers/mcSystem/sources.mk
include $(SRC_DIR)/mcMd/analyzers/mdSystem/sources.mk
include $(SRC_DIR)/mcMd/analyzers/simulation/sources.mk
include $(SRC_DIR)/mcMd/analyzers/mutable/sources.mk

mcMd_analyzers_=\
    $(mcMd_analyzers_util_) \
    $(mcMd_analyzers_system_) \
    $(mcMd_analyzers_mcSystem_) \
    $(mcMd_analyzers_mdSystem_) \
    $(mcMd_analyzers_simulation_) \
    $(mcMd_analyzers_mutable_) \
    mcMd/analyzers/Analyzer.cpp \
    mcMd/analyzers/AnalyzerManager.cpp 

ifdef MCMD_PERTURB
include $(SRC_DIR)/mcMd/analyzers/perturb/sources.mk
mcMd_analyzers_+=$(mcMd_analyzers_perturb_)
endif

mcMd_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_))
mcMd_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_:.cpp=.o))

