mcMd_diagnostics_mcSystem_SRCS=\
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McBondEnergyAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McDiagnosticFactory.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McEnergyAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McEnergyOutput.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McPairEnergyAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McPressureAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McStressAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McIntraBondStressAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McIntraBondTensorAutoCorr.cpp 

ifdef INTER_EXTERNAL
mcMd_diagnostics_mcSystem_SRCS+=\
    $(SRC_DIR)/mcMd/diagnostics/mcSystem/McExternalEnergyAverage.cpp 
endif

mcMd_diagnostics_mcSystem_OBJS=$(mcMd_diagnostics_mcSystem_SRCS:.cpp=.o)

