mcMd_diagnostics_mdSystem_SRCS=\
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdDiagnosticFactory.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdEnergyOutput.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdKineticEnergyAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdPairEnergyCoefficients.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdPotentialEnergyAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdPressureAverage.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdIntraBondStressAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/mdSystem/MdIntraBondTensorAutoCorr.cpp 

mcMd_diagnostics_mdSystem_OBJS=$(mcMd_diagnostics_mdSystem_SRCS:.cpp=.o)

