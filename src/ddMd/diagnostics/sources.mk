ddMd_diagnostics_SRCS=\
     $(SRC_DIR)/ddMd/diagnostics/Diagnostic.cpp\
     $(SRC_DIR)/ddMd/diagnostics/DiagnosticManager.cpp\
     $(SRC_DIR)/ddMd/diagnostics/DiagnosticFactory.cpp\
     $(SRC_DIR)/ddMd/diagnostics/WriteConfig.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputEnergy.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputPressure.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputBoxdim.cpp \
     $(SRC_DIR)/ddMd/diagnostics/OutputTemperature.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputPairEnergies.cpp\
     $(SRC_DIR)/ddMd/diagnostics/StructureFactor.cpp\
     $(SRC_DIR)/ddMd/diagnostics/StructureFactorGrid.cpp\
     $(SRC_DIR)/ddMd/diagnostics/AsymmSF.cpp\
     $(SRC_DIR)/ddMd/diagnostics/AsymmSFGrid.cpp

ddMd_diagnostics_OBJS=$(ddMd_diagnostics_SRCS:.cpp=.o)
