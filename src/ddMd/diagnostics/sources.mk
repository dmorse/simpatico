ddMd_diagnostics_SRCS=\
     $(SRC_DIR)/ddMd/diagnostics/Diagnostic.cpp\
     $(SRC_DIR)/ddMd/diagnostics/DiagnosticManager.cpp\
     $(SRC_DIR)/ddMd/diagnostics/DiagnosticFactory.cpp\
     $(SRC_DIR)/ddMd/diagnostics/WriteConfig.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputEnergy.cpp\
     $(SRC_DIR)/ddMd/diagnostics/OutputPressure.cpp

ddMd_diagnostics_OBJS=$(ddMd_diagnostics_SRCS:.cpp=.o)
